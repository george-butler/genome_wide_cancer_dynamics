#!/usr/bin/env python3
"""
Filter miniprot GFF3 to:
  1) keep the best mRNA per Target gene (by Rank, then F1),
  2) remove same-strand overlapping mRNAs within loci (keep best by Rank, then F1),
  3) allow overlapping models on opposite strands,
  4) drop frameshift models,
  5) optionally write:
        - a TSV of mRNAs removed due to same-strand overlaps,
        - a TSV summarizing all candidate mRNAs (kept + removed).

F1 = harmonic mean of Identity and Positive:
    F1 = 2 * Identity * Positive / (Identity + Positive)

Typical usage (matching your AWK thresholds):
    python filter_miniprot_best_per_target_strandaware.py \
        miniprot_fixed.gff3 \
        miniprot.filtered.nonoverlap.gff3 \
        --min_identity 0.75 \
        --min_positive 0.75 \
        --removed_tsv overlapping_removed.tsv \
        --summary_tsv all_models_summary.tsv
"""

import sys
import argparse


def parse_attrs(attr_str):
    """Parse GFF3 attributes into a dict."""
    attrs = {}
    for part in attr_str.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
        else:
            attrs[part] = True
    return attrs


def extract_target_gene(target_val):
    """
    miniprot Target attribute is usually:
      Target=geneID start end [strand]
    We only want the first token (geneID).
    """
    if not target_val:
        return None
    return target_val.split()[0]


def pick_best(candidates):
    """
    Pick best record based on:
      1) lowest Rank
      2) highest F1 if Rank ties
    Each candidate is a dict with keys: 'rank', 'f1'.
    """
    best = None
    for rec in candidates:
        if best is None:
            best = rec
        else:
            if (rec["rank"] < best["rank"]) or (
                rec["rank"] == best["rank"] and rec["f1"] > best["f1"]
            ):
                best = rec
    return best


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Filter miniprot GFF3 to best per Target gene and remove same-strand "
            "overlapping loci using Rank and F1."
        )
    )
    ap.add_argument("gff_in", help="Input miniprot GFF3 (from miniprot)")
    ap.add_argument("gff_out", help="Output filtered non-overlapping GFF3")
    ap.add_argument(
        "--min_identity",
        type=float,
        default=0.0,
        help="Minimum Identity to keep (>=). Default: 0.0",
    )
    ap.add_argument(
        "--min_positive",
        type=float,
        default=0.0,
        help="Minimum Positive to keep (>). Default: 0.0",
    )
    ap.add_argument(
        "--removed_tsv",
        help="Optional TSV listing mRNAs removed due to same-strand overlaps.",
        default=None,
    )
    ap.add_argument(
        "--summary_tsv",
        help="Optional TSV summarizing ALL candidate mRNAs (kept + removed).",
        default=None,
    )
    args = ap.parse_args()

    # First pass: collect all mRNA models passing basic filters
    # Group by Target gene so we can pick best per Target.
    by_target = {}           # target_gene -> list of recs
    all_recs_by_mrna = {}    # mrna_id -> rec (for later lookup)

    with open(args.gff_in) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attr_str = cols

            if ftype != "mRNA":
                continue

            attrs = parse_attrs(attr_str)

            # Skip frameshift models entirely
            if "frameshift" in attrs:
                continue

            # Identity / Positive / Rank / Target / ID
            identity = float(attrs.get("Identity", "0"))
            positive = float(attrs.get("Positive", "0"))
            rank_str = attrs.get("Rank", "999999")
            mrna_id = attrs.get("ID")
            target_raw = attrs.get("Target", "")
            target_gene = extract_target_gene(target_raw)

            if mrna_id is None or target_gene is None:
                continue

            # Apply thresholds (like your AWK: $2 >= min_identity && $3 > min_positive)
            if identity < args.min_identity or positive <= args.min_positive:
                continue

            # Compute F1-like score
            if identity + positive > 0:
                f1 = (2.0 * identity * positive) / (identity + positive)
            else:
                f1 = 0.0

            try:
                rank = float(rank_str)
            except ValueError:
                rank = 999999.0

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            parent_gene = attrs.get("Parent")  # gene ID if present

            rec = {
                "seqid": seqid,
                "source": source,
                "type": ftype,
                "start": start_i,
                "end": end_i,
                "strand": strand,
                "attr_str": attr_str,
                "attrs": attrs,
                "identity": identity,
                "positive": positive,
                "rank": rank,
                "f1": f1,
                "mrna_id": mrna_id,
                "target_gene": target_gene,
                "parent_gene": parent_gene,
            }

            all_recs_by_mrna[mrna_id] = rec

            if target_gene not in by_target:
                by_target[target_gene] = []
            by_target[target_gene].append(rec)

    # For each Target gene, keep only the best mRNA (by Rank, then F1).
    # Record which ones were removed at this stage.
    best_per_target = []
    removed_by_target = {}  # mrna_id -> dict(reason, winner_mrna, winner_target)

    for tgt, recs in by_target.items():
        winner = pick_best(recs)
        best_per_target.append(winner)
        for r in recs:
            if r["mrna_id"] == winner["mrna_id"]:
                continue
            removed_by_target[r["mrna_id"]] = {
                "reason": "worse_than_other_for_same_target",
                "winner_mrna": winner["mrna_id"],
                "winner_target": tgt,
            }

    # Now enforce same-strand non-overlap among these winners:
    # Group best_per_target by (seqid, strand), then cluster by overlapping loci.
    # Within each overlapping cluster, keep best (Rank, then F1), remove others.
    groups = {}  # (seqid, strand) -> list of recs
    for rec in best_per_target:
        key = (rec["seqid"], rec["strand"])
        if key not in groups:
            groups[key] = []
        groups[key].append(rec)

    kept_mrna_ids = set()
    removed_overlap = {}  # mrna_id -> dict(reason, winner_mrna, winner_target)

    removed_fh = None
    if args.removed_tsv is not None:
        removed_fh = open(args.removed_tsv, "w")
        header_cols = [
            "TargetGene",
            "mRNA_ID",
            "Seqid",
            "Start",
            "End",
            "Strand",
            "Identity",
            "Positive",
            "Rank",
            "F1",
            "Removed_reason",
            "Kept_mRNA_ID",
            "Kept_TargetGene",
        ]
        removed_fh.write("\t".join(header_cols) + "\n")

    for key in groups:
        recs = groups[key]
        # sort by start, then end
        recs_sorted = sorted(recs, key=lambda r: (r["start"], r["end"]))

        cluster = []
        current_end = None

        def flush_cluster(cluster_list):
            if not cluster_list:
                return
            # pick best in this overlapping locus
            winner = pick_best(cluster_list)
            kept_mrna_ids.add(winner["mrna_id"])
            # all others in this cluster are removed due to same-strand collision
            for r in cluster_list:
                if r["mrna_id"] == winner["mrna_id"]:
                    continue
                removed_overlap[r["mrna_id"]] = {
                    "reason": "same_strand_collision",
                    "winner_mrna": winner["mrna_id"],
                    "winner_target": winner["target_gene"],
                }
                if removed_fh is not None:
                    row = [
                        r["target_gene"],
                        r["mrna_id"],
                        r["seqid"],
                        str(r["start"]),
                        str(r["end"]),
                        r["strand"],
                        "{:.6f}".format(r["identity"]),
                        "{:.6f}".format(r["positive"]),
                        "{:.6f}".format(r["rank"]),
                        "{:.6f}".format(r["f1"]),
                        "same_strand_collision",
                        winner["mrna_id"],
                        winner["target_gene"],
                    ]
                    removed_fh.write("\t".join(row) + "\n")

        for r in recs_sorted:
            if not cluster:
                cluster = [r]
                current_end = r["end"]
            else:
                # overlap if start <= current_end
                if r["start"] <= current_end:
                    cluster.append(r)
                    if r["end"] > current_end:
                        current_end = r["end"]
                else:
                    # no overlap with current cluster: flush it
                    flush_cluster(cluster)
                    cluster = [r]
                    current_end = r["end"]

        # flush last cluster
        flush_cluster(cluster)

    if removed_fh is not None:
        removed_fh.close()

    # Collect gene IDs (parents) for kept mRNAs
    keep_gene_ids = set()
    for mrna_id in kept_mrna_ids:
        rec = all_recs_by_mrna.get(mrna_id)
        if rec is None:
            continue
        parent = rec["parent_gene"]
        if parent:
            for g in parent.split(","):
                keep_gene_ids.add(g)

    # Summary TSV (all candidate mRNAs, kept + removed)
    if args.summary_tsv is not None:
        with open(args.summary_tsv, "w") as sfh:
            header = [
                "TargetGene",
                "mRNA_ID",
                "Seqid",
                "Start",
                "End",
                "Strand",
                "Identity",
                "Positive",
                "Rank",
                "F1",
                "Status",
                "Reason",
                "Winner_mRNA_ID",
                "Winner_TargetGene",
            ]
            sfh.write("\t".join(header) + "\n")

            for mrna_id, r in all_recs_by_mrna.items():
                if mrna_id in kept_mrna_ids:
                    status = "kept"
                    reason = "kept_nonoverlapping_best"
                    winner_mrna = ""
                    winner_target = ""
                elif mrna_id in removed_overlap:
                    status = "removed"
                    reason = removed_overlap[mrna_id]["reason"]
                    winner_mrna = removed_overlap[mrna_id]["winner_mrna"]
                    winner_target = removed_overlap[mrna_id]["winner_target"]
                elif mrna_id in removed_by_target:
                    status = "removed"
                    reason = removed_by_target[mrna_id]["reason"]
                    winner_mrna = removed_by_target[mrna_id]["winner_mrna"]
                    winner_target = removed_by_target[mrna_id]["winner_target"]
                else:
                    # Shouldn't happen, but just in case
                    status = "removed"
                    reason = "unknown"
                    winner_mrna = ""
                    winner_target = ""

                row = [
                    r["target_gene"],
                    r["mrna_id"],
                    r["seqid"],
                    str(r["start"]),
                    str(r["end"]),
                    r["strand"],
                    "{:.6f}".format(r["identity"]),
                    "{:.6f}".format(r["positive"]),
                    "{:.6f}".format(r["rank"]),
                    "{:.6f}".format(r["f1"]),
                    status,
                    reason,
                    winner_mrna,
                    winner_target,
                ]
                sfh.write("\t".join(row) + "\n")

    # Second pass: write filtered GFF3
    with open(args.gff_in) as fh, open(args.gff_out, "w") as out:
        for line in fh:
            if line.startswith("#") or not line.strip():
                out.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            ftype = cols[2]
            attrs = parse_attrs(cols[8])
            fid = attrs.get("ID")
            parent = attrs.get("Parent", "")

            if ftype == "gene":
                if fid in keep_gene_ids:
                    out.write(line)
                continue

            if ftype == "mRNA":
                if fid in kept_mrna_ids:
                    out.write(line)
                continue

            # child features (CDS, exon, UTR, etc.): keep if any Parent is a kept mRNA
            if parent:
                parent_ids = parent.split(",")
                keep_this = any(p in kept_mrna_ids for p in parent_ids)
                if keep_this:
                    out.write(line)
                continue

            # everything else is dropped


if __name__ == "__main__":
    main()
