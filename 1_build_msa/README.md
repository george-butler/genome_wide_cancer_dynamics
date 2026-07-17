# MSA Pipeline

## Installation

Create the conda environment:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate geneMSAbuilder
```

Install any additional Python dependencies:

```bash
pip install -r requirements.txt
```

### Install MACSE

The pipeline expects the MACSE repository to be located **next to the pipeline directory**, because `scripts/submit_MSA.sh` automatically looks for:

```text
../MACSE_V2_PIPELINES
```

relative to the `scripts/` directory.

Clone MACSE into the pipeline root directory:

```bash
git clone https://github.com/ranwez/MACSE_V2_PIPELINES.git
```

Expected directory structure:

```text
geneMSAbuilder/
├── batch_submit_miniprot.sh
├── batch_submit_MSA.sh
├── environment.yml
├── requirements.txt
├── scripts/
│   └── submit_MSA.sh
└── MACSE_V2_PIPELINES/
    └── OMM_MACSE/
        └── S_OMM_MACSE_V12.02.sh
```

The script automatically resolves:

```bash
LAUNCH_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MACSE_PATH="${LAUNCH_DIR}/../MACSE_V2_PIPELINES"
```

and expects the following file to exist:

```text
MACSE_V2_PIPELINES/OMM_MACSE/S_OMM_MACSE_V12.02.sh
```

## Step 1: Create Gene Batches

```bash
scripts/make_msa_gene_chunks.sh data/test_prot.fasta data/test 1000
```

This creates:

```text
data/test/
├── human_proteins.gene_sizes.tsv
├── human_proteins.sorted_by_size.gene_list.txt
├── MSA_1/
├── MSA_2/
└── ...
```

## Step 2: Run Miniprot

```bash
./batch_submit_miniprot.sh
```

## Step 3: Run MSA Jobs

After the Miniprot jobs complete:

```bash
./batch_submit_MSA.sh data/test
```

This automatically finds all `MSA_*` directories and submits a SLURM array job for each batch.

## Workflow

```bash
conda activate geneMSAbuilder

scripts/make_msa_gene_chunks.sh data/test_prot.fasta data/test 1000

./batch_submit_miniprot.sh

./batch_submit_MSA.sh data/test
```

