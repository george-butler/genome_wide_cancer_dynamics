#!/usr/bin/env bash


GFF=$1
SPECIES=$2
GFF_d=$( dirname ${GFF} )

echo -e "#sp\tn\tno_start\tno_stop\tboth\tinframe" > "${GFF_d}/psuedogene_report.txt" 
awk -v sp="$SPECIES" '
BEGIN{
  FS="\t"; OFS="\t";
  n=0; no_start=0; no_stop=0; both=0; inframe=0;
}
function finalize_seq(s){
  if (s=="") return;

  n++;

  start_ok = (substr(s,1,1)=="M");
  stop_ok  = (substr(s,length(s),1)=="*");

  if (!start_ok && !stop_ok) both++;
  else if (!start_ok) no_start++;
  else if (!stop_ok) no_stop++;

  # inframe stop: any * not at last char
  if (index(substr(s,1,length(s)-1),"*")>0) inframe++;
}
(/^>/){
  finalize_seq(seq);
  seq="";
  next;
}
{
  gsub(/[ \t\r\n]/,"",$0);
  seq = seq $0;
}
END{
  finalize_seq(seq);
  print sp, n, no_start, no_stop, both, inframe;
}
' "$GFF" >> "${GFF_d}/psuedogene_report.txt" 

