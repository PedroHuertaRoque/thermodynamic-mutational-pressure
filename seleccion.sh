#!/usr/bin/env bash
set -euo pipefail

IN="basefreq_strand_longQ20.tsv"
OUT="alt_dominante_strand_filteredQ20.tsv"

MIN_DP=50000
MIN_AF=0.0005
MIN_MINOR_STRAND_FRAC=0.2

echo -e "sample\tpos\tref_base\talt_base\tDP\talt_count\tAF_alt\talt_fwd\talt_rev\tminor_frac" > "$OUT"

awk -F'\t' -v OFS="\t" \
  -v MIN_DP="$MIN_DP" -v MIN_AF="$MIN_AF" -v MIN_MINOR_STRAND_FRAC="$MIN_MINOR_STRAND_FRAC" '
NR==1{
  # opcional: podrías chequear nombres de columnas aquí
  next
}
{
  sample=$1; pos=$2; ref=toupper($3)

  Af=$4+0; Cf=$5+0; Gf=$6+0; Tf=$7+0
  Ar=$8+0; Cr=$9+0; Gr=$10+0; Tr=$11+0
  DP=$12+0

  if (DP < MIN_DP) next

  # Totales por base
  A=Af+Ar; C=Cf+Cr; G=Gf+Gr; T=Tf+Tr

  # Elegir ALT por conteo total (equivalente a mayor AF porque DP fijo)
  alt=""; altc=-1; altf=0; altr=0

  if (ref!="A" && A>altc) {alt="A"; altc=A; altf=Af; altr=Ar}
  if (ref!="C" && C>altc) {alt="C"; altc=C; altf=Cf; altr=Cr}
  if (ref!="G" && G>altc) {alt="G"; altc=G; altf=Gf; altr=Gr}
  if (ref!="T" && T>altc) {alt="T"; altc=T; altf=Tf; altr=Tr}

  if (alt=="" || altc<=0) next

  AF = altc / DP
  if (AF < MIN_AF) next

  minor = (altf < altr ? altf : altr)
  minor_frac = minor / altc
  if (minor_frac < MIN_MINOR_STRAND_FRAC) next

  printf "%s\t%s\t%s\t%s\t%d\t%d\t%.10f\t%d\t%d\t%.6f\n", \
    sample,pos,ref,alt,DP,altc,AF,altf,altr,minor_frac
}
' "$IN" >> "$OUT"

echo "✔ Listo: $OUT"


