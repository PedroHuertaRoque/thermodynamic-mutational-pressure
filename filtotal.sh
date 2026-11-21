#!/usr/bin/env bash
set -euo pipefail

THREADS=8
MINLEN=250
MAXLEN=500

OUTDIR="ygjV_250"
mkdir -p "$OUTDIR"

REPORT="${OUTDIR}/kept_report.tsv"
echo -e "bam\tout_bam\ttotal_alignments\tkept_alignments\tpct_kept" > "$REPORT"

shopt -s nullglob
bams=( *.bam )
if [[ ${#bams[@]} -eq 0 ]]; then
  echo "ERROR: No hay *.bam en $(pwd)"
  exit 1
fi

for bam in "${bams[@]}"; do
  base=${bam%.bam}
  out="${OUTDIR}/${base}_fil250.bam"

  total=$(samtools view -c "$bam")

  samtools view -@ "$THREADS" -h "$bam" \
  | awk -v MIN="$MINLEN" -v MAX="$MAXLEN" '
    function ref_span(cig,   n,op,rest,span,token){
      span=0; rest=cig
      while (match(rest, /^[0-9]+[A-Z=]/)) {
        token=substr(rest, RSTART, RLENGTH)
        n=substr(token, 1, length(token)-1)+0
        op=substr(token, length(token), 1)
        # Consume referencia (lo que cubre del gen)
        if (op=="M" || op=="=" || op=="X" || op=="D" || op=="N") span += n
        rest = substr(rest, RLENGTH+1)
      }
      return span
    }
    /^@/ {print; next}
    { s=ref_span($6); if(s>=MIN && s<=MAX) print }
  ' \
  | samtools view -@ "$THREADS" -b -o "$out" -

  samtools index -@ "$THREADS" "$out" || true

  kept=$(samtools view -c "$out")
  pct=$(awk -v k="$kept" -v t="$total" 'BEGIN{ if(t>0) printf "%.2f", (100*k/t); else print "0.00"}')

  echo -e "$bam\t$out\t$total\t$kept\t$pct" >> "$REPORT"
  echo "$bam: total=$total kept=$kept (${pct}%)"
done

echo "Listo. Reporte: $REPORT"
