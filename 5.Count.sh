#!/usr/bin/env bash
set -euo pipefail

# Mismo criterio de ordenamiento para sort y join
export LC_ALL=C

REF="ref.fasta"
POS="comunes.txt"
OUT="basefreq_strand_longQ20.tsv"

MIN_MAPQ=20
MIN_BASEQ=20
MAX_DEPTH=100000
THREADS=8

[ -f "$REF" ] || {
  echo "ERROR: No encuentro referencia: $REF"
  exit 1
}

[ -s "$POS" ] || {
  echo "ERROR: $POS no existe o está vacío"
  exit 1
}

[ -f "$REF.fai" ] || samtools faidx "$REF"

shopt -s nullglob
BAMS=( *.bam )

[[ ${#BAMS[@]} -gt 0 ]] || {
  echo "ERROR: No encontré *.bam"
  exit 1
}

echo -e "sample\tpos\tref_base\tA_fwd\tC_fwd\tG_fwd\tT_fwd\tA_rev\tC_rev\tG_rev\tT_rev\tDP\tfA\tfC\tfG\tfT\talt_fwd_frac" > "$OUT"

# Generar pileup, seleccionar posiciones comunes y contar bases
count_bases () {
  local bam="$1"

  samtools mpileup \
    -aa \
    -d "$MAX_DEPTH" \
    -q "$MIN_MAPQ" \
    -Q "$MIN_BASEQ" \
    -f "$REF" \
    "$bam" \
  | awk 'NR==FNR {pos[$1]; next} ($2 in pos)' "$POS" - \
  | awk 'BEGIN{OFS="\t"}
    function clean(b,   i,c,n,num,rest){
      gsub(/\^./,"",b)
      gsub(/\$/,"",b)

      while (match(b, /[+-][0-9]+/)) {
        rest = substr(b, RSTART+1)
        num=""

        for (i=1;i<=length(rest);i++){
          c=substr(rest,i,1)
          if(c~/[0-9]/) num=num c
          else break
        }

        n = num + 0
        b = substr(b,1,RSTART-1) substr(b, RSTART+1+length(num)+n)
      }

      gsub(/\*/,"",b)
      return b
    }
    {
      pos=$2
      ref=toupper($3)
      bases=clean($5)

      A=C=G=T=0

      for(i=1;i<=length(bases);i++){
        c=substr(bases,i,1)

        if(c=="." || c==",") c=ref
        c=toupper(c)

        if(c=="A") A++
        else if(c=="C") C++
        else if(c=="G") G++
        else if(c=="T") T++
      }

      print pos, ref, A, C, G, T
    }'
}

for BAM in "${BAMS[@]}"; do
  sample="${BAM%.bam}"

  [ -f "${BAM}.bai" ] || samtools index -@ "$THREADS" "$BAM"

  # Forward: lecturas SIN FLAG 16
  # Orden lexicográfico por posición, compatible con join
  count_bases <(
    samtools view -@ "$THREADS" -h -F 16 "$BAM" \
    | samtools view -@ "$THREADS" -b -
  ) | sort -k1,1 > "${sample}.fwd.tmp"

  # Reverse: lecturas CON FLAG 16
  # Orden lexicográfico por posición, compatible con join
  count_bases <(
    samtools view -@ "$THREADS" -h -f 16 "$BAM" \
    | samtools view -@ "$THREADS" -b -
  ) | sort -k1,1 > "${sample}.rev.tmp"

  # Unir conteos forward y reverse por posición
  join -t $'\t' -a1 -a2 -e 0 \
    -o '0,1.2,1.3,1.4,1.5,1.6,2.3,2.4,2.5,2.6' \
    "${sample}.fwd.tmp" \
    "${sample}.rev.tmp" \
  | awk -v OFS="\t" -v S="$sample" '
      {
        pos=$1
        ref=$2

        Af=$3+0; Cf=$4+0; Gf=$5+0; Tf=$6+0
        Ar=$7+0; Cr=$8+0; Gr=$9+0; Tr=$10+0

        A=Af+Ar
        C=Cf+Cr
        G=Gf+Gr
        T=Tf+Tr
        DP=A+C+G+T

        fA=(DP>0?A/DP:0)
        fC=(DP>0?C/DP:0)
        fG=(DP>0?G/DP:0)
        fT=(DP>0?T/DP:0)

        fwd=Af+Cf+Gf+Tf
        alt_fwd_frac=(DP>0?fwd/DP:0)

        print S,pos,ref,Af,Cf,Gf,Tf,Ar,Cr,Gr,Tr,DP,fA,fC,fG,fT,alt_fwd_frac
      }
    ' >> "$OUT"

  rm -f "${sample}.fwd.tmp" "${sample}.rev.tmp"
  echo "Procesado: $BAM"
done

echo "✔ Listo: $OUT"
