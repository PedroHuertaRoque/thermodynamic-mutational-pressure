#!/bin/bash

# =========================
# CONFIGURACIÓN
# =========================
REF="ref.fasta"
THREADS=8
MIN_MAPQ=20

SUMMARY_RAW="mapping_stats_raw.tsv"
SUMMARY_FILTERED="mapping_stats_filtered.tsv"

# =========================
# VERIFICAR REFERENCIA
# =========================
[ -f "$REF" ] || { echo "ERROR: No se encuentra $REF"; exit 1; }

# =========================
# ENCABEZADOS
# =========================
echo -e "sample\ttotal_reads\tmapped\tprimary\tsecondary\tsupplementary\tunmapped\tpct_mapped" > "$SUMMARY_RAW"
echo -e "sample\tprimary_filtered" > "$SUMMARY_FILTERED"

# =========================
# PROCESAR FASTQ
# =========================
for fq in pass_barcode*.fastq; do
    [ -e "$fq" ] || { echo "No se encontraron FASTQ"; exit 1; }

    sample="${fq%.fastq}"

    bam_raw="${sample}.bam"
    bam_filt="${sample}_gene.bam"

    echo "▶ Mapeando $fq"

    # 1️⃣ MAPEO
    minimap2 -t $THREADS -ax map-ont "$REF" "$fq" \
        | samtools sort -@ $THREADS -o "$bam_raw"

    samtools index "$bam_raw"

    # =========================
    # 2️⃣ ESTADÍSTICAS ANTES DEL FILTRO
    # =========================
    total=$(samtools view -c "$bam_raw")
    unmapped=$(samtools view -c -f 0x4 "$bam_raw")
    mapped=$(samtools view -c -F 0x4 "$bam_raw")
    primary=$(samtools view -c -F 0x904 "$bam_raw")
    secondary=$(samtools view -c -f 0x100 "$bam_raw")
    supplementary=$(samtools view -c -f 0x800 "$bam_raw")

    pct=$(awk -v m=$mapped -v t=$total \
        'BEGIN{ if(t>0) printf "%.2f", (m/t)*100; else print "0.00" }')

    echo -e "$sample\t$total\t$mapped\t$primary\t$secondary\t$supplementary\t$unmapped\t$pct" \
        >> "$SUMMARY_RAW"

    # =========================
    # 3️⃣ FILTRO (SOLO PRIMARIOS BUENOS)
    # =========================
    samtools view -b \
        -F 0x4 -F 0x900 -q $MIN_MAPQ \
        "$bam_raw" > "$bam_filt"

    samtools index "$bam_filt"

    primary_filt=$(samtools view -c "$bam_filt")
    echo -e "$sample\t$primary_filt" >> "$SUMMARY_FILTERED"

    # =========================
    # 4️⃣ LIMPIEZA
    # =========================
    rm "$bam_raw" "$bam_raw.bai"

    echo "✓ Hecho $sample"
    echo
done

echo "✔ Estadísticas crudas:    $SUMMARY_RAW"
echo "✔ Estadísticas filtradas: $SUMMARY_FILTERED"


