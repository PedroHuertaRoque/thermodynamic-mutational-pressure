#!/bin/bash

# Crear carpeta "depth" si no existe
mkdir -p depth

# Recorre todos los BAM filtrados en la misma carpeta
for bam in *.bam; do

    # Verifica que existan archivos BAM
    [ -e "$bam" ] || { echo "No se encontraron BAM filtrados."; exit 1; }

    # Nombre del archivo de salida dentro de depth/
    base=$(basename "${bam%.bam}")
    depth="depth/${base}.depth"

    echo "Generando depth para $bam → $depth"

    # Calcular cobertura por posición (-a incluye posiciones sin cobertura)
    samtools depth -a "$bam" > "$depth"
done

echo "✓ Todos los depth fueron guardados en la carpeta 'depth/'"
