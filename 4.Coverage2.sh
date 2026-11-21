#!/bin/bash
set -euo pipefail

# Usar el mismo criterio de ordenamiento en sort y comm
export LC_ALL=C

echo "Ejecutando desde: $(pwd)"
echo

# Seleccionar archivos de entrada
shopt -s nullglob
archivos=( *.depth )

# Verificar que existan archivos
if (( ${#archivos[@]} == 0 )); then
    echo "ERROR: no se encontraron archivos barcode*_filtered.depth"
    exit 1
fi

# Crear carpeta temporal y eliminarla al terminar
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

primero=1

for depth in "${archivos[@]}"; do
    echo "Procesando $depth"

    # Extraer posiciones con cobertura > 0
    # Ordenar como texto, como requiere comm
    awk '$3 > 0 {print $2}' "$depth" \
        | sort -u > "$tmpdir/actual.txt"

    if (( primero == 1 )); then
        cp "$tmpdir/actual.txt" "$tmpdir/comunes.txt"
        primero=0
    else
        # Conservar solo posiciones presentes en ambas listas
        comm --check-order -12 \
            "$tmpdir/comunes.txt" \
            "$tmpdir/actual.txt" \
            > "$tmpdir/interseccion.txt"

        mv "$tmpdir/interseccion.txt" "$tmpdir/comunes.txt"
    fi
done

# Ordenar numéricamente solo al finalizar las comparaciones
sort -n "$tmpdir/comunes.txt" > comunes.txt

echo
echo "Archivos procesados: ${#archivos[@]}"
echo "Posiciones con cobertura > 0 en TODOS los archivos:"
wc -l < comunes.txt

if [[ ! -s comunes.txt ]]; then
    echo "AVISO: no hay posiciones con cobertura > 0 en todos los archivos."
fi

