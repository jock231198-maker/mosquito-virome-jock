#!/bin/bash
# ---------------------------------------------------------------------------
# compress.sh — concatena los 2 lanes de cada muestra y comprime a .fastq.gz
#
#   ENTRADA   /lustre/.../transcriptome/MERIDA/
#             M1_S1_L001_R1_001.fastq   +   M1_S1_L002_R1_001.fastq
#
#   SALIDA    /lustre/.../transcriptome/MERIDA_cat/     (carpeta hermana)
#             M1_S1_R1_001.fastq.gz
#
# NO BORRA NADA. Al terminar te dice exactamente que borrar a mano.
#
# Uso:
#   ./compress.sh              # usa las rutas de config.sh
#   ./compress.sh --dry        # ensena que haria, sin escribir
#   ./compress.sh <entrada> <salida>
#
# Para muchas muestras, mandalo a la cola (es I/O + gzip):
#   bsub -o "$LOGS_DIR/compress.%J.log" -e "$LOGS_DIR/compress.%J.err" \
#        -q normal -n 8 -M 4000 \
#        -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
#        "$PWD/compress.sh"
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

DRY=0; ARGS=()
for a in "$@"; do
  case "$a" in
    --dry) DRY=1 ;;
    -*) echo "Opcion desconocida: $a  (solo --dry)"; exit 1 ;;
    *) ARGS+=("$a") ;;
  esac
done

indir="${ARGS[0]:-$DATA_DIR}"
outdir="${ARGS[1]:-${DATA_DIR}_cat}"      # carpeta hermana de MERIDA

# Cada muestra debe tener 2 lanes (L001 y L002). Una con menos = lane perdido.
EXPECTED_LANES="${EXPECTED_LANES:-2}"
THREADS="${LSB_DJOB_NUMPROC:-4}"

if command -v pigz >/dev/null 2>&1; then
  COMPRESS=(pigz -p "$THREADS" -c); CNAME="pigz -p $THREADS"
else
  COMPRESS=(gzip -c);               CNAME="gzip"
fi

echo "==========================================================="
echo "compress.sh   $(date)"
echo "Entrada  : $indir"
echo "Salida   : $outdir"
echo "Compresor: $CNAME    Lanes esperados: $EXPECTED_LANES"
[[ $DRY -eq 1 ]] && echo ">>> --dry: no se escribe nada <<<"
echo "==========================================================="
echo

[[ -d "$indir" ]] || { echo "ERROR: no existe $indir"; exit 1; }

# Detecta gzip por los magic bytes, no por la extension
is_gz() { [[ "$(head -c2 "$1" 2>/dev/null | od -An -tx1 | tr -d ' \n')" == "1f8b" ]]; }
emit()  { if is_gz "$1"; then zcat "$1"; else cat "$1"; fi; }
count_reads() {
  local n
  if is_gz "$1"; then n=$(zcat "$1" | wc -l); else n=$(wc -l < "$1"); fi
  echo $(( n / 4 ))
}

mapfile -t SAMPLES < <(
  find "$indir" -maxdepth 1 -type f \
       \( -name "*_L0*_R1_*.fastq" -o -name "*_L0*_R1_*.fastq.gz" \) -printf '%f\n' \
  | cut -d_ -f1-2 | sort -u
)
[[ ${#SAMPLES[@]} -gt 0 ]] || {
  echo "ERROR: no hay ficheros *_L0??_R1_*.fastq[.gz] en $indir"; exit 1; }

echo "Muestras detectadas: ${#SAMPLES[@]}"
echo
[[ $DRY -eq 0 ]] && mkdir -p "$outdir"

OK=0; FAIL=0

for sample in "${SAMPLES[@]}"; do
  echo "== $sample"
  for R in R1 R2; do
    mapfile -t lanes < <(
      find "$indir" -maxdepth 1 -type f \
           \( -name "${sample}_L0*_${R}_*.fastq" -o -name "${sample}_L0*_${R}_*.fastq.gz" \) | sort
    )
    out="$outdir/${sample}_${R}_001.fastq.gz"

    if [[ ${#lanes[@]} -eq 0 ]]; then
      echo "  [ERROR] $R: sin lanes"; FAIL=$((FAIL+1)); continue
    fi

    echo "  $R: ${#lanes[@]} lanes -> $(basename "$out")"

    # Un lane de menos son datos que faltan, no un detalle: la muestra queda a
    # media profundidad y nada aguas abajo te lo va a decir.
    if [[ "$EXPECTED_LANES" -gt 0 && ${#lanes[@]} -ne "$EXPECTED_LANES" ]]; then
      echo "  [ERROR] se esperaban $EXPECTED_LANES lanes y hay ${#lanes[@]}:"
      for l in "${lanes[@]}"; do echo "            $(basename "$l")"; done
      echo "          revisa la copia antes de seguir  (EXPECTED_LANES=0 para ignorar)"
      FAIL=$((FAIL+1)); continue
    fi

    [[ $DRY -eq 1 ]] && { for l in "${lanes[@]}"; do echo "       $(basename "$l")"; done; continue; }

    if [[ -s "$out" ]]; then
      echo "  [skip ] ya existe (borralo si quieres rehacerlo)"; continue
    fi

    in_reads=0
    for l in "${lanes[@]}"; do in_reads=$(( in_reads + $(count_reads "$l") )); done

    # Escribe a .tmp y renombra: un job cortado a medias nunca deja un .gz
    # truncado con pinta de bueno.
    tmp="${out}.tmp"
    if ! { for l in "${lanes[@]}"; do emit "$l"; done; } | "${COMPRESS[@]}" > "$tmp"; then
      echo "  [ERROR] fallo al comprimir"; rm -f "$tmp"; FAIL=$((FAIL+1)); continue
    fi

    out_reads=$(count_reads "$tmp")
    if [[ "$out_reads" -ne "$in_reads" ]]; then
      echo "  [ERROR] verificacion: entrada=$in_reads salida=$out_reads"
      rm -f "$tmp"; FAIL=$((FAIL+1)); continue
    fi
    if ! gzip -t "$tmp" 2>/dev/null; then
      echo "  [ERROR] el .gz no pasa gzip -t"; rm -f "$tmp"; FAIL=$((FAIL+1)); continue
    fi

    mv "$tmp" "$out"
    echo "  [ok   ] $in_reads reads verificadas"
    OK=$((OK+1))
  done
done

echo
echo "==========================================================="
echo "Creados OK : $OK"
echo "Fallos     : $FAIL"
echo "Salida     : $outdir"

if [[ $DRY -eq 1 ]]; then
  echo
  echo "--dry: no se escribio nada."
  exit 0
fi

if [[ $FAIL -gt 0 ]]; then
  echo
  echo "HAY FALLOS. Revisalos antes de borrar nada de $indir."
  echo "==========================================================="
  exit 1
fi

echo
echo "Todo verificado. Los originales siguen intactos en:"
echo "    $indir"
echo
echo "Para liberar espacio, cuando lo tengas claro, a mano:"
echo "    rm $indir/*_L0*.fastq*"
echo
echo "Ocupan: $(du -sh "$indir" 2>/dev/null | cut -f1)   |   nuevos: $(du -sh "$outdir" 2>/dev/null | cut -f1)"
echo
echo "Siguiente paso:  DATA_DIR=$outdir ./submit.sh fastqc"
echo "(o cambia DATA_DIR en config.sh para que apunte ahi de forma permanente)"
echo "==========================================================="
