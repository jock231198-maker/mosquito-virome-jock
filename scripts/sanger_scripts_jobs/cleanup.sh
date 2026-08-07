#!/bin/bash
# ---------------------------------------------------------------------------
# cleanup.sh — libera cuota en Lustre borrando intermedios del pipeline.
#
# POR DEFECTO NO BORRA NADA: lista qué borraría y cuánto espacio liberaría.
# Solo borra de verdad si le pasas --force.
#
# Uso:
#   ./cleanup.sh                      # dry-run de TODOS los pasos
#   ./cleanup.sh trimmed mapped       # dry-run solo de esos pasos
#   ./cleanup.sh --force trimmed      # BORRA los intermedios de trimming
#   ./cleanup.sh --list               # muestra los pasos disponibles
#
# Regla de oro: nunca toca $DATA_DIR (crudos) ni $RESULTS_DIR (resultados).
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

FORCE=0
STEPS=()

for arg in "$@"; do
  case "$arg" in
    --force) FORCE=1 ;;
    --list)
      echo "Pasos disponibles:"
      echo "  trimmed     unpaired de Trimmomatic (los paired se conservan)"
      echo "  mapped      BAM del mapeo al huesped (tras extraer unmapped)"
      echo "  sam         cualquier .sam suelto que ya tenga su .bam"
      echo "  startmp     carpetas _STARtmp y BAM sin ordenar de STAR"
      echo "  unmapped    FASTQ de no-mapeados (SOLO tras ensamblar)"
      echo "  spadestmp   carpetas tmp/ y K*/ internas de SPAdes"
      echo "  logs        logs de LSF de mas de 30 dias"
      exit 0 ;;
    -*) echo "Opcion desconocida: $arg (usa --force o --list)"; exit 1 ;;
    *) STEPS+=("$arg") ;;
  esac
done

# Sin pasos explicitos -> todos los seguros (unmapped NO entra por defecto)
if [[ ${#STEPS[@]} -eq 0 ]]; then
  STEPS=(trimmed mapped sam startmp spadestmp logs)
fi

TOTAL_BYTES=0

# purge <descripcion> <find args...>
purge() {
  local desc="$1"; shift
  local files bytes count
  files=$(find "$@" 2>/dev/null || true)
  if [[ -z "$files" ]]; then
    printf '  [ok]   %s: nada que borrar\n' "$desc"
    return
  fi
  count=$(printf '%s\n' "$files" | wc -l)
  bytes=$(printf '%s\n' "$files" | xargs -d '\n' du -sb 2>/dev/null | awk '{s+=$1} END{print s+0}')
  TOTAL_BYTES=$((TOTAL_BYTES + bytes))
  printf '  [%s] %s: %s elementos, %s\n' \
    "$([[ $FORCE -eq 1 ]] && echo BORRA || echo dry )" \
    "$desc" "$count" "$(numfmt --to=iec-i --suffix=B "$bytes" 2>/dev/null || echo "${bytes}B")"
  if [[ $FORCE -eq 1 ]]; then
    printf '%s\n' "$files" | xargs -d '\n' rm -rf --
  else
    printf '%s\n' "$files" | head -5 | sed 's/^/         /'
    if [[ $count -gt 5 ]]; then
      printf '         ... y %s mas\n' "$((count - 5))"
    fi
  fi
  return 0
}

echo "SCRATCH:  $SCRATCH"
echo "RESULTS:  $RESULTS_DIR  (intocable)"
echo "DATA:     $DATA_DIR     (intocable)"
[[ $FORCE -eq 1 ]] && echo ">>> MODO --force: SE VA A BORRAR <<<" || echo ">>> dry-run: no se borra nada <<<"
echo

for step in "${STEPS[@]}"; do
  echo "== $step"
  case "$step" in
    trimmed)
      # Los unpaired casi nunca se usan aguas abajo y pesan mucho
      purge "unpaired de Trimmomatic" \
        "$SCRATCH/trimmed" -type f -name "*_unpaired.fastq.gz"
      ;;
    mapped)
      # Solo si ya existen los FASTQ de no-mapeados
      if [[ -d "$SCRATCH/unmapped_fastq" ]] && \
         compgen -G "$SCRATCH/unmapped_fastq/*_unmapped_R1.fastq.gz" >/dev/null; then
        purge "BAM del mapeo al huesped" \
          "$SCRATCH/mapped" -type f -name "*.bam"
      else
        echo "  [skip] no hay unmapped_fastq todavia: NO se borran los BAM"
      fi
      ;;
    sam)
      # Borra .sam solo si existe su .bam correspondiente
      local_found=0
      while IFS= read -r sam; do
        [[ -z "$sam" ]] && continue
        base=$(basename "$sam" .sam)
        if find "$SCRATCH" -type f -name "${base}.bam" | grep -q .; then
          local_found=1
          size=$(du -sb "$sam" | awk '{print $1}')
          TOTAL_BYTES=$((TOTAL_BYTES + size))
          if [[ $FORCE -eq 1 ]]; then rm -f -- "$sam"; echo "  [BORRA] $sam"
          else echo "  [dry ] $sam (tiene ${base}.bam)"; fi
        fi
      done < <(find "$SCRATCH" -type f -name "*.sam" 2>/dev/null || true)
      if [[ $local_found -eq 0 ]]; then
        echo "  [ok]   ningun .sam con .bam correspondiente"
      fi
      ;;
    startmp)
      purge "carpetas _STARtmp" \
        "$SCRATCH/aligned" -type d -name "*_STARtmp"
      purge "STAR progress/SJ intermedios" \
        "$SCRATCH/aligned" -type f \( -name "*_Log.progress.out" -o -name "*_Log.out" \)
      ;;
    unmapped)
      # PELIGROSO: solo tras confirmar que los ensamblajes estan completos
      if compgen -G "$SCRATCH/spades/*/contigs.fasta" >/dev/null; then
        purge "FASTQ de no-mapeados (ya ensamblados)" \
          "$SCRATCH/unmapped_fastq" -type f -name "*.fastq.gz"
      else
        echo "  [skip] no hay contigs.fasta: NO se borran los unmapped"
      fi
      ;;
    spadestmp)
      purge "tmp/ y K*/ internos de SPAdes" \
        "$SCRATCH/spades" -maxdepth 2 -type d \( -name "tmp" -o -name "K[0-9]*" -o -name "misc" \)
      purge "before_rr / dataset intermedios de SPAdes" \
        "$SCRATCH/spades" -type f \( -name "before_rr.fasta" -o -name "*.yaml" \)
      ;;
    logs)
      purge "logs de LSF > 30 dias" \
        "$LOGS_DIR" -type f \( -name "*.log" -o -name "*.err" \) -mtime +30
      ;;
    *)
      echo "  paso desconocido: $step (usa --list)"
      ;;
  esac
  echo
done

echo "-----------------------------------------------------------"
printf 'Total %s: %s\n' \
  "$([[ $FORCE -eq 1 ]] && echo liberado || echo liberable)" \
  "$(numfmt --to=iec-i --suffix=B "$TOTAL_BYTES" 2>/dev/null || echo "${TOTAL_BYTES}B")"
if [[ $FORCE -eq 0 ]]; then
  echo "Repite con --force para borrar de verdad."
fi
echo
echo "Cuota actual:"
lfs quota -h "$SCRATCH" 2>/dev/null || echo "  (lfs quota no disponible aqui)"
