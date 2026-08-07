#!/bin/bash
# ---------------------------------------------------------------------------
# run_multiqc.sh — OBSOLETO: ahora es un alias de multi_qc.sh
#
# Antes había dos scripts de MultiQC (uno para la Mac, otro para el farm).
# Se han unificado en multi_qc.sh, que lee las rutas de config.sh.
# Este fichero se mantiene solo para no romper comandos antiguos.
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "[aviso] run_multiqc.sh esta obsoleto; usa multi_qc.sh" >&2
exec "$SCRIPT_DIR/multi_qc.sh" "$@"
