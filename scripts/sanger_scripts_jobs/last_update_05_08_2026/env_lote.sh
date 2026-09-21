#!/bin/bash
# ---------------------------------------------------------------------------
# env_lote.sh — entorno por LOTE, un arbol por origen (como MERIDA)
#
#   source env_lote.sh MEXICO_INSP
#   source env_lote.sh UGANDA_UVRI
#
# Cada lote tiene su work/ y su results/ propios, exactamente igual que MERIDA:
#
#   $USER_ROOT/work/MERIDA/        22 muestras, campo, Yucatan     (ya existe)
#   $USER_ROOT/work/MEXICO_INSP/    3 muestras, colonias INSP
#   $USER_ROOT/work/UGANDA_UVRI/    3 muestras, Uganda
#
# POR QUE EL NOMBRE LLEVA PAIS *Y* ORIGEN
#   Un `work/MEXICO` a secas chocaria con MERIDA, que tambien es Mexico, y el
#   dia que llegue otro lote mexicano de campo el nombre ya estaria gastado.
#   Pais primero para que ordene por pais, origen despues para que no colisione.
#   Si prefieres otra cosa, es UNA linea: el case de abajo.
#
# ESTE FICHERO ASIGNA SIN ${VAR:-...}
#   Por eso puedes cambiar de lote en la MISMA terminal sin arrastrar el valor
#   viejo. Es justo la trampa de la §13, desactivada aqui a proposito.
# ---------------------------------------------------------------------------

LOTE="${1:-}"

USER_ROOT="/lustre/scratch126/tol/teams/lawniczak/users/jm79"
VECTOR_EVE="/lustre/scratch127/tol/projects/vector_eve/data/MosquitoTranscriptome"

case "$LOTE" in
  MEXICO_INSP)
    LOTE_DATA="$VECTOR_EVE/LANCIS_Aaa_Aas_piloto/INSP_Aaa_Labcolonies"
    LOTE_KIND="crudo"          # entra por polyG
    LOTE_GLOB="*_R1_001.fastq.gz"
    LOTE_PREFIX=""             # IgualaC2_S2 ya es un nombre sano
    LOTE_TAXON="Aedes aegypti aegypti (colonias de laboratorio, INSP)"
    ;;
  UGANDA_UVRI)
    LOTE_DATA="$VECTOR_EVE/UVRI_Aaf"
    LOTE_KIND="recortado"      # entra por bowtie2, ya paso Trim Galore
    LOTE_GLOB="*_val_1.fq.gz"
    LOTE_PREFIX="UVRI"         # 2_S3 empieza por digito: se arregla ahora
    LOTE_TAXON="Aedes aegypti formosus (Uganda)"
    ;;
  MERIDA)
    LOTE_DATA="$USER_ROOT/raw_data/transcriptome/MERIDA"
    LOTE_KIND="crudo"
    LOTE_GLOB="*_R1_001.fastq.gz"
    LOTE_PREFIX=""
    LOTE_TAXON="Aedes aegypti (campo, Merida, Yucatan)"
    ;;
  *)
    echo "Uso: source env_lote.sh <LOTE>" >&2
    echo "Lotes: MEXICO_INSP  UGANDA_UVRI  MERIDA" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac

# --- Asignacion incondicional: cambiar de lote en la misma terminal es seguro
export SCRATCH="$USER_ROOT/work/$LOTE"
export RESULTS_DIR="$USER_ROOT/results/$LOTE"
export LOGS_DIR="$SCRATCH/logs"
export RESULT_FROM="$LOTE"
export DATA_DIR="$LOTE_DATA"
export LOTE LOTE_DATA LOTE_KIND LOTE_GLOB LOTE_PREFIX LOTE_TAXON

export S="${S:-$HOME/jkvirome/jm79/mosquito-virome-jock/scripts/sanger_scripts_jobs/last_update_05_08_2026}"
if [[ ! -f "$S/config.sh" ]]; then
    echo "ERROR: no encuentro $S/config.sh — exporta S=<ruta> antes" >&2
    return 1 2>/dev/null || exit 1
fi

source "$S/config.sh"
make_dirs

cat <<EOF
[lote] $LOTE — $LOTE_TAXON
  datos ($LOTE_KIND) = $LOTE_DATA
  SCRATCH            = $SCRATCH
  RESULTS_DIR        = $RESULTS_DIR
  RESULT_FROM        = $RESULT_FROM
  BT2_INDEX          = $BT2_INDEX
EOF

[[ -n "${ENV_FASTP:-}" ]] || echo "  !! ENV_FASTP vacio en config.sh — polyg_worker.sh activaria un entorno sin nombre"
[[ -d "$LOTE_DATA" ]]     || echo "  !! no existe $LOTE_DATA (¿scratch127 montado en este nodo?)"
