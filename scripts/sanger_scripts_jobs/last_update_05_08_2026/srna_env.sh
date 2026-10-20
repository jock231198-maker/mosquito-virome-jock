#!/bin/bash
# ---------------------------------------------------------------------------
# srna_env.sh — entorno del proyecto de SMALL RNA (se sourcea, no se ejecuta)
#
#   source "$S/srna_env.sh"              # lote por defecto: MERIDA_SRNA
#   source "$S/srna_env.sh" MERIDA_SRNA
#
# Mismo patron que env_lote.sh: un arbol work/ y results/ propio, separado del
# RNA largo. Nada se escribe en work/MERIDA ni en results/MERIDA.
#
#   $USER_ROOT/work/MERIDA_SRNA/       intermedios (borrables)
#   $USER_ROOT/results/MERIDA_SRNA/    lo que se conserva
#
# ASIGNACION SIN ${VAR:-...} en lo que define el arbol (la trampa de la §13):
# cambiar de lote en la misma terminal no arrastra valores viejos.
# Lo que SI se puede sobreescribir desde fuera: SRNA_DATA, SRNA_ADAPTER,
# ENV_SRNA (se exportan antes de sourcear).
# ---------------------------------------------------------------------------

SRNA_LOTE="${1:-MERIDA_SRNA}"

USER_ROOT="/lustre/scratch126/tol/teams/lawniczak/users/jm79"

export S="${S:-$HOME/jkvirome/jm79/mosquito-virome-jock/scripts/sanger_scripts_jobs/last_update_05_08_2026}"
if [[ ! -f "$S/config.sh" ]]; then
    echo "ERROR: no encuentro $S/config.sh — exporta S=<ruta> antes" >&2
    return 1 2>/dev/null || exit 1
fi
source "$S/config.sh"           # REFS_DIR, BT2_INDEX, activate_env, make_dirs...

case "$SRNA_LOTE" in
  MERIDA_SRNA)
    # Crudos de small RNA (solo lectura). Mexico1/ y Mexico2/ contienen las librerias de MERIDA.
    SRNA_DATA="${SRNA_DATA:-/lustre/scratch126/tol/teams/lawniczak/projects/vector_eve/data/MosquitoSmallRNA}"
    SRNA_GLOB="${SRNA_GLOB:-*.fq.gz}"
    SRNA_LONG_LOTE="MERIDA"     # lote de RNA largo con el que se empareja
    ;;
  *)
    echo "Uso: source srna_env.sh [MERIDA_SRNA]" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac

# --- Arbol (incondicional) --------------------------------------------------
export SCRATCH="$USER_ROOT/work/$SRNA_LOTE"
export RESULTS_DIR="$USER_ROOT/results/$SRNA_LOTE"
export LOGS_DIR="$SCRATCH/logs"
export RESULT_FROM="$SRNA_LOTE"
export SRNA_LOTE SRNA_DATA SRNA_GLOB SRNA_LONG_LOTE

# Samplesheet: la fuente de verdad de que libreria es que muestra MERIDA
export SRNA_SHEET="$RESULTS_DIR/srna_samplesheet.tsv"

# --- Libreria: Illumina TruSeq Small RNA ------------------------------------
# RA3 (3'): TGGAATTCTCGGGTGCCAAGG. Sin UMI ni random-mers: no hay que cortar 5'.
export SRNA_ADAPTER="${SRNA_ADAPTER:-TGGAATTCTCGGGTGCCAAGG}"

# Ventanas de tamano (iguales a Lau et al. 2026 para poder comparar):
#   siRNA 18-23 nt (pico esperado en 21), piRNA 24-35 nt
export SRNA_MINLEN=15 SRNA_MAXLEN=45            # lo que sobrevive al trimming
export SRNA_SI_MIN=18 SRNA_SI_MAX=23
export SRNA_PI_MIN=24 SRNA_PI_MAX=35

# --- Referencias bowtie1 (las construye srna_build_refs.sh) -----------------
export SRNA_REFS="$REFS_DIR/srna"
export BT1_GENOME="$SRNA_REFS/bt1_genome/AaegL5"
export BT1_RRNA="$SRNA_REFS/bt1_rrna/rrna"
export BT1_VIRAL="$SRNA_REFS/bt1_viral/viral_targets"

# --- Entorno conda: uno solo para todo small RNA ----------------------------
#   conda create -n srna -c conda-forge -c bioconda \
#       'cutadapt>=4.4' 'bowtie=1.3.1' samtools seqkit multiqc pysam python=3.11
export ENV_SRNA="${ENV_SRNA:-srna}"

make_dirs
mkdir -p "$SRNA_REFS"

cat <<EOT
[srna] $SRNA_LOTE  (empareja con RNA largo: $SRNA_LONG_LOTE)
  crudos      = $SRNA_DATA  ($SRNA_GLOB)
  SCRATCH     = $SCRATCH
  RESULTS_DIR = $RESULTS_DIR
  samplesheet = $SRNA_SHEET
  adaptador 3'= $SRNA_ADAPTER
  env conda   = $ENV_SRNA
EOT
[[ -d "$SRNA_DATA" ]] || echo "  !! no existe $SRNA_DATA — exporta SRNA_DATA=<ruta> y vuelve a sourcear"
