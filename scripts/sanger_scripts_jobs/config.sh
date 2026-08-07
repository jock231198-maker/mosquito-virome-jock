#!/bin/bash
# ---------------------------------------------------------------------------
# config.sh — rutas y entornos centrales del pipeline (farm de Sanger)
#
# Todos los scripts hacen:  source "$(dirname "$0")/config.sh"
# Así cambias una ruta en un solo sitio.
#
# Uso manual desde un nodo de login:
#   source ~/jkvirome/jm79/mosquito-virome-jock/sanger_scripts_jobs/config.sh
#   echo "$SCRATCH"
# ---------------------------------------------------------------------------

# --- Raíces -----------------------------------------------------------------
USER_ROOT="/lustre/scratch126/tol/teams/lawniczak/users/jm79"

# Datos crudos iniciales (SOLO LECTURA — nunca escribir aquí)
DATA_DIR="${DATA_DIR:-$USER_ROOT/raw_data/transcriptome/MERIDA}"

# Área de trabajo: intermedios, temporales, todo lo borrable
SCRATCH="${SCRATCH:-$USER_ROOT/work/MERIDA}"

# Resultados que se conservan (uno por paso)
RESULTS_DIR="${RESULTS_DIR:-$USER_ROOT/results/MERIDA}"

# Logs de LSF
LOGS_DIR="${LOGS_DIR:-$SCRATCH/logs}"

# Etiqueta del lote de datos (se usa como subcarpeta en varias salidas)
RESULT_FROM="${RESULT_FROM:-MERIDA}"

# --- Referencias y bases de datos -------------------------------------------
REFS_DIR="${REFS_DIR:-$USER_ROOT/refs}"
ADAPTERS="${ADAPTERS:-$REFS_DIR/adapters/TruSeq3-PE-2.fa}"
BT2_INDEX="${BT2_INDEX:-$REFS_DIR/bt2_index/aedes_aegypti_index}"
STAR_INDEX="${STAR_INDEX:-$REFS_DIR/star_index}"
CHECKV_DB="${CHECKV_DB:-$REFS_DIR/checkv-db-v1.5}"
GENOMAD_DB="${GENOMAD_DB:-$REFS_DIR/genomad_db}"
DIAMOND_DB="${DIAMOND_DB:-$REFS_DIR/diamond_db/viral_refseq_proteins.dmnd}"

# --- Entornos conda (nombres tal cual aparecen en `conda env list`) ----------
# Cámbialos aquí conforme vayas creando cada entorno en el farm.
ENV_FASTQC="${ENV_FASTQC:-fastqc_0.12.1}"
ENV_MULTIQC="${ENV_MULTIQC:-multiqc_3.15}"
ENV_TRIMMOMATIC="${ENV_TRIMMOMATIC:-trimmomatic}"
ENV_BOWTIE2="${ENV_BOWTIE2:-bowtie2}"      # debe incluir bowtie2 Y samtools
ENV_SAMTOOLS="${ENV_SAMTOOLS:-$ENV_BOWTIE2}"
ENV_STAR="${ENV_STAR:-star_env}"
ENV_SPADES="${ENV_SPADES:-spades}"
ENV_QUAST="${ENV_QUAST:-quast}"
ENV_COBRA="${ENV_COBRA:-cobra}"
ENV_CHECKV="${ENV_CHECKV:-checkv}"
ENV_GENOMAD="${ENV_GENOMAD:-genomad}"
ENV_DIAMOND="${ENV_DIAMOND:-diamond}"

# --- Helper: cargar conda y activar un entorno ------------------------------
# El módulo `conda` del farm YA define la función conda (setenv CONDA_EXE +
# set-function conda), así que NO hace falta `eval "$(conda shell.bash hook)"`.
activate_env() {
    local env_name="${1:?activate_env necesita el nombre del entorno}"
    # En jobs no interactivos el sistema de módulos puede no estar inicializado
    if ! command -v module >/dev/null 2>&1; then
        [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
    fi
    module load conda
    conda activate "$env_name"
    echo "[env] $env_name activado ($(command -v conda))"
}

# --- Helper: crear el árbol de directorios ----------------------------------
make_dirs() {
    mkdir -p "$SCRATCH" "$RESULTS_DIR" "$LOGS_DIR"
}
