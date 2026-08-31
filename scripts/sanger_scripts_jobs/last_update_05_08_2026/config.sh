#!/bin/bash
# ---------------------------------------------------------------------------
# config.sh — rutas y entornos centrales del pipeline (farm de Sanger)
#
# Todos los scripts hacen:  source "$(dirname "$0")/config.sh"
# Así cambias una ruta en un solo sitio.
#
# Uso manual desde un nodo de login:
#   source ~/jkvirome/jm79/mosquito-virome-jock/scripts/sanger_scripts_jobs/last_update_05_08_2026config.sh
#   echo "$SCRATCH"
# ---------------------------------------------------------------------------

# --- LSF ---------------------------------------------------------------------
export LSB_DEFAULT_USERGROUP="${LSB_DEFAULT_USERGROUP:-team222}"

# --- Raíces -----------------------------------------------------------------
SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
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
# --- Entrada del pipeline (lanes ya concatenados) ---
FASTQ_DIR="${FASTQ_DIR:-$SCRATCH/cat_fastq}"

# --- Referencias y bases de datos -------------------------------------------
REFS_DIR="${REFS_DIR:-$USER_ROOT/refs}"
GENOME_FA="${GENOME_FA:-$REFS_DIR/genome/GCF_002204515.2_AaegL5.0_genomic.fa}"
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
ENV_TRIMMOMATIC="${ENV_TRIMMOMATIC:-trimmomatic_0.41}"
ENV_BOWTIE2="${ENV_BOWTIE2:-bowtie2_2.5.5}"      # debe incluir bowtie2 Y samtools
ENV_SAMTOOLS="${ENV_SAMTOOLS:-$ENV_BOWTIE2}"
ENV_STAR="${ENV_STAR:-star_2.7.11b}"
ENV_SPADES="${ENV_SPADES:-spades_4.3.0}"
ENV_QUAST="${ENV_QUAST:-quast_5.3.0}"
ENV_CDHIT=cdhit_4.8.1
ENV_COBRA="${ENV_COBRA:-cobra}"
ENV_CHECKV="${ENV_CHECKV:-checkv}"
ENV_GENOMAD="${ENV_GENOMAD:-genomad}"
ENV_DIAMOND="${ENV_DIAMOND:-diamond}"
ENV_FASTP="${ENV_FASTP:-fastp_1.3.6}"
ENV_MEGAHIT=megahit_1.2.9
ENV_TRINITY=trinity_2.15.2
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
    mkdir -p "$SCRATCH"/{cat_fastq,nopolyg,trimmed/"$RESULT_FROM",mapped/"$RESULT_FROM",unmapped_fastq,spades,aligned,logs}
    mkdir -p "$RESULTS_DIR"/{fastqc_raw,multiqc,polyg_reports,mapping_stats,unmapped_counts,quast,genomad,checkv,diamond,qc_control}
}
