#!/bin/bash
# ---------------------------------------------------------------------------
# setup_farm_dirs.sh
# Crea el esqueleto de directorios recomendado para la farm.
#
# Uso:
#   ./setup_farm_dirs.sh <scratch_root>
#
# Ejemplo (pon TU filesystem y TU grupo):
#   ./setup_farm_dirs.sh /lustre/scratch127/teamYY/jk
#
# Crea:
#   - estructura de software/scripts en tu $HOME (con backup)
#   - el árbol de datos/cómputo en tu scratch de Lustre (sin backup)
# ---------------------------------------------------------------------------
set -euo pipefail

scratch_root="${1:?Falta scratch_root, p.ej. /lustre/scratch127/teamYY/jk}"

# --- $HOME: scripts y software (pequeño, permanente, con backup) ----------
home_base="${HOME}/jkvirome"
mkdir -p "${home_base}/scripts"
echo "HOME listo: ${home_base}/scripts"

# --- Lustre scratch: todo el cómputo y datos intermedios ------------------
proj="${scratch_root}/jkviromedata"
mkdir -p "${proj}/raw_fastq"        # copia de trabajo de los crudos
mkdir -p "${proj}/cat_fastq"        # salida de cat_samples.sh
mkdir -p "${proj}/fastqc_results"   # salida de FastQC
mkdir -p "${proj}/multiqc_report"   # salida de MultiQC
mkdir -p "${proj}/logs"             # .out / .err de LSF
echo "Scratch listo: ${proj}"

cat <<EOF

Estructura creada.

  Scripts (en \$HOME, con backup):
    ${home_base}/scripts/

  Datos y cómputo (en Lustre, SIN backup):
    ${proj}/

Sugerencia: copia tus .sh a ${home_base}/scripts/ y trabaja siempre
desde ${proj}. Recuerda revisar tu cuota con:
    lfs quota -h ${scratch_root}
EOF
