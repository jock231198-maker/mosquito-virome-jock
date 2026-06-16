#!/bin/bash
# spades_assembly.sh
SECONDS=0
eval "$(conda shell.bash hook)"
conda activate spades

# run like this: ./spades_assembly.sh "/directory/to/my/fastqs/" "experiment_name" "rnavirus"
# SPAdes modes available: rnavirus, rna, meta, metaviral, plasmid, bio, corona, isolate

input_dir="$1"
result_from="$2"
spades_mode="${3:-rnavirus}"   # rnavirus por defecto si no se especifica

outdir="/Users/JK/Desktop/spades/Assembly/${result_from}/"
mkdir -p "$outdir"

# Mapa de modos válidos de SPAdes
declare -A SPADES_FLAGS=(
  [rnavirus]="--rnaviral"
  [rna]="--rna"
  [meta]="--meta"
  [metaviral]="--metaviral"
  [plasmid]="--plasmid"
  [bio]="--bio"
  [corona]="--corona"
  [isolate]="--isolate"
)

# Validar que el modo existe
if [[ -z "${SPADES_FLAGS[$spades_mode]}" ]]; then
  echo "ERROR: Modo '$spades_mode' no reconocido."
  echo "Modos disponibles: ${!SPADES_FLAGS[@]}"
  exit 1
fi

MODE_FLAG="${SPADES_FLAGS[$spades_mode]}"
echo "=== SPAdes mode: $spades_mode ($MODE_FLAG) ==="

find "$input_dir" -type f -name "*_unmapped_R1.fastq.gz" | while read -r r1; do

  sample=$(basename "$r1" _unmapped_R1.fastq.gz)
  r2="${input_dir}/${sample}_unmapped_R2.fastq.gz"

  if [[ ! -f "$r2" ]]; then
    echo "Missing R2 for $sample — skipping."
    continue
  fi

  sample_outdir="${outdir}/${sample}/"
  mkdir -p "$sample_outdir"

  echo ""
  echo "=== Assembling $sample ==="

  spades.py \
    $MODE_FLAG \
    -1 "$r1" \
    -2 "$r2" \
    -o "$sample_outdir" \
    --threads 4 \
    --memory 16

  echo "  Done: $sample"
  echo "  Output: $sample_outdir"

done

duration=$SECONDS
echo ""
echo "Process finished in $(($duration / 3600)) hour(s), $(($duration % 3600 / 60)) minute(s) and $(($duration % 60)) second(s)."
conda deactivate


# Filtrar contigs >=1000 bp del ensamble original
seqkit seq -m 1000 \
    /Users/JK/Desktop/spades/Assembly/trimm_cat_fastq_unmapped/M1_S1/contigs.fasta \ # ruta al contigs.fasta original checar directorio de salida del spades_assembly.sh
    > /Users/JK/Desktop/spades/contigs_1000bp.fasta

grep -c ">" /Users/JK/Desktop/spades/contigs_1000bp.fasta