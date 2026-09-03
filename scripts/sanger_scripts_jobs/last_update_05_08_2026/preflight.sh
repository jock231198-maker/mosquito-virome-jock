#!/bin/bash
# ---------------------------------------------------------------------------
# preflight.sh — comprueba TODO antes de un bsub e imprime el comando listo.
#
# Nace de cuatro envios perdidos en una noche, ninguno por culpa del worker:
#   - falta de corchetes [1-N]  -> LSF pone LSB_JOBINDEX=0 y sed revienta
#   - bsub lanzado desde otra carpeta -> $PWD apuntaba mal (exit 127)
#   - variable de ruta vacia -> "/worker.sh" (exit 127)
#   - lista de muestras vacia porque su fichero fuente no existia (exit 1)
#
# Los cuatro se veian antes de enviar. Esto los mira por ti.
#
# Uso:
#   ./preflight.sh                # lista los pasos disponibles
#   ./preflight.sh fastqc         # comprueba e imprime el bsub
#   ./preflight.sh spades rna     # los pasos con variante aceptan un 2o argumento
#
# NO lanza nada. Imprime el comando para que lo copies tras leerlo.
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
: "${SCRIPTS_DIR:=$SCRIPT_DIR}"

ERR=0
ok()   { printf '  \033[32m[ok   ]\033[0m %s\n' "$*"; }
bad()  { printf '  \033[31m[FALTA]\033[0m %s\n' "$*"; ERR=$((ERR+1)); }
warn() { printf '  [warn ] %s\n' "$*"; }

# --- Comprobadores ----------------------------------------------------------
chk_worker() {
  local w="$SCRIPTS_DIR/$1"
  [[ -f "$w" ]] || { bad "no existe: $w"; return; }
  [[ -x "$w" ]] || { bad "sin permiso de ejecucion (chmod +x): $w"; return; }
  if head -1 "$w" | grep -q $'\r'; then bad "tiene CRLF (sed -i 's/\\r\$//'): $w"; return; fi
  bash -n "$w" 2>/dev/null || { bad "error de sintaxis: $w"; return; }
  ok "worker: $w"
}

chk_list() {   # chk_list <fichero> <esperado|->
  local f="$1" want="${2:--}" n
  [[ -f "$f" ]] || { bad "no existe la lista: $f"; N=0; return; }
  n=$(grep -cve '^[[:space:]]*$' "$f")
  N=$n
  (( n > 0 )) || { bad "lista VACIA: $f"; return; }
  if [[ "$want" != "-" && "$n" != "$want" ]]; then
    warn "lista con $n lineas, se esperaban $want: $f"
  else
    ok "lista: $f  ($n elementos)"
  fi
}

chk_dir() {
  local d="$1" pat="${2:-}" n
  [[ -d "$d" ]] || { bad "no existe el directorio: $d"; return; }
  if [[ -n "$pat" ]]; then
    n=$(find "$d" -maxdepth 2 -name "$pat" 2>/dev/null | wc -l)
    (( n > 0 )) || { bad "sin ficheros '$pat' en $d"; return; }
    ok "entrada: $d  ($n x $pat)"
  else
    ok "directorio: $d"
  fi
}

chk_env() {
  local e="$1"
  [[ -n "$e" ]] || { bad "variable de entorno conda vacia"; return; }
  if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$e"; then
    ok "entorno conda: $e"
  else
    bad "el entorno conda '$e' no existe (revisa config.sh)"
  fi
}

chk_file() { [[ -s "$1" ]] && ok "$2: $1" || bad "$2 no existe o vacio: $1"; }

chk_writable() {
  local d="$1"
  mkdir -p "$d" 2>/dev/null && [[ -w "$d" ]] && ok "salida escribible: $d" \
    || bad "no se puede escribir en: $d"
}

chk_queue() {
  local q="$1" line
  line=$(bqueues "$q" 2>/dev/null | awk 'NR==2{printf "PEND %s / RUN %s", $9, $10}')
  [[ -n "$line" ]] && ok "cola $q  ($line)" || warn "no se pudo consultar la cola $q"
}

# ---------------------------------------------------------------------------
usage() {
  cat <<'EOF'
Uso: ./preflight.sh <paso> [variante]

Pasos:
  cat        concatenado de lanes (compress.sh)      no es array
  fastqc     FastQC                  [carpeta]       cat_fastq | nopolyg | trimmed | unpaired
  polyg      recorte de polyG (fastp)
  trim       Trimmomatic
  smrdb      base + indice de SortMeRNA  [variante]  fast | default | sensitive
  rrna       depleción de rRNA (SortMeRNA) [variante]
  btmap      mapeo al huesped (bowtie2, sin BAM)
  btmapbam   mapeo al huesped (bowtie2, con BAM)
  star       mapeo al huesped (STAR corregido)
  spades     ensamblaje            [modo]           rnavirus | rna | meta | metaviral
  megahit    ensamblaje            [preset]         "" | meta-sensitive | meta-large
  trinity    ensamblaje            [sufijo lista]   usa $SCRATCH/trin_<sufijo>.txt
  quastcmp   QUAST comparando ensambladores (1 corrida por muestra)
  quast      QC de ensamblaje (1 corrida por ensamblaje)
  genomad    clasificacion viral
  checkv     calidad de contigs virales
  diamond    blastx

Recursos MEDIDOS en este proyecto (no estimados):
  FastQC       279 MB / 274 s    -n 1  -M 1000
  fastp       1210 MB / 180 s    -n 4  -M 2000
  Trimmomatic 3190 MB / 1103 s   -n 4  -M 6000   (JAVA_MEM=4g)
  bowtie2     1611 MB / 6959 s   -n 16 -M 6000
  STAR             — / 1400 s    -n 8  -M 20000  (indice 13.1 GB en RAM)
  SPAdes      2077 MB / 1703 s   -n 4  -M 8000
  MEGAHIT      558 MB /  277 s   -n 2  -M 2000   (mediana 300 MB / ~180 s)
  Trinity      754 MB / 6596 s   -n 8  -M 4000   (M11_S3, 1.26M pares)

  MEGAHIT es 20-50x mas barato que SPAdes para un numero de contigs comparable.
  Trinity cuesta POCA MEMORIA (754 MB medidos) pero mucho TIEMPO: 1h50m para
  1.26M pares. El limite es la cola, no la RAM: pedir 12 GB fue 16x de mas.
EOF
}

step="${1:-}"; variant="${2:-}"
[[ -n "$step" ]] || { usage; exit 0; }

echo "==========================================================="
echo "preflight: $step ${variant:+($variant)}"
echo "==========================================================="
echo
ok "SCRIPTS_DIR = $SCRIPTS_DIR"
[[ -d "$SCRIPTS_DIR" ]] || bad "SCRIPTS_DIR no es un directorio"

BSUB=""
case "$step" in

  cat)
    chk_worker compress.sh
    chk_dir "$DATA_DIR" "*_R1_*.fastq.gz"
    chk_writable "$SCRATCH/cat_fastq"
    chk_queue normal
    BSUB="bsub -J compress \\
     -o \"$LOGS_DIR/compress.%J.log\" -e \"$LOGS_DIR/compress.%J.err\" \\
     -q normal -n 8 -M 4000 \\
     -R \"select[mem>4000] rusage[mem=4000] span[hosts=1]\" \\
     \"$SCRIPTS_DIR/compress.sh\""
    ;;

  fastqc)
    src="${variant:-cat_fastq}"
    case "$src" in
      cat_fastq) indir="$SCRATCH/cat_fastq";              list="$SCRATCH/filelist.txt";          out="$RESULTS_DIR/fastqc_raw";      pat="*.fastq.gz" ;;
      nopolyg)   indir="$SCRATCH/nopolyg";                list="$SCRATCH/filelist_nopolyg.txt";  out="$RESULTS_DIR/fastqc_nopolyg";  pat="*.fastq.gz" ;;
      trimmed)   indir="$SCRATCH/trimmed/$RESULT_FROM";   list="$SCRATCH/filelist_trimmed.txt";  out="$RESULTS_DIR/fastqc_trimmed";  pat="*_001_paired.fastq.gz" ;;
      unpaired)  indir="$SCRATCH/trimmed/$RESULT_FROM";   list="$SCRATCH/filelist_unpaired.txt"; out="$RESULTS_DIR/fastqc_unpaired"; pat="*_unpaired.fastq.gz" ;;
      *) echo "variante desconocida: $src  (cat_fastq|nopolyg|trimmed|unpaired)"; exit 1 ;;
    esac
    chk_worker fastqc_worker.sh
    chk_dir "$indir" "$pat"
    chk_list "$list"
    chk_env "$ENV_FASTQC"
    chk_writable "$out"
    chk_queue normal
    echo
    echo "  Si la lista no existe o esta desfasada:"
    echo "    find $indir -type f -name \"$pat\" | sort > $list"
    BSUB="bsub -J \"fq${src}[1-$N]%20\" \\
     -o \"$LOGS_DIR/fq${src}.%J.%I.log\" -e \"$LOGS_DIR/fq${src}.%J.%I.err\" \\
     -q normal -n 1 -M 1000 \\
     -R \"select[mem>1000] rusage[mem=1000] span[hosts=1]\" \\
     \"$SCRIPTS_DIR/fastqc_worker.sh $list $out\""
    ;;

  polyg)
    chk_worker polyg_worker.sh
    chk_dir "$SCRATCH/cat_fastq" "*_R1_001.fastq.gz"
    chk_list "$SCRATCH/samples.txt" 22
    chk_env "$ENV_FASTP"
    chk_writable "$SCRATCH/nopolyg"
    chk_queue normal
    BSUB="bsub -J \"polyg[1-$N]%10\" \\
     -o \"$LOGS_DIR/polyg.%J.%I.log\" -e \"$LOGS_DIR/polyg.%J.%I.err\" \\
     -q normal -n 4 -M 2000 \\
     -R \"select[mem>2000] rusage[mem=2000] span[hosts=1]\" \\
     \"$SCRIPTS_DIR/polyg_worker.sh $SCRATCH/samples.txt $SCRATCH/cat_fastq\""
    ;;

  trim)
    chk_worker trimmomatic_worker.sh
    chk_dir "$SCRATCH/nopolyg" "*_R1_001.fastq.gz"
    chk_list "$SCRATCH/samples.txt" 22
    chk_env "$ENV_TRIMMOMATIC"
    chk_file "$ADAPTERS" "adaptadores"
    grep -q 'PrefixPE/1' "$ADAPTERS" 2>/dev/null \
      && ok "adaptadores con PrefixPE/1 y /2 (modo palindromico)" \
      || warn "el FASTA de adaptadores no tiene PrefixPE/1: cogiste TruSeq3-PE.fa en vez del -2?"
    chk_writable "$SCRATCH/trimmed/$RESULT_FROM"
    chk_queue normal
    BSUB="bsub -J \"trim[1-$N]%10\" \\
     -o \"$LOGS_DIR/trim.%J.%I.log\" -e \"$LOGS_DIR/trim.%J.%I.err\" \\
     -q normal -n 4 -M 6000 \\
     -R \"select[mem>6000] rusage[mem=6000] span[hosts=1]\" \\
     \"JAVA_MEM=4g $SCRIPTS_DIR/trimmomatic_worker.sh $SCRATCH/samples.txt $SCRATCH/nopolyg\""
    ;;

  btmap)
    chk_worker bowtie_map_worker_fast.sh
    chk_dir "$SCRATCH/trimmed/$RESULT_FROM" "*_R1_001_paired.fastq.gz"
    chk_list "$SCRATCH/samples.txt" 22
    chk_env "$ENV_BOWTIE2"
    for e in 1 2 3 4 rev.1 rev.2; do
      [[ -s "${BT2_INDEX}.${e}.bt2" || -s "${BT2_INDEX}.${e}.bt2l" ]] \
        || { bad "falta el indice bowtie2: ${BT2_INDEX}.${e}.bt2"; break; }
    done
    [[ $ERR -eq 0 ]] && ok "indice bowtie2: $BT2_INDEX (6 ficheros)"
    chk_writable "$SCRATCH/unmapped_fastq"
    chk_queue normal
    BSUB="bsub -J \"btmap[1-$N]%10\" \\
     -o \"$LOGS_DIR/btmap.%J.%I.log\" -e \"$LOGS_DIR/btmap.%J.%I.err\" \\
     -q normal -n 16 -M 6000 \\
     -R \"select[mem>6000] rusage[mem=6000] span[hosts=1]\" \\
     \"$SCRIPTS_DIR/bowtie_map_worker_fast.sh $SCRATCH/samples.txt $SCRATCH/trimmed/$RESULT_FROM\""
    ;;

  btmapbam)
    chk_worker bowtie_map_worker.sh
    chk_dir "$SCRATCH/trimmed/$RESULT_FROM" "*_R1_001_paired.fastq.gz"
    chk_list "$SCRATCH/samples.txt" 22
    chk_env "$ENV_BOWTIE2"
    conda run -n "$ENV_BOWTIE2" which samtools >/dev/null 2>&1 \
      && ok "samtools presente en $ENV_BOWTIE2" \
      || warn "no encuentro samtools en $ENV_BOWTIE2 (el worker encadena bowtie2 | samtools)"
    chk_writable "$SCRATCH/mapped_bam/$RESULT_FROM"
    chk_queue normal
    BSUB="bsub -J \"btmapbam[1-$N]\" \\
     -o \"$LOGS_DIR/btmapbam.%J.%I.log\" -e \"$LOGS_DIR/btmapbam.%J.%I.err\" \\
     -q normal -n 16 -M 6000 \\
     -R \"select[mem>6000] rusage[mem=6000] span[hosts=1]\" \\
     \"STATS_DIR=$RESULTS_DIR/mapping_stats_bam $SCRIPTS_DIR/bowtie_map_worker.sh \\
      $SCRATCH/samples.txt $SCRATCH/trimmed/$RESULT_FROM $BT2_INDEX $SCRATCH/mapped_bam\""
    ;;

  star)
    chk_worker star_align_worker_fast.sh
    chk_dir "$SCRATCH/trimmed/$RESULT_FROM" "*_R1_001_paired.fastq.gz"
    chk_list "$SCRATCH/samples.txt" 22
    chk_env "$ENV_STAR"
    chk_file "$STAR_INDEX/SA" "indice STAR (SA)"
    chk_file "$STAR_INDEX/Genome" "indice STAR (Genome)"
    chk_writable "$SCRATCH/unmapped_fastq_star"
    chk_queue normal
    warn "cada elemento carga $(du -sh "$STAR_INDEX" 2>/dev/null | cut -f1) de indice en RAM: no subas de %6"
    BSUB="bsub -J \"star[1-$N]%6\" \\
     -o \"$LOGS_DIR/star.%J.%I.log\" -e \"$LOGS_DIR/star.%J.%I.err\" \\
     -q normal -n 8 -M 20000 \\
     -R \"select[mem>20000] rusage[mem=20000] span[hosts=1]\" \\
     \"$SCRIPTS_DIR/star_align_worker_fast.sh $SCRATCH/samples.txt $SCRATCH/trimmed/$RESULT_FROM\""
    ;;

  spades)
    mode="${variant:-rnavirus}"
    case "$mode" in
      rnavirus) out="$SCRATCH/spades";           stats="$RESULTS_DIR/assembly_stats" ;;
      *)        out="$SCRATCH/spades_$mode";     stats="$RESULTS_DIR/assembly_stats_$mode" ;;
    esac
    chk_worker assembly_spades.sh
    chk_dir "$SCRATCH/unmapped_fastq" "*_unmapped_R1.fastq.gz"
    chk_list "$SCRATCH/asm_samples.txt" 22
    chk_env "$ENV_SPADES"
    chk_writable "$out"
    chk_writable "$stats"
    chk_queue normal
    echo
    echo "  Si la lista no existe:"
    echo "    ls \$SCRATCH/unmapped_fastq/*_unmapped_R1.fastq.gz | xargs -n1 basename \\"
    echo "      | sed 's/_unmapped_R1.fastq.gz//' | sort > \$SCRATCH/asm_samples.txt"
    BSUB="bsub -J \"spades${mode}[1-$N]%22\" \\
     -o \"$LOGS_DIR/spades${mode}.%J.%I.log\" -e \"$LOGS_DIR/spades${mode}.%J.%I.err\" \\
     -q normal -n 4 -M 8000 \\
     -R \"select[mem>8000] rusage[mem=8000] span[hosts=1]\" \\
     \"SPADES_MEM=6 ASM_STATS_DIR=$stats $SCRIPTS_DIR/assembly_spades.sh \\
      $SCRATCH/asm_samples.txt $SCRATCH/unmapped_fastq $mode $out\""
    ;;

  smrdb)
    v="${variant:-default}"
    chk_worker build_sortmerna_db.sh
    chk_env "$ENV_SORTMERNA"
    chk_writable "${SMR_DIR:-$REFS_DIR/sortmerna}"
    chk_queue normal
    echo
    echo "  No es un array. Necesita RED para descargar database.tar.gz (~2 GB)."
    BSUB="bsub -J smrdb \\
     -o \"$LOGS_DIR/smrdb.%J.log\" -e \"$LOGS_DIR/smrdb.%J.err\" \\
     -q normal -n 4 -M 8000 \\
     -R \"select[mem>8000] rusage[mem=8000] span[hosts=1]\" \\
     \"$SCRIPTS_DIR/build_sortmerna_db.sh $v\""
    ;;

  rrna)
    v="${variant:-default}"
    db="${SMR_DB:-$REFS_DIR/sortmerna/smr_v4.3_${v}_db.fasta}"
    ix="${SMR_IDX:-$REFS_DIR/sortmerna/idx_${v}}"
    chk_worker rrna_worker.sh
    chk_dir "$SCRATCH/unmapped_fastq" "*_unmapped_R1.fastq.gz"
    chk_list "$SCRATCH/asm_samples.txt" 22
    chk_env "$ENV_SORTMERNA"
    chk_file "$db" "base de rRNA"
    chk_dir "$ix"
    chk_writable "$SCRATCH/norrna"
    chk_writable "$RESULTS_DIR/rrna_stats"
    chk_queue long
    echo
    warn "El indice se usa en SOLO LECTURA. Si no existe, corre antes:  ./preflight.sh smrdb"
    warn "SortMeRNA es lento: cola 'long'. Mide con un smoke antes de las 22."
    BSUB="bsub -J \"rrna[1-$N]%8\" \\
     -o \"$LOGS_DIR/rrna.%J.%I.log\" -e \"$LOGS_DIR/rrna.%J.%I.err\" \\
     -q long -n 8 -M 16000 \\
     -R \"select[mem>16000] rusage[mem=16000] span[hosts=1]\" \\
     \"SMR_DB=$db SMR_IDX=$ix $SCRIPTS_DIR/rrna_worker.sh \\
      $SCRATCH/asm_samples.txt $SCRATCH/unmapped_fastq\""
    ;;

  megahit)
    pre="${variant:-}"
    out="$SCRATCH/megahit${pre:+_$pre}"
    stats="$RESULTS_DIR/assembly_stats_megahit"
    chk_worker assembly_megahit.sh
    chk_dir "$SCRATCH/unmapped_fastq" "*_unmapped_R1.fastq.gz"
    chk_list "$SCRATCH/asm_samples.txt" 22
    chk_env "$ENV_MEGAHIT"
    chk_writable "$out"
    chk_writable "$stats"
    chk_queue normal
    echo
    echo "  Recuerda: -m se le pasa en BYTES absolutos dentro del worker."
    echo "  Una fraccion haria que MEGAHIT leyera la RAM del NODO, no la del job."
    echo
    echo "  Recursos MEDIDOS en las 22 (31 ago): max 558 MB, 88-277 s por muestra."
    echo "  Con -n 2 se evita el 'Affinity resource requirement cannot be met'"
    echo "  que aparecio con -n 4, y MEGAHIT apenas pierde velocidad."
    BSUB="bsub -J \"megahit${pre}[1-$N]\" \\
     -o \"$LOGS_DIR/megahit${pre}.%J.%I.log\" -e \"$LOGS_DIR/megahit${pre}.%J.%I.err\" \\
     -q normal -n 2 -M 2000 \\
     -R \"select[mem>2000] rusage[mem=2000] span[hosts=1]\" \\
     \"MEGAHIT_MEM=1 $SCRIPTS_DIR/assembly_megahit.sh \\
      $SCRATCH/asm_samples.txt $SCRATCH/unmapped_fastq '$pre' $out\""
    ;;

  trinity)
    lst="${variant:+$SCRATCH/trin_${variant}.txt}"
    lst="${lst:-$SCRATCH/asm_samples.txt}"
    chk_worker assembly_trinity.sh
    chk_dir "$SCRATCH/unmapped_fastq" "*_unmapped_R1.fastq.gz"
    chk_list "$lst"
    chk_env "$ENV_TRINITY"
    chk_writable "$SCRATCH/trinity"
    chk_writable "$RESULTS_DIR/assembly_stats_trinity"
    chk_queue long
    echo
    warn "MEDIDO: 754 MB y 1h50m para 1.26M pares. El cuello es el TIEMPO, no la RAM."
    warn "Trinity son HORAS por muestra, no minutos: cola 'long', no 'normal'."
    warn "Chrysalis crea decenas de miles de ficheros. El worker trabaja en disco"
    warn "     local y usa --full_cleanup; no lo desactives sin mirar la cuota de inodes."
    echo "  Cuota actual de inodes:"
    lfs quota -h -u "$USER" "$SCRATCH" 2>/dev/null | sed 's/^/    /' || echo "    (no disponible)"
    BSUB="bsub -J \"trinity[1-$N]%6\" \\
     -o \"$LOGS_DIR/trinity.%J.%I.log\" -e \"$LOGS_DIR/trinity.%J.%I.err\" \\
     -q long -n 8 -M 4000 \\
     -R \"select[mem>4000] rusage[mem=4000] span[hosts=1]\" \\
     \"TRINITY_MEM=3 $SCRIPTS_DIR/assembly_trinity.sh \\
      $lst $SCRATCH/unmapped_fastq\""
    ;;

  quastcmp)
    chk_worker quast_compare_worker.sh
    chk_list "$SCRATCH/asm_samples.txt" 22
    chk_env "$ENV_QUAST"
    chk_writable "$RESULTS_DIR/quast_compare"
    chk_queue normal
    echo
    echo "  Ensamblajes que encontrara para comparar:"
    for d in "$SCRATCH"/spades "$SCRATCH"/spades_* "$SCRATCH"/megahit "$SCRATCH"/megahit_* \
             "$SCRATCH"/trinity "$SCRATCH"/trinity_* "$SCRATCH"/union; do
      [[ -d "$d" ]] || continue
      nn=$(find "$d" -mindepth 2 -maxdepth 2 -name contigs.fasta -size +0 2>/dev/null | wc -l)
      (( nn > 0 )) && printf '    %-28s %3d muestras\n' "$(basename "$d")" "$nn"
    done
    warn "Sin -r (genoma de referencia): QUAST da contiguidad, NO errores de ensamblaje."
    BSUB="bsub -J \"quastcmp[1-$N]\" \\
     -o \"$LOGS_DIR/quastcmp.%J.%I.log\" -e \"$LOGS_DIR/quastcmp.%J.%I.err\" \\
     -q normal -n 2 -M 8000 \\
     -R \"select[mem>8000] rusage[mem=8000] span[hosts=1]\" \\
     \"$SCRIPTS_DIR/quast_compare_worker.sh $SCRATCH/asm_samples.txt\""
    ;;

  quast)
    chk_worker quast_worker.sh
    chk_dir "$SCRATCH/spades" "contigs.fasta"
    chk_list "$SCRATCH/assembly_list.txt"
    chk_env "$ENV_QUAST"
    chk_writable "$RESULTS_DIR/quast"
    echo
    echo "  Si la lista no existe:"
    echo "    find \$SCRATCH/spades -type f -name contigs.fasta | sort > \$SCRATCH/assembly_list.txt"
    BSUB="bsub -J \"quast[1-$N]%10\" \\
     -o \"$LOGS_DIR/quast.%J.%I.log\" -e \"$LOGS_DIR/quast.%J.%I.err\" \\
     -q normal -n 8 -M 8000 \\
     -R \"select[mem>8000] rusage[mem=8000] span[hosts=1]\" \\
     \"$SCRIPTS_DIR/quast_worker.sh $SCRATCH/assembly_list.txt\""
    ;;

  genomad|checkv)
    chk_worker "${step}_worker.sh"
    chk_dir "$SCRATCH/spades" "*_final.fasta"
    chk_list "$SCRATCH/fasta_list.txt"
    [[ "$step" == genomad ]] && { chk_env "$ENV_GENOMAD"; chk_dir "$GENOMAD_DB"; nn=8; MM=32000; pct=4; }
    [[ "$step" == checkv  ]] && { chk_env "$ENV_CHECKV";  chk_dir "$CHECKV_DB";  nn=8; MM=16000; pct=6; }
    chk_writable "$RESULTS_DIR/$step"
    echo
    echo "  Si la lista no existe:"
    echo "    find \$SCRATCH/spades -type f -name '*_final.fasta' | sort > \$SCRATCH/fasta_list.txt"
    warn "OJO: check_step.sh deriva la clave como 'M1_S1' pero estos workers nombran"
    warn "     la salida con basename .fasta = 'M1_S1_final'. No cuadran (pendiente)."
    BSUB="bsub -J \"${step}[1-$N]%${pct}\" \\
     -o \"$LOGS_DIR/${step}.%J.%I.log\" -e \"$LOGS_DIR/${step}.%J.%I.err\" \\
     -q normal -n ${nn} -M ${MM} \\
     -R \"select[mem>${MM}] rusage[mem=${MM}] span[hosts=1]\" \\
     \"$SCRIPTS_DIR/${step}_worker.sh $SCRATCH/fasta_list.txt\""
    ;;

  diamond)
    chk_worker diamond_blastx_worker.sh
    chk_dir "$SCRATCH/spades" "contigs_1000bp.fasta"
    chk_list "$SCRATCH/query_list.txt"
    chk_env "$ENV_DIAMOND"
    chk_file "$DIAMOND_DB" "base DIAMOND"
    chk_writable "$RESULTS_DIR/diamond"
    echo
    echo "  Si la lista no existe:"
    echo "    find \$SCRATCH/spades -type f -name contigs_1000bp.fasta | sort > \$SCRATCH/query_list.txt"
    BSUB="bsub -J \"dmnd[1-$N]%6\" \\
     -o \"$LOGS_DIR/dmnd.%J.%I.log\" -e \"$LOGS_DIR/dmnd.%J.%I.err\" \\
     -q normal -n 8 -M 16000 \\
     -R \"select[mem>16000] rusage[mem=16000] span[hosts=1]\" \\
     \"$SCRIPTS_DIR/diamond_blastx_worker.sh $SCRATCH/query_list.txt\""
    ;;

  *) echo "Paso desconocido: $step"; echo; usage; exit 1 ;;
esac

# --- Cola limpia ------------------------------------------------------------
echo
running=$(bjobs -w 2>/dev/null | grep -vc "^JOBID" || true)
if [[ "${running:-0}" -gt 0 ]]; then
  warn "tienes $running trabajos en cola. Confirma que el paso anterior TERMINO:"
  bjobs -w 2>/dev/null | head -5
else
  ok "no hay trabajos en cola"
fi

# --- Resultado --------------------------------------------------------------
echo
echo "==========================================================="
if (( ERR > 0 )); then
  echo "$ERR COMPROBACIONES FALLIDAS — no envies nada hasta resolverlas."
  echo "==========================================================="
  exit 1
fi
echo "Todo listo. Comando (N=${N:-?}):"
echo "==========================================================="
echo
echo "$BSUB"
echo
echo "-----------------------------------------------------------"
echo "LEE la cadena de arriba antes de copiarla. Ninguna ruta debe"
echo "empezar por '/worker' ni quedar coja: eso es una variable vacia."
echo "-----------------------------------------------------------"