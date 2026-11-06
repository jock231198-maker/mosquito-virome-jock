#!/usr/bin/env bash
# =============================================================================
# isv_sanity.sh -- Chequeo de sanidad del viroma MERIDA contra ISVs conocidos
#
#   Dos pistas INDEPENDIENTES contra el MISMO set de referencias:
#     (1) READS   -> bowtie2 --very-sensitive-local -> breadth/depth por genoma
#     (2) CONTIGS -> blastn (dc-megablast) -> identidad y fraccion reconstruida
#   Se cruzan al final. NO se mezclan contigs y reads en un mismo indice.
#
#   Virus objetivo: PCLV (3 segmentos), HTV (segmentado), Aedes anphevirus.
#   Controles: CFAV (positivo esperado), TMV (negativo, mide ruido de mapeo),
#              y tu propio anfevirus ensamblado (auto-referencia).
#
# Subcomandos:
#   check      verifica binarios, rutas y muestras; no corre nada
#   refs       baja/arma el set de referencias + indices (correr en head node)
#   map        lanza el array LSF de mapeo de reads (22 muestras)
#   contigs    blastn de contigs/vOTUs contra las referencias
#   table      junta todo en TSVs finales
#
# Uso tipico:
#   ./isv_sanity.sh check
#   ./isv_sanity.sh refs
#   ./isv_sanity.sh map
#   ./isv_sanity.sh contigs
#   ./isv_sanity.sh table
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
[[ -f "${SCRIPT_DIR}/config.sh" ]] && source "${SCRIPT_DIR}/config.sh"

SCRATCH_ROOT="${SCRATCH_ROOT:-/lustre/scratch126/tol/teams/lawniczak/users/jm79}"
REFS_DIR="${REFS_DIR:-${SCRATCH_ROOT}/refs}"
SCRATCH="${SCRATCH:-${SCRATCH_ROOT}/work/MERIDA}"
RESULTS_DIR="${RESULTS_DIR:-${SCRATCH_ROOT}/results/MERIDA}"
LOGS_DIR="${LOGS_DIR:-${SCRATCH}/logs}"

# --- Salidas propias de este chequeo ---------------------------------------
ISV_DIR="${REFS_DIR}/isv_sanity"                 # referencias + indices
OUT_DIR="${RESULTS_DIR}/isv_sanity"              # resultados
BAM_DIR="${OUT_DIR}/bam"
COV_DIR="${OUT_DIR}/coverage"
BLAST_DIR="${OUT_DIR}/blast"
mkdir -p "${ISV_DIR}" "${BAM_DIR}" "${COV_DIR}" "${BLAST_DIR}" "${LOGS_DIR}"

# --- Entradas ---------------------------------------------------------------
# OJO: elige la fuente de reads a conciencia (ver nota EVE mas abajo).
#   unmapped_fastq  = reads tras depletar host  -> quita EVEs transcritos
#   trimmed         = reads crudos recortados   -> incluye EVEs + virus
READ_SET="${READ_SET:-unmapped}"                 # unmapped | trimmed
READS_UNMAPPED_DIR="${READS_UNMAPPED_DIR:-${SCRATCH}/unmapped_fastq}"
READS_TRIMMED_DIR="${READS_TRIMMED_DIR:-${SCRATCH}/trimmed}"
SAMPLES_FILE="${SAMPLES_FILE:-${SCRIPT_DIR}/samples.txt}"   # una muestra por linea

# Contigs a interrogar (edita si cambian las rutas)
CONTIGS_UNION="${CONTIGS_UNION:-${RESULTS_DIR}/derep/union_min500.fasta}"
CONTIGS_VOTU="${CONTIGS_VOTU:-${RESULTS_DIR}/votus/votus_95_85.fasta}"

# Tu propio anfevirus completo (auto-referencia). Deja vacio para omitirlo.
OWN_ANPHE_FASTA="${OWN_ANPHE_FASTA:-}"

# --- Parametros -------------------------------------------------------------
THREADS="${THREADS:-8}"
MEM_MB="${MEM_MB:-8000}"
QUEUE="${QUEUE:-normal}"
BT2_PRESET="${BT2_PRESET:---very-sensitive-local}"
MIN_ALN_ID="${MIN_ALN_ID:-0.85}"   # identidad minima por read (filtro NM/alnlen)
JOBNAME="${JOBNAME:-isv_sanity}"

RED=$'\033[31m'; GRN=$'\033[32m'; YLW=$'\033[33m'; RST=$'\033[0m'
log()  { printf '%s[%s]%s %s\n' "${GRN}" "$(date +%H:%M:%S)" "${RST}" "$*"; }
warn() { printf '%s[WARN]%s %s\n' "${YLW}" "${RST}" "$*" >&2; }
die()  { printf '%s[ERROR]%s %s\n' "${RED}" "${RST}" "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# MANIFEST de referencias:  label <TAB> consulta Entrez <TAB> n_segmentos_esperados
# n_segmentos_esperados = 0 significa "no se cuantos, solo reporta"
# ---------------------------------------------------------------------------
MANIFEST="${ISV_DIR}/refs_manifest.tsv"
write_manifest() {
  cat > "${MANIFEST}" <<'EOF'
PCLV	"Phasi Charoen-like virus"[Organism] AND srcdb_refseq[PROP]	3
HTV	"Humaita-Tubiacanga virus"[Organism] AND complete[Title]	0
AeAV	"Aedes anphevirus"[Organism] AND srcdb_refseq[PROP]	1
CFAV	"Cell fusing agent virus"[Organism] AND srcdb_refseq[PROP]	1
NEGCTRL	"Tobacco mosaic virus"[Organism] AND srcdb_refseq[PROP]	1
EOF
  log "MANIFEST escrito en ${MANIFEST} -- revisalo antes de seguir."
}

# ---------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------
cmd_check() {
  local fail=0
  log "Binarios:"
  for b in bowtie2 bowtie2-build samtools blastn makeblastdb; do
    if command -v "$b" >/dev/null 2>&1; then
      printf '  %-14s %s\n' "$b" "$(command -v "$b")"
    else
      printf '  %-14s %sFALTA%s\n' "$b" "${RED}" "${RST}"; fail=1
    fi
  done
  for b in esearch efetch seqkit bedtools; do
    command -v "$b" >/dev/null 2>&1 \
      && printf '  %-14s %s\n' "$b" "$(command -v "$b")" \
      || printf '  %-14s %sausente (opcional)%s\n' "$b" "${YLW}" "${RST}"
  done

  log "Rutas:"
  for p in "${REFS_DIR}" "${SCRATCH}" "${RESULTS_DIR}"; do
    [[ -d "$p" ]] && printf '  OK    %s\n' "$p" || { printf '  %sFALTA%s %s\n' "${RED}" "${RST}" "$p"; fail=1; }
  done

  local rdir; rdir="$(reads_dir)"
  [[ -d "${rdir}" ]] && printf '  OK    %s (READ_SET=%s)\n' "${rdir}" "${READ_SET}" \
                     || { printf '  %sFALTA%s %s\n' "${RED}" "${RST}" "${rdir}"; fail=1; }

  log "Contigs:"
  for f in "${CONTIGS_UNION}" "${CONTIGS_VOTU}"; do
    [[ -s "$f" ]] && printf '  OK    %s (%s seqs)\n' "$f" "$(grep -c '^>' "$f" || echo '?')" \
                  || printf '  %sFALTA%s %s\n' "${YLW}" "${RST}" "$f"
  done

  log "Muestras:"
  if [[ -s "${SAMPLES_FILE}" ]]; then
    printf '  %s muestras en %s\n' "$(grep -cve '^\s*$' "${SAMPLES_FILE}")" "${SAMPLES_FILE}"
    local missing=0
    while read -r s; do
      [[ -z "$s" ]] && continue
      read_pair "$s" >/dev/null 2>&1 || { printf '  %sSIN FASTQ%s %s\n' "${RED}" "${RST}" "$s"; missing=1; }
    done < "${SAMPLES_FILE}"
    [[ $missing -eq 0 ]] && printf '  todas con FASTQ R1/R2\n' || fail=1
  else
    printf '  %sFALTA%s %s -- crealo con una muestra por linea\n' "${RED}" "${RST}" "${SAMPLES_FILE}"; fail=1
  fi

  log "Referencias:"
  if [[ -s "${ISV_DIR}/isv_refs.fasta" ]]; then
    printf '  OK    %s\n' "${ISV_DIR}/isv_refs.fasta"
    grep '^>' "${ISV_DIR}/isv_refs.fasta" | sed 's/^/        /'
  else
    printf '  %sno construidas%s -- corre: %s refs\n' "${YLW}" "${RST}" "$0"
  fi

  [[ $fail -eq 0 ]] && log "check OK" || die "check con fallos (ver arriba)"
}

reads_dir() {
  case "${READ_SET}" in
    unmapped) echo "${READS_UNMAPPED_DIR}" ;;
    trimmed)  echo "${READS_TRIMMED_DIR}" ;;
    *) die "READ_SET debe ser 'unmapped' o 'trimmed'" ;;
  esac
}

# Devuelve "R1<TAB>R2" probando los patrones usuales del pipeline
read_pair() {
  local s="$1" d; d="$(reads_dir)"
  local r1 r2
  for pat in "${s}_1.fastq.gz:${s}_2.fastq.gz" \
             "${s}_R1.fastq.gz:${s}_R2.fastq.gz" \
             "${s}_unmapped_1.fastq.gz:${s}_unmapped_2.fastq.gz" \
             "${s}_val_1.fq.gz:${s}_val_2.fq.gz" \
             "${s}_1_val_1.fq.gz:${s}_2_val_2.fq.gz"; do
    r1="${d}/${pat%%:*}"; r2="${d}/${pat##*:}"
    if [[ -s "$r1" && -s "$r2" ]]; then printf '%s\t%s\n' "$r1" "$r2"; return 0; fi
  done
  return 1
}

# ---------------------------------------------------------------------------
# refs -- bajar y armar el set de referencias (CORRER EN HEAD NODE: usa red)
# ---------------------------------------------------------------------------
cmd_refs() {
  command -v esearch >/dev/null 2>&1 || die "entrez-direct no disponible: conda install -c bioconda entrez-direct"
  [[ -s "${MANIFEST}" ]] || write_manifest

  local out="${ISV_DIR}/isv_refs.fasta"
  : > "${out}.tmp"

  while IFS=$'\t' read -r label query nexp; do
    [[ -z "${label}" || "${label}" == \#* ]] && continue
    log "Bajando ${label} ..."
    local tmp="${ISV_DIR}/.${label}.fa"
    if ! esearch -db nuccore -query "${query}" < /dev/null | efetch -format fasta > "${tmp}" 2>/dev/null; then
      warn "${label}: la consulta fallo. Revisa el MANIFEST o baja el FASTA a mano a ${tmp}"
      continue
    fi
    local n; n="$(grep -c '^>' "${tmp}" || true)"
    if [[ "${n}" -eq 0 ]]; then
      warn "${label}: 0 secuencias. Ajusta la consulta en ${MANIFEST}"
      rm -f "${tmp}"; continue
    fi
    if [[ "${nexp}" -gt 0 && "${n}" -ne "${nexp}" ]]; then
      warn "${label}: esperaba ${nexp} secuencias, bajaron ${n}. Verifica antes de confiar."
    fi
    # Renombra: >LABEL__ACCESION  (descripcion original preservada como comentario)
    awk -v L="${label}" '
      /^>/ { acc=$1; sub(/^>/,"",acc); sub(/\..*$/,"",acc);
             desc=$0; sub(/^>[^ ]+ */,"",desc);
             print ">" L "__" acc " " desc; next } { print }' "${tmp}" >> "${out}.tmp"
    log "  ${label}: ${n} secuencia(s)"
    rm -f "${tmp}"
  done < "${MANIFEST}"

  # Auto-referencia: tu propio anfevirus ensamblado
  if [[ -n "${OWN_ANPHE_FASTA}" && -s "${OWN_ANPHE_FASTA}" ]]; then
    awk '/^>/ { n++; print ">AeAV_MERIDA__own" n; next } { print }' "${OWN_ANPHE_FASTA}" >> "${out}.tmp"
    log "  AeAV_MERIDA: auto-referencia anadida desde ${OWN_ANPHE_FASTA}"
  else
    warn "OWN_ANPHE_FASTA vacio: no se incluye tu anfevirus como auto-referencia"
  fi

  mv "${out}.tmp" "${out}"

  log "Composicion final del set de referencias:"
  if command -v seqkit >/dev/null 2>&1; then
    seqkit fx2tab -nl "${out}" | awk -F'\t' '{printf "  %-28s %s bp\n", $1, $NF}'
  else
    awk '/^>/{if(n)printf "  %-28s %s bp\n",h,n; h=substr($1,2); n=0; next}{n+=length($0)}
         END{if(n)printf "  %-28s %s bp\n",h,n}' "${out}"
  fi

  log "Construyendo indices ..."
  bowtie2-build --threads "${THREADS}" "${out}" "${ISV_DIR}/isv_refs" > "${LOGS_DIR}/bowtie2-build.log" 2>&1
  makeblastdb -in "${out}" -dbtype nucl -out "${ISV_DIR}/isv_refs_blastdb" > "${LOGS_DIR}/makeblastdb.log" 2>&1
  samtools faidx "${out}"
  log "refs listo en ${ISV_DIR}"
  warn "VERIFICA las longitudes de arriba contra la literatura antes de interpretar nada."
}

# ---------------------------------------------------------------------------
# map -- array LSF sobre las muestras
# ---------------------------------------------------------------------------
cmd_map() {
  [[ -s "${ISV_DIR}/isv_refs.1.bt2" ]] || die "indice bowtie2 ausente; corre '$0 refs'"
  [[ -s "${SAMPLES_FILE}" ]] || die "falta ${SAMPLES_FILE}"
  local n; n="$(grep -cve '^\s*$' "${SAMPLES_FILE}")"

  local runner="${ISV_DIR}/_map_one.sh"
  cat > "${runner}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
SAMPLE=\$(sed -n "\${LSB_JOBINDEX}p" "${SAMPLES_FILE}")
[[ -n "\${SAMPLE}" ]] || exit 0
source "${SCRIPT_DIR}/\$(basename "${BASH_SOURCE[0]}")" --lib 2>/dev/null || true
READ_SET="${READ_SET}"
EOF
  # Inyecta las funciones necesarias en el runner (evita depender de --lib)
  declare -f reads_dir read_pair >> "${runner}"
  cat >> "${runner}" <<EOF
READS_UNMAPPED_DIR="${READS_UNMAPPED_DIR}"
READS_TRIMMED_DIR="${READS_TRIMMED_DIR}"
PAIR=\$(read_pair "\${SAMPLE}") || { echo "sin FASTQ para \${SAMPLE}" >&2; exit 1; }
R1=\$(echo "\${PAIR}" | cut -f1); R2=\$(echo "\${PAIR}" | cut -f2)
BAM="${BAM_DIR}/\${SAMPLE}.${READ_SET}.bam"

bowtie2 ${BT2_PRESET} -p ${THREADS} -x "${ISV_DIR}/isv_refs" \\
        -1 "\${R1}" -2 "\${R2}" --no-unal \\
        2> "${LOGS_DIR}/bt2_\${SAMPLE}.${READ_SET}.log" \\
  | samtools view -@ 2 -b -F 0x900 - \\
  | samtools sort -@ 2 -o "\${BAM}" -
samtools index "\${BAM}"

# Filtro de identidad por read (descarta alineamientos espurios/divergentes)
samtools view -@ 2 -h "\${BAM}" \\
  | awk -v MINID=${MIN_ALN_ID} '
      /^@/ { print; next }
      {
        nm = -1; for (i=12; i<=NF; i++) if (\$i ~ /^NM:i:/) { split(\$i,a,":"); nm=a[3] }
        aln = 0; c = \$6; while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
          seg = substr(c, RSTART, RLENGTH); op = substr(seg, length(seg));
          len = substr(seg, 1, length(seg)-1);
          if (op ~ /[M=XI]/) aln += len;
          c = substr(c, RLENGTH+1);
        }
        if (nm < 0 || aln == 0) next
        if ((aln - nm) / aln >= MINID) print
      }' \\
  | samtools view -@ 2 -b -o "${BAM_DIR}/\${SAMPLE}.${READ_SET}.id.bam" -
samtools index "${BAM_DIR}/\${SAMPLE}.${READ_SET}.id.bam"

for TAG in "" ".id"; do
  B="${BAM_DIR}/\${SAMPLE}.${READ_SET}\${TAG}.bam"
  samtools coverage "\${B}" > "${COV_DIR}/\${SAMPLE}.${READ_SET}\${TAG}.coverage.tsv"
  samtools idxstats "\${B}" > "${COV_DIR}/\${SAMPLE}.${READ_SET}\${TAG}.idxstats.tsv"
  samtools depth -a -J "\${B}" \\
    | awk -v S="\${SAMPLE}" '{d[\$1]++; if(\$3>=10) t[\$1]++}
        END{for(r in d) printf "%s\\t%s\\t%.4f\\n", S, r, (t[r]+0)/d[r]}' \\
    > "${COV_DIR}/\${SAMPLE}.${READ_SET}\${TAG}.breadth10x.tsv"
done
echo "OK \${SAMPLE}"
EOF
  chmod +x "${runner}"

  log "Lanzando array LSF de ${n} muestras (READ_SET=${READ_SET})"
  bsub -J "${JOBNAME}[1-${n}]" -q "${QUEUE}" -n "${THREADS}" \
       -R "span[hosts=1] select[mem>${MEM_MB}] rusage[mem=${MEM_MB}]" -M "${MEM_MB}" \
       -o "${LOGS_DIR}/${JOBNAME}.%I.out" -e "${LOGS_DIR}/${JOBNAME}.%I.err" \
       "${runner}"
}

# ---------------------------------------------------------------------------
# contigs -- blastn de contigs/vOTUs contra las referencias
# ---------------------------------------------------------------------------
cmd_contigs() {
  [[ -s "${ISV_DIR}/isv_refs_blastdb.nin" ]] || die "blastdb ausente; corre '$0 refs'"
  local fmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs'

  for pair in "union:${CONTIGS_UNION}" "votu:${CONTIGS_VOTU}"; do
    local tag="${pair%%:*}" fa="${pair#*:}"
    [[ -s "${fa}" ]] || { warn "omito ${tag}: no existe ${fa}"; continue; }
    log "blastn (dc-megablast) ${tag} -> referencias"
    blastn -task dc-megablast -query "${fa}" -db "${ISV_DIR}/isv_refs_blastdb" \
           -outfmt "${fmt}" -evalue 1e-10 -num_threads "${THREADS}" \
           -max_target_seqs 50 > "${BLAST_DIR}/${tag}_vs_isv.tsv"

    # Fraccion de cada referencia reconstruida por contigs (union de HSPs)
    if command -v bedtools >/dev/null 2>&1; then
      awk -F'\t' '$3>=70 && $4>=100 {s=$9; e=$10; if(s>e){t=s;s=e;e=t} print $2"\t"(s-1)"\t"e}' \
          "${BLAST_DIR}/${tag}_vs_isv.tsv" \
        | sort -k1,1 -k2,2n | bedtools merge \
        | awk -F'\t' '{cov[$1]+=$3-$2} END{for(r in cov) printf "%s\t%d\n", r, cov[r]}' \
        | sort > "${BLAST_DIR}/${tag}_refcov.tmp"
      join -t $'\t' "${BLAST_DIR}/${tag}_refcov.tmp" \
           <(cut -f1,2 "${ISV_DIR}/isv_refs.fasta.fai" | sort) \
        | awk -F'\t' 'BEGIN{print "ref\tbases_cubiertas\tlargo_ref\tfraccion"}
                      {printf "%s\t%d\t%d\t%.4f\n", $1, $2, $3, $2/$3}' \
        > "${BLAST_DIR}/${tag}_refcov.tsv"
      rm -f "${BLAST_DIR}/${tag}_refcov.tmp"
      log "  -> ${BLAST_DIR}/${tag}_refcov.tsv"
    else
      warn "bedtools ausente: omito la fraccion reconstruida de ${tag}"
    fi
  done
}

# ---------------------------------------------------------------------------
# table -- consolidar
# ---------------------------------------------------------------------------
cmd_table() {
  local out="${OUT_DIR}/isv_sanity_reads.${READ_SET}.tsv"
  printf 'muestra\tfiltro\tref\tvirus\tlargo_ref\treads\tbreadth_1x\tmeandepth\tmeanmapq\n' > "${out}"
  shopt -s nullglob
  for f in "${COV_DIR}"/*.${READ_SET}*.coverage.tsv; do
    local base; base="$(basename "$f" .coverage.tsv)"
    local sample="${base%%.*}"
    local filt="raw"; [[ "${base}" == *".id" ]] && filt="id${MIN_ALN_ID}"
    awk -F'\t' -v S="${sample}" -v F="${filt}" 'NR>1 {
        split($1, p, "__");
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", S, F, $1, p[1], $3, $4, $6, $7, $9
      }' "$f" >> "${out}"
  done
  shopt -u nullglob
  log "Tabla de reads -> ${out}"
  [[ -s "${BLAST_DIR}/union_vs_isv.tsv" ]] && log "Tabla de contigs -> ${BLAST_DIR}/"

  log "Resumen rapido (breadth_1x medio por virus, filtro id):"
  awk -F'\t' 'NR>1 && $2!="raw" {s[$4]+=$7; n[$4]++}
              END{for(v in s) printf "  %-14s %6.1f%%\n", v, s[v]/n[v]}' "${out}" | sort
}

# ---------------------------------------------------------------------------
case "${1:-}" in
  check)    cmd_check ;;
  refs)     cmd_refs ;;
  map)      cmd_map ;;
  contigs)  cmd_contigs ;;
  table)    cmd_table ;;
  manifest) write_manifest ;;
  --lib)    : ;;
  *) sed -n '2,30p' "${BASH_SOURCE[0]}"; exit 1 ;;
esac
