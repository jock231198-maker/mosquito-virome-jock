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
#   check      verifica entornos conda, binarios, rutas y muestras
#   samples    genera samples.txt desde el directorio de reads
#   refs       baja/arma referencias + indices  (CORRER EN HEAD NODE: usa red)
#   map        lanza el array LSF de mapeo de reads
#   contigs    blastn de contigs/vOTUs contra las referencias
#   table      junta todo en TSVs finales
#
# Dependencias: solo lo que ya existe en el farm.
#   bowtie2 + samtools -> $ENV_BOWTIE2 ; blastn + makeblastdb -> $ENV_BLAST
#   curl (descarga NCBI E-utilities), awk, sort, join. Sin entrez-direct,
#   sin bedtools, sin seqkit.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/config.sh"

# --- Entornos conda ---------------------------------------------------------
# ENV_BOWTIE2 ya trae samtools (segun config.sh). ENV_BLAST hay que confirmarlo:
# ajusta aqui o en config.sh si blastn vive en otro entorno.
ENV_BLAST="${ENV_BLAST:-${ENV_VIRALDB:-viraldb}}"

# --- Salidas propias de este chequeo ---------------------------------------
ISV_DIR="${REFS_DIR}/isv_sanity"
OUT_DIR="${RESULTS_DIR}/isv_sanity"
BAM_DIR="${OUT_DIR}/bam"
COV_DIR="${OUT_DIR}/coverage"
BLAST_DIR="${OUT_DIR}/blast"
mkdir -p "${ISV_DIR}" "${BAM_DIR}" "${COV_DIR}" "${BLAST_DIR}" "${LOGS_DIR}"

# --- Entradas ---------------------------------------------------------------
#   unmapped = reads tras depletar host -> quita EVEs transcritos
#   trimmed  = reads recortados          -> incluye EVEs + virus
READ_SET="${READ_SET:-unmapped}"
READS_UNMAPPED_DIR="${READS_UNMAPPED_DIR:-${SCRATCH}/unmapped_fastq}"
READS_TRIMMED_DIR="${READS_TRIMMED_DIR:-${SCRATCH}/trimmed/${RESULT_FROM}}"
SAMPLES_FILE="${SAMPLES_FILE:-${ISV_DIR}/samples.txt}"

# Contigs a interrogar -- PENDIENTE: confirmar rutas reales del set _min500
CONTIGS_UNION="${CONTIGS_UNION:-}"
CONTIGS_VOTU="${CONTIGS_VOTU:-}"

# Tu propio anfevirus completo (auto-referencia). Vacio = omitir.
OWN_ANPHE_FASTA="${OWN_ANPHE_FASTA:-}"

# --- Parametros -------------------------------------------------------------
THREADS="${THREADS:-8}"
MEM_MB="${MEM_MB:-8000}"
QUEUE="${QUEUE:-normal}"
BT2_PRESET="${BT2_PRESET:---very-sensitive-local}"
MIN_ALN_ID="${MIN_ALN_ID:-0.85}"
JOBNAME="${JOBNAME:-isv_sanity}"

RED=$'\033[31m'; GRN=$'\033[32m'; YLW=$'\033[33m'; RST=$'\033[0m'
log()  { printf '%s[%s]%s %s\n' "${GRN}" "$(date +%H:%M:%S)" "${RST}" "$*"; }
warn() { printf '%s[WARN]%s %s\n' "${YLW}" "${RST}" "$*" >&2; }
die()  { printf '%s[ERROR]%s %s\n' "${RED}" "${RST}" "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# MANIFEST:  label <TAB> consulta Entrez <TAB> n_segmentos_esperados (0 = no se)
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
# Helpers de reads
# ---------------------------------------------------------------------------
reads_dir() {
  case "${READ_SET}" in
    unmapped) echo "${READS_UNMAPPED_DIR}" ;;
    trimmed)  echo "${READS_TRIMMED_DIR}" ;;
    *) die "READ_SET debe ser 'unmapped' o 'trimmed'" ;;
  esac
}

# Devuelve "R1<TAB>R2". Patron real del farm: M1_S1_unmapped_R1.fastq.gz
read_pair() {
  local s="$1" d r1 r2; d="$(reads_dir)"
  for pat in "${s}_unmapped_R1.fastq.gz:${s}_unmapped_R2.fastq.gz" \
             "${s}_R1.fastq.gz:${s}_R2.fastq.gz" \
             "${s}_1.fastq.gz:${s}_2.fastq.gz" \
             "${s}_val_1.fq.gz:${s}_val_2.fq.gz" \
             "${s}_R1_val_1.fq.gz:${s}_R2_val_2.fq.gz" \
             "${s}_trimmed_R1.fastq.gz:${s}_trimmed_R2.fastq.gz"; do
    r1="${d}/${pat%%:*}"; r2="${d}/${pat##*:}"
    if [[ -s "$r1" && -s "$r2" ]]; then printf '%s\t%s\n' "$r1" "$r2"; return 0; fi
  done
  return 1
}

cmd_samples() {
  local d; d="$(reads_dir)"
  [[ -d "$d" ]] || die "no existe ${d}"
  find "$d" -maxdepth 1 -name '*_R1.fastq.gz' -o -maxdepth 1 -name '*_R1.fq.gz' \
    | xargs -r -n1 basename \
    | sed -E 's/_(unmapped|trimmed)?_?R1\.f(ast)?q\.gz$//' \
    | sort -u > "${SAMPLES_FILE}"
  log "$(wc -l < "${SAMPLES_FILE}") muestras -> ${SAMPLES_FILE}"
  cat "${SAMPLES_FILE}" | sed 's/^/  /'
}

# ---------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------
cmd_check() {
  local fail=0

  log "Entorno ${ENV_BOWTIE2} (bowtie2 + samtools):"
  if ( activate_env "${ENV_BOWTIE2}" >/dev/null 2>&1
       command -v bowtie2 >/dev/null && command -v samtools >/dev/null ); then
    ( activate_env "${ENV_BOWTIE2}" >/dev/null 2>&1
      printf '  bowtie2    %s\n' "$(bowtie2 --version | head -1 | awk '{print $NF}')"
      printf '  samtools   %s\n' "$(samtools --version | head -1 | awk '{print $NF}')"
      samtools coverage 2>&1 | grep -q 'Usage' \
        || printf '  %ssamtools coverage no disponible (samtools muy viejo)%s\n' "${RED}" "${RST}" )
  else
    printf '  %sno se pudo activar o faltan binarios%s\n' "${RED}" "${RST}"; fail=1
  fi

  log "Entorno ${ENV_BLAST} (blastn + makeblastdb):"
  if ( activate_env "${ENV_BLAST}" >/dev/null 2>&1
       command -v blastn >/dev/null && command -v makeblastdb >/dev/null ); then
    ( activate_env "${ENV_BLAST}" >/dev/null 2>&1
      printf '  blastn     %s\n' "$(blastn -version | head -1 | awk '{print $NF}')" )
  else
    printf '  %sblastn no esta en %s -- busca con: conda env list; luego ENV_BLAST=<env>%s\n' \
           "${RED}" "${ENV_BLAST}" "${RST}"; fail=1
  fi

  log "Red (NCBI E-utilities):"
  if curl -sf --max-time 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi" >/dev/null; then
    printf '  OK    alcanzable desde este nodo\n'
  else
    printf '  %ssin salida a NCBI%s -- baja los FASTA a mano a %s\n' "${YLW}" "${RST}" "${ISV_DIR}"
  fi

  log "Rutas:"
  for p in "${REFS_DIR}" "${SCRATCH}" "${RESULTS_DIR}"; do
    [[ -d "$p" ]] && printf '  OK    %s\n' "$p" \
                  || { printf '  %sFALTA%s %s\n' "${RED}" "${RST}" "$p"; fail=1; }
  done
  local rdir; rdir="$(reads_dir)"
  [[ -d "${rdir}" ]] && printf '  OK    %s (READ_SET=%s)\n' "${rdir}" "${READ_SET}" \
                     || { printf '  %sFALTA%s %s\n' "${RED}" "${RST}" "${rdir}"; fail=1; }

  log "Contigs:"
  for f in "${CONTIGS_UNION}" "${CONTIGS_VOTU}"; do
    if [[ -z "$f" ]]; then
      printf '  %ssin definir%s -- exporta CONTIGS_UNION / CONTIGS_VOTU\n' "${YLW}" "${RST}"
    elif [[ -s "$f" ]]; then
      printf '  OK    %s (%s seqs)\n' "$f" "$(grep -c '^>' "$f")"
    else
      printf '  %sFALTA%s %s\n' "${RED}" "${RST}" "$f"; fail=1
    fi
  done

  log "Muestras:"
  if [[ -s "${SAMPLES_FILE}" ]]; then
    printf '  %s en %s\n' "$(grep -cve '^[[:space:]]*$' "${SAMPLES_FILE}")" "${SAMPLES_FILE}"
    local missing=0
    while read -r s; do
      [[ -z "$s" ]] && continue
      read_pair "$s" >/dev/null 2>&1 || { printf '  %sSIN FASTQ%s %s\n' "${RED}" "${RST}" "$s"; missing=1; }
    done < "${SAMPLES_FILE}"
    [[ $missing -eq 0 ]] && printf '  todas con R1/R2\n' || fail=1
  else
    printf '  %sfalta%s %s -- genera con: %s samples\n' "${YLW}" "${RST}" "${SAMPLES_FILE}" "$0"
  fi

  log "Referencias:"
  if [[ -s "${ISV_DIR}/isv_refs.fasta" ]]; then
    grep '^>' "${ISV_DIR}/isv_refs.fasta" | sed 's/^/  /'
  else
    printf '  %sno construidas%s -- corre: %s refs\n' "${YLW}" "${RST}" "$0"
  fi

  [[ $fail -eq 0 ]] && log "check OK" || die "check con fallos (ver arriba)"
}

# ---------------------------------------------------------------------------
# refs -- descarga via E-utilities con curl (sin entrez-direct)
# ---------------------------------------------------------------------------
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

ncbi_fetch() {   # $1 = consulta Entrez -> FASTA a stdout
  local q="$1" ids
  ids=$(curl -sf --max-time 60 -G "${EUTILS}/esearch.fcgi" \
          --data-urlencode "db=nuccore" \
          --data-urlencode "term=${q}" \
          --data-urlencode "retmax=20" \
        | grep -o '<Id>[0-9]*</Id>' | sed 's/<[^>]*>//g' | paste -sd,)
  [[ -z "${ids}" ]] && return 1
  sleep 1   # cortesia con NCBI
  curl -sf --max-time 120 -G "${EUTILS}/efetch.fcgi" \
       --data-urlencode "db=nuccore" \
       --data-urlencode "id=${ids}" \
       --data-urlencode "rettype=fasta" \
       --data-urlencode "retmode=text"
}

cmd_refs() {
  [[ -s "${MANIFEST}" ]] || write_manifest
  local out="${ISV_DIR}/isv_refs.fasta"
  : > "${out}.tmp"

  while IFS=$'\t' read -r label query nexp; do
    [[ -z "${label}" || "${label}" == \#* ]] && continue
    log "Bajando ${label} ..."
    local tmp="${ISV_DIR}/.${label}.fa"
    if ! ncbi_fetch "${query}" > "${tmp}"; then
      warn "${label}: consulta sin resultados o sin red. Baja el FASTA a mano a ${tmp} y repite."
      rm -f "${tmp}"; continue
    fi
    local n; n="$(grep -c '^>' "${tmp}" || true)"
    if [[ "${n}" -eq 0 ]]; then
      warn "${label}: 0 secuencias. Ajusta la consulta en ${MANIFEST}"; rm -f "${tmp}"; continue
    fi
    [[ "${nexp}" -gt 0 && "${n}" -ne "${nexp}" ]] \
      && warn "${label}: esperaba ${nexp} secuencias, bajaron ${n}. Verifica antes de confiar."
    awk -v L="${label}" '
      /^>/ { acc=$1; sub(/^>/,"",acc); sub(/\..*$/,"",acc);
             desc=$0; sub(/^>[^ ]+ */,"",desc);
             print ">" L "__" acc " " desc; next } { print }' "${tmp}" >> "${out}.tmp"
    log "  ${label}: ${n} secuencia(s)"
    rm -f "${tmp}"
  done < "${MANIFEST}"

  if [[ -n "${OWN_ANPHE_FASTA}" && -s "${OWN_ANPHE_FASTA}" ]]; then
    awk '/^>/ { n++; print ">AeAV_MERIDA__own" n; next } { print }' "${OWN_ANPHE_FASTA}" >> "${out}.tmp"
    log "  AeAV_MERIDA: auto-referencia desde ${OWN_ANPHE_FASTA}"
  else
    warn "OWN_ANPHE_FASTA vacio: sin auto-referencia de tu propio anfevirus"
  fi

  [[ -s "${out}.tmp" ]] || die "no se bajo ninguna referencia"
  mv "${out}.tmp" "${out}"

  log "Composicion del set (VERIFICA estas longitudes contra la literatura):"
  awk '/^>/{if(h)printf "  %-30s %7d bp\n",h,n; h=substr($1,2); n=0; next}
       {n+=length($0)} END{if(h)printf "  %-30s %7d bp\n",h,n}' "${out}"
  echo "  esperado: PCLV_L ~6.8kb, PCLV_M ~3.9kb, PCLV_S ~1.4kb, AeAV ~12.4-12.9kb"

  log "Construyendo indices ..."
  ( activate_env "${ENV_BOWTIE2}" >/dev/null
    bowtie2-build --threads "${THREADS}" "${out}" "${ISV_DIR}/isv_refs" > "${LOGS_DIR}/isv_bt2build.log" 2>&1
    samtools faidx "${out}" )
  ( activate_env "${ENV_BLAST}" >/dev/null
    makeblastdb -in "${out}" -dbtype nucl -out "${ISV_DIR}/isv_refs_blastdb" > "${LOGS_DIR}/isv_makeblastdb.log" 2>&1 )
  log "refs listo en ${ISV_DIR}"
}

# ---------------------------------------------------------------------------
# map -- array LSF
# ---------------------------------------------------------------------------
cmd_map() {
  [[ -s "${ISV_DIR}/isv_refs.1.bt2" ]] || die "indice bowtie2 ausente; corre '$0 refs'"
  [[ -s "${SAMPLES_FILE}" ]] || die "falta ${SAMPLES_FILE}; corre '$0 samples'"
  local n; n="$(grep -cve '^[[:space:]]*$' "${SAMPLES_FILE}")"
  local self="${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"
  local runner="${ISV_DIR}/_map_one.sh"

  cat > "${runner}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export READ_SET="${READ_SET}"
export SAMPLES_FILE="${SAMPLES_FILE}"
source "${self}" --lib
SAMPLE=\$(sed -n "\${LSB_JOBINDEX}p" "${SAMPLES_FILE}")
[[ -n "\${SAMPLE}" ]] || exit 0

activate_env "${ENV_BOWTIE2}"
PAIR=\$(read_pair "\${SAMPLE}") || { echo "sin FASTQ para \${SAMPLE}" >&2; exit 1; }
R1=\$(echo "\${PAIR}" | cut -f1); R2=\$(echo "\${PAIR}" | cut -f2)
BAM="${BAM_DIR}/\${SAMPLE}.${READ_SET}.bam"

bowtie2 ${BT2_PRESET} -p ${THREADS} -x "${ISV_DIR}/isv_refs" \\
        -1 "\${R1}" -2 "\${R2}" --no-unal \\
        2> "${LOGS_DIR}/bt2_\${SAMPLE}.${READ_SET}.log" \\
  | samtools view -@ 2 -b -F 0x900 - \\
  | samtools sort -@ 2 -o "\${BAM}" -
samtools index "\${BAM}"

# Filtro de identidad por read: descarta alineamientos locales espurios
samtools view -@ 2 -h "\${BAM}" \\
  | awk -v MINID=${MIN_ALN_ID} '
      /^@/ { print; next }
      { nm=-1; for(i=12;i<=NF;i++) if(\$i ~ /^NM:i:/){split(\$i,a,":"); nm=a[3]}
        aln=0; c=\$6;
        while (match(c, /^[0-9]+[MIDNSHP=X]/)) {
          seg=substr(c,RSTART,RLENGTH); op=substr(seg,length(seg));
          len=substr(seg,1,length(seg)-1);
          if (op ~ /[M=XI]/) aln+=len;
          c=substr(c,RLENGTH+1);
        }
        if (nm<0 || aln==0) next
        if ((aln-nm)/aln >= MINID) print }' \\
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

  log "Array LSF de ${n} muestras (READ_SET=${READ_SET})"
  bsub -J "${JOBNAME}[1-${n}]" -q "${QUEUE}" -n "${THREADS}" \
       -R "span[hosts=1] select[mem>${MEM_MB}] rusage[mem=${MEM_MB}]" -M "${MEM_MB}" \
       -o "${LOGS_DIR}/${JOBNAME}.%I.out" -e "${LOGS_DIR}/${JOBNAME}.%I.err" \
       "${runner}"
}

# ---------------------------------------------------------------------------
# contigs
# ---------------------------------------------------------------------------
cmd_contigs() {
  [[ -s "${ISV_DIR}/isv_refs_blastdb.nin" ]] || die "blastdb ausente; corre '$0 refs'"
  local fmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs'

  for pair in "union:${CONTIGS_UNION}" "votu:${CONTIGS_VOTU}"; do
    local tag="${pair%%:*}" fa="${pair#*:}"
    [[ -n "${fa}" && -s "${fa}" ]] || { warn "omito ${tag}: ruta vacia o inexistente"; continue; }
    log "blastn dc-megablast: ${tag} -> referencias"
    ( activate_env "${ENV_BLAST}" >/dev/null
      blastn -task dc-megablast -query "${fa}" -db "${ISV_DIR}/isv_refs_blastdb" \
             -outfmt "${fmt}" -evalue 1e-10 -num_threads "${THREADS}" \
             -max_target_seqs 50 > "${BLAST_DIR}/${tag}_vs_isv.tsv" )

    # Fraccion de cada referencia reconstruida por contigs (merge de HSPs en awk)
    awk -F'\t' '$3>=70 && $4>=100 {s=$9; e=$10; if(s>e){t=s;s=e;e=t} print $2"\t"(s-1)"\t"e}' \
        "${BLAST_DIR}/${tag}_vs_isv.tsv" \
      | sort -k1,1 -k2,2n \
      | awk -F'\t' '
          { if (chr != $1)        { if (chr!="") cov[chr]+=e-s; chr=$1; s=$2; e=$3 }
            else if ($2 > e)      { cov[chr]+=e-s; s=$2; e=$3 }
            else if ($3 > e)      { e=$3 } }
          END { if (chr!="") cov[chr]+=e-s; for(c in cov) printf "%s\t%d\n", c, cov[c] }' \
      | sort > "${BLAST_DIR}/.${tag}_refcov.tmp"

    join -t $'\t' "${BLAST_DIR}/.${tag}_refcov.tmp" \
         <(cut -f1,2 "${ISV_DIR}/isv_refs.fasta.fai" | sort) \
      | awk -F'\t' 'BEGIN{print "ref\tbases_cubiertas\tlargo_ref\tfraccion"}
                    {printf "%s\t%d\t%d\t%.4f\n", $1, $2, $3, $2/$3}' \
      > "${BLAST_DIR}/${tag}_refcov.tsv"
    rm -f "${BLAST_DIR}/.${tag}_refcov.tmp"
    log "  -> ${BLAST_DIR}/${tag}_refcov.tsv"
  done
}

# ---------------------------------------------------------------------------
# table
# ---------------------------------------------------------------------------
cmd_table() {
  local out="${OUT_DIR}/isv_sanity_reads.${READ_SET}.tsv"
  printf 'muestra\tfiltro\tref\tvirus\tlargo_ref\treads\tbreadth_1x\tmeandepth\tmeanmapq\n' > "${out}"
  shopt -s nullglob
  for f in "${COV_DIR}"/*."${READ_SET}"*.coverage.tsv; do
    local base sample filt
    base="$(basename "$f" .coverage.tsv)"
    sample="${base%%.*}"
    filt="raw"; [[ "${base}" == *".id" ]] && filt="id${MIN_ALN_ID}"
    awk -F'\t' -v S="${sample}" -v F="${filt}" 'NR>1 {
        split($1, p, "__");
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", S, F, $1, p[1], $3, $4, $6, $7, $9
      }' "$f" >> "${out}"
  done
  shopt -u nullglob
  log "Tabla de reads -> ${out}"

  log "Breadth_1x medio por virus (filtro de identidad):"
  awk -F'\t' 'NR>1 && $2!="raw" {s[$4]+=$7; n[$4]++}
              END{for(v in s) printf "  %-16s %6.1f%%\n", v, s[v]/n[v]}' "${out}" | sort
}

# ---------------------------------------------------------------------------
case "${1:-}" in
  check)    cmd_check ;;
  samples)  cmd_samples ;;
  refs)     cmd_refs ;;
  map)      cmd_map ;;
  contigs)  cmd_contigs ;;
  table)    cmd_table ;;
  manifest) write_manifest ;;
  --lib)    : ;;
  *) sed -n '2,25p' "${BASH_SOURCE[0]}"; exit 1 ;;
esac
