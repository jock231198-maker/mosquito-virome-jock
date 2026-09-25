#!/bin/bash
# ---------------------------------------------------------------------------
# build_viral_db.sh — base viral propia: RefSeq viral + GenBank con huésped artrópodo
#
# Reemplaza a viral_db_ftp.sh, viral_db_edirect_txid*.sh y merge_db.sh.
# Una sola descarga (NCBI Datasets) alimenta tres productos:
#
#   blast/viral_nt       blastn de contigs/vOTUs        (con taxid -> staxids/sscinames)
#   diamond/viral_prot   DIAMOND blastx propio          (con --taxonmap/nodes/names)
#   mapping/viral_nr95   bowtie2 de reads               (completos, CD-HIT 95/85)
#
# Uso (cada paso es un envío aparte, igual que build_downstream_dbs.sh):
#   ./build_viral_db.sh check      # herramientas, red, cuota.       EMPIEZA AQUI
#   ./build_viral_db.sh probe      # 1 genoma de prueba: valida el formato de Datasets
#   ./build_viral_db.sh download   # RED. RefSeq + GenBank(host) + taxdump + taxdb
#   ./build_viral_db.sh prepare    # descomprime, metadatos, dedup, mapas de taxid
#   ./build_viral_db.sh blastn     # makeblastdb + autocomprobacion
#   ./build_viral_db.sh diamond    # diamond makedb + autocomprobacion
#   ./build_viral_db.sh mapping    # CD-HIT-EST 95/85 + bowtie2-build
#   ./build_viral_db.sh summary    # conteos por fuente / completitud / familia
#
# Versionado: cada descarga va a $VIRALDB_ROOT/<AAAAMMDD>/ y el enlace
# $VIRALDB_ROOT/current apunta a la última. Los pasos posteriores usan 'current'
# salvo que exportes VIRALDB_VERSION=AAAAMMDD. Nada se sobrescribe entre versiones.
#
# Entorno conda (una vez):
#   module load conda
#   conda create -n viraldb -c conda-forge -c bioconda \
#         ncbi-datasets-cli seqkit blast cd-hit python=3.11 unzip
#   DIAMOND y bowtie2 salen de $ENV_DIAMOND y $ENV_BOWTIE2 (config.sh), para que
#   el .dmnd y el índice los lean las MISMAS versiones que usan tus workers.
#
# Envíos LSF (ejecutar desde la carpeta de scripts, con $PWD, no ./):
#   bsub -J vdbdl -o "$LOGS_DIR/vdbdl.%J.log" -e "$LOGS_DIR/vdbdl.%J.err" \
#        -q normal -n 2 -M 4000 -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
#        "$PWD/build_viral_db.sh download"
#   bsub -J vdbprep -o "$LOGS_DIR/vdbprep.%J.log" -e "$LOGS_DIR/vdbprep.%J.err" \
#        -q normal -n 4 -M 16000 -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "$PWD/build_viral_db.sh prepare"
#   (blastn y diamond: -n 8 -M 16000.  mapping: -n 16 -M 32000, cola long si el
#    CD-HIT pasa de 12 h. Son estimaciones: MEDIR y anotar en bsubs_reproducibilidad.md)
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
TOOLS="$SCRIPT_DIR/viral_db_tools.py"

THREADS="${LSB_DJOB_NUMPROC:-4}"
VIRALDB_ROOT="${VIRALDB_ROOT:-$REFS_DIR/viral_custom}"
ENV_VIRALDB="${ENV_VIRALDB:-viraldb}"

# --- Qué se descarga --------------------------------------------------------
HOST_TAXON="${HOST_TAXON:-Arthropoda}"   # --host de Datasets (nombre o taxid; 6656)
COMPLETE_ONLY="${COMPLETE_ONLY:-0}"      # 1 = solo genomas completos del set de host
INCLUDE="${INCLUDE:-genome,protein,annotation}"

# --- Filtros -----------------------------------------------------------------
MIN_LEN_NT="${MIN_LEN_NT:-200}"          # blastn: fuera fragmentos tipo cebador
MIN_LEN_MAP="${MIN_LEN_MAP:-1000}"       # mapeo
MAX_LEN_MAP="${MAX_LEN_MAP:-500000}"     # mapeo: fuera virus gigantes (CD-HIT MAX_SEQ)
CDHIT_C="${CDHIT_C:-0.95}"               # mismos umbrales que los vOTUs (95/85)
CDHIT_AS="${CDHIT_AS:-0.85}"
PROT_A2T="${PROT_A2T:-}"                 # opcional: prot.accession2taxid.FULL.gz para rellenar huecos

NCBI_FTP="https://ftp.ncbi.nlm.nih.gov"
step="${1:-}"

ok()  { printf '  \033[32m[ok   ]\033[0m %s\n' "$*"; }
bad() { printf '  \033[31m[FALLO]\033[0m %s\n' "$*"; }
hdr() { echo; echo "==========================================================="; \
        echo "$*"; echo "==========================================================="; }
die() { bad "$*"; exit 1; }
nseq() { grep -c '^>' "$1" 2>/dev/null || echo 0; }

# Directorio de la versión: 'download' la crea; el resto usa la existente.
resolve_version() {
  if [[ -n "${VIRALDB_VERSION:-}" ]]; then
    VDIR="$VIRALDB_ROOT/$VIRALDB_VERSION"
  elif [[ "$step" == "download" || "$step" == "probe" ]]; then
    VDIR="$VIRALDB_ROOT/$(date +%Y%m%d)"
  else
    [[ -L "$VIRALDB_ROOT/current" ]] || die "no existe $VIRALDB_ROOT/current: corre 'download' primero"
    VDIR="$(readlink -f "$VIRALDB_ROOT/current")"
  fi
  RAW="$VDIR/raw"; WORK="$VDIR/work"; TAX="$VDIR/taxonomy"
  echo "  version: $VDIR"
}

use_env() { [[ -n "$1" ]] || die "entorno vacio"; activate_env "$1" >/dev/null || die "no se pudo activar $1"; }

datasets_key() { [[ -n "${NCBI_API_KEY:-}" ]] && echo "--api-key $NCBI_API_KEY"; }

# Descarga con reintentos y verificación del zip. Escribe a .tmp y renombra:
# un job que muere a media no deja un zip "bueno" a medias.
dl_datasets() {
  local name="$1"; shift
  local zip="$RAW/$name.zip"
  if [[ -s "$zip" ]] && unzip -tq "$zip" >/dev/null 2>&1; then ok "ya esta: $zip"; return 0; fi
  local try wait=60
  for try in 1 2 3 4; do
    echo "  [$name] intento $try: datasets download virus genome $*"
    rm -f "$zip.tmp"
    # shellcheck disable=SC2046
    if datasets download virus genome "$@" --include "$INCLUDE" --no-progressbar \
         $(datasets_key) --filename "$zip.tmp" \
       && unzip -tq "$zip.tmp" >/dev/null 2>&1; then
      mv "$zip.tmp" "$zip"
      ok "$zip  ($(du -h "$zip" | cut -f1))"
      return 0
    fi
    bad "intento $try fallido; espero ${wait}s"; sleep "$wait"; wait=$((wait * 2))
  done
  return 1
}

dl_ftp_md5() {   # fichero del FTP de NCBI + su .md5
  local url="$1" f; f="$(basename "$url")"
  if [[ -s "$f" && -s "$f.md5" ]] && md5sum -c "$f.md5" >/dev/null 2>&1; then ok "ya esta: $f"; return 0; fi
  wget -q -c -O "$f" "$url" && wget -q -O "$f.md5" "$url.md5" || return 1
  md5sum -c "$f.md5" >/dev/null 2>&1 || { rm -f "$f"; return 1; }
  ok "$f (md5 ok)"
}

# Extrae del zip solo los ficheros que interesan, sin reventar inodes.
unpack() {
  local name="$1"; local out="$RAW/$name"
  mkdir -p "$out"
  unzip -oq -j "$RAW/$name.zip" 'ncbi_dataset/data/*' -d "$out" || die "no se pudo descomprimir $name.zip"
  for f in genomic.fna data_report.jsonl; do [[ -s "$out/$f" ]] || die "$name: falta $f"; done
  ok "$name: $(nseq "$out/genomic.fna") genomas, $(nseq "$out/protein.faa") proteinas"
}

# ---------------------------------------------------------------------------
case "$step" in

check)
  hdr "1. Herramientas"
  use_env "$ENV_VIRALDB"
  for t in datasets seqkit makeblastdb blastdbcmd blastn cd-hit-est python3 unzip wget md5sum; do
    command -v "$t" >/dev/null && ok "$t" || bad "$t no esta en $ENV_VIRALDB"
  done
  datasets --version 2>/dev/null
  [[ -s "$TOOLS" ]] && ok "viral_db_tools.py" || bad "falta $TOOLS (va junto a este script)"
  use_env "$ENV_DIAMOND";  command -v diamond >/dev/null && ok "$(diamond version)" || bad "diamond"
  use_env "$ENV_BOWTIE2";  command -v bowtie2-build >/dev/null && ok "bowtie2-build" || bad "bowtie2-build"

  hdr "2. Red hacia NCBI (este nodo)"
  for u in https://api.ncbi.nlm.nih.gov/datasets/v2/version "$NCBI_FTP/pub/taxonomy/"; do
    curl -s -o /dev/null -m 20 -w '%{http_code}' "$u" | grep -q '^2' && ok "$u" || bad "$u sin acceso"
  done
  [[ -n "${NCBI_API_KEY:-}" ]] && ok "NCBI_API_KEY definida" \
    || echo "  (sin NCBI_API_KEY: funciona, pero con limite de peticiones mas bajo)"

  hdr "3. Cuota"
  lfs quota -gh team222 /lustre/scratch126 2>/dev/null || echo "  (lfs quota no disponible)"
  echo "  Estimacion a CONFIRMAR con 'summary' tras la descarga: pocas decenas de GB"
  echo "  en total (zips + fasta + indices). Se pueden borrar raw/*/ tras 'prepare'."
  ;;

probe)
  hdr "Probe: un genoma pequeño para ver el formato real de Datasets"
  resolve_version
  use_env "$ENV_VIRALDB"
  P="$VIRALDB_ROOT/_probe"; rm -rf "$P"; mkdir -p "$P"; cd "$P" || exit 1
  # NC_001477 = dengue 1 RefSeq: pequeño, anotado, con proteínas maduras
  datasets download virus genome accession NC_001477.1 --include "$INCLUDE" \
    --no-progressbar $(datasets_key) --filename probe.zip || die "fallo la descarga de prueba"
  unzip -oq probe.zip
  python3 "$TOOLS" probe "$P"
  echo
  echo "  Si todo dice OK y 'virus.taxId' tiene un numero, lanza 'download'."
  echo "  Si algo dice NO ESTA / NO ENCONTRADO, pega esta salida antes de seguir."
  ;;

download)
  hdr "Descarga (RED)"
  resolve_version
  mkdir -p "$RAW" "$TAX"
  use_env "$ENV_VIRALDB"
  datasets --version

  CO=(); [[ "$COMPLETE_ONLY" == 1 ]] && CO=(--complete-only)
  dl_datasets refseq taxon Viruses --refseq                       || die "fallo refseq"
  dl_datasets host   taxon Viruses --host "$HOST_TAXON" "${CO[@]}" || die "fallo host=$HOST_TAXON"

  cd "$TAX" || exit 1
  # taxdump FRESCO: los taxids de Datasets son de hoy; un taxdump viejo (el del
  # nr de feb-2024) no conoce los binomiales ICTV nuevos.
  dl_ftp_md5 "$NCBI_FTP/pub/taxonomy/taxdump.tar.gz" || die "fallo taxdump"
  tar -xzf taxdump.tar.gz nodes.dmp names.dmp merged.dmp
  dl_ftp_md5 "$NCBI_FTP/blast/db/taxdb.tar.gz"        || die "fallo taxdb"

  ln -sfn "$(basename "$VDIR")" "$VIRALDB_ROOT/current"
  {
    echo "fecha_descarga: $(date -Iseconds)"
    echo "datasets: $(datasets --version 2>&1)"
    echo "refseq: datasets download virus genome taxon Viruses --refseq --include $INCLUDE"
    echo "host:   datasets download virus genome taxon Viruses --host $HOST_TAXON ${CO[*]} --include $INCLUDE"
    echo "taxdump/taxdb: $NCBI_FTP (md5 verificado)"
  } > "$VDIR/MANIFEST.txt"
  ok "current -> $VDIR"
  ;;

prepare)
  hdr "Preparacion: metadatos, dedup, mapas de taxid"
  resolve_version
  use_env "$ENV_VIRALDB"
  mkdir -p "$WORK"; cd "$WORK" || exit 1
  unpack refseq
  unpack host

  echo "== metadatos (RefSeq primero: gana en duplicados) =="
  python3 "$TOOLS" meta --taxdump "$TAX" --out metadata.tsv.tmp \
      "refseq=$RAW/refseq/data_report.jsonl" "host=$RAW/host/data_report.jsonl" \
    && mv metadata.tsv.tmp metadata.tsv || die "fallo meta"

  echo "== nucleotidos =="
  # rmdup -n: misma accesion en los dos sets.  rmdup -s: misma secuencia con otra
  # accesion (las NC_ son copias de registros GenBank). Se queda la primera = RefSeq.
  cat "$RAW/refseq/genomic.fna" "$RAW/host/genomic.fna" \
    | seqkit seq -m "$MIN_LEN_NT" -u -g \
    | seqkit rmdup -n 2>/dev/null \
    | seqkit rmdup -s -D dup_nt.tsv -o viral_nt.fa.tmp 2> rmdup_nt.log \
    && mv viral_nt.fa.tmp viral_nt.fa || die "fallo dedup nt"
  cat rmdup_nt.log
  ok "viral_nt.fa: $(nseq viral_nt.fa) secuencias"
  python3 "$TOOLS" ntmap --meta metadata.tsv --fasta viral_nt.fa --out viral_nt.taxidmap || die "fallo ntmap"

  echo "== proteinas =="
  FAA=(); ANN=()
  for s in refseq host; do
    [[ -s "$RAW/$s/protein.faa" ]] && FAA+=("$RAW/$s/protein.faa")
    [[ -s "$RAW/$s/annotation_report.jsonl" ]] && ANN+=("$RAW/$s/annotation_report.jsonl")
  done
  (( ${#FAA[@]} )) || die "no hay protein.faa: ¿INCLUDE sin 'protein'?"
  cat "${FAA[@]}" \
    | seqkit rmdup -n 2>/dev/null \
    | seqkit rmdup -s -D dup_prot.tsv -o viral_prot.faa.tmp 2> rmdup_prot.log \
    && mv viral_prot.faa.tmp viral_prot.faa || die "fallo dedup prot"
  cat rmdup_prot.log
  ok "viral_prot.faa: $(nseq viral_prot.faa) proteinas"
  python3 "$TOOLS" protmap --meta metadata.tsv --faa viral_prot.faa --annot "${ANN[@]}" \
      --out prot2taxid.gz --missing prot_sin_taxid.txt || die "fallo protmap"

  if [[ -s prot_sin_taxid.txt && -n "$PROT_A2T" && -s "$PROT_A2T" ]]; then
    echo "== rellenando $(wc -l < prot_sin_taxid.txt) proteinas con $PROT_A2T =="
    zcat "$PROT_A2T" | awk -F'\t' 'NR==FNR{m[$1]=1; next} ($1 in m){print $1"\t"$2}' \
        prot_sin_taxid.txt - > extra.tsv
    ok "rellenadas: $(wc -l < extra.tsv)"
    { zcat prot2taxid.gz; cat extra.tsv; } | gzip > prot2taxid.gz.tmp && mv prot2taxid.gz.tmp prot2taxid.gz
  fi

  python3 "$TOOLS" summary --meta metadata.tsv --fasta viral_nt.fa > summary_nt.txt
  head -20 summary_nt.txt
  echo "prepare: $(date -Iseconds) | nt=$(nseq viral_nt.fa) prot=$(nseq viral_prot.faa)" >> "$VDIR/MANIFEST.txt"
  echo
  echo "  Si todo cuadra, puedes liberar espacio:  rm -r $RAW/refseq $RAW/host"
  echo "  (los zips se quedan: son la copia fiel de lo descargado)"
  ;;

blastn)
  hdr "BLAST nucleotidos"
  resolve_version
  use_env "$ENV_VIRALDB"
  [[ -s "$WORK/viral_nt.fa" ]] || die "falta viral_nt.fa: corre 'prepare'"
  B="$VDIR/blast"; mkdir -p "$B"; cd "$B" || exit 1
  tar -xzf "$TAX/taxdb.tar.gz" -C "$B"
  makeblastdb -in "$WORK/viral_nt.fa" -dbtype nucl -parse_seqids \
      -taxid_map "$WORK/viral_nt.taxidmap" -blastdb_version 5 \
      -title "viral_nt RefSeq+host:$HOST_TAXON $(basename "$VDIR")" \
      -out "$B/viral_nt" || die "fallo makeblastdb"
  export BLASTDB="$B"
  blastdbcmd -db viral_nt -info | head -8

  echo "== autocomprobacion: 3 secuencias contra la base, deben encontrarse a si mismas con taxid =="
  # las 3 primeras (RefSeq) y las 2 ultimas (set de huesped)
  { seqkit range -r 1:3 "$WORK/viral_nt.fa"; seqkit range -r -2:-1 "$WORK/viral_nt.fa"; } > selftest.fa
  blastn -task megablast -query selftest.fa -db viral_nt -max_target_seqs 1 -num_threads "$THREADS" \
         -outfmt "6 qacc sacc pident length staxids sscinames" | sort -u -k1,1 | tee selftest.tsv
  # sacc sale sin version (NC_001477), qacc con ella (NC_001477.1)
  awk '{q=$1; sub(/\.[0-9]+$/,"",q)} q!=$2 || $5=="" || $5=="0"{e=1} END{exit (e || NR<5)}' selftest.tsv \
    && ok "autocomprobacion" || bad "la autocomprobacion no cuadra (¿taxid vacio o hit ajeno?)"
  echo
  echo "  Uso:  export BLASTDB=$B"
  echo "        blastn -task megablast -db viral_nt -query vOTUs.fa -max_target_seqs 50 \\"
  echo "          -outfmt '6 std qlen slen staxids sscinames' -num_threads 8"
  ;;

diamond)
  hdr "DIAMOND proteinas"
  resolve_version
  [[ -s "$WORK/viral_prot.faa" ]] || die "falta viral_prot.faa: corre 'prepare'"
  D="$VDIR/diamond"; mkdir -p "$D"
  use_env "$ENV_DIAMOND"; diamond version
  # TRAMPA: desde 2025 el taxdump de NCBI usa rangos nuevos (realm, domain,
  # acellular root...) que DIAMOND anteriores a esos cambios no conocen:
  # "Error: Invalid taxonomic rank: realm". Se reintenta con una copia de
  # nodes.dmp donde el rango que DIAMOND rechaza pasa a "no rank". Solo cambia la
  # etiqueta del rango, no el arbol: taxids, LCA y nombres siguen igual.
  NODES="$TAX/nodes.dmp"; cp "$NODES" "$D/nodes_diamond.dmp"; built=0; prev=""
  for i in $(seq 1 10); do
    diamond makedb --in "$WORK/viral_prot.faa" -d "$D/viral_prot" \
        --taxonmap "$WORK/prot2taxid.gz" --taxonnodes "$NODES" --taxonnames "$TAX/names.dmp" \
        --threads "$THREADS" 2> "$D/makedb.err" && { built=1; break; }
    r=$(sed -n 's/.*Invalid taxonomic rank: //p' "$D/makedb.err" | head -1)
    [[ -n "$r" && "$r" != "$prev" ]] || { cat "$D/makedb.err"; die "fallo makedb"; }
    prev="$r"
    echo "  DIAMOND no conoce el rango '$r' -> 'no rank' en nodes_diamond.dmp"
    awk -v r="$r" 'BEGIN{FS="\t[|]\t"; OFS="\t|\t"} $3==r{$3="no rank"} {print}' "$D/nodes_diamond.dmp" > "$D/nodes.tmp" \
      && mv "$D/nodes.tmp" "$D/nodes_diamond.dmp"
    NODES="$D/nodes_diamond.dmp"
  done
  (( built )) || { cat "$D/makedb.err"; die "fallo makedb"; }
  diamond dbinfo -d "$D/viral_prot.dmnd"

  echo "== autocomprobacion =="
  awk '/^>/{n++} n>3{exit} {print}' "$WORK/viral_prot.faa" > "$D/selftest.faa"
  diamond blastp -q "$D/selftest.faa" -d "$D/viral_prot.dmnd" -k 1 --quiet --threads "$THREADS" \
      -f 6 qseqid sseqid pident staxids sscinames | tee "$D/selftest.tsv"
  awk -F'\t' '$1!=$2 || $4=="" || $4=="0"{e=1} END{exit (e || NR<3)}' "$D/selftest.tsv" \
    && ok "autocomprobacion" || bad "la autocomprobacion no cuadra"
  echo
  echo "  Recuerda (build_downstream_dbs.sh): una base SOLO viral asigna a todo contig"
  echo "  su mejor hit viral. Sirve para DESCUBRIR; para DESCARTAR sigue haciendo falta nr."
  ;;

mapping)
  hdr "Referencia de mapeo: completos, CD-HIT-EST ${CDHIT_C}/${CDHIT_AS}, bowtie2"
  resolve_version
  use_env "$ENV_VIRALDB"
  [[ -s "$WORK/viral_nt.fa" ]] || die "falta viral_nt.fa: corre 'prepare'"
  M="$VDIR/mapping"; mkdir -p "$M"; cd "$M" || exit 1
  python3 "$TOOLS" select --meta "$WORK/metadata.tsv" --fasta "$WORK/viral_nt.fa" --out ids_map.txt
  seqkit grep -f ids_map.txt "$WORK/viral_nt.fa" \
    | seqkit seq -g -m "$MIN_LEN_MAP" -M "$MAX_LEN_MAP" > input_map.fa
  ok "entrada: $(nseq input_map.fa) secuencias"

  if [[ ! -s viral_nr95.fa ]]; then
    # -G 0 + -aS: cobertura sobre la secuencia corta (equivalente a los vOTUs).
    # -n 10 exige -c >= 0.95.  -M 0 = sin tope de RAM interno (el tope lo pone LSF).
    cd-hit-est -i input_map.fa -o viral_nr95.fa.tmp -c "$CDHIT_C" -aS "$CDHIT_AS" \
               -G 0 -n 10 -d 0 -M 0 -T "$THREADS" > cdhit.log 2>&1 || { tail cdhit.log; die "fallo cd-hit-est"; }
    mv viral_nr95.fa.tmp.clstr viral_nr95.fa.clstr; mv viral_nr95.fa.tmp viral_nr95.fa
  fi
  python3 "$TOOLS" clstr viral_nr95.fa.clstr --out viral_nr95_clusters.tsv
  ok "representantes: $(nseq viral_nr95.fa)"
  seqkit fx2tab -n -i -l viral_nr95.fa > viral_nr95_lengths.tsv

  use_env "$ENV_BOWTIE2"
  bowtie2-build --threads "$THREADS" viral_nr95.fa viral_nr95 > bt2build.log 2>&1 \
    || { tail bt2build.log; die "fallo bowtie2-build"; }
  ok "indice: $M/viral_nr95.*.bt2"
  echo
  echo "  viral_nr95_clusters.tsv dice que accesiones colapsan en cada representante:"
  echo "  una lectura que cae en un representante 'es' de cualquiera de sus miembros."
  ;;

summary)
  resolve_version
  use_env "$ENV_VIRALDB"
  cat "$VDIR/MANIFEST.txt"; echo
  python3 "$TOOLS" summary --meta "$WORK/metadata.tsv" --fasta "$WORK/viral_nt.fa"
  echo; du -sh "$VDIR"/* 2>/dev/null
  ;;

*)
  sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
  ;;
esac
