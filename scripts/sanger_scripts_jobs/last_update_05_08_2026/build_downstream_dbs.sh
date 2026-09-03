#!/bin/bash
# ---------------------------------------------------------------------------
# build_downstream_dbs.sh — bases de datos de geNomad, CheckV y DIAMOND.
#
# Los cuatro pasos posteriores al ensamblaje (quast, genomad, checkv, diamond)
# estan bloqueados por estas bases. Son horas de descarga que NO dependen de que
# los ensamblajes terminen: se lanzan en paralelo.
#
# Uso:
#   ./build_downstream_dbs.sh check      # cuota, copias compartidas. HAZ ESTO 1o
#   ./build_downstream_dbs.sh genomad    # ~3 GB
#   ./build_downstream_dbs.sh checkv     # ~7 GB
#   ./build_downstream_dbs.sh rvdb       # ~5 GB   base viral curada
#   ./build_downstream_dbs.sh nr         # ~250 GB  ¡lee el aviso de 'check'!
#
# CADA PASO ES UN ENVIO APARTE. Nada de 'all': si nr falla a las seis horas no
# quieres arrastrar a los demas.
#
# POR QUE nr Y rvdb, Y NO SOLO UNA
#   Una base solo viral asigna a TODO contig su mejor hit viral, aunque el contig
#   sea 28S de mosquito: no tiene con que competir. nr permite DESCARTAR, que en
#   este proyecto -donde el rRNA es el problema- importa mas que acertar.
#   rvdb da resultados pronto; nr los confirma.
#
# Envio (necesitan RED; van bien en la cola normal):
#   bsub -J dbgenomad -o "$LOGS_DIR/dbgenomad.%J.log" -e "$LOGS_DIR/dbgenomad.%J.err" \
#        -q normal -n 4 -M 8000 \
#        -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "$SCRIPTS_DIR/build_downstream_dbs.sh genomad"
#
#   # nr necesita mas de todo:
#   bsub -J dbnr -o "$LOGS_DIR/dbnr.%J.log" -e "$LOGS_DIR/dbnr.%J.err" \
#        -q long -n 8 -M 64000 \
#        -R "select[mem>64000] rusage[mem=64000] span[hosts=1]" \
#        "$SCRIPTS_DIR/build_downstream_dbs.sh nr"
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

THREADS="${LSB_DJOB_NUMPROC:-4}"
DB_ROOT="${DB_ROOT:-$REFS_DIR}"
step="${1:-}"

ok()  { printf '  \033[32m[ok   ]\033[0m %s\n' "$*"; }
bad() { printf '  \033[31m[FALLO]\033[0m %s\n' "$*"; }
hdr() { echo; echo "==========================================================="; \
        echo "$*"; echo "==========================================================="; }

# ---------------------------------------------------------------------------
case "$step" in

check)
  hdr "1. Cuota de disco e inodes"
  lfs quota -h -u "$USER" "$SCRATCH" 2>/dev/null || echo "  (lfs quota no disponible)"
  echo
  echo "  nr ocupa ~250 GB descargado + ~200 GB indexado, y durante el makedb"
  echo "  conviven los dos. Si no tienes ~500 GB libres, NO lances 'nr':"
  echo "  pide ampliacion de cuota o usa una copia compartida."

  hdr "2. Copias compartidas en la farm (antes de descargar nada)"
  # OJO: NO recorrer /lustre entero. Es un sistema de ficheros de petabytes y un
  # find a profundidad 6 tarda horas. Solo se miran rutas plausibles, poco
  # profundas y con timeout.
  cands=(
    /lustre/scratch*/tol/teams/*/refs
    /lustre/scratch*/blastdb
    /data/blastdb /data/blastdbs
    /nfs/*/blastdb
    /software/db /software/blastdb
    "$REFS_DIR"
  )
  found=0
  for d in "${cands[@]}"; do
    [[ -d "$d" ]] || continue
    while IFS= read -r hit; do
      echo "    $hit"; found=1
    done < <(timeout 30 find "$d" -maxdepth 3 \
               \( -name "nr.dmnd" -o -name "genomad_db" -o -name "checkv-db-*" \) 2>/dev/null | head -5)
  done
  (( found )) || echo "    (nada en las rutas habituales)"
  echo
  echo "  Esta busqueda es deliberadamente corta. Si quieres saber de verdad si"
  echo "  el equipo ya tiene nr montado, pregunta en el helpdesk de Sanger o mira:"
  echo "    module avail 2>&1 | grep -i -E 'blast|diamond'"
  echo "    ls -d /lustre/scratch*/tol/teams/lawniczak/*/refs* 2>/dev/null"
  echo
  echo "  Si aparece un nr.dmnd compartido, comprueba la version de formato:"
  echo "    diamond dbinfo -d <ruta>.dmnd"
  echo "  DIAMOND rompe compatibilidad entre versiones mayores. Si no la lee,"
  echo "  la copia compartida no te sirve y hay que reconstruir."

  hdr "3. Entornos conda"
  for e in "${ENV_GENOMAD:-}" "${ENV_CHECKV:-}" "${ENV_DIAMOND:-}" "${ENV_QUAST:-}"; do
    [[ -z "$e" ]] && { bad "variable de entorno vacia en config.sh"; continue; }
    conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$e" \
      && ok "$e" || bad "$e no existe"
  done
  ;;

genomad)
  hdr "geNomad — base de datos"
  out="$DB_ROOT/genomad"
  mkdir -p "$out"
  activate_env "$ENV_GENOMAD"
  command -v genomad >/dev/null || { bad "genomad no esta en el PATH"; exit 1; }
  genomad --version

  if [[ -d "$out/genomad_db" ]] && [[ -n "$(ls -A "$out/genomad_db" 2>/dev/null)" ]]; then
    ok "ya existe: $out/genomad_db"
  else
    # download-database crea el subdirectorio genomad_db/ dentro del destino
    genomad download-database "$out" || { bad "fallo la descarga"; exit 1; }
  fi

  n=$(find "$out/genomad_db" -type f | wc -l)
  sz=$(du -sh "$out/genomad_db" | cut -f1)
  (( n > 5 )) && ok "$out/genomad_db  ($n ficheros, $sz)" || { bad "solo $n ficheros"; exit 1; }
  echo
  echo "  En config.sh:  GENOMAD_DB=$out/genomad_db"
  ;;

checkv)
  hdr "CheckV — base de datos"
  out="$DB_ROOT/checkv"
  mkdir -p "$out"
  activate_env "$ENV_CHECKV"
  command -v checkv >/dev/null || { bad "checkv no esta en el PATH"; exit 1; }

  db=$(find "$out" -maxdepth 1 -type d -name "checkv-db-*" | sort | tail -1)
  if [[ -n "$db" && -d "$db/genome_db" ]]; then
    ok "ya existe: $db"
  else
    checkv download_database "$out" || { bad "fallo la descarga"; exit 1; }
    db=$(find "$out" -maxdepth 1 -type d -name "checkv-db-*" | sort | tail -1)
  fi
  [[ -n "$db" ]] || { bad "no aparece ningun checkv-db-*"; ls -la "$out"; exit 1; }

  # TRAMPA: el tarball trae un .dmnd preconstruido. Si tu DIAMOND es de otra
  # version mayor, no lo puede leer y CheckV falla a mitad, no al empezar.
  if [[ -s "$db/genome_db/checkv_reps.dmnd" ]]; then
    activate_env "$ENV_DIAMOND"
    if diamond dbinfo -d "$db/genome_db/checkv_reps.dmnd" >/dev/null 2>&1; then
      ok "el .dmnd de CheckV es compatible con tu DIAMOND"
    else
      echo "  el .dmnd viene de otra version de DIAMOND. Reconstruyendo..."
      diamond makedb --in "$db/genome_db/checkv_reps.faa" \
                     -d "$db/genome_db/checkv_reps" --threads "$THREADS" \
        && ok "reconstruido" || bad "no se pudo reconstruir"
    fi
  fi

  sz=$(du -sh "$db" | cut -f1)
  ok "$db  ($sz)"
  echo
  echo "  En config.sh:  CHECKV_DB=$db"
  ;;

rvdb)
  hdr "RVDB-prot — base viral curada para DIAMOND"
  out="$DB_ROOT/diamond"; mkdir -p "$out"
  # OJO: la version sube cada pocos meses. Comprueba la actual en
  #      https://rvdb-prot.pasteur.fr/  y pasala con RVDB_URL=...
  RVDB_URL="${RVDB_URL:-}"
  [[ -n "$RVDB_URL" ]] || {
    bad "define RVDB_URL con el fichero .fasta.xz actual de https://rvdb-prot.pasteur.fr/"
    echo "  ej: RVDB_URL=https://rvdb-prot.pasteur.fr/files/U-RVDBvXX.0-prot.fasta.xz"
    exit 1; }

  fa="$out/$(basename "${RVDB_URL%.xz}")"
  if [[ -s "$out/rvdb.dmnd" ]]; then ok "ya existe: $out/rvdb.dmnd"; exit 0; fi
  [[ -s "$fa" ]] || { wget -O "$fa.xz" "$RVDB_URL" && unxz "$fa.xz"; }
  [[ -s "$fa" ]] || { bad "no se descargo el fasta"; exit 1; }
  ok "fasta: $fa  ($(grep -c '^>' "$fa") proteinas)"

  activate_env "$ENV_DIAMOND"
  diamond makedb --in "$fa" -d "$out/rvdb" --threads "$THREADS" \
    || { bad "fallo makedb"; exit 1; }
  ok "$out/rvdb.dmnd  ($(du -sh "$out/rvdb.dmnd" | cut -f1))"
  echo
  echo "  En config.sh:  DIAMOND_DB_RVDB=$out/rvdb.dmnd"
  ;;

nr)
  hdr "NCBI nr — base completa para DIAMOND (LARGO)"
  out="$DB_ROOT/diamond"; mkdir -p "$out"
  cd "$out" || exit 1

  if [[ -s "$out/nr.dmnd" ]]; then ok "ya existe: $out/nr.dmnd"; exit 0; fi

  echo "== descargando nr.gz y la taxonomia =="
  # SIN --taxonmap el resultado son accesiones sin taxonomia, que para asignar
  # virus no vale de nada. Los tres ficheros de taxonomia NO son opcionales.
  for u in \
    https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz \
    https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz \
    https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz ; do
    f=$(basename "$u")
    [[ -s "$f" ]] && { ok "ya esta: $f"; continue; }
    wget -c -O "$f" "$u" || { bad "fallo la descarga de $f"; exit 1; }
  done
  [[ -s nodes.dmp ]] || tar -xzf taxdump.tar.gz nodes.dmp names.dmp

  echo "== diamond makedb (esto es lo que tarda) =="
  activate_env "$ENV_DIAMOND"
  diamond --version
  diamond makedb \
      --in nr.gz \
      -d "$out/nr" \
      --taxonmap prot.accession2taxid.FULL.gz \
      --taxonnodes nodes.dmp \
      --taxonnames names.dmp \
      --threads "$THREADS" \
    || { bad "fallo makedb"; exit 1; }

  ok "$out/nr.dmnd  ($(du -sh "$out/nr.dmnd" | cut -f1))"
  echo
  echo "  Ya puedes borrar nr.gz (~90 GB) si andas justo de cuota:"
  echo "    rm $out/nr.gz"
  echo "  En config.sh:  DIAMOND_DB=$out/nr.dmnd"
  ;;

*)
  echo "Uso: ./build_downstream_dbs.sh <check|genomad|checkv|rvdb|nr>"
  echo
  echo "  check     cuota, copias compartidas y entornos.  EMPIEZA POR AQUI."
  echo "  genomad   ~3 GB    clasificacion viral"
  echo "  checkv    ~7 GB    calidad y completitud de contigs virales"
  echo "  rvdb      ~5 GB    DIAMOND contra base viral curada"
  echo "  nr        ~250 GB  DIAMOND contra nr, la que permite DESCARTAR"
  exit 0
  ;;
esac