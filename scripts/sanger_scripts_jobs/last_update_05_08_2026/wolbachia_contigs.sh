#!/usr/bin/env bash
# =============================================================================
# wolbachia_contigs.sh -- de donde viene la Wolbachia de M20/M21/M35/M60
#
# Kraken2 dice "Wolbachia pipientis" (taxid 955), que es el NODO DE ESPECIE:
# debajo de el cuelgan wMel, wAlbA, wAlbB y wPip. Asi que el reporte no
# tipifica la cepa y no distingue las cuatro explicaciones posibles
# (Ae. albopictus simpatrico, colonia wAlbB del LCB-UADY, Culex con wPip,
# o algo divergente).
#
# Esto lo resuelve con los contigs que ya estan ensamblados: blastn del union
# completo contra todos los genomas completos de Wolbachia en RefSeq, cada uno
# etiquetado con su cepa. El mejor hit tipifica; la amplitud de cobertura dice
# si el genoma esta de verdad o si son hits en regiones conservadas.
#
# Las 18 muestras sin senal entran en la misma corrida y sirven de control
# interno: si tambien dan contigs de Wolbachia, hay arrastre generalizado.
#
# Subcomandos:
#   check   verifica herramientas, rutas y red. No escribe nada.
#   refs    baja los genomas de Wolbachia y construye la base de blast.
#   run     blastn del union y resumen por muestra y por cepa.
#
#   bash wolbachia_contigs.sh check
#   bash wolbachia_contigs.sh refs
#   bash wolbachia_contigs.sh run
#
# En la farm conviene mandar 'run' a la cola:
#   bsub -q normal -n 4 -M 8000 -R "select[mem>8000] rusage[mem=8000]" \
#        -o logs/wolb.%J.out -e logs/wolb.%J.err \
#        "bash wolbachia_contigs.sh run"
# =============================================================================
set -euo pipefail

# ----------------------------------------------------------------- parametros
RAMA="${RAMA:-min500}"          # min500 (catalogo principal) o 1kb
IDENT_MIN="${IDENT_MIN:-90}"    # identidad minima para contar un hit
LEN_MIN="${LEN_MIN:-200}"       # longitud minima de alineamiento
CPUS="${CPUS:-4}"

BLAST_BIN=/software/treeoflife/conda/users/envs/team222/jm79/blast_2.17.0/bin
SUMMARY_URL=https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

# ------------------------------------------------------------------ contexto
# config.sh define REFS_DIR / RESULTS_DIR / LOGS_DIR, como el resto del pipeline
for c in ./config.sh ../config.sh "$HOME/jkvirome/jm79/mosquito-virome-jock/scripts/sanger_scripts_jobs/last_update_05_08_2026/config.sh"; do
  if [[ -f "$c" ]]; then source "$c"; break; fi
done
: "${RESULTS_DIR:=/lustre/scratch126/tol/teams/lawniczak/users/jm79/results/MERIDA}"
: "${REFS_DIR:=/lustre/scratch126/tol/teams/lawniczak/users/jm79/refs}"
: "${LOGS_DIR:=$RESULTS_DIR/../../logs}"

if [[ "$RAMA" == "min500" ]]; then
  UNION="$RESULTS_DIR/contig_annot_min500/contigs_union_MERIDA.fasta"
else
  UNION="$RESULTS_DIR/contig_annot/contigs_union_MERIDA.fasta"
fi

WDIR="$REFS_DIR/wolbachia"
GEN="$WDIR/genomas"
DB="$WDIR/wolbachia_completos"
REFTAB="$WDIR/refs.tsv"
OUT="$RESULTS_DIR/wolbachia_$RAMA"

export PATH="$BLAST_BIN:$PATH"

msg() { printf '  %s\n' "$*"; }

# -------------------------------------------------------------------- check
cmd_check() {
  echo "== herramientas =="
  local falta=0
  for t in blastn makeblastdb curl python3 awk; do
    if command -v "$t" >/dev/null 2>&1; then
      msg "ok    $t  ($(command -v "$t"))"
    else
      msg "FALTA $t"; falta=1
    fi
  done

  echo
  echo "== rutas =="
  if [[ -s "$UNION" ]]; then
    msg "ok    union ($RAMA): $UNION"
    msg "      contigs: $(grep -c '^>' "$UNION")"
    msg "      ejemplo de ID: $(grep -m1 '^>' "$UNION" | cut -c2-)"
  else
    msg "FALTA union: $UNION"; falta=1
  fi
  for d in "$REFS_DIR" "$RESULTS_DIR"; do
    [[ -d "$d" ]] && msg "ok    $d" || { msg "FALTA $d"; falta=1; }
  done

  echo
  echo "== red =="
  local code
  code=$(curl -s -o /dev/null -w '%{http_code}' --max-time 30 "$SUMMARY_URL" || echo 000)
  if [[ "$code" == "200" ]]; then
    msg "ok    RefSeq assembly_summary alcanzable"
  else
    msg "FALLA RefSeq devuelve $code -- sin red no se pueden bajar las referencias"; falta=1
  fi

  echo
  echo "== estado =="
  [[ -s "$REFTAB" ]] && msg "referencias ya preparadas: $(wc -l < "$REFTAB") secuencias" \
                     || msg "referencias sin preparar -- corre: refs"
  [[ -s "$OUT/hits.tsv" ]] && msg "ya hay hits en $OUT/hits.tsv" || true

  echo
  [[ $falta -eq 0 ]] && echo "check: todo listo" || { echo "check: faltan cosas"; return 1; }
}

# --------------------------------------------------------------------- refs
cmd_refs() {
  mkdir -p "$GEN"
  local sum="$WDIR/assembly_summary_wolbachia.tsv"

  if [[ ! -s "$sum" ]]; then
    msg "filtrando assembly_summary de RefSeq (el fichero es grande, va en streaming)"
    # 8 organism_name, 9 infraspecific_name, 12 assembly_level, 20 ftp_path
    curl -sS --max-time 1800 "$SUMMARY_URL" \
      | awk -F'\t' 'BEGIN{OFS="\t"}
          $8 ~ /Wolbachia/ && ($12=="Complete Genome" || $12=="Chromosome") && $20!="na" {
            print $1, $8, $9, $12, $20
          }' > "$sum"
  fi
  msg "ensamblados de Wolbachia completos: $(wc -l < "$sum")"
  [[ -s "$sum" ]] || { echo "ERROR: el filtro no devolvio nada"; return 1; }

  : > "$REFTAB"
  local n=0
  while IFS=$'\t' read -r acc org infra nivel ftp; do
    local base fna cepa
    base=$(basename "$ftp")
    fna="$GEN/${acc}.fna"
    # cepa: de infraspecific_name (strain=wAlbB) o del ultimo token del organismo
    cepa=$(printf '%s' "$infra" | sed 's/^strain=//')
    [[ -z "$cepa" || "$cepa" == "na" ]] && cepa=$(printf '%s' "$org" | awk '{print $NF}')
    cepa=$(printf '%s' "$cepa" | tr -c 'A-Za-z0-9._-' '_')

    if [[ ! -s "$fna" ]]; then
      curl -sS --max-time 600 "$ftp/${base}_genomic.fna.gz" | gunzip -c > "$fna" || {
        msg "aviso: no pude bajar $acc, lo salto"; rm -f "$fna"; continue; }
    fi
    # accesion de secuencia -> cepa -> longitud (cromosoma y plasmidos aparte)
    awk -v c="$cepa" 'BEGIN{OFS="\t"}
        /^>/ { if (id!="") print id, c, n; id=substr($1,2); n=0; next }
             { n += length($0) }
        END  { if (id!="") print id, c, n }' "$fna" >> "$REFTAB"
    n=$((n+1))
  done < "$sum"

  msg "genomas descargados: $n"
  msg "secuencias en la tabla: $(wc -l < "$REFTAB")"
  msg "cepas: $(cut -f2 "$REFTAB" | sort -u | tr '\n' ' ')"

  cat "$GEN"/*.fna > "$WDIR/wolbachia_completos.fasta"
  makeblastdb -in "$WDIR/wolbachia_completos.fasta" -dbtype nucl -out "$DB" >/dev/null
  msg "base de blast lista: $DB"
}

# ---------------------------------------------------------------------- run
cmd_run() {
  [[ -s "$REFTAB" ]] || { echo "ERROR: faltan las referencias, corre 'refs' antes"; return 1; }
  mkdir -p "$OUT"

  if [[ ! -s "$OUT/hits.tsv" ]]; then
    msg "blastn del union contra Wolbachia (dc-megablast, sensible a ~80-90 % id)"
    blastn -task dc-megablast \
           -query "$UNION" -db "$DB" \
           -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
           -evalue 1e-20 -max_target_seqs 50 -num_threads "$CPUS" \
           -out "$OUT/hits.tsv"
  else
    msg "reutilizando $OUT/hits.tsv (borralo para repetir el blast)"
  fi
  msg "lineas de hit: $(wc -l < "$OUT/hits.tsv")"

  python3 -I "$(dirname "$0")/wolb_resumen.py" \
      "$OUT/hits.tsv" "$REFTAB" "$IDENT_MIN" "$LEN_MIN" \
      | tee "$OUT/resumen.txt"

  echo
  msg "escrito: $OUT/hits.tsv y $OUT/resumen.txt"
}

case "${1:-}" in
  check) cmd_check ;;
  refs)  cmd_refs ;;
  run)   cmd_run ;;
  *) sed -n '2,33p' "$0"; exit 1 ;;
esac
