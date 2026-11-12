#!/bin/bash
# ---------------------------------------------------------------------------
# srna_build_refs.sh — indices bowtie1 para small RNA
#
#   ./srna_build_refs.sh check            # que hay y que falta (login)
#   ./srna_build_refs.sh genome           # AaegL5 completo          (~1-3 h)
#   ./srna_build_refs.sh rrna             # rDNA (refs/rrna/rrna_filtro.fasta)
#   ./srna_build_refs.sh viral            # vOTUs + extras (AEFE1, etc.)
#
# Envios:
#   source "$S/srna_env.sh"
#   bsub -J srnaref_genome -o "$LOGS_DIR/srnaref_genome.%J.log" -e "$LOGS_DIR/srnaref_genome.%J.err" \
#        -q long -n 8 -M 24000 -R "select[mem>24000] rusage[mem=24000] span[hosts=1]" \
#        "$S/srna_build_refs.sh genome"
#   for t in rrna viral; do
#     bsub -J srnaref_$t -o "$LOGS_DIR/srnaref_$t.%J.log" -e "$LOGS_DIR/srnaref_$t.%J.err" \
#          -q normal -n 2 -M 4000 -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
#          "$S/srna_build_refs.sh $t"
#   done
#   (recursos = estimaciones; medir con bhist -l y anotar en bsubs_reproducibilidad.md)
#
# POR QUE BOWTIE1 Y NO BOWTIE2
#   Para lecturas de 18-30 nt sin gaps, bowtie1 en modo -v (mismatches totales,
#   sin calidad) es el estandar de facto (Olmo 2023: v1.1.2 -v 0; Lau 2026: MSRG
#   sobre bowtie1). bowtie2 esta pensado para reads largos: su puntuacion por
#   semilla y su modo local no tienen sentido en 21 nt, y los multi-mapeos de
#   piRNA de TEs se controlan mejor con -a --best --strata / -M.
#
# POR QUE TODAS LAS vOTUs Y NO SOLO LAS VIRALES
#   Los candidatos a NIRVS de RESULTADO_votus_vs_nr.md (Flavi 9, Rhabdo 4,
#   Peribunya 11) estan entre los "celulares". El small RNA es justo la prueba
#   que los separa: un virus que replica da siRNA de 21 nt en las dos hebras;
#   un EVE da piRNA, sobre todo antisentido, y nada de 21 nt. Mapear a todas
#   de una vez cuesta lo mismo y los subconjuntos se filtran despues.
#
# DONDE ACABA EL FASTA DEL SUBCOMANDO `viral`
#   Junto a SU PROPIO indice: $(dirname "$BT1_VIRAL")/viral_targets.fa, con su
#   MANIFEST.txt al lado. Antes era una ruta fija en $SRNA_REFS/fasta/, asi que
#   construir el indice del catalogo de 500 bp PISABA el FASTA del de 1 kb y con
#   el su trazabilidad. Cambiando BT1_VIRAL se cambia de catalogo sin tocar nada
#   de lo anterior.
#
# Variables (exportar antes para cambiar):
#   GENOME_FA   FASTA de AaegL5. Si no existe, se reconstruye desde $BT2_INDEX
#               con bowtie2-inspect (el mismo truco que con los vOTUs).
#   RRNA_FA     por defecto $REFS_DIR/rrna/rrna_filtro.fasta
#   VOTUS_FA    por defecto $REFS_DIR/votu_index/votus_from_index.fasta
#               para el catalogo de 5,899 vOTUs (corte >=500 bp, el anclaje del
#               nivel B desde el 12 nov 2026):
#               $REFS_DIR/votu_catalog/MERIDA_cdhit95_85_G1g1_min500/votus.fa
#   VIRAL_EXTRA FASTA(s) extra separados por comas (AEFE1, PCLV, HTV, CFAV...)
#   BT1_VIRAL   prefijo del indice viral; decide tambien donde va el FASTA
# ---------------------------------------------------------------------------
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export S="${S:-$SCRIPT_DIR}"
source "$SCRIPT_DIR/srna_env.sh" "${SRNA_LOTE:-MERIDA_SRNA}" >/dev/null

cmd="${1:?uso: srna_build_refs.sh check|genome|rrna|viral}"
THREADS="${LSB_DJOB_NUMPROC:-2}"
GENOME_FA="${GENOME_FA:-$SRNA_REFS/fasta/AaegL5.fa}"
RRNA_FA="${RRNA_FA:-$REFS_DIR/rrna/rrna_filtro.fasta}"
VOTUS_FA="${VOTUS_FA:-$REFS_DIR/votu_index/votus_from_index.fasta}"
VIRAL_EXTRA="${VIRAL_EXTRA:-}"

idx_ok() { [[ -s "$1.1.ebwt" || -s "$1.1.ebwtl" ]] && [[ -s "$1.rev.1.ebwt" || -s "$1.rev.1.ebwtl" ]]; }

manifest() {   # manifest <prefijo> <fasta> <origenes...>
  local p="$1" fa="$2"; shift 2
  { echo "fecha: $(date -Iseconds)"; echo "bowtie: $(bowtie --version | head -1)"
    echo "fasta: $fa  ($(grep -c '^>' "$fa") secuencias, $(md5sum < "$fa" | cut -c1-12))"
    for o in "$@"; do echo "origen: $o"; done; } > "$(dirname "$p")/MANIFEST.txt"
}

build() {      # build <fasta> <prefijo>
  local fa="$1" p="$2" d; d=$(dirname "$p"); mkdir -p "$d"
  local tmp="$d/.tmp_$$"; mkdir -p "$tmp"
  # IDs unicos o bowtie mezcla conteos sin avisar
  local dup; dup=$(grep '^>' "$fa" | sed 's/^>//; s/[[:space:]].*//' | sort | uniq -d | head -3)
  [[ -z "$dup" ]] || { echo "ERROR: IDs repetidos en $fa: $dup"; exit 1; }
  echo "bowtie-build ($THREADS hilos): $fa -> $p"
  bowtie-build --threads "$THREADS" -q "$fa" "$tmp/$(basename "$p")"
  mv "$tmp"/* "$d/" && rmdir "$tmp"
  idx_ok "$p" || { echo "ERROR: indice incompleto en $p"; exit 1; }
  # autocomprobacion: 5 trozos de 22 nt de la propia referencia deben mapear
  local q="$d/.selftest.fa"
  awk '/^>/{if(seq!=""&&n<5&&length(seq)>200){print ">t"n; print substr(seq,101,22); n++} seq=""; next}{seq=seq $0}
       END{if(n<5&&length(seq)>200){print ">t"n; print substr(seq,101,22)}}' "$fa" > "$q"
  local nq hit; nq=$(grep -c '^>' "$q" || true)
  hit=$(bowtie -f -v 0 -k 1 --no-unal "$p" "$q" 2>/dev/null | wc -l)
  rm -f "$q"
  echo "  autocomprobacion: $hit/$nq trozos de 22 nt mapean sin mismatches"
  (( hit >= 1 )) || { echo "ERROR: la autocomprobacion no mapea nada"; exit 1; }
}

case "$cmd" in
  check)
    activate_env "$ENV_SRNA" >/dev/null 2>&1 || echo "  !! no se activa el entorno $ENV_SRNA"
    for t in cutadapt bowtie bowtie-build samtools python3; do
      printf '  %-13s %s\n' "$t" "$(command -v "$t" >/dev/null && "$t" --version 2>&1 | head -1 || echo 'NO ESTA')"
    done
    # seqkit no acepta --version, es 'seqkit version'
    printf '  %-13s %s\n' seqkit "$(command -v seqkit >/dev/null && seqkit version 2>&1 | head -1 || echo 'NO ESTA')"
    echo
    for p in "$BT1_GENOME" "$BT1_RRNA" "$BT1_VIRAL"; do
      if idx_ok "$p"; then echo "  [ok ] $p"; else echo "  [ -- ] $p"; fi
    done
    echo
    for f in "$GENOME_FA" "$RRNA_FA" "$VOTUS_FA"; do
      [[ -s "$f" ]] && echo "  [ok ] fasta $f" || echo "  [ -- ] fasta $f"
    done
    [[ -s "$GENOME_FA" ]] || echo "        (sin FASTA del genoma: se sacara de $BT2_INDEX con bowtie2-inspect)"
    echo
    echo "  el FASTA de 'viral' se escribira en: $(dirname "$BT1_VIRAL")/viral_targets.fa"
    [[ -n "$VIRAL_EXTRA" ]] && echo "  VIRAL_EXTRA = $VIRAL_EXTRA" \
      || echo "  !! VIRAL_EXTRA vacio: sin AEFE1 no hay control positivo de perfil EVE"
    df -h "$SRNA_REFS" 2>/dev/null | tail -1 | awk '{print "  espacio en refs: "$4" libres"}'
    ;;

  genome)
    if [[ ! -s "$GENOME_FA" ]]; then
      echo "Reconstruyendo FASTA de AaegL5 desde $BT2_INDEX"
      mkdir -p "$(dirname "$GENOME_FA")"
      activate_env "$ENV_BOWTIE2"
      bowtie2-inspect "$BT2_INDEX" > "$GENOME_FA.tmp" && mv "$GENOME_FA.tmp" "$GENOME_FA"
    fi
    activate_env "$ENV_SRNA"
    build "$GENOME_FA" "$BT1_GENOME"
    manifest "$BT1_GENOME" "$GENOME_FA" "$GENOME_FA" "BT2_INDEX=$BT2_INDEX"
    ;;

  rrna)
    [[ -s "$RRNA_FA" ]] || { echo "ERROR: no existe $RRNA_FA"; exit 1; }
    activate_env "$ENV_SRNA"
    build "$RRNA_FA" "$BT1_RRNA"
    manifest "$BT1_RRNA" "$RRNA_FA" "$RRNA_FA"
    ;;

  viral)
    [[ -s "$VOTUS_FA" ]] || { echo "ERROR: no existe $VOTUS_FA"; exit 1; }
    # GUARDIA: el nombre del directorio del indice no puede mentir sobre su contenido.
    # Sin esto, un BT1_VIRAL=...bt1_viral_min500 con el VOTUS_FA por defecto (las 748
    # vOTUs del corte >=1 kb) construiria el catalogo equivocado en silencio.
    if [[ "$BT1_VIRAL" == *min500* && "$VOTUS_FA" != *min500* ]]; then
      echo "ERROR: BT1_VIRAL apunta a min500 pero VOTUS_FA no:"
      echo "  BT1_VIRAL = $BT1_VIRAL"
      echo "  VOTUS_FA  = $VOTUS_FA"
      echo "  exporta VOTUS_FA=\$REFS_DIR/votu_catalog/MERIDA_cdhit95_85_G1g1_min500/votus.fa"
      exit 1
    fi
    activate_env "$ENV_SRNA"
    # el FASTA vive junto a SU indice: cambiar BT1_VIRAL cambia de catalogo sin
    # pisar el anterior (ver la cabecera)
    out="$(dirname "$BT1_VIRAL")/viral_targets.fa"; mkdir -p "$(dirname "$out")"
    srcs=("$VOTUS_FA")
    if [[ -n "$VIRAL_EXTRA" ]]; then IFS=',' read -ra ex <<< "$VIRAL_EXTRA"; srcs+=("${ex[@]}"); fi
    for f in "${srcs[@]}"; do [[ -s "$f" ]] || { echo "ERROR: no existe $f"; exit 1; }; done
    # cabeceras cortas (bowtie corta en el primer espacio igualmente) y en mayusculas
    cat "${srcs[@]}" | awk '/^>/{sub(/[[:space:]].*/,""); print; next}{print toupper($0)}' > "$out.tmp"
    mv "$out.tmp" "$out"
    build "$out" "$BT1_VIRAL"
    manifest "$BT1_VIRAL" "$out" "${srcs[@]}"
    ;;

  *) echo "subcomando desconocido: $cmd"; exit 2 ;;
esac
echo "Done: $cmd ($((SECONDS/60)) min)"
