#!/bin/bash
# ---------------------------------------------------------------------------
# votu_catalog.sh — catalogo de vOTUs: CD-HIT-EST 95/85 + indice bowtie2
#
# Deja por escrito el paso que se corrio a mano: 1,144 contigs virales de
# geNomad (union de 6,004) -> 748 vOTUs -> indice bowtie2 para mapear lecturas.
# NO toca refs/votu_index/ (el indice en uso): escribe en su propia carpeta, y
# 'verificar' comprueba si reproduce exactamente los 748 de ese indice.
#
# Uso:
#   ./votu_catalog.sh check       # entradas, herramientas, bsub listo     (login)
#   ./votu_catalog.sh run         # CD-HIT-EST + tablas + bowtie2-build    (job)
#   ./votu_catalog.sh verificar   # compara con refs/votu_index/votus      (login)
#
# Envio (las variables exportadas viajan al job; LSF copia el entorno):
#   bsub -J votucat -o "$LOGS_DIR/votucat.%J.log" -e "$LOGS_DIR/votucat.%J.err" \
#        -q normal -n 4 -M 4000 -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
#        "$SCRIPTS_DIR/votu_catalog.sh run"
#   (estimacion para ~1,100 contigs: minutos. Medir con bhist -l y anotar en
#    bsubs_reproducibilidad.md)
#
# DEFINICION DE vOTU
#   Umbral MIUViG (Roux et al. 2019, Nat Biotechnol): >=95 % de identidad
#   nucleotidica sobre >=85 % de la secuencia MAS CORTA. En CD-HIT-EST:
#     -c 0.95   identidad
#     -aS 0.85  cobertura de la secuencia corta (la que se agrupa)
#     -G 0      identidad local, sobre el alineamiento (no sobre toda la corta)
#     -n 10     palabra de 10 (exige -c >= 0.95; si bajas -c, baja -n)
#     -d 0      nombre completo en el .clstr (si no, lo corta a 20 caracteres)
#     -M 0      sin tope de RAM interno (el tope lo pone LSF)
#   El representante de cada vOTU es el contig MAS LARGO del grupo (CD-HIT
#   ordena por longitud). Mismos parametros que 'mapping' en build_viral_db.sh.
#
# Variables (exportar antes para cambiar):
#   VIRAL_FA      entrada; por defecto el *_virus.fna de geNomad sobre la union
#   VOTU_TAG      version = carpeta de salida; por defecto <lote>_cdhit95_85
#   CDHIT_C CDHIT_AS CDHIT_G CDHIT_N   parametros (0.95 0.85 0 10)
#   CDHIT_EXTRA   flags extra de cd-hit-est, tal cual
#   ENV_CDHIT     entorno con cd-hit-est (por defecto viraldb)
#   REF_INDEX     indice en uso con el que comparar (refs/votu_index/votus)
#   FORCE=1       rehacer aunque la salida ya exista
#
# SALIDA ($REFS_DIR/votu_catalog/$VOTU_TAG/):
#   votus.fa               representantes = referencia de mapeo (+ .fai)
#   votus.fa.clstr         salida cruda de CD-HIT
#   votu_miembros.tsv      un contig por fila: su vOTU, identidad, hebra
#   votu_resumen.tsv       un vOTU por fila: longitud, n.o de contigs y en que
#                          muestras se ENSAMBLO (prefijo <muestra>__ de la cabecera;
#                          NA si la cabecera no lo lleva). Ojo: eso no es presencia;
#                          la presencia sale del mapeo (cov_worker.sh).
#   votus.*.bt2            indice bowtie2
#   MANIFEST.txt           versiones, comando, md5 y conteos. Se escribe al final:
#                          si existe, el paso esta completo.
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SELF="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")"
source "$SCRIPT_DIR/config.sh"

cmd="${1:-ayuda}"
THREADS="${LSB_DJOB_NUMPROC:-2}"
LOTE="${RESULT_FROM:-MERIDA}"

VIRAL_FA="${VIRAL_FA:-$RESULTS_DIR/genomad_all/all_final/all_final_summary/all_final_virus.fna}"
VOTU_TAG="${VOTU_TAG:-${LOTE}_cdhit95_85}"
OUT="$REFS_DIR/votu_catalog/$VOTU_TAG"
CDHIT_C="${CDHIT_C:-0.95}"
CDHIT_AS="${CDHIT_AS:-0.85}"
CDHIT_G="${CDHIT_G:-0}"
CDHIT_N="${CDHIT_N:-10}"
CDHIT_EXTRA="${CDHIT_EXTRA:-}"
ENV_CDHIT="${ENV_CDHIT:-viraldb}"
REF_INDEX="${REF_INDEX:-$REFS_DIR/votu_index/votus}"
REF_FA="${REF_FA:-$REFS_DIR/votu_index/votus_from_index.fasta}"
EXPECT_IN="${EXPECT_IN:-1144}"
EXPECT_OUT="${EXPECT_OUT:-748}"
FORCE="${FORCE:-0}"

# --- helpers ----------------------------------------------------------------
die()  { echo "ERROR: $*" >&2; exit 1; }
ok()   { echo "  [ok] $*"; }
warn() { echo "  [!!] $*"; }
nseq() { awk '/^>/{n++} END{print n+0}' "$1" 2>/dev/null || echo 0; }
md5()  { md5sum < "$1" | cut -c1-12; }
ids()  { grep '^>' "$1" | sed 's/^>//; s/[[:space:]].*//'; }
idx_ok() { [[ -s "$1.1.bt2" || -s "$1.1.bt2l" ]] && [[ -s "$1.rev.1.bt2" || -s "$1.rev.1.bt2l" ]]; }

# id<TAB>secuencia, en una linea, mayusculas y no-ACGT -> N (asi lo guarda
# bowtie2, de modo que un FASTA y su bowtie2-inspect comparan igual).
norm() {
  awk '/^>/{if (id != "") print id "\t" s; id = substr($1, 2); s = ""; next}
       {s = s toupper($0)}
       END{if (id != "") print id "\t" s}' "$1" \
  | awk -F'\t' '{gsub(/[^ACGT]/, "N", $2); print $1 "\t" $2}' \
  | LC_ALL=C sort -t$'\t' -k1,1
}

# .clstr de CD-HIT -> un contig por fila
#   0	12345nt, >M1_S1__NODE_1... *
#   1	3000nt, >M7_S7__NODE_9... at +/98.50%          (o  at 1:3000:5:3004/+/98.50%)
parse_clstr() {
  awk -v OFS='\t' '
    function flush(   i, rep) {
      if (n == 0) return
      for (i = 1; i <= n; i++) if (es[i] == "si") rep = id[i]
      for (i = 1; i <= n; i++) print rep, id[i], len[i], es[i], hebra[i], ident[i], n
      n = 0
    }
    BEGIN { print "votu", "contig", "longitud", "es_representante", "hebra", "identidad", "n_contigs" }
    /^>Cluster/ { flush(); next }
    NF >= 3 {
      n++
      l = $2; sub(/nt,$/, "", l); len[n] = l
      x = $0; sub(/^[^>]*>/, "", x); sub(/\.\.\..*$/, "", x); id[n] = x
      if ($NF == "*") { es[n] = "si"; hebra[n] = "+"; ident[n] = "100.00" }
      else { es[n] = "no"; m = $NF; sub(/%$/, "", m); k = split(m, p, "/"); ident[n] = p[k]; hebra[n] = p[k-1] }
    }
    END { flush() }' "$1"
}

# votu_miembros.tsv -> un vOTU por fila, de mayor a menor longitud
resumen() {
  awk -F'\t' -v OFS='\t' '
    NR == 1 { next }
    {
      if ($4 == "si") L[$1] = $3
      N[$1] = $7
      i = index($2, "__"); s = (i > 1) ? substr($2, 1, i - 1) : "NA"
      if (!(($1, s) in visto)) { visto[$1, s] = 1; M[$1] = (M[$1] == "" ? s : M[$1] "," s); nm[$1]++ }
    }
    END { for (v in L) print v, L[v], N[v], nm[v], M[v] }' "$1" \
  | sort -t$'\t' -k2,2nr -k1,1 \
  | { printf 'votu\tlongitud\tn_contigs\tn_muestras_ensamblado\tmuestras\n'; cat; }
}

cdhit_cmd() {
  # shellcheck disable=SC2206   # CDHIT_EXTRA se parte en palabras a proposito
  CDHIT=(cd-hit-est -i "$VIRAL_FA" -o votus.fa.tmp -c "$CDHIT_C" -aS "$CDHIT_AS"
         -G "$CDHIT_G" -n "$CDHIT_N" -d 0 -M 0 -T "$THREADS" $CDHIT_EXTRA)
}

# ---------------------------------------------------------------------------
case "$cmd" in

check)
  echo "== votu_catalog.sh check   $(date '+%F %T')  $(hostname)   lote: $LOTE"
  [[ "$LOTE" == MERIDA ]] || warn "el lote activo no es MERIDA (source env_lote.sh MERIDA)"
  echo
  echo "-- entrada: $VIRAL_FA"
  if [[ -s "$VIRAL_FA" ]]; then
    n=$(nseq "$VIRAL_FA")
    if (( n == EXPECT_IN )); then ok "$n secuencias"; else warn "$n secuencias (esperadas $EXPECT_IN)"; fi
    dup=$(ids "$VIRAL_FA" | sort | uniq -d | head -3 | tr '\n' ' ')
    [[ -z "$dup" ]] && ok "IDs unicos" || warn "IDs repetidos: $dup"
    echo "       ejemplo de cabecera: $(grep -m1 '^>' "$VIRAL_FA" | cut -c1-80)"
  else
    warn "no existe. Candidatos:"
    find "$RESULTS_DIR" -maxdepth 4 -name '*_virus.fna' 2>/dev/null | sed 's/^/         /'
    echo "       exporta VIRAL_FA=<ruta>"
  fi
  echo
  echo "-- herramientas"
  ( activate_env "$ENV_CDHIT" >/dev/null 2>&1
    command -v cd-hit-est >/dev/null \
      && ok "cd-hit-est ($ENV_CDHIT): $(cd-hit-est -h 2>&1 | grep -m1 -o 'version [0-9.]*')" \
      || warn "cd-hit-est no esta en el entorno '$ENV_CDHIT'" )
  ( activate_env "$ENV_BOWTIE2" >/dev/null 2>&1
    for t in bowtie2 bowtie2-build bowtie2-inspect samtools; do
      command -v "$t" >/dev/null && ok "$t ($ENV_BOWTIE2)" || warn "$t no esta en '$ENV_BOWTIE2'"
    done )
  echo
  echo "-- indice en uso: $REF_INDEX"
  idx_ok "$REF_INDEX" && ok "indice completo" || warn "no se encuentra"
  [[ -s "$REF_FA" ]] && ok "$(basename "$REF_FA"): $(nseq "$REF_FA") secuencias" \
                     || warn "sin $REF_FA ('verificar' lo sacara con bowtie2-inspect)"
  echo
  echo "-- salida: $OUT"
  if [[ -s "$OUT/MANIFEST.txt" ]]; then
    ok "ya hecho ($(nseq "$OUT/votus.fa") vOTUs). 'run' no hara nada sin FORCE=1"
  else
    echo "       (vacia: lista para 'run')"
  fi
  echo
  echo "-- parametros: -c $CDHIT_C -aS $CDHIT_AS -G $CDHIT_G -n $CDHIT_N -d 0 -M 0 $CDHIT_EXTRA"
  echo
  echo "-- envio:"
  cat <<EOF
  bsub -J votucat -o "$LOGS_DIR/votucat.%J.log" -e "$LOGS_DIR/votucat.%J.err" \\
       -q normal -n 4 -M 4000 -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \\
       "$SELF run"
EOF
  ;;

run)
  echo "votu_catalog.sh run — $(date '+%F %T') en $(hostname), $THREADS hilos, lote $LOTE"
  [[ -s "$VIRAL_FA" ]] || die "no existe VIRAL_FA=$VIRAL_FA (corre 'check')"
  if [[ -s "$OUT/MANIFEST.txt" && "$FORCE" != 1 ]]; then
    echo "Ya existe $OUT/MANIFEST.txt; nada que hacer (FORCE=1 para rehacer)."
    exit 0
  fi
  mkdir -p "$OUT" && cd "$OUT" || die "no puedo usar $OUT"
  rm -f MANIFEST.txt

  n_in=$(nseq "$VIRAL_FA")
  dup=$(ids "$VIRAL_FA" | sort | uniq -d | head -3 | tr '\n' ' ')
  [[ -z "$dup" ]] || die "IDs repetidos en la entrada: $dup"
  echo "entrada: $n_in contigs ($VIRAL_FA)"

  # --- 1. CD-HIT-EST --------------------------------------------------------
  activate_env "$ENV_CDHIT"
  command -v cd-hit-est >/dev/null || die "cd-hit-est no esta en '$ENV_CDHIT'"
  CDHIT_VER=$(cd-hit-est -h 2>&1 | grep -m1 -o 'version [0-9.]*')
  cdhit_cmd
  echo "+ ${CDHIT[*]}"
  "${CDHIT[@]}" > cdhit.log 2>&1 || { tail -20 cdhit.log; die "fallo cd-hit-est"; }
  mv votus.fa.tmp.clstr votus.fa.clstr && mv votus.fa.tmp votus.fa || die "no se genero la salida de CD-HIT"
  n_out=$(nseq votus.fa)
  echo "CD-HIT: $n_in contigs -> $n_out vOTUs"
  (( n_out == EXPECT_OUT )) || warn "esperaba $EXPECT_OUT vOTUs; corre 'verificar' para ver en que difiere"

  # --- 2. tablas ------------------------------------------------------------
  parse_clstr votus.fa.clstr > votu_miembros.tsv
  resumen votu_miembros.tsv > votu_resumen.tsv
  n_tab=$(( $(wc -l < votu_resumen.tsv) - 1 ))
  n_mem=$(( $(wc -l < votu_miembros.tsv) - 1 ))
  (( n_tab == n_out )) || die "votu_resumen.tsv tiene $n_tab vOTUs y votus.fa $n_out"
  (( n_mem == n_in  )) || die "votu_miembros.tsv tiene $n_mem contigs y la entrada $n_in"
  read -r n_single n_max <<< "$(awk -F'\t' 'NR>1{if($3==1)s++; if($3>m)m=$3} END{print s+0, m+0}' votu_resumen.tsv)"
  echo "tablas: $n_single vOTUs de un solo contig; el mayor agrupa $n_max contigs"

  # --- 3. indice bowtie2 (atomico: se construye aparte y se mueve) ----------
  activate_env "$ENV_BOWTIE2"
  BT2_VER=$(bowtie2 --version 2>&1 | grep -m1 -o 'version [0-9.]*')
  tmp=".bt2_tmp_$$"; rm -rf "$tmp"; mkdir -p "$tmp"
  bowtie2-build --threads "$THREADS" votus.fa "$tmp/votus" > bt2build.log 2>&1 \
    || { tail -20 bt2build.log; die "fallo bowtie2-build"; }
  rm -f votus.*.bt2 votus.*.bt2l
  mv "$tmp"/votus.* . && rmdir "$tmp"
  idx_ok votus || die "indice incompleto"
  n_idx=$(bowtie2-inspect -n votus | wc -l)
  (( n_idx == n_out )) || die "el indice tiene $n_idx secuencias y votus.fa $n_out"
  samtools faidx votus.fa 2>/dev/null || warn "samtools faidx fallo (no es grave)"

  # autocomprobacion: 100 nt del centro de 5 vOTUs >= 1 kb deben volver a su vOTU
  awk '/^>/{if (id != "" && n < 5 && length(s) >= 1000) {print ">" id; print substr(s, int(length(s)/2), 100); n++}
            id = substr($1, 2); s = ""; next}
       {s = s $0}
       END{if (id != "" && n < 5 && length(s) >= 1000) {print ">" id; print substr(s, int(length(s)/2), 100)}}' \
       votus.fa > .selftest.fa
  nq=$(nseq .selftest.fa)
  hit=$(bowtie2 -f -x votus -U .selftest.fa --end-to-end --no-unal --no-hd -p "$THREADS" 2>/dev/null \
        | awk -F'\t' '$1 == $3 {n++} END{print n+0}')
  rm -f .selftest.fa
  echo "autocomprobacion: $hit/$nq trozos de 100 nt vuelven a su propio vOTU"
  (( hit >= 1 )) || die "la autocomprobacion no recupera nada"

  # --- 4. MANIFEST (ultimo: su existencia marca el paso como completo) ------
  {
    echo "paso:        catalogo de vOTUs (votu_catalog.sh run)"
    echo "fecha:       $(date -Iseconds)"
    echo "host/job:    $(hostname) / ${LSB_JOBID:-local}"
    echo "script:      $SELF (md5 $(md5 "$SELF"))"
    echo "lote:        $LOTE"
    echo "cd-hit-est:  $CDHIT_VER ($ENV_CDHIT)"
    echo "bowtie2:     $BT2_VER ($ENV_BOWTIE2)"
    echo "comando:     ${CDHIT[*]}"
    echo "criterio:    vOTU = >=${CDHIT_C} identidad sobre >=${CDHIT_AS} de la secuencia corta (MIUViG); representante = el mas largo"
    echo "entrada:     $VIRAL_FA ($n_in secuencias, md5 $(md5 "$VIRAL_FA"))"
    echo "vOTUs:       $n_out (votus.fa md5 $(md5 votus.fa))"
    echo "agrupacion:  $n_single singletons; mayor vOTU = $n_max contigs"
    echo "indice:      $OUT/votus ($n_idx secuencias)"
  } > MANIFEST.txt
  echo
  cat MANIFEST.txt
  echo
  echo "Done ($((SECONDS/60)) min). Siguiente: ./votu_catalog.sh verificar"
  ;;

verificar)
  new="$OUT/votus.fa"
  [[ -s "$new" ]] || die "no existe $new: corre 'run' primero"
  if [[ ! -s "$REF_FA" ]]; then
    idx_ok "$REF_INDEX" || die "no existe ni $REF_FA ni el indice $REF_INDEX"
    activate_env "$ENV_BOWTIE2" >/dev/null
    REF_FA="$OUT/verificacion/ref_from_index.fa"; mkdir -p "$(dirname "$REF_FA")"
    bowtie2-inspect "$REF_INDEX" > "$REF_FA" || die "fallo bowtie2-inspect"
  fi
  d="$OUT/verificacion"; mkdir -p "$d"
  norm "$REF_FA" > "$d/ref.tsv"
  norm "$new"    > "$d/nuevo.tsv"

  n_ref=$(wc -l < "$d/ref.tsv"); n_new=$(wc -l < "$d/nuevo.tsv")
  LC_ALL=C comm -23 <(cut -f1 "$d/ref.tsv") <(cut -f1 "$d/nuevo.tsv") > "$d/solo_en_indice.txt"
  LC_ALL=C comm -13 <(cut -f1 "$d/ref.tsv") <(cut -f1 "$d/nuevo.tsv") > "$d/solo_en_nuevo.txt"
  LC_ALL=C join -t$'\t' "$d/ref.tsv" "$d/nuevo.tsv" \
    | awk -F'\t' -v OFS='\t' '$2 != $3 {print $1, length($2), length($3)}' > "$d/secuencia_distinta.tsv"
  a=$(wc -l < "$d/solo_en_indice.txt"); b=$(wc -l < "$d/solo_en_nuevo.txt"); c=$(wc -l < "$d/secuencia_distinta.tsv")

  echo "== verificar: $VOTU_TAG frente al indice en uso"
  printf '  %-28s %6d   %s\n' "vOTUs en el indice en uso" "$n_ref" "$REF_FA"
  printf '  %-28s %6d   %s\n' "vOTUs regenerados"          "$n_new" "$new"
  printf '  %-28s %6d\n'      "solo en el indice"          "$a"
  printf '  %-28s %6d\n'      "solo en el regenerado"      "$b"
  printf '  %-28s %6d\n'      "mismo ID, otra secuencia"   "$c"
  rm -f "$d/ref.tsv" "$d/nuevo.tsv"
  if (( a + b + c == 0 )); then
    echo
    echo "IDENTICO: el indice en uso se reproduce con este script y estos parametros."
    echo "IDENTICO $(date -Iseconds) frente a $REF_FA" >> "$OUT/MANIFEST.txt"
    exit 0
  fi
  echo
  echo "DIFIEREN. Listas en $d/"
  echo "  - Si solo cambia el numero de vOTUs, prueba otra identidad: por ejemplo"
  echo "      CDHIT_G=1 VOTU_TAG=${LOTE}_cdhit95_85_G1 ./votu_catalog.sh run   (en un bsub)"
  echo "  - Si los IDs coinciden y cambian las secuencias, revisa que VIRAL_FA sea"
  echo "    exactamente el FASTA que se uso a mano (md5 en MANIFEST.txt)."
  exit 1
  ;;

ayuda|-h|--help)
  awk 'NR > 2 && /^# -{10}/ {exit} NR > 2 {sub(/^# ?/, ""); print}' "$SELF"
  ;;

*)
  echo "subcomando desconocido: $cmd (check | run | verificar | ayuda)"; exit 2
  ;;
esac
