#!/bin/bash
# ---------------------------------------------------------------------------
# virsorter2_union.sh — VirSorter2 2.2.4 sobre la union de contigs
#
# Segunda opinion frente a geNomad (1,144 virales de los 6,004 de la union).
# Incluye el grupo RNA, que VirSorter2 NO usa por defecto (solo dsDNAphage,ssDNA).
# Probado con 8seq.fa (29 oct): RNA-linear pasa de 0.500 (dsDNAphage, 0 hallmarks)
# a 1.000 en RNA (1 hallmark) al anadir el grupo.
#
# Uso:
#   ./virsorter2_union.sh check   # entrada, base, entorno, salida; imprime el bsub (login)
#   ./virsorter2_union.sh run     # virsorter run + copia a results + MANIFEST      (job)
#
# Envio (LSF copia el entorno; hace falta LSB_DEFAULT_USERGROUP=team222 -> start.sh):
#   bsub -J vs2 -o "$LOGS_DIR/vs2.%J.log" -e "$LOGS_DIR/vs2.%J.err" \
#        -q long -n 16 -M 16000 -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \
#        "$SCRIPTS_DIR/virsorter2_union.sh run"
#   (sin medir: la prueba de 8 secuencias / 1.6 Mb tardo ~14 min con 4 hilos.
#    Medir con bhist -l y anotar en bsubs_reproducibilidad.md)
#
# Base: $REFS_DIR/virsorter2_db (virsorter setup, 28 oct). Los sub-entornos conda
#   ya estan construidos en virsorter2_db/conda_envs, asi que el job no necesita
#   internet. Para que se construyeran hubo que quitar '- defaults' de
#   site-packages/virsorter/envs/*.yaml (el proxy bloquea repo.anaconda.com).
#
# Variables (exportar antes para cambiar):
#   VS2_IN        entrada; por defecto $RESULTS_DIR/contig_annot/contigs_union_MERIDA.fasta
#                 (la que construye diamond_contigs.sh; cabecera <muestra>__<ensamblador>__<contig>)
#   VS2_MINLEN    --min-length (1000)
#   VS2_TAG       nombre de la corrida = carpeta; por defecto <nombre de VS2_IN>_min<VS2_MINLEN>
#   VS2_MINSCORE  --min-score (0.5, el de VirSorter2)
#   VS2_GROUPS    --include-groups (dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae = todos)
#   VS2_EXTRA     flags extra de 'virsorter run', tal cual
#   VS2_DB        base ($REFS_DIR/virsorter2_db)
#   ENV_VS2       entorno conda (virsorter2_2.2.4)
#   EXPECT_IN     contigs esperados en la entrada (6004; 0 = no comprobar)
#   FORCE=1       rehacer aunque ya exista MANIFEST
#
# SALIDA
#   $SCRATCH/virsorter2/$VS2_TAG/        directorio de trabajo de VirSorter2 (iter-0/, log/)
#   $RESULTS_DIR/virsorter2/$VS2_TAG/    final-viral-score.tsv, final-viral-combined.fa,
#                                        final-viral-boundary.tsv, config.yaml, MANIFEST.txt
#   MANIFEST.txt se escribe al final: si existe, el paso esta completo.
#   Si el job muere, relanzar 'run' reanuda (snakemake --rerun-incomplete).
#   Si cambias parametros, cambia tambien VS2_TAG para no mezclar corridas.
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SELF="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")"
source "$SCRIPT_DIR/config.sh"

cmd="${1:-ayuda}"
THREADS="${LSB_DJOB_NUMPROC:-4}"
LOTE="${RESULT_FROM:-MERIDA}"

VS2_IN="${VS2_IN:-$RESULTS_DIR/contig_annot/contigs_union_MERIDA.fasta}"
VS2_MINLEN="${VS2_MINLEN:-1000}"
VS2_TAG="${VS2_TAG:-$(basename "${VS2_IN%.*}")_min${VS2_MINLEN}}"
VS2_MINSCORE="${VS2_MINSCORE:-0.5}"
VS2_GROUPS="${VS2_GROUPS:-dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae}"
VS2_EXTRA="${VS2_EXTRA:-}"
VS2_DB="${VS2_DB:-$REFS_DIR/virsorter2_db}"
ENV_VS2="${ENV_VS2:-virsorter2_2.2.4}"
EXPECT_IN="${EXPECT_IN:-6004}"
FORCE="${FORCE:-0}"

WORK="$SCRATCH/virsorter2/$VS2_TAG"
OUT="$RESULTS_DIR/virsorter2/$VS2_TAG"
FINALS=(final-viral-score.tsv final-viral-combined.fa final-viral-boundary.tsv config.yaml)

# --- helpers ----------------------------------------------------------------
die()  { echo "ERROR: $*" >&2; exit 1; }
ok()   { echo "  [ok] $*"; }
warn() { echo "  [!!] $*"; }
nseq() { awk '/^>/{n++} END{print n+0}' "$1" 2>/dev/null || echo 0; }
md5()  { md5sum < "$1" | cut -c1-12; }
minlen() { awk '/^>/{if (id) print l; id = 1; l = 0; next} {l += length($0)} END{if (id) print l}' "$1" | sort -n | head -1; }
db_ok() { [[ -d "$VS2_DB/group/RNA" && -d "$VS2_DB/hmm" && -d "$VS2_DB/conda_envs" ]]; }

vs2_cmd() {
  # shellcheck disable=SC2206   # VS2_EXTRA se parte en palabras a proposito
  VS2=(virsorter run -w "$WORK" -i "$VS2_IN" --db-dir "$VS2_DB"
       --include-groups "$VS2_GROUPS" --min-length "$VS2_MINLEN"
       --min-score "$VS2_MINSCORE" --keep-original-seq -j "$THREADS" $VS2_EXTRA all)
}

# final-viral-score.tsv -> conteos por max_score_group y por tipo (full/partial/lt2gene)
conteos() {
  awk -F'\t' '
    NR == 1 { for (i = 1; i <= NF; i++) if ($i == "max_score_group") c = i; next }
    { n++; g[$c]++
      t = $1; sub(/^.*\|\|/, "", t); sub(/^[0-9]+_/, "", t); tipo[t]++ }
    END {
      printf "total %d |", n
      for (k in g) printf " %s=%d", k, g[k]
      printf " |"
      for (k in tipo) printf " %s=%d", k, tipo[k]
      print ""
    }' "$1"
}

# ---------------------------------------------------------------------------
case "$cmd" in

check)
  echo "== virsorter2_union.sh check   $(date '+%F %T')  $(hostname)   lote: $LOTE"
  [[ "$LOTE" == MERIDA ]] || warn "el lote activo no es MERIDA (source env_lote.sh MERIDA)"
  [[ -n "${LSB_DEFAULT_USERGROUP:-}" ]] && ok "LSB_DEFAULT_USERGROUP=$LSB_DEFAULT_USERGROUP" \
                                        || warn "falta LSB_DEFAULT_USERGROUP (source ./start.sh)"
  echo
  echo "-- entrada: $VS2_IN"
  if [[ -s "$VS2_IN" ]]; then
    n=$(nseq "$VS2_IN")
    if (( EXPECT_IN == 0 )); then ok "$n secuencias (sin comprobar)"
    elif (( n == EXPECT_IN )); then ok "$n secuencias"
    else warn "$n secuencias (esperadas $EXPECT_IN)"; fi
    ml=$(minlen "$VS2_IN")
    (( ml >= VS2_MINLEN )) && ok "longitud minima $ml (>= --min-length $VS2_MINLEN: no se pierde nada)" \
                           || warn "longitud minima $ml: VirSorter2 descartara las < $VS2_MINLEN"
    dup=$(grep '^>' "$VS2_IN" | sed 's/^>//; s/[[:space:]].*//' | sort | uniq -d | head -3 | tr '\n' ' ')
    [[ -z "$dup" ]] && ok "IDs unicos" || warn "IDs repetidos: $dup"
    nsp=$(grep -c '^>.*[[:space:]]' "$VS2_IN")
    (( nsp == 0 )) || warn "$nsp cabeceras con espacios (VirSorter2 se queda con la primera palabra)"
    echo "       ejemplo de cabecera: $(grep -m1 '^>' "$VS2_IN" | cut -c1-80)"
  else
    warn "no existe. Exporta VS2_IN=<ruta>"
  fi
  echo
  echo "-- entorno: $ENV_VS2"
  ( activate_env "$ENV_VS2" >/dev/null 2>&1
    if command -v virsorter >/dev/null; then
      ok "virsorter: $(virsorter --version 2>&1 | tail -1)"
      command -v mamba >/dev/null && ok "mamba en el entorno" || warn "falta mamba en el entorno"
      d=$(python -c 'import virsorter, os; print(os.path.dirname(virsorter.__file__))' 2>/dev/null)
      if [[ -n "$d" ]]; then
        grep -lE '^\s*-\s*defaults\s*$' "$d"/envs/*.yaml >/dev/null 2>&1 \
          && warn "algun envs/*.yaml aun lista 'defaults' (el proxy lo bloquea)" \
          || ok "envs/*.yaml sin canal 'defaults'"
        cfg=$(grep -m1 '^DBDIR:' "$d/template-config.yaml" 2>/dev/null | awk '{print $2}')
        [[ "$cfg" == "$VS2_DB" ]] && ok "DBDIR del config = VS2_DB" \
                                  || warn "DBDIR del config es '$cfg' (se pasa --db-dir de todos modos)"
      fi
    else
      warn "virsorter no esta en '$ENV_VS2'"
    fi )
  echo
  echo "-- base: $VS2_DB"
  db_ok && ok "group/RNA, hmm/ y conda_envs/ presentes ($(du -sh "$VS2_DB" 2>/dev/null | cut -f1))" \
        || warn "base incompleta (falta group/RNA, hmm/ o conda_envs/)"
  echo
  echo "-- salida"
  echo "       trabajo:    $WORK"
  echo "       resultados: $OUT"
  if [[ -s "$OUT/MANIFEST.txt" ]]; then
    ok "ya hecho. 'run' no hara nada sin FORCE=1"
  elif [[ -d "$WORK" ]]; then
    warn "existe $WORK sin MANIFEST: 'run' reanudara esa corrida"
  else
    echo "       (vacia: lista para 'run')"
  fi
  [[ -d "$LOGS_DIR" ]] || warn "no existe LOGS_DIR=$LOGS_DIR"
  echo
  vs2_cmd
  echo "-- comando (con los hilos del job): ${VS2[*]}"
  echo
  echo "-- envio:"
  cat <<EOF
  bsub -J vs2 -o "$LOGS_DIR/vs2.%J.log" -e "$LOGS_DIR/vs2.%J.err" \\
       -q long -n 16 -M 16000 -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \\
       "$SELF run"
EOF
  ;;

run)
  echo "virsorter2_union.sh run — $(date '+%F %T') en $(hostname), $THREADS hilos, lote $LOTE"
  if [[ -s "$OUT/MANIFEST.txt" && "$FORCE" != 1 ]]; then
    echo "Ya existe $OUT/MANIFEST.txt; nada que hacer (FORCE=1 para rehacer)."
    exit 0
  fi
  [[ -s "$VS2_IN" ]] || die "no existe VS2_IN=$VS2_IN (corre 'check')"
  db_ok || die "base incompleta en $VS2_DB"
  n_in=$(nseq "$VS2_IN")
  if (( EXPECT_IN > 0 && n_in != EXPECT_IN )); then
    die "la entrada tiene $n_in secuencias y se esperaban $EXPECT_IN (EXPECT_IN=0 para no comprobar)"
  fi
  echo "entrada: $n_in contigs ($VS2_IN)"

  activate_env "$ENV_VS2"
  command -v virsorter >/dev/null || die "virsorter no esta en '$ENV_VS2'"
  VS2_VER=$(virsorter --version 2>&1 | tail -1)

  if [[ "$FORCE" == 1 ]]; then rm -rf "$WORK" "$OUT"; fi
  mkdir -p "$(dirname "$WORK")" || die "no puedo crear $(dirname "$WORK")"
  [[ -d "$WORK" ]] && echo "reanudando $WORK"

  vs2_cmd
  echo "+ ${VS2[*]}"
  "${VS2[@]}" || die "fallo virsorter run (mira $WORK/log/ y $LOGS_DIR/vs2.${LSB_JOBID:-}.err)"

  for f in "${FINALS[@]}"; do [[ -e "$WORK/$f" ]] || die "falta $WORK/$f"; done
  mkdir -p "$OUT" || die "no puedo crear $OUT"
  for f in "${FINALS[@]}"; do cp "$WORK/$f" "$OUT/$f.tmp" && mv "$OUT/$f.tmp" "$OUT/$f"; done

  n_vir=$(( $(wc -l < "$OUT/final-viral-score.tsv") - 1 ))
  resumen=$(conteos "$OUT/final-viral-score.tsv")
  echo "virales: $resumen"

  {
    echo "paso:        VirSorter2 sobre la union (virsorter2_union.sh run)"
    echo "fecha:       $(date -Iseconds)"
    echo "host/job:    $(hostname) / ${LSB_JOBID:-local}, $THREADS hilos"
    echo "script:      $SELF (md5 $(md5 "$SELF"))"
    echo "lote:        $LOTE"
    echo "virsorter:   $VS2_VER ($ENV_VS2)"
    echo "base:        $VS2_DB"
    echo "comando:     ${VS2[*]}"
    echo "entrada:     $VS2_IN ($n_in secuencias, md5 $(md5 "$VS2_IN"))"
    echo "virales:     $n_vir  [$resumen]"
    echo "trabajo:     $WORK"
    echo "nota:        el campo max_score_group NO es taxonomia fiable (aviso de VirSorter2)"
  } > "$OUT/MANIFEST.txt"
  echo
  cat "$OUT/MANIFEST.txt"
  echo
  echo "Done ($((SECONDS/60)) min)."
  ;;

ayuda|-h|--help)
  awk 'NR > 2 && /^# -{10}/ {exit} NR > 2 {sub(/^# ?/, ""); print}' "$SELF"
  ;;

*)
  echo "subcomando desconocido: $cmd (check | run | ayuda)"; exit 2
  ;;
esac
