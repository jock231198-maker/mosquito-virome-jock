#!/bin/bash
# ---------------------------------------------------------------------------
# votu_vclust.sh — vOTUs con Vclust (ANI y AF reales, MIUViG) y comparacion con CD-HIT
#
# Paso 1 de la metodologia v2. Agrupa los mismos 1,144 contigs virales de geNomad
# que dieron las 748 vOTUs de CD-HIT-EST, pero midiendo ANI y fraccion alineada
# (AF) por alineamiento (Zielezinski et al. 2025, Nat Methods), y compara las
# dos particiones contig a contig. NO toca nada de v1: ni refs/votu_index/, ni
# la carpeta de CD-HIT, ni datos_v3.
#
# Uso:
#   ./votu_vclust.sh check      # entradas, vclust, referencia CD-HIT, bsub listo  (login)
#   ./votu_vclust.sh run        # prefilter + align + cluster (cd-hit y leiden)     (job)
#   ./votu_vclust.sh comparar   # Vclust vs CD-HIT (y MMseqs2 si se da)             (login)
#
# Envio (LSF copia el entorno exportado):
#   bsub -J votuvc -o "$LOGS_DIR/votuvc.%J.log" -e "$LOGS_DIR/votuvc.%J.err" \
#        -q normal -n 4 -M 8000 -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "$SCRIPTS_DIR/votu_vclust.sh run"
#   (estimacion para ~1,100 contigs de 1-15 kb: minutos. Medir con bhist -l)
#
# Instalacion (una vez; solo conda-forge/bioconda, el proxy bloquea 'defaults'):
#   conda create -n vclust_1.3.1 -c conda-forge -c bioconda vclust=1.3.1 -y
#   (alternativa: pip install vclust==1.3.1 dentro de un entorno propio)
#
# DEFINICION DE vOTU (MIUViG, Roux et al. 2019)
#   >=95 % de ANI sobre >=85 % de la secuencia MAS CORTA.
#     vclust prefilter --min-ident 0.95   descarta pares que no pueden llegar al 95 %
#     vclust align                        ANI y AF por pares (ani.tsv, ambas direcciones)
#     vclust cluster --metric ani --ani 0.95 --qcov 0.85
#   En ani.tsv cada par aparece dos veces; en la fila cuya query es la secuencia
#   corta, qcov es justo la AF de la corta. Representante = el contig mas largo
#   del grupo (empates: ID menor), igual que CD-HIT.
#   Algoritmos: cd-hit (voraz de mayor a menor, comparable con v1) y leiden.
#
# Variables (exportar antes para cambiar):
#   VIRAL_FA       entrada; por defecto el *_virus.fna de geNomad sobre la union
#   VOTU_TAG       carpeta de salida; por defecto <lote>_vclust_ani95_af85
#   VC_IDENT VC_ANI VC_QCOV   0.95 0.95 0.85
#   VC_ALGOS       algoritmos de clustering ("cd-hit leiden")
#   VC_EXTRA_ALIGN flags extra de 'vclust align', tal cual
#   ENV_VCLUST     entorno con vclust (vclust_1.3.1)
#   CDHIT_CLSTR    .clstr de las 748 vOTUs de v1 (por defecto se busca)
#   MMSEQS_TSV     opcional: <prefijo>_cluster.tsv de mmseqs (rep<TAB>miembro)
#   REF_NAME       nombre de la particion de referencia en 'comparar' (v1_cdhit)
#   COMP_DIR       carpeta de 'comparar' ($OUT/comparar); cambiala para comparar
#                  contra otra referencia sin pisar la anterior
#   EXPECT_IN      1144      EXPECT_V1   748
#   FORCE=1        rehacer aunque ya exista MANIFEST
#
# SALIDA ($REFS_DIR/votu_catalog/$VOTU_TAG/):
#   fltr.txt, ani.tsv, ani.ids.tsv     prefiltro, ANI/AF por pares, IDs y longitudes
#   clusters_<algo>.tsv                salida cruda de 'vclust cluster'
#   votu_miembros_<algo>.tsv           un contig por fila: vOTU (= representante),
#                                      longitud, es_representante, n_contigs
#   votus_<algo>.fa                    representantes
#   MANIFEST.txt                       versiones, comandos, md5, conteos (al final)
#   comparar/                          resumen.txt, contig_votu_comparacion.tsv,
#                                      discordantes_<algo>.tsv, reps_<ref>_perdidos_<algo>.txt
# No construye indice bowtie2: eso se hace cuando se elija el catalogo v2.
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SELF="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")"
source "$SCRIPT_DIR/config.sh"

cmd="${1:-ayuda}"
THREADS="${LSB_DJOB_NUMPROC:-2}"
LOTE="${RESULT_FROM:-MERIDA}"

VIRAL_FA="${VIRAL_FA:-$RESULTS_DIR/genomad_all/all_final/all_final_summary/all_final_virus.fna}"
VOTU_TAG="${VOTU_TAG:-${LOTE}_vclust_ani95_af85}"
OUT="$REFS_DIR/votu_catalog/$VOTU_TAG"
VC_IDENT="${VC_IDENT:-0.95}"
VC_ANI="${VC_ANI:-0.95}"
VC_QCOV="${VC_QCOV:-0.85}"
VC_ALGOS="${VC_ALGOS:-cd-hit leiden}"
VC_EXTRA_ALIGN="${VC_EXTRA_ALIGN:-}"
ENV_VCLUST="${ENV_VCLUST:-vclust_1.3.1}"
MMSEQS_TSV="${MMSEQS_TSV:-}"
REF_NAME="${REF_NAME:-v1_cdhit}"
COMP_DIR="${COMP_DIR:-$OUT/comparar}"
EXPECT_IN="${EXPECT_IN:-1144}"
EXPECT_V1="${EXPECT_V1:-748}"
FORCE="${FORCE:-0}"

# --- helpers ----------------------------------------------------------------
die()  { echo "ERROR: $*" >&2; exit 1; }
ok()   { echo "  [ok] $*"; }
warn() { echo "  [!!] $*"; }
nseq() { awk '/^>/{n++} END{print n+0}' "$1" 2>/dev/null || echo 0; }
md5()  { md5sum < "$1" | cut -c1-12; }
ids()  { grep '^>' "$1" | sed 's/^>//; s/[[:space:]].*//'; }

# .clstr de v1: el que diga CDHIT_CLSTR, o el primero que exista
find_clstr() {
  if [[ -n "${CDHIT_CLSTR:-}" ]]; then echo "$CDHIT_CLSTR"; return; fi
  local c
  for c in "$RESULTS_DIR/qc_control/votus.clstr" \
           "$REFS_DIR/votu_catalog/${LOTE}_cdhit95_85/votus.fa.clstr" \
           "$REFS_DIR/votu_index/votus.clstr"; do
    [[ -s "$c" ]] && { echo "$c"; return; }
  done
  echo ""
}
n_clstr() { grep -c '^>Cluster' "$1" 2>/dev/null || echo 0; }

# id<TAB>longitud a partir del FASTA (primera palabra de la cabecera)
lens() {
  awk '/^>/{if (id != "") print id "\t" L; id = substr($1, 2); L = 0; next}
       {L += length($0)}
       END{if (id != "") print id "\t" L}' "$1"
}

# clusters_<algo>.tsv (object<TAB>cluster) + longitudes -> votu_miembros
# representante = el mas largo; empate -> ID menor (orden C)
miembros() {
  awk -F'\t' -v OFS='\t' '
    NR == FNR { L[$1] = $2; next }
    FNR == 1 && $2 !~ /^-?[0-9]+$/ { next }      # cabecera object/cluster
    {
      c = $2; m[c] = m[c] (m[c] == "" ? "" : SUBSEP) $1; n[c]++
      if (!(c in rep) || L[$1] > L[rep[c]] || (L[$1] == L[rep[c]] && $1 < rep[c])) rep[c] = $1
    }
    END {
      print "votu", "contig", "longitud", "es_representante", "n_contigs"
      for (c in m) { k = split(m[c], a, SUBSEP)
        for (i = 1; i <= k; i++) print rep[c], a[i], L[a[i]], (a[i] == rep[c] ? "si" : "no"), n[c] }
    }' "$1" "$2" | { IFS= read -r h; echo "$h"; LC_ALL=C sort -t$'\t' -k1,1 -k4,4r -k2,2; }
}

# extrae del FASTA las secuencias cuyos IDs estan en un fichero
extrae() {
  awk 'NR == FNR { w[$1] = 1; next }
       /^>/ { keep = (substr($1, 2) in w) }
       keep' "$1" "$2"
}

# ---------------------------------------------------------------------------
case "$cmd" in

check)
  echo "== votu_vclust.sh check   $(date '+%F %T')  $(hostname)   lote: $LOTE"
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
  echo "-- vclust ($ENV_VCLUST)"
  ( activate_env "$ENV_VCLUST" >/dev/null 2>&1
    if command -v vclust >/dev/null; then
      ok "vclust $(vclust --version 2>&1 | tail -1)"
      for f in "--min-ident" "--qcov" "--algorithm"; do
        vclust prefilter -h 2>&1 | grep -q -- "$f" || vclust cluster -h 2>&1 | grep -q -- "$f" \
          && ok "opcion $f presente" || warn "opcion $f no aparece en la ayuda"
      done
    else
      warn "vclust no esta en '$ENV_VCLUST'. Instalar:"
      echo "         conda create -n vclust_1.3.1 -c conda-forge -c bioconda vclust=1.3.1 -y"
    fi )
  echo
  cl=$(find_clstr)
  echo "-- referencia v1 (CD-HIT): ${cl:-no encontrada}"
  if [[ -n "$cl" && -s "$cl" ]]; then
    nc=$(n_clstr "$cl")
    (( nc == EXPECT_V1 )) && ok "$nc grupos" || warn "$nc grupos (esperados $EXPECT_V1)"
  else
    warn "sin .clstr: 'comparar' no podra hacer la comparacion contig a contig. Candidatos:"
    find "$RESULTS_DIR" "$REFS_DIR" -maxdepth 4 -name '*.clstr' 2>/dev/null | grep -i votu | sed 's/^/         /'
    echo "       exporta CDHIT_CLSTR=<ruta>"
  fi
  [[ -n "$MMSEQS_TSV" ]] && { [[ -s "$MMSEQS_TSV" ]] && ok "MMSEQS_TSV: $MMSEQS_TSV" || warn "no existe MMSEQS_TSV=$MMSEQS_TSV"; }
  echo
  echo "-- salida: $OUT"
  if [[ -s "$OUT/MANIFEST.txt" ]]; then
    ok "ya hecho. 'run' no hara nada sin FORCE=1"
  else
    echo "       (vacia: lista para 'run')"
  fi
  echo
  echo "-- parametros: prefilter --min-ident $VC_IDENT | cluster --metric ani --ani $VC_ANI --qcov $VC_QCOV | algoritmos: $VC_ALGOS"
  echo
  echo "-- envio:"
  cat <<EOF
  bsub -J votuvc -o "$LOGS_DIR/votuvc.%J.log" -e "$LOGS_DIR/votuvc.%J.err" \\
       -q normal -n 4 -M 8000 -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \\
       "$SELF run"
EOF
  ;;

run)
  SECONDS=0
  echo "votu_vclust.sh run — $(date '+%F %T') en $(hostname), $THREADS hilos, lote $LOTE"
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
  lens "$VIRAL_FA" > .lens.tsv

  activate_env "$ENV_VCLUST"
  command -v vclust >/dev/null || die "vclust no esta en '$ENV_VCLUST' (ver 'check')"
  VC_VER=$(vclust --version 2>&1 | tail -1)

  # --- 1. prefiltro y alineamiento (una vez) --------------------------------
  PRE=(vclust prefilter -i "$VIRAL_FA" -o fltr.txt.tmp --min-ident "$VC_IDENT" -t "$THREADS")
  echo "+ ${PRE[*]}"
  "${PRE[@]}" > prefilter.log 2>&1 || { tail -20 prefilter.log; die "fallo vclust prefilter"; }
  mv fltr.txt.tmp fltr.txt
  # shellcheck disable=SC2206
  ALN=(vclust align -i "$VIRAL_FA" -o ani.tsv --filter fltr.txt -t "$THREADS" $VC_EXTRA_ALIGN)
  echo "+ ${ALN[*]}"
  "${ALN[@]}" > align.log 2>&1 || { tail -20 align.log; die "fallo vclust align"; }
  [[ -s ani.ids.tsv ]] || die "vclust align no genero ani.ids.tsv"
  n_pairs=$(( $(wc -l < ani.tsv) - 1 ))
  n_ids=$(( $(wc -l < ani.ids.tsv) - 1 ))
  echo "alineamiento: $n_pairs filas (pares en ambas direcciones); $n_ids secuencias en ani.ids.tsv"

  # --- 2. clustering por algoritmo -----------------------------------------
  resumen_algos=""
  CLU_CMDS=()
  for algo in $VC_ALGOS; do
    CLU=(vclust cluster -i ani.tsv -o "clusters_${algo}.tsv.tmp" --ids ani.ids.tsv
         --algorithm "$algo" --metric ani --ani "$VC_ANI" --qcov "$VC_QCOV")
    echo "+ ${CLU[*]}"
    "${CLU[@]}" > "cluster_${algo}.log" 2>&1 || { tail -20 "cluster_${algo}.log"; die "fallo vclust cluster ($algo)"; }
    mv "clusters_${algo}.tsv.tmp" "clusters_${algo}.tsv"
    CLU_CMDS+=("${CLU[*]}")

    miembros .lens.tsv "clusters_${algo}.tsv" > "votu_miembros_${algo}.tsv"
    n_mem=$(( $(wc -l < "votu_miembros_${algo}.tsv") - 1 ))
    (( n_mem == n_in )) || die "$algo: votu_miembros tiene $n_mem contigs y la entrada $n_in"
    awk -F'\t' 'NR > 1 && $4 == "si" {print $2}' "votu_miembros_${algo}.tsv" > .reps
    extrae .reps "$VIRAL_FA" > "votus_${algo}.fa"
    n_out=$(nseq "votus_${algo}.fa")
    n_reps=$(wc -l < .reps)
    (( n_out == n_reps )) || die "$algo: $n_reps representantes pero solo $n_out en el FASTA (los IDs de vclust no coinciden con las cabeceras)"
    read -r n_single n_max <<< "$(awk -F'\t' 'NR > 1 && $4 == "si" {if ($5 == 1) s++; if ($5 > m) m = $5} END{print s+0, m+0}' "votu_miembros_${algo}.tsv")"
    echo "$algo: $n_in contigs -> $n_out vOTUs ($n_single singletons; mayor = $n_max contigs)"
    resumen_algos+="  $algo: $n_out vOTUs (votus_${algo}.fa md5 $(md5 "votus_${algo}.fa")); $n_single singletons; mayor = $n_max contigs"$'\n'
  done
  rm -f .reps .lens.tsv

  # --- 3. MANIFEST (ultimo: su existencia marca el paso como completo) ------
  {
    echo "paso:        vOTUs con Vclust (votu_vclust.sh run)"
    echo "fecha:       $(date -Iseconds)"
    echo "host/job:    $(hostname) / ${LSB_JOBID:-local}, $THREADS hilos"
    echo "script:      $SELF (md5 $(md5 "$SELF"))"
    echo "lote:        $LOTE"
    echo "vclust:      $VC_VER ($ENV_VCLUST)"
    echo "prefilter:   ${PRE[*]}"
    echo "align:       ${ALN[*]}"
    for c in "${CLU_CMDS[@]}"; do echo "cluster:     $c"; done
    echo "criterio:    vOTU = ANI >= $VC_ANI sobre AF >= $VC_QCOV de la secuencia corta (MIUViG); representante = el mas largo"
    echo "entrada:     $VIRAL_FA ($n_in secuencias, md5 $(md5 "$VIRAL_FA"))"
    echo "pares:       $n_pairs filas en ani.tsv"
    echo "vOTUs:"
    printf '%s' "$resumen_algos"
  } > MANIFEST.txt
  echo
  cat MANIFEST.txt
  echo
  echo "Done ($((SECONDS/60)) min). Siguiente: ./votu_vclust.sh comparar"
  ;;

comparar)
  [[ -s "$OUT/MANIFEST.txt" ]] || die "no hay corrida completa en $OUT (corre 'run')"
  cl=$(find_clstr)
  [[ -n "$cl" && -s "$cl" ]] || die "no encuentro el .clstr de v1; exporta CDHIT_CLSTR=<ruta>"
  d="$COMP_DIR"; mkdir -p "$d"
  echo "referencia: $REF_NAME = $cl"
  parts=("$REF_NAME=clstr:$cl")
  for algo in $VC_ALGOS; do
    [[ -s "$OUT/votu_miembros_${algo}.tsv" ]] && parts+=("vclust_${algo}=miembros:$OUT/votu_miembros_${algo}.tsv")
  done
  [[ -n "$MMSEQS_TSV" && -s "$MMSEQS_TSV" ]] && parts+=("mmseqs=pares:$MMSEQS_TSV")

  python3 - "$OUT/ani.tsv" "$OUT/ani.ids.tsv" "$d" "${parts[@]}" <<'PY' | tee "$d/resumen.txt"
import sys, os, csv
from collections import defaultdict
ani_tsv, ids_tsv, outd, *specs = sys.argv[1:]

def read_tsv(p):
    with open(p) as f:
        return list(csv.reader(f, delimiter="\t"))

# longitudes (ani.ids.tsv: id seq_len no_parts)
L = {r[0]: int(r[1]) for r in read_tsv(ids_tsv)[1:]}

def part_clstr(p):
    m, rep, cur, members = {}, {}, None, []
    def flush():
        if not members: return
        r = next((x for x, star in members if star), members[0][0])
        for x, _ in members: m[x] = r
    for line in open(p):
        if line.startswith(">Cluster"):
            flush(); members = []; continue
        if ">" in line:
            x = line.split(">", 1)[1].split("...", 1)[0].split()[0]
            members.append((x, line.rstrip().endswith("*")))
    flush()
    return m

def part_miembros(p):
    rows = read_tsv(p)
    h = rows[0]; iv, ic = h.index("votu"), h.index("contig")
    return {r[ic]: r[iv] for r in rows[1:]}

def part_pares(p):          # mmseqs _cluster.tsv: rep<TAB>member
    return {r[1]: r[0] for r in read_tsv(p) if len(r) >= 2}

readers = {"clstr": part_clstr, "miembros": part_miembros, "pares": part_pares}
P = {}
for s in specs:
    name, rest = s.split("=", 1)
    fmt, path = rest.split(":", 1)
    P[name] = readers[fmt](path)

names = list(P)
ref = names[0]
common = set(P[ref])
for n in names[1:]:
    common &= set(P[n])

def comb2(k): return k * (k - 1) // 2

def ari(a, b, keys):
    cont = defaultdict(int); ca = defaultdict(int); cb = defaultdict(int)
    for k in keys:
        cont[(a[k], b[k])] += 1; ca[a[k]] += 1; cb[b[k]] += 1
    s_ij = sum(comb2(v) for v in cont.values())
    s_a = sum(comb2(v) for v in ca.values()); s_b = sum(comb2(v) for v in cb.values())
    n2 = comb2(len(keys))
    exp = s_a * s_b / n2 if n2 else 0
    mx = (s_a + s_b) / 2
    return 1.0 if mx == exp else (s_ij - exp) / (mx - exp)

def groups(p, keys):
    g = defaultdict(set)
    for k in keys: g[p[k]].add(k)
    return g

print("== Comparacion de particiones (contigs en comun: %d)" % len(common))
for n in names:
    extra = set(P[n]) - common
    if extra:
        print("  [!!] %s tiene %d contigs que no estan en todas las particiones (se ignoran)" % (n, len(extra)))
print()
print("%-16s %8s %8s %11s %8s" % ("particion", "contigs", "vOTUs", "singletons", "mayor"))
for n in names:
    g = groups(P[n], common)
    sizes = [len(v) for v in g.values()]
    print("%-16s %8d %8d %11d %8d" % (n, len(common), len(g), sum(1 for x in sizes if x == 1), max(sizes)))
print()

# ANI/AF entre pares (fila con query = la secuencia corta; si empatan, la mejor)
pair = {}
with open(ani_tsv) as f:
    rd = csv.reader(f, delimiter="\t"); h = next(rd)
    iq, ir, ia, ic = h.index("query"), h.index("reference"), h.index("ani"), h.index("qcov")
    for r in rd:
        q, t = r[iq], r[ir]
        if L.get(q, 0) > L.get(t, 0): continue
        key = frozenset((q, t)); val = (float(r[ia]), float(r[ic]))
        if key not in pair or val > pair[key]: pair[key] = val
def af(a, b):
    if a == b: return ("1.000", "1.000")
    v = pair.get(frozenset((a, b)))
    return ("%.3f" % v[0], "%.3f" % v[1]) if v else ("<pref", "<pref")

gref = groups(P[ref], common)
reps_ref = set(gref)
for n in names[1:]:
    g = groups(P[n], common)
    same = sum(1 for s in map(frozenset, gref.values()) if s in set(map(frozenset, g.values())))
    split = sum(1 for s in gref.values() if len({P[n][k] for k in s}) > 1)
    merge = sum(1 for s in g.values() if len({P[ref][k] for k in s}) > 1)
    shared_reps = len(reps_ref & set(g))
    print("-- %s frente a %s" % (n, ref))
    print("   ARI                                  %.4f" % ari(P[ref], P[n], common))
    print("   vOTUs identicas (mismos contigs)     %d de %d" % (same, len(gref)))
    print("   vOTUs de %s partidas en varias     %d" % (ref, split))
    print("   vOTUs de %s que juntan varias de %s  %d" % (n, ref, merge))
    print("   representantes de %s que siguen siendolo  %d de %d" % (ref, shared_reps, len(gref)))
    lost = sorted(reps_ref - set(g))
    with open(os.path.join(outd, "reps_%s_perdidos_%s.txt" % (ref, n)), "w") as f:
        for r in lost: f.write("%s\t%d\t%s\n" % (r, L.get(r, 0), P[n][r]))
    disc = [k for k in sorted(common)
            if len({P[n][x] for x in gref[P[ref][k]]}) > 1 or len({P[ref][x] for x in g[P[n][k]]}) > 1]
    with open(os.path.join(outd, "discordantes_%s.tsv" % n), "w") as f:
        f.write("contig\tlongitud\tvotu_%s\tvotu_%s\tani_a_rep_%s\taf_a_rep_%s\tani_a_rep_%s\taf_a_rep_%s\n"
                % (ref, n, ref, ref, n, n))
        for k in disc:
            a1, f1 = af(k, P[ref][k]); a2, f2 = af(k, P[n][k])
            f.write("\t".join([k, str(L.get(k, 0)), P[ref][k], P[n][k], a1, f1, a2, f2]) + "\n")
    print("   contigs en vOTUs discordantes        %d  (discordantes_%s.tsv)" % (len(disc), n))
    print()

with open(os.path.join(outd, "contig_votu_comparacion.tsv"), "w") as f:
    f.write("contig\tlongitud\t" + "\t".join("votu_" + n for n in names) + "\n")
    for k in sorted(common):
        f.write("\t".join([k, str(L.get(k, 0))] + [P[n][k] for n in names]) + "\n")
print("Tablas en " + outd)
print("ani/af '<pref' = el par no paso el prefiltro o no hubo alineamiento")
PY
  ;;

ayuda|-h|--help)
  awk 'NR > 2 && /^# -{10}/ {exit} NR > 2 {sub(/^# ?/, ""); print}' "$SELF"
  ;;

*)
  echo "subcomando desconocido: $cmd (check | run | comparar | ayuda)"; exit 2
  ;;
esac
