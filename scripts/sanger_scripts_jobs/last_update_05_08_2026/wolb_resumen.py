#!/usr/bin/env python3
"""
Resume el blastn de los contigs contra los genomas de Wolbachia.

Entradas:
  hits.tsv  blastn -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
  refs.tsv  accesion <TAB> cepa <TAB> longitud

Salidas por stdout: dos tablas.
  1. por muestra: contigs, bp alineado, mejor referencia, cepa, amplitud de cobertura
  2. por muestra y cepa: bp alineado, para ver coinfecciones A+B

La amplitud (breadth) es la fraccion del genoma de referencia cubierta por la
union de los intervalos alineados. Es la cifra que separa una infeccion real de
hits espurios: un genoma presente reparte cobertura por todo el cromosoma, los
hits en regiones conservadas se apilan en el operon ribosomal y tres genes.

    python3 -I wolb_resumen.py hits.tsv refs.tsv [ident_min] [len_min]
"""

import collections
import sys

if len(sys.argv) < 3:
    sys.exit(__doc__)

HITS, REFS = sys.argv[1], sys.argv[2]
IDENT_MIN = float(sys.argv[3]) if len(sys.argv) > 3 else 90.0
LEN_MIN = int(sys.argv[4]) if len(sys.argv) > 4 else 200

cepa, largo = {}, {}
with open(REFS) as fh:
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if len(p) >= 3 and p[0]:
            cepa[p[0]] = p[1]
            try:
                largo[p[0]] = int(p[2])
            except ValueError:
                largo[p[0]] = 0

# muestra -> ref -> lista de intervalos en la referencia
iv = collections.defaultdict(lambda: collections.defaultdict(list))
ctgs = collections.defaultdict(set)
ident = collections.defaultdict(float)
desconocidas = set()

n = 0
with open(HITS) as fh:
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 10:
            continue
        q, s = p[0], p[1]
        try:
            pid, alen = float(p[2]), int(p[3])
            ss, se = int(p[6]), int(p[7])
        except ValueError:
            continue
        if pid < IDENT_MIN or alen < LEN_MIN:
            continue
        # el ID de contig es muestra__ensamblador__nombre
        m = q.split("__")[0]
        if s not in cepa:
            desconocidas.add(s)
        iv[m][s].append((min(ss, se), max(ss, se)))
        ctgs[m].add(q)
        ident[m] = max(ident[m], pid)
        n += 1

if desconocidas:
    print(f"AVISO: {len(desconocidas)} referencias sin entrada en refs.tsv, "
          f"p.ej. {sorted(desconocidas)[:3]}", file=sys.stderr)


def union(ivs):
    """bp cubiertos por la union de intervalos cerrados [a, b]"""
    if not ivs:
        return 0
    ivs = sorted(ivs)
    tot = 0
    a, b = ivs[0]
    for x, y in ivs[1:]:
        if x > b + 1:
            tot += b - a + 1
            a, b = x, y
        else:
            b = max(b, y)
    tot += b - a + 1
    return tot


print(f"# hits usados: {n}  (identidad >= {IDENT_MIN} %, longitud >= {LEN_MIN} bp)")
print()
print("== por muestra ==")
print("muestra\tcontigs\tbp_union\tident_max\tmejor_ref\tcepa\tbreadth_pct")

filas = []
for m, refs in iv.items():
    # la mejor referencia es la que acumula mas bp en union
    mejor, mejor_bp = None, -1
    for s, ivs in refs.items():
        bp = union(ivs)
        if bp > mejor_bp:
            mejor, mejor_bp = s, bp
    L = largo.get(mejor, 0)
    breadth = round(100.0 * mejor_bp / L, 2) if L else float("nan")
    filas.append((m, len(ctgs[m]), mejor_bp, round(ident[m], 1),
                  mejor, cepa.get(mejor, "?"), breadth))

for f in sorted(filas, key=lambda r: -r[2]):
    print("\t".join(str(x) for x in f))

print()
print("== por muestra y cepa (para ver coinfecciones A+B) ==")
print("muestra\tcepa\tbp_union\trefs")
agg = collections.defaultdict(lambda: [0, 0])
for m, refs in iv.items():
    por_cepa = collections.defaultdict(list)
    for s, ivs in refs.items():
        por_cepa[cepa.get(s, "?")].append(union(ivs))
    for c, bps in por_cepa.items():
        agg[(m, c)] = [max(bps), len(bps)]
for (m, c), (bp, k) in sorted(agg.items(), key=lambda t: (t[0][0], -t[1][0])):
    print(f"{m}\t{c}\t{bp}\t{k}")
