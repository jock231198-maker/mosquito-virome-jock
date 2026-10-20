#!/usr/bin/env python3
"""
srna_b0_frecuencias.py — compara librerias de ARN pequeno y de RNA total por sus frecuencias
alelicas en el genoma mitocondrial (cada pool de mosquitos tiene su propia mezcla de haplotipos).

Uso:
  samtools mpileup ... bam1 bam2 ... | python3 srna_b0_frecuencias.py <nombres.txt> <salida_dir> \
        [--min-dp 30] [--min-af 0.05] [--parejas M2=M2_S2,M11=M11_S3,...,M41=M48_S5]

  --parejas: pareja esperada de cada libreria de ARN pequeno. Las de control calibran la prueba;
             la hipotesis (M41=M48_S5) se juzga contra ellas.

  nombres.txt: una linea por BAM, en el MISMO orden que en el mpileup:  <nombre> <tipo>
               tipo = srna | rnaseq

Salidas:
  frecuencias_sitios_variables.tsv   sitio x muestra: frecuencia del alelo no-referencia (NA si poca profundidad)
  distancias.tsv                      todas las parejas: distancia media |diferencia de frecuencia| y n de sitios
  emparejamiento.tsv                  para cada libreria de ARN pequeno, los pools de RNA total ordenados por distancia
"""
import sys, os, re, itertools
from collections import defaultdict

args = sys.argv[1:]
if len(args) < 2:
    sys.exit(__doc__)
def opt(flag, default):
    return type(default)(args[args.index(flag) + 1]) if flag in args else default
MIN_DP = opt('--min-dp', 30)
MIN_AF = opt('--min-af', 0.05)
PAREJAS = dict(x.split('=') for x in opt('--parejas', '').split(',') if '=' in x)
nombres, tipos = [], []
for l in open(args[0]):
    if l.strip():
        n, t = l.split()[:2]
        nombres.append(n); tipos.append(t)
out = args[1]
os.makedirs(out, exist_ok=True)
ns = len(nombres)

indel = re.compile(r'[+-](\d+)')
def contar(bases, ref):
    """Cuenta A C G T en una columna de mpileup (. y , = referencia; ignora indels, ^X, $, *)."""
    c = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    i, n = 0, len(bases)
    while i < n:
        ch = bases[i]
        if ch == '^':
            i += 2; continue
        if ch in '+-':
            m = indel.match(bases, i)
            i = m.end() + int(m.group(1)); continue
        if ch in '.,':
            if ref in c: c[ref] += 1
        else:
            u = ch.upper()
            if u in c: c[u] += 1
        i += 1
    return c

# frecuencia del alelo no-referencia por sitio y muestra
af = {}            # pos -> [af o None]
for linea in sys.stdin:
    f = linea.rstrip('\n').split('\t')
    pos, ref = int(f[1]), f[2].upper()
    fila = []
    for k in range(ns):
        bases = f[4 + 3 * k] if len(f) > 4 + 3 * k else ''
        c = contar(bases, ref)
        dp = sum(c.values())
        fila.append(None if dp < MIN_DP or ref not in c else 1 - c[ref] / dp)
    # sitio variable: alguna muestra con frecuencia >= MIN_AF y otra (o la misma) con <= 1-MIN_AF
    vals = [v for v in fila if v is not None]
    if len(vals) >= 2 and max(vals) >= MIN_AF and min(vals) <= 1 - MIN_AF and (max(vals) - min(vals)) >= MIN_AF:
        af[pos] = fila

with open(os.path.join(out, 'frecuencias_sitios_variables.tsv'), 'w') as fo:
    fo.write('pos\t' + '\t'.join(nombres) + '\n')
    for p in sorted(af):
        fo.write(f'{p}\t' + '\t'.join('NA' if v is None else f'{v:.4f}' for v in af[p]) + '\n')

def dist(i, j):
    d = [abs(af[p][i] - af[p][j]) for p in af if af[p][i] is not None and af[p][j] is not None]
    return (sum(d) / len(d), len(d)) if d else (None, 0)

D = {}
with open(os.path.join(out, 'distancias.tsv'), 'w') as fo:
    fo.write('muestra_1\ttipo_1\tmuestra_2\ttipo_2\tdistancia\tn_sitios\n')
    for i, j in itertools.combinations(range(ns), 2):
        d, n = dist(i, j)
        D[(i, j)] = D[(j, i)] = (d, n)
        fo.write(f'{nombres[i]}\t{tipos[i]}\t{nombres[j]}\t{tipos[j]}\t{"NA" if d is None else f"{d:.4f}"}\t{n}\n')

srna = [i for i in range(ns) if tipos[i] == 'srna']
rna = [i for i in range(ns) if tipos[i] == 'rnaseq']
print(f'Sitios variables (profundidad >= {MIN_DP}, frecuencia >= {MIN_AF}): {len(af)}')
print()
print('Libreria de ARN pequeno -> pools de RNA total, del mas parecido al menos (distancia, n sitios)')
with open(os.path.join(out, 'emparejamiento.tsv'), 'w') as fo:
    fo.write('srna\trango\trnaseq\tdistancia\tn_sitios\n')
    for i in srna:
        orden = sorted((D[(i, j)][0] if D[(i, j)][0] is not None else 9, j) for j in rna)
        txt = []
        for r, (d, j) in enumerate(orden, 1):
            n = D[(i, j)][1]
            fo.write(f'{nombres[i]}\t{r}\t{nombres[j]}\t{"NA" if d == 9 else f"{d:.4f}"}\t{n}\n')
            txt.append(f'{nombres[j]} ({"NA" if d == 9 else f"{d:.3f}"}, {n})')
        print(f'  {nombres[i]:10s} ' + '  >  '.join(txt))

# ------------------------------------------------------------------ calibracion
if PAREJAS:
    idx = {n: i for i, n in enumerate(nombres)}
    print()
    print('Calibracion: distancia a la pareja esperada frente al mejor pool alternativo')
    print(f'  {"srna":6s} {"pareja":9s} {"rango":>5s} {"d_pareja":>9s} {"d_otro":>8s} {"margen":>7s}  otro mas parecido')
    filas = []
    with open(os.path.join(out, 'calibracion.tsv'), 'w') as fo:
        fo.write('srna\tpareja_esperada\trango\td_pareja\tmejor_otro\td_mejor_otro\tmargen\n')
        for s, r in PAREJAS.items():
            if s not in idx or r not in idx:
                continue
            i, j = idx[s], idx[r]
            dp = D[(i, j)][0]
            otros = sorted((D[(i, k)][0], nombres[k]) for k in rna if k != j and D[(i, k)][0] is not None)
            rango = 1 + sum(1 for d, _ in otros if dp is not None and d < dp)
            d_o, n_o = otros[0] if otros else (None, '-')
            m = None if dp is None or d_o is None else d_o - dp
            filas.append((s, r, rango, dp, m))
            fo.write(f'{s}\t{r}\t{rango}\t{dp:.4f}\t{n_o}\t{d_o:.4f}\t{m:.4f}\n')
            print(f'  {s:6s} {r:9s} {rango:5d} {dp:9.3f} {d_o:8.3f} {m:7.3f}  {n_o}')
    # resumen: controles = todas las parejas menos la ultima indicada (la hipotesis va al final)
    if len(filas) >= 2:
        ctrl, hip = filas[:-1], filas[-1]
        ok_ctrl = all(f[2] == 1 for f in ctrl)
        dmax = max(f[3] for f in ctrl); mmin = min(f[4] for f in ctrl)
        print()
        print(f'  Controles en primer lugar: {sum(f[2] == 1 for f in ctrl)}/{len(ctrl)}; '
              f'distancia a su pareja hasta {dmax:.3f}; margen minimo {mmin:.3f}')
        s, r, rango, dp, m = hip
        print(f'  Hipotesis {s} = {r}: rango {rango}, distancia {dp:.3f}, margen {m:.3f}')
        if not ok_ctrl:
            print('  -> La prueba NO esta calibrada: algun control no reconoce a su propio pool. No concluir.')
        elif rango == 1 and dp <= dmax * 1.25 and m > 0:
            fuerte = m >= mmin
            print(f'  -> Compatible con el intercambio ({"margen como el de los controles" if fuerte else "pero con margen menor que los controles: revisar"}).')
        else:
            print('  -> No apoya el intercambio.')
