#!/usr/bin/env python3
"""
srna_b0b_resumen.py — resume asignaciones.tsv: para cada libreria de ARN pequeno, que fraccion de sus
lecturas informativas apunta a cada pool, y calibra con las parejas de control.

Uso: python3 srna_b0b_resumen.py <asignaciones.tsv> <salida_dir> --parejas M2=M2_S2,...,M41=M48_S5
"""
import sys, os
from collections import defaultdict

asig, out = sys.argv[1], sys.argv[2]
args = sys.argv[3:]
PAREJAS = {}
if '--parejas' in args:
    PAREJAS = dict(x.split('=') for x in args[args.index('--parejas') + 1].split(',') if '=' in x)

n = defaultdict(lambda: defaultdict(int))        # lib -> pool -> lecturas (todos los virus)
nv = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))   # lib -> virus -> pool -> lecturas
pools, libs = set(), []
for l in open(asig):
    if not l.strip() or l.startswith('libreria'):
        continue
    lib, v, p, c = l.rstrip('\n').split('\t')
    n[lib][p] += int(c); nv[lib][v][p] += int(c); pools.add(p)
    if lib not in libs:
        libs.append(lib)
pools = sorted(pools)

with open(os.path.join(out, 'fracciones.tsv'), 'w') as fo:
    fo.write('libreria\tvirus\tpool\tlecturas\tfraccion\n')
    for lib in libs:
        for v in sorted(nv[lib]):
            tot = sum(nv[lib][v].values())
            for p in pools:
                c = nv[lib][v][p]
                fo.write(f'{lib}\t{v}\t{p}\t{c}\t{c / tot:.4f}\n')
        tot = sum(n[lib].values())
        for p in pools:
            fo.write(f'{lib}\tTODOS\t{p}\t{n[lib][p]}\t{(n[lib][p] / tot if tot else 0):.4f}\n')

print('Fraccion de lecturas informativas asignadas a cada pool (todos los virus)')
print('  ' + f'{"libreria":9s}{"n":>8s}  ' + ''.join(f'{p:>9s}' for p in pools))
frac = {}
for lib in libs:
    tot = sum(n[lib].values())
    frac[lib] = {p: (n[lib][p] / tot if tot else 0) for p in pools}
    print('  ' + f'{lib:9s}{tot:8d}  ' + ''.join(f'{frac[lib][p]:9.2f}' for p in pools))

print()
print('Por virus (pool con mas lecturas, fraccion, n informativas)')
for lib in libs:
    partes = []
    for v in sorted(nv[lib]):
        tot = sum(nv[lib][v].values())
        p, c = max(nv[lib][v].items(), key=lambda x: x[1])
        partes.append(f'{v}: {p} {c / tot:.2f} (n={tot})')
    print(f'  {lib:9s} ' + ' | '.join(partes))

if PAREJAS:
    print()
    print('Calibracion: fraccion en la pareja esperada frente al mejor pool alternativo')
    print(f'  {"srna":6s} {"pareja":9s} {"n":>7s} {"f_pareja":>9s} {"f_otro":>8s}  otro')
    filas = []
    with open(os.path.join(out, 'calibracion.tsv'), 'w') as fo:
        fo.write('srna\tpareja\tn_informativas\tf_pareja\tmejor_otro\tf_mejor_otro\n')
        for s, r in PAREJAS.items():
            if s not in frac or r not in pools:
                print(f'  {s:6s} {r:9s}  sin lecturas informativas'); continue
            tot = sum(n[s].values())
            otros = sorted(((f, p) for p, f in frac[s].items() if p != r), reverse=True)
            fo_, po = otros[0] if otros else (0, '-')
            filas.append((s, r, tot, frac[s][r], fo_))
            fo.write(f'{s}\t{r}\t{tot}\t{frac[s][r]:.4f}\t{po}\t{fo_:.4f}\n')
            print(f'  {s:6s} {r:9s} {tot:7d} {frac[s][r]:9.2f} {fo_:8.2f}  {po}')
    if len(filas) >= 2:
        ctrl, hip = filas[:-1], filas[-1]
        ok = [f for f in ctrl if f[3] > f[4] and f[2] >= 50]
        print()
        print(f'  Controles que reconocen su propio pool (con >= 50 lecturas informativas): {len(ok)}/{len(ctrl)}')
        s, r, tot, fp, fo_ = hip
        print(f'  Hipotesis {s} = {r}: {fp:.2f} de {tot} lecturas informativas (mejor alternativo {fo_:.2f})')
        if len(ok) < len(ctrl) - 1 or len(ok) < 2:
            print('  -> La prueba NO esta calibrada. No concluir.')
        elif tot >= 50 and fp > fo_:
            print('  -> Apoya el intercambio: la libreria "M41" lleva las variantes virales de M48_S5.')
        else:
            print('  -> No apoya el intercambio.')
