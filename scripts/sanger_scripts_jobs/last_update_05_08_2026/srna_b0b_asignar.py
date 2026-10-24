#!/usr/bin/env python3
"""
srna_b0b_asignar.py — asigna lecturas de ARN pequeno al pool de RNA total cuya version de un virus
coinciden mejor (B0b: intercambio M41/M48 con variantes virales).

Entrada por stdin: SAM de bowtie1 con TODAS las alineaciones de cada lectura (-a -v 2), contra los
contigs virales de varios pools. Cada contig se llama  <pool>__<virus>__<contig original>.

Una lectura es INFORMATIVA si:
  - alinea (<= 2 errores) con algun contig de TODOS los pools que tienen ese virus
    -> la region esta ensamblada en todos; asi un ensamblaje incompleto no atrae lecturas
  - el minimo de errores se alcanza en UN solo pool -> esa lectura distingue pools
Se asigna a ese pool.

Uso:
  bowtie ... -S | python3 srna_b0b_asignar.py <libreria> <pools_por_virus.tsv> >> asignaciones.tsv
  pools_por_virus.tsv: virus <TAB> pool   (los pools con contigs de cada virus)
Salida (stdout, sin cabecera): libreria  virus  pool  lecturas_asignadas
Por stderr: lecturas alineadas, informativas y descartadas.
"""
import sys, re
from collections import defaultdict

lib, tabla = sys.argv[1], sys.argv[2]
pools_de = defaultdict(set)
for l in open(tabla):
    if l.strip():
        v, p = l.rstrip('\n').split('\t')[:2]
        pools_de[v].add(p)

nm_re = re.compile(r'\tNM:i:(\d+)')
cuenta = defaultdict(int)
stats = defaultdict(int)

def cerrar(alins):
    """alins: lista de (pool, virus, nm) de una lectura."""
    if not alins:
        return
    stats['alineadas'] += 1
    por_virus = defaultdict(dict)                 # virus -> pool -> minimo nm
    for p, v, nm in alins:
        d = por_virus[v]
        d[p] = min(nm, d.get(p, 99))
    if len(por_virus) > 1:                         # cae en dos virus distintos: ambigua
        stats['varios_virus'] += 1
        return
    v, d = next(iter(por_virus.items()))
    if set(d) != pools_de[v]:
        stats['region_no_en_todos'] += 1
        return
    mejor = min(d.values())
    ganadores = [p for p, nm in d.items() if nm == mejor]
    if len(ganadores) != 1:
        stats['empate'] += 1
        return
    stats['informativas'] += 1
    cuenta[(v, ganadores[0])] += 1

actual, alins = None, []
cerradas = set()
for linea in sys.stdin:
    if linea.startswith('@'):
        continue
    f = linea.split('\t', 3)
    q, flag, ref = f[0], int(f[1]), f[2]
    if flag & 4 or ref == '*':
        continue
    if q != actual:
        cerrar(alins)
        if actual is not None:
            cerradas.add(actual)
        if q in cerradas:
            stats['no_contiguas'] += 1          # alineaciones de una lectura separadas: no deberia pasar
        actual, alins = q, []
    pool, virus = ref.split('__')[:2]
    m = nm_re.search(linea)
    alins.append((pool, virus, int(m.group(1)) if m else 0))
cerrar(alins)

for (v, p), n in sorted(cuenta.items()):
    print(f'{lib}\t{v}\t{p}\t{n}')
print(f'  {lib}: alineadas {stats["alineadas"]}, informativas {stats["informativas"]}, '
      f'empate {stats["empate"]}, region no ensamblada en todos {stats["region_no_en_todos"]}, '
      f'en varios virus {stats["varios_virus"]}', file=sys.stderr)
if stats['no_contiguas']:
    print(f'  !! {lib}: {stats["no_contiguas"]} lecturas con alineaciones no contiguas en el SAM; '
          f'el conteo puede estar sesgado (correr con -p 1)', file=sys.stderr)
