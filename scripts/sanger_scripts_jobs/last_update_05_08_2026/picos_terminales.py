#!/usr/bin/env python3
"""
Los picos terminales de 100,000x: de donde salen y que cifra reportar.

EL PROBLEMA
  Verdadero, Chaq-like y Guadeloupe tienen profundidad media inflada por un
  pico en un extremo del contig. Los CV son 3.16, 3.48 y 1.86 -- los tres mas
  altos del catalogo, en los tres virus menos dudosos. Mientras no se resuelva
  no se pueden dar cifras de profundidad media de esos tres.

QUE HACE
  1. Localiza el pico de cada vOTU y dice si es terminal.
  2. Da la cifra que SI se puede reportar: mediana y media recortada, que no
     dependen de dos ventanas de 50 bp.
  3. Extrae la secuencia del pico a un FASTA, para pasarle barrnap y ver si es
     rRNA -- la hipotesis principal, porque las librerias no se depletaron.

USO
  python3 -I picos_terminales.py perfiles.tsv contigs.fasta ids.txt SALIDA/

  perfiles   qc_control_min500/perfil_cobertura/perfiles_cobertura.tsv
             columnas: muestra, votu, inicio, fin, profundidad
  contigs    contig_annot_min500/contigs_union_MERIDA.fasta
  ids        un votu por linea (p.ej. /tmp/tres.ids)
  SALIDA     se crean picos.fasta, extremos.fasta y picos_resumen.tsv
"""

import collections
import os
import statistics
import sys

if len(sys.argv) < 5:
    sys.exit(__doc__)

PERF, FASTA, IDS, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
os.makedirs(OUT, exist_ok=True)

MARGEN = 200      # bp desde un extremo para llamar "terminal" a una ventana
FACTOR = 20       # una ventana por encima de FACTOR x la mediana es un pico
FLANCO = 150      # bp a cada lado del pico que se extraen

quiero = {ln.strip() for ln in open(IDS) if ln.strip()}
print(f"vOTUs pedidas: {len(quiero)}")

# ------------------------------------------------- 1. secuencias y longitudes
seq = {}
nombre, trozos = None, []
with open(FASTA) as fh:
    for ln in fh:
        if ln.startswith(">"):
            if nombre in quiero:
                seq[nombre] = "".join(trozos)
            nombre, trozos = ln[1:].split()[0], []
        elif nombre in quiero:
            trozos.append(ln.strip())
if nombre in quiero:
    seq[nombre] = "".join(trozos)
print(f"secuencias encontradas: {len(seq)}")
faltan = quiero - set(seq)
if faltan:
    print(f"AVISO: sin secuencia para {len(faltan)}, p.ej. {sorted(faltan)[:3]}")

# --------------------------------------------------------- 2. los perfiles
# (muestra, votu) -> [(inicio, fin, profundidad), ...]
perf = collections.defaultdict(list)
with open(PERF) as fh:
    cab = fh.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(cab)}
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < len(cab) or f[ix["votu"]] not in quiero:
            continue
        try:
            perf[(f[ix["muestra"]], f[ix["votu"]])].append(
                (int(f[ix["inicio"]]), int(f[ix["fin"]]), float(f[ix["profundidad"]])))
        except ValueError:
            continue
print(f"perfiles (muestra x vOTU): {len(perf)}")

# ------------------------------------- 3. por vOTU, la muestra de mas media
mejor = {}
for (m, v), vs in perf.items():
    media = statistics.fmean(p for _, _, p in vs)
    if v not in mejor or media > mejor[v][1]:
        mejor[v] = (m, media, vs)

filas, picos_fa, extremos_fa = [], [], []
for v, (m, media, vs) in sorted(mejor.items(), key=lambda t: -t[1][1]):
    vs.sort()
    largo = len(seq.get(v, "")) or max(f for _, f, _ in vs)
    profs = [p for _, _, p in vs]
    mediana = statistics.median(profs)
    corte = mediana * FACTOR

    picos = [(i, f, p) for i, f, p in vs if p > corte]
    limpio = [p for _, _, p in vs if p <= corte]
    media_rec = statistics.fmean(limpio) if limpio else float("nan")

    i0, f0, p0 = max(vs, key=lambda t: t[2])
    desde_5 = i0
    desde_3 = largo - f0
    terminal = min(desde_5, desde_3) <= MARGEN
    extremo = "5'" if desde_5 <= desde_3 else "3'"

    filas.append((v, m, largo, len(vs), round(media, 1), round(mediana, 1),
                  round(media_rec, 1), round(media / media_rec, 1) if media_rec == media_rec and media_rec else "",
                  f"{i0}-{f0}", round(p0, 0), len(picos),
                  "si" if terminal else "no", extremo if terminal else "",
                  min(desde_5, desde_3)))

    s = seq.get(v, "")
    if s and picos:
        a = max(0, min(i for i, _, _ in picos) - 1 - FLANCO)
        b = min(len(s), max(f for _, f, _ in picos) + FLANCO)
        picos_fa.append(f">{v}__pico_{a + 1}_{b} largo={largo} prof={int(p0)}\n{s[a:b]}\n")
    if s:
        extremos_fa.append(f">{v}__5p_300 largo={largo}\n{s[:300]}\n")
        extremos_fa.append(f">{v}__3p_300 largo={largo}\n{s[-300:]}\n")

CAB = ["votu", "muestra", "largo", "ventanas", "media", "mediana", "media_recortada",
       "factor_inflado", "ventana_pico", "prof_pico", "n_ventanas_pico",
       "terminal", "extremo", "bp_al_extremo"]

print("\n== por vOTU, en su muestra de mayor profundidad media ==")
print("   media_recortada excluye las ventanas por encima de "
      f"{FACTOR}x la mediana; factor_inflado = media / media_recortada")
anchos = [52, 9, 7, 9, 11, 9, 16, 15, 14, 11, 8, 9, 8, 14]
print("".join(c[:w - 1].ljust(w) for c, w in zip(CAB, anchos)))
for f in filas:
    print("".join(str(x)[:w - 1].ljust(w) for x, w in zip(f, anchos)))

with open(os.path.join(OUT, "picos_resumen.tsv"), "w") as fh:
    fh.write("\t".join(CAB) + "\n")
    for f in filas:
        fh.write("\t".join(str(x) for x in f) + "\n")
with open(os.path.join(OUT, "picos.fasta"), "w") as fh:
    fh.writelines(picos_fa)
with open(os.path.join(OUT, "extremos.fasta"), "w") as fh:
    fh.writelines(extremos_fa)

# ------------------------------------------------- 4. composicion del pico
print("\n== composicion de la secuencia del pico ==")
print("   homopolimero = la racha mas larga de una sola base")
print(f"{'votu':<52}{'largo':<8}{'GC%':<8}{'homopol':<10}{'base':<7}{'N%':<7}")
for linea in picos_fa:
    cab_fa, s = linea.split("\n")[0], linea.split("\n")[1]
    v = cab_fa[1:].split("__pico")[0]
    s = s.upper()
    if not s:
        continue
    gc = 100.0 * sum(c in "GC" for c in s) / len(s)
    n = 100.0 * s.count("N") / len(s)
    mejor_racha, base, racha, prev = 0, "", 0, ""
    for c in s:
        racha = racha + 1 if c == prev else 1
        prev = c
        if racha > mejor_racha:
            mejor_racha, base = racha, c
    print(f"{v[:51]:<52}{len(s):<8}{gc:<8.1f}{mejor_racha:<10}{base:<7}{n:<7.1f}")

print(f"\nescritos {OUT}/picos_resumen.tsv, {OUT}/picos.fasta y {OUT}/extremos.fasta")
print("siguiente paso: barrnap --kingdom bac sobre picos.fasta y extremos.fasta")
