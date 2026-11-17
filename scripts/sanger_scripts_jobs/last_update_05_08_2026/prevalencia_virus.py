#!/usr/bin/env python3
"""
prevalencia_virus.py — prevalencia por especie SOLO sobre las vOTUs que el
arbitraje acepta como virus, a partir de TABLA_FINAL_votus_arb.tsv.

Por que existe: tabla_final.py escribe prevalencia_por_especie.tsv ANTES del
arbitraje y de la integracion, asi que reporta nombres de especie de vOTUs que
integrar_arbitro.py descarta despues (EVE_genoma, celular_otro, huesped). Por
eso esa tabla lista "Tick-borne encephalitis virus" en 20 de 22 muestras.

Agrupa por nombre_virus y une las muestras de mapeo de todas las vOTUs de cada
especie, en vez de quedarse con el maximo de una sola.

Uso:
  python3 prevalencia_virus.py TABLA_FINAL_votus_arb.tsv [SALIDA.tsv]
"""
import csv
import sys
from collections import defaultdict

VERDAD = {"1", "si", "sí", "sí", "true", "yes", "y", "v"}
N_MUESTRAS = 22


def truthy(x):
    return (x or "").strip().lower() in VERDAD


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    ent = sys.argv[1]
    sal = sys.argv[2] if len(sys.argv) > 2 else "prevalencia_virus.tsv"

    with open(ent, newline="") as fh:
        lector = csv.DictReader(fh, delimiter="\t")
        if "es_virus" not in (lector.fieldnames or []):
            sys.exit(f"ERROR: {ent} no tiene columna 'es_virus'. "
                     f"Columnas: {lector.fieldnames}")
        todas = list(lector)

    filas = [f for f in todas if truthy(f.get("es_virus"))]
    if not filas:
        vistos = sorted({(f.get("es_virus") or "").strip() for f in todas})
        sys.exit(f"ERROR: ninguna fila con es_virus verdadero. "
                 f"Valores presentes: {vistos}")

    mapeo = defaultdict(set)
    asm = defaultdict(int)
    votus = defaultdict(int)
    prof = defaultdict(float)
    grupo, clases = {}, defaultdict(set)

    for f in filas:
        nombre = (f.get("nombre_virus") or "").strip() or "(sin nombre)"
        votus[nombre] += 1
        grupo[nombre] = (f.get("grupo_virus") or "").strip()
        clases[nombre].add((f.get("clase_final") or "").strip())
        for m in (f.get("muestras_mapeo") or "").replace(";", ",").split(","):
            m = m.strip()
            if m:
                mapeo[nombre].add(m)
        for col, dest in (("n_muestras_asm", asm), ("profundidad_media", prof)):
            try:
                dest[nombre] = max(dest[nombre], float(f.get(col) or 0))
            except ValueError:
                pass

    orden = sorted(votus, key=lambda n: (-len(mapeo[n]), -asm[n], n))

    with open(sal, "w") as out:
        out.write("nombre_virus\tgrupo\tvOTUs\tn_mapeo\tn_asm\t"
                  "prof_max\tclase_final\tmuestras_mapeo\n")
        for n in orden:
            out.write(f"{n}\t{grupo[n]}\t{votus[n]}\t{len(mapeo[n])}\t"
                      f"{int(asm[n])}\t{prof[n]:.1f}\t"
                      f"{','.join(sorted(c for c in clases[n] if c))}\t"
                      f"{','.join(sorted(mapeo[n]))}\n")

    print(f"vOTUs con es_virus: {len(filas)} de {len(todas)}")
    print(f"especies: {len(votus)}")
    print(f"escrito: {sal}\n")
    for n in orden:
        nm, na = len(mapeo[n]), int(asm[n])
        aviso = "  <-- revisar: mapea mucho, ensambla poco" \
            if (nm >= 18 and na <= 3) else ""
        print(f"  {n:<45} {nm:>2}/{N_MUESTRAS} mapeo  "
              f"{na:>2}/{N_MUESTRAS} ensamblaje  {votus[n]} vOTUs{aviso}")


if __name__ == "__main__":
    main()
