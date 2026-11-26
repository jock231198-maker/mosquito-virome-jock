#!/usr/bin/env python3
"""
El pie de la tabla restringida decia que el hueco entre las dos columnas de
pools era la contribucion endogena. Con los datos delante no hay hueco: en las
seis especies la vOTU replicativa llega sola al total. Eso es un resultado, no
una ausencia, y el pie debe decirlo en vez de describir un hueco que no esta.

Tiene que estar en el mismo directorio que figuras.qmd, despues de
parche_tasas_votu.py.
No escribe nada si el aserto falla. Deja figuras.qmd.bak6.

    python3 parche_pie_coinciden.py
"""

import os
import shutil
import sys

QMD = "figuras.qmd"

if not os.path.exists(QMD):
    sys.exit(f"ERROR: no encuentro {QMD} en {os.getcwd()}")

with open(QMD, "rb") as fh:
    crudo = fh.read()
crlf = b"\r\n" in crudo
nl = "\r\n" if crlf else "\n"
texto = crudo.decode("utf-8").replace("\r\n", "\n")

VIEJO = '''      "above; the second counts only the pools where a vOTU the small RNA calls replicating reaches 70 % ",
      "coverage, the same threshold. The gap between the two columns is the endogenous contribution, which ",
      "a species-level count cannot separate, and the rates here are computed from the second. The ",'''

NUEVO = '''      "above; the second counts only the pools where a vOTU the small RNA calls replicating reaches 70 % ",
      "coverage, the same threshold. Where the two columns agree, every positive pool contains the ",
      "replicating sequence itself at that coverage, so the rate does not rest on the endogenous copies; ",
      "a gap between them would be the endogenous contribution, which a species-level count cannot ",
      "separate. The rates below are computed from the second column. One caveat survives either way: the ",
      "small RNA verdict is per vOTU across libraries, not per pool, so the rate still assumes that ",
      "detection in a pool means infection in that pool. The ",'''

n = texto.count(VIEJO)
assert n == 1, (f"pie de la tabla restringida: esperaba 1 coincidencia, encontre {n}. "
                f"Corriste antes parche_tasas_votu.py?")
texto = texto.replace(VIEJO, NUEVO)

shutil.copy2(QMD, QMD + ".bak6")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print("  pie reescrito: las columnas coinciden y eso es el resultado")
print("  anadida la salvedad de que el veredicto es por vOTU y no por pool")
print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  copia de seguridad: {QMD}.bak6")
print("ok")
