#!/usr/bin/env python3
"""
Sustituye el chunk g-srna-virus de srna.qmd por la version de S2_nuevo.qmd,
y anade library(stringr) al setup si falta.

Los dos ficheros tienen que estar en el mismo directorio que srna.qmd.
No escribe nada si cualquiera de los asertos falla.

    python3 parche_S2.py
"""

import os
import shutil
import sys

QMD = "srna.qmd"
NUEVO = "S2_nuevo.qmd"

for f in (QMD, NUEVO):
    if not os.path.exists(f):
        sys.exit(f"ERROR: no encuentro {f} en {os.getcwd()}")

# --- leer preservando el final de linea original
with open(QMD, "rb") as fh:
    crudo = fh.read()
crlf = b"\r\n" in crudo
nl = "\r\n" if crlf else "\n"
texto = crudo.decode("utf-8")
lineas = texto.replace("\r\n", "\n").split("\n")

with open(NUEVO, encoding="utf-8") as fh:
    nuevas = fh.read().replace("\r\n", "\n").rstrip("\n").split("\n")

# --- localizar el chunk por sus delimitadores, no por su contenido
inicios = [i for i, l in enumerate(lineas)
           if l.strip().startswith("```{r g-srna-virus")]
assert len(inicios) == 1, f"esperaba 1 chunk g-srna-virus, encontre {len(inicios)}"
i = inicios[0]

cierres = [j for j in range(i + 1, len(lineas)) if lineas[j].strip() == "```"]
assert cierres, "no encuentro el cierre ``` del chunk"
j = cierres[0]

viejo_n = j - i + 1
assert 30 < viejo_n < 90, f"el chunk viejo mide {viejo_n} lineas, no cuadra"
assert any("slice_max(rpm_max" in l for l in lineas[i:j]), \
    "el chunk viejo no tiene el slice_max que vengo a quitar; ya estaba parcheado?"
assert any("guardar_srna(pS2" in l for l in lineas[i:j]), \
    "el chunk viejo no guarda pS2; no es el chunk que creo"

# --- stringr en el setup (lo pide str_wrap del pie)
if "library(stringr)" not in texto:
    cand = [k for k, l in enumerate(lineas) if "library(patchwork)" in l]
    assert len(cand) == 1, f"esperaba 1 linea con library(patchwork), hay {len(cand)}"
    k = cand[0]
    lineas[k] = lineas[k].rstrip() + "; library(stringr)"
    print(f"  setup: {lineas[k].strip()}")
else:
    print("  setup: stringr ya estaba")

lineas[i:j + 1] = nuevas

shutil.copy2(QMD, QMD + ".bak")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write("\n".join(lineas))

print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  chunk sustituido: {viejo_n} lineas -> {len(nuevas)}")
print(f"  copia de seguridad: {QMD}.bak")
print("ok")
