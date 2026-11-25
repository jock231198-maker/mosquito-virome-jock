#!/usr/bin/env python3
"""
Deja P1b solo con el titulo y la nota del orden de librerias: fuera el
subtitulo de las ventanas de tamano y fuera el pie entero.

Tiene que estar en el mismo directorio que srna_metavir_v1.qmd.
No escribe nada si el aserto falla. Deja srna_metavir_v1.qmd.bak.

    python3 parche_P1b.py
"""

import os
import shutil
import sys

QMD = "srna_metavir_v1.qmd"

if not os.path.exists(QMD):
    sys.exit(f"ERROR: no encuentro {QMD} en {os.getcwd()}")

with open(QMD, "rb") as fh:
    crudo = fh.read()
crlf = b"\r\n" in crudo
nl = "\r\n" if crlf else "\n"
texto = crudo.decode("utf-8").replace("\r\n", "\n")

VIEJO = (
    '  labs(title = "Size composition of the 18–35 nt reads",\n'
    '       subtitle = paste0("In the pipeline\'s own size windows, which differ from the ",\n'
    '                         "ones used in level B (18–23 and 24–35 nt).\\n", sub_orden),\n'
    '       x = NULL, y = "of the 18–35 nt reads", fill = NULL,\n'
    '       caption = str_wrap(paste0(\n'
    '         "\'Other\' is what falls inside 18–35 nt but outside both named bands: ",\n'
    '         "18–19 and 31–35 nt. Two cautions. The 22 nt miRNA class — the most ",\n'
    '         "abundant small RNA in a whole mosquito — sits inside the siRNA band, so ",\n'
    '         "that band is not pure siRNA. And no rRNA was removed, so every band still ",\n'
    '         "carries rRNA fragments."), width = 120)) +\n'
)

NUEVO = (
    '  labs(title = "Size composition of the 18–35 nt reads",\n'
    '       subtitle = sub_orden,\n'
    '       x = NULL, y = "of the 18–35 nt reads", fill = NULL) +\n'
)

n = texto.count(VIEJO)
assert n == 1, f"esperaba 1 bloque labs() de P1b sin parchear, encontre {n}"

texto = texto.replace(VIEJO, NUEVO)

shutil.copy2(QMD, QMD + ".bak")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  copia de seguridad: {QMD}.bak")
print("ok")
