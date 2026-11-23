#!/usr/bin/env python3
"""
Quita "the provider" de los textos visibles de las figuras P y deja el nombre
del pipeline. De paso corrige el pie de P2, que mencionaba dos referencias
marcadas como artefacto que el filtro de >=70 % ya deja fuera, y lo sustituye
por la limitacion que si aplica.

Tiene que estar en el mismo directorio que srna_proveedor.qmd.
No escribe nada si cualquier aserto falla.

    python3 parche_P_nombres.py
"""

import os
import shutil
import sys

QMD = "srna_proveedor.qmd"
PIPE = "smallRNA-Metavir"

if not os.path.exists(QMD):
    sys.exit(f"ERROR: no encuentro {QMD} en {os.getcwd()}")

with open(QMD, "rb") as fh:
    crudo = fh.read()
crlf = b"\r\n" in crudo
nl = "\r\n" if crlf else "\n"
texto = crudo.decode("utf-8").replace("\r\n", "\n")

R = []

# --- titulo del documento
R.append((
    'title: "MERIDA — small RNA: the provider\'s own analysis"',
    'title: "MERIDA — small RNA: de novo virome analysis (smallRNA-Metavir)"',
))

# --- comentario del orden de librerias
R.append((
    "# Orden de librerias: por la fraccion de piRNA que mide EL PROVEEDOR, no la\n"
    "# nuestra. Son ventanas distintas (24-30 frente a 24-35), asi que el orden\n"
    "# deberia parecerse al de S1 pero no tiene por que ser identico.",
    "# Orden de librerias: por la fraccion de piRNA que mide smallRNA-Metavir, no\n"
    "# la nuestra. Son ventanas distintas (24-30 frente a 24-35), asi que el orden\n"
    "# deberia parecerse al de S1 pero no tiene por que ser identico.",
))

# --- P1, comentario del panel B
R.append((
    "# --- B. composicion siRNA / piRNA SEGUN EL PROVEEDOR.\n"
    "# Es una medida independiente de la de S1, calculada sobre otras ventanas\n"
    "# (20-23 y 24-30 frente a 18-23 y 24-35) y por otro pipeline.",
    "# --- B. composicion siRNA / piRNA en las ventanas de smallRNA-Metavir.\n"
    "# Es una medida independiente de la de S1, calculada sobre otras ventanas\n"
    "# (20-23 y 24-30 frente a 18-23 y 24-35) y por otro pipeline. El tramo\n"
    "# 'other' es lo que cae en 18-35 pero fuera de las dos bandas: 18-19 y\n"
    "# 31-35 nt. Y el miRNA de 22 nt cae DENTRO de la banda de siRNA.",
))

# --- P1, subtitulo del panel B
R.append((
    'labs(subtitle = "B · Size composition, as measured by the provider",',
    'labs(subtitle = "B · Size composition, in the pipeline\'s size windows",',
))

# --- P1, subtitulo general
R.append((
    '    subtitle = paste0("The provider\'s own pipeline: no host filter, bacterial filter, ",\n'
    '                      "generic mosquito viral reference.\\nLibraries ordered by their ",\n'
    '                      "piRNA fraction."),',
    '    subtitle = paste0("No host filter, bacterial filter, generic mosquito viral ",\n'
    '                      "reference.\\nLibraries ordered by their piRNA fraction."),',
))

# --- P2, encabezado de seccion
R.append((
    "## P2 · The provider's taxonomy",
    "## P2 · Taxonomy from the small RNA",
))

# --- P2, comentario del filtro
R.append((
    "# comparables. El filtro propio del proveedor es mas laxo: >=30 % y >=100\n"
    "# lecturas.",
    "# comparables. El filtro por defecto del pipeline es mas laxo: >=30 % y >=100\n"
    "# lecturas.",
))

# --- P2, titulo
R.append((
    'labs(title = "Viruses detected by the small RNA, on the provider\'s own reference set",',
    'labs(title = "Viruses detected by the small RNA, on the smallRNA-Metavir reference set",',
))

# --- P2, pie: la frase de los artefactos sobraba (el filtro de 70 % ya los
#     deja fuera) y faltaba la limitacion que si aplica
R.append((
    '       caption = paste0("Independent of the vOTU catalogue: the provider maps against a ",\n'
    '                        "generic mosquito viral reference set. Two references are flagged ",\n'
    '                        "by the provider as artefacts (≈50 % coverage in every library).")) +',
    '       caption = str_wrap(paste0("Independent of the vOTU catalogue: smallRNA-Metavir ",\n'
    '                        "maps against a generic mosquito viral reference set. This figure ",\n'
    '                        "shows detection, not interpretation: the pipeline applies no host ",\n'
    '                        "filter and no screen against AaegL5, so a circulating virus and an ",\n'
    '                        "endogenous element look the same here."), width = 120)) +',
))

# --- P3, subtitulo
R.append((
    '    subtitle = paste0("Classified by the provider\'s model from their small RNA ",\n'
    '                      "profile alone, with no reference.\\nThese are the de novo ",\n'
    '                      "discovery candidates."),',
    '    subtitle = paste0("Classified by the pipeline\'s model from the small RNA ",\n'
    '                      "profile alone, with no reference.\\nThese are the de novo ",\n'
    '                      "discovery candidates."),',
))

# --- P3, pie
R.append((
    '    caption = paste0("A hit means BLASTN or DIAMOND against the provider\'s databases. ",',
    '    caption = paste0("A hit means BLASTN or DIAMOND against the ",\n'
    '                     "smallRNA-Metavir databases. ",',
))

# --- str_wrap lo pide el pie nuevo de P2
if "library(stringr)" not in texto:
    R.append((
        "  library(scales); library(forcats); library(patchwork); library(knitr)",
        "  library(scales); library(forcats); library(patchwork); library(knitr)\n"
        "  library(stringr)",
    ))

for i, (viejo, nuevo) in enumerate(R, start=1):
    n = texto.count(viejo)
    assert n == 1, f"reemplazo {i}: esperaba 1 coincidencia, encontre {n}\n---\n{viejo[:140]}"
    texto = texto.replace(viejo, nuevo)

quedan = texto.lower().count("provider") + texto.lower().count("proveedor")
shutil.copy2(QMD, QMD + ".bak")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print(f"  {len(R)} reemplazos aplicados")
print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  quedan {quedan} menciones de provider/proveedor en el fichero")
print(f"  copia de seguridad: {QMD}.bak")
print("ok")
