#!/usr/bin/env python3
"""
Quita las menciones a "the provider" que quedan en los textos visibles de las
figuras P y deja el nombre del pipeline. Ajustado a la version del fichero del
23 nov, con str_wrap ya puesto en los pies de P2 y P3.

Corrige ademas el pie de P2: la frase de las dos referencias marcadas como
artefacto sobra, porque con el filtro de >=70 % nunca llegan a la figura, y en
su lugar va la limitacion que si aplica.

    python3 parche_P_nombres2.py

No escribe nada si cualquier aserto falla. Deja srna_proveedor.qmd.bak.
"""

import os
import shutil
import sys

QMD = "srna_proveedor.qmd"

if not os.path.exists(QMD):
    sys.exit(f"ERROR: no encuentro {QMD} en {os.getcwd()}")

with open(QMD, "rb") as fh:
    crudo = fh.read()
crlf = b"\r\n" in crudo
nl = "\r\n" if crlf else "\n"
texto = crudo.decode("utf-8").replace("\r\n", "\n")

R = [
    # 1. titulo del documento
    ('title: "MERIDA — small RNA: the provider\'s own analysis"',
     'title: "MERIDA — small RNA: de novo virome analysis (smallRNA-Metavir)"'),

    # 2. P1, subtitulo general
    ('    subtitle = paste0("The provider\'s own pipeline: no host filter, bacterial filter, ",\n'
     '                      "generic mosquito viral reference.\\nLibraries ordered by their ",\n'
     '                      "piRNA fraction."),',
     '    subtitle = paste0("No host filter, bacterial filter, generic mosquito viral ",\n'
     '                      "reference.\\nLibraries ordered by their piRNA fraction."),'),

    # 3. P2, encabezado de seccion
    ("## P2 · The provider's taxonomy",
     "## P2 · Taxonomy from the small RNA"),

    # 4. P2, titulo de la figura
    ('  labs(title = "Viruses detected by the small RNA, on the provider\'s own reference set",',
     '  labs(title = "Viruses detected by the small RNA, on the smallRNA-Metavir reference set",'),

    # 5. P2, pie: fuera la frase de los artefactos, dentro la limitacion real
    ('         "Independent of the vOTU catalogue: the provider maps against a ",\n'
     '         "generic mosquito viral reference set. Two references are flagged ",\n'
     '         "by the provider as artefacts (≈50 % coverage in every library)."),',
     '         "Independent of the vOTU catalogue: smallRNA-Metavir maps against a ",\n'
     '         "generic mosquito viral reference set and infers replicating-vs-endogenous ",\n'
     '         "from the small RNA profile alone, which is what the dots show. It applies ",\n'
     '         "no host filter and no screen against AaegL5, so it cannot say whether a ",\n'
     '         "sequence also sits in the mosquito genome."),'),

    # 6. P3, subtitulo
    ('    subtitle = paste0("Classified by the provider\'s model from their small RNA ",',
     '    subtitle = paste0("Classified by the pipeline\'s model from the small RNA ",'),

    # 7. P3, pie
    ('      "A hit means BLASTN or DIAMOND against the provider\'s databases. ",',
     '      "A hit means BLASTN or DIAMOND against the smallRNA-Metavir databases. ",'),
]

for i, (viejo, nuevo) in enumerate(R, start=1):
    n = texto.count(viejo)
    assert n == 1, (f"reemplazo {i}: esperaba 1 coincidencia, encontre {n}\n"
                    f"--- texto buscado ---\n{viejo}")
    texto = texto.replace(viejo, nuevo)

quedan = texto.lower().count("provider")
shutil.copy2(QMD, QMD + ".bak")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print(f"  {len(R)} reemplazos aplicados")
print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  quedan {quedan} menciones de 'provider' (deberia ser 0)")
print(f"  copia de seguridad: {QMD}.bak")
print("ok")
