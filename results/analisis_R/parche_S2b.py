#!/usr/bin/env python3
"""
Dos arreglos en S2:
  1. El eje x llegaba a 30 cuando solo hay 26 librerias. Cortes explicitos.
  2. Las especies con replicacion quedaban abajo. Con facet_grid el primer
     nivel es el panel de arriba, asi que se invierte el orden y el resultado
     queda donde se empieza a leer. Dentro de cada especie el eje discreto ya
     iba de mayor a menor RPM, asi que ahora los dos ordenes coinciden.

Tiene que estar en el mismo directorio que srna.qmd.
No escribe nada si cualquier aserto falla. Deja srna.qmd.bak.

    python3 parche_S2b.py
"""

import os
import shutil
import sys

QMD = "srna.qmd"

if not os.path.exists(QMD):
    sys.exit(f"ERROR: no encuentro {QMD} en {os.getcwd()}")

with open(QMD, "rb") as fh:
    crudo = fh.read()
crlf = b"\r\n" in crudo
nl = "\r\n" if crlf else "\n"
texto = crudo.decode("utf-8").replace("\r\n", "\n")

R = [
    # 1. comentario del orden
    ("# Orden. Especies: primero las que no tienen ninguna vOTU replicativa, luego\n"
     "# por RPM, de modo que las que replican quedan arriba. Dentro de cada especie,\n"
     "# por RPM.",
     "# Orden. Las especies con alguna vOTU replicativa van primero y, dentro de\n"
     "# ellas, por RPM descendente: con facet_grid el primer nivel es el panel de\n"
     "# arriba, asi que el resultado queda donde se empieza a leer. Dentro de cada\n"
     "# especie el eje discreto va al reves, de modo que el mayor RPM tambien queda\n"
     "# arriba y los dos ordenes coinciden."),

    # 2. el orden en si
    ('            rpm_esp = max(rpm_max), .groups = "drop") |>\n'
     '  arrange(rep_esp, rpm_esp)',
     '            rpm_esp = max(rpm_max), .groups = "drop") |>\n'
     '  arrange(desc(rep_esp), desc(rpm_esp))'),

    # 3. cortes del eje x: no hay 30 librerias
    ('  scale_x_continuous(expand = expansion(c(0, 0.18))) +',
     '  scale_x_continuous(breaks = c(0, 10, 20, n_lib),\n'
     '                     expand = expansion(c(0, 0.18))) +'),
]

for i, (viejo, nuevo) in enumerate(R, start=1):
    n = texto.count(viejo)
    assert n == 1, (f"reemplazo {i}: esperaba 1 coincidencia, encontre {n}\n"
                    f"--- texto buscado ---\n{viejo}")
    texto = texto.replace(viejo, nuevo)

shutil.copy2(QMD, QMD + ".bak")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print(f"  {len(R)} reemplazos aplicados")
print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  copia de seguridad: {QMD}.bak")
print("ok")
