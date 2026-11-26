#!/usr/bin/env python3
"""
Extiende el pie de la figura 08 con la limitacion de rRNA, medida.

POR QUE NO UN UMBRAL DE ABUNDANCIA
  La Wolbachia estaba al 16 % de las lecturas no-hospedero en M20. Ningun
  umbral razonable la quita: es de las senales mas fuertes de la figura. El
  fallo no fue de abundancia sino de naturaleza de la secuencia -- rRNA
  conservado que Kraken2 asigna a un genero concreto, y eso pasa a cualquier
  nivel. Un umbral darIa una falsa sensacion de control y ademas tiraria
  taxones de poca abundancia que si pueden ser reales.

QUE HACE
  El chunk YA quitaba Homo sapiens y Plasmodium por ser rRNA residual de
  mosquito asignado al eucariota mas cercano. Esto extiende ese mismo
  razonamiento a las bacterias y le pone la verificacion detras: la senal
  bacteriana mas fuerte de la figura se comprobo a nivel de contig y resulto
  ser rRNA 23S conservado.

  La figura NO se toca: se quedan las mismas barras, los mismos taxones y el
  mismo orden. Solo crece el pie, y con el la altura.

Tiene que estar en el mismo directorio que figuras.qmd.
No escribe nada si algun aserto falla. Deja figuras.qmd.bak7.

    python3 parche_08_rrna.py
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

CAMBIOS = [
    (
        "el pie, con la verificacion de Wolbachia",
        '''                                   paste(etiqueta_muestra(pocos), collapse = ", ")) else "", "."),
         width = 105)) +''',
        '''                                   paste(etiqueta_muestra(pocos), collapse = ", ")) else "", ". ",
         "The libraries were not rRNA-depleted, so these assignments rest largely on ribosomal sequence, which is ",
         "conserved across taxa: the strongest bacterial signal here (Wolbachia, 16 % of non-host reads in one pool) ",
         "was verified at contig level and proved to be conserved 23S rRNA, not that organism. Read this panel as the ",
         "relative composition of ribosomal signal; a specific taxon needs contig-level verification before it is ",
         "claimed."),
         width = 105)) +''',
    ),
    (
        "altura del chunk, para que quepa el pie",
        "```{r g-kraken, fig.height = 5.5}",
        "```{r g-kraken, fig.height = 6.2}",
    ),
    (
        "altura del fichero guardado",
        'guardar(p8, "08_kraken2", alto = 140)',
        'guardar(p8, "08_kraken2", alto = 158)',
    ),
]

for rotulo, viejo, nuevo in CAMBIOS:
    n = texto.count(viejo)
    assert n == 1, (f"{rotulo}: esperaba 1 coincidencia, encontre {n}. "
                    f"El chunk g-kraken no esta como esperaba.")

for rotulo, viejo, nuevo in CAMBIOS:
    texto = texto.replace(viejo, nuevo)
    print(f"  ok  {rotulo}")

shutil.copy2(QMD, QMD + ".bak7")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  copia de seguridad: {QMD}.bak7")
print("ok")
