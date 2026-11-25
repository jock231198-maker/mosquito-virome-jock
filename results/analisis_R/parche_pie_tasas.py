#!/usr/bin/env python3
"""
Arregla el pie de la tabla restringida. Decia que todo lo que queda fuera "es
una deteccion, no una infeccion demostrada", y eso mete en el mismo saco dos
cosas distintas: las especies donde el ARN pequeno dice activamente que es
endogeno, y aquellas donde no pudo ver nada porque no llegan a 10 rpm. La
primera es una afirmacion; la segunda es ausencia de datos.

Anade ademas el reparto en cifras, calculado del propio dato para que no se
quede viejo, y la advertencia de que las especies marcadas Replicating + EVE
llevan una tasa que es un limite superior: el recuento de pools viene del ARN
largo y no distingue cual de las dos formas se detecto.

Tiene que estar en el mismo directorio que figuras.qmd, despues de
parche_tasas.py.
No escribe nada si el aserto falla. Deja figuras.qmd.bak4.

    python3 parche_pie_tasas.py
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

VIEJO = '''if (!is.null(veredicto)) {
  tasas |>
    filter(replica > 0) |>
    select(-replica) |>
    kable(caption = paste0(
      "The same rates, restricted to the species with small RNA evidence of replication in at least one ",
      "library. These are the figures to quote as infection rates; everything left out of this table is a ",
      "detection, not a demonstrated infection."))
} else {'''

NUEVO = '''if (!is.null(veredicto)) {
  # El reparto sale del dato y no escrito a mano, para que no se quede viejo.
  n_rep <- sum(grepl("Replicating", tasas$`Small RNA`))
  n_eve <- sum(tasas$`Small RNA` == "EVE only")
  n_sin <- sum(tasas$`Small RNA` == "No evidence")
  n_mix <- sum(tasas$`Small RNA` == "Replicating + EVE")
  tasas |>
    filter(replica > 0) |>
    select(-replica) |>
    kable(caption = paste0(
      "The same rates, restricted to the ", n_rep, " species with small RNA evidence of replication in at ",
      "least one library. The ", nrow(tasas) - n_rep, " species left out are not one group but two, and ",
      "they carry different weight: in ", n_eve, " the small RNA actively reports an endogenous profile, ",
      "while ", n_sin, " stay below 10 rpm in every library, so it cannot judge them at all. Only the ",
      "first is evidence against infection; the second is absence of data. And in the ", n_mix,
      " species labelled 'Replicating + EVE' the pool count still comes from the long RNA, which does not ",
      "distinguish which of the two forms was detected, so those rates are upper bounds."))
} else {'''

n = texto.count(VIEJO)
assert n == 1, (f"pie de la tabla restringida: esperaba 1 coincidencia, encontre {n}. "
                f"Corriste antes parche_tasas.py?")
texto = texto.replace(VIEJO, NUEVO)

shutil.copy2(QMD, QMD + ".bak4")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print("  pie de la tabla restringida reescrito, con el reparto calculado del dato")
print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  copia de seguridad: {QMD}.bak4")
print("ok")
