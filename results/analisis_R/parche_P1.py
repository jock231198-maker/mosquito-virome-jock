#!/usr/bin/env python3
"""
Anade el conteo de contigs con hit viral encima de las barras del panel C de P1.

Tiene que estar en el mismo directorio que srna_proveedor.qmd.
No escribe nada si el aserto falla.

    python3 parche_P1.py
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

VIEJO = """pC <- ggplot(ctg, aes(lib, n, fill = clase)) +
  geom_col(width = 0.78, colour = SUPERFICIE, linewidth = 0.2) +
  scale_fill_manual(values = setNames(c(PAL[["azul"]], PAL[["violeta"]], GRIS), CT)) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) +
"""

NUEVO = """# El numero sobre la barra es el conteo de contigs con hit viral. En la pila
# es una franja de pocos pixeles y no hay forma de leerlo. El total sale de
# `ctg` y no de `mu` para que la etiqueta quede exactamente en lo alto de lo
# que se dibuja, aunque alguna columna del proveedor venga vacia.
tot_ctg <- ctg |>
  group_by(lib) |>
  summarise(total = sum(n, na.rm = TRUE), .groups = "drop") |>
  left_join(select(mu, lib, virales = contigs_virales), by = "lib")

pC <- ggplot(ctg, aes(lib, n, fill = clase)) +
  geom_col(width = 0.78, colour = SUPERFICIE, linewidth = 0.2) +
  geom_text(data = tot_ctg, aes(lib, total, label = virales),
            inherit.aes = FALSE, vjust = -0.7, size = 2.1,
            colour = PAL[["azul"]], fontface = "bold") +
  scale_fill_manual(values = setNames(c(PAL[["azul"]], PAL[["violeta"]], GRIS), CT)) +
  scale_y_continuous(expand = expansion(c(0, 0.12))) +
"""

assert texto.count(VIEJO) == 1, \
    f"esperaba 1 bloque pC sin parchear, encontre {texto.count(VIEJO)}"
assert "tot_ctg" not in texto, "ya estaba parcheado: tot_ctg existe"

texto = texto.replace(VIEJO, NUEVO)

shutil.copy2(QMD, QMD + ".bak")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  copia de seguridad: {QMD}.bak")
print("ok")
