#!/usr/bin/env python3
"""
Anade el mapa de calor de los 50 taxones mas abundantes de Kraken2/Bracken.

POR QUE
  La figura 08 pasa a suplementario, y una barra apilada con 7 taxones mas
  "Other" esconde casi todo lo que Kraken2 reporta. En suplementario, donde no
  compite por espacio con las figuras virales, cabe el panel completo.

COMO SE CONSTRUYE
  Sobre `br`, el mismo objeto que alimenta la figura 08: solo microbiota, ya
  renormalizado, sin los pools por debajo de MIN_LECTURAS. Y reutiliza
  `levels(kr$etiqueta)`, de modo que el orden de pools es EL MISMO que en la
  barra apilada y las dos figuras se pueden leer una contra otra.

  Depende de que el chunk g-kraken se haya evaluado antes. En un render de
  arriba a abajo eso se cumple; si alguna vez se ejecuta suelto, hay que correr
  g-kraken primero.

Tiene que estar en el mismo directorio que figuras.qmd, despues de
parche_08_rrna.py.
No escribe nada si el aserto falla. Deja figuras.qmd.bak8.

    python3 parche_08b_heatmap.py
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

ANCLA = '''guardar(p8, "08_kraken2", alto = 158)
p8
```

# Data used'''

NUEVO = '''guardar(p8, "08_kraken2", alto = 158)
p8
```

```{r g-kraken-heat, fig.height = 10.2}
# Mapa de calor de los 50 taxones mas abundantes. Va a suplementario junto con
# la figura 08. Se construye sobre `br` -- el mismo objeto renormalizado del
# chunk anterior -- y reutiliza levels(kr$etiqueta), asi que el orden de pools
# es identico al de la barra apilada y las dos se leen una contra otra.
TOP_HEAT <- 50

top_heat <- br |>
  group_by(name) |>
  summarise(media = sum(fraction_total_reads) / N_MUESTRAS, .groups = "drop") |>
  slice_max(media, n = TOP_HEAT) |>
  arrange(desc(media))

heat <- br |>
  filter(name %in% top_heat$name) |>
  select(muestra, name, frac = fraction_total_reads) |>
  complete(muestra, name, fill = list(frac = 0)) |>
  mutate(frac = if_else(frac <= 0, NA_real_, frac),
         etiqueta = factor(etiqueta_muestra(muestra), levels = levels(kr$etiqueta)),
         taxon = factor(name, levels = rev(top_heat$name)))

p8b <- ggplot(heat, aes(etiqueta, taxon, fill = frac)) +
  geom_tile(colour = SUPERFICIE, linewidth = 0.3) +
  scale_fill_viridis_c(trans = "log10", direction = -1,
                       labels = label_percent(accuracy = 0.1),
                       na.value = "#EDEDE8",
                       name = "fraction of\\nclassified\\nreads") +
  labs(title = paste0("Kraken2/Bracken: the ", TOP_HEAT, " most abundant taxa per pool"),
       subtitle = "Same filtering, re-normalisation and pool order as the stacked-bar figure; log colour scale",
       x = NULL, y = NULL,
       caption = str_wrap(paste0(
         "Grey cells: taxon not reported in that pool. Bacterial and eukaryotic taxa only; Homo sapiens, ",
         "Plasmodium and viral taxa were removed before re-normalising, as in the stacked-bar figure. The same ",
         "caveat applies, and is why both panels are supplementary: the libraries were not rRNA-depleted, so these ",
         "assignments rest largely on conserved ribosomal sequence. Contig-level verification of the strongest ",
         "bacterial signal (Wolbachia) showed it to be conserved 23S rRNA and not that organism. The taxa at the ",
         "top of this panel are, in large part, those whose ribosomal sequence is best represented in the reference ",
         "database. No individual taxon here should be claimed without contig-level verification."),
         width = 150)) +
  theme(axis.text.y = element_text(size = 5.5, face = "italic"),
        axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        legend.key.height = grid::unit(14, "mm"))
guardar(p8b, "08b_kraken2_heatmap", alto = 260)
p8b
```

# Data used'''

n = texto.count(ANCLA)
assert n == 1, (f"cierre del chunk g-kraken: esperaba 1 coincidencia, encontre {n}. "
                f"Corriste antes parche_08_rrna.py?")
texto = texto.replace(ANCLA, NUEVO)

shutil.copy2(QMD, QMD + ".bak8")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print("  anadido el chunk g-kraken-heat tras g-kraken")
print("  escribe figuras/08b_kraken2_heatmap")
print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  copia de seguridad: {QMD}.bak8")
print("ok")
