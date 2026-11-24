#!/usr/bin/env python3
"""
Tres arreglos en figuras.qmd, de la revision del 24 nov:

  02_assembly  el subtitulo decia "contigs >= 1,000 bp" y los datos son de la
               rama de 500. Ahora el umbral se lee del propio dato, asi que no
               puede volver a quedarse viejo. Y el eje de longitudes gana un
               corte en 500, que antes empezaba en 1,000 con datos por debajo.

  06b          59 taxones en fig.height = 5: las etiquetas se pisaban y la
               figura era ilegible. Pasa a la misma maqueta que 06_prevalence
               (fig.height 10, alto 230). Corte del eje en 22, no en 20. Y pie,
               que no tenia: la discordancia no se interpreta sola.

  07_richness  pie, que tampoco tenia. La riqueza no esta rarefizada y los
               pools van de 10 a 75 M de pares crudos.

Tiene que estar en el mismo directorio que figuras.qmd.
No escribe nada si cualquier aserto falla. Deja figuras.qmd.bak.

    python3 parche_figuras_rev.py
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

R = [
    # ---------------------------------------------------------- 02_assembly
    ('                  subtitle = paste0(format(nrow(contigs), big.mark = ","), " contigs ≥ 1,000 bp from ",\n'
     '                                    n_distinct(contigs$ens), " assemblers, dereplicated per sample"))',
     '                  subtitle = paste0(format(nrow(contigs), big.mark = ","), " contigs ≥ ",\n'
     '                                    format(min(contigs$longitud), big.mark = ","), " bp from ",\n'
     '                                    n_distinct(contigs$ens), " assemblers, dereplicated per sample"))'),

    ('  scale_x_log10(labels = label_number(big.mark = ","), breaks = c(1e3, 3e3, 1e4, 3e4)) +',
     '  scale_x_log10(labels = label_number(big.mark = ","), breaks = c(500, 1e3, 3e3, 1e4, 3e4)) +'),

    # ----------------------------------------------------------------- 06b
    ("```{r g-assembly-vs-mapping, eval = HAY_VOTUS, fig.height = 5}",
     "```{r g-assembly-vs-mapping, eval = HAY_VOTUS, fig.height = 10}"),

    ('  scale_x_continuous(limits = c(0, N_MUESTRAS), breaks = seq(0, N_MUESTRAS, 5)) +\n'
     '  labs(title = "Assembly versus read mapping",\n'
     '       subtitle = "Pools in which each virus was assembled versus pools in which it is detected by read mapping",\n'
     '       x = paste0("pools (of ", N_MUESTRAS, ")"), y = NULL) +\n'
     '  theme(panel.grid.major.y = element_blank())\n'
     'guardar(p6b, "06b_assembly_vs_mapping", alto = 110)',

     '  scale_x_continuous(limits = c(0, N_MUESTRAS),\n'
     '                     breaks = sort(unique(c(seq(0, N_MUESTRAS, 5), N_MUESTRAS)))) +\n'
     '  labs(title = "Assembly versus read mapping",\n'
     '       subtitle = "Pools in which each virus was assembled versus pools in which it is detected by read mapping",\n'
     '       x = paste0("pools (of ", N_MUESTRAS, ")"), y = NULL,\n'
     '       caption = paste0(\n'
     '         "A virus detected by read mapping in many pools but assembled in only one or two is the signature of an\\n",\n'
     '         "endogenous element or a cellular hit: every mosquito carries the locus, so every pool recruits reads over\\n",\n'
     '         "it, but it only assembles where coverage happens to suffice. The screen against AaegL5 catches the same\\n",\n'
     '         "cases more directly, which is why this is a review signal and not a criterion of its own.")) +\n'
     '  theme(panel.grid.major.y = element_blank(),\n'
     '        axis.text.y = element_text(size = 6.5))\n'
     'guardar(p6b, "06b_assembly_vs_mapping", alto = 230)'),

    # ---------------------------------------------------------- 07_richness
    ('  labs(title = "Arthropod virus species detected per pool",\n'
     '       subtitle = "Colour: collection date (2024). Phages excluded", x = NULL, y = "species") +',

     '  labs(title = "Arthropod virus species detected per pool",\n'
     '       subtitle = "Colour: collection date (2024). Phages excluded", x = NULL, y = "species",\n'
     '       caption = paste0(\n'
     '         "Richness is not rarefied. The pools range from 10 to 75 million raw read pairs and differ in their\\n",\n'
     '         "non-host fraction, so part of the spread between pools is sequencing effort rather than biology.")) +'),
]

for i, (viejo, nuevo) in enumerate(R, start=1):
    n = texto.count(viejo)
    assert n == 1, (f"reemplazo {i}: esperaba 1 coincidencia, encontre {n}\n"
                    f"--- texto buscado ---\n{viejo}")
    texto = texto.replace(viejo, nuevo)

assert "1,000 bp from" not in texto, "quedo una mencion al umbral de 1,000 bp"

shutil.copy2(QMD, QMD + ".bak")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print(f"  {len(R)} reemplazos aplicados")
print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  copia de seguridad: {QMD}.bak")
print("ok")
