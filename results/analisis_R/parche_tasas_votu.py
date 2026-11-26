#!/usr/bin/env python3
"""
La tabla restringida recalcula la tasa sobre las vOTUs replicativas y no sobre
la especie entera.

El problema: los 19/22 de Verdadero vienen de la prevalencia de la especie, que
no distingue si el pool dio positivo por el virus o por su copia endogena. La
tasa resultante es un limite superior.

El arreglo: contar los pools solo sobre las vOTUs que el nivel B llama
replicativas. Se hace sobre `cob`, que conserva la vOTU y `pct_cubierto`;
`pres` no sirve porque presencia_especies() ya agrega por especie. Mismo
umbral de 70 % que esa funcion, para que las dos cifras sean comparables.

La tabla pasa a llevar las dos columnas de pools, la de la especie y la de las
vOTUs replicativas, de modo que la diferencia entre ambas -- la contribucion
endogena -- se ve directamente.

Tiene que estar en el mismo directorio que figuras.qmd, despues de
parche_pie_tasas.py.
No escribe nada si el aserto falla. Deja figuras.qmd.bak5.

    python3 parche_tasas_votu.py
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

VIEJO = '''  tasas |>
    filter(replica > 0) |>
    select(-replica) |>
    kable(caption = paste0(
      "The same rates, restricted to the ", n_rep, " species with small RNA evidence of replication in at ",
      "least one library. The ", nrow(tasas) - n_rep, " species left out are not one group but two, and ",
      "they carry different weight: in ", n_eve, " the small RNA actively reports an endogenous profile, ",
      "while ", n_sin, " stay below 10 rpm in every library, so it cannot judge them at all. Only the ",
      "first is evidence against infection; the second is absence of data. And in the ", n_mix,
      " species labelled 'Replicating + EVE' the pool count still comes from the long RNA, which does not ",
      "distinguish which of the two forms was detected, so those rates are upper bounds."))'''

NUEVO = '''  # Recuento de pools sobre las vOTUs REPLICATIVAS y no sobre la especie
  # entera. Los 19/22 de Verdadero incluyen pools donde solo estaba la copia
  # endogena, y la especie no los separa. Se cuenta sobre `cob`, que conserva
  # la vOTU y pct_cubierto; `pres` no vale porque presencia_especies() ya
  # agrega por especie. Mismo umbral de 70 % que esa funcion.
  rep_votus <- votus |>
    filter(es_virus == "si") |>
    mutate(especie_fig = nombre_virus(pick(everything()))) |>
    select(votu, especie_fig) |>
    inner_join(select(nb, votu, consenso), by = "votu") |>
    filter(consenso %in% c("replica", "MIXTO_replica_y_EVE")) |>
    select(votu, especie_fig)

  prev_rep <- cob |>
    inner_join(rep_votus, by = "votu") |>
    filter(pct_cubierto >= 70) |>
    distinct(especie_fig, muestra) |>
    count(especie_fig, name = "pools_rep")

  tasas |>
    filter(replica > 0) |>
    select(Species, `Positive pools`, `Small RNA`) |>
    left_join(rename(prev_rep, Species = especie_fig), by = "Species") |>
    mutate(pools_rep = coalesce(pools_rep, 0L)) |>
    rowwise() |>
    mutate(t = list(tasa_infeccion(pools_rep, N_MUESTRAS, K_POOL))) |>
    ungroup() |>
    unnest(t) |>
    transmute(Species,
              `Pools (species)` = `Positive pools`,
              `Pools (replicating vOTUs)` = paste0(pools_pos, "/", N_MUESTRAS),
              `MIR (‰)` = round(mir, 1),
              `MLE (‰)` = round(mle, 1),
              `95% CI` = paste0(round(mle_inf, 1), "–", round(mle_sup, 1)),
              `Small RNA`) |>
    kable(caption = paste0(
      "The ", n_rep, " species with small RNA evidence of replication, with the rate recomputed over their ",
      "replicating vOTUs alone. The first pool column is the species-level count that feeds the table ",
      "above; the second counts only the pools where a vOTU the small RNA calls replicating reaches 70 % ",
      "coverage, the same threshold. The gap between the two columns is the endogenous contribution, which ",
      "a species-level count cannot separate, and the rates here are computed from the second. The ",
      nrow(tasas) - n_rep, " species left out are not one group but two, and they carry different weight: ",
      "in ", n_eve, " the small RNA actively reports an endogenous profile, while ", n_sin, " stay below ",
      "10 rpm in every library, so it cannot judge them at all. Only the first is evidence against ",
      "infection; the second is absence of data."))'''

n = texto.count(VIEJO)
assert n == 1, (f"cuerpo de la tabla restringida: esperaba 1 coincidencia, encontre {n}. "
                f"Corriste antes parche_pie_tasas.py?")
texto = texto.replace(VIEJO, NUEVO)

shutil.copy2(QMD, QMD + ".bak5")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print("  la tabla restringida recalcula la tasa sobre las vOTUs replicativas")
print("  y muestra las dos columnas de pools para que se vea la diferencia")
print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  copia de seguridad: {QMD}.bak5")
print("ok")
