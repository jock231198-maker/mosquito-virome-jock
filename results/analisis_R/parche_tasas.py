#!/usr/bin/env python3
"""
Dos cosas en figuras.qmd:

  table-rate  La tabla de tasas de infeccion incluye taxones que son elementos
              endogenos. Chaq-like, el primero por MLE, tiene un EVE al 0.882
              de fraccion en AaegL5 y el ARN pequeno no le encuentra una sola
              libreria replicativa. Una tasa de infeccion sobre un elemento
              endogeno no significa nada: lo lleva el genoma de todos los
              mosquitos. Se anade el veredicto del nivel B como columna y una
              segunda tabla restringida a lo que replica, que es la que se cita.

  g-kraken    El pie ya llevaba str_wrap pero a width = 125, que no cabe en
              200 mm y se salia por la derecha. Baja a 105.

El veredicto se une POR vOTU y no por nombre: el ARN largo escribe
"Hattula totivirus 3-like (59% aa)" y el pequeno "afin a Hattula totivirus 3",
que no casan, y ademas el largo distingue cuatro Hattula donde el pequeno ve
una. La vOTU es la unica clave comun.

Si no existe datos_srna_b_v1 el documento sigue renderizando: la columna no se
anade y la segunda tabla avisa en vez de fallar.

Tiene que estar en el mismo directorio que figuras.qmd.
No escribe nada si cualquier aserto falla. Deja figuras.qmd.bak3.

    python3 parche_tasas.py
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

# ------------------------------------------------------- 1. pie de kraken
VIEJO_W = "         width = 125)) +"
NUEVO_W = "         width = 105)) +"
n = texto.count(VIEJO_W)
assert n == 1, f"width del pie de kraken: esperaba 1 coincidencia, encontre {n}"
texto = texto.replace(VIEJO_W, NUEVO_W)

# ------------------------------------------------------------ 2. table-rate
VIEJO = '''prev |>
  filter(grupo == GRUPOS_VIRUS[["artropodo"]]) |>
  arrange(desc(pools)) |>
  rowwise() |>
  mutate(t = list(tasa_infeccion(pools, N_MUESTRAS, K_POOL))) |>
  ungroup() |>
  select(especie_fig, t) |>
  unnest(t) |>
  transmute(Species = especie_fig,
            `Positive pools` = paste0(pools_pos, "/", N_MUESTRAS),
            `MIR (‰)` = round(mir, 1),
            `MLE (‰)` = if_else(pools_pos == N_MUESTRAS, NA_real_, round(mle, 1)),
            `95% CI` = if_else(pools_pos == N_MUESTRAS, "all pools positive",
                                paste0(round(mle_inf, 1), "–", round(mle_sup, 1)))) |>
  kable(caption = paste0("Individual infection rate per 1,000 females (pools of ", K_POOL,
                         "). MLE after Chiang and Reeves with a transformed Clopper-Pearson CI; not estimable when all pools are positive."))'''

NUEVO = '''# El veredicto del ARN pequeno se une POR vOTU, no por nombre: el ARN largo
# escribe "Hattula totivirus 3-like (59% aa)" y el pequeno "afin a Hattula
# totivirus 3", que no casan, y ademas el largo distingue cuatro Hattula donde
# el pequeno ve una. nombre_virus() trabaja por vOTU, asi que el mapa sale de
# ahi y la agregacion se hace despues.
DATOS_B    <- "datos_srna_b_v1"
HAY_NIVELB <- file.exists(file.path(DATOS_B, "nivelB_consenso.tsv"))

veredicto <- NULL
if (HAY_NIVELB) {
  nb <- readr::read_tsv(file.path(DATOS_B, "nivelB_consenso.tsv"),
                        show_col_types = FALSE)
  veredicto <- votus |>
    filter(es_virus == "si") |>
    mutate(especie_fig = nombre_virus(pick(everything()))) |>
    select(votu, especie_fig) |>
    inner_join(select(nb, votu, consenso), by = "votu") |>
    group_by(especie_fig) |>
    summarise(replica = sum(consenso %in% c("replica", "MIXTO_replica_y_EVE")),
              eve     = sum(consenso == "EVE_NIRVS"),
              .groups = "drop") |>
    mutate(`Small RNA` = case_when(
      replica > 0 & eve > 0 ~ "Replicating + EVE",
      replica > 0           ~ "Replicating",
      eve > 0               ~ "EVE only",
      TRUE                  ~ "No evidence"))
}

tasas <- prev |>
  filter(grupo == GRUPOS_VIRUS[["artropodo"]]) |>
  arrange(desc(pools)) |>
  rowwise() |>
  mutate(t = list(tasa_infeccion(pools, N_MUESTRAS, K_POOL))) |>
  ungroup() |>
  select(especie_fig, t) |>
  unnest(t) |>
  transmute(Species = especie_fig,
            `Positive pools` = paste0(pools_pos, "/", N_MUESTRAS),
            `MIR (‰)` = round(mir, 1),
            `MLE (‰)` = if_else(pools_pos == N_MUESTRAS, NA_real_, round(mle, 1)),
            `95% CI` = if_else(pools_pos == N_MUESTRAS, "all pools positive",
                                paste0(round(mle_inf, 1), "–", round(mle_sup, 1))))

if (!is.null(veredicto)) {
  tasas <- tasas |>
    left_join(select(veredicto, Species = especie_fig, `Small RNA`, replica),
              by = "Species") |>
    mutate(`Small RNA` = coalesce(`Small RNA`, "Not tested"))
}

tasas |>
  select(-any_of("replica")) |>
  kable(caption = paste0(
    "Individual infection rate per 1,000 females (pools of ", K_POOL,
    "). MLE after Chiang and Reeves with a transformed Clopper-Pearson CI; not estimable when all ",
    "pools are positive. The last column is what the small RNA says about that species' vOTUs. A rate ",
    "computed over an endogenous element is not an infection rate: every mosquito carries the locus, ",
    "so the pool is positive whether or not anything is infected."))'''

n = texto.count(VIEJO)
assert n == 1, f"cuerpo de table-rate: esperaba 1 coincidencia, encontre {n}"
texto = texto.replace(VIEJO, NUEVO)

# ------------------------------------------- 3. segunda tabla, restringida
ANCLA = '''pools are positive. The last column is what the small RNA says about that species' vOTUs. A rate ",
    "computed over an endogenous element is not an infection rate: every mosquito carries the locus, ",
    "so the pool is positive whether or not anything is infected."))
```'''

EXTRA = ANCLA + '''

```{r table-rate-replica, eval = HAY_VOTUS}
if (!is.null(veredicto)) {
  tasas |>
    filter(replica > 0) |>
    select(-replica) |>
    kable(caption = paste0(
      "The same rates, restricted to the species with small RNA evidence of replication in at least one ",
      "library. These are the figures to quote as infection rates; everything left out of this table is a ",
      "detection, not a demonstrated infection."))
} else {
  cat("Level B data not found in `", DATOS_B, "`: the restricted table is not built.\\n", sep = "")
}
```'''

n = texto.count(ANCLA)
assert n == 1, f"ancla de la segunda tabla: esperaba 1 coincidencia, encontre {n}"
texto = texto.replace(ANCLA, EXTRA)

shutil.copy2(QMD, QMD + ".bak3")
with open(QMD, "w", encoding="utf-8", newline=nl) as fh:
    fh.write(texto)

print("  pie de kraken: width 125 -> 105")
print("  table-rate: columna Small RNA anadida")
print("  tabla nueva: table-rate-replica")
print(f"  finales de linea: {'CRLF' if crlf else 'LF'} (preservados)")
print(f"  copia de seguridad: {QMD}.bak3")
print("ok")
