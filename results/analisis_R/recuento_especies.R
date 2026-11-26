#!/usr/bin/env Rscript
# Recuento de virus distintos, resolviendo los nombres contra las secuencias.
#
# El problema de contar por nombre (especie_fig) es que falla en las dos
# direcciones a la vez, y las dos estan medidas ya por blastn todos contra
# todos:
#
#   - un nombre repartido en varias secuencias  -> fragmentacion del ensamblado
#     (o segmentos de un virus segmentado bajo un solo nombre)
#   - una secuencia bajo varios nombres         -> nombres de mas, casi siempre
#     por el sufijo (NN% aa), que cambia entre fragmentos del mismo genoma
#
# Asi que ni los nombres ni los grupos de secuencia sirven solos. La unidad
# correcta es la componente conexa de la union de las dos relaciones: dos vOTUs
# son el mismo virus si comparten grupo de secuencia O comparten nombre. Se
# resuelve por propagacion de etiquetas, que converge sin paquetes extra.
#
# Encima van las fusiones que blastn no puede ver: segmentos de un mismo virus
# sin similitud mutua. Esas van a mano y declaradas en MANUAL.
#
# Correr desde results/analisis_R con datos_srna_b_v1/grupos.tsv en su sitio.
#
#     Rscript recuento_especies.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
})

DATOS <- "datos_v4"
source("R/leer.R")
source("R/tema.R")

votus <- leer_votus()
cob   <- leer_cobertura()
gr    <- read_tsv("datos_srna_b_v1/grupos.tsv", show_col_types = FALSE)

# grupo_virus es columna en TABLA_FINAL y puede ser tambien funcion en R/.
# No lo doy por hecho: si la funcion existe la uso, si no leo la columna.
grupo_de <- if (exists("grupo_virus", mode = "function")) {
  function(d) grupo_virus(d)
} else {
  function(d) d$grupo_virus
}

v <- votus |>
  filter(es_virus == "si") |>
  mutate(especie_fig = nombre_virus(pick(everything())),
         grupo_v     = grupo_de(pick(everything()))) |>
  left_join(select(gr, votu, grupo), by = "votu") |>
  mutate(grupo = coalesce(grupo, paste0("solo_", votu)))

# presencia con el mismo criterio que presencia_especies().
# La longitud la pone cob y no v: son la misma cifra (el largo del contig) y
# tomarla de los dos lados deja longitud.x / longitud.y y ninguna longitud.
pres_v <- cob |>
  filter(pct_cubierto >= 70) |>
  inner_join(select(v, votu, especie_fig, grupo, grupo_v),
             by = "votu")

art <- pres_v |> filter(grupo_v == GRUPOS_VIRUS[["artropodo"]])

# las que se quedan fuera por no llegar a 70 % en ningun pool
fuera <- v |>
  filter(grupo_v == GRUPOS_VIRUS[["artropodo"]],
         !votu %in% art$votu)

cat("\n=============================================================\n")
cat("vOTUs de artropodo, en el catalogo:  ", nrow(filter(v, grupo_v == GRUPOS_VIRUS[["artropodo"]])), "\n")
cat("  de ellas, presentes a >=70 %:      ", n_distinct(art$votu), "\n")
cat("  sin llegar a 70 % en ningun pool:  ", nrow(fuera), "\n")
cat("especie_fig, que es lo que cuenta hoy:", n_distinct(art$especie_fig), "\n")
cat("grupos de secuencia:                  ", n_distinct(art$grupo), "\n")
cat("=============================================================\n")

if (nrow(fuera)) {
  cat("\n== las que no llegan a 70 % ==\n")
  fuera |>
    transmute(votu, longitud, especie_fig) |>
    arrange(desc(longitud)) |>
    as.data.frame() |>
    print(row.names = FALSE)
}

# --- fusiones que blastn no puede ver: segmentos sin similitud mutua.
# Verdadero virus: RNA1 (RdRp, hallmark), RNA2 y RNA3 tipo chaq. Los tres
# grupos de secuencia estan confirmados; el RNA3 lleva el nombre Chaq-like.
MANUAL <- list(
  c("Verdadero virus", "Chaq-like virus")
)

mapa <- art |> distinct(especie_fig) |> mutate(fusion = especie_fig)
for (g in MANUAL) {
  hit <- mapa$especie_fig %in% g
  presentes <- g[g %in% mapa$especie_fig]
  if (sum(hit) > 1) {
    mapa$fusion[hit] <- g[1]
    cat("\nfusion manual aplicada: ", paste(presentes, collapse = " + "), "\n")
  } else {
    cat("\nfusion manual NO aplicada, falta alguno de los dos nombres.",
        "\n  pedidos:   ", paste(g, collapse = " + "),
        "\n  presentes: ", if (length(presentes)) paste(presentes, collapse = " + ") else "ninguno",
        "\n")
  }
}

art2 <- art |> left_join(mapa, by = "especie_fig")

# --- componentes conexas sobre (grupo de secuencia, nombre fusionado).
# Propagacion de etiquetas: cada vOTU arranca con su propia etiqueta y se
# queda con la minima de su grupo y la minima de su nombre, hasta que nada
# cambia. Eso une a la vez la fragmentacion y los nombres repetidos, que es
# justo lo que ninguna de las dos relaciones hace sola.
ent <- art2 |> distinct(votu, grupo, fusion, longitud)
ent$comp <- seq_len(nrow(ent))
vueltas <- 0
repeat {
  antes <- ent$comp
  ent <- ent |> group_by(grupo)  |> mutate(comp = min(comp)) |> ungroup()
  ent <- ent |> group_by(fusion) |> mutate(comp = min(comp)) |> ungroup()
  vueltas <- vueltas + 1
  if (identical(ent$comp, antes) || vueltas > 50) break
}

n_fig   <- n_distinct(art$especie_fig)
n_grp   <- n_distinct(art$grupo)
n_final <- n_distinct(ent$comp)

# nombre representativo de cada entidad: se prefiere un nombre sin sufijo de
# identidad, es decir una especie nombrada, y a falta de eso el nombre con mas
# vOTUs detras
rep <- ent |>
  count(comp, fusion, name = "n") |>
  mutate(nombrada = !grepl("% aa\\)$", fusion)) |>
  arrange(comp, desc(nombrada), desc(n), fusion) |>
  group_by(comp) |>
  summarise(Entidad  = first(fusion),
            nombres  = n_distinct(fusion),
            nombrada = first(nombrada),
            .groups  = "drop")

pools <- art2 |>
  left_join(select(ent, votu, comp), by = "votu") |>
  group_by(comp) |>
  summarise(pools = n_distinct(muestra), .groups = "drop")

# OJO con kb_suma: hay que sumar UNA longitud por grupo de secuencia, no una
# por vOTU. Dentro de un grupo las vOTUs son copias de lo mismo (lo que cd-hit
# no fusiono por quedarse bajo el 85 % de cobertura mutua), asi que sumarlas
# cuenta el mismo trozo de genoma varias veces. Verdadero lo deja claro: 7
# vOTUs en 3 grupos daban 12.6 kb para un partitivirus de ~5 kb en total.
# De cada grupo se toma el contig mas largo, que es el representante.
por_grupo <- ent |>
  group_by(comp, grupo) |>
  summarise(len = max(longitud), .groups = "drop")

res <- por_grupo |>
  group_by(comp) |>
  summarise(grupos   = n(),
            kb_suma  = round(sum(len) / 1000, 1),
            kb_mayor = round(max(len) / 1000, 1),
            piezas_bajo_1kb = sum(len < 1000),
            .groups  = "drop") |>
  left_join(count(ent, comp, name = "votus"), by = "comp") |>
  left_join(rep,   by = "comp") |>
  left_join(pools, by = "comp")

cat("\n=============================================================\n")
cat("RECUENTO\n")
cat("  por nombre (lo que dice la figura hoy) ....", n_fig, "\n")
cat("  por grupo de secuencia ....................", n_grp, "\n")
cat("  virus distintos (componentes conexas) .....", n_final, "\n")
cat("    de ellos, especies nombradas ............", sum(res$nombrada), "\n")
cat("    parientes (-like / -related) ............", sum(!res$nombrada), "\n")
cat(sprintf("  (convergio en %d vueltas)\n", vueltas))
cat("=============================================================\n")

cat("\n== NOMBRES DE MAS: entidades que juntan varios nombres ==\n")
cat("   cada nombre extra es una especie que la figura cuenta dos veces\n")
mezcla <- ent |>
  left_join(select(res, comp, Entidad), by = "comp") |>
  distinct(comp, Entidad, fusion) |>
  group_by(comp, Entidad) |>
  filter(n() > 1) |>
  summarise(nombres = n(),
            cuales  = paste(sort(unique(fusion)), collapse = "  +  "),
            .groups = "drop") |>
  arrange(desc(nombres))
if (nrow(mezcla)) {
  cat(sprintf("   total de nombres de mas: %d\n\n", sum(mezcla$nombres) - nrow(mezcla)))
  print(as.data.frame(select(mezcla, nombres, cuales)), row.names = FALSE)
} else cat("  ninguno\n")

cat("\n== FRAGMENTACION: entidades repartidas en varias secuencias ==\n")
cat("   kb_suma es la suma de todas las piezas, kb_mayor la pieza mas larga.\n")
cat("   Si kb_suma se parece al genoma esperado, son pedazos de uno solo;\n")
cat("   si lo pasa de largo, puede que sean virus distintos con el mismo hit.\n\n")
frag <- res |>
  filter(grupos > 1) |>
  transmute(Entidad = substr(Entidad, 1, 42),
            grupos, votus, pools, kb_suma, kb_mayor, bajo_1kb = piezas_bajo_1kb) |>
  arrange(desc(grupos))
if (nrow(frag)) print(as.data.frame(frag), row.names = FALSE) else cat("  ninguna\n")

cat("\n== las", n_final, "entidades finales ==\n")
res |>
  transmute(Entidad = substr(Entidad, 1, 44),
            nombres, grupos, votus, pools, kb_suma) |>
  arrange(desc(pools), Entidad) |>
  as.data.frame() |>
  print(row.names = FALSE)

write_tsv(select(res, Entidad, nombres, grupos, votus, pools,
                 kb_suma, kb_mayor, piezas_bajo_1kb),
          "recuento_entidades.tsv")
write_tsv(left_join(select(ent, votu, grupo, fusion, longitud, comp),
                    select(res, comp, Entidad), by = "comp"),
          "recuento_votu_entidad.tsv")
cat("\nescritos recuento_entidades.tsv y recuento_votu_entidad.tsv\n")
