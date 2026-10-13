# cobertura.R — lectores y utilidades del informe de cobertura (cobertura.qmd).
# Se carga despues de leer.R y tema.R y no modifica nada de ellos.

if (!exists("COBERTURA")) COBERTURA <- "datos_cobertura_v1"

# Tres tipos de vOTU para comparar la forma de la cobertura. Los fagos van con los
# parientes: son dos vOTUs y un cuarto color no pasaria la prueba de daltonismo
# en un diagrama de dispersion.
TIPOS <- c("Identified virus", "Related virus (-like)", "EVE")
COL_TIPO <- setNames(unname(PAL[c("azul", "aqua", "naranja")]), TIPOS)

tipo_votu <- function(votus) {
  case_when(
    votus$clase_final == "identificado"                   ~ TIPOS[1],
    votus$es_virus == "si"                                ~ TIPOS[2],
    votus$clase_final %in% c("EVE_genoma", "EVE_candidato") ~ TIPOS[3],
    TRUE                                                  ~ NA_character_
  )
}

# Nombre de las vOTUs EVE: no tienen nombre_virus (no son virus), asi que se
# toma el del hit viral (blastn; si no, nr) y se marca como EVE. Si varias EVE
# comparten nombre se numeran de la mas larga a la mas corta.
nombre_eve <- function(votus) {
  hit <- coalesce(na_if(votus$especie, ""), na_if(votus$nr_nombre, ""), "unnamed")
  tibble(votu = votus$votu, hit = hit, longitud = votus$longitud) |>
    group_by(hit) |>
    mutate(n = n(), i = rank(-longitud, ties.method = "first")) |>
    ungroup() |>
    mutate(nombre = paste0("EVE · ", hit, if_else(n > 1, paste0(" #", i), ""))) |>
    pull(nombre)
}

leer_perfiles <- function(dir = COBERTURA) {
  f <- file.path(dir, "perfiles_cobertura.tsv.gz")
  if (!file.exists(f)) stop("No existe ", f, call. = FALSE)
  read_tsv(f, show_col_types = FALSE, progress = FALSE, col_types = "cciid") |>
    comprobar(c("muestra", "votu", "inicio", "fin", "profundidad"), f)
}

hay_perfiles <- function(dir = COBERTURA) file.exists(file.path(dir, "perfiles_cobertura.tsv.gz"))

# Las figuras de cobertura van a su propia carpeta. guardar() escribe en
# "figuras/<nombre>", asi que se le pasa una ruta relativa que sale de ahi.
guardar_cob <- function(p, nombre, ...) {
  dir.create("figuras_cobertura", showWarnings = FALSE)
  guardar(p, file.path("..", "figuras_cobertura", nombre), ...)
}

# Uniformidad de un perfil: coeficiente de variacion de la profundidad entre
# ventanas y fraccion de ventanas con al menos 1x. Un genoma cubierto de punta a
# punta tiene un CV moderado y casi todas las ventanas >= 1x; una cobertura a
# trozos tiene un CV alto y muchas ventanas vacias.
uniformidad <- function(perfiles) {
  perfiles |>
    group_by(muestra, votu) |>
    summarise(prof_media = mean(profundidad),
              cv = if (mean(profundidad) > 0) sd(profundidad) / mean(profundidad) else NA_real_,
              frac_ventanas_1x = mean(profundidad >= 1),
              .groups = "drop")
}
