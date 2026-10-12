# leer.R — lectura de las tablas congeladas (datos_vN).
# Cada lector comprueba que estan las columnas que usan las figuras: si un
# script de la farm cambia su salida, el error sale aqui y no en una figura rara.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

if (!exists("DATOS")) DATOS <- "datos_v1"

ruta_datos <- function(f) {
  r <- file.path(DATOS, f)
  if (!file.exists(r)) stop("No existe ", r, call. = FALSE)
  r
}

comprobar <- function(d, cols, f) {
  faltan <- setdiff(cols, names(d))
  if (length(faltan)) stop(f, ": faltan columnas ", paste(faltan, collapse = ", "), call. = FALSE)
  d
}

leer_tsv <- function(f, cols = NULL, ...) {
  d <- read_tsv(ruta_datos(f), show_col_types = FALSE, progress = FALSE, comment = "#", ...)
  if (!is.null(cols)) comprobar(d, cols, f) else d
}

existe <- function(f) file.exists(file.path(DATOS, f))

# "M11_S3" -> "M11"
etiqueta_muestra <- function(x) sub("_S[0-9]+$", "", x)

ENSAMBLADORES <- c(
  trinity   = "Trinity",
  rnaviral  = "SPAdes rnaviral",
  rna       = "rnaSPAdes",
  megahit   = "MEGAHIT",
  metaviral = "SPAdes metaviral"
)

# id de contig "<muestra>__<ensamblador>__<resto>"
partir_id <- function(id) {
  p <- strsplit(id, "__", fixed = TRUE)
  tibble(
    muestra     = vapply(p, `[`, "", 1),
    ensamblador = vapply(p, `[`, "", 2)
  )
}

# ---------------------------------------------------------------- tablas

leer_metadata <- function(f = "metadata_MERIDA.tsv") {
  read_tsv(f, show_col_types = FALSE, col_types = cols(fecha = col_date(), .default = col_guess())) |>
    comprobar(c("muestra", "fecha", "individuos_pool"), f)
}

leer_lecturas <- function() {
  leer_tsv("read_accounting.tsv",
           c("sample", "raw_pairs", "trimmed_pairs", "host_pct", "unmapped_pairs"),
           na = c("", "NA")) |>
    rename(muestra = sample)
}

leer_contigs <- function() {
  leer_tsv("contigs_union.tsv", c("muestra", "ensamblador", "contig", "longitud"))
}

leer_checkv <- function() {
  d <- leer_tsv("checkv_viral.tsv", c("contig_id", "contig_length", "checkv_quality", "completeness"))
  bind_cols(d, partir_id(d$contig_id))
}

leer_bracken <- function(dir = "kraken2_nohost") {
  fs <- list.files(ruta_datos(dir), pattern = "\\.bracken$", full.names = TRUE)
  if (!length(fs)) stop("Sin ficheros .bracken en ", dir, call. = FALSE)
  bind_rows(lapply(fs, function(f) {
    read_tsv(f, show_col_types = FALSE, progress = FALSE) |>
      comprobar(c("name", "taxonomy_lvl", "new_est_reads", "fraction_total_reads"), basename(f)) |>
      mutate(muestra = sub("\\.bracken$", "", basename(f)))
  }))
}

# votu_coverage.tsv no tiene cabecera (salida de cov_worker.sh)
leer_cobertura <- function() {
  read_tsv(ruta_datos("votu_coverage.tsv"), show_col_types = FALSE, progress = FALSE,
           col_names = c("muestra", "votu", "longitud", "bases_cubiertas", "pct_cubierto", "profundidad"))
}

leer_votus <- function() {
  leer_tsv("TABLA_FINAL_votus_arb.tsv",
           c("votu", "longitud", "clase", "especie", "familia", "n_muestras_asm",
             "nr_veredicto", "nr_nombre", "frac_en_genoma", "clase_final", "es_virus"),
           na = c("", "NA"))
}

leer_veredicto_contigs <- function() {
  leer_tsv("veredicto_contigs.tsv", c("votu", "veredicto"), na = c("", "NA")) |>
    rename(contig = votu)
}

# ---------------------------------------------------------------- derivadas

# Nombre de cada vOTU viral para las figuras. Solo "identificado" lleva el nombre
# de la especie tal cual; en pariente_debil y divergente el nombre es el del mejor
# hit, no una identificacion, y se marca como "afin a". Asi una vOTU cuyo pariente
# mas cercano es Oropouche virus no se lee como deteccion de OROV.
# Desde integrar_arbitro.py v3 el nombre lo decide el pipeline (columna
# nombre_virus: en los parientes, el mejor hit de proteina en nr). Esto es el
# respaldo para tablas anteriores.
nombre_virus <- function(votus) {
  if ("nombre_virus" %in% names(votus)) return(coalesce(na_if(votus$nombre_virus, ""), votus$votu))
  base <- coalesce(na_if(votus$especie, ""), na_if(votus$nr_nombre, ""), votus$votu)
  if_else(votus$clase_final %in% c("pariente_debil", "divergente"), paste0("afín a ", base), base)
}

GRUPOS_VIRUS <- c(artropodo = "Virus de artrópodos", fago = "Fagos (microbiota)")

# Fagos (virus de la microbiota bacteriana) frente a virus de artropodos.
grupo_virus <- function(votus) {
  if ("grupo_virus" %in% names(votus)) return(unname(GRUPOS_VIRUS[votus$grupo_virus]))
  nombres <- paste(votus$especie, votus$nr_nombre)
  if_else(grepl("phage|Leviviricetes|Caudoviricetes", nombres, ignore.case = TRUE),
          GRUPOS_VIRUS[["fago"]], GRUPOS_VIRUS[["artropodo"]])
}

# Presencia por especie y muestra: presente si alguna vOTU de esa especie
# cubre >= 70 % de su longitud. La profundidad se normaliza por millon de
# pares recortados (no por pares no-huesped: ver M36 y M74).
presencia_especies <- function(votus, cob, lecturas, umbral = 70) {
  v <- votus |>
    filter(es_virus == "si") |>
    mutate(especie_fig = nombre_virus(pick(everything())),
           grupo = grupo_virus(pick(everything())))
  cob |>
    inner_join(select(v, votu, especie_fig, grupo, familia), by = "votu") |>
    filter(pct_cubierto >= umbral) |>
    left_join(select(lecturas, muestra, trimmed_pairs), by = "muestra") |>
    mutate(prof_norm = profundidad / (trimmed_pairs / 1e6)) |>
    group_by(especie_fig, grupo, muestra) |>
    summarise(familia = first(na.omit(na_if(familia, "")), default = NA_character_),
              prof_norm = max(prof_norm), n_votus = n(), .groups = "drop")
}

# Tasa de infeccion individual a partir de pools del mismo tamano.
# MLE: p = 1 - (1 - x/m)^(1/k); IC transformando el Clopper-Pearson de la
# prevalencia de pools. MIR = x / (m * k). Ambas por 1,000 mosquitos.
tasa_infeccion <- function(x, m, k) {
  ci <- binom.test(x, m)$conf.int
  f <- function(p) 1000 * (1 - (1 - p)^(1 / k))
  tibble(pools_pos = x, pools = m, mir = 1000 * x / (m * k),
         mle = f(x / m), mle_inf = f(ci[1]), mle_sup = f(ci[2]))
}
