# R/hits.R — tablas y figura de los hits virales por contig (blastn + DIAMOND
# sobre los 6,004 contigs de la union, contra viral_custom).
#
# Fuente unica: datos_hits_vN/hits_por_contig.tsv. Todo lo demas (por virus,
# muestra x virus, grupos) se calcula aqui, para que no haya dos versiones de
# la misma cuenta.
#
# Las etiquetas que salen en tablas y figuras van en ingles y con la convencion
# "X-like" para parientes.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(stringr); library(ggplot2)
})

TOTAL_CONTIGS <- 6004   # contigs de union_all/all_final.fasta
N_MUESTRAS    <- 22
TOTAL_GENOMAD <- 1144   # contigs llamados virales por geNomad (all_final_virus_summary.tsv)

# ---- Curacion manual -------------------------------------------------------
# Evidencia de que es material del hospedero (cribado contra AaegL5 o rRNA ya
# demostrado). Estos no cuentan como virus aunque geNomad llame algun contig.
DERIVADO_HUESPED <- c(
  "Panayarvi Aedes partiti-like virus 2" = "EVE in AaegL5 (genome screen)",
  "Tick-borne encephalitis virus"        = "EVE in AaegL5 (genome screen)",
  "Hanko iflavirus 1"                    = "Host-derived: genome screen / rRNA"
)
# Sospechosos: pasan a "Unconfirmed" aunque geNomad llame algun contig (p. ej.
# Hattula, 2 de 130). La nota se muestra en la tabla.
NOTAS <- c(
  "Iflavirus sp."              = "Ubiquitous, no geNomad support: rRNA pattern (cf. Hanko)",
  "Hattula partiti-like virus" = "Ubiquitous, short nt hit: possible EVE",
  "Oropouche virus"            = "Short nt hit; protein closest to Culex mosquito virus 1",
  "Kamiti River virus"         = "Short nt hit; possible flavivirus EVE",
  "Hameenlinna phasivirus"     = "No protein hit in any database"
)

GRUPOS <- c(
  candidato              = "RNA virus families (candidates)",
  ADN_u_otro             = "DNA virus families",
  realm_desconocido      = "DNA virus families (no realm assigned)",
  fago                   = "Bacteriophages",
  retroelemento          = "Retroviruses / retroelements",
  sin_familia            = "Viral hit without family",
  debil                  = "Weak (blastn < 30% cover, no protein)",
  sin_hit                = "No viral hit"
)

# ---- Lectura ---------------------------------------------------------------
leer_hits <- function(dir) {
  f <- file.path(dir, "hits_por_contig.tsv")
  stopifnot(file.exists(f))
  d <- read_tsv(f, col_types = cols(.default = "c"), na = c("", "NA"))
  esperadas <- c("muestra","ensamblador","contig","longitud","categoria","nombre",
                 "nivel","fuente","familia","realm","grupo","en_genomad",
                 "bn_pident","bn_qcov","dm_pident","dm_qcov")
  falta <- setdiff(esperadas, names(d))
  if (length(falta)) stop("Faltan columnas en hits_por_contig.tsv: ", paste(falta, collapse = ", "))
  d |>
    mutate(across(c(longitud, bn_pident, bn_qcov, dm_pident, dm_qcov), as.numeric),
           en_genomad = en_genomad == "si")
}

orden_muestras <- function(x) str_sort(unique(x), numeric = TRUE)

# Nombre para presentar: especie (>= 90 % nt y >= 30 % cobertura) tal cual;
# pariente -> "X-like". Nombres genericos ("sp.") no se tocan, y uno que ya
# lleva "-like" pasa a "X relative" para no dar "partiti-like virus-like".
etiqueta <- function(nombre, nivel) {
  case_when(
    nivel == "especie"                 ~ nombre,
    str_detect(nombre, "\\bsp\\.?$")   ~ nombre,
    str_detect(nombre, "-like")        ~ paste(nombre, "relative"),
    TRUE                               ~ paste0(nombre, "-like")
  )
}

# ---- Tablas ----------------------------------------------------------------
tabla_grupos <- function(d) {
  cuenta <- d |>
    mutate(g = if_else(nivel == "debil", "debil", grupo)) |>
    count(g, en_genomad) |>
    pivot_wider(names_from = en_genomad, values_from = n, values_fill = 0) |>
    rename(no = `FALSE`, si = `TRUE`) |>
    mutate(total = si + no)
  resto <- TOTAL_CONTIGS - sum(cuenta$total)
  sin_hit <- tibble(g = "sin_hit", no = NA_integer_, si = TOTAL_GENOMAD - sum(cuenta$si), total = resto)
  bind_rows(cuenta, sin_hit) |>
    mutate(g = factor(g, levels = names(GRUPOS))) |>
    arrange(g) |>
    transmute(Group = unname(GRUPOS[as.character(g)]),
              Contigs = total,
              `% of 6,004` = round(100 * total / TOTAL_CONTIGS, 1),
              `Called viral by geNomad` = si)
}

# Por virus: solo familias de virus de ARN, sin los debiles.
tabla_virus <- function(d) {
  d |>
    filter(grupo == "candidato", nivel != "debil") |>
    mutate(virus = etiqueta(nombre, nivel)) |>
    group_by(virus, nombre, familia) |>
    summarise(n_muestras  = n_distinct(muestra),
              n_contigs   = n(),
              n_genomad   = sum(en_genomad),
              n_ensamb    = n_distinct(ensamblador),
              max_long    = max(longitud, na.rm = TRUE),
              max_pid_nt  = suppressWarnings(max(bn_pident, na.rm = TRUE)),
              max_pid_aa  = suppressWarnings(max(dm_pident, na.rm = TRUE)),
              .groups = "drop") |>
    mutate(across(c(max_pid_nt, max_pid_aa), ~ if_else(is.finite(.x), .x, NA_real_)),
           estado = case_when(nombre %in% names(DERIVADO_HUESPED) ~ "Host-derived",
                              nombre %in% names(NOTAS)            ~ "Unconfirmed",
                              n_genomad > 0                       ~ "Confirmed",
                              TRUE                                ~ "Unconfirmed"),
           nota = coalesce(DERIVADO_HUESPED[nombre], NOTAS[nombre], ""),
           estado = factor(estado, c("Confirmed", "Unconfirmed", "Host-derived"))) |>
    arrange(estado, desc(n_muestras), desc(n_contigs))
}

presentar_virus <- function(tv) {
  tv |>
    transmute(Virus = virus, Family = familia, Status = as.character(estado),
              Samples = sprintf("%d/%d", n_muestras, N_MUESTRAS),
              Contigs = n_contigs, `geNomad-supported` = n_genomad,
              `Max length (bp)` = max_long,
              `Max nt id. (%)` = round(max_pid_nt, 1),
              `Max aa id. (%)` = round(max_pid_aa, 1),
              Note = nota)
}

# Muestra x virus (n de contigs), formato ancho para la hoja de calculo.
matriz_muestras <- function(d, tv) {
  d |>
    filter(grupo == "candidato", nivel != "debil") |>
    mutate(virus = etiqueta(nombre, nivel)) |>
    count(virus, muestra) |>
    mutate(muestra = factor(muestra, orden_muestras(d$muestra))) |>
    arrange(muestra) |>
    pivot_wider(names_from = muestra, values_from = n, values_fill = 0) |>
    right_join(select(tv, virus, Family = familia, Status = estado), by = "virus") |>
    relocate(Virus = virus, Family, Status) |>
    arrange(Status)
}

# ---- Figura ----------------------------------------------------------------
COLORES_ESTADO <- c(Confirmed = "#2A9D8F", Unconfirmed = "#E9C46A", `Host-derived` = "#9AA0A6")

figura_heatmap <- function(d, tv) {
  m <- d |>
    filter(grupo == "candidato", nivel != "debil") |>
    mutate(virus = etiqueta(nombre, nivel)) |>
    count(virus, muestra) |>
    inner_join(select(tv, virus, estado), by = "virus")
  ord_v <- tv |> arrange(desc(estado), n_muestras, n_contigs) |> pull(virus)
  todas <- orden_muestras(d$muestra)
  m <- m |> mutate(virus = factor(virus, ord_v), muestra = factor(muestra, todas))
  ggplot(m, aes(muestra, virus, fill = estado)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = n), size = 2.6, colour = "grey15") +
    scale_fill_manual(values = COLORES_ESTADO, name = NULL, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL,
         title = "RNA virus hits in assembled contigs",
         subtitle = "Contigs per pool (blastn/DIAMOND vs. custom viral DB). Confirmed: called viral by geNomad, no host evidence") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank(), legend.position = "top",
          plot.title.position = "plot")
}
