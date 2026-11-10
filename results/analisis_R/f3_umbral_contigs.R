# =============================================================================
# f3_umbral_contigs.R
#
# F3 — Cuánto de cada referencia reconstruyen los contigs, en cuatro niveles:
#        dos umbrales de longitud (>=500 / >=1000 bp)
#        x dos etapas (union dereplicada / catalogo de vOTUs)
#
# Cinco referencias no se mueven ni un decimal. Solo el segmento L de PCLV
# cambia, y cambia dos veces: por el umbral de entrada y por la clusterizacion.
# Esa forma —cinco cuartetos planos y uno que se abre— es el argumento.
#
# Datos (copiados desde $RESULTS_DIR/isv_sanity_2refAeAV/blast/):
#   datos/union_min500_refcov.tsv
#   datos/union_min1000_refcov.tsv
#   datos/votu_min500_refcov.tsv     (catalogo de 5 899)
#   datos/votu_748_refcov.tsv        (catalogo de 748)
# =============================================================================

library(ggplot2)
library(dplyr)
library(readr)

DIR_DAT <- "datos"
DIR_FIG <- "figuras"
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------- lectura
NIVELES <- tribble(
  ~archivo,                    ~nivel,                                ~umbral,
  "union_min1000_refcov.tsv",  "Dereplicated union",                  "≥ 1000 bp",
  "votu_748_refcov.tsv",       "vOTU catalogue",                      "≥ 1000 bp",
  "union_min500_refcov.tsv",   "Dereplicated union",                  "≥ 500 bp",
  "votu_min500_refcov.tsv",    "vOTU catalogue",                      "≥ 500 bp"
)

leer <- function(archivo, nivel, umbral) {
  ruta <- file.path(DIR_DAT, archivo)
  if (!file.exists(ruta)) stop("falta ", ruta)
  read_tsv(ruta, show_col_types = FALSE) |>
    transmute(ref, pct = fraccion * 100, nivel = nivel, umbral = umbral)
}

dat <- NIVELES |>
  purrr::pmap_dfr(function(archivo, nivel, umbral) leer(archivo, nivel, umbral))

# ------------------------------------------------------------- etiquetas
ETIQ <- c(
  "AeAV_MERIDA__own1" = "Merida anphevirus (own reference)",
  "CFAV__NC_001564"   = "CFAV",
  "PCLV__NC_038263"   = "PCLV  S — nucleocapsid",
  "PCLV__NC_038261"   = "PCLV  M — glycoprotein",
  "PCLV__NC_038262"   = "PCLV  L — RdRp"
)

dat <- dat |>
  filter(ref %in% names(ETIQ)) |>          # descarta AeAV RML-12, ya retirado del set
  mutate(
    referencia = factor(ETIQ[ref], levels = unname(ETIQ)),   # PCLV L arriba del todo
    serie = factor(paste0(nivel, ", ", umbral),
                   levels = c("Dereplicated union, ≥ 1000 bp",
                              "vOTU catalogue, ≥ 1000 bp",
                              "Dereplicated union, ≥ 500 bp",
                              "vOTU catalogue, ≥ 500 bp"))
  )

# Gris para el umbral viejo, teal para el nuevo; tono claro para el catálogo.
PAL <- c("Dereplicated union, ≥ 1000 bp" = "#6E7F85",
         "vOTU catalogue, ≥ 1000 bp"     = "#B4C0C4",
         "Dereplicated union, ≥ 500 bp"  = "#028090",
         "vOTU catalogue, ≥ 500 bp"      = "#79BEC6")

# ----------------------------------------------------------------- figura
f3 <- ggplot(dat, aes(x = pct, y = referencia, fill = serie)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(width = 0.78),
            hjust = -0.15, size = 2.4, colour = "grey25") +
  scale_fill_manual(values = PAL, name = NULL,
                    guide = guide_legend(reverse = TRUE, nrow = 2, byrow = TRUE)) +
  scale_x_continuous(limits = c(0, 112), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, " %"),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title    = "How much of each reference the contigs reconstruct",
    subtitle = "Two length thresholds × two stages. Only the L segment of PCLV moves — and it loses ground twice.",
    x = "fraction of the reference covered by contig HSPs", y = NULL,
    caption = paste0(
      "HSPs with pident ≥ 70 and length ≥ 100 bp, merged per reference.\n",
      "Scope: the five contigs that tile PCLV L all come from M56_S7 and from a single assembler. ",
      "A convincing test case, not a general effect across the 22 samples.")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
    axis.text.y        = element_text(face = "bold", hjust = 0),
    axis.title.x       = element_text(size = 9.5, colour = "grey30"),
    plot.title         = element_text(size = 12.5, face = "bold"),
    plot.subtitle      = element_text(size = 9.5, colour = "grey35"),
    plot.caption       = element_text(size = 8, colour = "grey50", hjust = 0),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    legend.position    = "top",
    legend.text        = element_text(size = 9),
    legend.key.size    = unit(11, "pt")
  )

ggsave(file.path(DIR_FIG, "F3_umbral_contigs.png"), f3,
       width = 8.6, height = 4.6, dpi = 300, bg = "white")
ggsave(file.path(DIR_FIG, "F3_umbral_contigs.pdf"), f3,
       width = 8.6, height = 4.6, device = cairo_pdf)

message("escrito: ", DIR_FIG, "/F3_umbral_contigs.png / .pdf")

# ------------------------------------------------- la tabla, para el texto
dat |>
  select(referencia, serie, pct) |>
  tidyr::pivot_wider(names_from = serie, values_from = pct) |>
  print(width = Inf)
