#!/usr/bin/env Rscript
# hits_tablas.R — tablas para presentar y figura de los hits virales por contig.
#
# Uso (desde results/analisis_R/):
#   Rscript hits_tablas.R                 # usa datos_hits_v1
#   Rscript hits_tablas.R datos_hits_v2
#
# Sale:
#   tablas/hits_virales.xlsx      todas las tablas, una por hoja (+ hoja README)
#   tablas/hits_virales.pptx      una diapositiva por tabla + la figura
#   tablas/T1_*.png, T2*_*.png    las tablas como imagen, para pegar donde sea
#   figuras/10_hits_heatmap.{pdf,png}

args  <- commandArgs(trailingOnly = TRUE)
DATOS <- if (length(args)) args[1] else "datos_hits_v1"

suppressPackageStartupMessages({ library(flextable); library(officer); library(writexl) })
source("R/hits.R")
dir.create("tablas", showWarnings = FALSE); dir.create("figuras", showWarnings = FALSE)

d  <- leer_hits(DATOS)
tg <- tabla_grupos(d)
tv <- tabla_virus(d)
pv <- presentar_virus(tv)
mm <- matriz_muestras(d, tv)

stopifnot(sum(tg$Contigs) == TOTAL_CONTIGS)
message(sprintf("%s: %d contigs con hit, %d virus de ARN (%d confirmados)",
                DATOS, nrow(d), nrow(tv), sum(tv$estado == "Confirmed")))

# ---- Hoja de calculo -------------------------------------------------------
readme <- tibble::tibble(
  Item = c("Source", "Input", "Status: Confirmed", "Status: Unconfirmed", "Status: Host-derived",
           "Naming", "Weak hits", "DNA families"),
  Definition = c(
    "blastn (megablast + dc-megablast) and DIAMOND blastx of the 6,004 dereplicated contigs vs. custom viral DB (RefSeq viral + GenBank arthropod-host viruses, 2026-09-25)",
    paste0(DATOS, "/hits_por_contig.tsv (see MANIFEST.txt)"),
    "At least one contig called viral by geNomad and no host/curation flag",
    "No contig called viral by geNomad, or flagged in curation (see Note)",
    "Genome screen against AaegL5 or rRNA evidence (manual curation in R/hits.R)",
    "Species name if blastn >= 90% nt id. over >= 30% of the contig; otherwise 'X-like' (closest relative, nt if cover >= 30%, else protein)",
    "blastn < 30% cover without protein hit: counted but not named",
    "Excluded as giant/endogenous DNA viruses: Nimaviridae, Mimiviridae, Phycodnaviridae, Poxviridae, polydnaviruses, baculoviruses, etc."))
contigs_cand <- d |> dplyr::filter(grupo == "candidato", nivel != "debil") |>
  dplyr::mutate(virus = etiqueta(nombre, nivel)) |>
  dplyr::relocate(virus, .after = contig)
write_xlsx(list(README = readme, groups = tg, viruses = pv, sample_x_virus = mm,
                rna_virus_contigs = contigs_cand, all_contigs = d),
           "tablas/hits_virales.xlsx")

# ---- Tablas para presentar (flextable) -------------------------------------
estilo <- function(ft, titulo) {
  ft |>
    set_caption(titulo) |>
    theme_booktabs() |>
    bold(part = "header") |>
    fontsize(size = 9, part = "all") |>
    font(fontname = "Helvetica", part = "all") |>
    align(j = which(sapply(ft$body$dataset, is.numeric)), align = "right", part = "all") |>
    autofit()
}
ft1 <- flextable(tg) |>
  colformat_num(j = "Called viral by geNomad", na_str = "-") |>
  colformat_int(j = "Contigs", big.mark = ",") |>
  estilo("Table 1. Viral similarity hits among the 6,004 assembled contigs") |>
  bold(i = 1, part = "body") |>
  add_footer_lines("Best hit per contig vs. custom viral DB. DNA-virus and phage hits in an RNA-seq library are host or bacterial genes with viral homologues.")

tabla2 <- function(est, titulo) {
  x <- pv |> dplyr::filter(Status == est) |> dplyr::select(-Status)
  if (!nrow(x)) return(NULL)
  if (all(x$Note == "")) x <- dplyr::select(x, -Note)
  flextable(x) |>
    colformat_num(na_str = "-") |>
    italic(j = "Virus", part = "body") |>
    estilo(titulo)
}
ft2a <- tabla2("Confirmed",    "Table 2a. RNA viruses supported by geNomad")
ft2b <- tabla2("Unconfirmed",  "Table 2b. RNA virus hits without geNomad support (to verify)")
ft2c <- tabla2("Host-derived", "Table 2c. Hits attributed to host sequence (EVE / rRNA)")
fts  <- Filter(Negate(is.null), list(T1_grupos = ft1, T2a_confirmados = ft2a,
                                     T2b_sin_confirmar = ft2b, T2c_huesped = ft2c))
for (n in names(fts)) save_as_image(fts[[n]], path = file.path("tablas", paste0(n, ".png")), res = 200)

# ---- Figura ----------------------------------------------------------------
p  <- figura_heatmap(d, tv)
alto <- 1.6 + 0.22 * nrow(tv)
ggsave("figuras/10_hits_heatmap.pdf", p, width = 9, height = alto)
ggsave("figuras/10_hits_heatmap.png", p, width = 9, height = alto, dpi = 300, device = ragg::agg_png)

# ---- PowerPoint ------------------------------------------------------------
# Tablas como tablas nativas (editables en PowerPoint), ajustadas al ancho;
# la figura como imagen escalada sin deformar.
ANCHO <- 9.2; ALTO <- 5.6
encajar <- function(png) {
  dim <- dim(png::readPNG(png))            # alto, ancho (px)
  r <- min(ANCHO / dim[2], ALTO / dim[1])
  external_img(png, width = dim[2] * r, height = dim[1] * r)
}
TITULO <- c(T1_grupos = "Viral hits among assembled contigs",
            T2a_confirmados = "RNA viruses confirmed by geNomad",
            T2b_sin_confirmar = "RNA virus hits to verify",
            T2c_huesped = "Hits attributed to host (EVE / rRNA)")
ppt <- read_pptx()
for (n in names(fts)) {
  ft <- fit_to_width(fts[[n]], max_width = ANCHO)
  ppt <- add_slide(ppt, layout = "Title Only") |>
    ph_with(TITULO[[n]], location = ph_location_type("title")) |>
    ph_with(ft, location = ph_location(left = 0.4, top = 1.6))
}
ppt <- add_slide(ppt, layout = "Title Only") |>
  ph_with("RNA virus hits per pool", location = ph_location_type("title")) |>
  ph_with(encajar("figuras/10_hits_heatmap.png"),
          location = ph_location(left = 0.4, top = 1.6), use_loc_size = FALSE)
print(ppt, target = "tablas/hits_virales.pptx")

message("Listo: tablas/ (xlsx, pptx, png) y figuras/10_hits_heatmap.{pdf,png}")
