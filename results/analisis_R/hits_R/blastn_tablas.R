#!/usr/bin/env Rscript
# blastn_tablas.R — resultados de blastn de los contigs contra la base viral, en TSV,
# y heatmap de los N mejores hits.
#
# Uso (desde la carpeta donde estan R/ y datos_hits_vN/):
#   Rscript blastn_tablas.R                       # ultima datos_hits_vN, top 20
#   Rscript blastn_tablas.R datos_hits_v2 30      # otra version / otro top
#
# Sale en tablas/blastn/:
#   blastn_todos_los_hits.tsv     todos los HSP, con encabezado, muestra, familia, grupo
#   blastn_mejor_hit_contig.tsv   un hit por contig (el de mayor bitscore)
#   blastn_especies.tsv           resumen por especie (virus de ARN, >= 30 % cobertura)
#   blastn_topN.tsv               las N primeras de la anterior (las del heatmap)
#   blastn_topN_muestra.tsv       matriz especie x muestra (n de contigs)
# y figuras/11_blastn_topN_heatmap.{pdf,png}

args <- commandArgs(trailingOnly = TRUE)
vers <- sort(Sys.glob("datos_hits_v*"))
DATOS <- if (length(args) >= 1) args[1] else tail(vers[order(as.integer(sub(".*_v", "", vers)))], 1)
N_TOP <- if (length(args) >= 2) as.integer(args[2]) else 20L

source("R/blastn.R")
OUT <- file.path("tablas", "blastn"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create("figuras", showWarnings = FALSE)

b      <- leer_blastn(DATOS)
d_hits <- leer_hits(DATOS)
todos  <- tabla_blastn_todos(b, d_hits)
mejor  <- tabla_blastn_mejor(todos)
esp    <- tabla_blastn_especies(mejor)
top    <- slice_head(esp, n = N_TOP)

stopifnot(n_distinct(todos$contig) == nrow(mejor))
message(sprintf("%s: %s HSP, %d contigs con hit blastn, %d especies de virus de ARN nombrables",
                DATOS, format(nrow(todos), big.mark = ","), nrow(mejor), nrow(esp)))
print(count(mejor, group, sort = TRUE))

escribir <- function(x, f) write_tsv(x, file.path(OUT, f), na = "")
escribir(select(todos, -en_genomad), "blastn_todos_los_hits.tsv")
escribir(mejor,  "blastn_mejor_hit_contig.tsv")
escribir(esp,    "blastn_especies.tsv")
escribir(top,    sprintf("blastn_top%d.tsv", N_TOP))

matriz <- mejor |>
  filter(group == "candidato", qcovs >= QCOV_NOMBRE, species %in% top$species) |>
  count(species, muestra) |>
  mutate(muestra = factor(muestra, orden_muestras(mejor$muestra))) |> arrange(muestra) |>
  pivot_wider(names_from = muestra, values_from = n, values_fill = 0) |>
  right_join(select(top, virus, species, family, status), by = "species") |>
  relocate(virus, species, family, status) |>
  arrange(match(species, top$species))
escribir(matriz, sprintf("blastn_top%d_muestra.tsv", N_TOP))

p <- figura_blastn_top(mejor, esp, N_TOP)
alto <- 2 + 0.25 * nrow(top)
base <- sprintf("figuras/11_blastn_top%d_heatmap", N_TOP)
ggsave(paste0(base, ".pdf"), p, width = 10, height = alto)
ggsave(paste0(base, ".png"), p, width = 10, height = alto, dpi = 300, device = ragg::agg_png)

message("Listo: ", OUT, "/*.tsv y ", base, ".{pdf,png}")
