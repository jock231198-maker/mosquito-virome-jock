#!/usr/bin/env Rscript
# Prepara las capas base del mapa (se corre UNA vez; el resultado va al repo).
# Fuente: INEGI, Marco Geoestadístico 2023 (municipios), en GeoJSON WGS84 de
#   https://github.com/PhantomInsights/mexico-geojson  (carpeta 2023/states/)
# Uso: Rscript prep_geo.R <Yucatán.json> <Campeche.json> <Quintana Roo.json>

suppressPackageStartupMessages(library(sf))
sf_use_s2(FALSE)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
dir.create("geo", showWarnings = FALSE)

edos <- lapply(args, st_read, quiet = TRUE)
yuc  <- edos[[1]]

# 1) Municipios de la zona de colecta (recorte holgado alrededor de Mérida)
caja <- st_as_sfc(st_bbox(c(xmin = -90.05, xmax = -89.45,
                            ymin = 20.85,  ymax = 21.32), crs = 4326))
muni <- yuc[st_intersects(yuc, caja, sparse = FALSE)[, 1], c("CVEGEO", "NOMGEO")]
muni <- st_make_valid(muni)
st_write(muni, "geo/municipios_zona.gpkg", delete_dsn = TRUE, quiet = TRUE)

# 2) Contorno de los tres estados de la península (para el recuadro)
pen <- do.call(rbind, lapply(edos, function(e) {
  e <- st_make_valid(e)
  st_sf(estado = e$NOM_ENT[1], geometry = st_union(e))
}))
pen <- st_simplify(pen, dTolerance = 0.01, preserveTopology = TRUE)
st_write(pen, "geo/peninsula_estados.gpkg", delete_dsn = TRUE, quiet = TRUE)

cat("Listo:", nrow(muni), "municipios en geo/municipios_zona.gpkg;",
    nrow(pen), "estados en geo/peninsula_estados.gpkg\n")
