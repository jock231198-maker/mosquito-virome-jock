#!/usr/bin/env Rscript
# =============================================================================
# Mapa de sitios de colecta (MERIDA) + esfuerzo de colecta por especie y sexo
#
# Uso:
#   Rscript mapa_colecta.R Colecta_Aedes_MiguelGarcia.xlsx [salida_prefijo]
#
# Entradas (misma carpeta que este script):
#   localidades.tsv            coordenadas por localidad (columna `verificado`)
#   geo/municipios_zona.gpkg   municipios INEGI 2023 alrededor de Mérida
#   geo/peninsula_estados.gpkg contorno de Yucatán, Campeche y Quintana Roo
#   (las dos capas se generan con prep_geo.R)
#
# Salidas:
#   <prefijo>.pdf / <prefijo>.png   figura (A: mapa, B: individuos por sitio)
#   <prefijo>_colecta_limpia.tsv    tabla limpia muestra → localidad (para
#                                   cruzar después con el viroma)
#
# Paquetes: sf, ggplot2, dplyr, tidyr, readxl, readr, stringr, ggrepel,
#           ggspatial, patchwork. Opcional: rosm (teselas de fondo).
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(tidyr)
  library(readxl); library(readr); library(stringr)
  library(ggrepel); library(ggspatial); library(patchwork)
})
sf_use_s2(FALSE)

# ---- Parámetros -------------------------------------------------------------
# Especies a incluir (cambia aquí; p. ej. solo Ae. aegypti + Cx. quinquefasciatus)
ESPECIES <- c("Aedes aegypti", "Culex quinquefasciatus",
              "Aedes albopictus", "Aedes taeniorhynchus")
# Fondo de teselas (CARTO light vía OSM). Requiere internet y el paquete rosm.
USAR_TESELAS <- FALSE
CRS_UTM <- 32616   # WGS 84 / UTM zona 16N

# Paleta categórica en orden fijo (validada: CVD y visión normal PASS)
COL_ESP <- c("Aedes aegypti"          = "#2a78d6",
             "Culex quinquefasciatus" = "#eb6834",
             "Aedes albopictus"       = "#1baf7a",
             "Aedes taeniorhynchus"   = "#eda100")
TXT1 <- "#1f1f1d"; TXT2 <- "#5f5e5a"; GRIS <- "#b9b8b2"

# ---- Argumentos --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Uso: Rscript mapa_colecta.R <tabla.xlsx> [prefijo]")
xlsx    <- args[1]
prefijo <- if (length(args) >= 2) args[2] else "fig_mapa_colecta"

aqui <- tryCatch(dirname(normalizePath(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)))),
          error = function(e) ".")
ruta <- function(...) file.path(aqui, ...)

# ---- 1. Leer y limpiar la tabla de colecta ----------------------------------
raw <- read_excel(xlsx, col_types = "text")

fix_fecha <- function(x) {
  out <- rep(as.Date(NA), length(x))
  num <- suppressWarnings(as.numeric(x))
  es_serial <- !is.na(num)
  out[es_serial] <- as.Date(num[es_serial], origin = "1899-12-30")
  # texto dd/mm/aaaa (incluye el error "04/09/20224" → 2024)
  txt <- !es_serial & !is.na(x)
  t <- str_match(x[txt], "^(\\d{1,2})/(\\d{1,2})/(\\d{4,5})$")
  anio <- ifelse(nchar(t[, 4]) == 5,            # "20224": dígito de más
                 paste0(substr(t[, 4], 1, 2), substr(t[, 4], 4, 5)), t[, 4])
  out[txt] <- as.Date(sprintf("%s-%s-%s", anio, t[, 3], t[, 2]))
  out
}

col <- tibble(
  id        = raw$ID,
  fecha     = fix_fecha(raw$`FECHA COLECTA`),
  lugar     = str_squish(raw$`LUGAR DE COLECTA`),
  metodo    = str_squish(raw$`METODO DE COLECTA`),
  n         = as.numeric(raw$`No. INDIVIDUOS`),
  sexo      = str_to_lower(str_squish(raw$SEXO)),
  especie   = str_squish(raw$ESPECIE),
  rnaseq    = suppressWarnings(as.numeric(raw$`Transcriptoma LANCIS/INMEGEN`)) %in% 1,
  comentario = raw$Comentario
)

# Filas sin localidad = pools de secuenciación (PoolMxx) y submuestras
# (M32.1/M32.2): sus individuos ya están contados en la fila de campo de origen.
n_sin_lugar <- sum(is.na(col$lugar))
col <- filter(col, !is.na(lugar))

# IDs repetidos en la tabla original (M9, M10): se desambiguan con la localidad
col <- col %>%
  group_by(id) %>%
  mutate(id_unico = if (n() > 1) paste0(id, "_", str_sub(lugar, 1, 3)) else id) %>%
  ungroup()

# ---- 2. Coordenadas por localidad -------------------------------------------
loc <- read_tsv(ruta("localidades.tsv"), show_col_types = FALSE)
falta <- setdiff(unique(col$lugar), loc$lugar_tabla)
if (length(falta)) stop("Localidades sin coordenadas en localidades.tsv: ",
                        paste(falta, collapse = ", "))
col <- left_join(col, loc, by = c("lugar" = "lugar_tabla"))

write_tsv(select(col, id, id_unico, fecha, localidad, lat, lon, metodo,
                 especie, sexo, n, rnaseq, comentario),
          paste0(prefijo, "_colecta_limpia.tsv"), na = "")

# ---- 3. Resumen para la figura -----------------------------------------------
esp <- filter(col, especie %in% ESPECIES)

por_sitio <- esp %>%
  group_by(localidad, lat, lon, verificado) %>%
  summarise(n = sum(n, na.rm = TRUE), .groups = "drop")

por_esp_sexo <- esp %>%
  mutate(sexo = recode(sexo, hembra = "Females", macho = "Males")) %>%
  group_by(localidad, especie, sexo) %>%
  summarise(n = sum(n, na.rm = TRUE), .groups = "drop") %>%
  complete(localidad, especie, sexo, fill = list(n = 0))

orden_loc <- por_sitio %>% arrange(n) %>% pull(localidad)
esp_presentes <- ESPECIES[ESPECIES %in% esp$especie]
lab_esp <- setNames(
  lapply(esp_presentes, function(s) {
    p <- strsplit(s, " ")[[1]]
    abbr <- c(Aedes = "Ae.", Culex = "Cx.")[p[1]]
    bquote(italic(.(paste(abbr, p[2]))))
  }), esp_presentes)

# ---- 4. Panel A: mapa --------------------------------------------------------
muni <- st_read(ruta("geo/municipios_zona.gpkg"), quiet = TRUE) %>% st_transform(CRS_UTM)
pen  <- st_read(ruta("geo/peninsula_estados.gpkg"), quiet = TRUE)
pts  <- st_as_sf(por_sitio, coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
        st_transform(CRS_UTM)

bb <- st_bbox(st_buffer(st_as_sfc(st_bbox(pts)), 4500))
xy <- cbind(st_drop_geometry(pts), st_coordinates(pts))

hay_sin_verificar <- any(pts$verificado != "si")

pA <- ggplot()
if (USAR_TESELAS && requireNamespace("rosm", quietly = TRUE)) {
  pA <- pA + annotation_map_tile(type = "cartolight", zoomin = 0, progress = "none")
}
pA <- pA +
  geom_sf(data = muni, aes(fill = NOMGEO == "Mérida"),
          colour = GRIS, linewidth = 0.3, alpha = if (USAR_TESELAS) 0.35 else 1) +
  scale_fill_manual(values = c(`TRUE` = "#ecebe5", `FALSE` = "#fafaf8"), guide = "none") +
  # verificados: relleno oscuro; sin verificar: hueco con borde oscuro
  geom_point(data = filter(xy, verificado == "si"), aes(X, Y, size = n),
             shape = 21, fill = TXT1, colour = "white", stroke = 0.8) +
  geom_point(data = filter(xy, verificado != "si"), aes(X, Y, size = n),
             shape = 21, fill = "white", colour = TXT1, stroke = 1) +
  scale_size_area(max_size = 9, name = "Individuals",
                  breaks = c(25, 100, 250), limits = c(0, max(xy$n)),
                  guide = guide_legend(override.aes = list(fill = TXT1,
                                                           colour = "white"))) +
  geom_text_repel(data = xy, aes(X, Y, label = localidad), size = 3.1,
                  colour = TXT1, seed = 7, box.padding = 0.8, force = 6, max.iter = 20000,
                  point.size = sqrt(xy$n) * 0.9, point.padding = 0.3, min.segment.length = 0.2,
                  segment.colour = TXT2, segment.size = 0.3) +
  annotation_scale(location = "bl", width_hint = 0.25, style = "ticks",
                   line_col = TXT2, text_col = TXT2, text_cex = 0.7) +
  annotation_north_arrow(location = "tr", height = unit(0.9, "cm"),
                         width = unit(0.7, "cm"),
                         style = north_arrow_orienteering(text_size = 7,
                                   fill = c("white", TXT2), line_col = TXT2)) +
  coord_sf(crs = CRS_UTM, datum = 4326,
           xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")],
           expand = FALSE) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_line(colour = "#e6e5df", linewidth = 0.25),
        axis.text = element_text(colour = TXT2, size = 7),
        axis.title = element_blank(),
        panel.border = element_rect(colour = GRIS, fill = NA, linewidth = 0.4),
        legend.position = "bottom",
        legend.title = element_text(colour = TXT2, size = 8),
        legend.text = element_text(colour = TXT2, size = 7))

# Recuadro: península con el área del mapa
area <- st_as_sfc(bb) %>% st_transform(4326)
inset <- ggplot() +
  geom_sf(data = pen, fill = "#f2f1ec", colour = GRIS, linewidth = 0.25) +
  geom_sf(data = area, fill = NA, colour = TXT1, linewidth = 0.6) +
  coord_sf(xlim = c(-92.6, -86.6), ylim = c(17.7, 21.8), expand = FALSE) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = GRIS,
                                        linewidth = 0.4))

pA <- pA + inset_element(inset, left = 0.01, bottom = 0.52, right = 0.30,
                         top = 0.89, align_to = "panel",
                         ignore_tag = TRUE)

# ---- 5. Panel B: individuos por sitio, especie y sexo ------------------------
tot_sexo <- por_esp_sexo %>% group_by(localidad, sexo) %>%
  summarise(n = sum(n), .groups = "drop")

pB <- ggplot(por_esp_sexo,
             aes(n, factor(localidad, levels = orden_loc),
                 fill = factor(especie, levels = rev(esp_presentes)))) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.4) +
  geom_text(data = tot_sexo, inherit.aes = FALSE,
            aes(n, factor(localidad, levels = orden_loc), label = n),
            hjust = -0.25, size = 2.7, colour = TXT2) +
  facet_wrap(~sexo, nrow = 1) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = COL_ESP, breaks = esp_presentes,
                    labels = lab_esp, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Individuals collected", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "#e6e5df", linewidth = 0.25),
        axis.text = element_text(colour = TXT2),
        axis.text.y = element_text(colour = TXT1, size = 8.5),
        axis.title.x = element_text(colour = TXT2, size = 8),
        strip.text = element_text(colour = TXT1, face = "bold", hjust = 0),
        legend.position = "bottom",
        legend.text = element_text(colour = TXT2, size = 8),
        legend.key.size = unit(0.35, "cm")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# ---- 6. Componer y guardar ---------------------------------------------------
fechas <- range(col$fecha, na.rm = TRUE)
pie <- sprintf("Field collections %s – %s. Point size: individuals of the species shown. Shaded: Mérida municipality (INEGI 2023).",
               format(fechas[1], "%d %b %Y"), format(fechas[2], "%d %b %Y"))
if (hay_sin_verificar)
  pie <- paste(pie, "Hollow points: approximate coordinates (to be verified).")

tag <- function(p, t) wrap_elements(full = p) + labs(tag = t) +
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.tag.location = "plot")
fig <- (tag(pA, "A") | tag(pB, "B")) + plot_layout(widths = c(1.35, 1)) +
  plot_annotation(caption = pie,
                  theme = theme(plot.caption = element_text(colour = TXT2,
                                  size = 7.5, hjust = 0)))

ggsave(paste0(prefijo, ".pdf"), fig, width = 25, height = 12.5, units = "cm",
       device = if (capabilities("cairo")) cairo_pdf else "pdf")
ggsave(paste0(prefijo, ".png"), fig, width = 25, height = 12.5, units = "cm",
       dpi = 300, bg = "white")

# ---- 7. Avisos ----------------------------------------------------------------
message(sprintf("Filas de campo usadas: %d (se omitieron %d filas sin localidad: pools/submuestras)",
                nrow(col), n_sin_lugar))
message("Especies incluidas: ", paste(esp_presentes, collapse = ", "))
if (hay_sin_verificar)
  message("OJO - coordenadas sin verificar: ",
          paste(pts$localidad[pts$verificado != "si"], collapse = ", "),
          "  → edita lat/lon y pon verificado=si en localidades.tsv")
message("Escrito: ", prefijo, ".pdf/.png y ", prefijo, "_colecta_limpia.tsv")
