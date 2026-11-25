#!/usr/bin/env Rscript
# =============================================================================
# Mapa de sitios de colecta de las muestras del estudio (MERIDA):
#   Ae. aegypti (RNA total + ARN pequeño, y extras solo ARN pequeño)
#   Cx. quinquefasciatus (solo ARN pequeño)
#
# Uso:
#   Rscript mapa_colecta.R Colecta_Aedes_MiguelGarcia.xlsx [prefijo_salida]
#
# Entradas (misma carpeta que este script):
#   localidades.tsv            coordenadas por localidad (columna `verificado`)
#   geo/municipios_zona.gpkg   municipios INEGI 2023 alrededor de Mérida
#   geo/peninsula_estados.gpkg contorno de Yucatán, Campeche y Quintana Roo
#   (las dos capas se generan con prep_geo.R)
#
# Salidas:
#   <prefijo>.pdf / <prefijo>.png   figura (A: mapa, B: pools por sitio)
#   <prefijo>_muestras_estudio.tsv  muestra → localidad, coordenadas y tipo de
#                                   dato (para cruzar después con el viroma)
#
# Paquetes: sf, ggplot2, dplyr, tidyr, readxl, readr, stringr, ggrepel,
#           ggspatial, patchwork. Opcional: rosm (teselas de fondo).
# Instalar una vez (en R):
#   install.packages(c("sf", "ggplot2", "dplyr", "tidyr", "readxl", "readr",
#                      "stringr", "ggrepel", "ggspatial", "patchwork", "ragg"))
#   (en conda: r-ragg)
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(ggplot2); library(dplyr); library(tidyr)
  library(readxl); library(readr); library(stringr)
  library(ggrepel); library(ggspatial); library(patchwork)
})
sf_use_s2(FALSE)

# ---- Parámetros -------------------------------------------------------------
# Muestras del estudio (pools secuenciados). Edita aquí si cambia el conjunto.
AAE_TOTAL <- c("M1", "M2", "M7", "M9", "M11", "M20", "M21", "M22", "M35", "M36",
               "M48", "M55", "M56", "M57", "M58", "M59", "M60", "M62", "M71",
               "M72", "M74", "M79")                 # RNA total + ARN pequeño
AAE_SRNA  <- c("M41", "M45", "M61", "M78")           # solo ARN pequeño (extras)
CQU_SRNA  <- c("M16", "M32.1", "M32.2", "M43", "M65", "M69", "M77")  # solo ARN pequeño
# Pools con algunos individuos de otro sitio (comentarios de la tabla); se
# asignan al sitio mayoritario, que es el de su propia fila.
POOLS_MIXTOS <- c("M16", "M43", "M62", "M65", "M74")
# Fondo de teselas (CARTO light vía OSM). Requiere internet y el paquete rosm.
USAR_TESELAS <- FALSE
CRS_UTM <- 32616   # WGS 84 / UTM zona 16N

# Colores por grupo (azul = Ae. aegypti en dos tonos; naranja = Cx. quinquefasciatus)
GRUPOS  <- c("aae_total", "aae_srna", "cqu_srna")
COL_GRP <- c(aae_total = "#2a78d6", aae_srna = "#86b6ef", cqu_srna = "#eb6834")
LAB_GRP <- list(aae_total = quote(italic("Ae. aegypti") ~ "(RNA-seq + small RNA)"),
                aae_srna  = quote(italic("Ae. aegypti") ~ "(small RNA only)"),
                cqu_srna  = quote(italic("Cx. quinquefasciatus") ~ "(small RNA only)"))
TXT1 <- "#1f1f1d"; TXT2 <- "#5f5e5a"; GRIS <- "#b9b8b2"

# ---- Argumentos --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Uso: Rscript mapa_colecta.R <tabla.xlsx> [prefijo]")
xlsx    <- args[1]
prefijo <- if (length(args) >= 2) args[2] else "fig_mapa_estudio"

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

# ---- 3. Muestras del estudio --------------------------------------------------
estudio <- bind_rows(
  tibble(id = AAE_TOTAL, especie = "Aedes aegypti",          grupo = "aae_total"),
  tibble(id = AAE_SRNA,  especie = "Aedes aegypti",          grupo = "aae_srna"),
  tibble(id = CQU_SRNA,  especie = "Culex quinquefasciatus", grupo = "cqu_srna")
) %>%
  mutate(id_origen = sub("\\.\\d+$", "", id))   # M32.1 / M32.2 salen de M32

# Se casa por ID de origen Y especie (la tabla repite M9 y M10 con otra especie)
mues <- estudio %>%
  left_join(select(col, id, especie, localidad, lat, lon, verificado, fecha, metodo),
            by = c("id_origen" = "id", "especie"))
if (any(is.na(mues$localidad)) || anyDuplicated(mues$id))
  stop("No se pudo asignar localidad única a: ",
       paste(mues$id[is.na(mues$localidad) | duplicated(mues$id)], collapse = ", "))
mues <- mutate(mues, pool_mixto = id %in% POOLS_MIXTOS)

write_tsv(select(mues, id, especie, grupo, localidad, lat, lon, verificado,
                 fecha, metodo, pool_mixto),
          paste0(prefijo, "_muestras_estudio.tsv"), na = "")

por_sitio <- mues %>%
  group_by(localidad, lat, lon, verificado) %>%
  summarise(n = n(), .groups = "drop")
por_grupo <- mues %>% count(localidad, grupo) %>%
  complete(localidad, grupo = GRUPOS, fill = list(n = 0))
orden_loc <- por_sitio %>% arrange(n, localidad) %>% pull(localidad)

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
  scale_size_area(max_size = 9, name = "Pools",
                  breaks = c(1, 4, 8), limits = c(0, max(xy$n)),
                  guide = guide_legend(override.aes = list(fill = TXT1,
                                                           colour = "white"))) +
  geom_text_repel(data = xy, aes(X, Y, label = localidad), size = 3.1,
                  colour = TXT1, seed = 7, box.padding = 0.8, force = 6, max.iter = 20000,
                  point.size = 9 * sqrt(xy$n / max(xy$n)), point.padding = 0.4, min.segment.length = 0.2,
                  segment.colour = TXT2, segment.size = 0.3) +
  annotation_scale(location = "bl", width_hint = 0.25, style = "ticks",
                   line_col = TXT2, text_col = TXT2, text_cex = 0.7) +
  annotation_north_arrow(location = "tr", height = unit(0.9, "cm"),
                         width = unit(0.7, "cm"),
                         style = north_arrow_orienteering(text_size = 7,
                                   fill = c("white", TXT2), line_col = TXT2)) +
  scale_x_continuous(breaks = seq(-89.90, -89.60, by = 0.05)) +
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

# ---- 5. Panel B: pools por sitio y tipo de dato --------------------------------
tot <- count(mues, localidad)

pB <- ggplot(por_grupo,
             aes(n, factor(localidad, levels = orden_loc),
                 fill = factor(grupo, levels = rev(GRUPOS)))) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.4) +
  geom_text(data = tot, inherit.aes = FALSE,
            aes(n, factor(localidad, levels = orden_loc), label = n),
            hjust = -0.4, size = 2.9, colour = TXT2) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = COL_GRP, breaks = GRUPOS, labels = LAB_GRP,
                    name = NULL) +
  scale_x_continuous(breaks = seq(0, 20, 2),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Pools sequenced", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "#e6e5df", linewidth = 0.25),
        axis.text = element_text(colour = TXT2),
        axis.text.y = element_text(colour = TXT1, size = 8.5),
        axis.title.x = element_text(colour = TXT2, size = 8),
        legend.position = "bottom", legend.direction = "vertical",
        legend.text = element_text(colour = TXT2, size = 8),
        legend.key.size = unit(0.35, "cm"))

# ---- 6. Componer y guardar ---------------------------------------------------
fechas <- range(mues$fecha, na.rm = TRUE)
pie <- sprintf(paste0("Pools sequenced in this study (n = %d), collected %s to %s. ",
                      "Point size: pools per site. Shaded: Mérida municipality (INEGI 2023).\n",
                      "Pools %s include a few individuals from another site and are ",
                      "assigned to their majority site."),
               nrow(mues), format(fechas[1], "%d %b %Y"), format(fechas[2], "%d %b %Y"),
               paste(POOLS_MIXTOS, collapse = ", "))
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
       device = "pdf")
# PNG con ragg si está instalado: mismas fuentes y cursivas en Mac y Linux/WSL
ggsave(paste0(prefijo, ".png"), fig, width = 25, height = 12.5, units = "cm",
       dpi = 300, bg = "white",
       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png")

# ---- 7. Avisos ----------------------------------------------------------------
message(sprintf("Filas de campo usadas: %d (se omitieron %d filas sin localidad: pools/submuestras)",
                nrow(col), n_sin_lugar))
print(count(mues, localidad, grupo) %>% pivot_wider(names_from = grupo,
      values_from = n, values_fill = 0))
if (hay_sin_verificar)
  message("OJO - coordenadas sin verificar: ",
          paste(pts$localidad[pts$verificado != "si"], collapse = ", "),
          "  → edita lat/lon y pon verificado=si en localidades.tsv")
message("Escrito: ", prefijo, ".pdf/.png y ", prefijo, "_muestras_estudio.tsv")
