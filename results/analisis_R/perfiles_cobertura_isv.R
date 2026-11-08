# =============================================================================
# perfiles_cobertura_isv.R
#
# Figura de tres paneles a partir de isv_perfiles.sh:
#   (a) virus que se replica  -> cobertura pareja de extremo a extremo
#   (b) virus segmentado      -> los tres segmentos, tambien parejos
#   (c) repetido del hospedero-> un pico estrecho sobre fondo cero
#
# El panel (c) es la prueba visual de que HTV no esta en las muestras.
# Etiquetas en ingles, convencion "X-like" para virus emparentados.
#
# Datos esperados (congelados desde la farm):
#   datos/perfiles_w25.unmapped.tsv
#   datos/tamanos_libreria.unmapped.tsv   (opcional, para normalizar)
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ----------------------------------------------------------------- parametros
RUTA      <- "datos/perfiles_w25.unmapped.tsv"
SALIDA    <- "figuras/F2_perfiles_cobertura"
M_ANPHE   <- "M1_S1"     # positiva de anfevirus, 284x
M_PCLV    <- "M56_S7"    # la unica con PCLV
M_HTV     <- "M1_S1"     # la de mayor senal aparente de HTV

PAL <- c(anphevirus = "#028090",
         pclv       = "#C97B16",
         artefacto  = "#B04A3C")

tema_perfil <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    axis.title         = element_text(size = 9),
    plot.title         = element_text(size = 10, face = "bold"),
    plot.subtitle      = element_text(size = 8.5, colour = "grey35"),
    strip.text         = element_text(size = 8.5, face = "bold", hjust = 0),
    legend.position    = "none"
  )

# ------------------------------------------------------------------- lectura
perf <- readr::read_tsv(RUTA, show_col_types = FALSE)

# Escala log: 0 no es representable, se lleva a la mitad del minimo visible
piso <- 0.3
perf <- perf |> mutate(prof = pmax(prof_media, piso))

eje_y <- scale_y_log10(
  breaks = c(0.3, 1, 10, 100, 1000, 10000, 1e5),
  labels = c("0", "1", "10", "100", "1k", "10k", "100k")
)

# ------------------------------------------------- (a) anfevirus, referencia propia
a <- perf |>
  filter(muestra == M_ANPHE, ref == "AeAV_MERIDA__own1") |>
  ggplot(aes(pos_ini, prof)) +
  geom_area(fill = PAL[["anphevirus"]], alpha = 0.22) +
  geom_line(colour = PAL[["anphevirus"]], linewidth = 0.4) +
  eje_y +
  scale_x_continuous(labels = function(x) paste0(x / 1000, " kb")) +
  labs(
    title    = "a   Replicating virus — Merida anphevirus",
    subtitle = paste0(M_ANPHE, " · own 12,945 bp reference (CheckV complete) · even coverage end to end"),
    x = "genome position", y = "mean depth (25 bp windows)"
  ) +
  tema_perfil

# ------------------------------------------------------- (b) PCLV, tres segmentos
etiq_pclv <- c(
  PCLV__NC_038262 = "L — RdRp (6,806 bp)",
  PCLV__NC_038261 = "M — glycoprotein (3,935 bp)",
  PCLV__NC_038263 = "S — nucleocapsid (1,358 bp)"
)

b <- perf |>
  filter(muestra == M_PCLV, virus == "PCLV") |>
  mutate(seg = factor(etiq_pclv[ref], levels = unname(etiq_pclv))) |>
  ggplot(aes(pos_ini, prof)) +
  geom_area(fill = PAL[["pclv"]], alpha = 0.22) +
  geom_line(colour = PAL[["pclv"]], linewidth = 0.4) +
  facet_wrap(~ seg, nrow = 1, scales = "free_x") +
  eje_y +
  scale_x_continuous(labels = function(x) paste0(x / 1000, " kb")) +
  labs(
    title    = "b   Segmented virus — Phasi Charoen-like virus",
    subtitle = paste0(M_PCLV, " · the only positive sample of 22 · all three segments covered"),
    x = "segment position", y = "mean depth"
  ) +
  tema_perfil

# --------------------------------------------- (c) HTV: el artefacto del hospedero
c_df <- perf |> filter(muestra == M_HTV, ref == "HTV__PZ595460")

c <- ggplot(c_df, aes(pos_ini, prof)) +
  geom_area(fill = PAL[["artefacto"]], alpha = 0.22) +
  geom_line(colour = PAL[["artefacto"]], linewidth = 0.4) +
  annotate("segment", x = 1618, xend = 1672, y = max(c_df$prof) * 1.9,
           yend = max(c_df$prof) * 1.9, colour = "grey25", linewidth = 0.4) +
  annotate("text", x = 1645, y = max(c_df$prof) * 3.2, size = 2.9, colour = "grey25",
           label = "~55 bp · 86% identity to a\ntandem repeat on AaegL5 chr 1") +
  eje_y +
  scale_x_continuous(labels = function(x) paste0(x, " bp")) +
  labs(
    title    = "c   Host repeat, not a virus — HTV capsid reference",
    subtitle = paste0(M_HTV, " · up to 80,769 reads, all in one window · RdRp segment at zero in all 22 samples"),
    x = "reference position", y = "mean depth"
  ) +
  tema_perfil

# ------------------------------------------------------------------- montaje
fig <- a / b / c + plot_layout(heights = c(1, 1, 1))

dir.create(dirname(SALIDA), showWarnings = FALSE, recursive = TRUE)
ggsave(paste0(SALIDA, ".pdf"), fig, width = 7.5, height = 8.2, device = cairo_pdf)
ggsave(paste0(SALIDA, ".png"), fig, width = 7.5, height = 8.2, dpi = 300, bg = "white")

message("escrito: ", SALIDA, ".pdf / .png")

# =============================================================================
# Suplementaria: el anfevirus en las 22 muestras, para ver el gradiente de carga
# =============================================================================
sup <- perf |>
  filter(ref == "AeAV_MERIDA__own1") |>
  mutate(muestra = factor(muestra, levels = perf |>
           filter(ref == "AeAV_MERIDA__own1") |>
           group_by(muestra) |> summarise(m = mean(prof_media)) |>
           arrange(desc(m)) |> pull(muestra))) |>
  ggplot(aes(pos_ini, prof)) +
  geom_area(fill = PAL[["anphevirus"]], alpha = 0.25) +
  facet_wrap(~ muestra, ncol = 4) +
  eje_y +
  scale_x_continuous(labels = function(x) paste0(x / 1000, "kb")) +
  labs(
    title    = "Merida anphevirus coverage across all 22 samples",
    subtitle = "ordered by mean depth · uniform profiles are replicating infections, flat zero is true absence",
    x = "genome position", y = "mean depth (25 bp windows)"
  ) +
  tema_perfil + theme(strip.text = element_text(size = 7.5))

ggsave("figuras/S2_anfevirus_22_muestras.png", sup,
       width = 9, height = 10, dpi = 300, bg = "white")
message("escrito: figuras/S2_anfevirus_22_muestras.png")
