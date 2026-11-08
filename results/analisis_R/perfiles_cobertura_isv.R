# =============================================================================
# perfiles_cobertura_isv.R
#
# Perfiles de cobertura del chequeo de sanidad de ISVs, UNA FIGURA POR PANEL.
#
#   F2a  virus que se replica      -> cobertura pareja de extremo a extremo
#   F2b  virus segmentado          -> los tres segmentos de PCLV
#   F2c  repetido del hospedero    -> un pico estrecho sobre fondo cero
#   S2   el anfevirus en las 22 muestras
#
# F2c es la prueba visual de que HTV no esta en las muestras.
# Etiquetas en ingles, convencion "X-like" para virus emparentados.
#
# Datos (congelados desde la farm):
#   datos/perfiles_w25.unmapped.tsv
# =============================================================================

library(ggplot2)
library(dplyr)
library(readr)

# ----------------------------------------------------------------- parametros
RUTA    <- "datos/perfiles_w25.unmapped.tsv"
DIR_FIG <- "figuras"
M_ANPHE <- "M1_S1"     # positiva de anfevirus, 284x de media
M_PCLV  <- "M56_S7"    # la unica con PCLV
M_HTV   <- "M1_S1"     # la de mayor senal aparente de HTV

PAL <- c(anphevirus = "#028090",
         pclv       = "#C97B16",
         artefacto  = "#B04A3C")

PISO <- 0.3   # la profundidad 0 no existe en escala log: se dibuja en el suelo

dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

tema_perfil <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    panel.spacing.x    = unit(16, "pt"),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9),
    plot.title         = element_text(size = 12, face = "bold"),
    plot.subtitle      = element_text(size = 9.5, colour = "grey35"),
    plot.caption       = element_text(size = 8, colour = "grey50", hjust = 0),
    plot.title.position = "plot",
    strip.text         = element_text(size = 10, face = "bold", hjust = 0),
    legend.position    = "none"
  )

escala_log <- function(mx) {
  br <- c(PISO, 1, 10, 100, 1000, 1e4, 1e5)
  lb <- c("0", "1", "10", "100", "1k", "10k", "100k")
  keep <- br <= mx * 3
  scale_y_log10(breaks = br[keep], labels = lb[keep],
                expand = expansion(mult = c(0, 0.08)))
}

kb <- function(x) ifelse(x == 0, "0", paste0(x / 1000, " kb"))

# ------------------------------------------------------------------- lectura
perf <- read_tsv(RUTA, show_col_types = FALSE) |>
  mutate(prof = pmax(prof_media, PISO))

NOTA_PISO <- "Log scale; windows with zero depth are drawn at the axis floor."

# ============================================================ F2a  anfevirus
da <- perf |> filter(muestra == M_ANPHE, ref == "AeAV_MERIDA__own1")

fa <- ggplot(da, aes(pos_ini, prof)) +
  geom_ribbon(aes(ymin = PISO, ymax = prof), fill = PAL[["anphevirus"]], alpha = 0.20) +
  geom_line(colour = PAL[["anphevirus"]], linewidth = 0.45) +
  escala_log(max(da$prof)) +
  scale_x_continuous(labels = kb, expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title    = "Replicating virus — Merida anphevirus",
    subtitle = paste0(M_ANPHE, " · own 12,945 bp reference (CheckV complete) · even coverage end to end"),
    x = "genome position", y = "mean depth (25 bp windows)",
    caption = NOTA_PISO
  ) +
  tema_perfil

ggsave(file.path(DIR_FIG, "F2a_anfevirus.png"), fa,
       width = 8.5, height = 3.2, dpi = 300, bg = "white")
ggsave(file.path(DIR_FIG, "F2a_anfevirus.pdf"), fa,
       width = 8.5, height = 3.2, device = cairo_pdf)

# ================================================================ F2b  PCLV
etiq_pclv <- c(
  PCLV__NC_038262 = "L — RdRp (6,806 bp)",
  PCLV__NC_038261 = "M — glycoprotein (3,935 bp)",
  PCLV__NC_038263 = "S — nucleocapsid (1,358 bp)"
)

db <- perf |>
  filter(muestra == M_PCLV, virus == "PCLV") |>
  mutate(seg = factor(etiq_pclv[ref], levels = unname(etiq_pclv)))

fb <- ggplot(db, aes(pos_ini, prof)) +
  geom_ribbon(aes(ymin = PISO, ymax = prof), fill = PAL[["pclv"]], alpha = 0.20) +
  geom_line(colour = PAL[["pclv"]], linewidth = 0.45) +
  facet_wrap(~ seg, nrow = 1, scales = "free_x") +
  escala_log(max(db$prof)) +
  scale_x_continuous(labels = kb, breaks = scales::breaks_extended(4),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(
    title    = "Segmented virus — Phasi Charoen-like virus",
    subtitle = paste0(M_PCLV, " · the only positive sample of 22 · all three segments covered at low depth"),
    x = "segment position", y = "mean depth",
    caption = paste(NOTA_PISO,
                    "Short gaps in L are sampling dropouts at 1–17×, not absence.")
  ) +
  tema_perfil

ggsave(file.path(DIR_FIG, "F2b_pclv.png"), fb,
       width = 9, height = 3.4, dpi = 300, bg = "white")
ggsave(file.path(DIR_FIG, "F2b_pclv.pdf"), fb,
       width = 9, height = 3.4, device = cairo_pdf)

# ================================================================= F2c  HTV
dc  <- perf |> filter(muestra == M_HTV, ref == "HTV__PZ595460")
mxc <- max(dc$prof)
lx  <- 1560   # la anotacion termina aqui, a la izquierda del pico (empieza en ~1601)

fc <- ggplot(dc, aes(pos_ini, prof)) +
  geom_ribbon(aes(ymin = PISO, ymax = prof), fill = PAL[["artefacto"]], alpha = 0.20) +
  geom_line(colour = PAL[["artefacto"]], linewidth = 0.45) +
  annotate("segment", x = lx + 8, xend = 1608, y = mxc / 6, yend = mxc / 2.2,
           colour = "grey35", linewidth = 0.35,
           arrow = arrow(length = unit(5, "pt"), type = "closed")) +
  annotate("text", x = lx, y = mxc / 6, hjust = 1, vjust = 0.5,
           size = 3.1, colour = "grey25", lineheight = 1.1,
           label = "~55 bp window, 86% identity to a\ntandem repeat on AaegL5 chromosome 1") +
  escala_log(mxc) +
  scale_x_continuous(labels = function(x) paste0(x, " bp"),
                     limits = c(0, 1672),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(
    title    = "Host repeat, not a virus — HTV capsid reference",
    subtitle = paste0(M_HTV, " · up to 80,769 reads, all in one window · RdRp segment at zero in all 22 samples"),
    x = "reference position", y = "mean depth",
    caption = NOTA_PISO
  ) +
  tema_perfil

ggsave(file.path(DIR_FIG, "F2c_htv_artefacto.png"), fc,
       width = 8.5, height = 3.2, dpi = 300, bg = "white")
ggsave(file.path(DIR_FIG, "F2c_htv_artefacto.pdf"), fc,
       width = 8.5, height = 3.2, device = cairo_pdf)

# ===================================== S2  el anfevirus en las 22 muestras
orden <- perf |>
  filter(ref == "AeAV_MERIDA__own1") |>
  group_by(muestra) |>
  summarise(m = mean(prof_media), .groups = "drop") |>
  arrange(desc(m)) |>
  pull(muestra)

ds <- perf |>
  filter(ref == "AeAV_MERIDA__own1") |>
  mutate(muestra = factor(muestra, levels = orden))

fs <- ggplot(ds, aes(pos_ini, prof)) +
  geom_ribbon(aes(ymin = PISO, ymax = prof), fill = PAL[["anphevirus"]], alpha = 0.25) +
  geom_line(colour = PAL[["anphevirus"]], linewidth = 0.25) +
  facet_wrap(~ muestra, ncol = 4) +
  escala_log(max(ds$prof)) +
  scale_x_continuous(labels = function(x) paste0(x / 1000, "k"),
                     breaks = scales::breaks_extended(3)) +
  labs(
    title    = "Merida anphevirus coverage across all 22 samples",
    subtitle = "Ordered by mean depth. Uniform profiles are replicating infections; flat zero is true absence.",
    x = "genome position", y = "mean depth (25 bp windows)",
    caption = NOTA_PISO
  ) +
  tema_perfil +
  theme(strip.text = element_text(size = 8.5), axis.text = element_text(size = 7.5))

ggsave(file.path(DIR_FIG, "S2_anfevirus_22_muestras.png"), fs,
       width = 9.5, height = 11, dpi = 300, bg = "white")

message("escrito en ", DIR_FIG, "/ : F2a, F2b, F2c y S2")
