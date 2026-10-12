# tema.R — aspecto comun de todas las figuras.
#
# Paleta categorica validada para daltonismo (orden fijo: el color sigue a la
# entidad, nunca a su rango). Si hay mas de 8 categorias, se agrupan en "Otros".

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

PAL <- c(azul = "#2a78d6", naranja = "#eb6834", aqua = "#1baf7a", amarillo = "#eda100",
         magenta = "#e87ba4", verde = "#008300", violeta = "#6250d6", rojo = "#e34948")
GRIS        <- "#b9b8b2"   # "Otros", no determinado, fuera de foco
TINTA       <- "#0b0b0b"
TINTA_2     <- "#52514e"
REJILLA     <- "#e6e5e1"
SUPERFICIE  <- "#fcfcfb"

# rampa secuencial (un solo tono, claro -> oscuro)
AZUL_SEQ <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
# escalones ordinales: el mas claro sigue siendo visible sobre el fondo
AZUL_ORD <- c("#86b6ef", "#3987e5", "#256abf", "#184f95")

tema_viroma <- function(base = 10) {
  theme_minimal(base_size = base) +
    theme(
      plot.background    = element_rect(fill = SUPERFICIE, colour = NA),
      panel.background   = element_rect(fill = SUPERFICIE, colour = NA),
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = REJILLA, linewidth = 0.3),
      text               = element_text(colour = TINTA),
      axis.text          = element_text(colour = TINTA_2),
      axis.title         = element_text(colour = TINTA_2),
      plot.title         = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle      = element_text(colour = TINTA_2),
      plot.caption       = element_text(colour = TINTA_2, hjust = 0, size = rel(0.8)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      legend.position    = "top",
      legend.justification = "left",
      legend.title       = element_blank(),
      legend.key.size    = unit(3.5, "mm"),
      strip.text         = element_text(face = "bold", hjust = 0)
    )
}
theme_set(tema_viroma())

# Guarda PDF (para la tesis) y PNG (para presentaciones) con el mismo nombre.
# El PDF va por Cairo: el dispositivo pdf() por defecto no dibuja simbolos como >= o x.
DISPOSITIVO_PDF <- if (capabilities("cairo")) cairo_pdf else "pdf"

guardar <- function(p, nombre, ancho = 180, alto = 110) {
  dir.create("figuras", showWarnings = FALSE)
  ggsave(file.path("figuras", paste0(nombre, ".pdf")), p, width = ancho, height = alto, units = "mm",
         device = DISPOSITIVO_PDF)
  ggsave(file.path("figuras", paste0(nombre, ".png")), p, width = ancho, height = alto, units = "mm",
         dpi = 300, bg = SUPERFICIE)
  invisible(p)
}

millones <- label_number(scale = 1e-6, suffix = " M", accuracy = 1)
