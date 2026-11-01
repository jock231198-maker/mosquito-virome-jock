#!/usr/bin/env Rscript
# diamond_tablas.R — tablas y figuras de DIAMOND sobre los 6,004 contigs de la unión
#
# Uso (en casa, desde results/analisis_R/diamond_R):
#   Rscript diamond_tablas.R                      # usa la última datos_diamond_vN
#   Rscript diamond_tablas.R datos_diamond_v1     # una versión concreta
#
# Salida:
#   tablas/T1..T7_*.tsv  y  tablas/diamond_tables.xlsx (una hoja por tabla + README)
#   figuras/D1..D6_*.{pdf,png}
# Textos de salida en inglés, convención "X-like" para parientes.

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) setwd(dirname(normalizePath(f)))
})
source("R/diamond.R")
suppressPackageStartupMessages({ library(ggplot2); library(scales); library(forcats) })

args <- commandArgs(trailingOnly = TRUE)
DATOS <- if (length(args)) args[1] else ultima_version(".")
d  <- leer_datos(DATOS)
ct <- tabla_contigs(d)
cat_v <- catalogo(ct)
dir.create("tablas", showWarnings = FALSE); dir.create("figuras", showWarnings = FALSE)
N <- nrow(ct)
fmt <- function(x) format(x, big.mark = ",", trim = TRUE)
pct <- function(x, n) sprintf("%.1f%%", 100 * x / n)
version_txt <- basename(normalizePath(DATOS))

# =============================================================================
# Tema y paleta (paleta validada para daltonismo; el color va por la función)
# =============================================================================
TINTA <- "#0b0b0b"; TINTA2 <- "#52514e"; REJILLA <- "#e6e5e0"; SUPERFICIE <- "#fcfcfb"
COL_VER <- c("Viral" = "#2a78d6", "Arthropod" = "#eb6834", "Other cellular" = "#1baf7a",
             "Ambiguous" = "#eda100", "No hit in nr" = "#c9c8c2")
SECUENCIAL <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
GRIS <- "#b5b4ae"

theme_set(
  theme_minimal(base_size = 9) +
    theme(text = element_text(colour = TINTA),
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(colour = TINTA2, margin = margin(b = 6)),
          plot.caption = element_text(colour = TINTA2, size = 7, hjust = 0, margin = margin(t = 8)),
          plot.title.position = "plot", plot.caption.position = "plot",
          axis.text = element_text(colour = TINTA2), axis.title = element_text(colour = TINTA2),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = REJILLA, linewidth = 0.3),
          strip.text = element_text(face = "bold", hjust = 0, colour = TINTA),
          legend.position = "top", legend.justification = "left",
          legend.title = element_text(colour = TINTA2), legend.key.size = unit(9, "pt"),
          plot.background = element_rect(fill = SUPERFICIE, colour = NA),
          plot.margin = margin(8, 12, 8, 8))
)

guardar <- function(p, nombre, ancho = 180, alto = 100) {
  png <- file.path("figuras", paste0(nombre, ".png"))
  pdf <- file.path("figuras", paste0(nombre, ".pdf"))
  ggsave(png, p, width = ancho, height = alto, units = "mm", dpi = 300, bg = SUPERFICIE,
         device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png")
  ggsave(pdf, p, width = ancho, height = alto, units = "mm", bg = SUPERFICIE,
         device = if (capabilities("cairo")) cairo_pdf else "pdf")
  message("  ", png)
}
pie <- function(...) str_wrap(paste0(..., " Data: ", version_txt, "."), width = 150)

# =============================================================================
# TABLAS
# =============================================================================
message("== tablas")

# --- T1: corridas ------------------------------------------------------------
t1 <- BASES |>
  rowwise() |>
  mutate(query = paste0(fmt(N), " union contigs (>= 1,000 bp; 22 pools; 5 assemblers)"),
         hits = nrow(d[[base]]),
         contigs_with_hit = n_distinct(d[[base]]$qseqid),
         pct_with_hit = round(100 * contigs_with_hit / N, 1),
         max_evalue_observed = signif(max(d[[base]]$evalue), 3),
         viral_by_nr_arbiter = sum(ct$clase_ver[ct$contig %in% d[[base]]$qseqid] == "Viral")) |>
  ungroup() |>
  transmute(database = etiqueta, database_description = base_datos, query, parameters = parametros,
            hits, contigs_with_hit, pct_with_hit, viral_by_nr_arbiter, max_evalue_observed,
            source_file = origen)

# --- T2: veredictos ----------------------------------------------------------
t2 <- ct |>
  mutate(verdict = case_when(veredicto == "viral" ~ "Viral",
                             veredicto == "ambiguo" ~ "Ambiguous",
                             veredicto == "sin_hit_nr" ~ "No hit in nr",
                             TRUE ~ paste("Cellular:", grupo_celular))) |>
  count(verdict, name = "contigs") |>
  mutate(pct = round(100 * contigs / N, 1)) |>
  arrange(desc(contigs))

# --- T3: bases solo virales frente al árbitro --------------------------------
cob <- bind_rows(lapply(BASES$base, function(b) tibble(base = b, contig = unique(d[[b]]$qseqid)))) |>
  left_join(select(ct, contig, clase_ver), "contig") |>
  mutate(database = factor(BASES$etiqueta[match(base, BASES$base)], levels = BASES$etiqueta))
t3 <- cob |> count(database, clase_ver) |>
  pivot_wider(names_from = clase_ver, values_from = n, values_fill = 0) |>
  mutate(contigs_with_hit = rowSums(across(any_of(NIVELES_VER))),
         pct_viral = round(100 * Viral / contigs_with_hit, 1)) |>
  select(database, contigs_with_hit, any_of(NIVELES_VER), pct_viral)

# --- T4: catálogo de virus ---------------------------------------------------
t4 <- cat_v |>
  transmute(virus = etiqueta, type = tipo, identity_class = clase, family = familia, realm,
            best_hit_names, pools, contigs, aa_id_median, aa_id_max, qcov_median,
            longest_contig_bp, total_bp, pool_list)

# --- T5: por pool ------------------------------------------------------------
virus_pool <- ct |> filter(veredicto == "viral") |>
  left_join(select(cat_v, etiqueta, tipo), "etiqueta") |>
  group_by(muestra) |>
  summarise(viruses = n_distinct(etiqueta[tipo == "Virus"]),
            eve_candidates = n_distinct(etiqueta[tipo != "Virus" & tipo != "Bacteriophage"]),
            .groups = "drop")
t5 <- ct |> count(muestra, clase_ver) |>
  pivot_wider(names_from = clase_ver, values_from = n, values_fill = 0) |>
  mutate(contigs = rowSums(across(any_of(NIVELES_VER)))) |>
  left_join(virus_pool, "muestra") |>
  mutate(across(c(viruses, eve_candidates), ~ coalesce(.x, 0L))) |>
  arrange(orden_muestra(muestra)) |>
  select(pool = muestra, contigs, any_of(NIVELES_VER), viruses, eve_candidates)

# --- T6: contigs virales -----------------------------------------------------
t6 <- ct |> filter(veredicto == "viral") |>
  left_join(select(cat_v, etiqueta, tipo), "etiqueta") |>
  arrange(orden_muestra(muestra), desc(longitud)) |>
  transmute(contig, pool = muestra, assembler = as.character(ens), length_bp = longitud,
            virus = etiqueta, type = tipo, nr_best_viral_hit = nrv_nombre, family = familia,
            aa_identity = nrv_pident, query_cov_hsp = nrv_qcov, protein = nrv_proteina,
            bits_viral, bits_cellular = bits_celular, best_cellular_hit = celular_nombre,
            own_viraldb_hit = vdb_nombre, own_viraldb_aa_id = vdb_pident,
            rvdb_hit = rvdb_nombre, rvdb_aa_id = rvdb_pident)

# --- T7: correcciones aplicadas ----------------------------------------------
t7 <- CORREGIR_VEREDICTO |> rename(new_verdict = veredicto, reason = motivo) |>
  left_join(select(ct, contig, original_verdict = veredicto_original), "contig")

tablas <- list(T1_runs = t1, T2_verdicts = t2, T3_viraldb_vs_nr = t3, T4_virus_catalogue = t4,
               T5_per_pool = t5, T6_viral_contigs = t6, T7_corrections = t7)
for (n in names(tablas)) write_tsv(tablas[[n]], file.path("tablas", paste0(n, ".tsv")), na = "")

readme <- tibble(item = c(
  "Query", "Databases", "Arbiter", "Verdict classes", "Identity classes", "EVE candidate",
  "Bacteriophage", "Labels", "Corrections", "Counting", "Data version"),
  definition = c(
  paste0("All three DIAMOND blastx runs used the same ", fmt(N), " contigs: the dereplicated union of 5 assemblers (rnaviralSPAdes, rnaSPAdes, metaviralSPAdes, MEGAHIT, Trinity), >= 1,000 bp, 22 pools."),
  "Own viral DB (viral only, with taxonomy), RVDB v32.0 (viral only, no taxids) and nr (viral + cellular). See T1 for parameters; e-value cut-offs differ (1e-3 viral DBs, 1e-5 nr).",
  "arbiter_nr.py on the nr hits: viral if the best viral bitscore is >= 1.1 x the best cellular bitscore; ambiguous if within 10%; otherwise cellular (group from NCBI taxonomy). Hits without taxid are ignored.",
  "Viral, Arthropod, Other cellular (bacteria, other metazoa, fungi, plants, other eukaryotes, archaea), Ambiguous, No hit in nr.",
  paste0("Median aa identity of the best viral nr hit per virus: Known >= ", UMBRAL_CONOCIDO, "%, Related ", UMBRAL_PARIENTE, "-", UMBRAL_CONOCIDO, "%, Divergent < ", UMBRAL_PARIENTE, "%."),
  paste0("Virus in a family rich in endogenous viral elements of Ae. aegypti (", paste(FAMILIAS_EVE, collapse = ", "), ") with median aa identity < ", UMBRAL_EVE, "%. Not counted as a virus until screened against AaegL5."),
  "Best viral hit is a phage: viral, but a virus of bacteria, not of the mosquito.",
  "Species name when Known; 'X-like' otherwise. Hits that are the same virus are merged (AGRUPAR in R/diamond.R), e.g. the novel nudivirus.",
  "Manual verdict corrections listed in T7.",
  "Contigs are not viruses: one virus appears in several pools and assemblers. Use 'pools' for prevalence.",
  version_txt))
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(c(list(README = readme), tablas), file.path("tablas", "diamond_tables.xlsx"))
} else message("  (sin writexl: solo TSV)")
message("  tablas/ (", length(tablas), " TSV + xlsx)")

# =============================================================================
# FIGURAS
# =============================================================================
message("== figuras")

# --- D1: lo que encuentra cada base, según el árbitro de nr -------------------
d1 <- cob |> count(database, clase_ver) |>
  mutate(database = fct_rev(database))
tot1 <- d1 |> group_by(database) |>
  summarise(v = sum(n[clase_ver == "Viral"]), n = sum(n), .groups = "drop") |>
  mutate(lab = paste0(fmt(n), " contigs  |  ", fmt(v), " viral (", pct(v, n), ")"))
p1 <- ggplot(d1, aes(n, database, fill = clase_ver)) +
  geom_col(width = 0.62, colour = SUPERFICIE, linewidth = 0.5, position = position_stack(reverse = TRUE)) +
  geom_text(data = tot1, aes(n, database, label = lab), inherit.aes = FALSE,
            hjust = -0.04, size = 2.9, colour = TINTA) +
  scale_fill_manual(values = COL_VER, drop = FALSE, name = "Verdict of the nr arbiter") +
  scale_x_continuous(labels = label_comma(), expand = expansion(c(0, 0.5))) +
  labs(title = "Most hits in viral-only databases are not viruses",
       subtitle = paste0("Contigs with a DIAMOND hit in each database, out of ", fmt(N),
                         " union contigs, coloured by what nr says they are"),
       x = "contigs with a hit", y = NULL,
       caption = pie("Same query for all three runs. Own viral DB and RVDB: blastx --very-sensitive, e <= 1e-3. ",
                     "nr (jr46): blastx --sensitive -e 1e-5 -k 25, then the arbiter (viral if best viral bitscore >= 1.1 x best cellular).")) +
  theme(panel.grid.major.y = element_blank()) +
  guides(fill = guide_legend(nrow = 1))
guardar(p1, "D1_databases_vs_arbiter", alto = 72)

# --- D2: identidad del mejor hit en las bases virales, por veredicto ---------
d2 <- bind_rows(
  ct |> transmute(base = "Own viral DB", pid = vdb_pident, clase_ver),
  ct |> transmute(base = "RVDB v32.0", pid = rvdb_pident, clase_ver)) |>
  filter(!is.na(pid), clase_ver %in% c("Viral", "Arthropod", "Other cellular")) |>
  mutate(grupo = factor(if_else(clase_ver == "Viral", "Viral according to nr", "Cellular according to nr"),
                        levels = c("Viral according to nr", "Cellular according to nr")),
         base = factor(base, levels = c("Own viral DB", "RVDB v32.0")))
med2 <- d2 |> group_by(base, grupo) |>
  summarise(m = median(pid), n = n(), .groups = "drop") |>
  mutate(lab = paste0("n = ", fmt(n), "\nmedian ", round(m), "% aa"))
p2 <- ggplot(d2, aes(pid, fill = grupo)) +
  geom_histogram(binwidth = 2.5, boundary = 0, colour = SUPERFICIE, linewidth = 0.25) +
  geom_vline(data = med2, aes(xintercept = m), colour = TINTA2, linewidth = 0.35, linetype = "22") +
  geom_text(data = med2, aes(x = 2, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = 0, vjust = 1.3, size = 2.6, colour = TINTA, lineheight = 0.95) +
  facet_grid(grupo ~ base, scales = "free_y") +
  scale_fill_manual(values = c("Viral according to nr" = COL_VER[["Viral"]],
                               "Cellular according to nr" = COL_VER[["Arthropod"]]), guide = "none") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), labels = \(x) paste0(x, "%")) +
  labs(title = "Identity of viral-database hits for contigs nr calls viral or cellular",
       subtitle = "aa identity of the best hit in each viral-only database, split by the nr verdict",
       x = "aa identity of the best hit", y = "contigs",
       caption = pie("Ambiguous contigs and contigs without hit in nr are left out. Dashed line: median.")) +
  theme(strip.text.y = element_text(angle = 0))
guardar(p2, "D2_identity_by_verdict", alto = 95)

# --- D3: veredicto por ensamblador -------------------------------------------
d3 <- ct |> count(ens, clase_ver) |>
  group_by(ens) |> mutate(f = n / sum(n)) |> ungroup() |>
  mutate(ens = fct_rev(ens))
tot3 <- ct |> group_by(ens) |>
  summarise(n = n(), v = sum(clase_ver == "Viral"), .groups = "drop") |>
  mutate(ens = factor(ens, levels = levels(d3$ens)),
         lab = paste0(fmt(n), " contigs  |  ", fmt(v), " viral (", pct(v, n), ")"))
p3 <- ggplot(d3, aes(f, ens, fill = clase_ver)) +
  geom_col(width = 0.62, colour = SUPERFICIE, linewidth = 0.5, position = position_stack(reverse = TRUE)) +
  geom_text(data = tot3, aes(1, ens, label = lab), inherit.aes = FALSE,
            hjust = -0.04, size = 2.9, colour = TINTA) +
  scale_fill_manual(values = COL_VER, drop = FALSE, name = "Verdict of the nr arbiter") +
  scale_x_continuous(labels = label_percent(), expand = expansion(c(0, 0.55)),
                     breaks = seq(0, 1, 0.25)) +
  labs(title = "Where the viral contigs come from",
       subtitle = paste0("Verdict of the ", fmt(N), " union contigs, by assembler"),
       x = "share of the assembler's contigs", y = NULL,
       caption = pie("Contigs after dereplication of the five assemblies (CD-HIT-EST), >= 1,000 bp.")) +
  theme(panel.grid.major.y = element_blank()) +
  guides(fill = guide_legend(nrow = 1))
guardar(p3, "D3_verdict_by_assembler", alto = 75)

# --- D4: heatmap virus x pool ------------------------------------------------
POOLS <- unique(ct$muestra); POOLS <- POOLS[order(orden_muestra(POOLS))]
panel_de <- function(clase, tipo) {
  factor(if_else(tipo == "Virus", paste0(as.character(clase)), "EVE candidate"),
         levels = c(NIVELES_CLASE, "EVE candidate"))
}
cat4 <- cat_v |> filter(tipo != "Bacteriophage") |>
  mutate(panel = panel_de(clase, tipo),
         fila = paste0(etiqueta, "  (", pools, ")"))
d4 <- ct |> filter(veredicto == "viral", etiqueta %in% cat4$etiqueta) |>
  group_by(etiqueta, muestra) |> summarise(pid = max(nrv_pident), n = n(), .groups = "drop") |>
  left_join(select(cat4, etiqueta, panel, fila, pools), "etiqueta") |>
  mutate(muestra = factor(muestra, levels = POOLS))
orden4 <- cat4 |> arrange(panel, pools, contigs) |> pull(fila)
d4$fila <- factor(d4$fila, levels = unique(orden4))
n_fagos <- sum(cat_v$tipo == "Bacteriophage")
p4 <- ggplot(d4, aes(muestra, fila, fill = pid)) +
  geom_tile(colour = SUPERFICIE, linewidth = 0.6) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradientn(colours = SECUENCIAL, limits = c(20, 100), breaks = c(20, 40, 60, 80, 100),
                       labels = \(x) paste0(x, "%"), name = "best aa identity (nr)") +
  scale_x_discrete(drop = FALSE, position = "top") +
  labs(title = "Viruses found by DIAMOND against nr, per pool",
       subtitle = paste0(fmt(nrow(cat4)), " viruses from ", fmt(sum(cat4$contigs)),
                         " viral contigs; number in brackets = pools"),
       x = NULL, y = NULL,
       caption = pie("Known: median aa identity >= ", UMBRAL_CONOCIDO, "%; Related: ", UMBRAL_PARIENTE, "-",
                     UMBRAL_CONOCIDO, "%; Divergent: < ", UMBRAL_PARIENTE, "%. EVE candidate: family rich in ",
                     "endogenous viral elements of Ae. aegypti and < ", UMBRAL_EVE, "% aa; pending the AaegL5 screen. ",
                     if (n_fagos) paste0(n_fagos, " bacteriophage(s) not shown. ") else "")) +
  theme(panel.grid = element_blank(),
        axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 7),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 1),
        legend.key.width = unit(28, "pt"))
guardar(p4, "D4_virus_by_pool_heatmap", ancho = 190, alto = 40 + 4.2 * nrow(cat4))

# --- D5: identidad frente a longitud -----------------------------------------
d5 <- ct |> filter(veredicto == "viral")
top_fam <- d5 |> filter(familia != "Unassigned") |> count(familia, sort = TRUE) |> slice_head(n = 3) |> pull(familia)
d5 <- d5 |> mutate(fam = factor(if_else(familia %in% top_fam, familia, "Other / unassigned"),
                                levels = c(top_fam, "Other / unassigned")))
col5 <- setNames(c("#2a78d6", "#eb6834", "#1baf7a")[seq_along(top_fam)], top_fam)
col5[["Other / unassigned"]] <- GRIS
p5 <- ggplot(d5, aes(longitud, nrv_pident)) +
  annotate("rect", xmin = 900, xmax = 20000, ymin = UMBRAL_CONOCIDO, ymax = 100.8, fill = "#cde2fb", alpha = 0.45) +
  annotate("rect", xmin = 900, xmax = 20000, ymin = 15, ymax = UMBRAL_PARIENTE, fill = "#f1f0ec", alpha = 0.7) +
  annotate("text", x = Inf, y = UMBRAL_CONOCIDO - 1, label = "Known", hjust = 1.05, vjust = 1, size = 2.7, colour = TINTA2) +
  annotate("text", x = Inf, y = UMBRAL_PARIENTE - 1, label = "Divergent", hjust = 1.05, vjust = 1, size = 2.7, colour = TINTA2) +
  geom_point(aes(fill = fam), shape = 21, size = 2.3, stroke = 0.4, colour = SUPERFICIE, alpha = 0.95) +
  scale_x_log10(labels = label_comma(), breaks = c(1000, 2000, 5000, 10000, 15000), expand = c(0, 0)) +
  coord_cartesian(xlim = c(950, 18000)) +
  scale_y_continuous(limits = c(15, 101), breaks = seq(20, 100, 20), labels = \(x) paste0(x, "%")) +
  scale_fill_manual(values = col5, name = NULL) +
  labs(title = "Known viruses sit at ~100% aa identity; divergent contigs spread below 70%",
       subtitle = paste0(fmt(nrow(d5)), " viral contigs: aa identity of the best viral hit in nr against contig length"),
       x = "contig length (bp, log scale)", y = "aa identity (nr)",
       caption = pie("Families from the NCBI taxonomy of the best viral hit; three most frequent shown.")) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3)))
guardar(p5, "D5_identity_vs_length", alto = 95)

# --- D6: el nudivirus --------------------------------------------------------
d6 <- ct |> filter(veredicto == "viral", etiqueta == "Novel nudivirus") |>
  arrange(muestra, longitud) |>
  mutate(y = factor(contig, levels = contig),
         lab = paste0(nrv_proteina, "  |  ", round(nrv_pident), "% aa  (", nrv_nombre, ")"))
if (nrow(d6)) {
  m6 <- unique(d6$muestra); col6 <- setNames(c("#2a78d6", "#eb6834", "#1baf7a")[seq_along(m6)], m6)
  kb <- round(sum(d6$longitud) / 1000, 1)
  p6 <- ggplot(d6, aes(longitud / 1000, y, fill = muestra)) +
    geom_col(width = 0.7, colour = SUPERFICIE, linewidth = 0.3) +
    geom_text(aes(label = lab), hjust = -0.03, size = 2.4, colour = TINTA) +
    scale_fill_manual(values = col6, name = NULL) +
    scale_x_continuous(expand = expansion(c(0, 1.6)), breaks = seq(0, 16, 4)) +
    scale_y_discrete(labels = NULL) +
    facet_grid(muestra ~ ., scales = "free_y", space = "free_y") +
    labs(title = paste0("A divergent nudivirus in ", paste(m6, collapse = " and ")),
         subtitle = paste0(nrow(d6), " viral contigs, ", kb, " kb in total (upper bound: assemblers overlap); ",
                           "label = best viral hit in nr"),
         x = "contig length (kb)", y = NULL,
         caption = pie("Core nudivirus/baculovirus proteins at low aa identity point to a new nudivirus. ",
                       "Check cruce_votus/ for whether any contig reached the vOTU catalogue.")) +
    theme(panel.grid.major.y = element_blank(), axis.ticks.y = element_blank(),
          strip.text.y = element_text(angle = 0)) +
    guides(fill = "none")
  guardar(p6, "D6_nudivirus_contigs", alto = 30 + 4 * nrow(d6))
}

# =============================================================================
message("\n== resumen")
message("  contigs: ", fmt(N), " | viral: ", sum(ct$clase_ver == "Viral"),
        " | ambiguous: ", sum(ct$clase_ver == "Ambiguous"),
        " | virus en catálogo: ", sum(cat_v$tipo == "Virus"),
        " | candidatos a EVE: ", sum(str_starts(cat_v$tipo, "EVE")),
        " | fagos: ", n_fagos)
message("  correcciones aplicadas: ", sum(ct$veredicto != ct$veredicto_original))
