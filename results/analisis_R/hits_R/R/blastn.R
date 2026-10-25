# R/blastn.R — tablas y heatmap de los resultados de blastn (contigs de la union
# contra viral_custom/blast/viral_nt). Se apoya en R/hits.R (curacion, etiquetas).
#
# Entradas (datos_hits_vN/, desde v2):
#   blastn.tsv          salida cruda de annotate_votus.sh, sin encabezado: TODOS los HSP
#                       (megablast; dc-megablast solo para los que no dieron hit)
#   annot_por_votu.tsv  mejor hit blastn + DIAMOND por contig, con familia y realm
#   hits_por_contig.tsv grupo (ARN / ADN / fago ...) y en_genomad por contig

source("R/hits.R")

COLS_BLASTN <- c("contig", "accession", "pident", "aln_length", "qlen", "slen",
                 "evalue", "bitscore", "taxid", "species", "qcovs")
QCOV_NOMBRE <- 30   # cobertura minima de la query para dar nombre (igual que integrar_arbitro.py)
ID_ESPECIE  <- 90   # identidad nt minima para nombre de especie
FAGO  <- c("Steitzviridae","Fiersviridae","Solspiviridae","Duinviridae","Kyanoviridae",
           "Straboviridae","Autographiviridae","Microviridae","Inoviridae")
RETRO <- c("Retroviridae","Metaviridae","Belpaoviridae","Pseudoviridae")

partir_contig <- function(x) {
  p <- str_split_fixed(x, "__", 3)
  tibble(muestra = p[, 1], ensamblador = p[, 2])
}

leer_blastn <- function(dir) {
  f <- file.path(dir, "blastn.tsv"); a <- file.path(dir, "annot_por_votu.tsv")
  if (!file.exists(f) || !file.exists(a))
    stop("Faltan blastn.tsv / annot_por_votu.tsv en ", dir,
         ": congela con congelar_hits.sh v2 (la v1 no los trae)")
  crudo <- read_tsv(f, col_names = FALSE, col_types = cols(.default = "c"), na = c("", "NA"))
  if (ncol(crudo) != length(COLS_BLASTN))
    stop("blastn.tsv tiene ", ncol(crudo), " columnas; se esperaban ", length(COLS_BLASTN))
  names(crudo) <- COLS_BLASTN
  crudo <- crudo |>
    mutate(across(c(pident, aln_length, qlen, slen, evalue, bitscore, qcovs), as.numeric))
  annot <- read_tsv(a, col_types = cols(.default = "c"), na = c("", "NA", "-"))
  list(crudo = crudo, annot = annot)
}

# familia -> realm, a partir de las propias anotaciones (como en el awk de la farm)
mapa_realm <- function(annot) {
  annot |> filter(!is.na(familia), !is.na(realm)) |> distinct(familia, realm) |>
    group_by(familia) |> slice(1) |> ungroup() |> tibble::deframe()
}

clasificar <- function(fam, especie, realm_de) {
  r <- unname(realm_de[fam])
  case_when(
    is.na(fam)                                           ~ "sin_familia",
    fam %in% FAGO | str_detect(especie, "[Pp]hage")      ~ "fago",
    fam %in% RETRO                                       ~ "retroelemento",
    str_detect(coalesce(r, ""), "Riboviria") | fam == "Parvoviridae" ~ "candidato",
    is.na(r)                                             ~ "realm_desconocido",
    TRUE                                                 ~ "ADN_u_otro")
}

# Todos los HSP, con encabezado, muestra, familia y grupo.
tabla_blastn_todos <- function(b, d_hits) {
  realm_de <- mapa_realm(b$annot)
  fam_tax  <- b$annot |> filter(!is.na(bn_taxid), !is.na(bn_familia)) |>
    distinct(taxid = bn_taxid, family = bn_familia) |> group_by(taxid) |> slice(1) |> ungroup()
  gen <- select(d_hits, contig, en_genomad)
  b$crudo |>
    bind_cols(partir_contig(b$crudo$contig)) |>
    left_join(fam_tax, by = "taxid") |>
    left_join(gen, by = "contig") |>
    mutate(group = clasificar(family, species, realm_de),
           en_genomad = coalesce(en_genomad, FALSE)) |>
    group_by(contig) |> mutate(hit_rank = rank(-bitscore, ties.method = "first")) |> ungroup() |>
    relocate(muestra, ensamblador, contig, hit_rank) |>
    arrange(str_rank(muestra, numeric = TRUE), contig, hit_rank)
}

# Mejor hit por contig (el de annot_por_votu.tsv, mayor bitscore).
tabla_blastn_mejor <- function(todos) {
  todos |>
    filter(hit_rank == 1) |>
    mutate(naming = case_when(qcovs < QCOV_NOMBRE ~ "short hit (< 30% cover): not named",
                              pident >= ID_ESPECIE ~ "species",
                              TRUE                 ~ "relative (-like)"),
           virus = case_when(qcovs < QCOV_NOMBRE ~ NA_character_,
                             pident >= ID_ESPECIE ~ species,
                             TRUE                 ~ etiqueta(species, "pariente"))) |>
    select(muestra, ensamblador, contig, contig_length = qlen, accession, species, virus,
           family, group, pident, qcovs, aln_length, slen, evalue, bitscore, taxid,
           naming, en_genomad)
}

estado_de <- function(especie, n_genomad) {
  factor(case_when(especie %in% names(DERIVADO_HUESPED) ~ "Host-derived",
                   especie %in% names(NOTAS)            ~ "Unconfirmed",
                   n_genomad > 0                        ~ "Confirmed",
                   TRUE                                 ~ "Unconfirmed"),
         c("Confirmed", "Unconfirmed", "Host-derived"))
}

# Resumen por especie de blastn: solo virus de ARN con hit nombrable (>= 30 % cobertura).
tabla_blastn_especies <- function(mejor) {
  mejor |>
    filter(group == "candidato", qcovs >= QCOV_NOMBRE) |>
    group_by(species, family) |>
    summarise(n_samples = n_distinct(muestra), n_contigs = n(),
              n_genomad = sum(en_genomad),
              max_pident = max(pident), median_pident = median(pident),
              max_qcovs = max(qcovs), max_contig_length = max(contig_length),
              best_accession = accession[which.max(bitscore)], .groups = "drop") |>
    mutate(virus  = if_else(max_pident >= ID_ESPECIE, species, etiqueta(species, "pariente")),
           status = estado_de(species, n_genomad),
           note   = coalesce(DERIVADO_HUESPED[species], NOTAS[species], "")) |>
    relocate(virus, species, family, status) |>
    arrange(desc(n_samples), desc(n_contigs), desc(max_pident))
}

COLORES_ESTADO_TXT <- c(Confirmed = "#1F7A6F", Unconfirmed = "#9C6B00", `Host-derived` = "#6B7075")

figura_blastn_top <- function(mejor, esp, n_top = 20) {
  top <- slice_head(esp, n = n_top)
  m <- mejor |>
    filter(group == "candidato", qcovs >= QCOV_NOMBRE, species %in% top$species) |>
    group_by(species, muestra) |>
    summarise(pid = max(pident), n = n(), .groups = "drop") |>
    left_join(select(top, species, virus, status), by = "species")
  todas <- orden_muestras(mejor$muestra)
  m <- m |> mutate(virus = factor(virus, rev(top$virus)), muestra = factor(muestra, todas))
  ggplot(m, aes(muestra, virus, fill = pid)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = n, colour = pid > 85), size = 2.5, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
    scale_fill_viridis_c(name = "Max nt identity (%)", option = "mako", direction = -1,
                         limits = c(min(60, floor(min(m$pid))), 100)) +
    scale_x_discrete(drop = FALSE) +
    facet_grid(rows = vars(status), scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL,
         title = sprintf("Top %d blastn hits in assembled contigs", nrow(top)),
         subtitle = sprintf("RNA virus families, best hit per contig with >= %d%% query cover. Colour: best nt identity in the pool; number: contigs", QCOV_NOMBRE)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank(), legend.position = "top",
          strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
          plot.title.position = "plot")
}
