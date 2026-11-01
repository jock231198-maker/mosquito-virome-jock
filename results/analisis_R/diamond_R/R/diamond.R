# R/diamond.R — lectura y lógica de las tablas DIAMOND sobre los 6,004 contigs
# de la unión. Lo usa diamond_tablas.R; aquí no se dibuja nada.
#
# Las tres corridas son sobre los MISMOS 6,004 contigs (>= 1,000 pb, 22 pools,
# 5 ensambladores); solo cambian la base y la sensibilidad (ver BASES).
# El veredicto viral/celular es el del árbitro de nr (arbiter_nr.py), no se
# recalcula aquí: se lee de veredicto_contigs.tsv.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
})

# =============================================================================
# Listas editables (curación). Cambiar aquí, no en el código de abajo.
# =============================================================================

# Veredictos corregidos a mano. El árbitro ignora los hits sin taxid; en este
# contig el mejor hit sin taxid es una proteína de Aedes que gana al viral.
CORREGIR_VEREDICTO <- tribble(
  ~contig, ~veredicto, ~motivo,
  "M1_S1__rnaviral__NODE_255_length_1217_cov_3.691781", "celular_Artropodo",
  "best hit without taxid is an Aedes protein (133 vs 99 bits): host gene, NIRVS candidate"
)

# Mejores hits que son, en realidad, el mismo virus: se agrupan bajo una etiqueta.
# Los contigs del nudivirus de M71_S19 reparten su mejor hit entre varios
# nudivirus (y un granulovirus) según la proteína.
AGRUPAR <- c(
  # mismo virus con dos nombres en NCBI
  "Aedes aegypti anphevirus"     = "Aedes anphevirus",
  "Chaq virus"                   = "Chaq-like virus",
  # nudivirus de M71_S19 (+ M7_S15)
  "Tomelloso virus"              = "Novel nudivirus",
  "Oryctes rhinoceros nudivirus" = "Novel nudivirus",
  "Kallithea virus"              = "Novel nudivirus",
  "Esparto virus"                = "Novel nudivirus",
  "Mauternbach virus"            = "Novel nudivirus",
  "XiangYun nudi-like virus"     = "Novel nudivirus",
  "Drosophila innubila nudivirus"= "Novel nudivirus",
  "Apis mellifera nudivirus"     = "Novel nudivirus",
  "Agrotis segetum granulovirus" = "Novel nudivirus"
)

# Familias con muchos elementos virales integrados (EVE/NIRVS) en Ae. aegypti.
# Un virus de estas familias con identidad aa < UMBRAL_EVE se marca como
# candidato a EVE hasta pasar el cribado contra AaegL5 (cribado_genoma.sh).
FAMILIAS_EVE <- c("Flaviviridae", "Rhabdoviridae", "Chuviridae", "Phenuiviridae",
                  "Partitiviridae", "Orthototiviridae", "Totiviridae", "Xinmoviridae",
                  "Sedoreoviridae", "Spinareoviridae", "Phasmaviridae", "Orthomyxoviridae")
UMBRAL_EVE <- 90

# Virus de bacterias: se cuentan como virales pero no son virus del mosquito.
PATRON_FAGO <- "phage|Caudoviricetes|Andhaxevirus|Microviridae|Inoviridae|Leviviricetes"

# Clases por identidad aa del mejor hit viral en nr (mediana por especie)
UMBRAL_CONOCIDO  <- 95   # >= : virus ya descrito, se usa su nombre
UMBRAL_PARIENTE  <- 70   # entre los dos: pariente ("X-like"); por debajo: divergente

# =============================================================================
# Descripción de las corridas (va a la tabla T1 y a los pies de figura)
# =============================================================================
BASES <- tribble(
  ~base, ~etiqueta, ~origen, ~base_datos, ~parametros, ~e_max,
  "viraldb", "Own viral DB",
  "results/MERIDA/diamond_TODOS_vcustom.tsv",
  "RefSeq viral + GenBank viruses with arthropod host (refs/viral_custom/20260925, 680,459 proteins). Viral only; NCBI taxonomy",
  "blastx --very-sensitive -e 1e-3 -k 25", 1e-3,
  "rvdb", "RVDB v32.0",
  "results/MERIDA/diamond_TODOS_rvdb.tsv",
  "Reference Viral DataBase v32.0 (U-RVDBv32.0-prot_unique). Viral only, includes retroelements and giant viruses; no taxids",
  "blastx --very-sensitive, e-value <= 1e-3", 1e-3,
  "nr", "nr (jr46)",
  "results/MERIDA/contig_annot/contigs_union_MERIDA_nr_jr46/diamond.tsv",
  "NCBI nr, nr_tax.dmnd shared by jr46 (707,028,945 sequences; same content as the Feb 2024 nr). Viral + cellular; NCBI taxonomy",
  "blastx --sensitive -e 1e-5 -k 25; then arbiter (viral if best viral bitscore >= 1.1 x best cellular)", 1e-5
)

ENSAMBLADORES <- c(rnaviral = "rnaviralSPAdes", megahit = "MEGAHIT", metaviral = "metaviralSPAdes",
                   rna = "rnaSPAdes", trinity = "Trinity")

NIVELES_VER <- c("Viral", "Arthropod", "Other cellular", "Ambiguous", "No hit in nr")
GRUPOS_CEL <- c(Artropodo = "Arthropod", Metazoo_otro = "Other metazoan", Bacteria = "Bacteria",
                Hongo = "Fungi", Planta = "Plant", Eucariota_otro = "Other eukaryote",
                Arquea = "Archaea", otro = "Other")
NIVELES_CLASE <- c("Known", "Related", "Divergent")

# =============================================================================
# Lectura
# =============================================================================
NUM <- c("pident", "length", "qlen", "slen", "qcovhsp", "scovhsp", "evalue", "bitscore")

ultima_version <- function(base = ".") {
  d <- list.files(base, pattern = "^datos_diamond_v[0-9]+$")
  if (!length(d)) stop("no hay ninguna carpeta datos_diamond_vN en ", normalizePath(base))
  file.path(base, d[order(as.integer(sub("datos_diamond_v", "", d)), decreasing = TRUE)][1])
}

comprobar <- function(df, cols, que) {
  falta <- setdiff(cols, names(df))
  if (length(falta)) stop(que, ": faltan columnas ", paste(falta, collapse = ", "))
  df
}

con_origen <- function(df, col = "qseqid") {
  p <- str_split_fixed(df[[col]], "__", 3)
  df$muestra <- p[, 1]; df$ensamblador <- p[, 2]
  df
}

leer_hits <- function(f) {
  read_tsv(f, col_types = cols(.default = col_character()), quote = "", progress = FALSE) |>
    mutate(across(any_of(NUM), as.numeric)) |>
    con_origen()
}

leer_datos <- function(dir) {
  message("datos: ", normalizePath(dir))
  d <- list(
    viraldb = leer_hits(file.path(dir, "hits_viraldb.tsv.gz")),
    rvdb    = leer_hits(file.path(dir, "hits_rvdb.tsv.gz")),
    nr      = leer_hits(file.path(dir, "hits_nr.tsv.gz")),
    ver     = read_tsv(file.path(dir, "veredicto_contigs.tsv"), col_types = cols(.default = col_character()),
                       quote = "", progress = FALSE),
    ann     = read_tsv(file.path(dir, "annot_viraldb.tsv"), col_types = cols(.default = col_character()),
                       quote = "", progress = FALSE)
  )
  comprobar(d$viraldb, c("qseqid", "pident", "bitscore", "sscinames", "familia", "grupo"), "hits_viraldb")
  comprobar(d$rvdb, c("qseqid", "pident", "bitscore", "stitle"), "hits_rvdb")
  comprobar(d$nr, c("qseqid", "pident", "bitscore", "sscinames", "stitle", "grupo", "familia"), "hits_nr")
  comprobar(d$ver, c("votu", "longitud", "veredicto", "nr_grupo", "nr_nombre", "bits_mejor_viral", "bits_mejor_celular", "celular_nombre"), "veredicto")
  comprobar(d$ann, c("votu", "categoria", "familia", "realm"), "annot_viraldb")
  if (nrow(d$ver) != 6004) warning("veredicto con ", nrow(d$ver), " contigs (se esperaban 6004)")
  d$ver <- corregir_veredicto(d$ver)
  d
}

corregir_veredicto <- function(ver) {
  i <- match(CORREGIR_VEREDICTO$contig, ver$votu)
  if (anyNA(i)) warning("CORREGIR_VEREDICTO: contigs que no están: ",
                        paste(CORREGIR_VEREDICTO$contig[is.na(i)], collapse = ", "))
  ok <- !is.na(i)
  ver$veredicto_original <- ver$veredicto
  ver$veredicto[i[ok]] <- CORREGIR_VEREDICTO$veredicto[ok]
  ver
}

# =============================================================================
# Lógica
# =============================================================================
mejor_hit <- function(h) h |> arrange(qseqid, desc(bitscore), evalue) |> distinct(qseqid, .keep_all = TRUE)

organismo_rvdb <- function(stitle) str_match(stitle, "\\[([^\\[\\]]+)\\]\\s*$")[, 2]

# "YP_009553427.1 Ac68-like protein [Tomelloso virus]" -> "Ac68-like protein"
proteina <- function(stitle) {
  x <- str_remove(stitle, "^([A-Z]{1,4}_?[0-9]{4,}\\.[0-9]+\\s+)+")        # accesiones iniciales
  x <- str_remove_all(x, "\\s*\\[[^\\]]*\\]")
  str_squish(str_remove(x, "^(MAG|TPA):\\s*"))
}

clase_veredicto <- function(v) {
  factor(case_when(
    v == "viral" ~ "Viral",
    v == "ambiguo" ~ "Ambiguous",
    v == "celular_Artropodo" ~ "Arthropod",
    str_starts(v, "celular_") ~ "Other cellular",
    TRUE ~ "No hit in nr"), levels = NIVELES_VER)
}

clase_identidad <- function(pid) {
  factor(case_when(pid >= UMBRAL_CONOCIDO ~ "Known", pid >= UMBRAL_PARIENTE ~ "Related",
                   TRUE ~ "Divergent"), levels = NIVELES_CLASE)
}

# Nombre para figuras: el del virus si es conocido, "X-like" si es pariente o divergente
etiqueta_virus <- function(nombre, pid_med) {
  agr <- unname(AGRUPAR[nombre])
  like <- !is.na(pid_med) & pid_med < UMBRAL_CONOCIDO &
    !str_detect(nombre, "(\\bsp\\.|-like virus|-like)$")
  coalesce(agr, if_else(like, paste0(nombre, "-like"), nombre))
}

orden_muestra <- function(m) as.integer(str_extract(m, "(?<=^M)[0-9]+"))

# Una fila por contig (6,004): veredicto + mejor hit en cada base
tabla_contigs <- function(d) {
  v <- d$ver |>
    transmute(contig = votu, longitud = as.integer(longitud), veredicto, veredicto_original,
              bits_viral = suppressWarnings(as.numeric(bits_mejor_viral)),
              bits_celular = suppressWarnings(as.numeric(bits_mejor_celular)),
              nr_grupo, nr_mejor_nombre = nr_nombre, celular_nombre) |>
    con_origen("contig")
  vdb <- mejor_hit(d$viraldb) |>
    transmute(contig = qseqid, vdb_pident = pident, vdb_qcov = qcovhsp, vdb_nombre = sscinames,
              vdb_familia = na_if(familia, ""))
  rv <- mejor_hit(d$rvdb) |>
    transmute(contig = qseqid, rvdb_pident = pident, rvdb_qcov = qcovhsp, rvdb_nombre = organismo_rvdb(stitle))
  nv <- mejor_hit(filter(d$nr, grupo == "Virus")) |>
    transmute(contig = qseqid, nrv_pident = pident, nrv_qcov = qcovhsp, nrv_nombre = sscinames,
              nrv_familia = na_if(familia, ""), nrv_realm = na_if(realm, ""),
              nrv_proteina = proteina(stitle), nrv_titulo = stitle)
  an <- d$ann |>
    transmute(contig = votu, categoria_viraldb = categoria, familia_viraldb = na_if(familia, ""),
              realm_viraldb = na_if(realm, ""))
  ct <- v |> left_join(vdb, "contig") |> left_join(rv, "contig") |> left_join(nv, "contig") |>
    left_join(an, "contig") |>
    mutate(clase_ver = clase_veredicto(veredicto),
           grupo_celular = unname(GRUPOS_CEL[str_remove(veredicto, "^celular_")]),
           familia = coalesce(nrv_familia, familia_viraldb, "Unassigned"),
           fago = coalesce(str_detect(nrv_nombre, regex(PATRON_FAGO, ignore_case = TRUE)), FALSE),
           ens = factor(ENSAMBLADORES[ensamblador], levels = ENSAMBLADORES))
  if (nrow(ct) != nrow(d$ver)) stop("tabla_contigs duplicó filas")
  # clase y etiqueta POR ESPECIE (mediana de identidad de sus contigs virales)
  esp <- ct |> filter(veredicto == "viral") |> group_by(nrv_nombre) |>
    summarise(pid_esp = median(nrv_pident), .groups = "drop") |>
    mutate(etiqueta = etiqueta_virus(nrv_nombre, pid_esp))
  ct |> left_join(esp, "nrv_nombre") |>
    mutate(etiqueta = if_else(veredicto == "viral", etiqueta, NA_character_),
           pid_esp = if_else(veredicto == "viral", pid_esp, NA_real_))
}

moda <- function(x) { x <- x[!is.na(x) & x != "Unassigned"]; if (!length(x)) "Unassigned" else names(sort(table(x), decreasing = TRUE))[1] }

# Catálogo: una fila por virus (etiqueta) a partir de los contigs virales
catalogo <- function(ct) {
  ct |> filter(veredicto == "viral") |>
    group_by(etiqueta) |>
    summarise(familia = moda(familia),
              realm = moda(coalesce(nrv_realm, realm_viraldb)),
              best_hit_names = paste(sort(unique(nrv_nombre)), collapse = "; "),
              contigs = n(),
              pools = n_distinct(muestra),
              pool_list = paste(unique(muestra[order(orden_muestra(muestra))]), collapse = ", "),
              aa_id_median = round(median(nrv_pident), 1),
              aa_id_max = max(nrv_pident),
              qcov_median = round(median(nrv_qcov), 1),
              longest_contig_bp = max(longitud),
              total_bp = sum(longitud),
              fago = any(fago), .groups = "drop") |>
    mutate(clase = clase_identidad(aa_id_median),
           eve_candidato = !fago & aa_id_median < UMBRAL_EVE & familia %in% FAMILIAS_EVE,
           tipo = case_when(fago ~ "Bacteriophage",
                            eve_candidato ~ "EVE candidate (pending AaegL5 screen)",
                            TRUE ~ "Virus")) |>
    arrange(desc(pools), desc(contigs))
}
