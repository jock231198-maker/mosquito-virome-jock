#!/bin/bash
# recolectar_tablas.sh — reúne todas las tablas de identificación viral en un
# solo directorio, con un manifiesto que dice de dónde salió cada una.
#
# Uso:
#   ./recolectar_tablas.sh [destino]
#   destino por defecto: $RESULTS_DIR/tablas_identificacion
#
# No mueve nada: copia. Los originales se quedan donde están.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

DEST="${1:-$RESULTS_DIR/tablas_identificacion}"
mkdir -p "$DEST"
MAN="$DEST/MANIFIESTO.tsv"

printf 'destino\torigen\tlineas\tbytes\tque_es\n' > "$MAN"

# copiar(destino_relativo, ruta_origen, descripcion)
copiar () {
  local dst="$1" src="$2" desc="$3"
  if [[ -s "$src" ]]; then
    mkdir -p "$DEST/$(dirname "$dst")"
    cp -f "$src" "$DEST/$dst"
    printf '%s\t%s\t%s\t%s\t%s\n' "$dst" "$src" "$(wc -l < "$src")" "$(stat -c%s "$src")" "$desc" >> "$MAN"
    printf '  OK    %-52s %8s lineas\n' "$dst" "$(wc -l < "$src")"
  else
    printf '%s\t%s\tFALTA\tFALTA\t%s\n' "$dst" "$src" "$desc" >> "$MAN"
    printf '  FALTA %-52s %s\n' "$dst" "$src"
  fi
}

R="$RESULTS_DIR"
A="$R/votu_annot/votus_viraldb20260925"
B="$R/votu_annot/all_final_viraldb20260925"
G="$R/genomad_all/all_final"
Q="$R/qc_control"

echo "== 1. geNomad (entrada: union_all/all_final.fasta, 6004 contigs) =="
copiar genomad/virus_summary.tsv        "$G/all_final_summary/all_final_virus_summary.tsv"                       "1144 contigs virales: score, hallmarks, taxonomia"
copiar genomad/plasmid_summary.tsv      "$G/all_final_summary/all_final_plasmid_summary.tsv"                     "plasmidos (0 filas)"
copiar genomad/aggregated_scores.tsv    "$G/all_final_aggregated_classification/all_final_aggregated_classification.tsv" "scores de los 6004, virales y no virales"
copiar genomad/features.tsv             "$G/all_final_marker_classification/all_final_features.tsv"              "genes, hallmarks y rasgos por contig"
copiar genomad/genes.tsv                "$G/all_final_annotate/all_final_genes.tsv"                              "coordenadas y anotacion de cada gen"
copiar genomad/taxonomy.tsv             "$G/all_final_annotate/all_final_taxonomy.tsv"                           "asignacion taxonomica por contig"
copiar genomad/parametros_summary.json  "$G/all_final_summary/all_final_summary.json"                            "parametros efectivos de la corrida"

echo "== 2. CheckV =="
copiar checkv/contigs_quality_summary.tsv   "$R/checkv_contigs/all_final/quality_summary.tsv"                    "sobre los 6004 contigs (independiente de geNomad)"
copiar checkv/contigs_completeness.tsv      "$R/checkv_contigs/all_final/completeness.tsv"                       "completitud, 6004"
copiar checkv/contigs_contamination.tsv     "$R/checkv_contigs/all_final/contamination.tsv"                      "contaminacion, 6004"
copiar checkv/contigs_complete_genomes.tsv  "$R/checkv_contigs/all_final/complete_genomes.tsv"                   "genomas completos detectados"
copiar checkv/viral_quality_summary.tsv     "$R/checkv_viral/all_final_virus.fna/quality_summary.tsv"            "sobre los 1144 virales de geNomad"
copiar checkv/viral_completeness.tsv        "$R/checkv_viral/all_final_virus.fna/completeness.tsv"               "completitud, 1144"
copiar checkv/viral_contamination.tsv       "$R/checkv_viral/all_final_virus.fna/contamination.tsv"              "contaminacion, 1144"

echo "== 3. vOTUs y abundancia =="
copiar votus/votus.clstr            "$Q/votus.clstr"            "miembros de cada uno de los 748 clusters"
copiar votus/votus_master.tsv       "$Q/votus_master.tsv"       "tabla maestra: vOTU + score + muestras + CheckV + taxonomia"
copiar votus/contigs_union.tsv      "$Q/contigs_union.tsv"      "los 6004: muestra, ensamblador, longitud, cobertura"
copiar votus/contigs_union_cov.tsv  "$Q/contigs_union.tsv.cov"  "igual, con la cobertura de MEGAHIT recuperada"
copiar votus/votu_coverage.tsv      "$Q/votu_coverage.tsv"      "AMPLITUD por mapeo: % cubierto y profundidad (el bueno)"
copiar votus/votu_counts_long.tsv   "$Q/votu_counts_long.tsv"   "conteo de lecturas por vOTU (OJO: inflado por rRNA)"

echo "== 4. DIAMOND =="
copiar diamond/TODOS_vcustom.tsv     "$R/diamond_TODOS_vcustom.tsv"     "6004 contigs vs base viral propia (con sscinames)"
copiar diamond/TODOS_rvdb.tsv        "$R/diamond_TODOS_rvdb.tsv"        "6004 contigs vs RVDB v32 --very-sensitive"
copiar diamond/votus_rvdb_final.tsv  "$R/diamond_votus_rvdb_final.tsv"  "748 vOTUs vs RVDB, con qcovhsp y scovhsp"
copiar diamond/votus_rvdb_vs.tsv     "$R/diamond_votus_rvdb_vs.tsv"     "748 vOTUs vs RVDB (sin scovhsp, superado)"
copiar diamond/votus_rvdb.tsv        "$R/diamond_votus_rvdb.tsv"        "748 vOTUs vs RVDB, primera pasada (superado)"
copiar diamond/votus_nr.tsv          "$R/diamond_votus_nr.tsv"          "748 vOTUs vs nr (feb 2024), con sscinames"
copiar diamond/oscuras_nr.tsv        "$R/oscuras_nr.tsv"                "50 vOTUs sin hit vs nr --very-sensitive"
copiar diamond/candidata_nr.tsv      "$R/cand_nr.tsv"                   "candidata k119_1227 vs nr (0 hits)"

echo "== 5. Anotacion integrada contra la base viral propia =="
copiar anotacion/votus_annot.tsv     "$A/annot_por_votu.tsv"   "748 vOTUs: blastn + DIAMOND + familia/orden/realm"
copiar anotacion/votus_blastn.tsv    "$A/blastn.tsv"           "blastn crudo, 748 vOTUs"
copiar anotacion/votus_diamond.tsv   "$A/diamond.tsv"          "blastx crudo, 748 vOTUs"
copiar anotacion/votus_resumen.txt   "$A/resumen.txt"          "recuento por categoria y familia"
copiar anotacion/contigs_annot.tsv   "$B/annot_por_votu.tsv"   "6004 contigs: blastn + DIAMOND + linaje"
copiar anotacion/contigs_blastn.tsv  "$B/blastn.tsv"           "blastn crudo, 6004 contigs"
copiar anotacion/contigs_diamond.tsv "$B/diamond.tsv"          "blastx crudo, 6004 contigs"
copiar anotacion/contigs_resumen.txt "$B/resumen.txt"          "recuento por categoria y familia"

echo "== 6. Kraken2 / Bracken =="
for d in kraken2 kraken2_nohost; do
  if [[ -d "$R/$d" ]]; then
    mkdir -p "$DEST/kraken/$d"
    cp -f "$R/$d"/*.report "$R/$d"/*.bracken "$DEST/kraken/$d/" 2>/dev/null
    n=$(ls "$DEST/kraken/$d"/*.report 2>/dev/null | grep -vc bracken)
    printf 'kraken/%s/\t%s\t%s muestras\t-\t%s\n' "$d" "$R/$d" "$n" \
      "$([[ $d == kraken2 ]] && echo 'lecturas completas (artefacto Plasmodium)' || echo 'lecturas no-huesped, --confidence 0.1')" >> "$MAN"
    printf '  OK    kraken/%-45s %8s muestras\n' "$d/" "$n"
  else
    printf '  FALTA kraken/%s\n' "$d"
  fi
done

echo "== 7. Contabilidad de lecturas =="
copiar lecturas/read_accounting.tsv "$Q/read_accounting.tsv" "contabilidad de lecturas por etapa"
if [[ -d "$R/mapping_stats_host" ]]; then
  mkdir -p "$DEST/lecturas/flagstat_host"
  cp -f "$R/mapping_stats_host"/*_flagstat.txt "$DEST/lecturas/flagstat_host/" 2>/dev/null
  printf '  OK    lecturas/flagstat_host/                              %8s ficheros\n' \
    "$(ls "$DEST/lecturas/flagstat_host" 2>/dev/null | wc -l)"
fi

# ---- resumen -----------------------------------------------------------------
echo
echo "========================================================"
printf 'copiadas: %s   faltantes: %s\n' \
  "$(awk -F'\t' 'NR>1 && $3!="FALTA"' "$MAN" | wc -l)" \
  "$(awk -F'\t' 'NR>1 && $3=="FALTA"' "$MAN" | wc -l)"
echo "destino:  $DEST"
echo "tamanio:  $(du -sh "$DEST" | cut -f1)"
echo "manifiesto: $MAN"
echo "========================================================"
echo
echo "Lo que falte, revisar si ese paso llego a correr:"
awk -F'\t' 'NR>1 && $3=="FALTA"{printf "  %s\n", $2}' "$MAN"
