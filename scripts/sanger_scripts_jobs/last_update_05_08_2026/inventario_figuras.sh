#!/bin/bash
# inventario_figuras.sh — SOLO LECTURA. No escribe ni lanza nada.
#
# 1. Verifica que existan las entradas de las figuras en R (con conteos esperados)
# 2. Diagnostica los puntos abiertos: Panayarvi, parche de mapping_stats,
#    arbitro de nr, nudivirus de M71_S19, metadatos
#
# Uso:
#   cd "$SCRIPTS_DIR" && bash inventario_figuras.sh 2>&1 | tee "inventario_$(date +%F).txt"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh" >/dev/null 2>&1
export LC_ALL=C

Q="$RESULTS_DIR/qc_control"
G="$RESULTS_DIR/genomad_all/all_final/all_final_summary"

sec()   { printf '\n== %s ==\n' "$1"; }
corto() { local p="${1/#$RESULTS_DIR/R:}"; echo "${p/#$SCRATCH/S:}"; }

# chk <fichero> [lineas_esperadas] [nota]
chk() {
  local f="$1" exp="${2:-}" nota="${3:-}" p; p="$(corto "$f")"
  [[ -e "$f" ]] || { printf '  FALTA    %-55s %s\n' "$p" "$nota"; return; }
  [[ -s "$f" ]] || { printf '  VACIO    %-55s %s\n' "$p" "$nota"; return; }
  local n sz flag="OK"
  n=$(wc -l < "$f"); sz=$(du -h "$f" | cut -f1)
  [[ -n "$exp" && "$n" != "$exp" ]] && flag="DIF:$exp"
  printf '  %-8s %-55s %8s lin %6s  %s\n' "$flag" "$p" "$n" "$sz" "$nota"
}

# chkfa <fasta> [secuencias_esperadas] [nota]
chkfa() {
  local f="$1" exp="${2:-}" nota="${3:-}" p; p="$(corto "$f")"
  [[ -s "$f" ]] || { printf '  FALTA    %-55s %s\n' "$p" "$nota"; return; }
  local n flag="OK"; n=$(grep -c '^>' "$f")
  [[ -n "$exp" && "$n" != "$exp" ]] && flag="DIF:$exp"
  printf '  %-8s %-55s %8s seq %6s  %s\n' "$flag" "$p" "$n" "$(du -h "$f" | cut -f1)" "$nota"
}

echo "inventario  $(date '+%F %H:%M')  $(hostname)"
echo "R=$RESULTS_DIR"
echo "S=$SCRATCH"

# ---------------------------------------------------------------- muestras
sec "Muestras"
chk "$SCRATCH/asm_samples.txt" 22 "lista maestra"
[[ -s "$SCRATCH/asm_samples.txt" ]] && paste -sd' ' "$SCRATCH/asm_samples.txt" | fold -w 100 | sed 's/^/  /'

# ---------------------------------------------------------------- ensamblaje
sec "Ensamblaje y union"
printf '  directorios en S:/union: %s (esperado 22)\n' "$(ls -d "$SCRATCH"/union/*/ 2>/dev/null | wc -l)"
chk "$Q/contigs_union.tsv"     6005 "muestra/ensamblador/contig/longitud/cov"
chk "$Q/contigs_union.tsv.cov" 6005 "con cov_fix"
if [[ -s "$Q/contigs_union.tsv" ]]; then
  echo "  contigs por ensamblador:"
  awk -F'\t' 'NR>1{n[$2]++} END{for(a in n) printf "    %-12s %5d\n", a, n[a]}' "$Q/contigs_union.tsv" | sort -k2 -rn
fi

# ---------------------------------------------------------------- lecturas
sec "Contabilidad de lecturas (normalizacion)"
for d in "$RESULTS_DIR"/mapping_stats*; do
  [[ -d "$d" ]] && printf '  %-45s %4s ficheros\n' "$(corto "$d")" "$(ls "$d" | wc -l)"
done
find "$RESULTS_DIR" -maxdepth 3 -iname '*read_account*' 2>/dev/null | head -5 | while read -r f; do chk "$f"; done

# ---------------------------------------------------------------- geNomad / CheckV
sec "geNomad y CheckV"
chkfa "$G/all_final_virus.fna"         1144 "contigs virales"
chk   "$G/all_final_virus_summary.tsv" 1145
find "$RESULTS_DIR" -maxdepth 4 -name 'quality_summary.tsv' 2>/dev/null | while read -r f; do chk "$f"; done

# ---------------------------------------------------------------- vOTUs y mapeo
sec "vOTUs y mapeo"
chkfa "$Q/votus"            748 "representativas CD-HIT"
chk   "$Q/votus.clstr"      ""  "clusters"
chk   "$Q/votus_master.tsv" ""  ""
printf '  BAM en S:/mapped_votus/votus: %s (esperado 22)\n' "$(ls "$SCRATCH"/mapped_votus/votus/*.bam 2>/dev/null | wc -l)"
printf '  .cov en R:/qc_control/cov:    %s (esperado 22)\n' "$(ls "$Q"/cov/*.cov 2>/dev/null | wc -l)"
chk   "$Q/votu_coverage.tsv" 16456 "22 x 748"

# ---------------------------------------------------------------- anotacion
sec "Anotacion (DIAMOND, blastn, Kraken2/Bracken)"
find "$RESULTS_DIR" -maxdepth 3 \( -iname 'diamond*.tsv' -o -iname 'blastn*.tsv' -o -iname '*bracken*' \) \
     2>/dev/null | sort | while read -r f; do chk "$f"; done

# ---------------------------------------------------------------- tablas finales
sec "Tablas finales (entrada directa para R)"
chk "$Q/TABLA_FINAL_votus.tsv"        749 "cabecera + 748"
chk "$Q/prevalencia_por_especie.tsv"
chk "$Q/diversidad_por_muestra.tsv"   23  "cabecera + 22"
if [[ -s "$Q/TABLA_FINAL_votus.tsv" ]]; then
  echo "  columnas:"
  head -1 "$Q/TABLA_FINAL_votus.tsv" | tr '\t' '\n' | awk '{printf "    %2d %s\n", NR, $0}'
  echo "  reparto por clase:"
  awk -F'\t' 'NR>1{n[$5]++} END{for(c in n) printf "    %-16s %4d\n", c, n[c]}' "$Q/TABLA_FINAL_votus.tsv" | sort -k2 -rn
fi

# ---------------------------------------------------------------- Panayarvi
sec "PUNTO 1 — Panayarvi"
F="$Q/TABLA_FINAL_votus.tsv"
if [[ -s "$F" ]]; then
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[i]=$i; next}
              $7 ~ /Panayarvi/ {for(i=1;i<=NF;i++) if(i!=21) printf "  %-18s %s\n", h[i], $i; print "  --"}' "$F"
  V=$(awk -F'\t' 'NR>1 && $7 ~ /Panayarvi/ {print $1; exit}' "$F")
  if [[ -n "$V" ]]; then
    echo "  cobertura por muestra (vOTU $V), ordenada por profundidad:"
    printf '    %-10s %8s %12s\n' muestra "%cub" profundidad
    awk -F'\t' -v V="$V" '$2==V{printf "    %-10s %8.1f %12.1f\n", $1, $5, $6}' "$Q/votu_coverage.tsv" | sort -k3,3gr
    echo "  vOTU representativa: ¿es el partitivirus divergente M35_S31 (78.3% id)?"
    [[ "$V" == M35_S31__trinity__TRINITY_DN53_c0_g1_i1* ]] && echo "    si" || echo "    no — es otra vOTU, revisar"
    echo "  referencia — profundidad media de las otras Partitiviridae identificadas:"
    awk -F'\t' 'NR>1 && $5=="identificado" && $8 ~ /Partitiviridae/ && $7 !~ /Panayarvi/ \
                {printf "    %-30s n_map=%-3s prof=%s\n", substr($7,1,30), $18, $19}' "$F"
  fi
fi

# ---------------------------------------------------------------- parche
sec "PUNTO 2 — parche mapping_stats en bowtie_map_worker"
grep -n 'statsdir' "$SCRIPT_DIR"/bowtie_map_worker*.sh 2>/dev/null || echo "  no aparece 'statsdir' en ningun worker"

# ---------------------------------------------------------------- arbitro nr / jr46
sec "PUNTO 3 — arbitro de nr y corrida jr46 (hecha 7-9 oct)"
chk "$RESULTS_DIR/votus_arbitro.tsv"     "" "veredicto nr por vOTU"
chk "$RESULTS_DIR/votus_nr_arbiter.tsv"  "" ""
for d in "$RESULTS_DIR/votu_annot/votus_nr_jr46" \
         "$RESULTS_DIR/contig_annot/contigs_union_MERIDA_nr_jr46" \
         "$RESULTS_DIR/contig_annot/contigs_union_MERIDA_nr_jr46/cruce_votus" \
         "$RESULTS_DIR/votu_annot/votus_viraldb20260925"; do
  if [[ -d "$d" ]]; then printf '  OK       %-55s %4s ficheros\n' "$(corto "$d")" "$(ls "$d" | wc -l)"
  else printf '  FALTA    %s\n' "$(corto "$d")"; fi
done
echo "  ¿TABLA_FINAL ya lleva el veredicto del arbitro? (busca columna 'arb' o 'nr')"
cols=$(head -1 "$Q/TABLA_FINAL_votus.tsv" 2>/dev/null | tr '\t' '\n' | grep -in -E 'arb|nr')
if [[ -n "$cols" ]]; then echo "$cols" | sed 's/^/    /'
else echo "    no: el arbitro aun no esta integrado en la tabla final"; fi

# ---------------------------------------------------------------- nudivirus
sec "PUNTO 5 — nudivirus de M71_S19 (fuera del catalogo)"
N="$RESULTS_DIR/contig_annot/contigs_union_MERIDA_nr_jr46"
f=$(ls "$N"/*veredicto*.tsv "$N"/*arbit*.tsv 2>/dev/null | head -1)
if [[ -n "$f" ]]; then
  echo "  fichero: $(corto "$f")"
  grep -i 'nudi' "$f" | awk -F'\t' '{split($1,a,"__"); n[a[1]]++} END{for(s in n) printf "    %-10s %3d contigs con hit nudi\n", s, n[s]}' | sort -k2 -rn
else
  echo "  no encuentro la tabla de veredictos de contigs en $(corto "$N")"; ls "$N" 2>/dev/null | sed 's/^/    /'
fi

# ---------------------------------------------------------------- kraken
sec "Kraken2 / Bracken (sin huesped)"
for d in "$RESULTS_DIR/kraken2_nohost" "$RESULTS_DIR/kraken2"; do
  [[ -d "$d" ]] && printf '  %-45s %4s ficheros\n' "$(corto "$d")" "$(ls "$d" | wc -l)"
done

# ---------------------------------------------------------------- metadatos
sec "PUNTO 4 — metadatos existentes"
for base in "$SCRIPT_DIR/../.." "$RESULTS_DIR" "$SCRATCH"; do
  find "$base" -maxdepth 3 \( -iname '*metadat*' -o -iname '*sample_info*' -o -iname '*muestras*.tsv' -o -iname '*muestras*.csv' \) \
       2>/dev/null
done | sort -u | while read -r f; do chk "$f"; done

# ---------------------------------------------------------------- estado y entorno
sec "Repositorio y entorno"
git -C "$SCRIPT_DIR" status --short 2>/dev/null | head -15 | sed 's/^/  /'
git -C "$SCRIPT_DIR" log --oneline -3 2>/dev/null | sed 's/^/  /'
for t in Rscript quarto; do printf '  %-8s %s\n' "$t" "$(command -v $t || echo 'no en PATH')"; done

sec "Tamano a traer con rsync (tablas pequenas)"
du -ch "$Q"/TABLA_FINAL_votus.tsv "$Q"/prevalencia_por_especie.tsv "$Q"/diversidad_por_muestra.tsv \
       "$Q"/votu_coverage.tsv "$Q"/contigs_union.tsv "$Q"/votus.clstr 2>/dev/null | tail -1 | sed 's/^/  /'

echo; echo "fin"
