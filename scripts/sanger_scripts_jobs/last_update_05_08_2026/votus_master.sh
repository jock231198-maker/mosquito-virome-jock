#!/bin/bash
# votus_master.sh — tabla maestra de vOTUs (antes se hacia a mano)
# Une: votu_resumen.tsv (votu_catalog.sh) + virus_summary de geNomad + quality_summary de CheckV
# Uso: votus_master.sh <votu_resumen.tsv> <genomad_virus_summary.tsv> <checkv_quality_summary.tsv> <salida.tsv>
# Columnas: votu longitud score hallmarks n_contigs n_muestras muestras checkv_quality completeness taxonomia
set -euo pipefail
res=${1:?votu_resumen.tsv}; gen=${2:?virus_summary.tsv}; chk=${3:?quality_summary.tsv}; out=${4:?salida}
awk -F'\t' -v OFS='\t' '
  FILENAME==ARGV[1] { if (FNR==1) {for(i=1;i<=NF;i++) g[$i]=i; next}
                      id=$g["seq_name"]; sc[id]=$g["virus_score"]; hm[id]=$g["n_hallmarks"]; tx[id]=$g["taxonomy"]; next }
  FILENAME==ARGV[2] { if (FNR==1) {for(i=1;i<=NF;i++) k[$i]=i; next}
                      id=$k["contig_id"]; q[id]=$k["checkv_quality"]; c=$k["completeness"]; cp[id]=(c==""?"NA":c); next }
  FNR==1 { for(i=1;i<=NF;i++) r[$i]=i
           print "votu","longitud","score","hallmarks","n_contigs","n_muestras","muestras","checkv_quality","completeness","taxonomia"; next }
  { v=$r["votu"]; n=split($r["muestras"], a, ","); asort(a); m=a[1]; for(i=2;i<=n;i++) m=m "," a[i]
    print v, $r["longitud"], sc[v], hm[v], $r["n_contigs"], $r["n_muestras_ensamblado"], m, q[v], ((v in cp)?cp[v]:"NA"), tx[v] }
' "$gen" "$chk" "$res" > "$out"
echo "$(( $(wc -l < "$out") - 1 )) vOTUs -> $out"
