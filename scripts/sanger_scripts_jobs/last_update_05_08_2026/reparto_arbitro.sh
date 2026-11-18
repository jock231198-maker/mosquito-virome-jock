#!/bin/bash
# ---------------------------------------------------------------------------
# reparto_arbitro.sh — cruza las categorias de la anotacion viral con el
#                      veredicto del arbitro de nr.
#
# Responde las dos preguntas que quedan abiertas para Resultados:
#
#   1. Que familias virales NO producen ningun veredicto viral
#      -> son ruido sistematico de la anotacion por homologia
#   2. Que proporcion de cada categoria de la base viral propia
#      (conocido_nt, pariente_nt, divergente_solo_aa) resulta ser huesped
#      -> es la justificacion empirica del arbitraje
#
# Lee las columnas POR NOMBRE, asi que da igual el orden o el ancho del
# fichero. Sirve tanto para el fichero de contigs como para el de vOTUs.
#
# USO
#   ./reparto_arbitro.sh FICHERO_VEREDICTO.tsv [familia|categoria|ambos]
#
#   A=$RESULTS_DIR/contig_annot_min500/contigs_union_MERIDA_nr_jr46
#   ./reparto_arbitro.sh $A/veredicto_por_votu.tsv ambos
# ---------------------------------------------------------------------------
set -euo pipefail

F="${1:-}"
MODO="${2:-ambos}"
[[ -f "$F" ]] || { echo "uso: $0 FICHERO_VEREDICTO.tsv [familia|categoria|ambos]"; exit 1; }

reparto() {   # $1 = columna a agrupar, $2 = titulo
  awk -F'\t' -v COL="$1" -v TIT="$2" '
    NR == 1 {
      for (i = 1; i <= NF; i++) h[$i] = i
        # alias: TABLA_FINAL_votus_arb.tsv llama a estas columnas de otra
        # forma (categoria, familia, nr_veredicto). Asi el reparto corre
        # directo sobre ella, sin renombrar la cabecera antes.
        if (!("veredicto" in h) && ("nr_veredicto" in h)) h["veredicto"] = h["nr_veredicto"]
        if (!(COL in h)) sub(/_viral$/, "", COL)
      if (!(COL in h)) { printf("  (el fichero no tiene columna %s)\n", COL); exit 0 }
      if (!("veredicto" in h)) { print "  (el fichero no tiene columna veredicto)"; exit 0 }
      next
    }
    {
      g = $(h[COL]); v = $(h["veredicto"])
      if (g == "" || g == "NA") g = "(vacio)"
      tot[g]++
      if (v == "viral")                vir[g]++
      else if (v ~ /^celular_Artro/)   art[g]++
      else if (v ~ /^celular_/)        cel[g]++
      else if (v == "sin_hit_nr")      nohit[g]++
      else if (v == "ambiguo")         amb[g]++
      n++
    }
    END {
      if (!n) exit 0
      printf("\n== %s ==\n", TIT)
      printf("%-36s %8s %8s %8s %8s %8s %8s %7s\n",
             "", "total", "viral", "artrop.", "otro cel", "sin hit", "ambiguo", "% viral")
      for (k in tot)
        printf("%-36s %8d %8d %8d %8d %8d %8d %6.1f%%\n", substr(k,1,36), tot[k],
               vir[k]+0, art[k]+0, cel[k]+0, nohit[k]+0, amb[k]+0,
               100*(vir[k]+0)/tot[k])
    }' "$F" | { IFS= read -r a || exit 0; IFS= read -r b; IFS= read -r c
                printf '%s\n%s\n%s\n' "$a" "$b" "$c"; sort -k2 -rn; }

  echo
  echo "  --- las que NO producen ningun veredicto viral ---"
  local tmp; tmp="$(mktemp)"
  awk -F'\t' -v COL="$1" '
    NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i
      # alias: TABLA_FINAL_votus_arb.tsv llama a estas columnas de otra
      # forma (categoria, familia, nr_veredicto). Asi el reparto corre
      # directo sobre ella, sin renombrar la cabecera antes.
      if (!("veredicto" in h) && ("nr_veredicto" in h)) h["veredicto"] = h["nr_veredicto"]
      if (!(COL in h)) sub(/_viral$/, "", COL)
              if (!(COL in h) || !("veredicto" in h)) exit 0; next }
    { g = $(h[COL]); v = $(h["veredicto"])
      if (g == "" || g == "NA") g = "(vacio)"
      tot[g]++; if (v == "viral") vir[g]++ }
    END { n = 0
          for (k in tot) if (!(k in vir)) { printf("      %-36s %6d contigs, 0 virales\n", substr(k,1,36), tot[k]); n++ }
          if (!n) print "NOTA      (ninguna: todas aportan al menos un veredicto viral)"
          else printf("NOTA      -> %d categorias sin un solo veredicto viral\n", n) }' "$F" > "$tmp"
  grep -v '^NOTA' "$tmp" | sort -k2 -rn || true
  sed -n 's/^NOTA//p' "$tmp"
  rm -f "$tmp"
}

echo "fichero: $F"
echo "filas:   $(( $(wc -l < "$F") - 1 ))"

[[ "$MODO" == "familia"   || "$MODO" == "ambos" ]] && reparto "familia_viral"   "FAMILIA VIRAL x VEREDICTO DE nr"
[[ "$MODO" == "categoria" || "$MODO" == "ambos" ]] && reparto "categoria_viral" "CATEGORIA DE LA BASE VIRAL PROPIA x VEREDICTO DE nr"
exit 0
