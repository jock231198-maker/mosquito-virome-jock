#!/usr/bin/env bash
# =============================================================================
# anfevirus_especie.sh -- ¿El anfevirus de Merida es una especie nueva?
#
#   El criterio del ICTV para Xinmoviridae es por PROTEINA, no por nucleotido:
#     especies distintas dentro de un genero -> RdRp con <=66 % de identidad aa
#     generos distintos                      -> RdRp con <=60 % de identidad aa
#
#   Este script saca la L (RdRp) del contig de Merida, baja las L de todos los
#   Xinmoviridae de RefSeq, las alinea y calcula la matriz de identidad por
#   pares. El alineamiento sirve luego para el arbol en MEGA.
#
# Subcomandos:
#   check     entornos, proteinas de geNomad, contig objetivo
#   extraer   saca la proteina mas larga del contig de Merida (= la L)
#   refs      baja las L de Xinmoviridae desde NCBI (head node: usa red)
#   alinear   mafft + matriz de identidad por pares + veredicto
#
# Variables:
#   CONTIG_ID   nombre del contig del anfevirus en los FASTA del proyecto
#   MIN_L_AA    longitud minima para considerar una proteina como L (def. 1500)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/config.sh"

OUT_DIR="${RESULTS_DIR}/anfevirus_especie"
mkdir -p "${OUT_DIR}"

CONTIG_ID="${CONTIG_ID:-M1_S1__rnaviral__NODE_1_length_12945_cov_150.857553}"
MIN_L_AA="${MIN_L_AA:-1500}"
ENV_MAFFT="${ENV_MAFFT:-mafft_7.526}"
THREADS="${THREADS:-4}"

L_MERIDA="${OUT_DIR}/L_merida.faa"
L_REFS="${OUT_DIR}/L_xinmoviridae.faa"
L_TODAS="${OUT_DIR}/L_merida_y_referencias.faa"
L_ALN="${OUT_DIR}/L_align.faa"
MATRIZ="${OUT_DIR}/identidad_pares.tsv"

RED=$'\033[31m'; GRN=$'\033[32m'; YLW=$'\033[33m'; RST=$'\033[0m'
log()  { printf '%s[%s]%s %s\n' "${GRN}" "$(date +%H:%M:%S)" "${RST}" "$*"; }
warn() { printf '%s[WARN]%s %s\n' "${YLW}" "${RST}" "$*" >&2; }
die()  { printf '%s[ERROR]%s %s\n' "${RED}" "${RST}" "$*" >&2; exit 1; }

# Localiza el .faa de proteinas que dejo geNomad
proteinas_genomad() {
  local p
  for p in "${RESULTS_DIR}"/genomad_min500/all_final/*_annotate/*_proteins.faa \
           "${RESULTS_DIR}"/genomad_min500/*/*_annotate/*_proteins.faa \
           "${RESULTS_DIR}"/genomad/*/*_annotate/*_proteins.faa; do
    [[ -s "$p" ]] && { echo "$p"; return 0; }
  done
  return 1
}

# ---------------------------------------------------------------------------
cmd_check() {
  local fail=0

  log "Entorno ${ENV_MAFFT}:"
  if ( activate_env "${ENV_MAFFT}" >/dev/null 2>&1; command -v mafft >/dev/null ); then
    ( activate_env "${ENV_MAFFT}" >/dev/null 2>&1
      printf '  mafft      %s\n' "$(mafft --version 2>&1 | head -1)" )
  else
    printf '  %smafft no disponible%s\n' "${RED}" "${RST}"; fail=1
  fi

  log "Proteínas de geNomad:"
  local faa
  if faa="$(proteinas_genomad)"; then
    printf '  OK    %s  (%s proteínas)\n' "$faa" "$(grep -c '^>' "$faa")"
    local n; n="$(grep -c "^>${CONTIG_ID}" "$faa" || true)"
    if [[ "$n" -gt 0 ]]; then
      printf '  %s proteínas en el contig objetivo:\n' "$n"
      awk -v C="${CONTIG_ID}" '
        /^>/ { p = (index($0, ">" C) == 1); if (p) { h = substr($1,2) } ; l = 0; next }
        p { l += length($0); len[h] = l }
        END { for (k in len) printf "        %-70s %5d aa\n", k, len[k] }' "$faa" | sort -k2 -nr
    else
      printf '  %sel contig %s no aparece en ese .faa%s\n' "${RED}" "${CONTIG_ID}" "${RST}"
      printf '        prueba: grep -c "NODE_1_length_12945" %s\n' "$faa"
      fail=1
    fi
  else
    printf '  %sno encontré el *_proteins.faa de geNomad%s\n' "${RED}" "${RST}"
    printf '        búscalo con: find "%s" -name "*_proteins.faa" | head\n' "${RESULTS_DIR}"
    fail=1
  fi

  log "Red (NCBI):"
  curl -sf --max-time 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi" >/dev/null \
    && printf '  OK    alcanzable\n' \
    || printf '  %ssin salida a NCBI desde este nodo%s\n' "${YLW}" "${RST}"

  [[ $fail -eq 0 ]] && log "check OK" || die "check con fallos"
}

# ---------------------------------------------------------------------------
# extraer -- la L de Mononegavirales es con mucho la proteina mas larga (~2100 aa)
cmd_extraer() {
  local faa; faa="$(proteinas_genomad)" || die "no encontré las proteínas de geNomad; corre '$0 check'"

  awk -v C="${CONTIG_ID}" '
    /^>/ { p = (index($0, ">" C) == 1); if (p) h = substr($1,2); next }
    p { s[h] = s[h] $0 }
    END {
      best = ""; bl = 0
      for (k in s) if (length(s[k]) > bl) { bl = length(s[k]); best = k }
      if (best == "") exit 1
      printf ">L_MERIDA %s %d_aa\n%s\n", best, bl, s[best]
    }' "$faa" > "${L_MERIDA}" || die "no hallé proteínas del contig ${CONTIG_ID}"

  local la; la="$(awk 'NR>1{n+=length($0)} END{print n}' "${L_MERIDA}")"
  log "L de Mérida -> ${L_MERIDA}  (${la} aa)"
  head -1 "${L_MERIDA}" | sed 's/^/  /'

  if [[ "${la}" -lt "${MIN_L_AA}" ]]; then
    warn "solo ${la} aa. La L de Mononegavirales ronda las 2 100 aa."
    warn "Puede que prodigal partiera el ORF. Revisa las proteínas del contig con '$0 check'."
  fi
}

# ---------------------------------------------------------------------------
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

cmd_refs() {
  log "Buscando genomas RefSeq de Xinmoviridae ..."
  local ids
  ids=$(curl -sf --max-time 60 -G "${EUTILS}/esearch.fcgi" \
          --data-urlencode "db=nuccore" \
          --data-urlencode 'term="Xinmoviridae"[Organism] AND srcdb_refseq[PROP]' \
          --data-urlencode "retmax=100" \
        | grep -o '<Id>[0-9]*</Id>' | sed 's/<[^>]*>//g' | paste -sd,)
  [[ -n "${ids}" ]] || die "la consulta no devolvió nada; ¿hay red en este nodo?"
  log "  $(tr ',' '\n' <<< "${ids}" | wc -l) genomas"

  sleep 1
  log "Bajando las traducciones de sus CDS ..."
  curl -sf --max-time 180 -G "${EUTILS}/efetch.fcgi" \
       --data-urlencode "db=nuccore" \
       --data-urlencode "id=${ids}" \
       --data-urlencode "rettype=fasta_cds_aa" \
       --data-urlencode "retmode=text" > "${OUT_DIR}/.cds_aa.faa" \
    || die "falló la descarga de CDS"

  # Nos quedamos solo con las proteinas grandes: esas son las L (RdRp)
  awk -v MIN="${MIN_L_AA}" '
    /^>/ { if (h != "" && length(s) >= MIN) print h "\n" s; h = $0; s = ""; next }
    { s = s $0 }
    END { if (h != "" && length(s) >= MIN) print h "\n" s }' "${OUT_DIR}/.cds_aa.faa" \
    | awk '
        /^>/ { # nombre corto: organismo + accesion
               org = ""; acc = ""
               if (match($0, /\[protein=[^]]*\]/))  prot = substr($0, RSTART+9, RLENGTH-10)
               if (match($0, /lcl\|[A-Z0-9_.]+/))   acc  = substr($0, RSTART+4, RLENGTH-4)
               gsub(/[^A-Za-z0-9._-]/, "_", acc)
               print ">" acc; next }
        { print }' > "${L_REFS}"

  rm -f "${OUT_DIR}/.cds_aa.faa"
  log "L de referencia -> ${L_REFS}  ($(grep -c '^>' "${L_REFS}") secuencias ≥${MIN_L_AA} aa)"
  grep '^>' "${L_REFS}" | sed 's/^/  /' | head -30
  warn "Revisa la lista: deben ser L (RdRp), no glicoproteínas largas."
  warn "Comprueba que esté la de AeAV (NC_075274); es la comparación que importa."
}

# ---------------------------------------------------------------------------
cmd_alinear() {
  [[ -s "${L_MERIDA}" ]] || die "falta ${L_MERIDA}; corre '$0 extraer'"
  [[ -s "${L_REFS}"   ]] || die "falta ${L_REFS}; corre '$0 refs'"

  cat "${L_MERIDA}" "${L_REFS}" > "${L_TODAS}"
  log "$(grep -c '^>' "${L_TODAS}") secuencias a alinear"

  activate_env "${ENV_MAFFT}" >/dev/null
  mafft --auto --thread "${THREADS}" --quiet "${L_TODAS}" > "${L_ALN}" \
    || die "mafft falló"
  log "alineamiento -> ${L_ALN}"

  # Identidad por pares: posiciones identicas sobre posiciones comparables
  # (se omiten las columnas con hueco en cualquiera de las dos secuencias)
  awk '
    /^>/ { name = substr($1, 2); orden[++n] = name; next }
    { seq[name] = seq[name] toupper($0) }
    END {
      printf "seq_a\tseq_b\tidenticas\tcomparadas\tidentidad_pct\n"
      for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) {
        a = seq[orden[i]]; b = seq[orden[j]]
        L = length(a); id = 0; cmp = 0
        for (k = 1; k <= L; k++) {
          ca = substr(a, k, 1); cb = substr(b, k, 1)
          if (ca == "-" || cb == "-") continue
          cmp++
          if (ca == cb) id++
        }
        if (cmp > 0)
          printf "%s\t%s\t%d\t%d\t%.2f\n", orden[i], orden[j], id, cmp, 100 * id / cmp
      }
    }' "${L_ALN}" > "${MATRIZ}"
  log "matriz -> ${MATRIZ}"

  echo
  log "L de Mérida frente a cada referencia (de mayor a menor identidad):"
  awk -F'\t' 'NR>1 && ($1=="L_MERIDA" || $2=="L_MERIDA") {
        otra = ($1=="L_MERIDA") ? $2 : $1
        printf "  %-28s %7.2f %%\n", otra, $5
      }' "${MATRIZ}" | sort -k2 -nr

  echo
  local mejor
  mejor="$(awk -F'\t' 'NR>1 && ($1=="L_MERIDA" || $2=="L_MERIDA") {if($5>m) m=$5} END{printf "%.2f", m}' "${MATRIZ}")"
  printf '  Criterio ICTV (Xinmoviridae): especies distintas si la RdRp está en ≤66 %% de identidad aa.\n'
  printf '  Identidad máxima observada: %s %%\n' "${mejor}"
  if awk -v m="${mejor}" 'BEGIN{exit !(m > 66)}'; then
    printf '  %s=> por encima del umbral: misma especie que su pariente más cercano.%s\n' "${GRN}" "${RST}"
  else
    printf '  %s=> por debajo del umbral: candidata a especie nueva. Verifícalo con un árbol.%s\n' "${YLW}" "${RST}"
  fi
  echo
  printf '  El alineamiento %s se puede abrir en MEGA para el árbol.\n' "${L_ALN}"
}

# ---------------------------------------------------------------------------
case "${1:-}" in
  check)   cmd_check ;;
  extraer) cmd_extraer ;;
  refs)    cmd_refs ;;
  alinear) cmd_alinear ;;
  *) sed -n '2,24p' "${BASH_SOURCE[0]}"; exit 1 ;;
esac
