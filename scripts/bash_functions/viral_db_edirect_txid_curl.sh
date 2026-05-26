#!/bin/bash
# download_viral_nt_v6.sh — usa curl directo a NCBI eutils

EMAIL="tu@email.com"                          # ← cambia esto
DB="nuccore"
QUERY="Viruses[Organism]+AND+refseq[filter]"
OUTFILE="viral_nt.fasta"
BATCH=200
RETRIES=5
PROGRESS_FILE="progress.txt"
FAILED_FILE="failed_batches.txt"
BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# ─────────────────────────────────────────
# Paso 1: esearch vía curl → obtener WebEnv y total
# ─────────────────────────────────────────
echo "🔍 Consultando NCBI via curl..."

SEARCH=$(curl -s "${BASE}/esearch.fcgi" \
    --data-urlencode "db=${DB}" \
    --data-urlencode "term=Viruses[Organism] AND refseq[filter]" \
    --data-urlencode "usehistory=y" \
    --data-urlencode "retmax=0" \
    --data-urlencode "email=${EMAIL}")

echo "XML recibido:"
echo "$SEARCH"

TOTAL=$(echo  "$SEARCH" | grep -o '<Count>[^<]*'    | sed 's/<Count>//'    | tr -d '[:space:]')
WEBENV=$(echo "$SEARCH" | grep -o '<WebEnv>[^<]*'   | sed 's/<WebEnv>//'   | tr -d '[:space:]')
QKEY=$(echo   "$SEARCH" | grep -o '<QueryKey>[^<]*' | sed 's/<QueryKey>//' | tr -d '[:space:]')

echo ""
echo "Total : $TOTAL"
echo "WebEnv: $WEBENV"
echo "QKey  : $QKEY"

if [[ -z "$TOTAL" || -z "$WEBENV" || -z "$QKEY" ]]; then
    echo "❌ No se pudo parsear la respuesta."
    exit 1
fi

# ─────────────────────────────────────────
# Reanudar desde progress.txt si existe
# ─────────────────────────────────────────
START_FROM=0
if [[ -f "$PROGRESS_FILE" ]]; then
    START_FROM=$(cat "$PROGRESS_FILE")
    echo "▶️  Reanudando desde lote: $START_FROM"
fi

# ─────────────────────────────────────────
# Paso 2: efetch vía curl con retstart/retmax
# ─────────────────────────────────────────
for ((START=START_FROM; START<TOTAL; START+=BATCH)); do
    END=$((START + BATCH))
    echo "📥 Lote: $START - $END de $TOTAL"

    SUCCESS=false
    SLEEP=3

    for ((TRY=1; TRY<=RETRIES; TRY++)); do

        # ✅ curl directo — mismos parámetros que usaba nquire internamente
        RESULT=$(curl -s "${BASE}/efetch.fcgi" \
            --data-urlencode "db=${DB}"         \
            --data-urlencode "WebEnv=${WEBENV}" \
            --data-urlencode "query_key=${QKEY}"\
            --data-urlencode "retstart=${START}" \
            --data-urlencode "retmax=${BATCH}"  \
            --data-urlencode "rettype=fasta"    \
            --data-urlencode "retmode=text"     \
            --data-urlencode "email=${EMAIL}")

        if echo "$RESULT" | grep -q "^>"; then
            echo "$RESULT" >> "$OUTFILE"
            echo "$START"  >  "$PROGRESS_FILE"
            COUNT=$(echo "$RESULT" | grep -c "^>")
            echo "  ✅ $COUNT secuencias guardadas"
            SUCCESS=true
            break
        else
            echo "  ❌ Intento $TRY — Respuesta:"
            echo "$RESULT" | head -3
            echo "  ⏳ Esperando ${SLEEP}s..."
            sleep "$SLEEP"
            SLEEP=$((SLEEP * 2))
        fi
    done

    [[ "$SUCCESS" == false ]] && \
        echo "$START" >> "$FAILED_FILE" && \
        echo "  ⚠️  Lote $START en $FAILED_FILE"

    sleep 1   # respetar límite NCBI: máx 3 req/seg
    SLEEP=3
done

echo "✅ Descarga completa: $OUTFILE"
echo "Total secuencias: $(grep -c '^>' $OUTFILE)"