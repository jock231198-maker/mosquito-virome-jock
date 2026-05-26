#!/bin/bash
# download_viral_nt_v7.sh — sintaxis correcta para EDirect 22.4

DB="nuccore"
QUERY='Viruses[Organism] AND refseq[filter]'
OUTFILE="/Users/JK/Desktop/ncbi_database/viral_db_ftp_data/viral_nt.fasta"
BATCH=200
RETRIES=5
PROGRESS_FILE="/Users/JK/Desktop/ncbi_database/viral_db_ftp_data/progress.txt"
FAILED_FILE="/Users/JK/Desktop/ncbi_database/viral_db_ftp_data/failed_batches.txt"

# ─────────────────────────────────────────
# Paso 1: Obtener total
# ─────────────────────────────────────────
echo "🔍 Consultando NCBI..."

TOTAL=$(esearch -db "$DB" -query "$QUERY" | \
        grep -o '<Count>[^<]*</Count>'    | \
        sed 's/<[^>]*>//g'               | \
        tr -d '[:space:]')

echo "Total de secuencias: $TOTAL"
[[ -z "$TOTAL" ]] && echo "❌ Sin respuesta de NCBI" && exit 1

# ─────────────────────────────────────────
# Reanudar desde progress.txt si existe
# ─────────────────────────────────────────
START_FROM=1      # EDirect usa base 1, no base 0
if [[ -f "$PROGRESS_FILE" ]]; then
    START_FROM=$(cat "$PROGRESS_FILE")
    echo "▶️  Reanudando desde: $START_FROM"
fi

# ─────────────────────────────────────────
# Paso 2: Descarga con -start y -stop (sintaxis 22.4)
# ─────────────────────────────────────────
for ((START=START_FROM; START<=TOTAL; START+=BATCH)); do
    STOP=$((START + BATCH - 1))
    # No pasar del total
    [[ $STOP -gt $TOTAL ]] && STOP=$TOTAL

    echo "📥 Lote: $START - $STOP de $TOTAL"

    SUCCESS=false
    SLEEP=3

    for ((TRY=1; TRY<=RETRIES; TRY++)); do

        # ✅ Sintaxis exacta del help de EDirect 22.4
        RESULT=$(esearch -db "$DB" -query "$QUERY" | \
                 efetch -format fasta \
                        -start "$START" \
                        -stop  "$STOP"  \
                 2>&1)

        if echo "$RESULT" | grep -q "^>"; then
            echo "$RESULT" >> "$OUTFILE"
            # Guardar el inicio del SIGUIENTE lote
            echo "$((STOP + 1))" > "$PROGRESS_FILE"
            COUNT=$(echo "$RESULT" | grep -c "^>")
            echo "  ✅ $COUNT secuencias guardadas"
            SUCCESS=true
            break
        else
            echo "  ❌ Intento $TRY — Error:"
            echo "$RESULT" | head -3
            echo "  ⏳ Esperando ${SLEEP}s..."
            sleep "$SLEEP"
            SLEEP=$((SLEEP * 2))
        fi
    done

    [[ "$SUCCESS" == false ]] && \
        echo "$START" >> "$FAILED_FILE" && \
        echo "  ⚠️  Lote $START en $FAILED_FILE"

    sleep 1
    SLEEP=3
done

echo ""
echo "✅ Descarga completa: $OUTFILE"
echo "📊 Total secuencias: $(grep -c '^>' $OUTFILE)"