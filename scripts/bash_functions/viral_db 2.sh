#!/bin/bash
# download_viral_nt.sh — Descarga robusta con reintentos

EMAIL="tu@email.com"      # ← cambia esto (requerido por NCBI)
DB="nuccore"
QUERY='Viruses[Organism] AND refseq[filter]'   # RefSeq reduce volumen
OUTFILE="viral_nt.fasta"
BATCH=200          # lotes pequeños = menos fallos
RETRIES=5
SLEEP=3            # segundos entre lotes

# --- Paso 1: Obtener total de resultados ---
TOTAL=$(esearch -db "$DB" -query "$QUERY" | \
        grep "<Count>" | grep -o '[0-9]*')
echo "Total de secuencias: $TOTAL"

# --- Paso 2: Descarga por lotes con retry ---
for ((START=0; START<TOTAL; START+=BATCH)); do
    echo "Descargando lote: $START - $((START+BATCH))..."
    
    SUCCESS=false
    for ((TRY=1; TRY<=RETRIES; TRY++)); do
        RESULT=$(efetch \
            -db "$DB" \
            -query "$QUERY" \
            -format fasta \
            -start "$START" \
            -stop "$((START+BATCH-1))" \
            2>/dev/null)
        
        if [[ -n "$RESULT" ]]; then
            echo "$RESULT" >> "$OUTFILE"
            SUCCESS=true
            break
        else
            echo "  Intento $TRY fallido, esperando ${SLEEP}s..."
            sleep $SLEEP
            SLEEP=$((SLEEP * 2))   # backoff exponencial
        fi
    done
    
    if [[ "$SUCCESS" == false ]]; then
        echo "$START" >> failed_batches.txt
        echo "  ⚠️  Lote $START guardado en failed_batches.txt"
    fi
    
    sleep 1   # pausa entre lotes (política NCBI)
    SLEEP=3   # reset backoff
done

echo "Descarga completa: $OUTFILE"