#!/bin/bash
# Paso 1 — Descomprimir los archivos gbvrl
cd /Users/JK/Desktop/ncbi_database/viral_db

gunzip -k gbvrl*.seq.gz

# Paso 2 — Convertir formato .seq a FASTA
# Los archivos .seq están en formato GenBank, hay que convertirlos
for FILE in *.seq; do
    echo "Convirtiendo $FILE..."
    python3 -c "
from Bio import SeqIO
import sys
with open('${FILE}'.replace('.seq','.fasta'), 'w') as out:
    for rec in SeqIO.parse('$FILE', 'genbank'):
        out.write(f'>{rec.id} {rec.description}\n{str(rec.seq)}\n')
"
done

# Paso 3 — Combinar todo en un solo FASTA
cat /Users/JK/Desktop/ncbi_database/viral_db/*.fasta \
    /Users/JK/Desktop/ncbi_database/viral_db_ftp_data/viral_nt.fasta \
    > /Users/JK/Desktop/ncbi_database/viral_combined.fasta

# Verificar total de secuencias combinadas
grep -c "^>" /Users/JK/Desktop/ncbi_database/viral_combined.fasta