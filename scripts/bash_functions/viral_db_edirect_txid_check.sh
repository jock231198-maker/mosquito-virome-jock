# En otra terminal — monitorear el progreso en tiempo real
watch -n 30 "grep -c '^>' viral_nt.fasta"

# Ver cuánto pesa el archivo
watch -n 30 "ls -lh viral_nt.fasta"

# Ver el lote actual
watch -n 5 "cat progress.txt"

# Si se interrumpe — simplemente volver a ejecutar
# El script retoma desde progress.txt automáticamente
./download_viral_nt_v7.sh