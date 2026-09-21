# Chuleta de revisión — comandos rápidos por herramienta

Comandos para comprobar que un paso terminó bien y echar un vistazo a sus
resultados sin abrir ficheros enteros. Separados por herramienta, con cada
línea explicada.

Antes de cualquier bloque, **siempre**:

```bash
cd "$SCRIPTS_DIR" && source ./config.sh     # carga rutas: $SCRATCH, $RESULTS_DIR, $LOGS_DIR, $REFS_DIR...
export LC_ALL=C                              # ordenacion byte a byte; sin esto `join` falla en silencio
```

---

## 1. LSF — la cola y los trabajos

### ¿Qué hay corriendo?

```bash
bjobs -w                    # trabajos pendientes y en ejecucion; -w no recorta los nombres
bjobs -w | grep vmap        # solo los de un nombre
bjobs -a                    # incluye los terminados recientemente (DONE / EXIT)
bjobs -p                    # solo pendientes, con el MOTIVO de la espera
bjobs -l 234311 | head -30  # detalle completo de uno: recursos pedidos, nodo, motivo
```

`bjobs -p` es el que explica por qué algo no arranca: "Job array has reached
its running element limit" (el `%N`), "Affinity resource requirement cannot be
met" (pides demasiados núcleos en un solo nodo), o falta de memoria libre.

Columnas de `bjobs`: `JOBID USER STAT QUEUE FROM_HOST EXEC_HOST JOB_NAME
SUBMIT_TIME`. **`STAT` va antes que `JOB_NAME`**: un `grep "trinity.*RUN"` nunca
casa. Usa `grep trinity | awk '{print $3}'`.

### ¿Terminó bien?

```bash
grep -hE "Successfully|Exited|TERM_|Max Memory|Run time" "$LOGS_DIR"/gnall.*.log
```

| salida | significado |
|---|---|
| `Successfully completed.` | salió con código 0 |
| `Exited with exit code N` | falló; mira el `.err` |
| `TERM_MEMLIMIT` | LSF lo mató por pasar el `-M` |
| `TERM_RUNLIMIT` | se pasó del tiempo de la cola (`normal` = 720 min) |
| `Max Memory` | pico real — úsalo para el `-M` del próximo envío |
| `Run time` | tiempo de reloj |
| `CPU time` | tiempo de cómputo; si es ≪ `Run time × núcleos`, el trabajo esperaba disco |

Códigos de salida que han aparecido:

| código | causa típica |
|---|---|
| 1 | error del propio programa o del script |
| 2 | argumento o fichero de entrada mal |
| 127 | comando no encontrado (script ausente o entorno sin activar) |
| 132 | SIGILL: el binario usa instrucciones que el nodo no tiene |

### Ver el error

```bash
tail -20 "$LOGS_DIR"/smrdb.350740.err    # las ultimas 20 lineas de stderr
ls -lt "$LOGS_DIR" | head -20            # los logs mas recientes primero
```

### Actuar sobre trabajos

```bash
bkill 189069                  # matar uno (o un array entero por su JOBID)
bkill 189069[3]               # matar un elemento del array
bmod -n 4 -M 8000 189069      # cambiar recursos de uno PENDIENTE
bmod -Rn 189069               # borrar el -R antes de poner uno nuevo
```

**`bmod` no sustituye un `-R` existente** si no lo borras primero con `-Rn`:
cambia los núcleos pero deja la reserva de memoria vieja.

### Plantilla de envío

```bash
bsub -J "nombre[1-22]%8" \                      # array de 22, maximo 8 a la vez
     -o "$LOGS_DIR/nombre.%J.%I.log" \          # %J = jobid, %I = indice del array
     -e "$LOGS_DIR/nombre.%J.%I.err" \
     -q normal -n 4 -M 16000 \                  # cola, nucleos, memoria en MB
     -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" \  # todo en un nodo
     -env all \                                 # propaga las variables exportadas
     "$SCRIPTS_DIR/worker.sh arg1 arg2"
```

Sin `-env all`, variables como `ASM_DIRS`, `SPLITS` o `GENOMAD_EXTRA` **no
llegan al script**. Y los workers que exigen `LSB_JOBINDEX` necesitan array
aunque sea de un elemento: `-J "nombre[8]"`.

---

## 2. `check_step.sh` — verificación de arrays

```bash
./check_step.sh --list        # pasos disponibles
./check_step.sh trim          # un paso
./check_step.sh --all         # todo el pipeline de un vistazo
./check_step.sh btmap --exclude M41_S4,M12_S1
```

Comprueba cinco cosas: logs con `Exited`/`TERM_*`, elementos lanzados frente a
logs escritos, rastros de error en `.err`, muestras con entrada frente a
muestras con salida, y salidas vacías o sospechosamente pequeñas. Además
reporta memoria y tiempo reales.

Pasos: `cat fastqc polyg trim btmap unmap star spades quast genomad checkv diamond`.

**⚠ No fiarse de `genomad`, `checkv` ni `diamond`.** Esos tres bloques apuntan
todavía a `$SCRATCH/spades` con `*_final.fasta` y al esquema antiguo de un
trabajo por muestra. Las corridas reales fueron una sola sobre
`union_all/all_final.fasta`, y además hay un desajuste de claves (`M1_S1`
frente a `M1_S1_final`) que cuenta las muestras como fallidas. Para esos pasos,
usa los comandos de las secciones 7–10.

---

## 3. Directorios, espacio y ficheros

```bash
ls -la "$RESULTS_DIR"                          # contenido con tamaños y fechas
ls -lh "$REFS_DIR/diamond/"                    # -h: tamaños legibles (G, M)
ls -lt --time-style=+%d-%H:%M "$LOGS_DIR" | head   # ordenado por fecha, formato corto
ls -d "$SCRATCH"/spades* "$SCRATCH"/trinity*   # -d: el directorio, no su contenido
ls "$RESULTS_DIR"/kraken2/*.report | wc -l     # contar ficheros que casan
```

```bash
du -sh "$SCRATCH/mapped_votus/"                # cuanto ocupa un directorio
find "$RESULTS_DIR" -name "quality_summary.tsv" -printf '%s\t%p\n'   # buscar por nombre
find "$REFS_DIR/kraken2" -name "*.k2d"         # donde quedaron unos ficheros
stat -c '%y %n' "$SCRATCH"/union/*/ | sort     # fecha exacta de modificacion
```

**Cuota** — en Lustre la tuya sale `0k` (sin límite propio); la que manda es
la del grupo:

```bash
lfs quota -h -g team222 /lustre/scratch126 | tail -2
```

Columnas: `used quota limit ... files quota limit`. El límite de espacio del
grupo es 75 T.

---

## 4. FASTA y FASTQ

### Contar secuencias

```bash
grep -c '^>' contigs.fasta                     # numero de secuencias en un FASTA
cat "$SCRATCH"/union/*/*_final.fasta | grep -c '^>'   # total sobre varios ficheros
```

### Longitudes sin seqkit

```bash
awk '/^>/{ if(l>0){ if(l>=2500)a++; else b++ } l=0; next }
     { l+=length($0) }
     END{ if(l>0){ if(l>=2500)a++; else b++ }
          printf ">=2500: %d  <2500: %d\n", a, b }' fichero.fasta
```

- `/^>/{...; next}` — en cada cabecera, cierra la secuencia anterior y reinicia
- `{ l+=length($0) }` — suma la longitud de cada línea de secuencia
- `END{...}` — **no olvidar la última secuencia**, que no va seguida de cabecera

Funciona con FASTA multilínea, que es lo que falla si cuentas caracteres de una
sola línea.

### Longitud de lectura

```bash
zcat "$SCRATCH/trimmed/MERIDA/M11_S3_R1_001_paired.fastq.gz" | head -2 | tail -1 | tr -d '\n' | wc -c
```

- `head -2 | tail -1` — la segunda línea de un FASTQ es la secuencia
- `tr -d '\n'` — quita el salto final para que no cuente como carácter

Número de lecturas: `zcat fichero.fastq.gz | wc -l` dividido entre 4.

### Extraer secuencias por nombre

```bash
awk 'NR==FNR{ids[$0]=1;next} /^>/{keep=(substr($0,2) in ids)} keep' \
  lista_ids.txt entrada.fasta > salida.fasta
```

- `NR==FNR` — verdadero solo mientras se lee el primer fichero (la lista)
- `ids[$0]=1` — guarda cada nombre como clave
- `substr($0,2)` — la cabecera sin el `>`
- `keep` — imprime la cabecera y sus líneas de secuencia mientras sea verdadero

### seqkit

No está en ningún entorno conda. Hay módulo:

```bash
module load seqkit/2.8.2--h9ee0642_0
seqkit stats *.fasta                           # n, longitud total, min, media, max, N50
```

`module load` antepone al PATH y puede pisar binarios del entorno conda activo.

---

## 5. Tablas TSV con `awk`, `sort` y `join`

### Ver la cabecera numerada

```bash
head -1 tabla.tsv | tr '\t' '\n' | nl
```

Convierte las columnas en filas numeradas. Úsalo **siempre** antes de escribir
un `awk` por número de columna: así evitas sumar la columna equivocada (pasó con
`n_uscg` y `n_virus_hallmarks` en `features.tsv`).

### Localizar una columna por nombre

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="checkv_quality") q=i; next}
            {c[$q]++}
            END{for(k in c) printf "  %-18s %5d\n", k, c[k]}' quality_summary.tsv
```

- `NR==1{...}` — en la cabecera busca la columna por nombre y guarda su índice
- `{c[$q]++}` — cuenta cada valor distinto
- `END{...}` — imprime el recuento

Robusto frente a cambios de versión que reordenan columnas.

### Filtrar y contar

```bash
awk -F'\t' 'NR>1 && $7>=0.80' tabla.tsv | wc -l          # filas que cumplen una condicion
awk -F'\t' 'NR>1{print $9}' tabla.tsv | sort | uniq -c | sort -rn   # frecuencia de valores
```

### Mejor hit por consulta

```bash
sort -k1,1 -k10,10gr hits.tsv | awk -F'\t' '!s[$1]++'
```

- `sort -k1,1` — agrupa por la primera columna
- `-k10,10gr` — dentro de cada grupo, bitscore de mayor a menor (`g` numérico
  general, `r` inverso)
- `!s[$1]++` — imprime solo la primera vez que aparece cada consulta

### Unir dos tablas

```bash
export LC_ALL=C
join -t$'\t' -j1 <(sort -k1,1 a.tsv) <(sort -k1,1 b.tsv)
```

`join` exige las dos entradas ordenadas igual. Sin `LC_ALL=C` el `sort` usa
reglas de idioma y `join` avisa "input is not in sorted order" y **pierde
filas**.

### Comparar listas

```bash
comm -12 a.txt b.txt       # lineas en ambas
comm -23 a.txt b.txt       # solo en a
comm -13 a.txt b.txt       # solo en b
```

Mismo requisito: ambas ordenadas con `LC_ALL=C`.

### Leer bonito

```bash
column -t tabla.tsv | head -25                 # alinea columnas
cut -f1,2,7 tabla.tsv | head                   # solo algunas columnas
cut -c1-150 tabla.tsv                          # recortar lineas largas
```

---

## 6. CD-HIT — dereplicación y vOTUs

```bash
grep -c '^>' votus                             # numero de representantes (clusters)
grep -c '^>Cluster' votus.clstr                # lo mismo, desde el .clstr
```

Formato del `.clstr`:

```
>Cluster 0
0   12887nt, >M11_S3__rnaviral__NODE_1... at -/99.08%
3   13129nt, >M21_S29__megahit__k119_2559... *
```

- `*` marca el representante (el más largo)
- `at -/99.08%` — hebra (`+`/`-`) e identidad con el representante

Extraer los miembros de un cluster **sin romper con los puntos del nombre**:

```bash
sed -n '/^>Cluster 0$/,/^>Cluster 1$/p' votus.clstr \
  | sed -n 's/.*>\(.*\)\.\.\..*/\1/p'
```

El `\.\.\.` se ancla en los tres puntos que CD-HIT pone al final. Un patrón
`>[^.]+` corta en el primer punto del nombre (`cov_36.29`) y trunca las claves.

**CD-HIT agrupa, no mide.** Sus porcentajes fallan con contigs en hebras
opuestas o con distinto punto de inicio (el anfevirus de M1_S1 salía al 86.73% y
era 99.98% idéntico). Para afirmar identidad, BLAST o MAFFT.

---

## 7. geNomad

Salida principal en `genomad_all/all_final/all_final_summary/`:

```bash
G="$RESULTS_DIR/genomad_all/all_final"
S="$G/all_final_summary/all_final_virus_summary.tsv"

echo "virus: $(( $(wc -l < "$S") - 1 ))"      # -1 por la cabecera
head -1 "$S" | tr '\t' '\n' | nl              # columnas
```

Columnas de `virus_summary.tsv`: `1 seq_name  2 length  3 topology
4 coordinates  5 n_genes  6 genetic_code  7 virus_score  8 fdr  9 n_hallmarks
10 marker_enrichment  11 taxonomy`.

```bash
# reparto de scores
awk -F'\t' 'NR>1{ b=int($7*10)/10; n[b]++ } END{ for(x in n) printf "  %.1f: %d\n", x, n[x] }' "$S" | sort

# cuantos sin taxonomia
awk -F'\t' 'NR>1 && ($11=="" || $11=="Unclassified")' "$S" | wc -l

# top por longitud
awk -F'\t' 'NR>1{print $2"\t"$1"\t"$7"\t"$9"\t"$11}' "$S" | sort -rn | head
```

### Parámetros que usó de verdad

```bash
cat "$G"/*_summary/*_summary.json
```

Registra los filtros efectivos. Es la forma de confirmar que un parche llegó:
si `min_virus_hallmarks_short_seqs` sale `1` cuando pediste `0`, el
`GENOMAD_EXTRA` no se aplicó.

### Diagnóstico de un contig concreto

```bash
grep -F "NOMBRE" "$G/all_final_marker_classification/all_final_features.tsv" | cut -f1-8
grep -F "NOMBRE" "$G/all_final_annotate/all_final_genes.tsv"
```

En `features.tsv`: `2 n_genes  3 n_uscg  4 n_plasmid_hallmarks
5 n_virus_hallmarks`. **La columna 3 son marcadores bacterianos, no virales.**

`mmseqs2.tsv` vacío en `annotate/` significa cero hits contra la base de
marcadores.

---

## 8. CheckV

```bash
Q="$RESULTS_DIR/checkv_viral/all_final_virus.fna/quality_summary.tsv"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="checkv_quality")q=i;next}{c[$q]++}
            END{for(k in c) printf "  %-18s %5d\n", k, c[k]}' "$Q"
```

Ojo con la ruta: el worker nombra el directorio con el basename **completo**
del FASTA, extensión incluida (`all_final_virus.fna/`).

Columnas clave: `2 contig_length  6 viral_genes  7 host_genes
8 checkv_quality  10 completeness  12 contamination  13 kmer_freq  14 warnings`.

```bash
# los buenos
awk -F'\t' '$8=="Complete" || $8=="High-quality"{print $1, $2, $10, $13}' "$Q"
```

**`kmer_freq`**: 1.0 = cada tramo es único. 1.2 = un 20% duplicado. Descarta
contigs "completos" que en realidad repiten secuencia (el de Trinity de
15,488 bp).

---

## 9. samtools — mapeos

```bash
samtools flagstat -@ 4 muestra.bam            # resumen: total, mapeadas, pares
samtools idxstats muestra.sorted.bam          # lecturas por referencia (requiere indice)
samtools coverage muestra.sorted.bam          # amplitud y profundidad por referencia
```

`idxstats` y `coverage` necesitan el BAM **ordenado e indexado**:

```bash
samtools sort -@ 4 -m 1G -T "$tmp/sort" -o ordenado.bam entrada.bam
samtools index ordenado.bam
```

Pon el temporal (`-T`) en Lustre, no en `/tmp` del nodo: con BAM de varios GB
el `/tmp` se llena y el trabajo muere sin explicación clara.

Columnas de `samtools coverage`: `1 rname  3 endpos (longitud)  5 covbases
6 coverage (%)  7 meandepth`.

**Contar lecturas no es detectar.** Con "≥1 lectura" como umbral, las 748 vOTUs
salían en las 22 muestras. Presencia = amplitud (`coverage ≥ 70%`), no volumen.

### Porcentaje al huésped

```bash
grep " mapped (" muestra_flagstat.txt | grep -v primary
```

---

## 10. DIAMOND

### Verificar una base

```bash
diamond dbinfo -d "$REFS_DIR/diamond/rvdb.dmnd"
```

Da número de secuencias y letras. Si falla, la versión de DIAMOND no lee ese
`.dmnd`.

### Resumen de una corrida

```bash
D="$RESULTS_DIR/diamond_votus_rvdb_final.tsv"
echo "lineas: $(wc -l < "$D")"
echo "consultas con hit: $(cut -f1 "$D" | sort -u | wc -l)"
awk -F'\t' '{s+=$3; n++} END{printf "identidad media: %.1f%%\n", s/n}' "$D"
```

Las columnas dependen del `--outfmt` que pediste. En `diamond_votus_rvdb_final.tsv`:
`1 qseqid  2 sseqid  3 pident  4 length  5 qlen  6 slen  7 qcovhsp  8 scovhsp
9 evalue  10 bitscore  11 stitle`.

### Leer el organismo

```bash
# RVDB: la especie va entre corchetes al final de stitle
awk -F'\t' '{d=$11; sub(/.*\[/,"",d); sub(/\]$/,"",d); print d}' "$D" | sort | uniq -c | sort -rn | head

# nr: igual, pero cuidado con nombres de virus que contienen el del huesped
awk -F'\t' '$2 ~ /virus|phage|viridae/ {next} $2 ~ /Aedes|Culex/ {print $1}' best_nr.tsv
```

**`Aedes anphevirus` contiene "Aedes".** Un filtro por nombre de huésped sin
excluir antes los virus los cuenta como mosquito.

### Qué mirar en un hit

| señal | lectura |
|---|---|
| `pident ≥ 90` | especie descrita o muy cercana |
| `pident 30–50` | pariente lejano; el nombre es orientativo |
| `qcovhsp` alto | el contig entero es esa proteína |
| `scovhsp` alto | se recuperó la proteína de referencia completa |
| `slen < 150` | referencia corta: un `scov` alto no significa nada |
| `evalue > 1e-3` | ruido |

`qcovhsp` penaliza los genomas largos (el anfevirus de 13 kb daba 37–47%). Usa
`qcov ≥ 50` **o** `scov ≥ 50`.

---

## 11. Kraken2 y Bracken

```bash
K="$RESULTS_DIR/kraken2"
ls "$K"/*.report | grep -v bracken | wc -l    # reports de Kraken2 (sin los de Bracken)
```

Columnas del `.report`:
`1 %  2 lecturas_clado  3 lecturas_directas  4 rango  5 taxid  6 nombre`.
Con `--report-minimizer-data` hay dos columnas más en medio y el nombre pasa a
la 8.

Rangos: `U` sin clasificar, `R` raíz, `D` dominio, `P` filo, `C` clase,
`O` orden, `F` familia, `G` género, `S` especie.

```bash
# % sin clasificar
awk -F'\t' '$4=="U"{print $1}' "$K/M7_S15.report"

# top especies
awk -F'\t' '$4=="S"' "$K/M7_S15.report" | sort -t$'\t' -k2 -rn | head -15

# buscar un organismo
grep -i "wolbachia" "$K"/*.report | grep -v bracken
```

**PlusPF no incluye insectos.** Las lecturas de *Aedes* salen sin clasificar,
salvo las de rRNA, que comparten k-meros con el rRNA humano de la base y salen
como *Homo sapiens*. Ese porcentaje es una medida indirecta de rRNA, no
contaminación humana.

Bracken: `<sample>.bracken` es la tabla de abundancia reestimada;
`<sample>.bracken.report` es el árbol corregido.

---

## 12. BLAST y MAFFT

```bash
module load blast/2.15.0--pl5321h6f7f691_1     # o: conda activate blast_2.17.0

blastn -query a.fasta -subject b.fasta \
       -outfmt "6 qseqid sseqid pident length qstart qend sstart send"
```

`sstart > send` significa que alinea en la hebra opuesta. Dos bloques que suman
la longitud completa suelen ser un indel, no una permutación.

```bash
mafft --adjustdirection --auto --thread 4 entrada.fasta > salida.aln 2> mafft.log
grep -c '^>' salida.aln                        # cuantas secuencias alineadas
grep '^>_R_' salida.aln                        # las que MAFFT invirtio
```

---

## 13. Entornos

```bash
conda env list                                 # todos los entornos
conda env list | grep -i diamond               # uno concreto
conda activate diamond_2.2.6
diamond --version                              # confirmar que el binario responde
module avail 2>&1 | grep -i -E 'blast|seqkit'  # modulos de la farm
```

`2>&1` es necesario en `module avail` porque escribe en stderr.

---

## 14. Trampas que ya costaron tiempo

| síntoma | causa |
|---|---|
| `join` pierde filas | falta `export LC_ALL=C` |
| `grep "x.*RUN"` no encuentra nada | en `bjobs` el `STAT` va antes del nombre |
| `awk` con dos `END` no imprime | usar un solo `END`, o varios `awk` separados |
| fichero de 0 bytes tras descargar | URL caducada o mal escrita |
| `mapping_stats/` sobrescrito | el worker de bowtie2 usa siempre el mismo directorio |
| dos arrays escribiendo lo mismo | revisar `bjobs` antes de reenviar |
| una columna "de marcadores" que no lo es | numerar la cabecera con `nl` antes de filtrar |
| un parámetro que no llegó | confirmar en el `*_summary.json` o en el log |
| conteo de lecturas absurdo | rRNA en el catálogo, o umbral de ≥1 lectura |
| "virus" con nombre de huésped | excluir `virus|phage|viridae` antes de filtrar por organismo |
