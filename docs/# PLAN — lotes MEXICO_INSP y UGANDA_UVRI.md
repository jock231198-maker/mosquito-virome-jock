# PLAN — lotes MEXICO_INSP y UGANDA_UVRI hasta ensamblaje

Runbook para pasar 6 muestras nuevas por `QC → polyG → trimming → depleción de
huésped → SPAdes rnaviral`, reutilizando los workers de MERIDA sin tocarles una
línea.

Alcance decidido: **hasta ensamblaje**. El downstream se decide con los contigs
delante (§11).

**Estructura decidida: un árbol por origen, igual que MERIDA.** No hay árbol
compartido: cada lote tiene su `work/` y su `results/`, su `samples.txt` y su
propio recuento en `check_step.sh`. Se corre la misma cadena dos veces,
cambiando una línea.

---

## 1. Los tres lotes

| | MERIDA | MEXICO_INSP | UGANDA_UVRI |
|---|---|---|---|
| Muestras | 22 | 3 | 3 |
| | `M1_S1`… | `IgualaC2_S2`, `PalmiraC1_S1`, `XochiC3_S3` | `UVRI2_S3`, `UVRI3_S2`, `UVRI5T_S1` |
| Taxón | *Ae. aegypti* | *Ae. aegypti aegypti* | *Ae. aegypti formosus* |
| Contexto | campo, Yucatán | colonias de laboratorio (INSP) | Uganda |
| Datos | crudos | **crudos**, 1 lane | **ya recortados** (Trim Galore) |
| Entra por | polyG | Trimmomatic | bowtie2 |
| `work/` | `work/MERIDA` | `work/MEXICO_INSP` | `work/UGANDA_UVRI` |

### Por qué el nombre lleva país *y* origen

Un `work/MEXICO` a secas chocaría con MERIDA, que también es México, y el día
que llegue otro lote mexicano de campo el nombre ya estaría gastado. País
primero para que ordene por país, origen después para que no colisione. Si
prefieres otra cosa, es el `case` de `env_lote.sh`: una línea.

### Dos metadatos que conviene cerrar antes de escribir métodos

- **Iguala es Guerrero y Xochitepec es Morelos.** Las colonias del INSP suelen
  llevar el nombre de la localidad de origen aunque se mantengan en otra sede.
  "Morelos" describe dónde está el insectario, no de dónde salió cada colonia.
- **`LANCIS_Aaa_Aas_piloto` dice `Aaa_Aas`**, dos taxones, pero solo listaste
  `INSP_Aaa_Labcolonies`. Si hay un `Aas` al lado, mira qué es antes de decidir
  si entra.

---

## 2. Cinco avisos metodológicos, por delante

Ninguno bloquea el pipeline. Los cinco cambian cómo se interpreta el resultado,
y son más baratos de anotar ahora que de descubrir en la discusión.

**a) El índice de bowtie2 es AaegL5 (cepa LVP) y las UVRI son *formosus*.**
*Ae. aegypti formosus* es la forma africana ancestral; LVP es una cepa doméstica
de laboratorio. Más divergencia respecto a la referencia = **menos huésped
eliminado** = `host_pct` más bajo y `unmapped_fastq` más grande en UGANDA_UVRI
**sin que eso sea señal biológica**. Es el mismo sesgo que documentaste en
bowtie2-vs-STAR (§6 del HANDOFF), apuntando al otro lado.

Consecuencia: **el `host_pct` no es comparable entre lotes**, y el recuento de
contigs normalizado por pares no-mapeados tampoco del todo. Lo que sí compara es
qué virus aparecen. Si quieres el número duro, mapea un subconjunto de UVRI con
`--local` y mira cuánto sube el `host_pct`: esa diferencia es el sesgo.

**b) Las INSP son colonias de laboratorio.** Espera un viroma pobre: sin
exposición de campo persisten sobre todo virus específicos de insecto de
transmisión vertical (CFAV, PCLV y similares). **No es un fallo del experimento,
es el control.** Y es lo que le da valor a la comparación: si el anfevirus
Xinmoviridae de ~12.9 kb que tienes en 12 de 22 muestras de MERIDA aparece
también en las colonias, tienes un argumento de transmisión vertical; si no
aparece en ninguna, uno de adquisición ambiental.

**c) Las UVRI están recortadas con otra herramienta y otros parámetros.**
Trim Galore corta por defecto a `--length 20`; tu cadena usa `MINLEN:35`. Reads
más cortas mapean peor y ensamblan peor. Dos salidas honestas: re-trimar desde
crudos si aparecen, o **homogeneizar el mínimo de longitud** (paso 2b, cinco
minutos). Lo que no vale es comparar rendimiento de ensamblaje sin mencionarlo.

**d) `2_S3`, `3_S2`, `5T_S1` empiezan por dígito.** `make.names()` en R les mete
una X delante y awk los coerciona a número. El staging los renombra a `UVRI2_S3`,
`UVRI3_S2`, `UVRI5T_S1` y deja la correspondencia en `staging_map.tsv`.

**e) Los dos lotes son MiSeq, y son MUCHO menos profundos que MERIDA.**
Medido en el staging (21 sep):

| lote | instrumento | flowcell | pares por muestra |
|---|---|---|---|
| MERIDA | 2 colores (el polyG era la sobrerrepresentada nº1 en 23/23) | — | 8–64 M, media ~43 M |
| MEXICO_INSP | MiSeq `M07836` | `000000000-K34R3` | 4.4 / 5.8 / 12.0 M |
| UGANDA_UVRI | MiSeq `M02853` | `000000000-L8TL6` | 8.9 / 10.7 / 11.6 M |

Cambia dos cosas, las dos de interpretación:

- **La ausencia de un virus en estos lotes vale poco.** Con 3-10× menos
  profundidad, no detectar algo es tan compatible con "no está" como con "está
  por debajo del límite de detección". La presencia sí vale; la ausencia hay que
  escribirla con cuidado, sobre todo en el argumento de las colonias del aviso b.
- **El recuento de contigs no se compara crudo con MERIDA.** Normalizado por
  millón de pares no-mapeados, como ya hiciste en la §8, sí.

Y una técnica: MiSeq da reads más largas que NovaSeq. Si el FastQC las devuelve
a 250-300 bp, hay que mirar los k-meros de SPAdes antes del paso 6 — los valores
por defecto se eligen a partir de la longitud de read, y `spades_parametros.md`
está escrito para las de MERIDA. Tres muestras de INSP comparten flowcell: si
aparece efecto de lote entre ellas, ahí está el sospechoso.

---

## 3. Cómo se cambia de lote

```bash
source "$S/env_lote.sh" MEXICO_INSP     # o UGANDA_UVRI, o MERIDA
```

`env_lote.sh` asigna `SCRATCH`, `RESULTS_DIR`, `LOGS_DIR`, `RESULT_FROM` y
`DATA_DIR` **sin `${VAR:-...}`**, a propósito: así puedes cambiar de lote en la
misma terminal sin arrastrar el valor viejo. Es la trampa de la §13, desactivada
para este caso concreto.

El árbol que queda:

```
$USER_ROOT/work/MEXICO_INSP/
    cat_fastq/            symlinks a las 3 crudas
    trimmed/MEXICO_INSP/  Trimmomatic (sin paso de fastp: ver paso 2)
    unmapped_fastq/       <- entrada de SPAdes
    spades/<sample>/
$USER_ROOT/work/UGANDA_UVRI/
    trimmed/UGANDA_UVRI/  symlinks a las val_* (no pasan por fastp ni Trimmomatic)
    unmapped_fastq/
    spades/<sample>/
```

Misma forma que `work/MERIDA/trimmed/MERIDA/`. Cada `check_step.sh <paso>` cuenta
**3/3** dentro de su árbol, sin mezclas ni exclusiones.

> **Lo que cuesta separar**: SPAdes son dos arrays de 3 en vez de uno de 6, y
> cuando llegue el downstream habrá que construir la lista de FASTA cruzando los
> dos árboles. Eso último no es problema: `genomad_worker.sh`, `checkv_worker.sh`
> y `quast_worker.sh` ya toman la lista y el `outbase` como argumentos.

---

## Paso 0 — acceso a scratch127 y despliegue

Los datos viven en **otro proyecto (`vector_eve`) y otro filesystem
(`scratch127`)**. Que se lean desde un nodo de cómputo, no solo desde el de
login:

```bash
V=/lustre/scratch127/tol/projects/vector_eve/data/MosquitoTranscriptome
ls -l "$V/INSP_Aaa_Labcolonies" "$V/UVRI_Aaf" 2>&1 | head -20

bsub -J "acceso[1-1]" -o /tmp/acceso.%J.log -q normal -M 500 \
     -R "select[mem>500] rusage[mem=500]" \
     "ls -l $V/UVRI_Aaf && gzip -t $V/UVRI_Aaf/2_S3_L001_R1_001_val_1.fq.gz && echo LEGIBLE_DESDE_COMPUTO"
```

Si el nodo de cómputo no monta `lus127`, todo lo demás falla en el segundo 1 con
el `.err` vacío. Compruébalo aquí, no dentro de un array de 12.

Despliegue:

```bash
S=~/jkvirome/jm79/mosquito-virome-jock/scripts/sanger_scripts_jobs/last_update_05_08_2026
cp env_lote.sh prep_lote.sh dup_unmapped.sh "$S/"
chmod +x "$S"/{prep_lote.sh,dup_unmapped.sh}
sed -i 's/\r$//' "$S"/{env_lote.sh,prep_lote.sh,dup_unmapped.sh}     # el exit 126

grep -n 'ENV_FASTP' "$S/config.sh" || echo "!! config.sh NO define ENV_FASTP y polyg_worker.sh lo usa"
```

Ese `grep` importa: `polyg_worker.sh` activa `$ENV_FASTP`, y la copia de
`config.sh` que hay en el proyecto **no define esa variable**. Si en la farm
tampoco está, `activate_env ""` y el paso 2 muere sin decir por qué.
`env_lote.sh` avisa también al sourcearse.

Staging, un lote cada vez:

```bash
source "$S/env_lote.sh" MEXICO_INSP  && "$S/prep_lote.sh" --count
source "$S/env_lote.sh" UGANDA_UVRI  && "$S/prep_lote.sh" --count
```

`prep_lote.sh` no copia nada: crea symlinks con el nombre exacto que esperan los
workers, escribe `$SCRATCH/samples.txt` (mismo nombre y formato que MERIDA),
imprime **el instrumento de secuenciación de cada fichero** (paso 1) y verifica
que R1 y R2 llevan el mismo ID en el primer read — el chequeo de la §5 de
`check_inputs.sh`, el único que detecta un cruce de parejas. Sale con código 1
si algo no cuadra.

---

## Paso 1 — FastQC

### La decisión del polyG: RESUELTA (21 sep), no hace falta fastp

`prep_lote.sh` leyó el ID del instrumento de la cabecera del primer read:
**`M07836` en INSP y `M02853` en UVRI, los dos MiSeq.** Química de 4 colores: el
MiSeq llama base con dos imágenes por ciclo y una G no es "ausencia de señal",
así que **el artefacto de polyG no puede producirse**. No es que Trim Galore lo
quitara: es que nunca estuvo.

Consecuencias:

- **El paso 2 (fastp) se salta en los dos lotes.** En MERIDA hacía falta porque
  aquello era química de dos colores (el polyG era la secuencia
  sobrerrepresentada nº1 en 23/23 muestras). Aquí sería un no-op que duplica
  2.5 GB y añade un paso que puede fallar.
- Trimmomatic toma su entrada de `cat_fastq/`, no de `nopolyg/`.
- La evidencia para métodos la da el FastQC de abajo: cero polyG en
  *Overrepresented sequences*. Es la misma prueba que daría el JSON de fastp.

Queda abierta solo la mitad del paso 2b: **si las UVRI traen reads por debajo de
35 bp**, hay que igualarlas al `MINLEN:35` de tu Trimmomatic. Eso lo decide la
*Sequence Length Distribution* del FastQC.

### Correr FastQC, por lote

```bash
source "$S/env_lote.sh" MEXICO_INSP
find "$SCRATCH/cat_fastq" -type f -name "*.fastq.gz" | sort > "$SCRATCH/filelist.txt"
N=$(wc -l < "$SCRATCH/filelist.txt"); echo "N=$N"        # 6

bsub -J "fqc[1-$N]%6" \
     -o "$LOGS_DIR/fqc.%J.%I.log" -e "$LOGS_DIR/fqc.%J.%I.err" \
     -q normal -n 4 -M 1000 -R "select[mem>1000] rusage[mem=1000] span[hosts=1]" \
     "SCRATCH=$SCRATCH RESULTS_DIR=$RESULTS_DIR $S/fastqc_worker.sh $SCRATCH/filelist.txt $RESULTS_DIR/fastqc_entrada"
```

Para UGANDA_UVRI, lo mismo con `find "$SCRATCH/trimmed/$RESULT_FROM"`.

(`-M 1000`: en MERIDA FastQC consumió 279 MB y se pidieron 2000. Pedir de más
bloquea la cola tanto como pedir de menos.)

Qué mirar en el MultiQC, y qué decide cada cosa:

| Qué | Lote | Si falla |
|---|---|---|
| *Overrepresented sequences* sin polyG | los dos | confirma el diagnóstico de MiSeq; si apareciera polyG, replantear |
| *Sequence Length Distribution* con masa < 35 bp | UVRI | paso 2b (`--length_required 35`) |
| Longitud de read (¿150? ¿250? ¿300?) | los dos | si ≥ 250, revisar los k-meros de SPAdes antes del paso 6 |
| *Adapter content* | UVRI | si falla, Trim Galore no hizo su trabajo: replantear el lote |
| *Overrepresented sequences*, las que no son polyG | los dos | candidatas a rRNA, igual que en la §9 de MERIDA |

Criterio, escrito antes de ver el dato para no racionalizar después:

- **longitud mínima ≥ 35** → UVRI entra directo a bowtie2 (paso 4).
- **masa por debajo de 35 bp** → paso 2b.

---

## Paso 2 — polyG: NO SE CORRE

Ver la cabecera del paso 1. Los dos lotes son MiSeq y el artefacto no puede
existir. Trimmomatic toma su entrada de `cat_fastq/`.

Para material y métodos: *"polyG trimming no se aplicó a estos lotes: ambos se
secuenciaron en MiSeq (química de 4 colores), donde el artefacto no se produce.
En MERIDA, secuenciado con química de 2 colores, sí fue necesario."*

### Paso 2b — solo si el FastQC lo pide (UGANDA_UVRI)

fastp haciendo dos cosas y nada más: quitar polyG y tirar lo que quede por debajo
de 35 bp, para igualar el `MINLEN:35` de Trimmomatic. Ni adaptadores ni calidad:
eso ya lo hizo Trim Galore y pisarlo sería recortar dos veces.

```bash
source "$S/env_lote.sh" UGANDA_UVRI
IN="$SCRATCH/trimmed/$RESULT_FROM"; OUT="$SCRATCH/trimmed/${RESULT_FROM}_len35"
mkdir -p "$OUT" "$RESULTS_DIR/polyg_reports"

while read -r s; do
  bsub -J "uvrifix_$s" -o "$LOGS_DIR/uvrifix.$s.%J.log" -e "$LOGS_DIR/uvrifix.$s.%J.err" \
       -q normal -n 4 -M 4000 -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
       "module load conda && conda activate $ENV_FASTP && fastp \
          -i $IN/${s}_R1_001_paired.fastq.gz -I $IN/${s}_R2_001_paired.fastq.gz \
          -o $OUT/${s}_R1_001_paired.fastq.gz -O $OUT/${s}_R2_001_paired.fastq.gz \
          --trim_poly_g --poly_g_min_len 8 --length_required 35 \
          --disable_adapter_trimming --disable_quality_filtering \
          --thread 4 --compression 6 \
          --json $RESULTS_DIR/polyg_reports/${s}.json --html $RESULTS_DIR/polyg_reports/${s}.html"
done < "$SCRATCH/samples.txt"
```

Si corres esto, el `indir` del paso 4 para UVRI pasa a ser `trimmed/UGANDA_UVRI_len35`.

> Los nombres de salida **acaban en `.gz`** a propósito: fastp decide si comprime
> mirando la extensión, y un nombre sin `.gz` deja FASTQ en claro que luego no
> pasa `gzip -t`.

---

## Paso 3 — trimming (solo MEXICO_INSP)

```bash
source "$S/env_lote.sh" MEXICO_INSP
N=$(wc -l < "$SCRATCH/samples.txt")

bsub -J "trim[1-$N]%3" \
     -o "$LOGS_DIR/trim.%J.%I.log" -e "$LOGS_DIR/trim.%J.%I.err" \
     -q normal -n 4 -M 12000 -R "select[mem>12000] rusage[mem=12000] span[hosts=1]" \
     "SCRATCH=$SCRATCH RESULTS_DIR=$RESULTS_DIR JAVA_MEM=8g $S/trimmomatic_worker.sh $SCRATCH/samples.txt $SCRATCH/cat_fastq"
```

**La entrada es `cat_fastq`, no `nopolyg`**, porque el paso 2 no se corre. No
hacen falta los argumentos 3 y 4: el worker usa `$SCRATCH/trimmed` y
`$RESULT_FROM`, que ya apuntan a este lote.

Mismos parámetros que MERIDA (`ILLUMINACLIP:2:30:10 LEADING:3 TRAILING:3
SLIDINGWINDOW:4:25 MINLEN:35`) — que es justo lo que hace comparables INSP y las
22. `JAVA_MEM` por debajo del `-M`, siempre. MERIDA usó 20000/16g; 12000/8g sobra
para 3 muestras y entra antes en una cola al 3.8%. Si se queja del heap, vuelve a
20000/16g.

---

## Paso 4 — depleción de huésped (los dos lotes, por separado)

```bash
"$S/preflight.sh" btmap

# --- MEXICO_INSP
source "$S/env_lote.sh" MEXICO_INSP
"$S/check_step.sh" trim
N=$(wc -l < "$SCRATCH/samples.txt")
bsub -J "btmap[1-$N]%3" \
     -o "$LOGS_DIR/btmap.%J.%I.log" -e "$LOGS_DIR/btmap.%J.%I.err" \
     -q normal -n 16 -M 12000 -R "select[mem>12000] rusage[mem=12000] span[hosts=1]" \
     "SCRATCH=$SCRATCH RESULTS_DIR=$RESULTS_DIR $S/bowtie_map_worker_fast.sh $SCRATCH/samples.txt $SCRATCH/trimmed/$RESULT_FROM"

# --- UGANDA_UVRI   (…_len35 si corriste el paso 2b)
source "$S/env_lote.sh" UGANDA_UVRI
N=$(wc -l < "$SCRATCH/samples.txt")
bsub -J "btmap[1-$N]%3" \
     -o "$LOGS_DIR/btmap.%J.%I.log" -e "$LOGS_DIR/btmap.%J.%I.err" \
     -q normal -n 16 -M 12000 -R "select[mem>12000] rusage[mem=12000] span[hosts=1]" \
     "SCRATCH=$SCRATCH RESULTS_DIR=$RESULTS_DIR $S/bowtie_map_worker_fast.sh $SCRATCH/samples.txt $SCRATCH/trimmed/$RESULT_FROM"
```

**Ventaja de haber separado los árboles**: estos dos arrays escriben en
`unmapped_fastq/` y `.btmap_stage/` distintos, así que se pueden lanzar a la vez
sin repetir la colisión que mató los arrays 329577 y 330925.

> El paso que hay que comprobar después es **`unmap`, no `btmap`**: el worker
> *fast* usa `--un-conc-gz` y `-S /dev/null`, no escribe BAM, y el nombre de
> salida es justo el que busca `check_step.sh unmap`. `check_step.sh btmap`
> dirá 0/3 y estará en lo cierto.

Lo primero que mirar al terminar, con el aviso 2a en la cabeza:

```bash
grep -H "overall alignment rate" "$RESULTS_DIR"/mapping_stats/*_bowtie2_summary.txt
```

**Si UVRI sale sistemáticamente por debajo de INSP, la primera hipótesis es la
divergencia *formosus* frente a LVP, no la biología.**

---

## Paso 5 — duplicación en `unmapped_fastq`

Por lote:

```bash
source "$S/env_lote.sh" MEXICO_INSP && "$S/check_step.sh" unmap && "$S/dup_unmapped.sh"
source "$S/env_lote.sh" UGANDA_UVRI && "$S/check_step.sh" unmap && "$S/dup_unmapped.sh"
```

Réplica del conteo de la §9: secuencias únicas en las primeras 400,000 reads.
Referencia MERIDA: cinco muestras entre 19.5% y 37.0%, el resto entre 37.3% y
68.1%. Por debajo de ~37% es la zona donde el rRNA 28S residual copaba la
librería y el ensamblaje se venía abajo pese a sobrar material.

El script escribe también `top_secuencias_unmapped.fasta` con las 3 secuencias
más frecuentes de cada muestra — el fichero que alimenta el test de SortMeRNA de
la §11 cuando retomes esa rama.

**Lee el `%` de únicas junto a la longitud media que imprime el script.** Reads
más cortas colisionan más: si UVRI sale baja pero también sale corta, parte de
esa duplicación es aritmética, no biología.

---

## Paso 6 — SPAdes `--rnaviral`

Dos arrays de 3. En MERIDA `--rnaviral` dio **mediana 1913 MB y máximo 2077 MB**,
mediana 663 s y máximo 1703 s, sobre muestras de hasta 11.2 M pares. Los 32 GB
que se pidieron iban quince veces sobrados.

```bash
for L in MEXICO_INSP UGANDA_UVRI; do
  source "$S/env_lote.sh" "$L"
  "$S/check_step.sh" unmap || { echo "$L: unmap incompleto, no lanzo"; continue; }
  cp "$SCRATCH/samples.txt" "$SCRATCH/asm_samples.txt"
  N=$(wc -l < "$SCRATCH/asm_samples.txt")
  echo "SPADES_MEM=8 $S/assembly_spades.sh $SCRATCH/asm_samples.txt $SCRATCH/unmapped_fastq"
  bsub -J "spades[1-$N]%3" \
       -o "$LOGS_DIR/spades.%J.%I.log" -e "$LOGS_DIR/spades.%J.%I.err" \
       -q normal -n 4 -M 12000 -R "select[mem>12000] rusage[mem=12000] span[hosts=1]" \
       "SCRATCH=$SCRATCH RESULTS_DIR=$RESULTS_DIR SPADES_MEM=8 $S/assembly_spades.sh $SCRATCH/asm_samples.txt $SCRATCH/unmapped_fastq"
done
```

Ese `echo` antes del `bsub` es el hábito que caza los cuatro fallos de arranque
de cero segundos: imprime la cadena exacta con las variables ya expandidas.

`SPADES_MEM` siempre por debajo del `-M`: si SPAdes topa con su propio límite se
para y lo dice; si topa con el de LSF, lo matan sin explicación.

Salida esperada por muestra: `contigs.fasta`, `contigs_1000bp.fasta`,
`<sample>_final.fasta` y `$RESULTS_DIR/assembly_stats/<sample>_assembly.tsv`.

---

## 10. Puertas de verificación

La regla de la §17: **confirmar que el paso anterior terminó, no que sus ficheros
existan.** Todo con el lote sourceado; `check_step.sh` hereda `$SCRATCH` y cuenta
solo ese árbol.

| Antes de | Correr | MEXICO_INSP | UGANDA_UVRI |
|---|---|---|---|
| paso 1 | `prep_lote.sh` (sale 1 si falla) | 6 symlinks, 3 parejas OK | ídem |
| paso 3 | `check_inputs.sh --deep` | layout sano en `cat_fastq` | n/a |
| paso 4 | `check_step.sh trim` | 3/3 | n/a (entran ya recortadas) |
| paso 5 | `check_step.sh unmap` | 3/3 | 3/3 |
| paso 6 | `check_step.sh spades` | 3/3 | 3/3 |

`check_step.sh btmap` dará 0/3 en los dos: el worker *fast* no escribe BAM. Es
correcto, no es un fallo.

`read_accounting.sh` sobre UGANDA_UVRI dará `trimmed_pairs = NA`: no hay log de
Trimmomatic porque no pasaron por él. También correcto.

---

## 11. Lo que este plan NO hace, y por qué

- **Sin depleción de rRNA.** SortMeRNA sigue con el SIGILL sin resolver, y estas
  6 van por el mismo camino que las 22 para que sean comparables. El paso 5 mide
  si el problema se repite; si se repite, entra en la misma cola de espera que
  MERIDA y se reensambla todo junto, no estas por separado.
- **Sin `rna` ni `metaviral`.** La rama principal es `rnaviral`. Añadir modos
  antes de saber cuántos contigs largos salen es trabajo a ciegas.
- **Sin downstream.** geNomad, CheckV y los vOTUs se deciden con los contigs
  delante. La pregunta de fondo entonces será si las 6 entran en un **único
  catálogo de vOTUs junto a las 22** —lo que permite decir "este virus está en
  Mérida y en Kampala pero no en la colonia"— o si van aparte. La primera opción
  obliga a rehacer derep + geNomad + CheckV + CD-HIT sobre la unión y a remapear
  las 28. Es coste real; conviene presupuestarlo antes, no después. Que los
  árboles estén separados no lo estorba: los workers de aguas abajo toman la
  lista de FASTA y el `outbase` como argumentos.

---

## 12. Registro de decisiones (para material y métodos)

1. Un árbol de trabajo y resultados por origen: `MERIDA`, `MEXICO_INSP`,
   `UGANDA_UVRI`. Sin directorios compartidos entre lotes.
2. Staging por symlink, sin copia: los crudos de `vector_eve` (scratch127) no se
   duplican ni se tocan.
3. MEXICO_INSP: Trimmomatic `ILLUMINACLIP:2:30:10 LEADING:3 TRAILING:3
   SLIDINGWINDOW:4:25 MINLEN:35` → bowtie2 vs AaegL5 → SPAdes `--rnaviral`,
   mismos parámetros que MERIDA.
3b. Sin recorte de polyG en ninguno de los dos lotes: ambos se secuenciaron en
   MiSeq (4 colores, `M07836` y `M02853`), donde el artefacto no se produce.
   MERIDA sí lo necesitó por ser química de 2 colores.
3c. Profundidad muy inferior a MERIDA (4.4–12.0 M pares por muestra frente a una
   media de ~43 M): la ausencia de un virus en estos lotes no es evidencia
   fuerte de ausencia.
4. UGANDA_UVRI: se aceptan las reads recortadas con Trim Galore por el proveedor;
   la asimetría de recorte queda documentada y, si el paso 1 lo indica, se mitiga
   igualando el mínimo de longitud a 35 bp.
5. Muestras UVRI renombradas `2_S3 → UVRI2_S3` etc.; correspondencia en
   `results/UGANDA_UVRI/qc_control/staging_map.tsv`.
6. Sin depleción de rRNA, por consistencia con las 22 de MERIDA.
7. El `host_pct` no se compara entre lotes: la divergencia *formosus* / LVP lo
   sesga a la baja en UGANDA_UVRI.