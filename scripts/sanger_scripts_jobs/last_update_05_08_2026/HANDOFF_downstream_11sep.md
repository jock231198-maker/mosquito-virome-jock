# Handoff — pasos 4 a 9 del downstream (4–11 sep 2026)

Continuación de `PLAN_downstream.md`, que quedó con los pasos 4 y 5 en "smoke
pendiente" y 6–10 "por escribir". Este documento recoge lo hecho desde entonces,
los números medidos y las decisiones que hay que poder defender.

---

## Estado por paso

| # | paso | estado |
|---|---|---|
| 1 | Dereplicación por muestra | **HECHO** — rehecha con Trinity, sin `meta`: **6,004** contigs |
| 2 | Listas por conjunto | HECHO |
| 3 | QUAST comparativo | HECHO 22/22 |
| 4 | geNomad | **HECHO** — 1,144 virales de 6,004, con controles |
| 5 | CheckV | **HECHO** — las dos vías (sobre contigs y sobre salida viral) |
| 6 | Filtro de calidad | **criterio fijado** (ver abajo), pendiente de limpiar rRNA |
| 7 | Catálogo vOTU | **HECHO** — 748 vOTUs |
| 8 | Mapeo de vuelta | BAM y cuentas hechas, **resultado contaminado por rRNA** |
| 9 | DIAMOND | **HECHO** — 20 vOTUs firmes, 6 especies |
| 10 | Tabla final | pendiente |

---

## Paso 1 rehecho — la unión cambió

La unión original (4,667 contigs) se hizo **antes de que Trinity terminara** y
**antes de decidir sacar `meta`**. Se rehízo el 7 de septiembre.

```bash
rm -rf "$SCRATCH/union"
export ASM_DIRS="$SCRATCH/spades $SCRATCH/spades_rna $SCRATCH/spades_metaviral $SCRATCH/megahit $SCRATCH/trinity"
bsub ... -env all "$SCRIPTS_DIR/derep_contigs.sh $SCRATCH/asm_samples.txt $SCRATCH/union"
```

Dos cosas imprescindibles ahí: el `rm -rf` (el script sale con "Ya existe... No se
rehace") y el `-env all` (sin él no ve `ASM_DIRS`, hace autodescubrimiento y vuelve
a meter `spades_meta` y `trinity_rescate`).

**Trinity no añade, desplaza.** CD-HIT elige el contig más largo de cada grupo
como representante:

| ensamblador | antes | ahora | ≥2.5 kb | ≥5 kb |
|---|---:|---:|---:|---:|
| trinity | — | **1,964** | 157 | **28** |
| rnaviral | 1,892 | 1,882 | 88 | 22 |
| rna | 1,519 | 1,438 | 160 | 8 |
| megahit | 827 | **711** | 81 | 16 |
| metaviral | 12 | 9 | 3 | 2 |
| meta | 417 | — (excluido) | | |

Trinity ensambló versiones más largas de 116 contigs que antes representaba
MEGAHIT y de 81 que representaba `rna`. **Lidera en ≥5 kb**, por delante de
`rnaviral` y MEGAHIT. Esto matiza la conclusión de agosto sobre ensambladores:
MEGAHIT era el mejor descubridor *de los cinco disponibles entonces*.

`metaviral` merece mención aparte: 9 contigs en toda la unión y 2 de ellos ≥5 kb
(22% de acierto frente al 2.2% de MEGAHIT). Aporta poco volumen y mucho valor.

### Trampa: dos arrays de derep solapados

Se enviaron por error dos arrays (329577 y 330925) que coincidían en los índices
7, 14 y 18 — dos trabajos escribiendo el mismo `$SCRATCH/union/<sample>`. El
directorio de trabajo intermedio es único por `LSB_JOBID`, pero el movimiento
final al destino colisiona. Se mataron los dos y se rehizo entero con `%8`.
Verificación posterior: 22 directorios, sin cabeceras duplicadas, todos con marca
de tiempo dentro del mismo minuto.

---

## Paso 4 — geNomad

### El coste real: `--splits`, no el número de contigs

| corrida | entrada | splits | memoria | tiempo |
|---|---|---:|---:|---:|
| smoke M36_S17 | 16 contigs | 8 | 3,569 MB | **25,735 s (7 h 09 m)** |
| control | 20 genomas | 1 | 17,548 MB | 413 s |
| **corrida final** | **6,004 contigs** | 1 | 17,794 MB | **382 s** |

Siete horas para 16 contigs y seis minutos para 6,004. El `CPU time` del smoke
fue de 982 s sobre 25,735 s de reloj: **3.8% de eficiencia**. Las 8 pasadas de
`--splits` releen 1.4 GB de base desde Lustre cada vez; el coste es la base, no
la consulta.

**Cambios aplicados a `genomad_worker.sh`:**
- línea 29: `SPLITS="${SPLITS:-8}"` → `"${SPLITS:-1}"`
- línea 41: añadido `${GENOMAD_EXTRA:-}` al final del `genomad end-to-end`

El `-M 32000` que traía el worker nunca estuvo medido. Con `--splits 1` el pico
real es ~17.8 GB; con 8 GB muere (`TERM_MEMLIMIT`). **24 GB es el número.**

### Los controles — por qué se quitó el filtro de hallmarks

geNomad aplica `--min-virus-hallmarks-short-seqs 1` por defecto: toda secuencia
**menor de 2,500 bp** necesita al menos un gen *hallmark* viral, tenga el score
que tenga. En la unión, **4,300 de 4,667 contigs (92%) caen bajo ese umbral**.

Tres controles, todos con `checkv_reps.fna` y el genoma de *Aedes* como material:

**1. Control positivo (20 genomas de referencia, 5–50 kb).** 20/20 llamados
virus. Sesgo: el muestreo cogió 20 secuencias consecutivas del fichero, todas
*Priklausovirales*. Sirve para verificar que la maquinaria funciona, no como
medida de sensibilidad. **6 de los 20 tenían `hall=0`** — virus de referencia
completos sin un solo hallmark reconocible.

**2. Experimento de fragmentación (30 genomas diversos + sus 699 fragmentos de
1,500 bp).** La misma secuencia, entera y troceada:

```
                 fragmentos recuperados    genomas con ≥1 fragmento
con filtro           117 / 699  (16.7%)          29 / 30
sin filtro           380 / 699  (54.4%)          30 / 30
                     ───────────────────
                     ganancia 3.25×
```

Fragmentos con ≥1 hallmark viral: **139/699 (19.9%)**. El filtro exige algo que
solo tiene una quinta parte de los fragmentos de virus *de referencia*.

**3. Control negativo (300 fragmentos de 1,500 bp de AaegL5).** 7/300 llamados
virus = **2.3% de falsos positivos**, todos con scores entre 0.704 y 0.754. Es
una cota superior: el control corrió sin calibración, y la calibración elimina
justamente esa franja. Matiz: los 300 salieron del principio del fichero, zona
rica en repeticiones y transposones.

**Parámetros de la corrida final:**

```
--min-virus-hallmarks-short-seqs 0 --enable-score-calibration
```

Justificación en métodos: el filtro por defecto está calibrado para fagos y virus
de ADN; sobre virus ARN fragmentados descarta el 83% de las detecciones posibles,
medido sobre genomas de referencia conocidos.

### Resultado

6,004 contigs → **1,144 virus, 0 plásmidos** (coherente: biblioteca de ARN).

- por ensamblador: rnaviral 444, trinity 348, rna 200, megahit 149, metaviral 3
- por longitud: 79 ≥2,500 bp, 1,065 <2,500 bp
- **1,008 de 1,144 sin taxonomía asignada** (88%)

Todos los scores caen en 0.8–1.0 pese a `min_score 0.7`: la calibración filtró la
franja baja vía `max_fdr 0.1`, que es exactamente donde estaban los 7 falsos
positivos del control negativo.

---

## Paso 5 — CheckV

Dos corridas, y las dos hacen falta:

```
checkv_contigs/all_final/          sobre los 6,004   → red de seguridad independiente
checkv_viral/all_final_virus.fna/  sobre los 1,144   → encadenado tras geNomad
```

| calidad | 6,004 | 1,144 |
|---|---:|---:|
| Complete | 3 | 3 |
| High-quality | 13 | 13 |
| Medium-quality | 21 | 21 |
| Low-quality | 329 | 150 |
| Not-determined | 5,638 | 957 |

**Convergencia total: los 16 mejores de CheckV están los 16 entre los 1,144 de
geNomad**, y 15 con completitud 100%. Las dos herramientas se confirman.

El encadenado solo no basta: si geNomad devuelve 0 virus (como en el smoke de
M36_S17), CheckV no tendría nada que evaluar. Correr también sobre los contigs
completos da la lectura independiente.

`kmer_freq` es la columna que descarta duplicaciones: el contig de Trinity de
15,488 bp de M1_S1 salía Complete al 100% pero con `kmer_freq 1.2` — un 20%
duplicado. 15,488 × 0.8 ≈ 12,400, la longitud real. El representante bueno es el
de ~12.9–13.0 kb.

---

## El anfevirus — y una lección sobre CD-HIT

Doce contigs de ~12.9 kb, completitud 100%, en once muestras. CD-HIT los agrupó
al 95% dando identidades de 95.5–99.5%, con M1_S1 aparentemente divergente al
95.48% y dos contigs suyos separados en clusters propios al 86.73%.

**Todo eso era artefacto.** BLASTn sobre los mismos contigs:

```
metaviral vs megahit    99.98%   aln = 11,290 + 1,872
rnaviral  vs megahit    99.99%   aln =  8,353 + 4,711
rnaviral  vs metaviral  99.98%   aln =  6,628 + 6,464
```

Son 99.98–100% idénticos. CD-HIT falla porque compara palabras en orden colineal
y estos contigs venían en **hebras opuestas**; los dos bloques de alineamiento
son un indel corto alrededor de la posición ~7,722, no una permutación.

**Regla práctica: CD-HIT para agrupar rápido, BLAST o MAFFT para afirmar
identidad.** No usar porcentajes de `.clstr` como medida de divergencia.

### Diversidad real (MAFFT `--adjustdirection`, 12 secuencias)

- alineamiento: 13,562 columnas; **core sin gaps: 12,302**
- sitios variables en el core: **286 (2.32%)**
- distancia por pares: **0.06% – 0.93%**

Menos del 1% entre los extremos: una sola cepa con variación intrapoblacional.

**Estructura:** dos grupos separados de forma consistente (intra ≤0.67%, inter
≥0.74%):

- **A**: M2_S2, M55_S6, M7_S15, M11_S3, M60_S11, M48_S5
- **B**: M1_S1, M62_S18, M20_S25, M21_S29, M22_S30
- **M72_S20**: intermedia (0.51–0.63% a todos)

Un análisis de ventanas de 1,500 bp dio afinidad alternante para M72_S20
(A-A-A-B-A-B-B-A-A), pero **no es evidencia suficiente de recombinación**: con
divergencias del 0.5% hablamos de 7–8 diferencias por ventana. Para afirmarlo
harían falta RDP4, 3SEQ o GARD. Redactar como observación, no como hallazgo.

**Lo que falta cruzar: los metadatos de colecta.** Si A y B corresponden a sitios
o fechas distintas, la estructura se explica sola.

---

## Paso 7 — catálogo de vOTUs

```bash
cd-hit-est -i all_final_virus.fna -o votus -c 0.95 -aS 0.85 -n 10 -d 0 -T 4 -M 8000 -g 1
```

Criterio MIUViG (95% identidad sobre 85% del contig más corto). **1,144 contigs →
748 vOTUs.**

Distribución de prevalencia por ensamblaje:

```
en  1 muestra:  606 vOTUs  (81%)
en  2:  71    en 3: 32    en 4: 15    en 5: 7
en 10+: 12 vOTUs
```

Cola larga clásica de viroma: núcleo pequeño compartido, mayoría de raros.

**M1_S1 no es excepcionalmente diverso.** Tiene 214 vOTUs de 1,043 contigs
(20.5%), pero la fracción viral es 18–20% en todo el lote. Lo que tiene es el
doble de contigs que la siguiente muestra.

---

## Paso 8 — mapeo de vuelta (⚠ resultado contaminado)

Índice bowtie2 sobre las 748 representativas; array de 22 con
`bowtie_map_worker.sh`. 22 BAM, **67 GB** (el worker guarda también las no
alineadas). Cuentas con `samtools idxstats`.

**Los números no son creíbles.** M57_S8 da 104 M de lecturas mapeadas a 1.2 Mb de
vOTUs cuando el 93.76% de su librería mapea a *Aedes*. Y las 748 vOTUs aparecen
en las 22 muestras, lo cual es umbral de detección, no biología.

Las 15 vOTUs más abundantes tienen coberturas de ensamblaje `cov_5569`,
`cov_7318`, `cov_13092`, `cov_3912`, `cov_5284`... **es rRNA que geNomad
clasificó como viral y que se arrastró hasta el catálogo.** La primera son 172 M
de lecturas en 1,684 bp.

**Pendiente:** `samtools coverage` en vez de `idxstats` (amplitud, no volumen), y
presencia = ≥70% del genoma cubierto a ≥1×. Y sacar el rRNA del catálogo: la
firma es cobertura de ensamblaje >1,000× sin hit viral canónico.

### Trampa: `mapping_stats/` se sobrescribió

`bowtie_map_worker.sh` escribe siempre en `$RESULTS_DIR/mapping_stats`,
independientemente del conjunto. El mapeo a vOTUs pisó los ficheros del mapeo al
huésped. Recuperados regenerando `flagstat` desde `$SCRATCH/mapped_bam/MERIDA/`
(los BAM seguían ahí, no en `mapped/`).

**Parche recomendado (sin aplicar aún):** `statsdir="$RESULTS_DIR/mapping_stats_$result_from"`.

---

## Paso 9 — DIAMOND

### Las dos bases

| | tamaño | secuencias | tiempo (748 vOTUs) |
|---|---:|---:|---:|
| `nr.dmnd` | 350 GB | — | 59 min |
| `rvdb.dmnd` (U-RVDB v32 unique) | 357 MB | 783,103 | 8 s / 164 s (`--very-sensitive`) |

**Dos problemas con `nr`:**
1. Construida **sin taxonomía** (faltaron `--taxonmap`, `--taxonnodes`,
   `--taxonnames`). Las columnas `staxids`/`sscinames` salen vacías; hay que leer
   el organismo de `stitle` entre corchetes.
2. El `nr.gz` es de **febrero de 2024** — el script lo encontró y no descargó.
   Para virus de mosquito, dos años y medio importan.

RVDB v32.0 es de este año y específica de virus. **Es la base principal, no la
secundaria.** URL correcta (ojo al guion bajo):
`https://rvdb-prot.pasteur.fr/files/U-RVDBv32.0-prot_unique.fasta.xz`

### Parámetros

`--very-sensitive --max-target-seqs 25 --evalue 1e-3` sube de **253 a 333** vOTUs
con hit (+32%) por 164 s. El modo por defecto está pensado para homología
cercana; el objetivo aquí es lo contrario.

`--outfmt` debe incluir **`qcovhsp` y `scovhsp`**. Sin ellas no se distingue un
hit que cubre el contig entero de uno sobre 40 aminoácidos.

### El filtrado — advertencia metodológica

**Estos filtros no son de ningún programa: son `awk` sobre la salida de DIAMOND.**
Hay que declararlos como criterio propio. Los pasos:

1. Excluir descripciones con `reverse transcriptase|Gag-Pol|retrotranspos|integrase`
2. Exigir proteína viral canónica: `RdRp|RNA-dependent RNA polymerase|capsid|nucleoprotein|glycoprotein|polymerase|nucleocapsid|coat protein`
3. `length ≥ 100` aa y `slen ≥ 150` aa
4. `pident ≥ 90`
5. `qcovhsp ≥ 50` **O** `scovhsp ≥ 50`

**Puntos frágiles:** los pasos 1–2 son búsqueda por palabras clave sobre texto
libre. Una proteína viral llamada "structural protein" u "ORF1" no se captura; un
retrotransposón anotado como "polyprotein" sí pasa.

**Alternativa citable, pendiente:** CAT/BAT (taxonomía por consenso de genes), o
reconstruir `nr.dmnd` con taxonomía y usar el LCA de DIAMOND (`--outfmt 102`).

Por qué el paso 5 es un **O** y no un **Y**: `qcovhsp` penaliza los genomas
largos. El anfevirus de 13 kb tiene `qcov` 37–47% porque la RdRP de referencia
cubre solo una parte del contig — el criterio inicial descartaba el mejor genoma
del catálogo. Con `scovhsp` entra con 90–100%.

**Descartes documentados:**
- *Human papillomavirus* (5 vOTUs): `qcov` 7–10%, bitscore 57–65. Homología
  espuria de la L1.
- *Hanko totivirus 4* (4 vOTUs): `slen = 46` aa. Cubrir entero un péptido de 46
  aminoácidos no identifica nada. Es lo que motivó el filtro `slen ≥ 150`.

### Catálogo firme — 20 vOTUs, 6 especies

| especie | muestras (por ensamblaje) | evidencia |
|---|---:|---|
| **Verdadero virus** | **19 / 22** | RdRP + cápside, 99.8–100% id |
| **Aedes anphevirus** | **14 / 22** | RdRP, 99.3–99.4%, `scov` 90–100% |
| *Aedes aegypti totivirus* | 3 / 22 | RdRP + cápside |
| *Netjeret virus* | 2 / 22 | RdRP + cápside |
| *Phasivirus phasiense* | 1 / 22 | RdRP + glicoproteína ×2 + nucleocápside (4 segmentos) |
| *San Gabriel mononegavirus* | 1 / 22 | nucleoproteína + glicoproteína |

Recuperar RdRP **y** proteínas estructurales del mismo virus es la confirmación
más fuerte posible sin aislamiento. *Phasivirus* con cuatro genes a 96–99% en
M56_S7 está prácticamente completo — es un bunyaviral segmentado, así que los
cuatro contigs son sus segmentos.

**El virus dominante es *Verdadero virus*, no el anfevirus.** La vOTU más
prevalente del catálogo (`M20_S25__megahit__k119_2619`, 19 muestras) es su RdRP
al 99.8%.

Nota: esta prevalencia es "en cuántas muestras se ensambló". La prevalencia real
la da el mapeo, pendiente de limpiar.

### Contexto que no es viral

Entre los mejores hits de `nr` aparecen **Wuchereria bancrofti (33), Brugia timori
(11) y Vibrio cholerae (12)** — filarias transmitidas por mosquito y un patógeno
bacteriano. Es la caja "hospedadores y patógenos" del diagrama del pipeline, sin
explotar. Justifica Kraken2 + Bracken.

---

## Contabilidad de lecturas — el hallazgo transversal

Porcentaje de mapeo al genoma de *Aedes*, frente a contigs ensamblados:

```
M60_S11   98.85% huésped   45.7M lecturas →   91 contigs   (1.05M no-huésped)
M79_S28   98.29%                          →  106
M11_S3    97.61%                          →  106
M48_S5    96.83%                          →   94
─────────────────────────────────────────────────────────
M72_S20   61.16%                          →  158          (17.9M no-huésped)
M74_S21   40.53%                          →  122          (17.0M)
M36_S17   40.44%           18.0M          →   25          (21.5M)
```

**La profundidad no predice nada.** M62_S18 con 8.1 M de lecturas produce los
mismos contigs que M60_S11 con 45.7 M. Y las tres muestras con más lecturas
no-huésped del lote están entre las que menos ensamblan: M36_S17 tiene **el doble
de material no-huésped que M1_S1 y produce 25 contigs frente a 1,043**.

Razón formal: el grafo de De Bruijn no crece con el número de lecturas. Copias de
la misma molécula generan los mismos k-meros y no añaden nodos, solo suben el
contador. El tamaño del grafo lo fija la diversidad, no la profundidad.

**Dos regímenes distintos:**
- **alto huésped (92–99%)**: pocas lecturas útiles, pero diversas
- **bajo huésped (40–61%)**: muchísimas lecturas no-huésped, casi todas rRNA

Ninguno se arregla secuenciando más. El primero pide depleción de huésped en el
protocolo; el segundo, depleción de rRNA.

---

## Rama de rRNA — bloqueada

Tres intentos fallidos con SortMeRNA:

1. `exit 127` — `build_sortmerna_db.sh` no estaba en la farm (solo en el paquete
   local) y no existía entorno conda
2. `exit 132` (SIGILL) — binario de bioconda con instrucciones no soportadas por
   el nodo. La farm mezcla Intel_Platinum (102 nodos) y EPYC7713 (45). Paliado
   con `-R "select[... model==Intel_Platinum]"`
3. `exit 1` tras 590 s — la descarga funcionó (236 MB), falla el indexado. Sin
   diagnosticar

Alternativas si se retoma: `sortmerna` 6.0.2 o 7.0.0 de conda-forge (revisar que
los flags `--ref/--reads/--workdir/--idx-dir/--kvdb` no hayan cambiado), o mapear
con bowtie2 contra secuencias de rRNA de *Aedes*, que usa infraestructura ya
probada.

**No bloquea los pasos 6–10.** Y conviene esperar al `samtools coverage` del paso
8 para saber cuánto rRNA hay realmente en el catálogo antes de invertir más.

---

## Entorno — cosas que costaron tiempo

- **`seqkit` no está en ningún entorno conda.** Hay módulo: `seqkit/2.8.2--h9ee0642_0`.
  El paso 4 del PLAN está escrito con `seqkit grep`. Recomendado:
  `conda install -n cdhit_4.8.1 -c bioconda seqkit`
- **`blastn` tampoco**; módulo `blast/2.15.0--pl5321h6f7f691_1`, y existe un
  entorno `blast_2.17.0`
- **`LC_ALL=C` es obligatorio** antes de `sort` + `join`, o `join` falla en
  silencio con "input is not in sorted order"
- **`bmod` no sustituye un `-R` existente** sin borrarlo antes (`bmod -Rn`).
  Cambió los tasks pero mantuvo `select[mem>32000]`
- **Cuota del grupo `team222`: 66.85 T de 75 T.** `nr.dmnd` (350 GB) + `nr.gz`
  (187 GB) + BAM de vOTUs (67 GB) son tuyos. Los BAM se pueden borrar tras
  extraer las cuentas

---

## Lo siguiente

1. **Limpiar el rRNA del catálogo** y rehacer el paso 8 con `samtools coverage`
2. **Paso 10**: tabla final uniendo `votus_master.tsv` + abundancia + DIAMOND
3. **Kraken2 + Bracken** — la rama de patógenos, sin empezar
4. **Script reproducible** (`build_votu_catalog.sh`): hoy todo el análisis son
   comandos sueltos en el historial
5. **BLASTn contra `nt_core`** sobre las 20 vOTUs firmes: diría si son la misma
   cepa o parientes
6. Actualizar el diagrama del pipeline — no incluye geNomad ni el catálogo de
   vOTUs, que acabaron siendo el núcleo
