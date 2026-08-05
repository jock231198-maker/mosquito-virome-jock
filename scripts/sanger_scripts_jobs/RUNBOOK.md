# RUNBOOK — cómo correr el pipeline en la farm

Todo se lanza desde un **nodo de login** (`farm22-head1`). Los scripts no
hacen cómputo: solo lo mandan a LSF con `bsub`.

---

## 0. Instalación (una sola vez)

```bash
# Carpeta de scripts
mkdir -p ~/jkvirome/jm79/mosquito-virome-jock/sanger_scripts_jobs
cd ~/jkvirome/jm79/mosquito-virome-jock/sanger_scripts_jobs

# ...copia aquí los .sh...

chmod +x *.sh
dos2unix *.sh 2>/dev/null   # por si los editaste en Windows/Mac
```

`dos2unix` no es paranoia: un `\r` al final del shebang da el error
`bad interpreter: /bin/bash^M`, que es críptico y te hace perder media hora.

### Crear los entornos conda

Ya tienes `fastqc_0.12.1` y `multiqc_3.15`. Los demás:

```bash
module load conda

conda create -y -n trimmomatic_0.39 -c bioconda -c conda-forge trimmomatic=0.39
conda create -y -n bowtie2_2.5.4    -c bioconda -c conda-forge bowtie2=2.5.4 samtools=1.21
conda create -y -n spades_4.0.0     -c bioconda -c conda-forge spades=4.0.0
conda create -y -n quast_5.2.0      -c bioconda -c conda-forge quast=5.2.0
conda create -y -n star_2.7.11b     -c bioconda -c conda-forge star=2.7.11b
conda create -y -n checkv_1.0.3     -c bioconda -c conda-forge checkv=1.0.3
conda create -y -n genomad_1.11.0   -c bioconda -c conda-forge genomad=1.11.0
conda create -y -n diamond_2.1.10   -c bioconda -c conda-forge diamond=2.1.10
conda create -y -n cobra_1.2.3      -c bioconda -c conda-forge cobra-meta=1.2.3
```

Nota: `bowtie2_2.5.4` incluye **samtools** a propósito — varios workers lo
necesitan y así no duplicas entornos.

Después edita `config.sh` con los nombres reales:

```bash
ENV_TRIMMOMATIC="${ENV_TRIMMOMATIC:-trimmomatic_0.39}"
ENV_BOWTIE2="${ENV_BOWTIE2:-bowtie2_2.5.4}"
# ...etc
```

### Crear el árbol de directorios

```bash
source ~/jkvirome/jm79/mosquito-virome-jock/sanger_scripts_jobs/config.sh
make_dirs
echo "$DATA_DIR"; ls "$DATA_DIR" | head
```

---

## 1. Prueba de humo (HAZLA ANTES DE NADA)

Antes de lanzar 200 trabajos, comprueba que conda funciona **dentro** de un job.
Los nodos de ejecución no tienen tu `.bashrc` ni tu sesión interactiva.

```bash
cd ~/jkvirome/jm79/mosquito-virome-jock/sanger_scripts_jobs

bsub -o /tmp/smoke.%J.log -e /tmp/smoke.%J.err -q normal -M 1000 \
     -R "select[mem>1000] rusage[mem=1000]" \
     'module load conda && conda activate fastqc_0.12.1 && fastqc --version && which fastqc'

bjobs -w          # espera a que termine
cat /tmp/smoke.*.log
```

Debe imprimir la versión de FastQC. Si falla con `module: command not found`,
`config.sh` ya lo cubre (sourcea `/etc/profile.d/modules.sh`), pero conviene
saberlo.

Segunda prueba: un solo elemento del array, simulando LSF a mano en el login node.

```bash
source ./config.sh
find "$DATA_DIR" -type f -name "*.fastq.gz" | sort > "$SCRATCH/filelist.txt"
wc -l < "$SCRATCH/filelist.txt"

LSB_JOBINDEX=1 ./fastqc_worker.sh "$SCRATCH/filelist.txt"
```

Si esto funciona, el array funcionará.

---

## 2. Pipeline paso a paso

```bash
cd ~/jkvirome/jm79/mosquito-virome-jock/sanger_scripts_jobs
source ./config.sh
S="$PWD"      # atajo: ruta a los scripts
```

### Paso 1 — FastQC sobre los crudos

```bash
find "$DATA_DIR" -type f -name "*.fastq.gz" | sort > "$SCRATCH/filelist.txt"
N=$(wc -l < "$SCRATCH/filelist.txt")
echo "$N archivos"

bsub -J "fastqc[1-$N]%10" \
     -o "$LOGS_DIR/fastqc.%J.%I.log" -e "$LOGS_DIR/fastqc.%J.%I.err" \
     -q normal -n 4 -M 2000 \
     -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" \
     "$S/fastqc_worker.sh $SCRATCH/filelist.txt"
```

### Paso 2 — MultiQC (espera al array anterior)

`-w "done(fastqc)"` encadena: no arranca hasta que **todos** los elementos
del array `fastqc` terminen bien.

```bash
bsub -J multiqc -w "done(fastqc)" \
     -o "$LOGS_DIR/multiqc.%J.log" -e "$LOGS_DIR/multiqc.%J.err" \
     -q normal -M 4000 -R "select[mem>4000] rusage[mem=4000]" \
     "$S/multi_qc.sh"
```

Mira el HTML antes de seguir: te dice qué adaptadores tienes y si hace falta
ajustar los parámetros de Trimmomatic.

### Paso 3 — Trimming

```bash
ls "$DATA_DIR"/*_R1_001.fastq.gz | xargs -n1 basename | cut -d_ -f1-2 \
    | sort -u > "$SCRATCH/samples.txt"
N=$(wc -l < "$SCRATCH/samples.txt")

bsub -J "trim[1-$N]%10" \
     -o "$LOGS_DIR/trim.%J.%I.log" -e "$LOGS_DIR/trim.%J.%I.err" \
     -q normal -n 4 -M 20000 \
     -R "select[mem>20000] rusage[mem=20000] span[hosts=1]" \
     "$S/trimmomatic_worker.sh $SCRATCH/samples.txt $DATA_DIR"
```

Comprueba que `cut -d_ -f1-2` te da los nombres correctos **antes** de lanzar:

```bash
head -3 "$SCRATCH/samples.txt"
```

### Paso 4 — Mapeo al huésped

```bash
bsub -J "btmap[1-$N]%10" -w "done(trim)" \
     -o "$LOGS_DIR/btmap.%J.%I.log" -e "$LOGS_DIR/btmap.%J.%I.err" \
     -q normal -n 4 -M 8000 \
     -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
     "$S/bowtie_map_worker.sh $SCRATCH/samples.txt $SCRATCH/trimmed/$RESULT_FROM"
```

### Paso 5 — Extraer no-mapeados

Este array se alimenta de una lista que **todavía no existe** al momento de
enviar. Por eso se lanza en dos tiempos: primero un job que genera la lista,
luego el array. Lo más simple es esperar y lanzarlo a mano:

```bash
# cuando btmap haya terminado:
find "$SCRATCH/mapped" -type f -name "*.bam" | sort > "$SCRATCH/bam_list.txt"
N2=$(wc -l < "$SCRATCH/bam_list.txt")

bsub -J "unmap[1-$N2]%10" \
     -o "$LOGS_DIR/unmap.%J.%I.log" -e "$LOGS_DIR/unmap.%J.%I.err" \
     -q normal -n 4 -M 8000 \
     -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
     "$S/extract_unmapped_worker.sh $SCRATCH/bam_list.txt"
```

### Paso 6 en adelante

Ensamblaje (SPAdes) → QUAST → geNomad / CheckV / DIAMOND. Mismo patrón:
generar la lista, contar, `bsub` con array.

---

## 3. Vigilar los trabajos

```bash
bjobs -w              # en curso
bjobs -p              # pendientes + POR QUÉ están pendientes
bjobs -A              # resumen del array (cuántos DONE / EXIT)
bhist -l <jobid>      # histórico, incluye memoria máxima usada
bkill <jobid>         # matar
bkill "<jobid>[5]"    # matar solo el elemento 5
```

Cuando un array termina, revisa que no haya fallos:

```bash
grep -l "Exited" "$LOGS_DIR"/fastqc.*.log | wc -l
```

**Truco de memoria**: `bhist -l <jobid> | grep -i "max mem"` te dice cuánta RAM
usó de verdad. Si pediste 32 GB y usó 6, baja el `-M` y tus trabajos entrarán
antes en cola.

---

## 4. Limpieza

Cada vez que termines una etapa:

```bash
"$S/cleanup.sh"                    # mira qué se puede borrar
"$S/cleanup.sh" --force trimmed    # borra los unpaired
lfs quota -h "$SCRATCH"            # comprueba la cuota
```

---

## 5. Guardar lo importante

`results/` está en **Lustre scratch, sin backup**. Cuando tengas resultados
definitivos, cópialos a un sitio con backup:

```bash
rsync -av "$RESULTS_DIR/multiqc/" /nfs/teamYY/results/MERIDA/multiqc/
rsync -av "$RESULTS_DIR/genomad/" /nfs/teamYY/results/MERIDA/genomad/
```

Los nodos de ejecución **no ven `/warehouse`**, así que estas copias se hacen
desde el nodo de login, nunca dentro de un job.

---

## Errores frecuentes

| Síntoma | Causa |
|---|---|
| `bad interpreter: /bin/bash^M` | Finales de línea Windows → `dos2unix *.sh` |
| `conda: command not found` en el job | Falta `module load conda` (lo hace `config.sh`) |
| `TERM_MEMLIMIT` en el log | Pediste poco `-M`. Mira `bhist -l` y sube |
| El array procesa el fichero equivocado | Cambió el orden de `find` → falta el `sort` |
| Job pendiente eternamente | `-M` demasiado alto, o `%N` alto y no hay hueco. `bjobs -p` lo explica |
| `No such file or directory` en el worker | Lanzaste `bsub` con ruta relativa desde otro `cwd`. Usa `$S/script.sh` |
