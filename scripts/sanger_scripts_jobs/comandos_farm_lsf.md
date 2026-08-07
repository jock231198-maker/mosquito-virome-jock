# Comandos LSF del farm (Sanger) — chuleta

El farm usa **LSF** (Platform LSF). Todos los comandos empiezan por `b`.

---

## 1. Enviar trabajos (`bsub`)

| Comando | Qué hace |
|---|---|
| `bsub script.sh` | Envía un trabajo a la cola por defecto (`normal`) |
| `bsub -o out.%J -e err.%J script.sh` | Ficheros de salida/error (`%J` = jobID, `%I` = índice del array) |
| `bsub -J nombre ...` | Da nombre al trabajo |
| `bsub -q long ...` | Elige cola (ver tabla de colas abajo) |
| `bsub -n 4 -R "span[hosts=1]" ...` | 4 CPUs, todas en la misma máquina |
| `bsub -M4000 -R "select[mem>4000] rusage[mem=4000]" ...` | Reserva 4 GB de RAM (límite y `rusage` deben coincidir) |
| `bsub -R "select[tmp>30000]" ...` | Nodos con >30 GB en `/tmp` |
| `bsub -J "job[1-100]%10" -o out.%J.%I script.sh` | Job array: 100 tareas, máximo 10 simultáneas |
| `bsub -G <grupo> ...` | Imputa el trabajo a tu grupo LSF (`show_my_lsf_groups` lista los tuyos) |
| `bsub -w "done(<jobid>)" ...` | Dependencia: espera a que otro trabajo termine bien |
| `bsub -Is -G <grupo> bash` | Sesión **interactiva** en un nodo de cómputo |
| `bsub < script.sh` | Envía un script con directivas `#BSUB` embebidas |
| `bsub -gpu "mode=shared:j_exclusive=no:gmem=6000:num=1" ...` | Pide 1 GPU con 6 GB de memoria |

> Nota: el `-M` va en MB. Si el trabajo excede la memoria pedida, LSF lo mata (`TERM_MEMLIMIT`).

---

## 2. Monitorizar trabajos

| Comando | Qué hace |
|---|---|
| `bjobs` | Lista tus trabajos: `PEND`, `RUN`, `DONE`, `EXIT` |
| `bjobs -w` | Igual pero sin truncar los campos |
| `bjobs -l <jobid>` | Detalle completo: recursos, host de ejecución, motivo de espera |
| `bjobs -p` | Solo pendientes, con el motivo |
| `bjobs -p1 <jobid>` | Motivo **principal** por el que un trabajo está pendiente |
| `bjobs -d` | Trabajos terminados (solo la última hora) |
| `bjobs -u all` | Trabajos de todo el mundo |
| `bjobs -A <jobid>` | Resumen de un job array |
| `bhist -l <jobid>` | Historial de trabajos más antiguos |
| `bacct <jobid>` | Contabilidad final: CPU, memoria máxima y tiempo usados |
| `bqueues` | Estado de las colas y número de trabajos |
| `bqueues -l <cola>` | Detalle de una cola |
| `bqueues -r` | Tu prioridad (fairshare) |
| `bhosts -w` | Estado de los nodos |
| `bhosts -l <host>` | Detalle de un nodo (memoria libre, carga...) |
| `bhosts -gpu` | Nodos con GPU y su ocupación |
| `bhosts -s` | Recursos compartidos / tokens de licencia |
| `lsload` / `lshosts` | Carga y modelos de máquina disponibles |

---

## 3. Controlar trabajos

| Comando | Qué hace |
|---|---|
| `bkill <jobid>` | Mata un trabajo |
| `bkill 0` | Mata **todos** tus trabajos |
| `bkill "<jobid>[3]"` | Mata un elemento concreto de un array |
| `bstop <jobid>` | Suspende |
| `bresume <jobid>` | Reanuda |
| `bmod -M8000 <jobid>` | Modifica recursos de un trabajo aún pendiente |
| `btop <jobid>` | Sube el trabajo al principio de tu cola |
| `bbot <jobid>` | Lo manda al final |
| `brequeue <jobid>` | Vuelve a encolar un trabajo |

---

## 4. Colas más habituales

| Cola | Uso típico |
|---|---|
| `small` | Trabajos muy cortos y ligeros |
| `normal` | Por defecto (~12 h) |
| `long` | Hasta ~48 h |
| `week` | Hasta 7 días |
| `basement` | Muy largos, prioridad baja |
| `hugemem` / `teramem` | Mucha RAM |
| `parallel` | Trabajos multinodo |
| `transfer` | Movimiento de datos |
| `gpu-normal` / `gpu-huge` / `gpu-basement` | Trabajos con GPU |
| `yesterday` | Prioridad alta, uso restringido |

`bqueues` te dice si el farm está abierto (`Open:Active`) y los límites de cada cola.

---

## 5. Plantilla de script con directivas `#BSUB`

```bash
#!/bin/bash
#BSUB -J mi_job
#BSUB -o logs/mi_job-%J.out
#BSUB -e logs/mi_job-%J.err
#BSUB -q normal
#BSUB -n 4
#BSUB -M 20000
#BSUB -R "select[mem>20000] rusage[mem=20000] span[hosts=1]"

set -euo pipefail
mi_programa --threads 4 --input datos.fastq.gz --output resultados/
```

Se envía con:

```bash
bsub < mi_script.sh
```

### Plantilla de job array

```bash
#!/bin/bash
#BSUB -J "muestras[1-96]%20"
#BSUB -o logs/muestra-%J.%I.out
#BSUB -e logs/muestra-%J.%I.err
#BSUB -q long
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000]"

SAMPLE=$(sed -n "${LSB_JOBINDEX}p" lista_muestras.txt)
mi_programa --sample "$SAMPLE"
```

Variables útiles dentro del trabajo: `$LSB_JOBID`, `$LSB_JOBINDEX`, `$LSB_JOBNAME`, `$LSB_HOSTS`.

---

## 6. Diagnóstico rápido

| Síntoma | Qué mirar |
|---|---|
| El trabajo no arranca | `bjobs -p1 <jobid>` — suele ser memoria/CPU pedidos imposibles o cola llena |
| Se murió sin explicación | Fichero `-e`, luego `bacct <jobid>` o `bhist -l <jobid>` (busca `TERM_MEMLIMIT` / `TERM_RUNLIMIT`) |
| Saber cuánta memoria usó de verdad | `bacct <jobid>` → `MAX MEM` (útil para ajustar `-M` la próxima vez) |
| Nada se ejecuta en todo el farm | `bqueues` — comprueba que la cola esté `Open:Active` |
| Ver el nodo donde corre | `bjobs -l <jobid>` → `EXEC_HOST`, luego `ssh` a ese nodo para `top` |
