# jkvirome — análisis de transcriptoma en la farm de Sanger

Notas de organización y flujo de trabajo. Guarda este archivo en
`~/jkvirome/README.md`.

## Dónde va cada cosa

La regla de oro: cada filesystem de la farm tiene un propósito distinto.
No mezclarlos evita el 90% de los problemas.

| Sitio | Qué va aquí | Backup |
|-------|-------------|--------|
| `$HOME/jkvirome/` | scripts, notas, entornos conda | Sí |
| `/lustre/scratchXXX/teamYY/jk/jkviromedata/` | TODO el cómputo: crudos de trabajo, intermedios, resultados, logs | **No** |
| `/nfs/teamYY` o `/warehouse/teamYY` | copia maestra de crudos y resultados finales comprimidos | Sí |

Notas importantes:

- Los nodos de **ejecución no pueden ver `/warehouse`** — el dato tiene que
  estar en Lustre cuando corre el trabajo.
- No poner rutas de Lustre en `.bashrc` ni en el `PATH`: hace la sesión
  vulnerable a caídas transitorias. Instalar software en `$HOME` o `/software`.
- Revisar la cuota a menudo: `lfs quota -h /lustre/scratchXXX`.

## Estructura de directorios

    ~/jkvirome/
    
    └─/git/mosquito-virome-jock
      ├── scripts    bash R jupyter py.3.14 scripts
      ├── results    
      ├── docs
      ├── sanger_scripts_jobs    Scripts made for sanger server (cat_samples.sh, fastqc_worker.sh, run_multiqc.sh)
    
    /lustre/scratchXXX/teamYY/jk/jkviromedata/
    ├── raw_fastq/               copia de trabajo de los crudos
    ├── cat_fastq/               salida de cat_samples.sh
    ├── fastqc_results/          salida de FastQC
    ├── multiqc_report/          salida de MultiQC
    └── logs/                    .out / .err de LSF

Mantengo los nombres de subcarpetas iguales a los de la Mac: así el mismo
script corre en ambos sitios cambiando solo la raíz (que pasa por argumento).

## Pipeline (de principio a fin)

Todo desde un **nodo de login**, situado en tu scratch:

    cd /lustre/scratchXXX/teamYY/jk/jkviromedata

    # 0. Traer crudos desde el área del equipo (login node, no /nfs para computar)
    cp /nfs/teamYY/runs/M57/*.fastq.gz raw_fastq/

    # 1. Concatenar lanes  (I/O puro; pocas muestras -> login node)
    ~/jkvirome/scripts/cat_samples.sh raw_fastq cat_fastq

    # 2. FastQC en paralelo con un job array
    find cat_fastq -type f -name "*.fastq.gz" | sort > filelist.txt
    N=$(wc -l < filelist.txt)
    bsub -J "fastqc[1-$N]%10" \
         -o "logs/fastqc.%J.%I.log" -e "logs/fastqc.%J.%I.err" \
         -q normal -n 4 -M 2000 \
         -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" \
         "~/jkvirome/scripts/fastqc_worker.sh filelist.txt fastqc_results"

    # 3. Agregar reportes con MultiQC (espera a que termine el array)
    bsub -o "logs/multiqc.%J.log" -e "logs/multiqc.%J.err" \
         -q normal -M 4000 -R "select[mem>4000] rusage[mem=4000]" \
         "~/jkvirome/scripts/run_multiqc.sh fastqc_results multiqc_report"

    # 4. Guardar resultados y limpiar
    cp multiqc_report/*.html /nfs/teamYY/results/
    # borra intermedios que no necesites para liberar cuota

El entorno conda `Fast_QC` contiene FastQC y MultiQC.

## Consejos de LSF

- Empieza pequeño: el `%10` del array limita a 10 a la vez. Sube solo
  después de medir. Nunca lances miles de golpe.
- Memoria por defecto = 100 MB; reserva más con `-M` + `rusage[mem=...]`.
- Diagnóstico: `bjobs -w` (en curso), `bjobs -p` (por qué pende),
  `bjobs -d` (terminados), `bhist` (históricos), `bkill <id>` (matar).
- Encadenar pasos automáticamente: `-w "done(<jobname>)"` hace que MultiQC
  espere a que termine el array de FastQC sin intervención manual.

## Acceso (Teleport vs SSH)

Te conectas por **Teleport** (terminal web), que ya gestiona el login —
no necesitas un `~/.ssh/config` para entrar. Lanza los trabajos con `bsub`
(corren en segundo plano en la farm), así si se cae la sesión web no pierdes
nada.

Un `~/.ssh/config` solo te hace falta si algún día transfieres datos
directamente con `scp`/`rsync` a través del gateway desde tu portátil.
Plantilla mínima en ese caso:

    Host sanger-gw
        HostName ssh.sanger.ac.uk
        User TU_USUARIO
        LocalForward 2222 farm5-login.internal.sanger.ac.uk:22

    # luego, con la conexión al gateway abierta:
    #   rsync -avz -e 'ssh -p 2222' ./datos localhost:/lustre/.../raw_fastq/
