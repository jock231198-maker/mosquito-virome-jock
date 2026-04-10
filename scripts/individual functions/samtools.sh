#!/bin/bash
# Ejemplo de comando para extraer solo las lecturas no mapeadas de un archivo BAM usando samtools
samtools view -b -f 4 mapped.bam > unmapped.bam

# Ejemplo de comando para convertir un archivo BAM a FASTQ usando samtools
samtools fastq unmapped.bam > unmapped.fastq
