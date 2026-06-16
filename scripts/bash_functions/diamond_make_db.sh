#!/bin/bash

cd /Users/JK/Desktop/diamond # diamond database directory

# Combinar ambos archivos
cat viral.1.protein.faa.gz viral.2.protein.faa.gz > viral_refseq.faa.gz

# Construir base de datos .dmnd
diamond makedb \
    --in viral_refseq.faa \
    --db viral_refseq_proteins \
    --threads 4