                                                                   #####PARA MANDAR LOS TRABAJOS AL SERVER######





#check_inputs
bsub -J chkcat -o "$LOGS_DIR/chkcat.%J.log" -e "$LOGS_DIR/chkcat.%J.err" \
     -q normal -n 8 -M 4000 -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
     "$PWD/check_inputs.sh --vs-raw --expect 23 $SCRATCH/cat_fastq"

#fastqc_worker
bsub -J "fastqc[1-$N]%20" \
     -o "$LOGS_DIR/fastqc.%J.%I.log" -e "$LOGS_DIR/fastqc.%J.%I.err" \
     -q normal -n 1 -M 1000 \
     -R "select[mem>1000] rusage[mem=1000] span[hosts=1]" \
     "$PWD/fastqc_worker.sh $SCRATCH/filelist.txt"