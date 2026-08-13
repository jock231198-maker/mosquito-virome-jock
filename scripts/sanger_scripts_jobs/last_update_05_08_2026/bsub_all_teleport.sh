                                                                  #####PARA MANDAR LOS TRABAJOS AL SERVER######


#rooting and directories stablishment 
source ~/jkvirome/jm79/mosquito-virome-jock/scripts/sanger_scripts_jobs/last_update_05_08_2026/config.sh
alias jkv='cd ~/jkvirome/jm79/mosquito-virome-jock/scripts/sanger_scripts_jobs/last_update_05_08_2026 && source ./config.sh'


#check_inputs
bsub -J chkcat -o "$LOGS_DIR/chkcat.%J.log" -e "$LOGS_DIR/chkcat.%J.err" \
     -q normal -n 8 -M 4000 -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
     "$PWD/check_inputs.sh --vs-raw --expect 23 $SCRATCH/cat_fastq"

#fastqc_worker
#make_array 
find "$SCRATCH/cat_fastq" -type f -name "*.fastq.gz" | sort > "$SCRATCH/filelist.txt" ##check that we have to change the find directory 
N=$(wc -l < "$SCRATCH/filelist.txt"); echo "N=$N"   #change the name of the txt to de find folder directory name
#bsub_fastqc
bsub -J "fastqc[1-$N]%20" \
     -o "$LOGS_DIR/fastqc.%J.%I.log" -e "$LOGS_DIR/fastqc.%J.%I.err" \
     -q normal -n 1 -M 1000 \
     -R "select[mem>1000] rusage[mem=1000] span[hosts=1]" \
     "$PWD/fastqc_worker.sh $SCRATCH/filelist.txt"  #stablishment of the array file txt and the output directory ->  Ex.: "$PWD/fastqc_worker.sh $SCRATCH/filelist_nopolyg.txt $RESULTS_DIR/fastqc_nopolyg"

./check_step.sh fastqc        # cuando termine el array

#multi_qc does not need a bsub
./multi_qc.sh "$RESULTS_DIR/<fastqc_directory>" "$RESULTS_DIR/multiqc" multiqc_name_of_the_fastqc_directory 2>&1 | tail -20
ls "$RESULTS_DIR/multiqc"/name_of_the_fastqc_directory.html
./multi_qc.sh "$RESULTS_DIR/fastqc_nopolyg" "$RESULTS_DIR/multiqc" multiqc_nopolyg 2>&1 | tail -20
ls "$RESULTS_DIR/multiqc"/multiqc_nopolyg.html
#polyg_worker
#make_array 
find "$SCRATCH/cat_fastq" -type f -name "*.fastq.gz" | sort > "$SCRATCH/filelist.txt"
N=$(wc -l < "$SCRATCH/filelist.txt"); echo "N=$N"   
#bsub_polygremove
bsub -J "polyg[1-22]%10" \
     -o "$LOGS_DIR/polyg.%J.%I.log" -e "$LOGS_DIR/polyg.%J.%I.err" \
     -q normal -n 4 -M 2000 \
     -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" \
     "$PWD/polyg_worker.sh $SCRATCH/samples.txt $SCRATCH/cat_fastq"

     ./check_step.sh polyg        # cuando termine el array

