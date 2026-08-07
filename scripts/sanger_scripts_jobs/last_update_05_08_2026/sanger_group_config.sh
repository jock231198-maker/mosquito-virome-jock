#Sanger group bsub
bugroup -w 2>/dev/null | grep -w jm79 | head
export LSB_DEFAULT_USERGROUP=team222
#La farm exige declarar el grupo de facturación. El tuyo es team222 — sale en el ls -lh (jm79 team222) y en la cuota (grp team222).


#prueba de humo para comprobar bsubs normalmente
source ./config.sh
bsub -o "$LOGS_DIR/smoke.%J.log" -e "$LOGS_DIR/smoke.%J.err" \
     -q normal -M 1000 -R "select[mem>1000] rusage[mem=1000]" \
     'module load conda && conda activate fastqc_0.12.1 && fastqc --version && which fastqc'