echo "$0 $@" >> history.txt
source settings.sh
set -e
mkdir -p data_metrics

for PROBLEM in ${PROBLEMS[@]}
do
  EPSILON=${EPSILON_MAP[$PROBLEM]}

  for ALGORITHM in ${ALGORITHMS[@]}
  do
    for SEED in ${SEEDS[@]}
    do
      if [ -z "$USEPBS" ]
      then
        for SAMPLE in ${SAMPLES[@]}
        do
          echo "Processing ${ALGORITHM}_${PROBLEM}_S${SEED}_P${SAMPLE}.metrics"
          ${TEST} rm -f data_metrics/${ALGORITHM}_${PROBLEM}_S${SEED}_P${SAMPLE}.metrics
          ${TEST} java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.ResultFileEvaluator -b ${PROBLEM} -e ${EPSILON} -i data_raw/${ALGORITHM}_${PROBLEM}_S${SEED}_P${SAMPLE}.data -r data_ref/${PROBLEM}.ref -o data_metrics/${ALGORITHM}_${PROBLEM}_S${SEED}_P${SAMPLE}.metrics
        done
      else
        if [ ! -t 1 ] && [ -z "${TEST}" ]
        then
          >&2 echo "Error: Use -t when piping output in USEPBS mode"
          exit -1
        fi

        NAME=EVAL_${ALGORITHM}_${PROBLEM}_S${SEED}
        SCRIPT="\
#!/bin/sh                
#SBATCH --job-name=${NAME}\n\
#SBATCH --nodes=1\n\
#SBATCH -t ${WALLTIME}\n\
#SBATCH -o data_mpi/${NAME}.out\n\
#SBATCH -e data_mpi/${NAME}.err\n\
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=16
cd \$SLURM_SUBMIT_DIR\n\
for SAMPLE in ${SAMPLES[@]}\n\
do\n\
  rm -f data_metrics/${ALGORITHM}_${PROBLEM}_S${SEED}_P\${SAMPLE}.metrics\n\
  java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.ResultFileEvaluator -b ${PROBLEM} -e ${EPSILON} -i data_raw/${ALGORITHM}_${PROBLEM}_S${SEED}_P\${SAMPLE}.data -r data_ref/${PROBLEM}.ref -o data_metrics/${ALGORITHM}_${PROBLEM}_S${SEED}_P\${SAMPLE}.metrics\n\
done"

        echo -e "$SCRIPT" > data_mpi/${NAME}.pbs
        ${TEST} sbatch data_mpi/${NAME}.pbs
      fi
    done
  done
done

