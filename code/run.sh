echo "$0 $@" >> history.txt
source settings.sh
set -e


for PROBLEM in ${PROBLEMS[@]}
do
  EPSILON=${EPSILON_MAP[$PROBLEM]}
  MAXNFE=${MAXNFE_MAP[$PROBLEM]}
  
  for SAMPLE in ${SAMPLES[@]}
  do
    for ALGORITHM in ${ALGORITHMS[@]}
    do
      if [ "${ALGORITHM}" == "Borg" ]
      then
        PARAMETERS="-e ${EPSILON}"
  # -x \"restartMode=3;probabilityMode=2\""
      elif [ "${ALGORITHM}" == "NSGAIII" ]
      then
        PARAMETERS="-x \"sbx.swap=false\""
      else
        PARAMETERS=""
      fi
  
      if [ "${ALGORITHM}" == "NSGAIII" ] || [ "${ALGORITHM}" == "RVEA" ]
      then
        PARAM_FILE="${ALGORITHM}_${PROBLEM}_Params.txt"
        SAMPLES_FILE="${ALGORITHM}_${PROBLEM}_Samples.txt"
      else
        PARAM_FILE="${ALGORITHM}_Params.txt"
        SAMPLES_FILE="${ALGORITHM}_Samples_.txt"
      fi
  
      for SEED in ${SEEDS[@]}
      do
        if [ -z "$USEPBS" ]
        then
          ${TEST} rm -f data_raw/${ALGORITHM}_${PROBLEM}_S${SEED}_P${SAMPLE}.data
          ${TEST} java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.DetailedEvaluator -p ./algorithm_params/${PARAM_FILE} -i ./${SAMPLES_FILE} -b ${PROBLEM} -a ${ALGORITHM} -s ${SEED} -f 5000 -x maxEvaluations=${MAXNFE} ${PARAMETERS} -o data_raw/${ALGORITHM}_${PROBLEM}_S${SEED}_P${SAMPLE}.data
        else
          if [ ! -t 1 ] && [ -z "${TEST}" ]
          then
            >&2 echo "Error: Use -t when piping output in USEPBS mode"
            exit -1
          fi
  
          NAME=RUN_${ALGORITHM}_${PROBLEM}_S${SEED}_P${SAMPLE}
          SCRIPT="\
#!/bin/sh                
#SBATCH --job-name=${NAME}\n\
#SBATCH --nodes=1\n\
#SBATCH -t ${WALLTIME}\n\
#SBATCH -o data_mpi/${NAME}.out\n\
#SBATCH -e data_mpi/${NAME}.err\n\
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
cd \$SLURM_SUBMIT_DIR\n\
rm -f data_raw/${ALGORITHM}_${PROBLEM}_S${SEED}_P${SAMPLE}.data\n\
java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.DetailedEvaluator -p ./algorithm_params/${PARAM_FILE} -i ./${SAMPLES_FILE} -b ${PROBLEM} -a ${ALGORITHM} -f 5000 -s ${SEED} -x maxEvaluations=${MAXNFE} ${PARAMETERS} -o data_raw/${ALGORITHM}_${PROBLEM}_S${SEED}_P${SAMPLE}.data"
  
          echo -e "$SCRIPT" > data_mpi/${NAME}.pbs
          ${TEST} sbatch data_mpi/${NAME}.pbs
        fi
      done
    done
  done
done

