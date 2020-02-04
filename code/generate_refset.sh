echo "$0 $@" >> history.txt
source settings.sh
set -e
mkdir -p data_ref

for PROBLEM in ${PROBLEMS[@]}
do
  IDS=""
  EPSILON=${EPSILON_MAP[$PROBLEM]}

  for ALGORITHM in ${ALGORITHMS[@]}
  do
    if [ -z "$USEPBS" ]
    then
      for SAMPLE in ${SAMPLES[@]}
      do
        echo "Generating ${ALGORITHM}_${PROBLEM}_P${SAMPLE}.set"
        ${TEST} java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.ResultFileMerger -b ${PROBLEM} -o data_ref/${ALGORITHM}_${PROBLEM}_P${SAMPLE}.set -e ${EPSILON} data_raw/${ALGORITHM}_${PROBLEM}_S*_P${SAMPLE}.data
      done
      ${TEST} java ${JAVA_ARGS} org.moeaframework.util.ReferenceSetMerger -o data_ref/${ALGORITHM}_${PROBLEM}.set -e ${EPSILON} data_ref/${ALGORITHM}_${PROBLEM}_*.set > /dev/null
    else
      if [ ! -t 1 ] && [ -z "${TEST}" ]
      then
        >&2 echo "Error: Use -t when piping output in USEPBS mode"
        exit -1
      fi

      NAME=SET_${ALGORITHM}_${PROBLEM}
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
java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.ResultFileMerger -b ${PROBLEM} -o data_ref/${ALGORITHM}_${PROBLEM}_P\${SAMPLE}.set -e ${EPSILON} data_raw/${ALGORITHM}_${PROBLEM}_S*_P\${SAMPLE}.data\n\
done\n\
java ${JAVA_ARGS} org.moeaframework.util.ReferenceSetMerger -o data_ref/${ALGORITHM}_${PROBLEM}.set -e ${EPSILON} data_ref/${ALGORITHM}_${PROBLEM}_*.set > /dev/null"

      echo -e "$SCRIPT" > data_mpi/${NAME}.pbs

      if [ -z "$TEST" ]
      then
        ID=$(sbatch data_mpi/${NAME}.pbs)
        echo $ID

        # Create array of ids so we can put a job hold on the final script
        ID=$(echo $ID | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}')
        IDS="$IDS:$ID"
      else
        ${TEST} sbatch data_mpi/${NAME}.pbs
      fi
    fi
  done

  if [ -z "$USEPBS" ]
  then
    echo "Generating ${PROBLEM}.ref"
    ${TEST} java ${JAVA_ARGS} org.moeaframework.util.ReferenceSetMerger -o data_ref/${PROBLEM}.ref -e ${EPSILON} data_ref/*_${PROBLEM}_*.set > /dev/null
  else
    NAME=REF_${PROBLEM}
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
java ${JAVA_ARGS} org.moeaframework.util.ReferenceSetMerger -o data_ref/${PROBLEM}.ref -e ${EPSILON} data_ref/*_${PROBLEM}_*.set > /dev/null"

    echo -e "$SCRIPT" > data_mpi/${NAME}.pbs
    ${TEST} sbatch --dependency=afterok$IDS data_mpi/${NAME}.pbs
  fi
done

