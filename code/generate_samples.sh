echo "$0 $@" >> history.txt
source settings.sh
set -e

for ALGORITHM in ${ALGORITHMS[@]}
do
  if [ "${ALGORITHM}" == "NSGAIII" ] || [ "${ALGORITHM}" == "RVEA" ]
  then
    for PROBLEM in ${PROBLEMS[@]}
    do
      ${TEST} java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.SampleGenerator -m latin -n ${NSAMPLES} -p ${ALGORITHM}_${PROBLEM}_Params.txt -o ${ALGORITHM}_${PROBLEM}_Samples.txt
    done
  else 
    ${TEST} java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.SampleGenerator -m latin -n ${NSAMPLES} -p ${ALGORITHM}_Params.txt -o ${ALGORITHM}_Samples.txt
  fi
done

