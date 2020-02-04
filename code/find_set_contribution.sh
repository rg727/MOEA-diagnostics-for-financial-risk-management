echo "$0 $@" >> history.txt
source settings.sh

for PROBLEM in ${PROBLEMS[@]}
do
  for ALGORITHM in ${ALGORITHMS[@]}
  do
    for SAMPLE in ${SAMPLES[@]}
    do
      ${TEST} java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.SetContribution -e ${EPSILON_MAP[$PROBLEM]} -r ./data_ref/${PROBLEM}.reference  ./data_ref/Borg_${PROBLEM}.set ./data_ref/MOEAD_${PROBLEM}.set ./data_ref/NSGAII_${PROBLEM}.set ./data_ref/NSGAIII_${PROBLEM}.set ./data_ref/RVEA_${PROBLEM}.set > ./data_ref/${PROBLEM}_set_contribution.txt
    done
  done
done

