echo "$0 $@" >> history.txt
source settings.sh

for PROBLEM in ${PROBLEMS[@]}
do
  for ALGORITHM in ${ALGORITHMS[@]}
  do
    for SAMPLE in ${SAMPLES[@]}
    do
     java ${JAVA_ARGS}  org.moeaframework.analysis.sensitivity.SimpleStatistics --mode average --output ./average_metrics/${ALGORITHM}_${PROBLEM}_P${SAMPLE}.average ./data_metrics_new/${ALGORITHM}_${PROBLEM}_S*_P${SAMPLE}.txt

    done
  done
done

