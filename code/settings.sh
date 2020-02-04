ALGORITHMS=( Borg NSGAII NSGAIII RVEA MOEAD )
NSAMPLES=100  
NSEEDS=25
PROBLEMS=( portfolio )

USEPBS=true
SAMPLES=( $(seq 1 ${NSAMPLES}) )
SEEDS=( $(seq 1 ${NSEEDS}) )
#JAVA_ARGS="-cp \"$(echo lib/*.jar | tr ' ' ':'):.\" -Xmx512m"
JAVA_ARGS="-cp MOEAFramework-2.12-Demo.jar"
WALLTIME=10:00:00

while [[ $# > 0 ]]
do
  key="$1"

  case $key in
    -s|--seed)
      if [[ $2 =~ .*\-.* ]]
      then
        SEEDS=$(eval echo {${2//-/..}})
      else
        SEEDS=( $2 )
      fi
      shift
      ;;
    -a|--algorithm)
      ALGORITHMS=( $2 )
      shift
      ;;
    -b|--problem)
      PROBLEMS=( $2 )
      shift
      ;;
    -t|--test)
      TEST=echo
      ;;
    *)
      echo "Unknown option $1"
      ;;
  esac
  shift
done

declare -A EPSILON_MAP
EPSILON_MAP["portfolio"]=0.05,0.1,0.05,0.1


declare -A MAXNFE_MAP

MAXNFE_MAP["portfolio"]=200000
