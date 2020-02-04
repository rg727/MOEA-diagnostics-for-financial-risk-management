
setwd("C:/Users/Rohini/Box Sync/MOEAFramework_Group/metrics/Cut_Files/")


SEEDS = 10
PARAMS = 10

i=1
j=1

initial_row_of_zeros=numeric(PARAMS)
df = read.delim(paste0("Borg_lake_S",j,"_P",i,".txt"), sep=" ") 
nrow=dim(df)[1]

seed_values = data.frame(matrix(0, ncol = SEEDS, nrow = nrow))
final_param_consolidation = data.frame(matrix(0, ncol = PARAMS, nrow = nrow))

for (i in 1:PARAMS){
  for (j in 1:SEEDS){
    df = read.delim(paste0("Borg_lake_S",j,"_P",i,".txt"), sep=" ") 
    names(df)[1]<-"Hypervolume"
    seed_values[,j] = (df)[1]
  }
  final_param_consolidation[,i]=rowMeans(seed_values)
}


hypervolume_Borg=rbind(final_param_consolidation)
filename='hypervolume_Borg.txt'
write.table(hypervolume_Borg, file = filename, sep = " ", row.names = FALSE, col.names=TRUE, quote=FALSE)

#####################################################################################

rm(list=ls())

setwd("C:/Users/Rohini/Box Sync/MOEAFramework_Group/metrics/Cut_Files/")


SEEDS = 10
PARAMS = 10

i=1
j=1

initial_row_of_zeros=numeric(PARAMS)
df = read.delim(paste0("NSGAII_lake_S",j,"_P",i,".txt"), sep=" ") 
nrow=dim(df)[1]

seed_values = data.frame(matrix(0, ncol = SEEDS, nrow = nrow))
final_param_consolidation = data.frame(matrix(0, ncol = PARAMS, nrow = nrow))

for (i in 1:PARAMS){
  for (j in 1:SEEDS){
    df = read.delim(paste0("NSGAII_lake_S",j,"_P",i,".txt"), sep=" ") 
    names(df)[1]<-"Hypervolume"
    seed_values[,j] = (df)[1]
  }
  final_param_consolidation[,i]=rowMeans(seed_values)
}


hypervolume_NSGAII=rbind(final_param_consolidation)
filename='hypervolume_NSGAII.txt'
write.table(hypervolume_NSGAII, file = filename, sep = " ", row.names = FALSE, col.names=TRUE, quote=FALSE)

