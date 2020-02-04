#The purpose of this script is to take the parameter files and the consoliated hypervolume
#and make them into a single text file that can be read by the python file to make control maps
setwd("C:/Users/Rohini/Box Sync/MOEAFramework_Group/metrics/Cut_Files/")

PARAMS=10

#First import the hypervolume file for Borg

df = read.delim("hypervolume_Borg.txt", sep=" ") 

#Concatenate values into one long vector 

z=data.frame(unlist(df))

#import the population values 

df_pop=read.delim("Borg_Samples.txt", sep=" ",header=F) 

x=df_pop[,1]

x_final=data.frame(rep(x,each=dim(df)[1]))

#Make a column for NFE 

NFE=seq(from = 1000, to = 10000, length.out = dim(df)[1])

NFE_final=data.frame(rep(NFE,PARAMS))

final_table=data.frame(x_final,NFE_final,z)

colnames(final_table)=c("x","y","z")

filename='FINAL_TABLE_BORG.txt'
write.table(final_table, file = filename, sep = " ", row.names = FALSE, col.names=TRUE, quote=FALSE)

###################################################################
rm(list=ls())
setwd("C:/Users/Rohini/Box Sync/MOEAFramework_Group/metrics/Cut_Files/")


PARAMS=10

#First import the hypervolume file for NSGAII

df = read.delim("hypervolume_NSGAII.txt", sep=" ") 

#Concatenate values into one long vector 

z=data.frame(unlist(df))

#import the population values 

df_pop=read.delim("NSGAII_Samples.txt", sep=" ",header=F) 

x=df_pop[,1]

x_final=data.frame(rep(x,each=dim(df)[1]))

#Make a column for NFE 

NFE=seq(from = 1000, to = 10000, length.out = dim(df)[1])

NFE_final=data.frame(rep(NFE,PARAMS))

final_table=data.frame(x_final,NFE_final,z)

colnames(final_table)=c("x","y","z")

filename='FINAL_TABLE_NSGAII.txt'
write.table(final_table, file = filename, sep = " ", row.names = FALSE, col.names=TRUE, quote=FALSE)

#############################################
