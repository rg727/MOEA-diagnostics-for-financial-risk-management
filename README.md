# MOEA_diagnostics_for_financial_risk_management
Repository for Gupta, Hamilton, Reed, &amp; Characklis, 2020, "Can Modern Multi-Objective Evolutionary Algorithms Discover High-Dimensional Financial Risk Portfolio Tradeoffs for Snow-Dominated Water-Energy Systems?" (in prep)

This repository contains all code needed to run the diagnostics experiments using MOEAFramework and generate data. The raw data are exluded due to size constraints, but the data that is ultimately needed to create figures is provided in the respective folders of ./figures. This experiment was run on THCUBE Cluster at Cornell University. 

Note: Two missing files are:

1) Borg.jar that resides in ./code/lib, which is a Java version of the Borg Multi-Objective Evolutionary Algorithm (http://borgmoea.org/), a licensed software. For the .jar file, please contact: rg727@cornell.edu. Alternatively, the diagnostics can be performed without Borg, using any of the standard algorithms in the MOEAFramework: http://moeaframework.org/javadoc/org/moeaframework/algorithm/StandardAlgorithms.html 

2) HHSamp11122018.txt which can be found at: https://drive.google.com/file/d/1sKbrKztGQgz6AWbVi3yLJCGWBqCBSIfd/view?usp=sharing and goes into ./code (produced by Andrew L. Hamilton from https://github.com/ahamilton144/hamilton-2020-managing-financial-risk-tradeoffs-for-hydropower)

## Contents ##
* Code: Contains all code necessary to perform the diagnostics experiments, create reference sets, and metrics
* Data Processing: Contains all code necessary to process output for figure generation 
* Figures: Contains all code and corresponding data to create figures

## Steps to Run for Generating Data (./code) ## 

1) Pairing The Problem with MOEAFramework: The relevant files are: 
    * portfolio.cpp : c++ implementation of the test case (test case written by Andrew L. Hamilton) 
    * portfolio.java : java file with problem specifications
    * global.properties: sets global properties for MOEAFramework (notably the name of the problem) 
    * Makefile : compiles relevant libraries into a single executable called "portfolio"
    * MOEAFramework-2.12-Demo.jar : MOEAFramework source code combined into a .jar file
    
 -type "make" to create the executable  
 -type "javac -classpath MOEAFramework-2.12-Demo.jar portfolio.java (this will create portfolio.class file that MOEAFramework can read)
 
 2) Problem Specifications: The relevant files are: 
    * settings.sh : Specify algorithms, samples, seeds, NFE, epsilons, etc (Written by Dave Hadka) 
    
-type "./settings.sh" to update problem settings 
 
 3) Algorithm Sampling: The relevant files are: 
    * generate_samples.sh : generates Latin hypercube sample (of size specified in settings.sh) of algorithm parameterization from  parameter files in ./algorithm_params (written by Dave Hadka)
-type "./generate_samples.sh" to generate the text files with samples 
(note that the default paramaterizations can be found in ./default_parameters and the whole experiment can be rerun just using the defaults)

 4) Running the Optimization: The relevant files are: 
    * run.sh : performs optimization of portfolio problem and outputs runtime files for each seed and parameterization in the data_raw folder. The relevant libraries that will be loaded for this are housed in ./lib. The files in data_mpi have information on if each run is sucessful. This is a good first place to look for troubleshooting. (Written by Dave Hadka)
 -type "./run.sh"   
 
 5) Generation the Reference Set: 
    * generate_ref_set.sh : generates a reference set for each algorithm (across all seeds and parameterizations) and across all algorithms to create an overall problem reference set stored in data_ref. (Written by Dave Hadka)
  -type "./generate_ref_set.sh"
  
6)  Determine Reference Set Contribution: 
    Copy MOEAFramework .jar file into data_ref folder  
    -Add a "#" at the end of individual algorithm reference sets  
    - type "java -cp MOEAFramework-2.12-Demo.jar org.moeaframework.analysis.sensitivity.SetContribution -e 0.05,0.1,0.5,0.1 -r portfolio.ref Borg_portfolio.set NSGAII_portfolio.set MOEAD_portfolio.set RVEA_portfolio.set NSGAIII_portfolio.set > portfolio_set_contribution.txt". The result is a text file with percentages of algorithm contribution
    
7) Generation of Metrics: 
   * evaluate.sh : Creates metrics file for each seed and parameterization for each algorithm. (Written by Dave Hadka)
  -type "./evaluate.sh" and metrics files will be created in ./data_metrics 

8) Calculating Overall Problem Hypervolume: The relevant files are: 

   * HypervolumeEval.java (and the corresponding .class file that is created from it (-type "javac -classpath MOEAFramework-2.12-   Demo.jar HypervolumeEval.java) (written by Dave Hadka)
-type "java -cp MOEAFramework-2.12-Demo.jar HypervolumeEval ./data_ref/portfolio.ref >> lake_ref.hypervolume
 
## Steps to Run for Processing Data (./data_processing) ##

1) Cutting Metrics to a Uniform Length:  
    -drag your metrics files to ./data_processing/metrics and change all extensions to .txt (command prompt command: ren *.metrics .*txt)  
    -Use Cutting_Script.R to find the maximum number of rows that are common among all seeds and creates new metric files in the folder,
Cut_Files  

2) More Processing for Figures  
  -Use Seed_Merge.R: to create a text file with the average hypervolume for each
parameterization for each algorithm (i.e. hypervolume_Borg.txt)  
  -Add all algorithm samples.txt files to the folder   
  -Use Make_Final_Table.R to takes the population values from the sample file and the hypervolume values for each parameterization and puts them into a final matrix in a form to be accepted by the control map code.
  -To get data ready for attainment plots, drag the "cut" metrics back into a new directory, data_metrics_new in ./code.  
  -Use ./average_metrics.sh to obtain a set of average metrics across seeds for each algorithm  
  -Concatenate the files (using Borg as an example, but do for all algorithms): "cat Borg_portfolio_*.average >> Borg_Concatenate_Average_Metrics.txt
  
## Steps to Run for Generating Figures ##

All necessary code and data needed for figure generation is provided in each folder of ./figures. Note that for the random seed analysis (Figure 11), the experiment must be re-run, but using only the default parameters, not a LHS. 


    
 
 
