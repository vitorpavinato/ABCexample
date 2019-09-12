#!/bin/R

# Settings 
seed <- 1234
set.seed(seed,"Mersenne-Twister")
sample_size           <- 10     # number of diploid individuals
number_of_simulations <- 2

# Write headers for reference table files (summary statistics and parameters)
#sfs_head <- paste0("SFS",1:(2*sample_size-1))
#write(paste(c("pi S TD",sfs_head),collapse=" "),file="results/sumstats.txt")
write("simID pi S TD",file="results/sumstats.txt")
write("simID seed2 N mu bs gamma_m gamma_k theta_1 theta_2",file="results/params.txt")

for (sim in 1:number_of_simulations){
  SEED       <- round(runif(3,0,2^32))
  MU         <- 10^runif(1,-9,-7)
  N          <- round(10^runif(1,log10(SampleSize),2))
  BS         <- rbeta(1,1.5,5)
  GAMMA_M    <- -10^runif(1,-3,1)
  GAMMA_K    <- runif(1,0.001,1)
  G          <- N*100
  
  # run SLiM to generate initial state of simulation
  SLiM_command_1 <- paste0("./bin/slim -s ",SEED[1],
                           " -d N=",N,
                           " -d MU=",MU,
                           " -d BS=",BS,
                           " -d GAMMA_M=",GAMMA_M,
                           " -d GAMMA_K=",GAMMA_K,
                           " -d G=",G,
                           " ./src/SLiM_initial_state_generation.txt")
  system(SLiM_command_1)
  
  # run SLiM inference model
  SLiM_command_2 <- paste0("./bin/slim -s ",SEED[2],
                           " -d simID=",SEED[1],
                           " -d N=",N,
                           " -d MU=",MU,
                           " -d BS=",BS,
                           " -d GAMMA_M=",GAMMA_M,
                           " -d GAMMA_K=",GAMMA_K,
                           " -d G=",G,
                           " -d SAMPLE_SIZE=",sample_size,
                           " ./src/SLiM_inference_model.txt")
  system(SLiM_command_2)
  
  
  # run: pyslim for addition of neutral mutations
  #      tskit for calculating summary statistics (export to vcf alternatively)
  
  NEUTRAL_MU <- MU*(1-BS)
  python_command <- paste("python3 ./src/addNeutralMutation.py",
                          SEED[1],
                          sample_size,
                          N,
                          NEUTRAL_MU, #mutation rate//
                          SEED[3])
  
  system(python_command)
}



