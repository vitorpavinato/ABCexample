#!/usr/bin/python

# Add neutral mutations to tree simulated with SLiM

import msprime, pyslim, sys, random

simID = int(sys.argv[1])
SAMPLE_SIZE = int(sys.argv[2])
N = int(sys.argv[3])
MU = float(sys.argv[4])
SEED = int(sys.argv[5])

#sample=random.sample(population=range(1,N), k=SAMPLE_SIZE)

ts = pyslim.load("/tmp/slim_%i.trees" % simID)
mutated = pyslim.SlimTreeSequence(msprime.mutate(ts, rate=MU, random_seed=SEED, keep=True))
#mutated_sample = mutated.samples()[sample]

#with open("results/sim_data.vcf", "w") as vcf_file:
#    mutated.write_vcf(vcf_file,individuals=range(1, SAMPLE_SIZE+1))


pi = mutated.diversity(span_normalise=False)
S = int(mutated.segregating_sites(span_normalise=False))
TD = mutated.Tajimas_D()
#sfs = mutated.allele_frequency_spectrum(sample_sets=mutated_sample,span_normalise=False)


ss_file = open("results/sumstats.txt","a")
ss_file.write("%i %i %.10f %i %.10f\n" % (simID,SEED,pi,S,TD))
ss_file.close() 



