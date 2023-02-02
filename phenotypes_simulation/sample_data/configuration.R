#dataset folder
root = '/home/nelson/research/genetic_inheritance/phenotypes_simulation/sample_data'

#markers matrix input file, expected in 0/1/2, markers x samples
SNP_infile = file.path(root, 'SNPs.csv.gz')

#--- output files
#simulated phenotypes
pheno_outfile = file.path(root, 'simulation_1', 'simulated_phenotype.csv')
#SNP matrix without the causative effect (we are always saving a .gz file)
SNP_outfile = file.path(root, 'simulation_1', 'SNPs_nonCausative.csv.gz') 
#text file with list of causative SNPs
SNPlist_outfile = file.path(root, 'simulation_1', 'SNPs_causative.list.txt') 
#parameter used for the simulation, just for checks
params_outfile = file.path(root, 'simulation_1', 'usedpars.txt') 

#--- simulation parameters
#overall mean
pheno_mean = 0
#overall heritability
pheno_h2 = 0.7
#coefficient of variation of QTN
QTN_cv = 0.1
#number of QTN in the main effect (must be even, because epistasis
#works in pairs of SNPs)
QTN_num = 4
#mean absolute effect of QTN, one element per simulated scenario
#actual final effect has 50% chance of being negative,
#i.e. being multiplied by -1
QTN_additive_mean            = c(10)
QTN_dominance_mean           = c( 0)
QTN_additive_additive_mean   = c( 0)
QTN_additive_dominance_mean  = c( 0)
QTN_dominance_additive_mean  = c( 0)
QTN_dominance_dominance_mean = c( 0)

#--- utility
#tmp folder: for SimPhe interface I'll need to write there genos and parameters
#files, and then read everything again...
#The following code uses the standard tempdir() function to get a suitable location
#but you may want to move it somewhere else
tmp_folder = tempdir(check = TRUE)
SNP_file_simphe = file.path(tmp_folder, 'genotypes.txt')
param_file_simphe = file.path(tmp_folder, 'params.txt')
