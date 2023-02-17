#this list will contain all the configurations
conf = list()

#dataset folder
conf$root = '/home/nelson/research/genetic_inheritance/phenotype_simulation_v2/test_data'

#markers matrix input file, expected in 0/1/2, markers x samples
conf$SNP_infile = file.path(conf$root, 'SNPs.csv.gz')

#--- output files
#simulated phenotypes
conf$pheno_outfile = file.path(conf$root, 'simulation_1', 'simulated_phenotype.csv')
#SNP matrix without the causative effect (we are always saving a .gz file)
conf$SNP_outfile = file.path(conf$root, 'simulation_1', 'SNPs_nonCausative.csv.gz') 
#text file with list of causative SNPs
conf$SNPlist_outfile = file.path(conf$root, 'simulation_1', 'SNPs_causative.list.txt') 
#parameter used for the simulation, just for checks
conf$params_outfile = file.path(conf$root, 'simulation_1', 'usedpars.txt') 

#--- simulation parameters
#overall mean
conf$pheno_mean = 0
#number of QTN in the main effect (must be even, because epistasis
#works in pairs of SNPs)
conf$QTN_num = 2

#overall heritability
conf$pheno_h2 = 0.7
#coefficient of variation of QTN
conf$QTN_cv = 0.1
#flip allele effect, put to TRUE if you want the sign of the effect of the 
#counted alleles to be randomly assigned, so that a QTN can have either positive
#or negative effect. Put to FALSE to deactivate the functionality (the counted
#allele will always have a positive effect)
conf$QTN_random_sign = TRUE

#intensity of dominance, set to zero to have completely additive effect, 
#set to one to have completely dominant (heterozygous equal to the dominant
#homozygous). Intermediate values create intermediate conditions
conf$QTN_dominance_fraction = 0

#intensity of additive-additive epistasis
conf$QTN_epistasis_additive_overexpression = -0.25
#--- utility (not really interesting for end user)
#tmp folder: for SimPhe interface I'll need to write there genos and parameters
#files, and then read everything again...
#The following code uses the standard tempdir() function to get a suitable location
#but you may want to move it somewhere else
conf$tmp_folder = tempdir(check = TRUE)
conf$SNP_file_simphe = file.path(conf$tmp_folder, 'genotypes.txt')
conf$param_file_simphe = file.path(conf$tmp_folder, 'params.txt')

#test mode: when TRUE all randomness is removed, heritability goes to one,
#cv goes to zero, no random flip of effect sign, causative SNP are the 
#first conf$QTN_num taken in order (and not picked 
#at random), and epistatic pairings are fixed (first and second SNPs, third and 
#fourth SNPs, and so forth)
conf$test_mode = TRUE