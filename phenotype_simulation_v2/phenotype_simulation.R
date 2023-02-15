#install.packages('SimPhe')
library(SimPhe)

#steps to success
# - decide the number of markers for the main effects (additive, dominant)
# - those same markers will be used, in pairs, for epistasis
# - decide trait heritability (influences added final noise)
# - decide overall mean (not strictly enforced by simphe)
# - decide a plan for effect sizes

# INPUT PARAMETERS --------------------------------------------------------

#a number of parameters are read from the configuration file and are thus 
#expected to be defined there. If something doesn't work the code crashes :)
#Please see the sample configuration file for the list and explanation of
#the expected parameters

#loading the parameters from the first argument
args = commandArgs(trailingOnly = TRUE)

#we have a secret default config, just for testing purposes
if (length(args) == 0){
  test_config_file = '/home/nelson/research/genetic_inheritance/phenotype_simulation_v2/test_data/configuration.R'
  args[1] = test_config_file
  writeLines(paste('No config file provided, using default:', test_config_file))
}
source(args[1])

#Below this point is just code, no need to actually read it if you are using
#the script as a black box

# SUPPORT FUNCTIONS -------------------------------------------------------
pick_main_effects = function(SNP_causative, QTN_dominance_fraction, cv, QTN_random_sign){
  #generating effects: additive and dominant
  A_eff = rnorm_posneg_effect(n = length(SNP_causative), mean = 100, cv = cv, random_flip = QTN_random_sign)
  D_eff = A_eff * QTN_dominance_fraction
  #D_eff = rnorm_posneg_effect(n = QTN, mean = D_mean, cv = cv, random_flip = QTN_random_sign)
  #A_eff = A_eff + D_eff
  
  #putting all together
  return(data.frame(
    SNP = SNP_causative,
    additive = A_eff,
    dominance = D_eff
  ))
}

#picks from the normal distribution with specification of mean and coefficient
#of variation. Optionally, results have 50% chance to have the sign flipped
rnorm_posneg_effect = function(n, mean, cv, random_flip){
  res = rnorm(n = n, mean = mean, sd = mean * cv)
  if (random_flip){
    sel = runif(n = n) > 0.5
    res[sel] = res[sel] * -1
  }
  return(res)
}

#SNP with epistasis effect are the same as those with the main effect,
#picked in pairs
#pick_epistasis_effects = function(SNP_names, AA_mean, AD_mean, DA_mean, DD_mean, cv){
pick_epistasis_effects = function(main_effects, QTN_epistasis_additive_overexpression){
  warning('Stub function')
  #epistatic pairs are determined by alternating SNPs with main effects: first
  #and second SNPs, third and fourth SNPs, and so on
  SNP_names = main_effects$SNP
  #SNPA_selected = SNP_names[seq(from=1, to=nrow(main_effects), by=2)]
  #SNPB_selected = SNP_names[seq(from=2, to=nrow(main_effects), by=2)]
  SNPA_selected = seq(from=1, to=nrow(main_effects), by=2)
  SNPB_selected = seq(from=2, to=nrow(main_effects), by=2)
  
  #generating effects: combinations of additive and dominant
  AA_eff = (main_effects$additive[SNPA_selected] + main_effects$additive[SNPB_selected]) * QTN_epistasis_additive_overexpression
  #AD_eff = rnorm_posneg_effect(n = QTN_epi, mean = AD_mean, cv = cv)
  #DA_eff = rnorm_posneg_effect(n = QTN_epi, mean = DA_mean, cv = cv)
  #DD_eff = rnorm_posneg_effect(n = QTN_epi, mean = DD_mean, cv = cv)
  
  #putting all together
  return(data.frame(
    SNPA = SNP_names[SNPA_selected],
    SNPB = SNP_names[SNPB_selected],
    additive_additive = AA_eff,
    additive_dominance = 0,
    dominance_additive = 0,
    dominance_dominance = 0
  ))
}

# FOLDER PREP -------------------------------------------------------------
dir.create(dirname(conf$pheno_outfile), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(conf$SNP_outfile), recursive = TRUE, showWarnings = FALSE)

# GENOTYPES PREP ----------------------------------------------------------
#loading genotypes, bringing them to the SimPhe format
writeLines(paste('Loading SNP data:', conf$SNP_infile))
SNP = read.csv(conf$SNP_infile, row.names = 1, stringsAsFactors = FALSE)
SNP = data.frame(t(SNP))
writeLines('DONE')

#fixed SNP create NAs in SimPhe, we cannot accept them
SNP_min = apply(SNP, MARGIN = 2, FUN = min)
SNP_max = apply(SNP, MARGIN = 2, FUN = max)
if (any(SNP_max == SNP_min)){
  stop('Some SNP markers are fixed, SimPhe would crash with this kind of data. Remove the fixed SNP are try again.')
}
rm(SNP_min, SNP_max)

# PHENOTYPE SIMULATION ----------------------------------------------------
#are we in test mode? if so, some parameters should be overwritten and the
#causative SNPs should not picked at random
if (conf$test_mode){
  conf$pheno_h2 = 1
  conf$QTN_cv = 0
  conf$QTN_random_sign = FALSE
  SNP_causative =  colnames(SNP)[1:conf$QTN_num]
}else{
  #not in test mode, picking causative SNPs at random
  SNP_causative =  sample(x = colnames(SNP), size = conf$QTN_num, replace = FALSE)
}

#epoch for this run (so that we can have multiple repetitions of the
#current config and still have unique phenotype names)
epoch = as.integer(Sys.time())

#number of main QTN must be even
if (conf$QTN_num %% 2 == 1){
  stop(paste('Number of QTN must be even, currently QTN_num =', conf$QTN_num))
}

#if there are already phenotypes, we are going to add to the table
pheno_values = NULL
if(file.exists(conf$pheno_outfile)){
  pheno_values = read.csv(conf$pheno_outfile, stringsAsFactors = FALSE, row.names = 1)
}

#building a string that will be used to name this phenotype and for interface purposes
phenoname = paste(sep='', 'simphe',
                  '_mean', conf$pheno_mean,
                  '_hSquare', conf$pheno_h2,
                  '_cv', conf$QTN_cv,
                  '_QTN', conf$QTN_num,
                  '_RandomSign', conf$QTN_random_sign,
                  '_Dfraction', conf$QTN_dominance_fraction,
                  '_EpiAddOver', conf$QTN_epistasis_additive_overexpression,
                  '_testMode', conf$test_mode,
                  '_epoch', epoch
)
writeLines(paste(sep='', ' - ', phenoname))

#-------- building the param file - preface
fp = file(conf$param_file_simphe, open = 'w')
param = c(
  #----mean
  '[P1mean]',
  'mean',
  paste(conf$pheno_mean),
  ''
)
writeLines(param, con = fp)
close(fp)

#-------- building the param file - main part
fp = file(conf$param_file_simphe, open = 'a')
param = c(
  '[P1main]',
  'SNP additive dominance'
)
writeLines(param, con = fp)
close(fp)

main_effects = pick_main_effects(SNP_causative = SNP_causative,
                                 QTN_dominance_fraction = conf$QTN_dominance_fraction,
                                 QTN_random_sign = conf$QTN_random_sign,
                                 cv = conf$QTN_cv)
write.table(main_effects, file = conf$param_file_simphe, sep = '\t', row.names = FALSE, quote = FALSE, append = TRUE, col.names = FALSE)

#-------- building the param file - epistasis part
fp = file(conf$param_file_simphe, open = 'a')
param = c(
  '',
  '',
  '[P1epistasis]',
  'SNPA SNPB additive_additive additive_dominance dominance_additive dominance_dominance'
)
writeLines(param, con = fp)
close(fp)
epi_effects = pick_epistasis_effects(main_effects = main_effects, QTN_epistasis_additive_overexpression = conf$QTN_epistasis_additive_overexpression)

write.table(epi_effects, file = conf$param_file_simphe, sep = '\t', row.names = FALSE, quote = FALSE, append = TRUE, col.names = FALSE)

#-------- building the param file - heritability part
fp = file(conf$param_file_simphe, open = 'a')
param = c(
  '',
  '[P1heritability]',
  'heritability',
  conf$pheno_h2
)
writeLines(param, con = fp)
close(fp)

#------- invoking simphe
#here we need a safeguard: sometimes sim.phe fails with warning:
#In rnorm(nrow(geno), sd = sqrt(exp.noise.var)) : NAs produced
#and returns NAs
maximum_cycles = 3
while(TRUE){
  phenotypes = sim.phe(
    sim.pars = conf$param_file_simphe, 
    fgeno = SNP, 
    fusepar = conf$params_outfile, 
    ftype = "snp.head", 
    seed = floor(runif(1) * 100000), #this is so each invocation is different 
    fwrite = FALSE, 
    pattern = "[[:alpha:]]+",
    genetic.model	 = 'epistasis')
  
  #if data is good let's break the cycle
  if (sum(is.na(phenotypes[,1])) == 0){
    break
  }
  
  #but we are not doing this more than 100 times
  maximum_cycles = maximum_cycles - 1
  if (maximum_cycles == 0){
    stop('sim.phe keeps returning NAs, stopping')
  }
  
  #if we get here we have to cycle again and try another round of generation
  writeLines('sim.phe returned NAs, trying again')
}

#putting phenotype with already present ones
newcolnames = c('sample', colnames(pheno_values), phenoname)
if (is.null(pheno_values)){
  pheno_values = phenotypes
}else{
  pheno_values$sample = NULL
  pheno_values = cbind(pheno_values, phenotypes)  
}
pheno_values = as_tibble(cbind(data.frame(sample = rownames(SNP)), pheno_values))
colnames(pheno_values) = newcolnames

#saving, overwriting at each iteration
write.csv(pheno_values, conf$pheno_outfile, row.names = FALSE)

# GENOTYPES EXCLUDING CAUSATIVE SNPS --------------------------------------
#we need to save the SNP matrix removing all the selected causative SNPs
noncausative = setdiff(colnames(SNP), main_effects$SNP)

#are we going to output a compressed
fp = gzfile(description = conf$SNP_outfile, open = 'w')
write.csv(x = t(SNP[,noncausative, drop = FALSE]), file = fp)
close(fp)

#we also write the list of causative SNPs, just to be sure
writeLines(main_effects$SNP, con = conf$SNPlist_outfile)