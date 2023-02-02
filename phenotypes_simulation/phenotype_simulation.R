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
source(args[1])

#Below this point is just code, no need to actually read it if you are using
#the script as a black box

# SUPPORT FUNCTIONS -------------------------------------------------------
pick_main_effects = function(SNP_names, QTN, A_mean, D_mean, cv){
  #sampling the markers with an effect
  SNP_selected = sample(x = SNP_names, size = QTN, replace = FALSE)

  #generating effects: additive and dominant
  A_eff = rnorm_posneg_effect(n = QTN, mean = A_mean, cv = cv)
  D_eff = rnorm_posneg_effect(n = QTN, mean = D_mean, cv = cv)
  
  #putting all together
  return(data.frame(
    SNP = SNP_selected,
    additive = A_eff,
    dominance = D_eff
  ))
}

#picks from the normal distribution with specification of mean and coefficient
#of variation, plus results have 50% chance to have the sign flipped
rnorm_posneg_effect = function(n, mean, cv){
  res = rnorm(n = n, mean = mean, sd = mean * cv)
  sel = runif(n = n) > 0.5
  res[sel] = res[sel] * -1
  return(res)
}

#SNP with epistasis effect are the same as those with the main effect,
#picked in pairs
pick_epistasis_effects = function(SNP_names, AA_mean, AD_mean, DA_mean, DD_mean, cv){
  QTN_epi = length(SNP_names) / 2
  SNPA_selected = SNP_names[1:QTN_epi]
  SNPB_selected = SNP_names[(QTN_epi + 1):(2 * QTN_epi)]
  
  #generating effects: combinations of additive and dominant
  AA_eff = rnorm_posneg_effect(n = QTN_epi, mean = AA_mean, cv = cv)
  AD_eff = rnorm_posneg_effect(n = QTN_epi, mean = AD_mean, cv = cv)
  DA_eff = rnorm_posneg_effect(n = QTN_epi, mean = DA_mean, cv = cv)
  DD_eff = rnorm_posneg_effect(n = QTN_epi, mean = DD_mean, cv = cv)
  
  #putting all together
  return(data.frame(
    SNPA = SNPA_selected,
    SNPB = SNPB_selected,
    additive_additive = AA_eff,
    additive_dominance = AD_eff,
    dominance_additive = DA_eff,
    dominance_dominance = DD_eff
  ))
}

# FOLDER PREP -------------------------------------------------------------
dir.create(dirname(pheno_outfile), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(SNP_outfile), recursive = TRUE, showWarnings = FALSE)

# GENOTYPES PREP ----------------------------------------------------------
#loading genotypes, bringing them to the SimPhe format
writeLines('Preparing SNP data...')
SNP = read.csv(SNP_infile, row.names = 1, stringsAsFactors = FALSE)
SNP = data.frame(t(SNP))
writeLines('DONE')

# PHENOTYPE SIMULATION ----------------------------------------------------
#epoch for this run (so that we can have multiple repetitions of the
#current config and still have unique phenotype names)
epoch = as.integer(Sys.time())

#number of main QTN must be even
if (QTN_num %% 2 == 1){
  stop('Number of QTN must be even')
}

#if there are already phenotypes, we are going to add to the table
pheno_values = NULL
if(file.exists(pheno_outfile)){
  pheno_values = read.csv(pheno_outfile, stringsAsFactors = FALSE)
}

#building a string that will be used to name this phenotype and for interface purposes
phenoname = paste(sep='', 'simphe',
                  '_mean', pheno_mean,
                  '_hSquare', pheno_h2,
                  '_cv', QTN_cv,
                  '_QTN', QTN_num,
                  '_A', QTN_additive_mean,
                  '_D', QTN_dominance_mean,
                  '_AA', QTN_additive_additive_mean,
                  '_AD', QTN_additive_dominance_mean,
                  '_DA', QTN_dominance_additive_mean,
                  '_DD', QTN_dominance_dominance_mean,
                  '_epoch', epoch
)
writeLines(paste(sep='', ' - ', phenoname))

#-------- building the param file - preface
fp = file(param_file_simphe, open = 'w')
param = c(
  #----mean
  '[P1mean]',
  'mean',
  paste(pheno_mean),
  ''
)
writeLines(param, con = fp)
close(fp)

#-------- building the param file - main part
fp = file(param_file_simphe, open = 'a')
param = c(
  '[P1main]',
  'SNP additive dominance'
)
writeLines(param, con = fp)
close(fp)

main_effects = pick_main_effects(SNP_names = colnames(SNP), 
                  QTN = QTN_num, 
                  A_mean = QTN_additive_mean,
                  D_mean = QTN_dominance_mean,
                  cv = QTN_cv
                  )
write.table(main_effects, file = param_file_simphe, sep = '\t', row.names = FALSE, quote = FALSE, append = TRUE, col.names = FALSE)

#-------- building the param file - epistasis part
fp = file(param_file_simphe, open = 'a')
param = c(
  '',
  '',
  '[P1epistasis]',
  'SNPA SNPB additive_additive additive_dominance dominance_additive dominance_dominance'
)
writeLines(param, con = fp)
close(fp)
epi_effects = pick_epistasis_effects(SNP_names = main_effects$SNP, 
                                 AA_mean = QTN_additive_additive_mean,
                                 AD_mean = QTN_additive_dominance_mean,
                                 DA_mean = QTN_dominance_additive_mean,
                                 DD_mean = QTN_dominance_dominance_mean,
                                 cv = QTN_cv
)
write.table(epi_effects, file = param_file_simphe, sep = '\t', row.names = FALSE, quote = FALSE, append = TRUE, col.names = FALSE)

#-------- building the param file - heritability part
fp = file(param_file_simphe, open = 'a')
param = c(
  '',
  '[P1heritability]',
  'heritability',
  pheno_h2
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
    sim.pars = param_file_simphe, 
    fgeno = SNP, 
    fusepar = params_outfile, 
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
newcolnames = c(colnames(pheno_values), phenoname)
if (is.null(pheno_values)){
  pheno_values = phenotypes
}else{
  pheno_values = cbind(pheno_values, phenotypes)  
}
colnames(pheno_values) = newcolnames

#saving, overwriting at each iteration
write.csv(pheno_values, pheno_outfile, row.names = FALSE)

# GENOTYPES EXCLUDING CAUSATIVE SNPS --------------------------------------
#we need to save the SNP matrix removing all the selected causative SNPs
noncausative = setdiff(colnames(SNP), main_effects$SNP)

#are we going to output a compressed
fp = gzfile(description = SNP_outfile, open = 'w')
write.csv(x = t(SNP[,noncausative]), file = fp)
close(fp)

#we also write the list of causative SNPs, just to be sure
writeLines(main_effects$SNP, con = SNPlist_outfile)