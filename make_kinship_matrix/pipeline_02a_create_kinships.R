

# INPUT CONFIGURATION MANAGEMENT ------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 1){
  #loading the parameters
  source(args[1])
}else{
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #the shared script library, where all kinships are declared  
  source('~/Documents/deep_learning_for_breeding/genetic_inheritance/make_kinship_matrix/kinships.R')
  
  #the common stats file with number of filtered SNPs for each dataset
  stats_infile = '~/Documents/deep_learning_for_breeding/deep_learning_for_breeding/datasets/filtering_stats.csv'
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    base_folder = '~/Documents/deep_learning_for_breeding/genetic_inheritance/data',
    kinship_algo = 'additive',
    MAF_min = 0.05,
    MAF_max = 0.06,
    force_overwrite = FALSE
  ))
  
  config = rbind(config, data.frame(
    base_folder = '~/Documents/deep_learning_for_breeding/genetic_inheritance/data',
    kinship_algo = 'dominance',
    MAF_min = 0.05,
    MAF_max = 0.06,
    force_overwrite = FALSE
  ))
}

# SETUP -------------------------------------------------------------------
#if not installed, use package remotes, function install_gitlab()
library("sommer")

# ACTUAL PIPELINE ---------------------------------------------------------
#stats on SNPs, to be updated if needed
stats = NULL
if (file.exists(stats_infile)){
  stats = read.csv(stats_infile, row.names = NULL, stringsAsFactors = FALSE)
}

#for each required operation
for (i in 1:nrow(config)){
  #current configuration
  curr = config[i,]
  print(curr)
  
  #building the target output file
  outfile = paste(sep='', curr$base_folder, '/kinships/kinship_', curr$kinship_algo, 
                  '_minMAF', curr$MAF_min, '_maxMAF', curr$MAF_max, '.csv.gz')
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  
  #should we proceed?
  if (file.exists(outfile) & curr$force_overwrite == FALSE){
    writeLines('Required configuration already present. To recreate set force_overwrite to TRUE')
    next()
  }
  
  #is it a brand new file? or are we overwriting an old one?
  if (file.exists(outfile) & curr$force_overwrite == TRUE){
    writeLines('Recreating kinship file')
  }else{
    writeLines('Creating new kinship file')
  }
  
  #loading imputed file
  infile = paste(sep='', curr$base_folder, '/SNPs.csv.gz')
  if(!file.exists(infile)){
    stop(paste('Missing expected file', infile))
  }
  writeLines(' - loading data')
  genos = read.csv(infile, row.names = 1)
  
  #filtering on MAF
  writeLines(' - filtering on MAF')
  genos_filtered = filter_on_MAF(genotypes = genos, min.MAF = curr$MAF_min, max.MAF = curr$MAF_max, ploidy = 2)
  
  #stats
  stats = rbind(stats, data.frame(
    dataset = basename(curr$base_folder),
    samples = ncol(genos),
    SNP_original = nrow(genos),
    MAF_min = curr$MAF_min,
    MAF_max = curr$MAF_max,
    SNP_filtered = nrow(genos_filtered)
  ))
  stats = unique(stats)
  write.csv(stats, stats_infile, row.names = FALSE)
  
  #computing kinship
  writeLines(' - computing kinship')
  kinship_function = paste(sep='', 'kinship_', curr$kinship_algo)
  kinship_args = list(SNP_matrix = genos_filtered)
  kinship_computed = do.call(what = kinship_function, args = kinship_args)
  
  #writing the output
  writeLines(' - writing outfile')
  fp = gzfile(outfile, open = 'w')
  write.csv(kinship_computed, file = fp)
  close(fp)
}
