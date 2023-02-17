#the shared script library, where all kinships are declared  
source('~/Documents/deep_learning_for_breeding/genetic_inheritance/make_kinship_matrix/kinships.R')

#the common stats file with number of filtered SNPs for each dataset
stats_infile = '~/Documents/deep_learning_for_breeding/datasets/filtering_stats.csv'

#actual elaboration config
config = NULL

#------- ADDITIVE
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'additive',
  MAF_min = 0.01,
  MAF_max = 0.5,
  force_overwrite = FALSE
))
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'additive',
  MAF_min = 0.05,
  MAF_max = 0.5,
  force_overwrite = FALSE
))
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'additive',
  MAF_min = 0.01,
  MAF_max = 0.05,
  force_overwrite = FALSE
))

#------- DOMINANCE
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'dominance',
  MAF_min = 0.01,
  MAF_max = 0.5,
  force_overwrite = FALSE
))
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'dominance',
  MAF_min = 0.05,
  MAF_max = 0.5,
  force_overwrite = FALSE
))
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'dominance',
  MAF_min = 0.01,
  MAF_max = 0.05,
  force_overwrite = FALSE
))

#------- EPISTASIS AA
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'epistasis_AA',
  MAF_min = 0.01,
  MAF_max = 0.5,
  force_overwrite = FALSE
))
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'epistasis_AA',
  MAF_min = 0.05,
  MAF_max = 0.5,
  force_overwrite = FALSE
))
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'epistasis_AA',
  MAF_min = 0.01,
  MAF_max = 0.05,
  force_overwrite = FALSE
))



#------- EPISTASIS AD
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'epistasis_AD',
  MAF_min = 0.01,
  MAF_max = 0.5,
  force_overwrite = FALSE
))
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'epistasis_AD',
  MAF_min = 0.05,
  MAF_max = 0.5,
  force_overwrite = FALSE
))
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'epistasis_AD',
  MAF_min = 0.01,
  MAF_max = 0.05,
  force_overwrite = FALSE
))



#------- EPISTASIS DD
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'epistasis_DD',
  MAF_min = 0.01,
  MAF_max = 0.5,
  force_overwrite = FALSE
))
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'epistasis_DD',
  MAF_min = 0.05,
  MAF_max = 0.5,
  force_overwrite = FALSE
))
config = rbind(config, data.frame(
  base_folder = '~/research/deep_learning_for_breeding/datasets/cattle/',
  kinship_algo = 'epistasis_DD',
  MAF_min = 0.01,
  MAF_max = 0.05,
  force_overwrite = FALSE
))

