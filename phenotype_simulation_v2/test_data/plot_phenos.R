#building a phenotype plot that works only for this specific configuration of
#test data

library(stringi)
library(ggplot2)
library(tidyr)
library(readr)
setwd('~/research/genetic_inheritance/phenotype_simulation_v2/test_data/simulation_1/')


# SUPPORT FUNCTIONS -------------------------------------------------------
simplify_phenoname = function(x){
  x = gsub(x = x, pattern = 'simphe_', replacement = '')
  x = gsub(x = x, pattern = 'mean0_', replacement = '')
  x = gsub(x = x, pattern = 'hSquare1_cv0_QTN2_RandomSignFALSE_', replacement = '')
  x = gsub(x = x, pattern = 'testModeTRUE_', replacement = '')
  x = gsub(x = x, pattern = '_epoch[0-9]*', replacement = '')
  return(x)
}

# ACTUAL SCRIPT -----------------------------------------------------------
#reading data in
phenos = read_csv(file = 'simulated_phenotype.csv', col_names = TRUE, col_types = NULL, show_col_types = FALSE)
colnames(phenos) = simplify_phenoname(colnames(phenos)) 
pieces = data.frame(stri_split_fixed(str = phenos$sample, pattern = '_', simplify = TRUE))
phenos$locus_1 = pieces$X2
phenos$locus_2 = pieces$X3
phenos$sample = NULL

#a form palatable to ggplot2
phenos = pivot_longer(data = phenos, 
                      cols = setdiff(colnames(phenos), c('locus_1', 'locus_2')), 
                      names_to = 'phenoConfig', values_to = 'phenoValue')

#forcing facet order
phenos$phenoConfig = factor(x = phenos$phenoConfig, levels = c(
  'Dfraction0_EpiAddOver0.25',
  'Dfraction0_EpiAddOver0.5',
  'Dfraction0_EpiAddOver0.75',
  'Dfraction0_EpiAddOver-0.25',
  'Dfraction0_EpiAddOver-0.5',
  'Dfraction0_EpiAddOver-0.75',
  'Dfraction0_EpiAddOver0'
))



#doing the plot
p = ggplot(phenos, aes(
  x = locus_1, 
  y=phenoValue, 
  group = locus_2, color = locus_2, shape = locus_2)) + 
  geom_point(size=5) + geom_line() + facet_wrap(facets = vars(phenoConfig))
print(p)

ggsave(filename = 'epistasis_additive.png', plot = p)
