#read the passed .csv.gz SNP file (one row per SNP, one column per sample, first row
#and column used as SNP/sample names, values in 0/1/2)
#creates a new filtered file removing fixed markers (fully het too) pluse 
#the SNPs that have a Minor Allele Frequency (MAF)
#equal or less than the passed value. If MAF=0 only fixed markers are removed



# INPUT MANAGEMENT --------------------------------------------------------
#loading the parameters from the first argument
args = commandArgs(trailingOnly = TRUE)

#we have a secret default config, just for testing purposes
if (length(args) != 2){
  stop('Usage: Rscript genotype_MAF_filtering.R <path/to/SNPs.csv.gz> <min_MAF_allowed>')
}

#naming arguments
infile = args[1]
MAF_threshold = as.numeric(args[2])

#checking arguments
if(!file.exists(infile)){
  stop(paste('SNP file does not exists:', infile))
}
if (!grepl('.csv.gz$', infile)){
  stop(paste('SNP file should be a .csv.gz'))
}
if (is.na(MAF_threshold)){
  stop(paste('<min_MAF_allowed> should be numeric between 0 and 1, instead I found', min_MAF_allowed))
}
if (MAF_threshold > 1 || MAF_threshold < 0){
  stop(paste('<min_MAF_allowed> should be numeric between 0 and 1, instead I found', min_MAF_allowed))
}

#reading data in
SNP = read.csv(infile, stringsAsFactors = FALSE, row.names = 1)

# FIXED MARKERS -----------------------------------------------------------
SNP_min = apply(SNP, MARGIN = 1, FUN = min)
SNP_max = apply(SNP, MARGIN = 1, FUN = max)
fixed = SNP_max == SNP_min
SNP = SNP[!fixed, ]

# MAF ---------------------------------------------------------------------
#computing MAF
MAF = rowSums(SNP) / (2*ncol(SNP))
MAF[MAF > 0.5] =  1 - MAF[MAF > 0.5]

#selecting markers to keep
good = MAF > MAF_threshold
SNP = SNP[good, ]

# OUTPUT ------------------------------------------------------------------
#building the filename for outfile
suffix = paste(sep='', '_filteredMAF', MAF_threshold, '.csv.gz')
outfile = gsub(infile, pattern = '.csv.gz', replacement = suffix)

#saving
fp = gzfile(outfile, open = 'w')
write.csv(SNP, file = fp)
close(fp)
writeLines(paste('File written to', outfile))

