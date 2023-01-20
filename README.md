# Mechanisms of genetic inheritance
Exploring different types of genetic inheritance: additivity, dominance and epistasis

This repository hosts code, workflows and ideas on the use of simulated and real phenotypic and genetic data for the study of genetic inheritance.
Different scenarios are simulated, considering additive. dominant and epistatic phenotypes. Genetic data are then used to predict and/or approximate the simulated phenotypes by accounting for additivity, dominance and espistasis.

## Operational steps

1. Use data from [Nazzicari & Biscarini 2022](https://www.nature.com/articles/s41598-022-24405-0) (~1,000 samples) to simulate phenotypes under scenarios of extreme epistasis, and add these to the existing fully additive and fully dominant simulated phenotypes
2. Use data from the $1,000$ [genomes project](https://www.internationalgenome.org/) to explore the same scenario on a different population (still diploid genomes, but human)
3. Use synthetic data from [HAPNEST](https://www.ebi.ac.uk/biostudies/studies/S-BSST936) (~1 million samples) to explore the effect of data scale on the results

[Infinium Global Screening Array](https://www.illumina.com/products/by-type/microarray-kits/infinium-global-screening.html) with $654,027$ SNP variants
