#all kinship functions should follow the template:
#kinship_SOMETHING = function(SNP_matrix){...}
#where "SOMETHING" is the name of the algorithm or some distinctive string
#and SNP_matrix is SNP matrix, markers as rows, samples as columns, values in 0/1/2
#(or more) in case of polyploids

library("sommer")
# library("madbito")

## ADDITIVE
kinship_AB = function(SNP_matrix){
  #using madbito's implementation, which requires a tranposition
  return(kinship.AB(t(SNP_matrix)))
}

# kinship_dominance = function(SNP_matrix){
#   #using madbito's implementation, which requires a tranposition
#   res = kinship.A_plus_D(t(SNP_matrix))
#   return(res$dominant)
# }

kinship_VR = function(SNP_matrix){
  #using madbito's implementation, which requires a tranposition
  return(kinship.VR(t(SNP_matrix)))
}

kinship_additive = function(SNP_matrix){
  #using sommer's implementation, which requires -1/0/1 coded SNPs
  #transposing the input matrix: from (m,n) to (n,m)
  #endelman=FALSE --> VanRaden
  M = as.matrix(t(SNP_matrix)) - 1
  res = A.mat(M, endelman = FALSE, min.MAF = 0)
  return(res)
}

## DOMINANCE
kinship_dominance = function(SNP_matrix) {
  #using sommer's implementation, which requires -1/0/1 coded SNPs
  #transposing the input matrix: from (m,n) to (n,m)
  M = as.matrix(t(SNP_matrix)) - 1
  res = D.mat(M, nishio = TRUE, min.MAF = 0)
  return(res)
}

## EPISTASIS
kinship_epistasis = function(SNP_matrix, epi_type) {
  #using sommer's implementation, which requires -1/0/1 coded SNPs
  #transposing the input matrix: from (m,n) to (n,m)
  M = as.matrix(t(SNP_matrix)) - 1
  res = E.mat(M, endelman = FALSE, nishio = TRUE, min.MAF = 0, type = epi_type)
  return(res)
}

kinship_epistasis_AA = function(SNP_matrix) {
  res = kinship_epistasis(SNP_matrix, "A#A")
  return(res)
}

kinship_epistasis_AD = function(SNP_matrix) {
  res = kinship_epistasis(SNP_matrix, "A#D")
  return(res)
}

kinship_epistasis_DD = function(SNP_matrix) {
  res = kinship_epistasis(SNP_matrix, "D#D")
  return(res)
}

