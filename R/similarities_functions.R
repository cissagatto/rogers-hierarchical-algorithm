##############################################################################
# Generate ROGERS Partitions                                                 #
# Copyright (C) 2022                                                         #
#                                                                            #
# This program is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by the      #
# Free Software Foundation, either version 3 of the License, or (at your     #
# option) any later version. This program is distributed in the hope that    #
# it will be useful, but WITHOUT ANY WARRANTY; without even the implied      #
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   #
# GNU General Public License for more details.                               #
#                                                                            #
# Elaine Cecilia Gatto | Prof. Dr. Ricardo Cerri | Prof. Dr. Mauri Ferrandin #
# Federal University of Sao Carlos (UFSCar: https://www2.ufscar.br/) Campus  #
# Sao Carlos Computer Department (DC: https://site.dc.ufscar.br/)            #
# Program of Post Graduation in Computer Science                             #
# (PPG-CC: http://ppgcc.dc.ufscar.br/) Bioinformatics and Machine Learning   #
# Group (BIOMAL: http://www.biomal.ufscar.br/)                               #
#                                                                            #
##############################################################################


########################################################################
# WORSKSPACE
########################################################################
FolderRoot = "~/Generate-Partitions-Rogers"
FolderScripts = "~/Generate-Partitions-Rogers/R"


#####################################################################
reorder_mat_cor <- function(matrix_correlation){
  #cat("\n\tReorder Matrix")
  dd <- as.dist((1-matrix_correlation)/2)
  hc <- hclust(dd)
  print(hc)
  matrix_correlation <- matrix_correlation[hc$order, hc$order]
  return(matrix_correlation)
  cat("\n")
  gc()
}

#####################################################################
get_upper_tri <- function(matrix_correlation){
  #cat("\n\tGet Upper Tri")
  matrix_correlation[lower.tri(matrix_correlation)]<- NA
  return(matrix_correlation)
  cat("\n")
  gc()
}



###############################################################################
#' Compute marginal probabilities
#'
#' @family MultiLabel Binary Measures Functions
#' @param labels The label space from the dataset
#' @param num.labels total number of labels from label space
#' @param  a frequency which label i and label j occur together
#' @param  b frequency which label i occur alone
#' @param  c frequency which label j occur alone
#' @param  d frequency which label i and label j not occur together
#' @return return the computed marginal probabilities for all labels
#'
###############################################################################
compute.marg.probs <- function(labels, num.labels, a, b, c, d){

  retorno = list()

  mab <- build.matrix.corr(num.labels, labels)
  mac <- build.matrix.corr(num.labels, labels)
  mad <- build.matrix.corr(num.labels, labels)
  mbc <- build.matrix.corr(num.labels, labels)
  mbd <- build.matrix.corr(num.labels, labels)
  mcd <- build.matrix.corr(num.labels, labels)
  mn <- build.matrix.corr(num.labels, labels)

  u = (num.labels*num.labels)
  pb <- progress_bar$new(total = u)

  #i = 1
  #j = 1
  for (i in 1:num.labels){
    for (j in 1:num.labels){

      x = a[i,j]
      y = b[i,j]
      mab[i,j] = compute.ab(x,y)
      #cat("\nAB ", mab[i,j])

      w = a[i,j]
      v = c[i,j]
      mac[i,j] = compute.ac(w,v)
      #cat("\nAC ", mac[i,j])

      e = a[i,j]
      f = d[i,j]
      mad[i,j] = compute.ad(e,f)
      #cat("\nAD ", mad[i,j])

      g = b[i,j]
      h = c[i,j]
      mbc[i,j] = compute.bc(g,h)
      #cat("\nBC ", mbc[i,j])

      k = b[i,j]
      l = d[i,j]
      mbd[i,j] = compute.bd(k,l)
      #cat("\nBD ", mbd[i,j])

      m = c[i,j]
      n = d[i,j]
      mcd[i,j] = compute.cd(m,n)
      #cat("\nCD ", mcd[i,j])

      o = a[i,j]
      p = b[i,j]
      q = c[i,j]
      r = d[i,j]
      mn[i,j] = compute.n(o,p,q,r)
      #cat("\nN ", mn[i,j])

      pb$tick()
      Sys.sleep(1/u)

      #i = i + 1
      gc()
    } # end intern for
    #j = j + 1
    gc()
  } # enf extern for

  retorno$mab = mab
  retorno$mac = mac
  retorno$mad = mad
  retorno$mbc = mbc
  retorno$mbd = mbd
  retorno$mcd = mcd
  retorno$mn = mn
  return(retorno)

  gc()
}



###############################################################################
# proporção de 1s que ambas as variáveis compartilham nas mesmas posições
# correspondências positivas entre x e y: x e y == 1
###############################################################################
compute.a <- function(x, y){
  return(sum(x == 1 & y == 1))
}

###############################################################################
# proporção de 0s na primeira variável e 1s na segunda variável nas mesmas posições
# x ausente: x == 0 e y == 1
###############################################################################
compute.b <- function(x, y){
  return(sum(x == 0 & y == 1))
}

###############################################################################
# proporção de 1s na primeira variável e 0s na segunda variável nas mesmas posições
# y ausente: x == 1 e y == 0
###############################################################################
compute.c <- function(x, y){
  return(sum(x == 1 & y == 0))
}

###############################################################################
# proporção de zeros que ambas as variáveis compartilham
# correspondências positivas entre x e y: x and y == 0
###############################################################################
compute.d <- function(x,y,m){
  return(sum(x == 0 & y == 0))
}

###############################################################################
# marginal probabilities
# p1 = a + b --> proporção de uns na primeira variável
###############################################################################
compute.ab <- function(a, b){
  return(a+b)
}

###############################################################################
# p2 = a + c --> proporção de uns na segunda variável
compute.ac <- function(a, c){
  return(a+c)
}

###############################################################################
# p3 = a + d --> diagonal (11) (00)
compute.ad <- function(a, d){
  return(a+d)
}

###############################################################################
# p4 = b + c --> diagonal (10) (01)
compute.bc <- function(b, c){
  return(b+c)
}

###############################################################################
# p5 = b + d --> proporção de zeros na segunda variável
compute.bd <- function(b, d){
  return(b+d)
}

###############################################################################
# p6 = c + d --> proporção de zeros na primeira variável
compute.cd <- function(c, d){
  return(c+d)
}

###############################################################################
compute.n <- function(a,b,c,d){
  return(a+b+c+d)
}




###############################################################################
#' Compute contingence table
#'
#' @family MultiLabel Binary Measures Functions
#' @param labels The label space from the dataset
#' @param num.labels total number of labels from label space
#' @return return the contingence table
#'
###############################################################################
compute.cont.table <- function(labels, num.labels){

  retorno = list()

  ma <- build.matrix.corr(num.labels, labels)
  mb <- build.matrix.corr(num.labels, labels)
  mc <- build.matrix.corr(num.labels, labels)
  md <- build.matrix.corr(num.labels, labels)

  u = (num.labels*num.labels)
  pb <- progress_bar$new(total = u)

  #i = 1
  #j = 1

  for (i in 1:num.labels){

    for (j in 1:num.labels){

      x = labels[,i]
      y = labels[,j]

      ma[i,j] = compute.a(x,y)
      mb[i,j] = compute.b(x,y)
      mc[i,j] = compute.c(x,y)
      md[i,j] = compute.d(x,y)

      pb$tick()
      Sys.sleep(1/u)
      gc()
    } # end intern for

    gc()
  } # enf extern for

  retorno$ma = ma
  retorno$mb = mb
  retorno$mc = mc
  retorno$md = md
  return(retorno)

  gc()
}


###############################################################################
#' Build the correlation matrix
#'
#' @family MultiLabel Binary Measures Functions
#' @param labels The label space from the dataset
#' @param num.labels total number of labels from label space
#' @return correlation matrix num.labels X num.labels
#'
###############################################################################
build.matrix.corr <- function(num.labels, labels){
  matrix.corr <- matrix(nrow=num.labels, ncol=num.labels, data=0)
  colnames(matrix.corr) <- colnames(labels)
  rownames(matrix.corr) <- colnames(labels)
  return(matrix.corr)
  gc()
}



###############################################################################
#' Compute all measures for categorial data
#'
#' @family MultiLabel Binary Measures Functions
#' @param labels  the label space from the dataset
#' @param num.labels total number of labels from label space
#' @param  a frequency which label i and label j occur together
#' @param  b frequency which label i occur alone
#' @param  c frequency which label j occur alone
#' @param  d frequency which label i and label j not occur together
#' @param  n n = a + b + c + d
#' @param name the name of the function
#' @param FUN the function
#' @return values computed for all binary measures
#' @examples
#' setwd(Folder)
#' dados = foreign::read.arff(bibtex.arff)
#' labels = data.frame(dados[,1837:1995])
#' names_labels = colnames(labels)
#' result <- compute.measure.2(labels, num.labels, a, b, c, d, n, name, FUN)
#'
###############################################################################
compute.measure.2 <- function(labels, num.labels, a, b, c, d, n, name, FUN){

  retorno = list()

  m <- build.matrix.corr(num.labels, labels)
  u = (num.labels*num.labels) # tamanho da matriz
  pb <- progress_bar$new(total = u) # barra de progresso

  for (i in 1:num.labels){
    for (j in 1:num.labels){
      x = as.numeric(a[i,j]) # a
      y = as.numeric(b[i,j]) # b
      w = as.numeric(c[i,j]) # c
      z = as.numeric(d[i,j]) # d
      k = as.numeric(n[i,j]) # n
      m[i,j] = FUN(x,y,w,z,k) # function
      pb$tick()
      Sys.sleep(1/u)
      gc()
    } # end intern for
    gc()
  } # enf extern for

  return(m)
  gc()
}


#############################################################################
#' Compute Rogers Tanimoto
#'
#' @family MultiLabel Binary Measures Functions
#' @param  a frequency which label i and label j occur together
#' @param  b frequency which label i occur alone
#' @param  c frequency which label j occur alone
#' @param  d frequency which label i and label j not occur together
#' @param  n a + b + c + d
#'
#############################################################################
rogers.tanimoto.e.1 <- function(l){
  d1 = l$a + l$d
  d2 = l$a + (2*(l$b+l$c)) + l$d
  d3 = d1/d2
  return(d3)
}

rogers.tanimoto.e <- function(a,b,c,d,n){
  d1 = a + d
  d2 = a + (2*(b+c)) + d
  d3 = d1/d2
  return(d3)
}


