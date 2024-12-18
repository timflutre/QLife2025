setwd("~/Desktop/qlife_2025/simu")

data <- read.table("test.vcf", comment.char = "", skip = 5, header = TRUE)

data <- read.table("pedigree/pedigree2.txt", col.names = c("ind", "sire", "dam", "sex", "aff", "gen", "pheno", "BV"))
xmean <- aggregate(data[c("pheno", "BV")], by = list(g = data$gen), mean)

plot(xmean$g, xmean$BV)

data <- read.table("freq/freq1.txt")
colnames(data) <- c("gen", "N", paste0("L", 1:(ncol(data)-2)))

hist(unlist(data[11,-c(1:2)]))



setwd("~/Desktop/qlife_2025/20241210_simu/genotype")

data <- read.table("origins_2.txt")



#########################################"
############################################"""""
## the process to simulate correlated traits is taken from AlphaSimR
nind <- 100
nqtl <- 10
ntrait <- 3
corrA <- matrix(c(1,-1, 0.5,-1,1,0.25,0.5,0.25,1), nrow = ntrait)

beta <- matrix(rnorm(nqtl*ntrait), ncol = ntrait) %*% transMat(corrA)

transMat <- function(corMat) {
    ## corMat = correlation matrix between the traits
    ## needs to be positive semi-definite
    eig = eigen(corMat, symmetric=TRUE)
    return(t(eig$vectors %*% (t(eig$vectors)*sqrt(pmax(eig$values, 0)))))
}

## this also gives me correlated traits, with the right correlation coefficients
## question: can the correlation change with time ?
## possibility of not having all the QTL shared between the different traits? 
## -> probably not... when doing the multiplication by transMat(corMat) the zeros disappear
geno <- matrix(sample(0:1, nqtl*nind, replace = TRUE), nrow = nind)

bv <- geno %*% beta



