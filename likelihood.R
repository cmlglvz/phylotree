library(tidyverse)
library(ape)
library(ade4)
library(adegenet)
library(phangorn)
library(Biostrings)
library(DECIPHER)

bingo <- fasta2DNAbin(file = "pmskrd.fasta", quiet = FALSE)
dna <- as.phyDat(bingo)
dm <- dist.ml(dna, "JC69", k = 4)
treeUPGMA <- upgma(dm)
treeNJ <- phangorn::NJ(dm)
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "unrooted", main = "NJ")
fun <- function(x) upgma(dist.ml(x))
bs_upgma <- bootstrap.phyDat(dna, fun)
plotBS(treeUPGMA, bs_upgma, main="UPGMA")
nfun <- function(x) phangorn::NJ(dist.ml(x))
bs_nj <- bootstrap.phyDat(dna, nfun)
plotBS(treeNJ, bs_nj, main = "NJ")
fit <- pml(treeNJ, data = dna, k = 4)
fitJC <- optim.pml(fit, TRUE)
logLik(fitJC)
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "NNI", control = pml.control(trace = 0))
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR
mt <- modelTest(dna)



