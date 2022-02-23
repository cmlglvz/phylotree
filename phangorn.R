library(tidyverse)
library(ape)
library(phangorn)

ASV <- read.FASTA(file = "C:/Users/Camilo/Dropbox/sasvsmmskd.fasta", type = "DNA")
pASV <- phyDat(ASV, type = "DNA", levels = NULL, return.index = TRUE)
#bASV <- adegenet::fasta2DNAbin(file = "C:/Users/Camilo/Dropbox/sasvsmmskd.fasta", quiet = FALSE)

dm <- dist.ml(pASV)
#dnm <- dist.dna(x = bASV, model = "TN93", variance = TRUE)

tNJ <- NJ(dm)
#tnNJ <- NJ(dnm)
#tUP <- upgma(dm)

#fun <- function(x) NJ(dist.ml(x))
#bs_NJ <- bootstrap.phyDat(pASV, fun)
#plotBS(tNJ, bs_NJ, main = "bootstrap NJ")

mt <- modelTest(pASV, tree = tNJ)
bfit <- as.data.frame(mt)
write.csv2(bfit, "C:/Users/Camilo/Desktop/bfit.csv")

fit <- pml(tNJ, data = pASV)
fitG <- optim.pml(fit, TRUE)
fitT <- update(fit, k = 4, inv = 0.26)
fitGTR <- optim.pml(fitT, 
                    model = "GTR", 
                    optNni = TRUE, 
                    optBf = TRUE, 
                    optGamma = TRUE, 
                    rearrangement = "stochastic", 
                    control = pml.control(trace = 0)
                    )
bs <- bootstrap.pml(fitGTR, bs = 1000, optNni = TRUE, control = pml.control(trace = 0))
plotBS(midpoint(fitGTR$tree), bs, p = 50, type = "p")


