library(tidyverse)
library(ape)
library(phangorn)
library(ips)
library(Biostrings)
library(DECIPHER)

edited <- readDNAMultipleAlignment("aleRAs.afa", format = "fasta")
eDNAStr <- as(edited, "DNAStringSet")
BrowseSeqs(eDNAStr, htmlFile = "aleras.html")
writeXStringSet(eDNAStr, filepath = "eAlRA.fa", format = "fasta")
autoMasked <- maskGaps(edited, min.fraction = 0.3, min.block.width = 4)
ATMskd <- as(autoMasked, "DNAStringSet")
writeXStringSet(ATMskd, filepath = "atmmskd.fa", format = "fasta")
BrowseSeqs(ATMskd, htmlFile = "atmskd.html")
alfabeto <- alphabetFrequency(autoMasked)
ConMat <- consensusMatrix(autoMasked, as.prob = TRUE, baseOnly = TRUE)[, 1:504]
nidea <- substr(consensusString(autoMasked), 219, 431)
consensusViews(autoMasked)
lsdist <- stringDist(ATMskd, method = "levenshtein")
hsdist <- stringDist(ATMskd, method = "hamming")
ls.clust <- hclust(lsdist, method = "single")
lc.clust <- hclust(lsdist, method = "complete")
hs.clust <- hclust(hsdist, method = "single")
hc.clust <- hclust(hsdist, method = "complete")

par(mar = c(5, 4, 4, 1.8) + 0.1)
plot(ls.clust, hang = -1, cex = 0.7, main = "Levenshtein distance with single hierarchical clustering")
dev.off()

dhc <- as.dendrogram(ls.clust)

par(mar = c(1, 1, 1, 7) + 0.1)
plot(dhc, 
     type = "rectangle", 
     main = "Levenshtein distance with single hierarchical clustering", 
     xlab = "height", 
     ylim = c(1, 504),
     horiz = TRUE)
dev.off()

eMasked <- ATMskd[-c(7, 11, 15, 19, 20, 25, 26, 31, 40, 43, 52, 58, 62, 64, 67, 75, 79, 82, 88, 95, 99, 100, 111, 113, 121, 123, 
                     125, 131, 132, 135, 136, 138, 142, 147, 150:153, 160, 165, 168, 173, 175, 176, 179, 181, 187, 193, 195, 197, 198, 
                     204, 208, 215:217, 221, 227, 229, 230, 236:238, 240, 245, 248, 257, 258, 260, 262, 264, 268, 274, 276, 278, 299, 
                     300, 307:310, 313, 321, 322, 325, 331, 333, 334, 340, 342:344, 357, 361:365, 367, 371, 380, 385, 407, 414:416, 
                     425, 426, 432, 434:436, 439, 440, 448, 449, 454, 461, 463, 470, 476, 479, 491, 494, 496, 503), ]
BrowseSeqs(eMasked, htmlFile = "emasked.html")
writeXStringSet(eMasked, filepath = "eMasked.fas", format = "fasta")

elsdist <- stringDist(ATMskd, method = "levenshtein")
ehsdist <- stringDist(ATMskd, method = "hamming")
els.clust <- hclust(elsdist, method = "single")
elc.clust <- hclust(elsdist, method = "complete")
ehs.clust <- hclust(ehsdist, method = "single")
ehc.clust <- hclust(ehsdist, method = "complete")

plot(els.clust, hang = -1, cex = 0.7)
plot(elc.clust, hang = -1, cex = 0.7)
plot(ehs.clust, hang = -1, cex = 0.7)
plot(ehc.clust, hang = -1, cex = 0.7)

















fast <- readDNAMultipleAlignment("eRAs.afa", format = "fasta")
fDNAStr <- as(fast, "DNAStringSet")
writeXStringSet(fDNAStr, filepath = "fRAMxd.fa", format = "fasta")
fAutoMasked <- maskGaps(fast, min.fraction = 0.3, min.block.width = 4)
fATMasked <- as(fAutoMasked, "DNAStringSet")
writeXStringSet(fATMasked, "fATMasked.fa", format = "fasta")
BrowseSeqs(fATMasked, htmlFile = "fatmasked.html")
BrowseSeqs(fDNAStr, htmlFile = "fdnastr.html")
defMasked <- maskGaps(fast, min.fraction = 0.5, min.block.width = 4)
dMasked <- as(defMasked, "DNAStringSet")
BrowseSeqs(dMasked, htmlFile = "dmasked.html")


muscl <- readDNAMultipleAlignment("SrccEasv.afa", format = "fasta")
xmuscl <- readDNAMultipleAlignment("RCC_ASV.afa", format = "fasta")
DNAStr <- as(muscl, "DNAStringSet")
xDNAStr <- as(xmuscl, "DNAStringSet")
writeXStringSet(DNAStr, filepath = "RAMxd.fa", format = "fasta")
writeXStringSet(xDNAStr, filepath = "xRAMxd.fa", format = "fasta")
autoMasked <- maskGaps(muscl, min.fraction = 0.5, min.block.width = 4) #default settings
xautoMasked <- maskGaps(xmuscl, min.fraction = 0.4, min.block.width = 4)
ATMasked <- as(autoMasked, "DNAStringSet")
xATMasked <- as(xautoMasked, "DNAStringSet")
BrowseSeqs(ATMasked, htmlFile = "automasked.html")
BrowseSeqs(xATMasked, htmlFile = "xautomasked.html")
writeXStringSet(xATMasked, filepath = "xATMasked.fa", format = "fasta")
MaskEx <- MaskAlignment(myXStringSet = DNAStr, 
                        type = "sequences", 
                        windowSize = 5, 
                        threshold = 1, 
                        maxFractionGaps = 0.2, 
                        includeTerminalGaps = FALSE, 
                        correction = FALSE, 
                        showPlot = TRUE)
vMaskEx <- MaskAlignment(myXStringSet = DNAStr, 
                        type = "sequences", 
                        windowSize = 5, 
                        threshold = 1, 
                        maxFractionGaps = 0.2, 
                        includeTerminalGaps = FALSE, 
                        correction = TRUE, 
                        showPlot = TRUE)
#Display only masked nucleotides covered by the mask
masked <- MaskEx
vmasked <- vMaskEx
colmask(masked, append = "replace", invert = TRUE) <- colmask(masked)
colmask(vmasked, append = "replace", invert = TRUE) <- colmask(vmasked)
masked <- as(masked, "DNAStringSet")
vmasked <- as(vmasked, "DNAStringSet")
BrowseSeqs(masked, htmlFile = "masked_first.html")
BrowseSeqs(vmasked, htmlFile = "vmasked.html")