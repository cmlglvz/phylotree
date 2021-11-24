library(tidyverse)
library(ape)
library(phangorn)
library(ips)
library(Biostrings)
library(DECIPHER)

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