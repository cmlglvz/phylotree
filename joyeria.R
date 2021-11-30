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
writeXStringSet(eMasked, filepath = "eMasked.fasta", format = "fasta")

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

alist <- c("Fragilaria_sp-_RCC5541_")
eMasked[alist]
emksd <- eMasked[-c(2, 3, 5, 8, 10, 13:14, 16:18, 21, 23, 32, 34, 40:41, 47, 52:53, 56, 59, 62, 64, 68, 71, 76, 77, 80, 83, 85, 87, 
                    96, 98, 105, 109, 111, 113, 116, 117, 122, 129, 131:135, 143, 147, 150:152, 154, 156, 158, 161:164, 168:170, 
                    172:173, 175, 184, 185, 191, 197, 202, 207, 212, 215, 217, 220:221, 226, 232, 238, 243, 246:247, 251, 254, 256, 
                    258:259, 266, 268:269, 273:276, 278, 281, 283, 286, 290:292,306, 308:309, 311:312, 314:315, 317:318, 321, 322, 
                    325, 328, 329, 335, 339:341, 344, 350, 353, 362, 368, 370:371, 375, 377)]
firstep <- eDNAStr[-c(7, 11, 15, 19, 20, 25, 26, 31, 40, 43, 52, 58, 62, 64, 67, 75, 79, 82, 88, 95, 99, 100, 111, 113, 121, 123, 
                      125, 131, 132, 135, 136, 138, 142, 147, 150:153, 160, 165, 168, 173, 175, 176, 179, 181, 187, 193, 195, 197, 198, 
                      204, 208, 215:217, 221, 227, 229, 230, 236:238, 240, 245, 248, 257, 258, 260, 262, 264, 268, 274, 276, 278, 299, 
                      300, 307:310, 313, 321, 322, 325, 331, 333, 334, 340, 342:344, 357, 361:365, 367, 371, 380, 385, 407, 414:416, 
                      425, 426, 432, 434:436, 439, 440, 448, 449, 454, 461, 463, 470, 476, 479, 491, 494, 496, 503), ]
BrowseSeqs(firstep, htmlFile = "firstep.html")
secondstep <- firstep[-c(2, 3, 5, 8, 10, 13:14, 16:18, 21, 23, 32, 34, 40:41, 47, 52:53, 56, 59, 62, 64, 68, 71, 76, 77, 80, 83, 85, 87, 
                         96, 98, 105, 109, 111, 113, 116, 117, 122, 129, 131:135, 143, 147, 150:152, 154, 156, 158, 161:164, 168:170, 
                         172:173, 175, 184, 185, 191, 197, 202, 207, 212, 215, 217, 220:221, 226, 232, 238, 243, 246:247, 251, 254, 256, 
                         258:259, 266, 268:269, 273:276, 278, 281, 283, 286, 290:292,306, 308:309, 311:312, 314:315, 317:318, 321, 322, 
                         325, 328, 329, 335, 339:341, 344, 350, 353, 362, 368, 370:371, 375, 377), ]
BrowseSeqs(secondstep, htmlFile = "secondstep.html")
writeXStringSet(secondstep, filepath = "aleRAs_edited.afa", format = "fasta")
custom <- readDNAMultipleAlignment(filepath = "aleRAs_edited.afa", format = "fasta")
mascara <- maskGaps(custom, min.fraction = 0.3, min.block.width = 4)
cMSKR <- as(mascara, "DNAStringSet")
writeXStringSet(cMSKR, filepath = "cmskr.fa", format = "fasta")
BrowseSeqs(cMSKR, htmlFile = "cmskr.html")

clsdist <- stringDist(cMSKR, method = "levenshtein")
chsdist <- stringDist(cMSKR, method = "hamming")
cls.clust <- hclust(clsdist, method = "single")
clc.clust <- hclust(clsdist, method = "complete")
chs.clust <- hclust(chsdist, method = "single")
chc.clust <- hclust(chsdist, method = "complete")

plot(cls.clust, hang = -1, cex = 0.7)
plot(clc.clust, hang = -1, cex = 0.7)
plot(chs.clust, hang = -1, cex = 0.7)
plot(chc.clust, hang = -1, cex = 0.7)

eras <- readDNAStringSet(filepath = "eRAs.fasta", format = "fasta", skip = 0, use.names = TRUE)
lceras <- firstep@ranges@NAMES
cRAs <- eras[lceras]
BrowseSeqs(cRAs)
lfceras <- secondstep@ranges@NAMES
cRAs <- cRAs[lfceras]
BrowseSeqs(cRAs, htmlFile = "cras.html")
writeXStringSet(x = cRAs, filepath = "cRAs.fas", format = "fasta")
#Se alineo con MUSCLEv5

praxis <- readDNAMultipleAlignment(filepath = "cRAs.afa", format = "fasta")
maskered <- maskGaps(praxis, min.fraction = 0.3, min.block.width = 4)
pMSKRD <- as(maskered, "DNAStringSet")
writeXStringSet(pMSKRD, filepath = "pmskrd.fas", format = "fasta")
BrowseSeqs(pMSKRD, htmlFile = "pmskrd.html")

plsdist <- stringDist(pMSKRD, method = "levenshtein")
phsdist <- stringDist(pMSKRD, method = "hamming")
pls.clust <- hclust(plsdist, method = "single")
plc.clust <- hclust(plsdist, method = "complete")
phs.clust <- hclust(phsdist, method = "single")
phc.clust <- hclust(phsdist, method = "complete")

plot(pls.clust, hang = -1, cex = 0.7)
plot(plc.clust, hang = -1, cex = 0.7)
plot(phs.clust, hang = -1, cex = 0.7)
plot(phc.clust, hang = -1, cex = 0.7)

Ali <- "pmskrd.fasta"
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(Ali, type = "FASTA", dbFile = dbConn, "")

# identify the sequences by their species
x <- dbGetQuery(dbConn, "select description from Seqs")$description
Add2DB(myData = data.frame(identifier = x, stringsAsFactors = FALSE), dbConn)

# form a consensus for each species
cons <- IdConsensus(dbConn, threshold = 0.3, minInformation = 0.1)

# calculate a maximum likelihood tree
d <- DistanceMatrix(cons, correction = "Jukes-Cantor", processors = NULL, verbose = TRUE)
dend <- IdClusters(d, 
                   method = "ML", 
                   showPlot = TRUE, 
                   type = "dendrogram", 
                   myXStringSet = cons, 
                   processors = NULL, 
                   verbose = TRUE)

hclust <- as.hclust(dend)
plot(hclust, cex = 0.7)
devolve <- as.dendrogram(hclust)
plot(devolve)

dbDisconnect(dbConn)





