library(tidyverse)
library(ape)
library(phangorn)
library(seqinr)
library(Biostrings)
library(ShortRead)
library(DECIPHER)
library(phylotools)
library(muscle)
library(tidytree)
library(treeio)
library(ggtree)
library(phylogram)
library(dendextend)


rcc_accession <- read.table("D:/CamiloG/RCCdatabase/accesion_rcc")
rcc_accession <- rcc_accession$V1
as.character(rcc_accession)
rcc_gen <- read.GenBank(rcc_accession, species.names = TRUE)
names_rcc <- data.frame(species = attr(rcc_gen, "species"), accs = names(rcc_gen))
head(names_rcc, n = 10)
names(rcc_gen) <- attr(rcc_gen, "species")
write.dna(rcc_gen, "rcc_seqs.fasta", format = "fasta")

RCC <- read.fasta(file = "rcc_seqs.fasta", clean_name = TRUE)
Neu <- names_rcc
Neu$accs <- paste("_", Neu$accs, sep = "")
Neu <- add_column(Neu, ID = with(Neu, paste0(species, accs)), .after = "species")
ID <- Neu[, 2]
RCC <- add_column(RCC, ID = ID, .after = "seq.name")
RCC <- RCC[, -1]
colnames(RCC) <- c("seq.name", "seq.text")
dat2fasta(RCC, outfile = "seqs_rcc.fasta")

Alignment <- readDNAStringSet("RCC_ASV.fasta", format = "fasta", use.names = TRUE)
Muscle <- muscle(Alignment, gapopen = -400, log = "log.txt", verbose = TRUE, maxhours = 24.0)
fullM <- as.matrix(Muscle)
Alphabet <- alphabetFrequency(Muscle, as.prob = TRUE, collapse = FALSE)
Consensus <- consensusMatrix(Muscle, as.prob = TRUE, baseOnly = TRUE)
write.phylip(Muscle, filepath = "C:/Users/Camilo/Desktop/Muscle.phy")
write.phylip(x = Muscle, filepath = "RCC_M.txt")
sdist <- stringDist(as(Muscle, "DNAStringSet"), method = "hamming")
clust <- hclust(sdist, method = "single")
pdf(file = "fullTree.pdf")
plot(clust)
dev.off()
mDNAStr <- as(Muscle, "DNAStringSet")
writeXStringSet(mDNAStr, file = "RCC_MAlign.fasta", format = "fasta")

fullM <- ape::as.alignment(fullM)
fullM <- ape::as.DNAbin.alignment(fullM)
image.DNAbin(fullM)



autoMasked <- maskGaps(Muscle, min.fraction = 0.5, min.block.width = 4)
partialM <- as.matrix(autoMasked)
write.csv2(partialM, file = "C:/Users/Camilo/Desktop/masked_muscle_alignment.csv")
pAlphabet <- alphabetFrequency(autoMasked)
pConView <- consensusViews(autoMasked)
sdistA <- stringDist(as(autoMasked, "DNAStringSet"), method = "hamming")
clustA <- hclust(sdistA, method = "single")
pdf(file = "MaskedTree.pdf")
plot(clustA)
dev.off()
write.phylip(autoMasked, filepath="RCC_aligned.txt")

objeto <- read.phyDat("RCC_M.txt", format = "interleaved")
mt <- modelTest(objeto)



mDNAStr <- as(Muscle, "DNAStringSet")
writeXStringSet(mDNAStr, file = "RCC_Muscle.fasta")
DMat <- DistanceMatrix(mDNAStr, correction = "Jukes-Cantor", processors = NULL, verbose = TRUE)
ddMuscle <- IdClusters(DMat, method = "ML", showPlot = TRUE, type = "dendrogram", myXStringSet = mDNAStr, processors = 6, verbose = TRUE)
ggplot(ddMuscle, 
       ladderize = TRUE, 
       branch.length = "none") + 
  geom_tree() + theme_tree() + 
  geom_tiplab(geom = "text", 
              size = 4) + 
  geom_treescale(x = 0, 
                 y = 1300, 
                 width = 3, 
                 offset = 1,
                 offset.label = 1,
                 color = "red")

