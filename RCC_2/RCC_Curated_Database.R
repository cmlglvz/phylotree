library(tidyverse) #Si
library(ape) #Si
library(phangorn) #Si
library(seqinr) #Si
library(Biostrings) #Si
library(ShortRead) #Si
library(DECIPHER) #Si
library(phylotools) #Si
library(muscle)
library(tidytree)
library(treeio)
library(ggtree)
library(phylogram)
library(dendextend)

curated_rcc <- read.csv2("RCC_2/curated_rcc.txt", 
                         header = TRUE, 
                         sep = ";", 
                         dec = ".", 
                         skip = 0)
write.csv2(curated_rcc, "RCC_2/curated_rcc.csv")
selected_rcc <- filter(curated_rcc, Keep == "Yes")
selected_rcc <- selected_rcc[-c(1110,1111),]
no_rcc <- filter(curated_rcc, Keep == "No")
rcc_accession <- read.table("RCC_2/accesion_rcc.txt")
rcc_accession <- rcc_accession$V1
as.character(rcc_accession)
print(all(rcc_accession%in%selected_rcc$Genbank_accession))
base::setdiff(rcc_accession,selected_rcc$Genbank_accession)

#Get sequences
rcc_gen <- read.GenBank(rcc_accession, species.names = TRUE)
#Another option
rcc_gen <- read.GenBank(selected_rcc$Genbank_accession, species.names = TRUE)
write.dna(rcc_gen, "RCC_2/og_rcc_gen.fasta", format = "fasta")
#Create a dataframe with sequences
names_rcc <- data.frame(species = attr(rcc_gen, "species"), Genbank_accession = names(rcc_gen)) #note that you get names from NCBI database and not from RCC
head(names_rcc, n = 10)
#We add "species" names from NCBI to downloaded sequences that are in a DNAbin object
ncbi <- rcc_gen
names(ncbi) <- attr(ncbi, "species")
#Export the sequences with its names to a fasta file
write.dna(ncbi, "RCC_2/ncbi_seqs.fasta", format = "fasta")
#Another option is to add the RCC code instead of NCBI names
rcc_codes <- data.frame(names_rcc$species,
                        paste(selected_rcc$RCC, selected_rcc$Genbank_accession, sep = "_")
                        )
write.csv2(rcc_codes, "RCC_2/rcc_codes.csv")
#I externally edited the file
rcc_codes_edited <- read.csv2("RCC_2/rcc_codes.csv", header = TRUE, sep = ";", dec = ".", skip = 0, row.names = 1)
#you can do it as you please
rename.fasta(infile = "RCC_2/ncbi_seqs.fasta", ref_table = rcc_codes, outfile = "RCC_2/rcc_seqs.fasta")

#Bit easier
names(rcc_gen) <- paste(rcc_codes_edited$names, rcc_codes_edited$codes, sep = "_")
write.dna(rcc_gen, "RCC_2/rcc_id_seqs.fasta", format = "fasta")
check <- read.fasta("RCC_2/rcc_id_seqs.fasta", clean_name = TRUE)
dat2fasta(check, outfile = "RCC_2/checked.fasta")
unified <- read.fasta("RCC_2/unified.fasta")
og_unified <- "RCC_2/unified.fasta"
unidos <- readDNAStringSet(og_unified)
unimil <- subseq(unidos, width = c(1:1000))
uni_align <- readDNAMultipleAlignment("RCC_2/uni_aligned.afa", format = "fasta")
uni <- as(uni_align, "DNAStringSet")
BrowseSeqs(uni, htmlFile = "RCC_2/multiple_align.html")
writeXStringSet(uni, filepath = "RCC_2/multi_uni.fas", format = "fasta")
autoMasked <- maskGaps(uni_align, min.fraction = 0.3, min.block.width = 4)
ATMskd <- as(autoMasked, "DNAStringSet")
writeXStringSet(ATMskd, filepath = "RCC_2/multi_masked.fasta", format = "fasta")
BrowseSeqs(ATMskd, htmlFile = "RCC_2/auto_masked.html")
########################################################################################################################################################

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

