library(tidyverse)
library(ape)
library(phangorn)
library(seqinr)
library(Biostrings)
library(ShortRead)
library(DECIPHER)
library(phylotools)

rcc_accession <- read.table("D:/CamiloG/RCCdatabase/accesion_rcc")
rcc_accession <- rcc_accession$V1
as.character(rcc_accession)
rcc_gen <- read.GenBank(rcc_accession, species.names = TRUE)
names_rcc <- data.frame(species = attr(rcc_gen, "species"), accs = names(rcc_gen))
head(names_rcc, n = 10)
names(rcc_gen) <- attr(rcc_gen, "species")
write.dna(rcc_gen, "rcc_seqs.fasta", format = "fasta")

Neu <- names_rcc
Neu$accs <- paste("_", Neu$accs, sep = "")
Neu <- add_column(Neu, ID = with(Neu, paste0(species, accs)), .after = "species")
ID <- Neu[, 2]
Actual <- get.fasta.name("/Users/Artemis/Dropbox/R/MEGA/RCC_MUSCLE.fas", clean_name = TRUE)
Muscle <- read.fasta(file = "/Users/Artemis/Dropbox/R/MEGA/RCC_MUSCLE.fas", clean_name = FALSE)
Muscle <- add_column(Muscle, ID = ID, .after = "seq.name")
Muscle <- Muscle[, -1]
colnames(Muscle) <- c("seq.name", "seq.text")
dat2fasta(Muscle, outfile = "/Users/Artemis/Dropbox/R/RCC.fasta")
write.csv2(Muscle, file = "/Users/Artemis/Dropbox/R/RCC_MUSCLE.csv")
RefNames <- tibble(Actual, ID)
rename.fasta(infile = "/Users/Artemis/Dropbox/R/MEGA/RCC_MUSCLE.fas", 
             ref_table = RefNames, 
             outfile = "/Users/Artemis/Dropbox/R/RCC.fasta")
edited <- read.fasta(file = "/Users/Artemis/Dropbox/R/RCC.fasta", clean_name = TRUE)




