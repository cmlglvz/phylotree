library(tidyverse)
library(ape)
library(seqinr)

rcc_accession <- read.table("D:/CamiloG/RCCdatabase/accesion_rcc")
rcc_accession <- rcc_accession$V1
as.character(rcc_accession)
rcc_gen <- read.GenBank(rcc_accession, species.names = TRUE)
names_rcc <- data.frame(species = attr(rcc_gen, "species"), accs = names(rcc_gen))
head(names_rcc, n = 10)
names(rcc_gen) <- attr(rcc_gen, "species")
write.dna(rcc_gen, "rcc_seqs.fasta", format = "fasta")
