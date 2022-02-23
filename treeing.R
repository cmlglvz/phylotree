library(tidyverse)
library(ape)
library(ade4)
library(adegenet)
library(phangorn)
library(Biostrings)
library(DECIPHER)
library(phytools)
library(phylogram)
library(viridis)
library(dendextend)
library(heatmaply)
library(htmlwidgets)
library(gtools)

ogtree <- read.tree("pmskrd_ogtree.nwk", skip = 0)
plot(ogtree, no.margin = TRUE)
bstree <- read.tree("pmskrd_bootstrapctree.nwk", skip = 0)
plot(bstree, no.margin = TRUE, main = "Bootstrap 500")

vtree <- phytools::read.newick("pmskrd_ogtree.nwk")
plot(vtree, no.margin = TRUE, main = "Original Tree by MEGA7")
vtoh <- phylogram::as.dendrogram(vtree)
c.dend <- vtoh %>% set("branches_lwd", 0.3) %>% ladderize()

joya <- readDNAMultipleAlignment("pmskrd.fasta", format = "fasta")
joyeria <- as(joya, "DNAStringSet")
lista <- joyeria@ranges@NAMES

oASVs <- read.csv2("/Users/Artemis/Dropbox/R/eAnalisis/oASVs.csv", 
                   sep = ";", dec = ".", header = TRUE, skip = 0)
rownames(oASVs) <- oASVs[, 1]
oASVs <- oASVs[, -1]