library(tidyverse)
library(viridis)
library(ape)
library(phangorn)
library(ips)
library(Biostrings)
library(DECIPHER)
library(dendextend)
library(heatmaply)
library(htmlwidgets)
library(gtools)
library(phyloseq)

sASVs <- readDNAMultipleAlignment("sASVs.afa", format = "fasta")
sDNAStr <- as(sASVs, "DNAStringSet")
BrowseSeqs(sDNAStr, htmlFile = "sasvs.html")
writeXStringSet(sDNAStr, filepath = "sasvs.fas", format = "fasta")
autoMasked <- maskGaps(sASVs, min.fraction = 0.3, min.block.width = 4)
ATMskd <- as(autoMasked, "DNAStringSet")
writeXStringSet(ATMskd, filepath = "sasvsmmskd.fas", format = "fasta")
BrowseSeqs(ATMskd, htmlFile = "atmskdsasvs.html")
alfabeto <- alphabetFrequency(autoMasked)
ConMat <- consensusMatrix(autoMasked, as.prob = TRUE, baseOnly = TRUE)[, 1:158]
nidea <- substr(consensusString(autoMasked), 1, 418)
consensusViews(autoMasked)
lsdist <- stringDist(ATMskd, method = "levenshtein")
hsdist <- stringDist(ATMskd, method = "hamming")
ls.clust <- hclust(lsdist, method = "single")
lc.clust <- hclust(lsdist, method = "complete")
hs.clust <- hclust(hsdist, method = "single")
hc.clust <- hclust(hsdist, method = "complete")
msasvs <- as.matrix(sASVs)
mOTU <- rownames(msasvs)

plot(ls.clust, hang = -1, cex = 0.7)
plot(lc.clust, hang = -1, cex = 0.7)
plot(hs.clust, hang = -1, cex = 0.7)
plot(hc.clust, hang = -1, cex = 0.7)

col_dend <- hc.clust %>% as.dendrogram() %>% set("branches_lwd", 0.3) %>% ladderize()
oASVs <- read.csv2("/Users/Artemis/Dropbox/R/eAnalisis/oASVs.csv", 
                     sep = ";", dec = ".", header = TRUE, skip = 0)
rownames(oASVs) <- oASVs[, 1]
oASVs <- oASVs[, -1]
APwATs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/APwATs.csv", 
                    sep = ";", dec = ".", header = TRUE, skip = 0)
Shared <- APwATs %>% filter(Cha == 1 & Fla == 1 & Hu == 1 & Pc == 1)
ShaSeq <- Shared[, 2]
ShaOTU <- Shared[, 1]
ShaASVs <- select(oASVs, all_of(ShaSeq))
dOTU <- mixedsort(mOTU, decreasing = FALSE)
xTXs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/xTXs.csv", 
                  sep = ";", dec = ".", header = TRUE, skip = 0)
ShaTXs <- xTXs %>% filter(Seq %in% all_of(ShaSeq))
nombres <- data.frame(ShaOTU, mOTU, dOTU)
colnames(ShaASVs) <- dOTU
eShaASVs <- ShaASVs[, mOTU]

hmsasvs<- heatmaply(normalize(eShaASVs), 
                   Colv = col_dend, 
                   Rowv = NA, 
                   main = "Shared ASVs across all sites ordered by clustering alignments based on their distance to each other",
                   margins = c(50, 50, 70, 0),
                   grid_gap = 1, 
                   width = 1920,
                   height = 1080, 
                   subplot_heights = c(0.35, 0.65)
                   )
saveWidget(hmsasvs, file = "hmsasvs.html")

shmsasvs<- heatmaply(percentize(eShaASVs), 
                     Colv = col_dend, 
                     Rowv = NA, 
                     main = "Shared ASVs across all sites ordered by clustering alignments based on their distance to each other", 
                     margins = c(50, 50, 70, 0), 
                     grid_gap = 1, 
                     width = 1920, 
                     height = 1080,
                     subplot_heights = c(0.35, 0.64)
                     )
saveWidget(shmsasvs, file = "shmsasvs.html")

gShaASVs <- as.data.frame(t(eShaASVs))
gShaASVs <- gShaASVs %>% mutate(Cha = rowSums(gShaASVs[1:12]), 
                    Fla = rowSums(gShaASVs[13:24]), 
                    Hu = rowSums(gShaASVs[25:36]),
                    Pc = rowSums(gShaASVs[37:41]), 
                    .before = "C1_2017.08")
gShaASVs <- gShaASVs[, -c(5:45)]

ghmsasvs<- heatmaply(percentize(t(gShaASVs)), 
                     Colv = col_dend, 
                     Rowv = NA, 
                     main = "Shared ASVs across all sites ordered by clustering alignments based on their distance to each other", 
                     margins = c(50, 50, 70, 0), 
                     grid_gap = 1, 
                     width = 1920, 
                     height = 1080,
                     subplot_heights = c(0.5, 0.5)
                     )
saveWidget(ghmsasvs, file = "ghmsasvs.html")

hmgsha <- heatmaply(normalize(t(gShaASVs)), 
                    Colv = col_dend, 
                    Rowv = NA, 
                    main = "Shared ASVs across all sites ordered by clustering alignments based on their distance to each other",
                    margins = c(50, 50, 70, 0),
                    grid_gap = 1, 
                    width = 1920,
                    height = 1080,
                    subplot_heights = c(0.5, 0.5)
                    )
saveWidget(hmgsha, file = "hmgsha.html")

#Mejoremos el dendrograma
Ali <- "sasvs.fas"
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
plot(hclust, hang = -1, cex = 0.7)
devolve <- as.dendrogram(hclust)
plot(devolve)

dbDisconnect(dbConn)

c.dend <- devolve %>% set("branches_lwd", 0.3) %>% ladderize()

HM.MLsASVs<- heatmaply(normalize(eShaASVs), 
                    Colv = c.dend, 
                    Rowv = NA, 
                    main = "Shared ASVs across all sites ordered by clustering alignments based on their maximum likelihood",
                    margins = c(50, 50, 70, 0),
                    grid_gap = 1, 
                    width = 1920,
                    height = 1080, 
                    subplot_heights = c(0.35, 0.65)
)
saveWidget(HM.MLsASVs, file = "ml_hmsasvs.html")

plasma <- heatmaply(normalize(eShaASVs), 
                       Colv = c.dend, 
                       Rowv = NA, 
                       main = "Shared ASVs across all sites ordered by clustering alignments based on their maximum likelihood",
                       margins = c(50, 50, 70, 0),
                       grid_gap = 1, 
                       width = 1920,
                       height = 1080, 
                       subplot_heights = c(0.35, 0.65), 
                    color = viridis(n = 256, 
                                    alpha = 1, 
                                    begin = 0, 
                                    end = 1, 
                                    option = "plasma")
)
saveWidget(plasma, file = "plasma_hmsasvs.html")

magma <- heatmaply(normalize(eShaASVs), 
                    Colv = c.dend, 
                    Rowv = NA, 
                    main = "Shared ASVs across all sites ordered by clustering alignments based on their maximum likelihood",
                    margins = c(50, 50, 70, 0),
                    grid_gap = 1, 
                    width = 1920,
                    height = 1080, 
                    subplot_heights = c(0.35, 0.65), 
                    color = viridis(n = 256, 
                                    alpha = 1, 
                                    begin = 0, 
                                    end = 1, 
                                    option = "magma")
)
saveWidget(magma, file = "magma_hmsasvs.html")

inferno <- heatmaply(normalize(eShaASVs), 
                   Colv = c.dend, 
                   Rowv = NA, 
                   main = "Shared ASVs across all sites ordered by clustering alignments based on their maximum likelihood",
                   margins = c(50, 50, 70, 0),
                   grid_gap = 1, 
                   width = 1920,
                   height = 1080, 
                   subplot_heights = c(0.35, 0.65), 
                   color = viridis(n = 256, 
                                   alpha = 1, 
                                   begin = 0, 
                                   end = 1, 
                                   option = "inferno")
)
saveWidget(inferno, file = "inferno_hmsasvs.html")

cividis <- heatmaply(normalize(eShaASVs), 
                   Colv = c.dend, 
                   Rowv = NA, 
                   main = "Shared ASVs across all sites ordered by clustering alignments based on their maximum likelihood",
                   margins = c(50, 50, 70, 0),
                   grid_gap = 1, 
                   width = 1920,
                   height = 1080, 
                   subplot_heights = c(0.35, 0.65), 
                   color = viridis(n = 256, 
                                   alpha = 1, 
                                   begin = 0, 
                                   end = 1, 
                                   option = "cividis")
)
saveWidget(cividis, file = "cividis_hmsasvs.html")

mako <- heatmaply(normalize(eShaASVs), 
                   Colv = c.dend, 
                   Rowv = NA, 
                   main = "Shared ASVs across all sites ordered by clustering alignments based on their maximum likelihood",
                   margins = c(50, 50, 70, 0),
                   grid_gap = 1, 
                   width = 1920,
                   height = 1080, 
                   subplot_heights = c(0.35, 0.65), 
                   color = viridis(n = 256, 
                                   alpha = 1, 
                                   begin = 0, 
                                   end = 1, 
                                   option = "mako")
)
saveWidget(mako, file = "mako_hmsasvs.html")

rocket <- heatmaply(normalize(eShaASVs), 
                   Colv = c.dend, 
                   Rowv = NA, 
                   main = "Shared ASVs across all sites ordered by clustering alignments based on their maximum likelihood",
                   margins = c(50, 50, 70, 0),
                   grid_gap = 1, 
                   width = 1920,
                   height = 1080, 
                   subplot_heights = c(0.35, 0.65), 
                   color = viridis(n = 256, 
                                   alpha = 1, 
                                   begin = 0, 
                                   end = 1, 
                                   option = "rocket")
)
saveWidget(rocket, file = "rocket_hmsasvs.html")

classic <- heatmaply(normalize(eShaASVs), 
                    Colv = c.dend, 
                    Rowv = NA, 
                    main = "Shared ASVs across all sites ordered by clustering alignments based on their maximum likelihood",
                    margins = c(50, 50, 70, 0),
                    grid_gap = 1, 
                    grid_color = "#dee2e6",
                    width = 1920,
                    height = 1080, 
                    subplot_heights = c(0.35, 0.65), 
                    color = heat.colors(n = 256, 
                                        alpha = 1, 
                                        rev = TRUE)
                    )
saveWidget(classic, file = "classic_hmsasvs.html")

HM.MLgSha <- heatmaply(normalize(t(gShaASVs)), 
                    Colv = c.dend, 
                    Rowv = NA, 
                    main = "Shared ASVs across all sites ordered by clustering alignments based on their maximum likelihood",
                    margins = c(50, 50, 70, 0),
                    grid_gap = 1, 
                    width = 1920,
                    height = 1080,
                    subplot_heights = c(0.5, 0.5), 
                    color = viridis(n = 256, 
                                    alpha = 1, 
                                    begin = 0, 
                                    end = 1, 
                                    option = "turbo")
)
saveWidget(HM.MLgSha, file = "turbo_hmgsha.html")

oHM.MLgSha <- heatmaply(normalize(t(gShaASVs)), 
                       Colv = c.dend, 
                       Rowv = NA, 
                       main = "Shared ASVs across all sites ordered by clustering alignments based on their maximum likelihood",
                       margins = c(50, 50, 70, 0),
                       grid_gap = 1, 
                       grid_color = "#dee2e6",
                       width = 1920,
                       height = 1080,
                       subplot_heights = c(0.5, 0.5), 
                       color = heat.colors(n = 256, 
                                           alpha = 1, 
                                           rev = TRUE)
)
saveWidget(oHM.MLgSha, file = "heat_hmgsha.html")

ShaASVs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/ShaASVs.csv", header = TRUE, sep = ";", dec = ".", row.names = 1, skip = 0)
asvtable <- as.matrix(ShaASVs)
ShaTXs <- read.csv2(file = "https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/ShaTXs.csv", header = TRUE, sep = ";", dec = ".", row.names = 1, skip = 0)
cnms <- sASVs@unmasked@ranges@NAMES
colnames(asvtable) <- ShaTXs$OTU
taxmat <- as.matrix(ShaTXs)
rownames(taxmat) <- ShaTXs$OTU
taxmat <- taxmat[, -c(1,2)]
OTU <- otu_table(object = asvtable, taxa_are_rows = FALSE)
TAX <- tax_table(object = taxmat)
physeq <- phyloseq(OTU, TAX)
plot_bar(physeq, fill = "Genus")













