library(tidyverse)
library(viridis)
library(ape)
library(phangorn)
library(Biostrings)
library(DECIPHER)
library(dendextend)
library(heatmaply)
library(htmlwidgets)
library(gtools)

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


noms <- sASVs@unmasked@ranges@NAMES
ornoms <- mixedsort(noms, decreasing = FALSE)
ornoms
comps <- Shared[, 1]
comps
secuencias <- Shared[, 2]
nuevo <- data.frame(comps, ornoms, secuencias)
nuevo
orig <- nuevo[match(noms, nuevo$ornoms),]
View(orig)
tShaASVs <- as.data.frame(t(ShaASVs))
View(tShaASVs)
tShaASVs <- tShaASVs %>% 
  mutate(ASV = comps, 
         F.Name = ornoms, 
         .before = "C1_2017.08")
tShaASVs <- tShaASVs[match(noms, tShaASVs$F.Name), ]
write.csv2(tShaASVs, file = "tShaASVs.csv")
nShaASVs <- read.csv2("tShaASVs.csv", header = TRUE, sep = ";", dec = ".", skip = 0)
rownames(nShaASVs) <- nShaASVs[, 1]
nShaASVs <- nShaASVs[, -c(1:3)]
cShaASVs <- nShaASVs %>% 
  mutate(Cha = rowSums(nShaASVs[1:12]), 
         Fla = rowSums(nShaASVs[13:24]), 
         Hu = rowSums(nShaASVs[25:36]), 
         Pc = rowSums(nShaASVs[37:41])
         )
cShaASVs <- cShaASVs[, -c(1:41)]

nViridis <- heatmaply(normalize(t(nShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across all sample clustered with maximum likelihood", 
                      margins = c(50, 50, 70, 0), 
                      grid_gap = 1, 
                      width = 1920, 
                      height = 1080, 
                      subplot_heights = c(0.35, 0.65), 
                      color = viridis(n = 256, 
                                      alpha = 1, 
                                      begin = 0, 
                                      end = 1, 
                                      option = "viridis")
                      )
saveWidget(nViridis, file = "viridis_hmshaml.html")

nMagma <- heatmaply(normalize(t(nShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across all samples clustered with maximum likelihood", 
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
saveWidget(nMagma, file = "magma_hmshaml.html")

nPlasma <- heatmaply(normalize(t(nShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across all samples clustered with maximum likelihood", 
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
saveWidget(nPlasma, file = "plasma_hmshaml.html")

nInferno <- heatmaply(normalize(t(nShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across all samples clustered with maximum likelihood", 
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
saveWidget(nInferno, file = "inferno_hmshaml.html")

nCividis <- heatmaply(normalize(t(nShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across all samples clustered with maximum likelihood", 
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
saveWidget(nCividis, file = "cividis_hmshaml.html")

nMako <- heatmaply(normalize(t(nShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across all samples clustered with maximum likelihood", 
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
saveWidget(nMako, file = "mako_hmshaml.html")

nRocket <- heatmaply(normalize(t(nShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across all samples clustered with maximum likelihood", 
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
saveWidget(nRocket, file = "rocket_hmshaml.html")

nHeat <- heatmaply(normalize(t(nShaASVs)), 
                   Colv = c.dend, 
                   Rowv = NA, 
                   main = "Shared ASVs across all samples clustered with maximum likelihood", 
                   margins = c(50, 50, 70, 0), 
                   grid_gap = 1, 
                   grid_color = "#DEE2E6", 
                   width = 1920, 
                   height = 1080, 
                   subplot_heights = c(0.35, 0.65), 
                   color = heat.colors(n = 256, 
                                       alpha = 1,
                                       rev = TRUE)
                   )

saveWidget(nHeat, file = "heat_hmshaml.html")

write.csv2(nShaASVs, file = "nShaASVs.csv")

gViridis <- heatmaply(normalize(t(cShaASVs)), 
                     Colv = c.dend, 
                     Rowv = NA, 
                     main = "Shared ASVs across the four sites clustered with maximum likelihood", 
                     margins = c(50, 50, 70, 0), 
                     grid_gap = 1, 
                     width = 1920, 
                     height = 1080, 
                     subplot_heights = c(0.48, 0.52), 
                     color = viridis(n = 256, 
                                     alpha = 1, 
                                     begin = 0, 
                                     end = 1, 
                                     option = "viridis")
)
saveWidget(gViridis, file = "viridis_grouped_hmshaml.html")

gMagma <- heatmaply(normalize(t(cShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across the four sites clustered with maximum likelihood", 
                      margins = c(50, 50, 70, 0), 
                      grid_gap = 1, 
                      width = 1920, 
                      height = 1080, 
                      subplot_heights = c(0.48, 0.52), 
                      color = viridis(n = 256, 
                                      alpha = 1, 
                                      begin = 0, 
                                      end = 1, 
                                      option = "magma")
)
saveWidget(gMagma, file = "magma_grouped_hmshaml.html")

gPlasma <- heatmaply(normalize(t(cShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across the four sites clustered with maximum likelihood", 
                      margins = c(50, 50, 70, 0), 
                      grid_gap = 1, 
                      width = 1920, 
                      height = 1080, 
                      subplot_heights = c(0.48, 0.52), 
                      color = viridis(n = 256, 
                                      alpha = 1, 
                                      begin = 0, 
                                      end = 1, 
                                      option = "plasma")
)
saveWidget(gPlasma, file = "plasma_grouped_hmshaml.html")

gInferno <- heatmaply(normalize(t(cShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across the four sites clustered with maximum likelihood", 
                      margins = c(50, 50, 70, 0), 
                      grid_gap = 1, 
                      width = 1920, 
                      height = 1080, 
                      subplot_heights = c(0.48, 0.52), 
                      color = viridis(n = 256, 
                                      alpha = 1, 
                                      begin = 0, 
                                      end = 1, 
                                      option = "inferno")
)
saveWidget(gInferno, file = "inferno_grouped_hmshaml.html")

gCividis <- heatmaply(normalize(t(cShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across the four sites clustered with maximum likelihood", 
                      margins = c(50, 50, 70, 0), 
                      grid_gap = 1, 
                      width = 1920, 
                      height = 1080, 
                      subplot_heights = c(0.48, 0.52), 
                      color = viridis(n = 256, 
                                      alpha = 1, 
                                      begin = 0, 
                                      end = 1, 
                                      option = "cividis")
)
saveWidget(gCividis, file = "cividis_grouped_hmshaml.html")

gMako <- heatmaply(normalize(t(cShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across the four sites clustered with maximum likelihood", 
                      margins = c(50, 50, 70, 0), 
                      grid_gap = 1, 
                      width = 1920, 
                      height = 1080, 
                      subplot_heights = c(0.48, 0.52), 
                      color = viridis(n = 256, 
                                      alpha = 1, 
                                      begin = 0, 
                                      end = 1, 
                                      option = "mako")
)
saveWidget(gMako, file = "mako_grouped_hmshaml.html")

gRocket <- heatmaply(normalize(t(cShaASVs)), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Shared ASVs across the four sites clustered with maximum likelihood", 
                      margins = c(50, 50, 70, 0), 
                      grid_gap = 1, 
                      width = 1920, 
                      height = 1080, 
                      subplot_heights = c(0.48, 0.52), 
                      color = viridis(n = 256, 
                                      alpha = 1, 
                                      begin = 0, 
                                      end = 1, 
                                      option = "rocket")
)
saveWidget(gRocket, file = "rocket_grouped_hmshaml.html")

gHeat <- heatmaply(normalize(t(cShaASVs)), 
                   Colv = c.dend, 
                   Rowv = NA, 
                   main = "Shared ASVs across the four sites clustered with maximum likelihood", 
                   margins = c(50, 50, 70, 0), 
                   grid_gap = 1, 
                   grid_color = "#DEE2E6", 
                   width = 1920, 
                   height = 1080, 
                   subplot_heights = c(0.48, 0.52), 
                   color = heat.colors(n = 256, 
                                       alpha = 1,
                                       rev = TRUE)
)
saveWidget(gHeat, file = "heat_grouped_hmshaml.html")






