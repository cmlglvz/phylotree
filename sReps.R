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

cASVs <- readDNAMultipleAlignment("sASVs_reps.afa", format = "fasta")
cDNAStr <- as(cASVs, "DNAStringSet")
BrowseSeqs(cDNAStr, htmlFile = "sasvs_reps.html")
writeXStringSet(cDNAStr, filepath = "sasvs_reps.fas", format = "fasta")
autoMasked <- maskGaps(cASVs, min.fraction = 0.3, min.block.width = 4)
cATMskd <- as(autoMasked, "DNAStringSet")
writeXStringSet(cATMskd, filepath = "sasvsmmskd_reps.fas", format = "fasta")
BrowseSeqs(cATMskd, htmlFile = "atmskd_sasvs_reps.html")
mcsasvs <- as.matrix(cASVs)
mcOTU <- rownames(mcsasvs)
dcOTU <- mixedsort(mcOTU, decreasing = FALSE)

Ali <- "sasvs_reps.fas"
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(Ali, type = "FASTA", dbFile = dbConn, "")

x <- dbGetQuery(dbConn, "select description from Seqs")$description
Add2DB(myData = data.frame(identifier = x, stringsAsFactors = FALSE), dbConn)

cons <- IdConsensus(dbConn, threshold = 0.3, minInformation = 0.1)

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
cShaASVs <- ShaASVs[, -c(1, 6, 9, 10, 16, 17, 18, 21, 26, 28:29, 34:35, 43, 45, 48:50, 54, 57, 60, 62:63, 66, 69:71, 73:74, 77:79, 
                         81, 83, 86, 92:93, 96, 98, 103, 106:108, 110, 112, 114:115, 120:122, 125, 127, 129, 131:134, 138, 142:145, 
                         148)]
cShared <- Shared[-c(1, 6, 9, 10, 16, 17, 18, 21, 26, 28:29, 34:35, 43, 45, 48:50, 54, 57, 60, 62:63, 66, 69:71, 73:74, 77:79, 
                     81, 83, 86, 92:93, 96, 98, 103, 106:108, 110, 112, 114:115, 120:122, 125, 127, 129, 131:134, 138, 142:145, 
                     148), ]
cShaOTU <- cShared[, 1]
cShaSeq <- cShared[, 2]
xTXs <- read.csv2("https://raw.githubusercontent.com/cmlglvz/datasets/master/Data/eAnalisis/xTXs.csv", 
                  sep = ";", dec = ".", header = TRUE, skip = 0)
rownames(xTXs) <- xTXs[, 2]
xTXs <- xTXs[, -1]
ShaTXs <- xTXs %>% filter(Seq %in% all_of(ShaSeq))
cShaTXs <- ShaTXs[-c(1, 6, 9, 10, 16, 17, 18, 21, 26, 28:29, 34:35, 43, 45, 48:50, 54, 57, 60, 62:63, 66, 69:71, 73:74, 77:79, 
                     81, 83, 86, 92:93, 96, 98, 103, 106:108, 110, 112, 114:115, 120:122, 125, 127, 129, 131:134, 138, 142:145, 
                     148), ]
cnombres <- data.frame(cShaOTU, mcOTU, dcOTU)
colnames(cShaASVs) <- dcOTU
cShaASVs <- cShaASVs[, mcOTU]

cHM <- heatmaply(normalize(cShaASVs), 
                 Colv = c.dend, 
                 Rowv = NA, 
                 margins = c(50, 50, 70, 0), 
                 grid_gap = 1,
                 width = 1920, 
                 height = 1080, 
                 subplot_heights = c(0.35, 0.65)
)
saveWidget(cHM, file = "ml_chmsasvs_reps.html")

cgShaASVs <- as.data.frame(t(cShaASVs))
cgShaASVs <- cgShaASVs %>% 
  mutate(Cha = rowSums(cgShaASVs[1:12]), 
         Fla = rowSums(cgShaASVs[13:24]), 
         Hu = rowSums(cgShaASVs[25:36]), 
         Pc = rowSums(cgShaASVs[37:41])
         )
cgShaASVs <- cgShaASVs[, -c(1:41)]
cgHM <- heatmaply(normalize(t(cgShaASVs)), 
                  Colv = c.dend, 
                  Rowv = NA, 
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
saveWidget(cgHM, file = "turbo_cghmgsha_reps.html")

cgHMt <- heatmaply(normalize(t(cgShaASVs)), 
                   Colv = c.dend, 
                   Rowv = NA,
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
saveWidget(cgHMt, file = "heat_cghmgsha_reps.html")

#Alternativa con generos pero sin dendrograma
Tax.sum <- function(OTU.Table, Tax.Table, Tax.lvl ){
  z <- NULL
  y <- NULL
  for (i in 1:length(unique(Tax.Table[colnames(OTU.Table),Tax.lvl]))) {
    if (length(OTU.Table[,which(Tax.Table[colnames(OTU.Table),Tax.lvl]==unique(Tax.Table[colnames(OTU.Table),Tax.lvl])[i])])!=length(rownames(OTU.Table))) {
      z <- which(Tax.Table[colnames(OTU.Table),Tax.lvl]==unique(Tax.Table[colnames(OTU.Table),Tax.lvl])[i])
      y <- cbind(y, apply(OTU.Table[,which(Tax.Table[colnames(OTU.Table),Tax.lvl]==unique(Tax.Table[colnames(OTU.Table),Tax.lvl])[i])], 1, function(x) sum(x)))
    } else { 
      y <- cbind(y, OTU.Table[,which(Tax.Table[colnames(OTU.Table),Tax.lvl]==unique(Tax.Table[colnames(OTU.Table),Tax.lvl])[i])])
    }
  }
  colnames(y) <- unique(Tax.Table[colnames(OTU.Table),Tax.lvl])
  invisible((y))
}

calor <- Tax.sum(ShaASVs, xTXs, 9)
cheat <- as.data.frame(t(calor))
cheat <- cheat %>% 
  mutate(Cha = rowSums(cheat[1:12]), 
         Fla = rowSums(cheat[13:24]), 
         Hu = rowSums(cheat[25:36]), 
         Pc = rowSums(cheat[37:41])
         )
cheat <- cheat[-15, -c(1:41)]

cheatmap <- heatmaply(normalize(t(cheat)), 
                      Colv = NA, 
                      Rowv = NA, 
                      margins = c(50, 50, 70, 0), 
                      grid_gap = 1, 
                      grid_color = "#DEE2E6", 
                      width = 1920, 
                      height = 1080, 
                      colors = heat.colors(n = 256, 
                                           alpha = 1, 
                                           rev = TRUE)
                      )
saveWidget(cheatmap, file = "cheatmap_reps.html")

vtree <- phytools::read.newick("pmskrd_ogtree.nwk")
tree <- force.ultrametric(vtree, method = "nnls")
plot(vtree, no.margin = TRUE, main = "Original Tree by MEGA7")
vtoh <- phylogram::as.dendrogram.phylo(tree)
hst <- as.hclust(vtoh)
plot(hclust, cex = 0.7)
devolve <- as.dendrogram(hclust)
plot(devolve)
c.dend <- vtoh %>% set("branches_lwd", 0.3) %>% ladderize()

joya <- readDNAMultipleAlignment("pmskrd.fasta", format = "fasta")
pMSKRD <- as(joya, "DNAStringSet")

corrected <- readDNAMultipleAlignment("sasvs_reps.fasta", format = "fasta")
jirillos <- corrected@unmasked@ranges@NAMES
ordillos <- mixedsort(jirillos, decreasing = FALSE) #no es necesario, solo para tener en orden las ASV

nShaASVs <- read.csv2(file = "/Users/Artemis/Documents/GitHub/datasets/Data/eAnalisis/nShaASVs.csv", header = TRUE, sep = ";", dec = ".", skip = 0)
rownames(nShaASVs) <- nShaASVs[, 1]
nShaASVs <- nShaASVs[, -1]
gsites <- nShaASVs %>% 
  mutate(Cha = rowSums(nShaASVs[1:12]), 
         Fla = rowSums(nShaASVs[13:24]), 
         Hu = rowSums(nShaASVs[25:36]), 
         Pc = rowSums(nShaASVs[37:41])
         )
gsites <- gsites[, -c(1:41)]
gsites <- as.data.frame(t(gsites))
nShaASVs <- as.data.frame(t(nShaASVs))
Reps <- select(nShaASVs, jirillos)

rViridis <- heatmaply(normalize(Reps), 
                 Colv = c.dend, 
                 Rowv = NA, 
                 main = "Representative shared ASVs across all samples clustered with maximum likelihood",
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
saveWidget(rViridis, file = "viridis_mlreps.html")

rHeat <- heatmaply(normalize(Reps), 
                   Colv = c.dend, 
                   Rowv = NA, 
                   main = "Representative shared ASVs across all samples clustered with maximum likelihood", 
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
saveWidget(rHeat, file = "heat_mlreps.html")

gReps <- select(gsites, jirillos)

gViridis <- heatmaply(normalize(gReps), 
                      Colv = c.dend, 
                      Rowv = NA, 
                      main = "Representative shared ASVs across all samples clustered with maximum likelihood",
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
saveWidget(gViridis, file = "viridis__grouped_mlreps.html")

gHeat <- heatmaply(normalize(gReps), 
                   Colv = c.dend, 
                   Rowv = NA, 
                   main = "Representative shared ASVs across all samples clustered with maximum likelihood", 
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
saveWidget(gHeat, file = "heat_grouped_mlreps.html")

cvirmap <- heatmaply(normalize(t(cheat)), 
                     Colv = NA, 
                     Rowv = NA, 
                     main = "PPE selected ASV shown as genus at the four sites", 
                     margins = c(50, 50, 70, 0), 
                     grid_gap = 1, 
                     width = 1920, 
                     height = 1080
                     )
saveWidget(cvirmap, file = "cvirmap_reps.html")

cheatmap <- heatmaply(normalize(t(cheat)), 
                      Colv = NA, 
                      Rowv = NA, 
                      main = "PPE selected ASV shown as genus at the four sites",
                      margins = c(50, 50, 70, 0), 
                      grid_gap = 1, 
                      grid_color = "#DEE2E6", 
                      width = 1920, 
                      height = 1080, 
                      colors = heat.colors(n = 256, 
                                           alpha = 1, 
                                           rev = TRUE)
)
saveWidget(cheatmap, file = "cheatmap_reps.html")
browseURL("cheatmap_reps.png")

gcheatmap <- ggheatmap(normalize(t(cheat)), 
                       Colv = NA, 
                       Rowv = NA, 
                       main = "PPE selected ASV group as genus for Cha, Fla, Hu and Pc", 
                       margins = c(50, 50, 70, 10), 
                       grid_gap = 1, 
                       grid_color = "#DEE2E6", 
                       widths = 1920, 
                       heights = 1080, 
                       colors = heat.colors(n =256, 
                                            alpha = 1, 
                                            rev = TRUE), 
                       hide_colorbar = FALSE)

pcheatmap <- heatmaply(normalize(t(cheat)), 
                      Colv = NA, 
                      Rowv = NA, 
                      main = "PPE selected ASV shown as genus at the four sites",
                      margins = c(50, 50, 70, 0), 
                      grid_gap = 1, 
                      grid_color = "#DEE2E6", 
                      width = 1920, 
                      height = 1080, 
                      colors = heat.colors(n = 256, 
                                           alpha = 1, 
                                           rev = TRUE), 
                      file = "cheatmap_reps.png")

pngHeat <- heatmaply(normalize(gReps), 
                   Colv = c.dend, 
                   Rowv = NA, 
                   main = "Representative shared ASVs across all samples clustered with maximum likelihood", 
                   margins = c(50, 50, 70, 0), 
                   grid_gap = 1, 
                   grid_color = "#DEE2E6", 
                   widths = 3840, 
                   heights = 2160, 
                   subplot_heights = c(0.35, 0.65), 
                   color = heat.colors(n = 256, 
                                       alpha = 1, 
                                       rev = TRUE), 
                   file = "pngheat.png"
)

ggHeat <- ggheatmap()







