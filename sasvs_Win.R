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
library(vegan)
library(BiodiversityR)
library(ggrepel)

sASVs <- readDNAMultipleAlignment("sASVs.afa", format = "fasta")
sDNAStr <- as(sASVs, "DNAStringSet")
BrowseSeqs(sDNAStr, htmlFile = "sasvs.html")
writeXStringSet(sDNAStr, filepath = "sasvs.fas", format = "fasta")
autoMasked <- maskGaps(sASVs, min.fraction = 0.3, min.block.width = 4)
ATMskd <- as(autoMasked, "DNAStringSet")
writeXStringSet(ATMskd, filepath = "sasvsmmskd.fasta", format = "fasta")
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

eShaTXs <- ShaTXs[, -c(1:3)]
rownames(eShaTXs) <- dOTU
eShaTXs <- eShaTXs[mOTU, ]
asvtable <- as.matrix(t(eShaASVs))
taxmat <- as.matrix(eShaTXs)
EnvMe <- read.csv2(file = "/Users/Artemis/Documents/GitHub/datasets/Data/eAnalisis/EnvMe.csv", header = TRUE, sep = ";", dec = ".", row.names = 1, skip = 0)
OTU <- otu_table(object = asvtable, taxa_are_rows = TRUE)
TAX <- tax_table(object = taxmat)
SAMPLE <- sample_data(object = EnvMe)
ASV <- read.FASTA(file = "C:/Users/Camilo/Dropbox/sasvsmmskd.fasta", type = "DNA")
pASV <- phyDat(ASV, type = "DNA", levels = NULL, return.index = TRUE)
dm <- dist.ml(pASV)
tNJ <- NJ(dm)
mt <- modelTest(pASV, tree = tNJ)
bfit <- as.data.frame(mt)
write.csv2(bfit, "C:/Users/Camilo/Desktop/bfit.csv")

fit <- pml(tNJ, data = pASV)
fitG <- optim.pml(fit, TRUE)
fitT <- update(fit, k = 4, inv = 0.26)
fitGTR <- optim.pml(fitT, 
                    model = "GTR", 
                    optNni = TRUE, 
                    optBf = TRUE, 
                    optGamma = TRUE, 
                    rearrangement = "stochastic", 
                    control = pml.control(trace = 0)
)
bs <- bootstrap.pml(fitGTR, bs = 1000, optNni = TRUE, control = pml.control(trace = 0))
plotBS(midpoint(fitGTR$tree), bs, p = 50, type = "p")

physq <- phyloseq(OTU, TAX, SAMPLE, phytools::midpoint.root(fitGTR$tree))

plot_tree(physq, color = "Site", label.tips = "taxa_names", ladderize = "left", plot.margin = 0.3)
plot_heatmap(physq)
plot_tree(physq, size = "abundance", color = "Site", label.tips = "taxa_names")

trank <- rankabundance(eShaASVs)
write.csv2(trank, "total_rankabundance.csv")
rankabunplot(xr = trank, scale = "proportion", addit = FALSE, scaledx = FALSE, specnames = c(1:10))

charank <- rankabundance(eShaASVs[c(1:12), ])
rankabunplot(charank, scale = "proportion", addit = FALSE, scaledx = FALSE, specnames = c(1:10))

flarank <- rankabundance(eShaASVs[c(13:24), ])
rankabunplot(flarank, scale = "proportion", addit = FALSE, scaledx = FALSE, specnames = c(1:10))

huarank <- rankabundance(eShaASVs[c(24:36)], )
rankabunplot(huarank, scale = "proportion", addit = FALSE, scaledx = FALSE, specnames = c(1:10))

pcrank <- rankabundance(eShaASVs[c(37:41), ])
rankabunplot(pcrank, scale = "proportion", addit = FALSE, scaledx = FALSE, specnames = c(1:10))

rcomp <- rankabuncomp(eShaASVs, EnvMe, factor = "Site", scale = "proportion", return.data = TRUE, legend = TRUE, specnames = c(1:10))
write.csv2(rcomp, "composite_rankabundance.csv")
compRA <- read.csv2("eCompRAbun.csv", header = TRUE, sep = ";", dec = ",", skip = 0, row.names = 1, fill = TRUE)

BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12, family="Arial"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

prcomp <- ggplot(data = compRA, aes(x = rank, y = proportion)) + 
  scale_x_continuous(expand = c(0, 1), sec.axis = dup_axis(labels = NULL, name = NULL)) +
  scale_y_continuous(expand = c(0, 2), sec.axis = dup_axis(labels = NULL, name = NULL)) +
  geom_line(aes(colour = Grouping), size = 1) +
  geom_point(aes(colour = Grouping), size = 3, alpha = 0.7) +
  geom_text_repel(data = subset(compRA, labelit == TRUE), 
                  aes(label = species), 
                  angle = 5, nudge_x = 2, nudge_y = 1, show.legend = FALSE, max.overlaps = Inf) +
  BioR.theme +
  scale_color_manual(values = c("#440154", "#1C3B74", "#20A387", "#FDE725")) +
  facet_wrap(~ Grouping) +
  labs(x = "Rank", y = "Abundance Proportion", colour = "Site")
prcomp

ggsave(prcomp, file = "composite_rankabundance.png", path = "/Users/Artemis/Documents/GitHub/phylotree/", width = 13, height = 13, dpi = 600)

ChlIASVs <- select(eShaASVs, all_of(lChl_I))
CURAbn <- rankabuncomp(ChlIASVs, EnvMe, factor = "Site", scale = "proportion", return.data = TRUE, legend = TRUE, specnames = c(1:10))
write.csv2(CURAbn, "Chl_I_RAbun.csv")
#edited externally
ChlI <- read.csv2("eChl_I_RAbun.csv", header = TRUE, sep = ";", dec = ",", row.names = 1, skip = 0, fill = TRUE)
ChlUno <- ggplot(data = ChlI, aes(x = rank, y = proportion)) + 
  scale_x_continuous(expand = c(0, 1), sec.axis = dup_axis(labels = NULL, name = NULL)) + 
  scale_y_continuous(expand = c(0, 2), sec.axis = dup_axis(labels = NULL, name = NULL)) + 
  geom_line(aes(colour = Grouping), size = 1) + 
  geom_point(aes(colour = Grouping), size = 3, alpha = 0.7) + 
  geom_text_repel(data = subset(ChlI, labelit == TRUE), 
                  aes(label = species), 
                  angle = 5, nudge_x = 2, nudge_y = 1, show.legend = FALSE, max.overlaps = Inf) + 
  BioR.theme + 
  scale_color_manual(values = c("#440154", "#1C3B74", "#20A387", "#FDE725")) + 
  facet_wrap(~ Grouping) + 
  labs(x = "Rank", y = "Proportional Abundance", colour = "Site")
ggsave(ChlUno, file = "Chl_I_CRAbun.png", path = "/Users/Artemis/Documents/GitHub/phylotree/", width = 13, height = 13, dpi = 600)

ChlIIASVs <- select(eShaASVs, all_of(lChl_II))
CDRAbn <- rankabuncomp(ChlIIASVs, EnvMe, factor = "Site", scale = "proportion", return.data = TRUE, legend = TRUE, specnames = c(1:10))
write.csv2(CDRAbn, "Chl_II_RAbun.csv")
#edited externally
ChlII <- read.csv2("eChl_II_RAbun.csv", header = TRUE, sep = ";", dec = ",", row.names = 1, skip = 0, fill = TRUE)
ChlDos <- ggplot(data = ChlII, aes(x = rank, y = proportion)) + 
  scale_x_continuous(expand = c(0, 1), sec.axis = dup_axis(labels = NULL, name = NULL)) + 
  scale_y_continuous(expand = c(0, 2), sec.axis = dup_axis(labels = NULL, name = NULL)) + 
  geom_line(aes(colour = Grouping), size = 1) + 
  geom_point(aes(colour = Grouping), size = 3, alpha = 0.7) + 
  geom_text_repel(data = subset(ChlII, labelit == TRUE), 
                  aes(label = species), 
                  angle = 5, nudge_x = 2, nudge_y = 1, show.legend = FALSE, max.overlaps = Inf) + 
  BioR.theme + 
  scale_color_manual(values = c("#440154", "#1C3B74", "#20A387", "#FDE725")) + 
  facet_wrap(~ Grouping) + 
  labs(x = "Rank", y = "Proportional Abundance", colour = "Site")
ggsave(ChlDos, file = "Chl_II_CRAbun.png", path = "/Users/Artemis/Documents/GitHub/phylotree/", width = 13, height = 13, dpi = 600)

HapASVs <- select(eShaASVs, all_of(lHap))
HRAbn <- rankabuncomp(HapASVs, EnvMe, factor = "Site", scale = "proportion", return.data = TRUE, legend = TRUE, specnames = c(1:10))
write.csv2(HRAbn, "Hap_RAbun.csv")
#externally edited
Hap <- read.csv2("eHap_RAbun.csv", header = TRUE, sep = ";", dec = ",", row.names = 1, skip = 0, fill = TRUE)
pHap <- ggplot(data = Hap, aes(x = rank, y = proportion)) + 
  scale_x_continuous(expand = c(0, 1), sec.axis = dup_axis(labels = NULL, name = NULL)) + 
  scale_y_continuous(expand = c(0, 2), sec.axis = dup_axis(labels = NULL, name = NULL)) + 
  geom_line(aes(colour = Grouping), size = 1) + 
  geom_point(aes(colour = Grouping), size = 3, alpha = 0.7) + 
  geom_text_repel(data = subset(Hap, labelit == TRUE), 
                  aes(label = species), 
                  angle = 5, nudge_x = 2, nudge_y = 1, show.legend = FALSE, max.overlaps = Inf) + 
  BioR.theme + 
  scale_color_manual(values = c("#440154", "#1C3B74", "#20A387", "#FDE725")) + 
  facet_wrap(~ Grouping) + 
  labs(x = "Rank", y = "Proportional Abundance", colour = "Site")
ggsave(pHap, file = "Hap_CRAbun.png", path = "/Users/Artemis/Documents/GitHub/phylotree/", width = 13, height = 13, dpi = 600)

CKASVs <- select(eShaASVs, all_of(lCK))
CKRAbn <- rankabuncomp(CKASVs, EnvMe, factor = "Site", scale = "proportion", return.data = TRUE, legend = TRUE, specnames = c(1:10))
write.csv2(CKRAbn, "CK_RAbun.csv")
#Ex edt
CK <- read.csv2("eCK_RAbun.csv", header = TRUE, sep = ";", dec = ",", row.names = 1, skip = 0, fill = TRUE)
pCK <- ggplot(data = CK, aes(x = rank, y = proportion)) + 
  scale_x_continuous(expand = c(0, 1), sec.axis = dup_axis(labels = NULL, name = NULL)) + 
  scale_y_continuous(expand = c(0, 2), sec.axis = dup_axis(labels = NULL, name = NULL)) + 
  geom_line(aes(colour = Grouping), size = 1) + 
  geom_point(aes(colour = Grouping), size = 3, alpha = 0.7) + 
  geom_text_repel(data = subset(CK, labelit == TRUE), 
                  aes(label = species), 
                  angle = 5, nudge_x = 2, nudge_y = 1, show.legend = FALSE, max.overlaps = Inf) + 
  BioR.theme + 
  scale_color_manual(values = c("#440154", "#1C3B74", "#20A387", "#FDE725")) + 
  facet_wrap(~ Grouping) + 
  labs(x = "Rank", y = "Proportional Abundance", colour = "Site")
ggsave(pCK, file = "CK_CRAbun.png", path = "/Users/Artemis/Documents/GitHub/phylotree/", width = 15, height = 15, dpi = 600)

OchIASVs <- select(eShaASVs, all_of(lOch_I))
OURAbn <- rankabuncomp(OchIASVs, EnvMe, factor = "Site", scale = "proportion", return.data = TRUE, legend = TRUE, specnames = c(1:10))
write.csv2(OURAbn, "Och_I_RAbun.csv")
#ext edt
OchI <- read.csv2("eOch_I_RAbun.csv", header = TRUE, sep = ";", dec = ",", row.names = 1, skip = 0, fill = TRUE)
OchUno <- ggplot(data = OchI, aes(x = rank, y = proportion)) + 
  scale_x_continuous(expand = c(0, 1), sec.axis = dup_axis(labels = NULL, name = NULL)) + 
  scale_y_continuous(expand = c(0, 2), sec.axis = dup_axis(labels = NULL, name = NULL)) + 
  geom_line(aes(colour = Grouping), size = 1) + 
  geom_point(aes(colour = Grouping), size = 3, alpha = 0.7) + 
  geom_text_repel(data = subset(OchI, labelit == TRUE), 
                  aes(label = species), 
                  angle = 5, nudge_x = 2, nudge_y = 1, show.legend = FALSE, max.overlaps = Inf) + 
  BioR.theme + 
  scale_color_manual(values = c("#440154", "#1C3B74", "#20A387", "#FDE725")) + 
  facet_wrap(~ Grouping) + 
  labs(x = "Rank", y = "Proportional Abundance", colour = "Site")
ggsave(OchUno, file = "Och_I_CRAbun.png", path = "/Users/Artemis/Documents/GitHub/phylotree/", width = 15, height = 15, dpi = 600)

OchIIASVs <- select(eShaASVs, all_of(lOch_II))
ODRAbn <- rankabuncomp(OchIIASVs, EnvMe, factor = "Site", scale = "proportion", return.data = TRUE, legend = TRUE, specnames = c(1:10))
write.csv2(ODRAbn, "Och_II_RAbun.csv")
#ext edt
OchII <- read.csv2("eOch_II_RAbun.csv", header = TRUE, sep = ";", dec = ",", row.names = 1, skip = 0, fill = TRUE)
OchDos <- ggplot(data = OchII, aes(x = rank, y = proportion)) + 
  scale_x_continuous(expand = c(0, 1), sec.axis = dup_axis(labels = NULL, name = NULL)) + 
  scale_y_continuous(expand = c(0, 2), sec.axis = dup_axis(labels = NULL, name = NULL)) + 
  geom_line(aes(colour = Grouping), size = 1) + 
  geom_point(aes(colour = Grouping), size = 3, alpha = 0.7) + 
  geom_text_repel(data = subset(OchII, labelit == TRUE), 
                  aes(label = species), 
                  angle = 5, nudge_x = 2, nudge_y = 1, show.legend = FALSE, max.overlaps = Inf) + 
  BioR.theme + 
  scale_color_manual(values = c("#440154", "#1C3B74", "#20A387", "#FDE725")) + 
  facet_wrap(~ Grouping) + 
  labs(x = "Rank", y = "Proportional Abundance", colour = "Site")
ggsave(OchDos, file = "Och_II_CRAbun.png", path = "/Users/Artemis/Documents/GitHub/phylotree/", width = 13, height = 13, dpi = 600)



















