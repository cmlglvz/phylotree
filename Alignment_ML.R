library(tidyverse)
library(Biostrings)
library(DECIPHER)
library(treeio)
library(ggtree)



fasta <- "D:/Documents/GitHub/phylotree/seqs_rcc.fasta"
RCC <- readDNAStringSet(fasta)
RCC <- OrientNucleotides(RCC)
aligned <- AlignSeqs(RCC, guideTree = NULL, iterations = 8, refinements = 1, gapOpening = -400, verbose = TRUE)
BrowseSeqs(aligned, htmlFile = "RCC_Aligned.html", colorPatterns = TRUE, highlight = NA, colWidth = Inf)
writeXStringSet(aligned, filepath = "D:/Documents/GitHub/phylotree/RCC_Decipher.fasta", format = "fasta")

Ali <- "D:/Documents/GitHub/phylotree/RCC_Decipher.fasta"
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(Ali, type = "FASTA", dbFile = dbConn, "")

# identify the sequences by their species
x <- dbGetQuery(dbConn, "select description from Seqs")$description
Add2DB(myData = data.frame(identifier = x, stringsAsFactors = FALSE), dbConn)

# form a consensus for each species
cons <- IdConsensus(dbConn, threshold = 0.3, minInformation = 0.1)

# calculate a maximum likelihood tree
d <- DistanceMatrix(cons, correction = "Jukes-Cantor", processors = NULL, verbose = TRUE)
dend <- IdClusters(d, method = "ML", showPlot = TRUE, type = "dendrogram", myXStringSet = cons, processors = 4, verbose = TRUE)
dend <- dendrapply(dend, 
                   FUN = function(n) {
                     if(is.leaf(n)) 
                       attr(n, "label") <- as.expression(substitute(italic(leaf), 
                                                                    list(leaf = attr(n, "label"))))
                     n
                     }
                   )

v.dend <- IdClusters(d, method = "ML", showPlot = FALSE, type = "dendrogram", myXStringSet = cons, model = c("TN93", "TN93+G4"), processors = 5, verbose = TRUE)

# display the phylogenetic tree
p <- par(mar=c(5, 1, 4, 10), xpd = TRUE)
plot(dend, yaxt = "n", horiz = TRUE)
arrows(-0.1, 6, -0.2, 6, angle = 90, length = 0.05, code = 3)
text(-0.15, 6, "0.1", adj = c(0.5, -0.5))
par(p)

dbDisconnect(dbConn)

WriteDendrogram(dend, file = "D:/Documents/GitHub/phylotree/aligned_dendro_rcc.nwk", quoteLabels = TRUE, convertBlanks = FALSE, internalLabels = TRUE)
WriteDendrogram(v.dend, file = "D:/Documents/GitHub/phylotree/vdendro.nwk", quoteLabels = TRUE, convertBlanks = FALSE, internalLabels = TRUE)
write.dendrogram(dend, file = "D:/Documents/GitHub/phylotree/aligned_dendro_rcc_vd.nwk", edges = TRUE)


nwk <- read.newick(file = "D:/Documents/GitHub/phylotree/aligned_dendro_rcc.nwk", node.label = "label")
o.nwk <- read.tree(file = "D:/Documents/GitHub/phylotree/aligned_dendro_rcc.nwk")
vNWK <- read.newick(file = "D:/Documents/GitHub/phylotree/vdendro.nwk", node.label = "label")

ggplot(nwk, branch.length='none', aes(x, y)) + geom_tree() + layout_inward_circular(xlim = 200) + theme_tree2()
ggplot(nwk, branch.length='none') + geom_tree() + layout_inward_circular(xlim = 200) + theme_tree()
ggplot(nwk, layout = "fan", open.angle=120) + geom_tree() + theme_tree()
ggplot(nwk, ladderize=FALSE, aes(x, y)) + geom_tree() + theme_tree()
ggplot(o.nwk, aes(x, y)) + geom_tree() + theme_tree()

phyl <- as.phylo(dend)


ggplot(v.dend, ladderize = TRUE, branch.length = "none") + geom_tree() + theme_tree() + geom_tiplab(geom = "text", size = 4)
ggplot(v.dend, ladderize = FALSE, branch.length = "none") + geom_tree() + theme_tree() + geom_tiplab(geom = "text", size = 4)
ggplot(v.dend, 
       ladderize = FALSE, 
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
ggplot(v.dend, 
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
ggplot2::ggsave(plot = t, filename = "Tree.png", path = "D:/Documents/GitHub/phylotree/", dpi = 600)
