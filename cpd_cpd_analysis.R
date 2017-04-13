# informer set
# drug to drug distance analysis using stitch data

options(stringsAsFactors = FALSE)
setwd("/fh/fast/kemp_c/olga_files/projects/informer_set/")
require(synapseClient)
source("code/awesome_visualization.R")
source("code/get_numeric_pubchem_id.R")
require(pvclust) 
require(ComplexHeatmap)
require(circlize)

INF <- read.delim(synGet("syn8049193")@filePath, sep="\t") # informer set
INF$PUBCHEM_CID <- get_numeric_pubchem_id(INF$PUBCHEM_CID)
  
# TODO: the uncompressed files in synapse should 
# be reloaded compressed::takes too long to pull

# chem - chem 
TMP <- read.delim(synGet("syn8361056")@filePath, sep="\t", header=T)
D <- TMP
D$chemical1 <- get_numeric_pubchem_id(TMP$chemical1)
D$chemical2 <- get_numeric_pubchem_id(TMP$chemical2)

D2 <- subset(D, chemical1 %in% unique(INF$PUBCHEM_CID) & chemical2 %in% unique(INF$PUBCHEM_CID))
# final stich set dim = 229 x229
# score matrix
SM <- tapply(D2$textmining, list(r=D2$chemical1, c=D2$chemical2), FUN = mean)

# plot missing values
pdf("is_chem_chem_score.pdf")
colstr(SM)
dev.off()

# pre-sort by missing values and replot
sm_col <- colSums(is.na(SM))
sm_row <- rowSums(is.na(SM))

pdf("is_chem_chem_score_ordered.pdf")
colstr(SM[names(sort(sm_row)), names(sort(sm_col))])
dev.off()

png(file      = "is_chem_chem_score_ordered.png",
    width     = 900,
    height    = 600,
    units     = "px",
    res       = NA,
    pointsize = 14)
colstr(SM[names(sort(sm_row)), names(sort(sm_col))])
dev.off()

# set NAs to 0
SM[is.na(SM)] <- 0

bar <- unlist(lapply(1:ncol(SM), function(COL){foo <- SM[,COL]; return(length(foo[foo!=0]))}))
names(bar) <- colnames(SM)

INCLUDE <- names(bar[bar != 0])

DF <- as.data.frame(SM)
DF1 <- subset(DF[, INCLUDE], row.names(DF) %in% INCLUDE) # 226 x 226

SM_CLUSTERS <- pvclust(DF1, nboot=1000)
SM_CLUSTERS_T <- pvclust(t(DF1), nboot=1000)

plot(SM_CLUSTERS)
pvrect(SM_CLUSTERS, alpha=0.90)

ht <- Heatmap(DF1,
              #col = colorRamp2(c(-2, 0, 4), c("green", "white", "red")), #colorRampPalette(rev(brewer.pal(11, 'RdBu')))(100)
              cluster_columns = SM_CLUSTERS$hclust,
              #cluster_rows = SM_CLUSTERS$hclust,
              show_row_names = FALSE,
              show_column_names = FALSE,
              #row_dend_side = "right",
              #row_dend_width = unit(2.5, "cm"),
              column_dend_side = "top",
              column_dend_height = unit(2.5, "cm"),
              name="Stitch")

draw(ht)




CLUST_SM <- pvpick(as.data.frame(DF1), alpha = 0.5)$clusters

############

DF2 <- as.matrix(scale(DF1))

heatmap(as.matrix(DF1), scale='none')

pdf("is_chem_chem_score_clust.pdf")
res <- pheatmap(DF1, scale = 'none',
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         #kmeans_k = 50,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()

png(file      = "is_chem_chem_score_clust.png",
    width     = 900,
    height    = 600,
    units     = "px",
    res       = NA,
    pointsize = 14)

pheatmap(DF1, scale = 'none',
                color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                #kmeans_k = 50,
                show_rownames = FALSE,
                show_colnames = FALSE)
dev.off()

df1.clust <- cbind(DF1, cluster = cutree(res$tree_row, k = 3))
head(df1.clust)

clust1 <- subset(df1.clust, cluster == 1)
rownames(clust1)
clust2 <- subset(df1.clust, cluster == 2)

clust3 <- subset(df1.clust, cluster == 3)

dev.off()

