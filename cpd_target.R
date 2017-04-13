# informer set
# drug to drug distance analysis using stich data

options(stringsAsFactors = FALSE)
setwd("/fh/fast/kemp_c/olga_files/projects/informer_set/")
require(synapseClient)
source("code/get_numeric_pubchem_id.R")
require("gProfileR")
require('biomaRt')

#================================================
# (1) enrichment analysis for targets annotated 
#     via CTRP, DrugBank and Selleck

# read in integated annotation file
DAT <- read.delim(synGet("syn8647538")@filePath, sep="\t", header=T)

# aggregate targets so far
# include columns: 
# gene_symbol_of_protein_target (CTRP) = 131 compounds
# DB_Gene.Name (DrugBank) = 22 compounds
# selleck_Target (Selleck) = 22 compounds; targets are not gene names in most cases; excluded after all

nrow(subset(DAT, is.na(gene_symbol_of_protein_target) == T & is.na(DB_Gene.Name) == T)) # 166 remain without target information

TARG1 <- c(DAT$gene_symbol_of_protein_target,
           DAT$DB_Gene.Name)#,
           #DAT$selleck_Target)

TARG2 <- unique( unlist(lapply(TARG1, strsplit, ";")) )
TARG3 <- TARG2[!is.na(TARG2)] # total: 223

# pathway analysis with background library
# background: sigma+hutch druggable library
LIB <- read.delim(synGet("syn7104369")@filePath, sep="\t")
LIB_GENES <- unique(LIB$Updated_symbol) # N =  6,600

inboth <- intersect(TARG3, LIB_GENES) # 191

initTARG <- gprofiler(TARG3, 
                      custom_bg = LIB_GENES,
                      ordered_query=F,
                      correction_method="fdr",
                      max_p_value=0.1,
                      png_fn=NULL) # total: 

initPATH <- subset(initTARG, domain %in% c("keg", "rea")) # 349 
initGO <- subset(initTARG, domain %in% c("BP","CC","MF")) # 2158

# upload results to synapse
write.table(initPATH, file="results/gprofiler_init_targets_pathway_keg_rea.tsv", sep="\t", col.names=T, row.names=F, quote=F)
write.table(initGO, file="results/gprofiler_init_targets_go.tsv", sep="\t", col.names=T, row.names=F, quote=F)

file <- File("results/gprofiler_init_targets_pathway_keg_rea.tsv", parentId = "syn8367102")
file <- synStore(file)

file <- File("results/gprofiler_init_targets_go.tsv", parentId = "syn8367102")
file <- synStore(file)

#================================================
# (2) enrichment analysis for targets annotated 
#     by integrating with STITCH
INF <- read.delim(synGet("syn8049193")@filePath, sep="\t") # informer set
INF$pubchem_cid <- get_numeric_pubchem_id(INF$PUBCHEM_CID)

# detect duplicated pubch_ids and set these aside to have unique mapping;
# will add back at the end
TMP_REMOVE <- !duplicated(INF$pubchem_cid)
INF1 <- INF[TMP_REMOVE,]

# integrate stitch data
SCH <- read.delim(synGet("syn8360870")@filePath, sep="\t")
SCH$chemical <- get_numeric_pubchem_id(SCH$chemical)
SCH_SUB <- unique(subset(SCH, chemical %in% INF1$pubchem_cid )) # 54,927 pairs for 268 unique drugs
#SCH_SUB1 <- SCH_SUB[grep("9606.", SCH_SUB$protein),]
SCH_SUB$protein <- gsub("9606.","",SCH_SUB$protein)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
qq <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"), values=unique(SCH_SUB$protein), mart= mart)

SCH_SUB2 <- merge(SCH_SUB, qq, by.x = "protein", by.y = "ensembl_peptide_id", all.x = T) # 54,927

SCH_FIN <- unique(subset(SCH_SUB2[,c("chemical","protein","hgnc_symbol","combined_score")], is.na(hgnc_symbol) == F) )# 52,845; 12710 unique proteins

# gprofiler on all targets 
# inboth <- intersect(unique(SCH_FIN$hgnc_sym), LIB_GENES)

TARG <- gprofiler(unique(SCH_FIN$hgnc_symbol), 
                  #custom_bg=LIB_GENES, 
                  ordered_query=F,
                  correction_method="fdr",
                  max_p_value=0.1,
                  png_fn=NULL)

PATH <- subset(TARG, domain %in% c("keg", "rea")) # 615
GO <- subset(TARG, domain %in% c("BP","CC","MF")) # 4349

write.table(PATH, file="results/gprofiler_stitch_targets_pathway_keg_rea.tsv", sep="\t", col.names=T, row.names=F, quote=F)
write.table(GO, file="results/gprofiler_stitch_targets_go.tsv", sep="\t", col.names=T, row.names=F, quote=F)

file <- File("results/gprofiler_stitch_targets_pathway_keg_rea.tsv", parentId = "syn8367102")
file <- synStore(file)

file <- File("results/gprofiler_stitch_targets_go.tsv", parentId = "syn8367102")
file <- synStore(file)

# prep for cytoscape visualization

### (1) Drug and target, scored
write.table(SCH_FIN, file="results/strich_chem_protein_for_cytoscape.tsv", sep="\t", col.names=T, row.names=F, quote=F)
file <- File("results/strich_chem_protein_for_cytoscape.tsv", parentId = "syn8367102")
file <- synStore(file)

png("results/cpr_target_score_dn_stitch.png")
hist(SCH_FIN$combined_score,col="grey", 
     main="STITCH compound-target score distribution",
     xlab="score", ylab="count")
dev.off()
file <- File("results/cpr_target_score_dn_stitch.png", parentId = "syn8367102")
file <- synStore(file)

part <- subset(SCH_FIN, combined_score > 600)
write.table(part, file="results/strich_chem_protein_for_cytoscape_part.tsv", sep="\t", col.names=T, row.names=F, quote=F)
file <- File("results/strich_chem_protein_for_cytoscape_part.tsv", parentId = "syn8367102")
file <- synStore(file)





