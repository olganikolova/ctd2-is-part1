# informer set
# clinical status annotation

options(stringsAsFactors = FALSE)
setwd("Projects/CTDD/informer-set/")
require(webchem)
require(synapseClient)
synapseLogin()

INF <- read.delim(synGet("syn8049193")@filePath, sep="\t") # informer set
names(INF) <- tolower(names(INF))

# detect duplicated pubch_ids and set these aside to have unique mapping;
# will add back at the end
TMP_REMOVE <- !duplicated(INF$pubchem_cid)
INF1 <- INF[TMP_REMOVE,]

ANN <- read.delim(synGet("syn8049194")@filePath, sep = "\t") # clinical status annotation in CTRP

# Warning message:
#   In createS4ObjectFromList(elem, listElementType) :
#   Entry storageLocationId found in list but S3UploadDestination has no such slot.

# inf set annotated via CTRP cpd_status:
# 133 have clinical status
# 191 (188 unique) are missing
RES <- merge(INF1, ANN, by = "broad_cpd_id", all.x = TRUE)
RES$pubchem_cid <- gsub("CID:","", RES$pubchem_cid)

# retrieve available CAS identifiers from pub chem ids
cas <- cts_convert(query = RES$pubchem_cid, from = 'PubChem CID', to = 'CAS')

cas1 <- lapply(foo, function(STR){
  return(paste(STR, collapse=";"))
})

RES$webchem_CAS <- do.call(rbind, cas1)

# BRANCH 1 ===========================================================
# mapping between pubchem_ids and inchi aquired from pubchem id mapper
MAP <- read.delim(synGet("syn8049778")@filePath, sep="\t", header=F) 
names(MAP) <- c("inchi", "pub_chem_id")

# subset by the informer set
MAP1 <- unique(subset(MAP, pub_chem_id %in% RES$pubchem_cid)) # only 105 drugs could be mapped this way

# map pubchemid to inchi
RES1 <- merge(RES, MAP1, by.x = "pubchem_cid", by.y = "pub_chem_id", all.x = TRUE)

# mapping between drug bank ids and inChI keys from drug bank
# inchi keys not available for nutraceuticals!
VOC <- read.delim(synGet("syn8049174")@filePath, sep=",")
VOC1 <- unique(subset(VOC[,c("DrugBank.ID","CAS","Standard.InChI.Key")], Standard.InChI.Key != ""))
names(VOC1) <- gsub("CAS", "DB_CAS", names(VOC1))

RES2 <- merge(RES1, VOC1, by.x = "inchi", by.y = "Standard.InChI.Key", all.x = TRUE)
# 3081361 seems to have 2 different drug bank ids; leaving both in (total 319)

# synapse query to get different cpd_status groups from drug bank
q <- synQuery("select id, name from entity where entity.parentId=='syn8049148'")
qq <- subset(q, entity.name %in% grep("_all", q$entity.name, value=T))

# by cpd_status list of drug to targets
DAT <- lapply(qq$entity.id, function(SID){
  # read in data
  D <- read.delim(synGet(SID)@filePath, sep=",")
  # single drug ID per row
  DD <- lapply(1:nrow(D), function(ROW){
    X <- D[ROW,]
    STR <- unlist(strsplit(X$Drug.IDs,";"))
    Y <- as.matrix(X)
    Z <- as.data.frame(matrix(rep(Y, length(STR)), nrow=length(STR), byrow = T))
    names(Z) <- names(X)
    Z$Drug.IDs <- STR
    return(Z)
  })
  U <- as.data.frame(do.call(rbind, DD))
  # unique drug id per row (accumulate all drug targets in single row)
  L <- split(U, U$Drug.IDs)
  
  LL <- lapply(L, function(W){
    ID <- paste(W$ID, collapse=";")
    Name <- paste(W$Name, collapse=";") 
    Gene.Name <- paste(W$Gene.Name, collapse=";") 
    GeneCard.ID <- paste(W$GeneCard.ID, collapse=";") 
    GenAtlas.ID <- paste(W$GenAtlas.ID, collapse=";")    
    HGNC.ID <- paste(W$HGNC.ID, collapse=";")
    Species <- paste(W$Species, collapse=";")
    
    WW <- data.frame(ID = ID,
                     Name = Name,
                     Gene.Name = Gene.Name,
                     GeneCard.ID = GeneCard.ID,
                     GenAtlas.ID = GenAtlas.ID,
                     HGNC.ID = HGNC.ID,
                     Species = Species,
                     Drug.IDs = W$Drug.IDs[1])
    return(WW)
  })
  
  TMP <- as.data.frame(do.call(rbind, LL))
  
  return(TMP)
})
names(DAT) <- gsub("_all.csv","", qq$entity.name)

# summary of Drug Bank data per category
sapply(DAT, nrow)
# approved    experimental         illicit investigational       withdrawn 
# 789            4877             108            1007             133 

# aux drug ids list to assess overall overlap with ctd2-is
DAT2 <- lapply(DAT, function(D){
  TMP <- subset(D, Drug.IDs %in% RES2$DrugBank.ID)
  return(TMP$Drug.IDs)
})

do.call(sum, lapply(DAT2, length))
[1] 48

# overlapping categories
length(unique(unlist(DAT2)))
all <- unlist(DAT2)

# intersect(DAT2$investigational,DAT2$approved)
# [1] "DB00328" "DB00987" "DB01143"
# > intersect(DAT2$investigational,DAT2$withdrawn)
# [1] "DB00533" "DB01041"

# are any of DB-annotated drugs help resolve
# drugs currently missing CTRP annotation?
DAT3 <- lapply(1:length(DAT), function(IDX){
  D <- DAT[[IDX]]
  TMP <- subset(D, Drug.IDs %in% RES2$DrugBank.ID)
  names(TMP) <- paste("DB", names(TMP), sep="_")
  x <- merge(RES2, TMP, by.x = "DrugBank.ID", by.y = "DB_Drug.IDs")
  x$DB_cpd_status <- rep(names(DAT)[IDX], nrow(x))
  return(x)
})
names(DAT3) <- names(DAT)

RES3 <- as.data.frame(do.call(rbind, DAT3))
write.table(RES3, file="ctd2is_additional_cps_status_annot_inchi.txt", sep="\t",row.names=F)

nrow(subset(RES3, is.na(cpd_status) == T)) # total of 25 previously un-annotated compounds

L3 <- split(RES3, RES3$pubchem_cid)

# # BRANCH 2 ===========================================================
# # using CAS identifiers
# 
# INF_UN_CAS <- unique(unlist(cas))
# VOC2 <- unique(subset(VOC, CAS %in% INF_UN_CAS)) # only 107
# 
# L4 <- lapply(RES$webchem_CAS, function(STR){
#   L <- unlist(strsplit(STR,";"))
#   return(paste(subset(VOC2, CAS %in% L)$DrugBank.ID, collapse=";"))
# })
# names(L4) <- RES$webchem_CAS
# 
# RES4 <- as.vector(do.call(cbind, L4))
# names(RES4) <- names(L4)
# RES5 <- RES4[RES4!=""]
# RES6 <- data.frame(CAS = names(RES5),
#                    DB_ID =RES5)
# row.names(RES6) <- NULL
# 
# RES7 <- merge(subset(RES, webchem_CAS != ""), RES6, by.x = "webchem_CAS", by.y = "CAS")
# 
# DAT4 <- lapply(1:length(DAT), function(IDX){
#   D <- DAT[[IDX]]
#   IDS <- unlist(strsplit(RES6$DB_ID,";"))
#   TMP <- subset(D, Drug.IDs %in% IDS)
#   names(TMP) <- paste("DB", names(TMP), sep="_")
#   
#   
#   
#   x <- merge(RES7, TMP, by.x = "webchem_CAS", by.y = "DB_Drug.IDs")
#   x$DB_cpd_status <- rep(names(DAT)[IDX], nrow(x))
#   return(x)
# })
# names() <- names(DAT)
