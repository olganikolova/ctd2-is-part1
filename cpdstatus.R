# informer set
# clinical status annotation

options(stringsAsFactors = FALSE)
setwd("/fh/fast/kemp_c/olga_files/projects/informer_set/")
require(webchem)
require(synapseClient)
require(VennDiagram)

INF <- read.delim(synGet("syn8049193")@filePath, sep="\t") # informer set
names(INF) <- tolower(names(INF))

# detect duplicated pubch_ids and set these aside to have unique mapping;
# will add back at the end
TMP_REMOVE <- !duplicated(INF$pubchem_cid)
INF1 <- INF[TMP_REMOVE,]

ANN <- read.delim(synGet("syn8049194")@filePath, sep = "\t") # clinical status annotation in CTRP
# STR <- paste("CTRP", names(ANN), sep="_")
# names(ANN) <- STR

# Warning message:
#   In createS4ObjectFromList(elem, listElementType) :
#   Entry storageLocationId found in list but S3UploadDestination has no such slot.

# inf set annotated via CTRP cpd_status:
# 133 have clinical status
# 191 (188 unique) are missing
RES <- merge(INF1, ANN, by = "broad_cpd_id", all.x = TRUE)
RES$pubchem_cid <- gsub("CID:","", RES$pubchem_cid)

# BRANCH 1 ===========================================================
# mapping between pubchem_ids and inchi aquired from pubchem id mapper
# this is a 2-column map only
MAP <- read.delim(synGet("syn8049778")@filePath, sep="\t", header=F) 
names(MAP) <- c("inchi", "pub_chem_id")
MAP$pub_chem_id <- as.character(MAP$pub_chem_id)

# subset by the informer set
MAP1 <- unique(subset(MAP, pub_chem_id %in% RES$pubchem_cid)) # only 105 drugs could be mapped this way

# map pubchemid to inchi
RES1 <- merge(RES, MAP1, by.x = "pubchem_cid", by.y = "pub_chem_id", all.x = TRUE)

# mapping between drug bank ids and inChI keys from drug bank
# inchi keys not available for nutraceuticals!
VOC <- read.delim(synGet("syn8049174")@filePath, sep=",")
VOC1 <- unique(subset(VOC[,c("DrugBank.ID","Standard.InChI.Key")], Standard.InChI.Key != ""))

RES2 <- merge(RES1, VOC1, by.x = "inchi", by.y = "Standard.InChI.Key", all.x = TRUE)
# 3081361 seems to have 2 different drug bank ids; leaving both in (total 319)
# at the end, neither of those DrugBank ids ends up being status-annotated   <<======================= TODO
# so could remove ither one to keep the pubchem ids unique

# synapse query to get different cpd_status groups from drug bank
q <- synQuery("select id, name from entity where entity.parentId=='syn8049148'")
qq <- subset(q, entity.name %in% grep("_all", q$entity.name, value=T))

# by cpd_status list of drug to targets from DrugBank
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
#[1] 48 # total drugs from DrugBank also included in ctd2-is

# overlapping categories
length(unique(unlist(DAT2)))
all <- unlist(DAT2)

# intersect(DAT2$investigational,DAT2$approved)
# [1] "DB00328" "DB00987" "DB01143"
# > intersect(DAT2$investigational,DAT2$withdrawn)
# [1] "DB00533" "DB01041"

# do any of DB-annotated drugs help resolve
# drugs currently missing CTRP annotation?
DAT3 <- lapply(1:length(DAT), function(IDX){
  D <- DAT[[IDX]] # DrugBank compound category
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

# are there drugs with multiple clinical annotation classes?
#L3 <- split(RES3, RES3$pubchem_cid)
# 5 drugs appear to be
# we will aggregate using the closest to approval status

# record all ctd2-is with additional DB annotation
DAT4 <- lapply(1:length(DAT), function(IDX){
  D <- DAT[[IDX]] # DrugBank compound category
  TMP <- subset(D, Drug.IDs %in% RES2$DrugBank.ID)
  names(TMP) <- paste("DB", names(TMP), sep="_")
  TMP$DB_cpd_status <- rep(names(DAT)[IDX], nrow(TMP))
  return(TMP)
})
names(DAT4) <- names(DAT)

RES4 <- as.data.frame(do.call(rbind, DAT4))

REMOVE <- which(duplicated(RES4$DB_Drug.IDs)) # removes 5 compounds
RES5 <- RES4[-REMOVE,]

RES6 <- merge(RES2, RES5, by.x = "DrugBank.ID", by.y = "DB_Drug.IDs", all.x = T)

write.table(RES6, file="ctd2is_all_cps_status_annot_drugbank-inchi.txt", sep="\t",row.names=F)

##### parsing selleck to see if we can adda annotation
# read in fda selleck compounds
# SELLECK <- read.delim(synGet("syn8646135")@filePath, sep="\t", header=T) # syn8646135 # FDA: syn8638067
# SELLECK_MAP <- read.delim(synGet("syn8646137")@filePath, sep="\t", header=F) # syn8646137 # FDA: syn8642971

SELLECK <- read.delim(synGet("syn8638067")@filePath, sep="\t", header=T) # syn8646135 # FDA: syn8638067
SELLECK_MAP <- read.delim(synGet("syn8642971")@filePath, sep="\t", header=F) # syn8646137 # FDA: syn8642971

names(SELLECK_MAP) <- c("product_id","pubchem_id")
SELLECK_MAP$product_id <- gsub("_Selleck","",SELLECK_MAP$product_id)

SELLECK_MAPPED <- merge(SELLECK, SELLECK_MAP, by.x = "Catalog.Number", by.y = "product_id") # 2661 x 23
SELLECK_withPubchemIDs <- subset(SELLECK_MAPPED, is.na(pubchem_id) == FALSE) # 884 have pubchem id; 1777 do not


length(intersect(as.integer(RES$pubchem_cid), SELLECK_MAPPED$pubchem_id)) # ALL: 99; FDA: 69

# incorporate SELLECK results with results so far
SELLECK_withPubchemIDs$pubchem_id <- as.character(SELLECK_withPubchemIDs$pubchem_id) 
names(SELLECK_withPubchemIDs) <- paste("selleck", names(SELLECK_withPubchemIDs), sep="_")
FIN <- merge(RES3, SELLECK_withPubchemIDs, by.x = "pubchem_cid", by.y = "selleck_pubchem_id")
# write.table(FIN, file="ctd2is_additional_cps_status_annot_pubchme-inchi_selleck.txt",
#             sep="\t", row.names = F, col.names = T, quote = F)

write.table(FIN, file="ctd2is_additional_cps_status_annot_pubchme-inchi_selleck-fda.txt",
            sep="\t", row.names = F, col.names = T, quote = F)

FINAL <- merge(RES6, SELLECK_withPubchemIDs, by.x = "pubchem_cid", by.y = "selleck_pubchem_id", all.x = T)
write.table(FINAL, file="ctd2is_all_cps_status_annot_pubchme-inchi_selleck-fda.txt",
            sep="\t", row.names = F, col.names = T, quote = F)

# counting up
REMAIN1 <- subset(FINAL, is.na(cpd_status) == TRUE) # after CTRP: 188
REMAIN2 <- subset(FINAL, is.na(cpd_status) == TRUE & is.na(DB_cpd_status) == TRUE) # after CTRP and DrugBank: 166 (gain 22)
REMAIN3 <- subset(FINAL, is.na(cpd_status) == TRUE & is.na(DB_cpd_status) == TRUE & is.na(selleck_STATUS) == TRUE) # after CTRP and DrugBank and Selleck:144 (gain 22)

REMAIN4 <- subset(FINAL, is.na(DB_cpd_status) == TRUE) # after DrugBank: 276 (gain 43)
REMAIN5 <- subset(FINAL, is.na(selleck_STATUS) == TRUE) # after Selleck: 249 (gain 70)

length(intersect(REMAIN1$pubchem_cid, as.character(SELLECK_MAPPED$pubchem_id))) # ALL=47; FDA=30 new labeled compounds
length(intersect(REMAIN1$pubchem_cid, as.character(RES3$pubchem_cid))) 

m_ctrp <- unique(subset(FINAL, is.na(cpd_status) == F)$pubchem_cid) # 130
m_db <- unique(subset(FINAL, is.na(DB_cpd_status) == F)$pubchem_cid) # 43
m_sell <- unique(subset(FINAL, is.na(selleck_STATUS) == F)$pubchem_cid) # 69
ctd2 <- unique(FINAL$pubchem_cid)

length(intersect(m_ctrp,m_db)) # 21
length(intersect(m_ctrp, m_sell)) # 39
length(intersect(m_db, m_sell)) # 19
length(Reduce(intersect, list(m_ctrp, m_db, m_sell))) # 11

ids <- list(CTRP=m_ctrp,
            DrugBank=m_db,
            SelleckFDA=m_sell,
            CTD2=ctd2)

png("integrated_cpd_clinical_status.png")
grid.draw(venn.diagram(ids, filename=NULL, 
                       fill = c("red", "green", "blue", "yellow"), 
                       alpha = c(0.5, 0.5,0.5, 0.5)))
dev.off()
