## aux function to clean pubchem ids ##
source("/fh/fast/kemp_c/olga_files/projects/informer_set/code/rm_lead_zeros.R")
get_numeric_pubchem_id <- function(STR=NULL)
{
  if(is.null(STR)){
    stop("null argument in get_numeric_pubchem_id...\n")
  }
  TMP <- gsub("CID[m|s]*:*", "", STR)
  return(rm_lead_zeros(TMP))
}