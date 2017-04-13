## aux function to remove leading zeros ##
rm_lead_zeros <- function(STR=NULL)
{
  if(is.null(STR)){
    stop("null argument in rm_lead_zeros...\n")
  }
  return(gsub("(?<![0-9])0+", "", STR, perl = TRUE))
}