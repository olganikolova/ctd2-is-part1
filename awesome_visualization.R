# visualizing missing data and
# other issues in a data frame
# came from here:
# http://stackoverflow.com/questions/27545423/visual-structure-of-a-data-frame-locations-of-nas-and-much-more

# PACKAGES 
require(ggplot2)
require(RColorBrewer)
require(reshape2)

# Test if an object is empty (data.frame, matrix, vector)
is.empty = function (input) {
  df <- data.frame(input)
  (is.null(df) || nrow(df) == 0 || ncol(df) == 0 || NROW(df) == 0)
}

#  min/max normalization (R->[0;1]), (all columns must be numerical)
minmax <- function(data, ...) {
  .minmax = function(x) (x-min(x, ...))/(max(x, ...)-min(x, ...))
  # find constant columns, replaces with O.5:
  constant <- which(apply(data, 2, function(u) {min(u, ...)==max(u, ...)}))
  if(is.vector(data)) {
    res <- .minmax(data)
  } else {
    res <- apply(data, 2, .minmax)
  }
  res[, constant] <- 0.5
  return(res)
}

# MAIN function
colstr = function(input, size.max=500, export=FALSE) {
  data      <- as.data.frame(input)
  if (NCOL(data) == 1) {
    data    <- cbind(data, data)
    message("warning: input data is a vector")
  }
  miror     <- data # miror data.frame will contain a coulour coding for all cells
  wholeNA   <- which(sapply(miror, function(x) all(is.na(x))))
  whole0    <- which(sapply(miror, function(x) all(x==0)))
  numeric   <- which(sapply(data, is.numeric))
  character <- which(sapply(data, is.character))
  factor    <- which(sapply(data, is.factor))
  # characters to code
  miror[character] <- 12 
  # factor coding
  miror[factor] <- 11
  # min/max normalization, coerce it into 9 classes.
  if (!is.empty(numeric)) {miror[numeric] <- minmax(miror[numeric], na.rm=T)}
  miror[numeric] <- data.frame(lapply(miror[numeric], function(x) cut(x, breaks=9, labels=1:9))) # 9 classes numériques
  miror <- data.frame(lapply(miror, as.numeric))
  # Na coding
  miror[is.na(data)] <- 10
  miror[whole0]    <- 13
  # color palette vector
  mypalette <- c(brewer.pal(n=9, name="Blues"), "red", "green", "purple", "grey")
  colnames <- c(paste0((1:9)*10, "%"), "NA", "factor (lvls)", "character", "zero")
  # subset if too large
  couper <- nrow(miror) > size.max
  if (couper) miror <- head(miror, size.max)
  # plot
  g <- ggplot(data=melt(as.matrix(unname(miror)))) + 
    geom_tile(aes(x=Var2, y=Var1, fill=factor(value, levels=1:13))) +
    scale_fill_manual("legend", values=mypalette, labels=colnames, drop=FALSE) +
    ggtitle(paste("graphical structure of", deparse(substitute(input)), paste(dim(input), collapse="X"), ifelse(couper, "(truncated)", ""))) +
    xlab("columns of the dataframe") + ylab("rows of the dataframe") +
    geom_point(data=data.frame(x=0, y=1:NROW(input)), aes(x,y), alpha=1-all(row.names(input)==seq(1, NROW(input)))) +
    scale_y_reverse(limits=c(min(size.max, nrow(miror)), 0))
  if (!is.empty(factor)) {
    g <- g + geom_text(data=data.frame(x     = factor, 
                                       y     = round(runif(length(factor), 2, NROW(miror)-2)), 
                                       label = paste0("(", sapply(data[factor], function(x) length(levels(x))), ")")),
                       aes(x=x, y=y, label=label))
  }
  if (export) {png("colstr_output.png"); print(g); dev.off()}
  return(g)
}