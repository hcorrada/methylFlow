tab2gr <- function(tab, onlypos=FALSE) {
  if (onlypos) {
    gr <- GRanges(seqnames=tab$chr, ranges=IRanges(start=tab$pos, end=tab$pos))
    colsToKeep <- setdiff(colnames(tab), c("chr","pos"))
  } else {
    gr <- GRanges(seqnames=tab$chr, ranges=IRanges(start=tab$start, end=tab$end))      
    colsToKeep <- setdiff(colnames(tab), c("start","end","chr"))
  }
  mcols(gr) <- DataFrame(tab[,colsToKeep])
  gr
}

read.methylflow.dir <- function(dir, sampleName, verbose = TRUE, has.header=TRUE) {
  stopifnot(file.exists(dir, "components.tsv") && 
              file.exists(dir, "patterns.tsv") && 
              file.exists(dir, "regions.tsv") &&
              file.exists(dir, "cpgs.tsv"))

  if (has.header) {
    componentsTab <- read.delim(file.path(dir, "components.tsv"), stringsAsFactors=FALSE)
    patternsTab <- read.delim(file.path(dir, "patterns.tsv"), stringsAsFactors=FALSE)
    regionsTab <- read.delim(file.path(dir, "regions.tsv"), stringsAsFactors=FALSE)
    cpgTab <- read.delim(file.path(dir, "cpgs.tsv"), stringsAsFactors=FALSE)
  } else {
    componentsTab <- read.delim(file.path(dir,"components.tsv"), stringsAsFactors=FALSE, header=FALSE)
    colnames(componentsTab) <- c("chr", "start","end","cid","npatterns","total_coverage","total_flow")

    patternsTab <- read.delim(file.path(dir, "patterns.tsv"), stringsAsFactors=FALSE, header=FALSE)
    colnames(patternsTab) <- c("chr", "start","end","cid","pid","abundance","methylpat")
    
    regionsTab <- read.delim(file.path(dir, "regions.tsv"), stringsAsFactors=FALSE, header=FALSE)
    colnames(regionsTab) <- c("chr","start","end","cid","rid","raw_coverage","norm_coverage","exp_coverage","methylpat")
    
    cpgTab <- read.delim(file.path(dir, "cpgs.tsv"), stringsAsFactors=FALSE, header=FALSE)
    colnames(cpgTab) <- c("chr", "pos", "Cov", "Meth")
  }

  componentsGR <- tab2gr(componentsTab)
  patternsGR <- tab2gr(patternsTab)
  regionsGR <- tab2gr(regionsTab)
  cpgGR <- tab2gr(cpgTab, onlypos=TRUE)
  
  seqinfo <- seqinfo(componentsGR)
  new("MFDataSet", components=componentsGR, patterns=patternsGR, regions=regionsGR, cpgs=cpgGR, seqinfo=seqinfo)
}
