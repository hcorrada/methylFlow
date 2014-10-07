tab2gr <- function(tab) {
  gr <- GRanges(seqnames=tab$chr, ranges=IRanges(start=tab$start, end=tab$end))
  colsToKeep <- setdiff(colnames(tab), c("start","end","chr"))
  mcols(gr) <- DataFrame(tab[,colsToKeep])
  gr
}

read.methylflow.dir <- function(dir, sampleName, verbose = TRUE, has.header=TRUE) {
  stopifnot(file.exists(dir, "components.tsv") && file.exists(dir, "patterns.tsv") && file.exists(dir, "regions.tsv"))

  if (has.header) {
    componentsTab <- read.delim(file.path(dir, "components.tsv"), stringsAsFactors=FALSE)
    patternsTab <- read.delim(file.path(dir, "patterns.tsv"), stringsAsFactors=FALSE)
    regionsTab <- read.delim(file.path(dir, "regions.tsv"), stringsAsFactors=FALSE)
  } else {
    componentsTab <- read.delim(file.path(dir,"components.tsv"), stringsAsFactors=FALSE, header=FALSE)
    colnames(componentsTab) <- c("chr", "start","end","cid","npatterns","total_coverage","total_flow")

    patternsTab <- read.delim(file.path(dir, "patterns.tsv"), stringsAsFactors=FALSE, header=FALSE)
    colnames(patternsTab) <- c("chr", "start","end","cid","pid","abundance","methylpat")
    
    regionsTab <- read.delim(file.path(dir, "regions.tsv"), stringsAsFactors=FALSE, header=FALSE)
    colnames(regionsTab) <- c("chr","start","end","cid","rid","raw_coverage","norm_coverage","exp_coverage","methylpat")
  }

  componentsGR <- tab2gr(componentsTab)
  patternsGR <- tab2gr(patternsTab)
  regionsGR <- tab2gr(regionsTab)

  seqinfo <- seqinfo(componentsGR)
  new("MFDataSet", components=componentsGR, patterns=patternsGR, regions=regionsGR, seqinfo=seqinfo)
}
