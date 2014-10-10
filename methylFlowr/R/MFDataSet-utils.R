regions <- function(obj) obj@regions
patterns <- function(obj) obj@patterns
components <- function(obj) obj@components

setMethod("show", "MFDataSet", function(object) {
    cat("MFDataSet object:\n")
    cat("  ", length(components(object)), " components\n")
    cat("  ", length(patterns(object)), " patterns\n")
    cat("  ", length(regions(object)), " regions\n")
    show(seqinfo(object))
})

setMethod("seqinfo", "MFDataSet", function(x) x@seqinfo)
setMethod("seqinfo<-", "MFDataSet", function(x, new2old = NULL, force=FALSE, value) {
    x@seqinfo <- value
    x@regions <- callGeneric(regions(x), new2old, force, value)
    x@patterns <- callGeneric(patterns(x), new2old, force, value)
    x@components <- callGeneric(components(x), new2old, force, value)
    x
})

mfFilterBy <- function(obj,
                       minComponentCoverage=NULL,
                       minComponentWidth=NULL,
                       minNumberOfPatterns=NULL) {
    keep <- rep(TRUE, length(components(obj)))
    if (!is.null(minComponentCoverage)) {
        keep <- keep & (counts(obj, level="component") >= minComponentCoverage)
    }
    if (!is.null(minComponentWidth)) {
        keep <- keep & (width(components(obj)) >= minComponentWidth)
    }
    if (!is.null(minNumberOfPatterns)) {
        keep <- keep & (npatterns(obj, by.component=TRUE) >= minNumberOfPatterns)
    }
    componentIdsToKeep <- components(obj)$cid[keep]
    obj@components <- components(obj)[keep,]

    regionsToKeep <- regions(obj)$cid %in% componentIdsToKeep
    obj@regions <- regions(obj)[regionsToKeep,]

    patternsToKeep <- patterns(obj)$cid %in% componentIdsToKeep
    obj@patterns <- patterns(obj)[patternsToKeep,]
    obj
}

nregions <- function(obj, by.component=TRUE) {
    if (isTRUE(by.component)) {
        tab <- table(regions(obj)$cid)
        m <- match(names(tab), components(obj)$cid)
        return(as.integer(tab[m]))
    }
    length(regions(obj))
}

npatterns <- function(obj, by.component=TRUE) {
  if (isTRUE(by.component)) {
      return(components(obj)$npatterns)
  }
  length(patterns(obj))
}

counts <- function(obj, level=c("region","component"), kind=c("raw","normalized"))
  {
    level <- match.arg(level)
    kind <- match.arg(kind)
    if (level == "region") {
      return(switch(kind,
                    raw=regions(obj)$raw_coverage,
                    normalized=regions(obj)$norm_coverage))
    } else {
      return(switch(kind,
                    raw=components(obj)$total_coverage,
                    normalized=components(obj)$total_coverage))
    }
  }


.parseMethylpats <- function(x) {
    tmp <- strsplit(x, ",")
    ncpgs <- ifelse(x=="*", 0, sapply(tmp,length))
    tmp2 <- lapply(seq(along=x), function(i) {
        if (ncpgs[i] == 0) return(NULL)
        strsplit(tmp[[i]], ":")
    })
    
    locs <- lapply(tmp2, function(y) {
      if (is.null(y)) return(0)
      as.integer(sapply(y,"[",1))
    })

    meth <- lapply(tmp2, function(y) {
        if (is.null(y)) return("")
        sapply(y,"[",2)
    })
    list(ncpgs=ncpgs,locs=locs,meth=meth)
}

processMethylpats <- function(obj) {
  patternMethylPats <- .parseMethylpats(patterns(obj)$methylpat)
  obj@patterns$ncpgs <- patternMethylPats$ncpgs
  obj@patterns$locs <- patternMethylPats$locs
  obj@patterns$meth <- patternMethylPats$meth

  regionMethylPats <- .parseMethylpats(regions(obj)$methylpat)
  obj@regions$ncpgs <- regionMethylPats$ncpgs
  obj@regions$locs <- regionMethylPats$locs
  obj@regions$meth <- regionMethylPats$meth
  obj
}

ncpgs <- function(obj, level=c("region","pattern"), summary=c("none", "median", "min", "max")) {
  level <- match.arg(level)
  summary <- match.arg(summary)
  
  if (level == "region") {
      if (!is.null(regions(obj)$ncpgs)) {
          res <- regions(obj)$ncpgs
          cids <- regions(obj)$cid
      } else {
          return(NULL)
      }
  } else if (level == "pattern") {
      if (!is.null(patterns(obj)$ncpgs)) {
          res <- patterns(obj)$ncpgs
          cids <- patterns(obj)$cid
      } else {
          return(NULL)
      }
  }

  if (summary != "none") {
      fn <- get(summary)
      res <- tapply(res, cids, fn)
  }
  res
}

makeCpgGR <- function(obj, kind=c("raw", "normalized", "estimated")) {
    kind <- match.arg(kind)
    if (kind %in% c("raw","normalized")) {
        gr <- regions(obj)
    } else {
        gr <- patterns(obj)
    }

    keep <- gr$ncpgs > 0
    gr <- gr[keep,]
    locs <- rep(start(gr), gr$ncpgs) + unlist(gr$locs)

    if (kind == "raw") {
        cov <- gr$raw_coverage / width(gr) * 10
    } else if (kind == "normalized") {
        cov <- gr$norm_coverage / width(gr) * 10
    } else {
        cov <- gr$abundance
    }

    Cov <- rep(cov, gr$ncpgs)
    Meth <- 1*(unlist(gr$meth) == "M") * Cov

    Cov <- tapply(Cov, locs, sum)
    Meth <- tapply(Meth, locs, sum)
    Loc <- tapply(locs,locs, function(x) x[1])
    chr <- tapply(rep(as.character(seqnames(gr)), gr$ncpgs), locs, function(x) x[1])
    newGR <- GRanges(chr, IRanges(start=Loc, width=1),
                       Cov=Cov, Meth=Meth)
    newGR$Beta <- Meth / Cov
    names(newGR) <- NULL
    seqinfo(newGR) <- seqinfo(gr)
    newGR
}

getEntropyStats <- function(obj) {
    npatcov <- coverage(patterns(obj))

    covRegions <- as(npatcov, "GRanges")
    covRegions <- covRegions[covRegions$score>0,]

    olaps <- findOverlaps(patterns(obj), covRegions)

    newPatterns <- covRegions[subjectHits(olaps),]
    newPatterns$rid <- subjectHits(olaps)
    newPatterns$abundance <- patterns(obj)$abundance[queryHits(olaps)]
    totalFlow <- tapply(newPatterns$abundance, newPatterns$rid, sum)
    newPatterns$prop <- newPatterns$abundance / totalFlow[as.character(newPatterns$rid)]

    covRegions$gini <- 1 - tapply(newPatterns$prop^2, newPatterns$rid, sum)
    covRegions$entr <- tapply(-newPatterns$prop * log2(newPatterns$prop), newPatterns$rid, sum)
    covRegions$maxEntr <- log2(tapply(rep(1, length(newPatterns)), newPatterns$rid, sum))
    covRegions$normEntr <- covRegions$entr / covRegions$maxEntr
    seqinfo(covRegions) <- seqinfo(obj)
    covRegions
}

positionCoverage <- function(obj, kind=c("raw", "estimated")){
    kind <- match.arg(kind)
    if (kind == "raw") {
        gr <- regions(obj)
    } else {
        gr <- patterns(obj)
    }
    coverage(gr)
}

positionWeightedCoverage <- function(obj){
  regionGR <- regions(obj)
  keep <- regionGR[regionGR$exp_coverage >0]
  x <- IRanges(start = keep@ranges@start, width= keep@ranges@width)
  coverage(x, weight= keep$exp_coverage )
  }


componentEntropy <- function(obj) {
  pats <- patterns(obj)
  .ent <- function(x) {
    p <- x/sum(x)
    -sum(p*log(p))
   # sum(p*p)
  }
  tapply(pats$abundance,pats$cid,  .ent)
 
}

#(this doesn't work)
componentEpipolymorphism <- function(obj) {
 pats <- patterns(obj)
  .norm <- function(x) {
     zoo(x/sum(x))
    #   -sum(p*log(p))
    #sum(p*p)
  }
  z <- tapply(pats$abundance,pats$cid, .norm)
  lapply(z , function(x) {1 - rollapply(x, width=4, function(i) sum(i*i))})  
}


#(this dosn't work)
posistionEntropy <- function(obj){
  regionGR <- regions(obj)
  keep <- regionGR[regionGR$exp_coverage >0]
  pwc <- positionWeightedCoverage(obj)
  makerle <- function(gr){
    x<- IRanges(start= gr@ranges@start, width = gr@ranges@width)
    p <- coverage(x, weight= gr@exp_coverage)/pwc
    -p*log(p)
    }
 lapply(keep, makerle)
}


