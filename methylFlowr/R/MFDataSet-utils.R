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

mfFilterBy <- function(obj, minComponentCoverage=NULL, minComponentWidth=NULL) {
    keep <- rep(TRUE, length(components(obj)))
    if (!is.null(minComponentCoverage)) {
        keep <- keep & (counts(obj, level="component") >= minComponentCoverage)
    }
    if (!is.null(minComponentWidth)) {
        keep <- keep & (width(components(obj)) >= minComponentWidth)
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
  if (by.component)
    return(table(regions(obj)$cid))
  nrow(regions(obj))
}

npatterns <- function(obj, by.component=TRUE) {
  if (by.component)
    return(table(patterns(obj)$cid))
  nrow(patterns(obj))
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


processMethylpats <- function(obj) {
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

ncpgs <- function(obj, level=c("region","pattern")) {
  level <- match.arg(level)
  if (level == "region") {
  }
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

positionCoverage <- function(obj){
  regionGR <- regions(obj)
  keep <- regionGR[regionsGR$exp_coverage >0]
  x <- IRanges(keep@ranges@start, keep@ranges@width)
  coverage(x)
  
}

componentAvgMeth <- function(obj) {
  regionGR <- regions(obj)
  cind <- split(seq(len=length(regionGR)), regionGR$cid)
  sapply(cind, function(ii) {
    if (sum(regionGR$ncpgs[ii]>0) == 0)
      return(NA)
    
    tab <- lapply(ii[regionGR$ncpgs[ii]>0], function(j) cbind(start(regionGR)[j]+regionGR$locs[[j]]-1, 1*(regionGR$meth[[j]]=="M"),regionGR$raw_coverage[j]))
    tab <- Reduce(rbind, tab)
    tab <- aggregate(tab[,3], list(tab[,1],tab[,2]),sum)
    
    mtab <- cbind(tab[,1],tab[,2]*tab[,3])
    mtab <- aggregate(mtab[,2],list(mtab[,1]),sum)

    covtab <- cbind(tab[,1],tab[,3])
    covtab <- aggregate(covtab[,2],list(covtab[,1]),sum)

    mean(mtab[,2] / covtab[,2])
  })
}
