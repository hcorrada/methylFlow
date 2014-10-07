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
        npatterns <- npatterns(obj, by.component=TRUE)
        cids <- names(npatterns)[npatterns >= minNumberOfPatterns]
        keep <- keep & (components(obj)$cid %in% cids)
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
      tab <- table(patterns(obj)$cid)
      m <- match(names(tab), components(obj)$cid)
    return(as.integer(tab[m]))
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

makeCpgGR <- function(obj, kind=c("raw", "estimated")) {
    kind <- match.arg(kind)
    if (kind == "raw") {
        gr <- regions(obj)
    } else {
        gr <- patterns(obj)
    }

    keep <- gr$ncpgs > 0
    gr <- gr[keep,]
    locs <- rep(start(gr), gr$ncpgs) + unlist(gr$locs)

    if (kind == "raw") {
        cov <- gr$raw_coverage
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
                                          
positionCoverage <- function(obj){
  regionGR <- regions(obj)
  keep <- regionGR[regionsGR$exp_coverage >0]
  x <- IRanges(keep@ranges@start, keep@ranges@width)
  coverage(x)
  
}

positionWeightedCoverage <- function(obj){
  regionGR <- regions(obj)
  keep <- regionGR[regionGR$exp_coverage >0]
  x <- IRanges(start = keep@ranges@start, width= keep@ranges@width)
  coverage(x, weight= keep$exp_coverage )
  }

componentAvgMeth <- function(obj) {
  #obj = objs2[[2]]
  regionGR <- regions(obj)
  cind <- split(seq(len=length(regionGR)), regionGR$cid)
  sapply(cind, function(ii) {
    if (sum(regionGR$ncpgs[ii]>0) == 0)
      return(NA)
    ii = cind[[100]]
    ii
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

regionMethPrecentage <- function(obj) {
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
    
    mtab[,2] <- (mtab[,2] / covtab[,2])
    mtab
  })
  
}

patternMethPrecentage <- function(obj) {
  #obj = objs2[[2]]
  patternGR <- patterns(obj)
  cind <- split(seq(len=length(patternGR)), patternGR$cid)
  sapply(cind, function(ii) {
    ii = cind[[100]]
    ii
    if (sum(patternGR$ncpgs[ii]>0) == 0)
      return(NA)
    
    tab <- lapply(ii[patternGR$ncpgs[ii]>0], function(j) cbind(start(patternGR)[j]+patternGR$locs[[j]]-1, 1*(patternGR$meth[[j]]=="M"),patternGR$abundance[j]))
    tab <- Reduce(rbind, tab)
    tab <- aggregate(tab[,3], list(tab[,1],tab[,2]),sum)
    
    mtab <- cbind(tab[,1],tab[,2]*tab[,3])
    mtab <- aggregate(mtab[,2],list(mtab[,1]),sum)
    
    covtab <- cbind(tab[,1],tab[,3])
    covtab <- aggregate(covtab[,2],list(covtab[,1]),sum)
    
    mtab[,2] <- (mtab[,2] / covtab[,2])
    mtab
  })
}

methPercentages2gr <- function(obj){
  #obj= objs[[2]]
  chr = levels(seqnames(obj@regions))
  obj2 <- processMethylpats(obj)
  rmp <- regionMethPrecentage(obj2)
  pmp <- patternMethPrecentage(obj2)
  keep <- width(components(obj)) > 100
  tabR <- Reduce(rbind, rmp[keep])
  tabP <- Reduce(rbind, pmp[keep])
  tabCombined = na.omit(merge(tabP,tabR,by='Group.1'))
  gr <- GRanges(seqnames=rep(chr , nrow(tabCombined)),
                ranges=IRanges(start=tabCombined[,1], width=1),
                readPercentage = tabCombined[,3],
                patternPercentage = tabCombined[,2])
  
  #gr <- GRanges(seqnames=rep(chr, length(tabR$x[!is.na(tabR$x)])),
   #             ranges=IRanges(start=tabR$Group.1[!is.na(tabR$Group.1)], width=1),
    #            readPercentage = tabR$x[!is.na(tabR$x)],
     #           patternPercentage = tabP$x[!is.na(tabR$x)])
  gr
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


