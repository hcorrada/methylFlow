plotPatterns <- function(obj, plotRegion, ...) {
    keep <- components(obj) %over% plotRegion

    npats <- npatterns(obj, by.component=TRUE)[keep]
    cids <- components(obj)$cid[keep]

    npatterns <- max(npats)
    plot(0, type="n", xlim=c(start(plotRegion), end(plotRegion)+10), ylim=c(0,npatterns+1), ...)
    
    for (j in seq(along=cids)) {
        cid <- cids[j]
        patterns <- patterns(obj)[patterns(obj)$cid == cid,]

        o <- order(patterns$abundance)
        patterns <- patterns[o,]

        npatterns <- length(patterns)
        segments(min(start(patterns)), seq(along=patterns),
                 max(end(patterns)), seq(along=patterns))
	     
#        text(end(patterns), seq(along=patterns),
 #            labels=sprintf("%.3f", patterns$abundance), pos=4,cex=.8)

        scaledAbundance <- patterns$abundance / max(patterns$abundance)
        quantizedAbundance <- cut(scaledAbundance, breaks=seq(0,1,len=11))

        cols <- colorRampPalette(brewer.pal(5,"Blues"))(11)
        col <- cols[quantizedAbundance]

        rect(end(patterns)+10,seq(along=patterns)-.5,end(patterns) + 50,seq(along=patterns)+.5,
             col=col)
        text(end(patterns)+60, (npatterns+1)/2, label="abundance",cex=.7,srt=90)
#        points(end(patterns)+10, seq(along=patterns), pch=22, bg=col)
        
        locs <- rep(start(patterns), patterns$ncpgs) +
            unlist(patterns$locs)
        ylocs <- rep(seq(along=patterns), patterns$ncpgs)
        methyl <- unlist(patterns$meth)
        col  <- ifelse(methyl=="M", "black", "white")
        points(locs,ylocs,pch=21,bg=col)
    }
}
