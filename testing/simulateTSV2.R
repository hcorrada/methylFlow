nRegions <- 2

theLocs <- c(2,6,11,22,25,32,35,38,45)
regionWidth <- 50

nPatterns <- 3
nLocs <- length(theLocs)
cpgMat <- matrix("", nr=nPatterns, nc=nLocs)
rownames(cpgMat) <- c("red","green","blue")

cpgMat["red",] <- rep(c("black","white","black"), c(5,3,1))
cpgMat["green",] <- rep(c("black","white","black"), c(3,2,4))
cpgMat["blue",] <- rep(c("black"), c(9))

fragments <- c(2,4,6)
nFragments <- sum(fragments)

readsPerFragment <- 20
rlen <- 20

set.seed(1)
# sample starting points
readLocs <- sample(1:(regionWidth-rlen),size=readsPerFragment*nFragments,replace=TRUE)

# sample frag
frag <- sample(size=readsPerFragment*nFragments, x=rownames(cpgMat), p=fragments/nFragments, replace=TRUE)

# sample component
offset <- sample(c(0,100), length(readLocs), replace=TRUE)
o <- order(readLocs + offset,frag)
readLocs <- readLocs[o]
frag <- frag[o]
offset <- offset[o]


outdf <- data.frame(readid=paste0('aread',seq(along=readLocs)),
                    loc=readLocs + offset,
                    rlen=rlen,
                    strand='W')

cpgstring <- rep('*', length(readLocs))

for (i in seq(along=readLocs)) {
  curLoc <- readLocs[i]
  locOverlap <- theLocs >= curLoc & theLocs <= curLoc+rlen-1

  if (any(locOverlap)) {
    locOverlap <- which(locOverlap)
    cpgLoc <- theLocs[locOverlap]
    readLoc <- cpgLoc - curLoc + 1
    meth <- ifelse(cpgMat[frag[i], locOverlap] == 'white', 'U', 'M')
    cpgstring[i] <- paste(readLoc,meth,sep=":",collapse=',')
  }
}

outdf <- cbind(outdf, cpgstring=cpgstring, substring='*')
write.table(outdf,sep='\t', file='sim2.tsv',quote=FALSE,col.names=FALSE,row.names=FALSE)




