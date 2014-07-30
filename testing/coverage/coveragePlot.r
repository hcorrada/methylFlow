#!/usr/bin/env Rscript

### run ./coveragePlot.r par1 par2

### par1 = 0 > Auto-lambda
### par1 = 1 > Non-Auto - lambda is hard coded

### par2 = 0 > simple
### par2 = 1 > moderate
### par2 = 2 > Hard


data <- commandArgs(T)
print(data)
dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/"


##### reading files ##################
if (data[1] == "0"){
    if ( data[2] == "2"){
        print("Hard Setting Plot")
        coverageAvg <- read.table(paste(dir,"hard-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
        
        mcfCoverage <- read.table(paste(dir,"hard-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
        
        dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard-Auto/"
        
    }
    
    if ( data[2] == "1"){
        print("Moderate Setting Plot")
        coverageAvg <- read.table(paste(dir,"moderate-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
        
        mcfCoverage <- read.table(paste(dir,"moderate-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
        
        dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate-Auto/"
        
        
    }
    if ( data[2] == "0"){
        print("Simple Setting Plot")
        coverageAvg <- read.table(paste(dir,"simple-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
        
        mcfCoverage <- read.table(paste(dir,"simple-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
        
        dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple-Auto/"
        
        
    }
}

if (data[1] == "1"){
if ( data[2] == "2"){
    print("Hard Setting Plot")
    coverageAvg <- read.table(paste(dir,"hard/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
    
    mcfCoverage <- read.table(paste(dir,"hard/mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
    
    dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/hard/"
    
}

if ( data[2] == "1"){
    print("Moderate Setting Plot")
    coverageAvg <- read.table(paste(dir,"moderate/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
    
    mcfCoverage <- read.table(paste(dir,"moderate/mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
    
    dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/moderate/"
    
    
}
if ( data[2] == "0"){
    print("Simple Setting Plot")
    coverageAvg <- read.table(paste(dir,"simple/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
    
    mcfCoverage <- read.table(paste(dir,"simple/mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
    
    dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/coverage/simple/"
    
    
}
}


############## different plots for differnet coverage rate #################################################################


#### plot the abundance Error for different coverage rate

print("plot abundance Error vs coverage rate")
pdf(paste(dir,"abdVcoverage.pdf",sep=""))

# get the range for the x and y axis
xrange <- range(coverageAvg[,1])
yrange <- range(coverageAvg[,3])
ntrees <- length(unique(coverageAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="coverage",
ylab="Abundance Error" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(coverageAvg[,2]==j)
    lines(coverageAvg[sel, 1],
    coverageAvg[sel, 3],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(coverageAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()


#### plot the Methyl Call  Error for different coverage rate

print("plot methyl call Error vs coverage rate ")
pdf(paste(dir, "methylVcoverage.pdf",sep=""))


# get the range for the x and y axis
xrange <- range(coverageAvg[,1])
yrange <- range(coverageAvg[,4])
ntrees <- length(unique(coverageAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="coverage",
ylab="Methyl Call Error" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(coverageAvg[,2]==j)
    lines(coverageAvg[sel, 1],
    coverageAvg[sel, 4],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(coverageAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()

#### plot #TP for different coverage rate

print("plot TP vs coverage rate ")
pdf(paste(dir, "TPVcoverage.pdf",sep=""))


# get the range for the x and y axis
xrange <- range(coverageAvg[,1])
yrange <- range(coverageAvg[,5])
ntrees <- length(unique(coverageAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="coverage",
ylab="TP" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(coverageAvg[,2]==j)
    lines(coverageAvg[sel, 1],
    coverageAvg[sel, 5],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(coverageAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()




######## plot sensitivity  for different coverage  ###############
#### TP / (TP + FN)


print("plot sensitivity vs coverage")
pdf(paste(dir,"sensitivityVcoverage.pdf",sep=""))

sensitivity = coverageAvg[,5]/(coverageAvg[,5] + coverageAvg[,6])

# get the range for the x and y axis
xrange <- range(coverageAvg[,1])
yrange <- range(sensitivity)

yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
yrange[2] <- yy

yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
yrange[1] <- yy

ntrees <- length(unique(coverageAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="coverage",
ylab="sensitivity" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(coverageAvg[,2]==j)
    lines(coverageAvg[sel, 1],
    sensitivity[sel],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(xrange[1] + 1, xrange[1]+10, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(coverageAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()



######## plot Precision(Positive Predictive Rate)  for different Coverage ###############
##### TP / (TP + FP )


print("plot precision vs coverage")
pdf(paste(dir,"precisionVcoverage.pdf",sep=""))

precision = coverageAvg[,5]/(coverageAvg[,5] + coverageAvg[,7])

# get the range for the x and y axis
xrange <- range(coverageAvg[,1])
yrange <- range(precision)
ntrees <- length(unique(coverageAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="coverage",
ylab="precision" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(coverageAvg[,2]==j)
    lines(coverageAvg[sel, 1],
    precision[sel],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(coverageAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()


######## plot FDR false discovery rate  for Coverage ###############


print("plot false discovery rate vs coverage ")
pdf(paste(dir,"FDRVCoverage.pdf",sep=""))

FDR = coverageAvg[,7]/(coverageAvg[,5] + coverageAvg[,7])

# get the range for the x and y axis
xrange <- range(coverageAvg[,1])
yrange <- range(FDR)
ntrees <- length(unique(coverageAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="coverage",
ylab="FDR" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(coverageAvg[,2]==j)
    lines(coverageAvg[sel, 1],
    FDR[sel],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(coverageAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()


####### plot min cost flow error for differnet coverage ########


print("plot min cost flow error vs Coverage")
pdf(paste(dir,"mcfCoverage.pdf",sep=""))

agg = aggregate(mcfCoverage[,2], list(coverage = mcfCoverage[,1]), FUN =  mean)

plot(agg)

dev.off()


