#!/usr/bin/env Rscript

### ./readLengthPlot.r par1 par2

### par1 = 0 > Auto-lambda
### par1 = 1 > Non-Auto - lambda is hard coded

### par2 = 0 > simple
### par2 = 1 > moderate
### par2 = 2 > Hard


data <- commandArgs(T)
print(data)
dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/"


##### reading files ##################
if (data[1] == "0"){
if ( data[2] == "2"){
    print("Hard Setting Plot")
    readLengthAvg <- read.table(paste(dir,"hard-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfReadLength <- read.table(paste(dir,"hard-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/hard-Auto/"
    
}

if ( data[2] == "1"){
    print("Moderate Setting Plot")
    readLengthAvg <- read.table(paste(dir,"moderate-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfReadLength <- read.table(paste(dir,"moderate-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/moderate-Auto/"
    
    
}
if ( data[2] == "0"){
    print("Simple Setting Plot")
    readLengthAvg <- read.table(paste(dir,"simple-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    mcfReadLength <- read.table(paste(dir,"simple-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
    
    dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/simple-Auto/"
    
    
}
}

if (data[1] == "1"){
    if ( data[2] == "2"){
        print("Hard Setting Plot")
        readLengthAvg <- read.table(paste(dir,"hard/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        mcfReadLength <- read.table(paste(dir,"hard/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/hard/"
        
    }
    
    if ( data[2] == "1"){
        print("Moderate Setting Plot")
        readLengthAvg <- read.table(paste(dir,"moderate/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        mcfReadLength <- read.table(paste(dir,"moderate/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/moderate/"
        
        
    }
    if ( data[2] == "0"){
        print("Simple Setting Plot")
        readLengthAvg <- read.table(paste(dir,"simple/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        mcfReadLength <- read.table(paste(dir,"simple/mcf.txt",sep=""), sep="\t", row.names=NULL, header = TRUE)
        
        dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/readLength/simple/"
        
        
    }
}






############## different plots for differnet read Length #################################################################


#### plot the abundance Error for different read Length

print("plot abundance Error vs read Length")
pdf(paste(dir,"abdVreadLength.pdf",sep=""))

# get the range for the x and y axis
xrange <- range(readLengthAvg[,1])
yrange <- range(readLengthAvg[,3])
ntrees <- length(unique(readLengthAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="readLength",
ylab="Abundance Error" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(readLengthAvg[,2]==j)
    lines(readLengthAvg[sel, 1],
    readLengthAvg[sel, 3],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(readLengthAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()


#### plot the Methyl Call  Error for different read Length

print("plot methyl call Error vs readLength ")
pdf(paste(dir,"methylVreadLength.pdf",sep=""))


# get the range for the x and y axis
xrange <- range(readLengthAvg[,1])
yrange <- range(readLengthAvg[,4])
ntrees <- length(unique(readLengthAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="readLength",
ylab="Methyl Call Error" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(readLengthAvg[,2]==j)
    lines(readLengthAvg[sel, 1],
    readLengthAvg[sel, 4],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(readLengthAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()

#### plot #TP for different read Length

print("plot TP vs readLength ")
pdf(paste(dir,"TPVreadLength.pdf",sep=""))


# get the range for the x and y axis
xrange <- range(readLengthAvg[,1])
yrange <- range(readLengthAvg[,5])
ntrees <- length(unique(readLengthAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="readLength",
ylab="TP" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(readLengthAvg[,2]==j)
    lines(readLengthAvg[sel, 1],
    readLengthAvg[sel, 5],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(readLengthAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()




######## plot sensitivity  for different readLength  ###############
#### TP / (TP + FN)


print("plot sensitivity vs readLength")
pdf(paste(dir,"sensitivityVreadLength.pdf",sep=""))

sensitivity = readLengthAvg[,5]/(readLengthAvg[,5] + readLengthAvg[,6])

# get the range for the x and y axis
xrange <- range(readLengthAvg[,1])
yrange <- range(sensitivity)

yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
yrange[2] <- yy

yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
yrange[1] <- yy

print(yrange[1])
print(yrange[2])

ntrees <- length(unique(readLengthAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="readLength",
ylab="sensitivity" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(readLengthAvg[,2]==j)
    lines(readLengthAvg[sel, 1],
    sensitivity[sel],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(readLengthAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()



######## plot Precision(Positive Predictive Rate)  for different readLength ###############
##### TP / (TP + FP )


print("plot precision vs readLength")
pdf(paste(dir,"precisionVreadLength.pdf",sep=""))

precision = readLengthAvg[,5]/(readLengthAvg[,5] + readLengthAvg[,7])

# get the range for the x and y axis
xrange <- range(readLengthAvg[,1])
yrange <- range(precision)
ntrees <- length(unique(readLengthAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="readLength",
ylab="precision" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(readLengthAvg[,2]==j)
    lines(readLengthAvg[sel, 1],
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

txt <- unique(readLengthAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()


######## plot FDR false discovery rate  for readLength ###############


print("plot false discovery rate vs readLength ")
pdf(paste(dir,"FDRVreadLength.pdf",sep=""))

FDR = readLengthAvg[,7]/(readLengthAvg[,5] + readLengthAvg[,7])

# get the range for the x and y axis
xrange <- range(readLengthAvg[,1])
yrange <- range(FDR)
ntrees <- length(unique(readLengthAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="readLength",
ylab="FDR" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(readLengthAvg[,2]==j)
    lines(readLengthAvg[sel, 1],
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

txt <- unique(readLengthAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()


####### plot min cost flow error for differnet read length ########


print("plot min cost flow error vs ReadLength")
pdf(paste(dir,"mcfReadLength.pdf",sep=""))

agg = aggregate(mcfReadLength[,2], list(readLength = mcfReadLength[,1]), FUN =  mean)

plot(agg)

dev.off()

