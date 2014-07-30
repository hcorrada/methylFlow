#!/usr/bin/env Rscript

## run ./cpgPlot.r par1 par2

### par1 = 0 > Auto-lambda
### par1 = 1 > Non-Auto - lambda is hard coded

### par2 = 0 > simple
### par2 = 1 > moderate
### par2 = 2 > Hard

data <- commandArgs(T)
print(data)
dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/"


##### reading files ##################
if ( data[1] == 0){
    if ( data[2] == "2"){
        print("Hard Setting Plot")
        CpGAvg <- read.table(paste(dir,"hard-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
        
        mcfCpG <- read.table(paste(dir,"hard-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
        
        dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard-Auto/"
        
    }
    
    if ( data[2] == "1"){
        print("Moderate Setting Plot")
        CpGAvg <- read.table(paste(dir,"moderate-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
        
        mcfCpG <- read.table(paste(dir,"moderate-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
        
        dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate-Auto/"
        
        
    }
    if ( data[2] == "0"){
        print("Simple Setting Plot")
        CpGAvg <- read.table(paste(dir,"simple-Auto/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
        
        mcfCpG <- read.table(paste(dir,"simple-Auto/mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
        
        dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple-Auto/"
        
        
    }
}

if ( data[1] == 1){
if ( data[2] == "2"){
    print("Hard Setting Plot")
    CpGAvg <- read.table(paste(dir,"hard/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
    
    mcfCpG <- read.table(paste(dir,"hard/mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
    
    dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/hard/"

}

if ( data[2] == "1"){
    print("Moderate Setting Plot")
    CpGAvg <- read.table(paste(dir,"moderate/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
    
    mcfCpG <- read.table(paste(dir,"moderate/mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
    
    dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/moderate/"


}
if ( data[2] == "0"){
    print("Simple Setting Plot")
    CpGAvg <- read.table(paste(dir,"simple/evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
    
    mcfCpG <- read.table(paste(dir,"simple/mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)
    
    dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/cpg/simple/"


}
}




############## different plots for differnet # CPG #################################################################


#### plot the abundance Error for different number of CpG sites ############

print("plot abundance Error vs #CpG sites")
pdf(paste(dir,"abdVCpG.pdf",sep=""))

# get the range for the x and y axis
xrange <- range(CpGAvg[,1])
yrange <- range(CpGAvg[,3])
ntrees <- length(unique(CpGAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="CpG",
ylab="Abundance Error" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(CpGAvg[,2]==j)
    lines(CpGAvg[sel, 1],
    CpGAvg[sel, 3],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(CpGAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()


#### plot the Methyl Call  Error for different number of CpG sites #############

print("plot methyl call Error vs #CpG sites")
pdf(paste(dir,"methylVCpG.pdf",sep=""))


# get the range for the x and y axis
xrange <- range(CpGAvg[,1])
yrange <- range(CpGAvg[,4])
ntrees <- length(unique(CpGAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="CpG",
ylab="Methyl Call Error" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(CpGAvg[,2]==j)
    lines(CpGAvg[sel, 1],
    CpGAvg[sel, 4],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(CpGAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()

#### plot #TP for different number of CpG sites ##################

print("plot TP vs #CpG sites")
pdf(paste(dir,"TPVCpG.pdf",sep=""))


# get the range for the x and y axis
xrange <- range(CpGAvg[,1])
yrange <- range(CpGAvg[,5])
ntrees <- length(unique(CpGAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="CpG",
ylab="TP" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(CpGAvg[,2]==j)
    lines(CpGAvg[sel, 1],
    CpGAvg[sel, 5],
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(CpGAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()




######## plot sensitivity  for different number of CpG sites ###############
#### TP / (TP + FN)

print("plot sensitivity vs #CpG sites")
pdf(paste(dir,"sensitivityVCpG.pdf",sep=""))

sensitivity = CpGAvg[,5]/(CpGAvg[,5] + CpGAvg[,6])

# get the range for the x and y axis
xrange <- range(CpGAvg[,1])
yrange <- range(sensitivity)

yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
yrange[2] <- yy

yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
yrange[1] <- yy


ntrees <- length(unique(CpGAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="CpG",
ylab="sensitivity" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(CpGAvg[,2]==j)
    lines(CpGAvg[sel, 1],
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

txt <- unique(CpGAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()



######## plot Precision(Positive Predictive Rate)  for different number of CpG sites ###############
##### TP / (TP + FP )

print("plot precision vs #CpG sites")
pdf(paste(dir,"precisionVCpG.pdf",sep=""))

precision = rep(0, length(CpGAvg[,1]));

sel <- which(CpGAvg[,5] + CpGAvg[,7] != 0)
precision[sel] = CpGAvg[sel,5]/(CpGAvg[sel,5] + CpGAvg[sel,7])


# get the range for the x and y axis
xrange <- range(CpGAvg[,1])
yrange <- range(precision)
ntrees <- length(unique(CpGAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="CpG",
ylab="precision" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(CpGAvg[,2]==j)
    lines(CpGAvg[sel, 1],
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

txt <- unique(CpGAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()


######## plot FDR false discovery rate  for different number of CpG sites ###############


print("plot false discovery rate vs #CpG sites")
pdf(paste(dir,"FDRVCpG.pdf",sep=""))

FDR = CpGAvg[,7]/(CpGAvg[,5] + CpGAvg[,7])

# get the range for the x and y axis
xrange <- range(CpGAvg[,1])
yrange <- range(FDR)
ntrees <- length(unique(CpGAvg[,2]))
# set up the plot
plot(0, 0,
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = xrange,
xlab="CpG",
ylab="FDR" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    print(j)
    
    sel <- which(CpGAvg[,2]==j)
    lines(CpGAvg[sel, 1],
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

txt <- unique(CpGAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()



####### plot min cost flow error ########


print("plot min cost flow error vs #CpG sites")
pdf(paste(dir, "mcfCpG.pdf",sep=""))

agg = aggregate(mcfCpG[,2], list(numberofCpG = mcfCpG[,1]), FUN =  mean)

plot(agg)

dev.off()


