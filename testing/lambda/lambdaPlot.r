#!/usr/bin/env Rscript
data <- commandArgs(T)
print(data)
## data = 1 is in log scale
## data = 0 is in linear scale


dir <- "/cbcb/project-scratch/fdorri/Code/methylFlow/testing/lambda/"


    print("Hard Setting Plot")
    lambdaAvg <- read.table(paste(dir,"evalAvg.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)

    mcfLambda <- read.table(paste(dir,"mcf.txt",sep=""), sep="\t", row.names=NULL, header = FALSE)



if ( data == "1"){

############## different plots for differnet Lambda ################################################################


#### plot the abundance Error for different Lambda ############

print("plot abundance Error vs Lambda sites")
pdf(paste(dir,"abdVLambda.pdf",sep=""))

# get the range for the x and y axis
xrange <- range(lambdaAvg[,1])
yrange <- range(lambdaAvg[,3])
ntrees <- length(unique(lambdaAvg[,2]))
# set up the plot
plot(0, 0,
log ="x",
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
xlab="Lambda",
ylab="Abundance Error" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    #print(j)
    
    sel <- which(lambdaAvg[,2]==j)
    lines(lambdaAvg[sel, 1],
    lambdaAvg[sel, 3],
    log = "x",
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
log="x",
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(lambdaAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()


#### plot the Methyl Call  Error for different Lambda #############

print("plot methyl call Error vs Lambda")
pdf(paste(dir,"methylVLambda.pdf",sep=""))


# get the range for the x and y axis
xrange <- range(lambdaAvg[,1])
yrange <- range(lambdaAvg[,4])
ntrees <- length(unique(lambdaAvg[,2]))
# set up the plot
plot(0, 0,
log = "x",
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
xlab="Lambda",
ylab="Methyl Call Error" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    #print(j)
    
    sel <- which(lambdaAvg[,2]==j)
    lines(lambdaAvg[sel, 1],
    lambdaAvg[sel, 4],
    log = "x",
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
log = "x",
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(lambdaAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()

#### plot #TP for different Lambda ##################

print("plot TP vs Lambda")
pdf(paste(dir,"TPVLambda.pdf",sep=""))


# get the range for the x and y axis
xrange <- range(lambdaAvg[,1])
yrange <- range(lambdaAvg[,5])
ntrees <- length(unique(lambdaAvg[,2]))
# set up the plot
plot(0, 0,
log = "x",
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
xlab="Lambda",
ylab="TP" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    #print(j)
    sel <- which(lambdaAvg[,2]==j)
    lines(lambdaAvg[sel, 1],
    lambdaAvg[sel, 5],
    log = "x",
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
log = "x",
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(lambdaAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()




######## plot sensitivity  for different Lambda ###############
#### TP / (TP + FN)

print("plot sensitivity vs Lambda")
pdf(paste(dir,"sensitivityVLambda.pdf",sep=""))

sensitivity = lambdaAvg[,5]/(lambdaAvg[,5] + lambdaAvg[,6])

# get the range for the x and y axis
xrange <- range(lambdaAvg[,1])
yrange <- range(sensitivity)

yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
yrange[2] <- yy

yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
yrange[1] <- yy


ntrees <- length(unique(lambdaAvg[,2]))
# set up the plot
plot(0, 0,
log = "x",
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
xlab="Lambda",
ylab="sensitivity" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    #print(j)
    
    sel <- which(lambdaAvg[,2]==j)
    lines(lambdaAvg[sel, 1],
    sensitivity[sel],
    log = "x",
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
log = "x",
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(lambdaAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()



######## plot Precision(Positive Predictive Rate)  for different Lambda ###############
##### TP / (TP + FP )

print("plot precision vs Lambda")
pdf(paste(dir,"precisionVLambda.pdf",sep=""))

precision = rep(0, length(lambdaAvg[,1]));

sel <- which(lambdaAvg[,5] + lambdaAvg[,7] != 0)
precision[sel] = lambdaAvg[sel,5]/(lambdaAvg[sel,5] + lambdaAvg[sel,7])


# get the range for the x and y axis
xrange <- range(lambdaAvg[,1])
yrange <- range(precision)

yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
yrange[2] <- yy

yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
yrange[1] <- yy


ntrees <- length(unique(lambdaAvg[,2]))
# set up the plot
plot(0, 0,
log = "x",
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
xlab="Lambda",
ylab="precision" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    #print(j)
    
    sel <- which(lambdaAvg[,2]==j)
    lines(lambdaAvg[sel, 1],
    precision[sel],
    log = "x",
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
log = "x",
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(lambdaAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()


######## plot FDR false discovery rate  for different Lambda ###############


print("plot false discovery rate vs Lambda ")
pdf(paste(dir,"FDRVLambda.pdf",sep=""))

FDR = rep(0, length(lambdaAvg[,1]));

sel <- which(lambdaAvg[,5] + lambdaAvg[,7] != 0)
FDR[sel] = lambdaAvg[sel,7]/(lambdaAvg[sel,5] + lambdaAvg[sel,7])


# get the range for the x and y axis
xrange <- range(lambdaAvg[,1])
yrange <- range(FDR)

yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
yrange[2] <- yy

yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
yrange[1] <- yy


ntrees <- length(unique(lambdaAvg[,2]))
# set up the plot
plot(0, 0,
log = "x",
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
xlab="Lambda",
ylab="FDR" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    #print(j)
    
    sel <- which(lambdaAvg[,2]==j)
    lines(lambdaAvg[sel, 1],
    FDR[sel],
    log = "x",
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
log = "x",
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(lambdaAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()



####### plot min cost flow error for Lambda ########


print("plot min cost flow error vs Lambda sites")
pdf(paste(dir,"mcfLambda.pdf",sep=""))

agg = aggregate(mcfLambda[,2], list(Lambda = mcfLambda[,1]), FUN =  mean)
yrange <- range(agg$x)
xrange <- range(agg$Lambda)

yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
yrange[2] <- yy

yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
yrange[1] <- yy


ntrees <- length(unique(lambdaAvg[,2]))
# set up the plot
plot(0, 0,
log = "x",
pch = "",
ylim = c(yrange[1], yrange[2] + .1),
xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
xlab="Lambda",
ylab="min cost flow error" )

lines(agg$Lambda, agg$x,
log = "x",
pch = 16, cex=0.5)

##plot(agg$x, agg$Lambda, log='x')
dev.off()



#######

print("plot FP vs Lambda ")
pdf(paste(dir,"FPVLambda.pdf",sep=""))

FP = (lambdaAvg[,5]+ lambdaAvg[,7])

# get the range for the x and y axis
xrange <- range(lambdaAvg[,1])
yrange <- range(FP)
ntrees <- length(unique(lambdaAvg[,2]))
# set up the plot
plot(0, 0,
log = "x",
pch = "",
ylim = c(yrange[1], yrange[2] + 10),
xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
xlab="Lambda",
ylab="FP" )

colors <- rainbow(ntrees)
# add lines
for (i in 1:ntrees) {
    j = 0.02 * i + 0.02
    #print(j)
    
    sel <- which(lambdaAvg[,2]==j)
    lines(lambdaAvg[sel, 1],
    FP[sel],
    log = "x",
    col = colors[i])
}
# cex scale the size
#pch = 16 is circle
lx <- seq(30, 40, length.out=ntrees) + 5
ly <- rep(yrange[2] + .1, ntrees)
points(lx, ly,
log = "x",
col = colors[1:28],
pch = 16, cex=0.5)

txt <- unique(lambdaAvg[,2])
sel1 <- c(1, (1:7)*4)

text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)


dev.off()
}

if (data == 0 ){
    
    ############## different plots for differnet Lambda ################################################################
    
    
    #### plot the abundance Error for different Lambda ############
    
    print("plot abundance Error vs Lambda sites")
    pdf(paste(dir,"abdVLambda.pdf",sep=""))
    
    # get the range for the x and y axis
    xrange <- range(lambdaAvg[,1])
    yrange <- range(lambdaAvg[,3])
    ntrees <- length(unique(lambdaAvg[,2]))
    # set up the plot
    plot(0, 0,
    pch = "",
    ylim = c(yrange[1], yrange[2] + .1),
    xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
    xlab="Lambda",
    ylab="Abundance Error" )
    
    colors <- rainbow(ntrees)
    # add lines
    for (i in 1:ntrees) {
        j = 0.02 * i + 0.02
        #print(j)
        
        sel <- which(lambdaAvg[,2]==j)
        lines(lambdaAvg[sel, 1],
        lambdaAvg[sel, 3],
        col = colors[i])
    }
    # cex scale the size
    #pch = 16 is circle
    lx <- seq(30, 40, length.out=ntrees) + 5
    ly <- rep(yrange[2] + .1, ntrees)
    points(lx, ly,
    col = colors[1:28],
    pch = 16, cex=0.5)
    
    txt <- unique(lambdaAvg[,2])
    sel1 <- c(1, (1:7)*4)
    
    text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
    text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)
    
    
    dev.off()
    
    
    #### plot the Methyl Call  Error for different Lambda #############
    
    print("plot methyl call Error vs Lambda")
    pdf(paste(dir,"methylVLambda.pdf",sep=""))
    
    
    # get the range for the x and y axis
    xrange <- range(lambdaAvg[,1])
    yrange <- range(lambdaAvg[,4])
    ntrees <- length(unique(lambdaAvg[,2]))
    # set up the plot
    plot(0, 0,
    pch = "",
    ylim = c(yrange[1], yrange[2] + .1),
    xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
    xlab="Lambda",
    ylab="Methyl Call Error" )
    
    colors <- rainbow(ntrees)
    # add lines
    for (i in 1:ntrees) {
        j = 0.02 * i + 0.02
        #print(j)
        
        sel <- which(lambdaAvg[,2]==j)
        lines(lambdaAvg[sel, 1],
        lambdaAvg[sel, 4],
        col = colors[i])
    }
    # cex scale the size
    #pch = 16 is circle
    lx <- seq(30, 40, length.out=ntrees) + 5
    ly <- rep(yrange[2] + .1, ntrees)
    points(lx, ly,
    col = colors[1:28],
    pch = 16, cex=0.5)
    
    txt <- unique(lambdaAvg[,2])
    sel1 <- c(1, (1:7)*4)
    
    text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
    text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)
    
    
    dev.off()
    
    #### plot #TP for different Lambda ##################
    
    print("plot TP vs Lambda")
    pdf(paste(dir,"TPVLambda.pdf",sep=""))
    
    
    # get the range for the x and y axis
    xrange <- range(lambdaAvg[,1])
    yrange <- range(lambdaAvg[,5])
    ntrees <- length(unique(lambdaAvg[,2]))
    # set up the plot
    plot(0, 0,
    pch = "",
    ylim = c(yrange[1], yrange[2] + .1),
    xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
    xlab="Lambda",
    ylab="TP" )
    
    colors <- rainbow(ntrees)
    # add lines
    for (i in 1:ntrees) {
        j = 0.02 * i + 0.02
        #print(j)
        sel <- which(lambdaAvg[,2]==j)
        lines(lambdaAvg[sel, 1],
        lambdaAvg[sel, 5],
        col = colors[i])
    }
    # cex scale the size
    #pch = 16 is circle
    lx <- seq(30, 40, length.out=ntrees) + 5
    ly <- rep(yrange[2] + .1, ntrees)
    points(lx, ly,
    col = colors[1:28],
    pch = 16, cex=0.5)
    
    txt <- unique(lambdaAvg[,2])
    sel1 <- c(1, (1:7)*4)
    
    text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
    text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)
    
    
    dev.off()
    
    
    
    
    ######## plot sensitivity  for different Lambda ###############
    #### TP / (TP + FN)
    
    print("plot sensitivity vs Lambda")
    pdf(paste(dir,"sensitivityVLambda.pdf",sep=""))
    
    sensitivity = lambdaAvg[,5]/(lambdaAvg[,5] + lambdaAvg[,6])
    
    # get the range for the x and y axis
    xrange <- range(lambdaAvg[,1])
    yrange <- range(sensitivity)
    
    yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
    yrange[2] <- yy
    
    yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
    yrange[1] <- yy
    
    
    ntrees <- length(unique(lambdaAvg[,2]))
    # set up the plot
    plot(0, 0,
    pch = "",
    ylim = c(yrange[1], yrange[2] + .1),
    xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
    xlab="Lambda",
    ylab="sensitivity" )
    
    colors <- rainbow(ntrees)
    # add lines
    for (i in 1:ntrees) {
        j = 0.02 * i + 0.02
        #print(j)
        
        sel <- which(lambdaAvg[,2]==j)
        lines(lambdaAvg[sel, 1],
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
    
    txt <- unique(lambdaAvg[,2])
    sel1 <- c(1, (1:7)*4)
    
    text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
    text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)
    
    
    dev.off()
    
    
    
    ######## plot Precision(Positive Predictive Rate)  for different Lambda ###############
    ##### TP / (TP + FP )
    
    print("plot precision vs Lambda")
    pdf(paste(dir,"precisionVLambda.pdf",sep=""))
    
    precision = rep(0, length(lambdaAvg[,1]));
    
    sel <- which(lambdaAvg[,5] + lambdaAvg[,7] != 0)
    precision[sel] = lambdaAvg[sel,5]/(lambdaAvg[sel,5] + lambdaAvg[sel,7])
    
    
    # get the range for the x and y axis
    xrange <- range(lambdaAvg[,1])
    yrange <- range(precision)
    
    yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
    yrange[2] <- yy
    
    yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
    yrange[1] <- yy
    
    
    ntrees <- length(unique(lambdaAvg[,2]))
    # set up the plot
    plot(0, 0,
    pch = "",
    ylim = c(yrange[1], yrange[2] + .1),
    xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
    xlab="Lambda",
    ylab="precision" )
    
    colors <- rainbow(ntrees)
    # add lines
    for (i in 1:ntrees) {
        j = 0.02 * i + 0.02
        #print(j)
        
        sel <- which(lambdaAvg[,2]==j)
        lines(lambdaAvg[sel, 1],
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
    
    txt <- unique(lambdaAvg[,2])
    sel1 <- c(1, (1:7)*4)
    
    text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
    text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)
    
    
    dev.off()
    
    
    ######## plot FDR false discovery rate  for different Lambda ###############
    
    
    print("plot false discovery rate vs Lambda ")
    pdf(paste(dir,"FDRVLambda.pdf",sep=""))
    
    FDR = rep(0, length(lambdaAvg[,1]));
    
    sel <- which(lambdaAvg[,5] + lambdaAvg[,7] != 0)
    FDR[sel] = lambdaAvg[sel,7]/(lambdaAvg[sel,5] + lambdaAvg[sel,7])
    
    
    # get the range for the x and y axis
    xrange <- range(lambdaAvg[,1])
    yrange <- range(FDR)
    
    yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
    yrange[2] <- yy
    
    yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
    yrange[1] <- yy
    
    
    ntrees <- length(unique(lambdaAvg[,2]))
    # set up the plot
    plot(0, 0,
    pch = "",
    ylim = c(yrange[1], yrange[2] + .1),
    xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
    xlab="Lambda",
    ylab="FDR" )
    
    colors <- rainbow(ntrees)
    # add lines
    for (i in 1:ntrees) {
        j = 0.02 * i + 0.02
        #print(j)
        
        sel <- which(lambdaAvg[,2]==j)
        lines(lambdaAvg[sel, 1],
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
    
    txt <- unique(lambdaAvg[,2])
    sel1 <- c(1, (1:7)*4)
    
    text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
    text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)
    
    
    dev.off()
    
    
    
    ####### plot min cost flow error for Lambda ########
    
    
    print("plot min cost flow error vs Lambda sites")
    pdf(paste(dir,"mcfLambda.pdf",sep=""))
    
    agg = aggregate(mcfLambda[,2], list(Lambda = mcfLambda[,1]), FUN =  mean)
    yrange <- range(agg$x)
    xrange <- range(agg$Lambda)
    
    yy  <- ifelse(is.na(yrange[2]) , 0, yrange[2])
    yrange[2] <- yy
    
    yy  <- ifelse(is.na(yrange[1]) , 0, yrange[1])
    yrange[1] <- yy
    
    
    ntrees <- length(unique(lambdaAvg[,2]))
    # set up the plot
    plot(0, 0,
    pch = "",
    ylim = c(yrange[1], yrange[2] + .1),
    xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
    xlab="Lambda",
    ylab="min cost flow error" )
    
    lines(agg$Lambda, agg$x,
    pch = 16, cex=0.5)
    
    ##plot(agg$x, agg$Lambda, log='x')
    dev.off()
    
    
    
    #######
    
    print("plot FP vs Lambda ")
    pdf(paste(dir,"FPVLambda.pdf",sep=""))
    
    FP = (lambdaAvg[,5]+ lambdaAvg[,7])
    
    # get the range for the x and y axis
    xrange <- range(lambdaAvg[,1])
    yrange <- range(FP)
    ntrees <- length(unique(lambdaAvg[,2]))
    # set up the plot
    plot(0, 0,
    pch = "",
    ylim = c(yrange[1], yrange[2] + 10),
    xlim = c(xrange[1]+0.00001,xrange[2]+0.0001),
    xlab="Lambda",
    ylab="FP" )
    
    colors <- rainbow(ntrees)
    # add lines
    for (i in 1:ntrees) {
        j = 0.02 * i + 0.02
        #print(j)
        
        sel <- which(lambdaAvg[,2]==j)
        lines(lambdaAvg[sel, 1],
        FP[sel],
        col = colors[i])
    }
    # cex scale the size
    #pch = 16 is circle
    lx <- seq(30, 40, length.out=ntrees) + 5
    ly <- rep(yrange[2] + .1, ntrees)
    points(lx, ly,
    col = colors[1:28],
    pch = 16, cex=0.5)
    
    txt <- unique(lambdaAvg[,2])
    sel1 <- c(1, (1:7)*4)
    
    text(xrange[1], yrange[2]+.1, "Legend", cex=0.5, pos=4)
    text(lx[sel1], ly[sel1] - 0.02, txt[sel1], cex=0.5, srt=90)
    
    
    dev.off()
}
else{
    print( "enter 0 for linear 1 for log")
}

















