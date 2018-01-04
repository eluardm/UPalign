#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)

mat = read.table(paste(args[1],'plotTMscore.dat', sep=""), head=T)
svg(filename=paste(args[1],"plotTMscore.svg", sep=""),width=5,height=5,bg='transparent')

plot(x=mat$Level, y= mat$TM.scoreA, xlab="Level", ylab="TM-score", type="l", col="red", lwd=2, ylim=c(0,1))
lines(x=mat$Level, y= mat$TM.scoreB, xlab="Level", ylab="TM-score", col="blue", lwd=2)
legend("bottomleft", legend=c("Normalisé selon A", "Normalisé selon B"),col=c("red", "blue"), lwd=2, cex=0.8)