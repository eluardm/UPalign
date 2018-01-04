#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)

mat = NULL

for(i in seq(1,strtoi(args[1]))){

	name = paste(args[3],"plotAligned",i,".dat", sep="")
	matTMP = read.table(name, head=T)
	mat = rbind(mat,matTMP)
}

svg(filename=paste(args[3],"plotAlign.svg"),width=15,height=5,bg='transparent')
plot(x=mat$Residus, y=mat$Level, pch=15, xlim=c(0,strtoi(args[2])), col=mat$NumPU, xlab="Residus", ylab="Level", xaxt="n")
xtick<-seq(0, strtoi(args[2]), by=10)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick, par("usr")[3], labels = xtick, pos = 1, xpd = TRUE)

invisible(dev.off())
