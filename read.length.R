#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

data=read.csv(args[1])

print (args[1])




pdf(args[2],width=6,height=4)
ggplot(data, aes(x=l, color=flag)) +geom_histogram(fill="white")+facet_grid(flag ~ .)


#+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="none")+labs(x = "Read length")+theme(axis.text=element_text(size=15,face="bold"),axis.title=element_text(size=15))

#ggplot(data, aes(x = data$l, y = data$flag, fill=data$f)) + geom_density_ridges(scale = 2) +

#theme(axis.text.x = element_text(angle=90)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(legend.position="none")+labs(x = "Read length")+theme(axis.text=element_text(size=15,face="bold"),axis.title=element_text(size=15))
dev.off()
