# Description: Plotting CNVs across the chromosome
#
# Author:      Mohamood Adhil
# Date:        1st Feb 2018
# help:        Rscript <R-script> -h (or) --help
# TODO:        
# #################################################################################
############## 
#library
library(argparse)
library(ggplot2)
options("scipen"=100, "digits"=4)
##############

parser <- ArgumentParser()
parser <- ArgumentParser(description='This function is for plotting CNVs across the chromosome')
parser$add_argument("-int", dest="intf", help=".alldat file from output of cnvlocator script. Multiple files can be given with comma seperator (Mandatory)", required = TRUE)
parser$add_argument("-lfc", dest="lfc", type="integer",help="Log fold change threshold value Default = 0.6", default=0.6)
parser$add_argument("-mn", dest="mn", type="integer", help="Minimum percentage of reads atleast in one sample test or reference Default = 0.05", default=0.05)
parser$add_argument("-ulc", dest="ulc", type="integer", help="Upper limit cutoff for the graph Default = 3", default=3)
parser$add_argument("-llc", dest="llc", type="integer", help="Lower limit cutoff for the graph, value should be below 0, Default= -3", default=-3)
parser$add_argument("-cloc", dest="cloc", default="NA", help="Plot for particular chromosome location (format eg chr;start;end) default = plot All chromosome")
parser$add_argument("-tgs", dest="tgs", required=TRUE, help="Prefix tag for the plot. Multiple tags can be given with comma seperator, number of tags should be same as the number of file names (Mandatory)")
parser$add_argument("-dir", dest="dir", required=TRUE, help="Output directory path (Mandatory)")
args <- parser$parse_args()

#############################

samples = strsplit(args$intf,",")[[1]]
lfc = as.numeric(args$lfc)
mn = as.numeric(args$mn)
if (args$cloc != "NA"){
  clocval = strsplit(args$cloc,";")[[1]]
  cloc = "av"
} else {
  cloc = "NA"
}
tgs = strsplit(args$tgs,",")[[1]]
outdir = args$dir
ulc = args$ulc
llc = args$llc
###############################

for (i in 1:length(samples)){
data = read.delim(file = samples[i], sep="\t", h=T)
data = data[-which(data$tnorm <= mn & data$rnorm <= mn),]
data = data[,c(1,2,10)]
data$chr = as.vector(data$chr)
data$order = c(1:dim(data)[1])
if (cloc != "NA"){
data = data[which(data$chr == clocval[1] & data$start >= as.numeric(clocval[2]) & data$start <= as.numeric(clocval[3])),]
}
data$Type = "NO"
data[which(data$lfc >= 0.6),"Type"] = "Amplification"
data[which(data$lfc <= -0.6),"Type"] = "Deletion"
data$Type = factor(data$Type)
data$size = 0.50
data$size = as.numeric(data$size)
te = aggregate(order ~ chr, data, max)
te = te[order(te$order),]
rownames(te) = c(1:dim(te)[1])
level <- as.character(unique(data$chr))
colour = 7
data$colour <- as.character(match(data$chr, level)%%colour)
pcolors = c("orangered1","mediumpurple","brown2","mediumorchid1","palevioletred1","darkseagreen4","darkorange")
for (j in 1:colour){
  data[which(data$colour == j-1),"colour"] = pcolors[j]
}
data$opacity = 0.80
data[which(data$Type != "NO"),"size"] = 1
if (cloc != "NA"){
  xtag = paste("chromosome",";",clocval[1],"_",clocval[2],"-",clocval[3],sep="")
  p <- ggplot(data,aes(x = start, y = lfc, colour = colour)) + geom_point(alpha = data$opacity,size=data$size,colour=data$colour) + ylim(llc,ulc) +
    theme(panel.grid.minor = element_blank()) + theme(legend.position="none") + xlab(xtag) + ylab("Log2FoldChange") + geom_hline(yintercept = c(-lfc,lfc)) + theme_bw()
} else {
p <- ggplot(data,aes(x = order, y = lfc, colour = colour)) + geom_point(alpha = data$opacity,size=data$size,colour=data$colour) + ylim(llc,ulc) + scale_x_continuous(breaks=c(1,te$order[-17]+1),labels=te$chr) +
  theme(panel.grid.minor = element_blank()) + theme(legend.position="none") + xlab ("Chromosome") + ylab("Log2FoldChange") + geom_hline(yintercept = c(-lfc,lfc)) + theme_bw()
}
ggsave(file = paste(outdir,"/",tgs[i],"_cnvplot.jpeg", sep=""), p, height = 8, width = 16, dpi = 600)
}
