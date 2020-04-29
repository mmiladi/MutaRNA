

rm(list=ls())

# supported output extensions
outExt = data.frame(
  mode = c("png","pdf","svg"),
  regex = c("^.+\\.[pP][nN][gG]$","^.+\\.[pP][dD][fF]$","^.+\\.[sS][vV][gG]$")
)

#### dependencies ###

suppressPackageStartupMessages({
  library(stringr)
  library(ggplot2)
  library(ggforce)
  library(ggthemes)
})

#### input ###

args = commandArgs(trailingOnly=TRUE)
#args = c("tmp.ps","tmp.png")
# check arguments
if (length(args)!=2) {
  stop(paste("\nCall Rscript --vanilla dot2ps2arcplot.R <viennaDot.ps> <output.[",paste(outExt$mode,collapse="|"),"]>\n",sep=""), call.=FALSE)
} 

# distribute input arguments
psFile = args[1]
outFile = args[2]

#### parsing ###

# parse output mode via file ending
outMode = NA
for (i in 1:nrow(outExt) ) {
  if (1==sum(grep(outExt$regex[i],outFile))) {
    outMode = outExt$mode[i]
    break;
  } 
}
if (is.na(outMode)) {
  stop(paste("output mode (defined by outFile extension) not supported. Has to be one of [",outExt$mode,"]",sep=""))
}

# parse dotplot PS file
if (!file.exists(psFile)) {
  stop(paste("Cannot open input file",psFile,"for reading!"))
}
FileInput = readLines(psFile) 

# parse sequence
seq = gsub("\\\\","",paste(FileInput[grep("^[acgutnACGUTN\\\\\\s]+$",FileInput)],collapse=""))
seqLength = str_length(seq)

# parse base pair probabilities
bpProb = data.frame(lapply(data.frame(matrix((unlist(lapply(FileInput[grep("^(\\d+.+ubox)$",FileInput)],strsplit, " "))),ncol=4,byrow = T)[,1:3], stringsAsFactors=F), as.numeric))
colnames(bpProb)=c("i","j","prob")
# compute additional values
bpProb[,"r"] = (bpProb$j - bpProb$i)/2
bpProb[,"c"] = bpProb$i + bpProb$r
# sort by probability
bpProb <- bpProb[sort.list(bpProb$prob,decreasing = F),]

# set output theme settings
old <- theme_set(theme_clean())
theme_update(axis.title.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),
             axis.line.y=element_blank(),
             legend.background = element_rect(fill="white",color="white"))

if (outMode=="pdf") {
  pdf(outFile, width= seqLength/100*6, height = max(bpProb$r)/100*5)
} else if (outMode=="png") {
  png(outFile, units="in", res=300, width= seqLength/100*6, height = max(bpProb$r)/100*5)
} else if (outMode=="svg") {
  svg(outFile, bg="white", width= seqLength/100*6, height = max(bpProb$r)/100*5)
}

# plot 
(p <- ggplot(bpProb) +
  geom_arc(aes(x0 = c, y0 = 0, r=r,
               start = -pi/2, end = pi/2,
               col = prob), size=0.6) +
  scale_colour_gradient(low = "yellow", high = "red") +
    ylab(NULL) +
  xlim(1,seqLength) +
  labs(col="Pr(bp)",x="sequence position (nt)",y="") +
  annotate("text", x=1:seqLength,y= -1, label=unlist(strsplit(seq,"")),size=1.5)
)
seqChar = unlist(strsplit(seq,""))
for( i in 1)

# close output file
dev.off()

