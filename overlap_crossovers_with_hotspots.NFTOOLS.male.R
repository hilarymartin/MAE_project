#source("http://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)
#QTR and FC - hg18
#NTR hg19
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

out.prefix = args[2]
recomb.file=args[1]
if(out.prefix %in% c("FC","QTR610","QTR370")){
  hotspot.file="hotspots_from_Anjali.hg18.plus_1.bed"
} else {
  hotspot.file  = "hotspots_from_Anjali.lifted_to_hg19.bed"
}

hg18.hotspots=read.delim(hotspot.file,header=F,stringsAsFactors=F)

hotspots=GRanges(seqnames=hg18.hotspots[,1],ranges=IRanges(start=hg18.hotspots[,2],end=hg18.hotspots[,3]))



duoHMM.all=NULL

for(i in 1:22){
d=read.delim(gsub(".1.",paste0(".",i,"."),recomb.file,fixed=T),header=T,stringsAsFactors=F)
d$chr=i
duoHMM.all=rbind(duoHMM.all,d)
}

###keep only females
duoHMM.all=duoHMM.all[duoHMM.all$sex=="male",]

duoHMM.all$chr=paste0("chr",duoHMM.all$chr)
duoHMM.all$chr=gsub("k","",duoHMM.all$chr,perl=T)
duoHMM.all$chr=gsub(".","",duoHMM.all$chr,fixed=T)

print(unique(duoHMM.all$chr))
#sizes=c(max(duoHMM.all$right-duoHMM.all$left) + 1,100000,50000,30000,20000,10000,5000)
sizes=c(max(duoHMM.all$right-duoHMM.all$left) + 1,60000,50000,40000,30000,20000)
#sizes=c(30000)
if(sum((duoHMM.all$right-duoHMM.all$left)<0) >0){
cat("ERROR: some of the crossovers in ",recomb.file," have negative length\n")
}
overlap.stats=NULL
for(size in sizes){
print(size)
duoHMM.resolved=duoHMM.all[duoHMM.all$right-duoHMM.all$left < size,]
crossovers=GRanges(seqnames=duoHMM.resolved$chr,ranges=IRanges(start=duoHMM.resolved$left,end=duoHMM.resolved$right))
n.true.overlaps=countOverlaps(crossovers,hotspots)
total.true.overlapping=sum(n.true.overlaps >0)
prop.true.overlapping=total.true.overlapping/nrow(duoHMM.resolved)

### calculate proportion by chance
#permutations
true.crossovers=duoHMM.resolved
prop.overlapping.samples=c()

for(i in 1:100){
#for(i in 1:10){
permuted.crossovers=true.crossovers
permuted.crossovers$shift=rnorm(n=nrow(permuted.crossovers),mean=0,sd=200000)
permuted=GRanges(seqnames=permuted.crossovers$chr,ranges=IRanges(start=permuted.crossovers$left+permuted.crossovers$shift,end=permuted.crossovers$right + permuted.crossovers$shift))
n.overlaps=countOverlaps(permuted,hotspots)
total.overlapping=sum(n.overlaps >0)
prop.overlapping=total.overlapping/nrow(permuted.crossovers)
prop.overlapping.samples=c(prop.overlapping.samples,prop.overlapping)
}
mean.prop.overlapping=mean(prop.overlapping.samples)

#Bleazard's method
actual.prop.in.hotspots=(prop.true.overlapping-mean.prop.overlapping)/(1-mean.prop.overlapping)

#Coop's method - evaluate likelihood under a grid of alphas, and pick most likely
alpha.grid=seq(0,1,0.001)
log.likelihood = sapply(alpha.grid,function(alpha){total.true.overlapping *log(alpha+(1-alpha)*mean.prop.overlapping) + (nrow(duoHMM.resolved)-total.true.overlapping)*log(1-alpha-(1-alpha)*mean.prop.overlapping)})
MLE.alpha = alpha.grid[which.max(log.likelihood)]
conf.interval=which(abs(log.likelihood-max(log.likelihood))<2)
conf.int.lower=alpha.grid[conf.interval[1]]
conf.int.upper=alpha.grid[conf.interval[length(conf.interval)]]

overlap.stats=rbind(overlap.stats,c(nrow(duoHMM.resolved),nrow(duoHMM.resolved)/nrow(duoHMM.all),total.true.overlapping,prop.true.overlapping,mean.prop.overlapping,actual.prop.in.hotspots,MLE.alpha,conf.int.lower,conf.int.upper))
}

overlap.stats=data.frame(overlap.stats)
colnames(overlap.stats)=c("n.crossovers","prop.all.crossovers","n.overlapping.hotspots","prop.overlapping.hotspots","prop.overlapping.hotspots.by.chance","alpha.30kb.Bleazard","alpha.30kb.Coop","alpha.30kb.Coop.95pc.lower","alpha.30kb.Coop.95pc.upper")
#overlap.stats$resolution=c("all","100kb","50kb","30kb","20kb","10kb","5kb")
overlap.stats$resolution=sizes


write.table(overlap.stats,paste0("hotspot_overlap/",out.prefix,".NFTOOLS.male.overlap_with_hotspots.txt"),quote=F,sep="\t")
