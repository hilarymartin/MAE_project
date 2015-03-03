cols=c("red","darkgreen","purple","blue","pink")
argv <- commandArgs(TRUE)

threshold=0.5
thr <- 1E6 #recombinations within this distance of one another are filtered

male.lengths <- c(1.9512,1.8955,1.6071,1.4654,1.512,1.3762,1.2835,1.0794,1.1725,1.3389,1.0936,1.3554,1.0131,0.9462,1.0257,1.081,1.0856,0.9862,0.9264,0.7472,0.4731,0.4896)
female.lengths <- c(3.4541,3.2541,2.7564,2.5906,2.6019,2.4159,2.3033,2.0994,1.982,2.1813,1.9553,2.0664,1.5588,1.4236,1.5496,1.4962,1.6153,1.4257,1.2682,1.2197,0.764,0.8276)
malesize <- sum(male.lengths)
femalesize <- sum(female.lengths)
male.thr=thr * femalesize/malesize


cohort.list=read.delim("cohort_list.txt",header=F,stringsAsFactors=F)

snp.filter=c()
for(j in 1:nrow(cohort.list)){
  mydir=cohort.list[j,1]
  cohort=cohort.list[j,2]
  prefix=cohort.list[j,3]
print(cohort)

rall <- list()
#recomb.file=argv[1]
recomb.file=paste0(mydir,"/duoHMM_results/",prefix,".more_stringent.generr_removed.recombinations.txt")
#map.file=argv[4]
map.file=paste0(mydir,"/input_for_SHAPEIT2/",prefix,".more_stringent.generr_removed.map")

total.SNPs=c()
total.dist=c()

maps=list()
                                        #read in crossovers and map files
for(i in 1:22) {
        rall[[i]] <- read.table(gsub("1.",paste(i,".",sep=""),recomb.file,fixed=T),header=TRUE,as.is=TRUE,colClasses=c("character","character","integer","integer","numeric"))
        rall[[i]]$chr <- i
        map<- read.table(gsub("1.",paste(i,".",sep=""),map.file,fixed=T),header=FALSE,as.is=TRUE,colClasses=c("integer","character","numeric","integer"))
        colnames(map)=c("chr","SNP","cM","pos")
        map$SNP.n = 1:nrow(map)
        maps[[i]]=map
#work out the number of the start and end SNPs
        rall[[i]]$START.SNP.n = map[match(rall[[i]]$START,map$pos),"SNP.n"]
        rall[[i]]$END.SNP.n = map[match(rall[[i]]$END,map$pos),"SNP.n"]
        total.SNPs=c(total.SNPs,nrow(map))
        total.dist=c(total.dist,map[nrow(map),4]-map[1,4])
    }

mean.snp.space=sum(as.numeric(total.SNPs))/sum(as.numeric(total.dist))
snp.filter=c(snp.filter,mean.snp.space)

}
names(snp.filter) = cohort.list[,3]

n.snps.threshold = thr * snp.filter

write.table(n.snps.threshold,"qqplots_duoHMM_v3/threshold_used_for_filtering_double_crossovers_within_X_SNPs.by_cohort.txt",quote=F,sep="\t")
