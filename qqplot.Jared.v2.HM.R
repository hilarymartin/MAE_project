#cols <- paste0(c("#E41A1C","#377EB8","#4DAF4A","#984EA3"),"CC")
cols=c("red","darkgreen","purple","blue","pink")
argv <- commandArgs(TRUE)
threshold=0.5
thr <- 1E6 #recombinations within this distance of one another are filtered

male.lengths <- c(1.9512,1.8955,1.6071,1.4654,1.512,1.3762,1.2835,1.0794,1.1725,1.3389,1.0936,1.3554,1.0131,0.9462,1.0257,1.081,1.0856,0.9862,0.9264,0.7472,0.4731,0.4896)
female.lengths <- c(3.4541,3.2541,2.7564,2.5906,2.6019,2.4159,2.3033,2.0994,1.982,2.1813,1.9553,2.0664,1.5588,1.4236,1.5496,1.4962,1.6153,1.4257,1.2682,1.2197,0.764,0.8276)
malesize <- sum(male.lengths)
femalesize <- sum(female.lengths)
male.thr=thr * femalesize/malesize
#argv=c("NTR/duoHMM_results/NTR_v2_minus_MZs.1.post_MERLIN.max_0.01_missing.generr_removed.recombinations.txt",
#    "NTR/input_for_SHAPEIT2/NTR_v2_minus_MZs.1.post_MERLIN.max_0.01_missing.generr_removed.fam",
#    "qqplots_duoHMM_v2/NTR_v2.post_MERLIN.max_0.01_missing.generr_removed","NTR/input_for_SHAPEIT2/NTR_v2_minus_MZs.1.post_MERLIN.max_0.01_missing.generr_removed.map","NTR_v2.max_0.01_missing",
#    "NTR/duoHMM_results/NTR_v2.post_MERLIN.max_0.05_missing.generr_removed_sample_information.txt") 
#argv=c("NTR/duoHMM_results/NTR_v2_minus_MZs.1.post_MERLIN.max_0.05_missing.generr_removed.recombinations.txt","NTR/input_for_SHAPEIT2/NTR_v2_minus_MZs.1.post_MERLIN.max_0.05_missing.generr_removed.fam",
#    "qqplots_duoHMM_v2/NTR_v2.post_MERLIN.generr_removed","NTR/input_for_SHAPEIT2/NTR_v2_minus_MZs.1.post_MERLIN.max_0.05_missing.generr_removed.map","NTR_v2",
#    "NTR/duoHMM_results/NTR_v2.post_MERLIN.max_0.05_missing.generr_removed_sample_information.txt")
#argv=c("FC/duoHMM_results/FC.1.post_MERLIN.generr_removed.recombinations.txt","FC/input_for_SHAPEIT2/FC.1.post_MERLIN.generr_removed.max_0.05_missing.fam",
#    "/home/hilary/maternal_age_recombination/qqplots_duoHMM_v2/FC.post_MERLIN.generr_removed","FC/input_for_SHAPEIT2/FC.1.post_MERLIN.generr_removed.map","FC",
#    "FC/duoHMM_results/FC.post_MERLIN.generr_removed_sample_information.txt")
#argv=c("GPC/duoHMM_results/GPC.1.post_MERLIN.generr_removed.recombinations.txt","GPC/input_for_SHAPEIT2/GPC.chr1.no_mendel_errors.post_MERLIN.generr_removed.fam",
#    "/home/hilary/maternal_age_recombination/qqplots_duoHMM_v2/GPC.post_MERLIN.generr_removed","GPC/input_for_SHAPEIT2/GPC.chr1.no_mendel_errors.post_MERLIN.generr_removed.map","GPC",
#    "GPC/duoHMM_results/GPC.post_MERLIN.generr_removed_sample_information.txt")
#argv=c("QTR/duoHMM_results/QTR_minus_MZs.610K.1.post_MERLIN.max_0.05_missing.generr_removed.recombinations.txt","QTR/input_for_SHAPEIT2/QTR_minus_MZs.610K.1.post_MERLIN.max_0.05_missing.generr_removed.fam",
##    "qqplots_duoHMM_v2/QTR610.post_MERLIN.generr_removed","QTR/input_for_SHAPEIT2/QTR_minus_MZs.610K.1.post_MERLIN.max_0.05_missing.generr_removed.map","QTR610",
 #   "QTR/duoHMM_results/QTR.610K.post_MERLIN.max_0.05_missing.generr_removed_sample_information.txt")
#argv=c("FVG/duoHMM_results/FVG.1.generr_removed.recombinations.txt","FVG/clean-fvg-flip_b37-related.fam","/home/hilary/maternal_age_recombination/qqplots_duoHMM_v2/FVG.post_MERLIN.generr_removed",
#    "FVG/input_for_SHAPEIT2/FVG.chr1.generr_removed.map","FVG","FVG/duoHMM_results/FVG.generr_removed_sample_information.txt")
#argv=c("ORCADES/duoHMM_results/ORCADES_raw_merged.1.max_0.05_missing.generr_removed.recombinations.txt","ORCADES/input_for_SHAPEIT2/ORCADES_raw_merged.1.no_mendel_errors.max_0.05_missing.generr_removed.fam",
#    "/home/hilary/maternal_age_recombination/qqplots_duoHMM_v2/ORCADES.post_MERLIN.generr_removed","ORCADES/input_for_SHAPEIT2/ORCADES_raw_merged.1.no_mendel_errors.max_0.05_missing.generr_removed.map","ORCADES",
#    "ORCADES/duoHMM_results/ORCADES_raw_merged.max_0.05_missing.generr_removed_sample_information.txt")
#argv=c("VIS_KORCULA/duoHMM_results/VIS.1.post_MERLIN.generr_removed.recombinations.txt","VIS_KORCULA/input_for_SHAPEIT2/VIS.chr1.no_mendel_errors.post_MERLIN.generr_removed.fam",
#    "qqplots_duoHMM_v2/VIS.post_MERLIN.generr_removed","VIS_KORCULA/input_for_SHAPEIT2/VIS.chr1.no_mendel_errors.post_MERLIN.generr_removed.map","VIS",
#    "VIS_KORCULA/duoHMM_results/VIS.post_MERLIN.generr_removed_sample_information.txt")
#argv=c("VIS_KORCULA/duoHMM_results/KORCULA.1.post_MERLIN.generr_removed.recombinations.txt","VIS_KORCULA/input_for_SHAPEIT2/KORCULA.chr1.no_mendel_errors.post_MERLIN.generr_removed.fam",
#    "/home/hilary/maternal_age_recombination/qqplots_duoHMM_v2/KORCULA.post_MERLIN.generr_removed","VIS_KORCULA/input_for_SHAPEIT2/KORCULA.chr1.no_mendel_errors.post_MERLIN.generr_removed.map","KORCULA",
#    "VIS_KORCULA/duoHMM_results/KORCULA.post_MERLIN.generr_removed_sample_information.txt")

cohort=argv[5]
#load(paste0(argv[3],".multiple_filters.RData"))
qqr <- function(k1,rate,col1,makeplot=TRUE,title,limits,xlowerlimit,...) {
  k1 <- sort(k1)
  n <- length(k1)
  e <- qpois(1:n/(n+1),rate)
  l <- c(min(k1),min(max(k1),rate+6*sqrt(rate)))
  etmp <- e[k1>=l[2]]
  ret <- cbind(jitter(etmp),l[2])
  
  if(makeplot) {
     # plot(e,k1,type='n',xlab="",ylab="",ylim=l,main=title,...)
       plot(e,k1,type='n',xlab="",ylab="",ylim=limits,xlim=c(xlowerlimit,max(limits)),main=title,...)
      abline(0,1,col=1,lty=2);grid()
      points(jitter(e[k1<l[2]]),jitter((k1[k1<l[2]])),col=col1,pch=16)
      points(ret,col=col1,pch=17)
  }

  return(  cbind(e,k1))
}

tabulateRecombination <- function(r,keep) {
    rec <- data.frame(r[!duplicated(r[,c('CHILD','PARENT',keep)]),c("CHILD","PARENT",keep)],stringsAsFactors=FALSE)
    colnames(rec)[3]="keep"
    rec$duo <- paste(rec$CHILD,rec$PARENT,sep="-")
    rec$sex <- as.factor(c("Male","Female")[as.integer(fam[match(rec$PARENT,fam$V2),5])] )
#    rec$informative <- fam[match(rec$PARENT,fam$V2),]$grandparent | fam[match(rec$PARENT,fam$V2),]$nkid>2
    rec$nkid <-  fam[match(rec$PARENT,fam$V2),]$nkid
    rec$nrec <- table(r$duo)[rec$duo]
    rec=cbind(rec,r[match(rec$duo,r$duo),c("informative","informative.2gen","informative.2gen.2parents","informative.2gen.1parent","noninf.2kids.2gen.2parents","noninf.2kids.2gen.1parent",
        "noninf.1kid.2gen.2parents","noninf.1kid.2gen.1parent","informative.3gen","informative.3gen.2parents", "informative.3gen.1parent")])
    rec <- subset(rec,keep)
    print(tapply(rec$nrec,rec$sex,mean))
    rec
}

#if(FALSE){

fam.file=argv[2]
#fam <- read.table(fam.file,colClasses="character",header=F)
sample.info.file=argv[6]

###integrate into data.frame r and print summaries for different family types
#sample.info=read.delim(sample.info.file,header=T,colClasses=c(rep("character",5),rep("integer",17-5),"logical"))
fam=read.delim(sample.info.file,header=T,colClasses=c(rep("character",5),rep("integer",4),rep("character",4),rep("integer",4),"logical"))
colnames(fam)[1:5]=paste0("V",1:5)
 #change names because it will interfer with code later on
fam[,2]=gsub("-","_",fam[,2])
fam[,3]=gsub("-","_",fam[,3])
fam[,4]=gsub("-","_",fam[,4])

#n kids with same parents
fam$nkid <- sapply(fam$V2,function(x) sum(x==fam$V3)+sum(x==fam$V4))
#mother or father is genotyped
fam$grandparent <- fam$V3 %in% fam$V2 | fam$V4 %in% fam$V2
#fam=cbind(fam,sample.info[match(fam[,2],sample.info[,2]),-c(1:5)])
rall <- list()
recomb.file=argv[1]
map.file=argv[4]
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

#load(paste0(argv[3],".RData"))
                                        #threshold for number of SNPs between adjacent crossovers - this corrects for SNP density, so we are not biased towards removing double crossovers in regions with denser SNPs
n.snps.threshold= round(thr * mean.snp.space)
n.snps.threshold.male= round(male.thr * mean.snp.space)



#keep those with prob > threshold
r <- subset(do.call("rbind",rall),PROB_RECOMBINATION>threshold)
r$CHILD=gsub("-","_",r$CHILD)
r$PARENT=gsub("-","_",r$PARENT)

remove.outliers=function(r,cohort){
    #remove families in FC cohort due to MZ twins/duplicates (76), and because of low maternal counts, reason unknown (52)
    if(cohort=="FC"){
        r=r[!r[,2] %in% fam[fam[,1] %in% c("52", "76"),2],]
    }
    #remove families FVG to due to pedigree erros
    if(cohort=="FVG"){
      r=r[!r[,2] %in% fam[fam[,1] %in% c("69", "148"),2],]
    }
    #remove families in VB: 143 and 64 each contain a cryptic pair of MZ twins or duplicate samples
    if(cohort=="VB"){
        r=r[!r[,2] %in% fam[fam[,1] %in% c("64", "143"),2],]
    }
    #remove families from NTR: 10689 seems to have three times the same sample; 11880 is strange (maybe too much missingness in mother)
    if(cohort=="NTR"){
       r=r[!r[,2] %in% fam[fam[,1] %in% c("10689", "11880"),2],]
    }
    if(cohort=="GPC"){
       #remove dodgy family 188 (mother APP5117398, father APP5117401 or    APP5117402_dummyfather from GPC; think relationships are mis-specified
        r=r[r$PARENT!="APP5117398" & r$PARENT !="APP5117401",]
    }
    return(r)
}

r=remove.outliers(r,cohort)


r$duo <- paste0(r$CHILD,"-",r$PARENT)
r$sex <- as.factor(c("Male","Female")[as.integer(fam[match(r$PARENT,fam$V2),5])] )
#if the parent's mother and/or father is genotyped, or if the parent has > 2 kids
r$informative <- fam[match(r$PARENT,fam$V2),]$grandparent | fam[match(r$PARENT,fam$V2),]$nkid>2

#at least 3 kids with same mother; only two generations
r$informative.2gen <- !fam[match(r$PARENT,fam$V2),]$grandparent & fam[match(r$PARENT,fam$V2),]$nkid>2

#at least 3 kids with same mother and father; only two generations
r$informative.2gen.2parents <- !fam[match(r$PARENT,fam$V2),]$grandparent & fam[match(r$PARENT,fam$V2),]$nkid>2 & !r$CHILD %in% fam[fam[,3]=="0"|fam[,4]=="0",2] &
    r$CHILD %in% fam[fam$total_same_parents==fam$total_same_mother,2]

#at least 3 kids with same mother but father missing for all; only two generations
r$informative.2gen.1parent <- !fam[match(r$PARENT,fam$V2),]$grandparent & fam[match(r$PARENT,fam$V2),]$nkid>2 & r$CHILD %in% fam[fam[,3]=="0"|fam[,4]=="0",2] &
    (r$CHILD %in% fam[fam$total_same_parents==fam$total_same_mother,2]|r$CHILD %in% fam[fam$total_same_parents==fam$total_same_father,2])

#2 kids with same mother and father; only two generations
r$noninf.2kids.2gen.2parents <- !fam[match(r$PARENT,fam$V2),]$grandparent & fam[match(r$PARENT,fam$V2),]$nkid==2 & !r$CHILD %in% fam[fam[,3]=="0"|fam[,4]=="0",2] &
    r$CHILD %in% fam[fam$total_same_parents==fam$total_same_mother,2]
#2 kids with same mother but father missing for all; only two generations
r$noninf.2kids.2gen.1parent <- !fam[match(r$PARENT,fam$V2),]$grandparent & fam[match(r$PARENT,fam$V2),]$nkid==2 & r$CHILD %in% fam[fam[,3]=="0"|fam[,4]=="0",2] &
    (r$CHILD %in% fam[fam$total_same_parents==fam$total_same_mother,2]|r$CHILD %in% fam[fam$total_same_parents==fam$total_same_father,2])

#1 kid with same mother and father; only two generations
r$noninf.1kid.2gen.2parents <- !fam[match(r$PARENT,fam$V2),]$grandparent & fam[match(r$PARENT,fam$V2),]$nkid==1 & !r$CHILD %in% fam[fam[,3]=="0"|fam[,4]=="0",2] &
    r$CHILD %in% fam[fam$total_same_parents==fam$total_same_mother,2]

#1 kid with same mother but father missing for all; only two generations
r$noninf.1kid.2gen.1parent <- !fam[match(r$PARENT,fam$V2),]$grandparent & fam[match(r$PARENT,fam$V2),]$nkid==1 & r$CHILD %in% fam[fam[,3]=="0"|fam[,4]=="0",2] &
    (r$CHILD %in% fam[fam$total_same_parents==fam$total_same_mother,2]|r$CHILD %in% fam[fam$total_same_parents==fam$total_same_father,2])

#1 or both maternal grandparents
r$informative.3gen = fam[match(r$PARENT,fam$V2),]$grandparent 
#1 or both maternal grandparents, both parents same and genotyped
r$informative.3gen.2parents = fam[match(r$PARENT,fam$V2),]$grandparent & !r$CHILD %in% fam[fam[,3]=="0"|fam[,4]=="0",2] & r$CHILD %in% fam[fam$total_same_parents==fam$total_same_mother,2]
#1 or both maternal grandparents, father missing
r$informative.3gen.1parent = fam[match(r$PARENT,fam$V2),]$grandparent & r$CHILD %in% fam[fam[,3]=="0"|fam[,4]=="0",2] &
    (r$CHILD %in% fam[fam$total_same_parents==fam$total_same_mother,2]|r$CHILD %in% fam[fam$total_same_parents==fam$total_same_father,2])

r$dummy=rep(TRUE,nrow(r))
#calculate distance to nearest crossover
r$diff <- c(NA,tail(r$END,-1) - head(r$START,-1))
#set distance to NA if the adjacent rows are not the same chromosome, child or parent
r$diff[c(FALSE,tail(r$chr,-1) != head(r$chr,-1)) | c(FALSE,tail(r$CHILD,-1) != head(r$CHILD,-1)) | c(FALSE,tail(r$PARENT,-1) != head(r$PARENT,-1))] <- NA
#remove double crossovers within the specified distance
r.flt <- r[-c(which(r$diff<thr),which(r$diff<thr)-1),]

r.double.xovers.within.N.Mb <- r[c(which(r$diff<thr),which(r$diff<thr)-1),]

#calculate distance to nearest crossover, in number of SNPs
r$SNP.diff <- c(NA,tail(r$END.SNP.n,-1) - head(r$START.SNP.n,-1))
#set distance to NA if the adjacent rows are not the same chromosome, child or parent
r$SNP.diff[c(FALSE,tail(r$chr,-1) != head(r$chr,-1)) | c(FALSE,tail(r$CHILD,-1) != head(r$CHILD,-1)) | c(FALSE,tail(r$PARENT,-1) != head(r$PARENT,-1))] <- NA
#remove double crossovers within the specified distance
r.flt1 <- r[-c(which(r$SNP.diff<n.snps.threshold),which(r$SNP.diff<n.snps.threshold)-1),]
#r.male=r[r$sex=="Male",]
#r.female=r[r$sex=="Female",]
#r.male.flt1=r.male[-c(which(r.male$SNP.diff<n.snps.threshold.male),which(r.male$SNP.diff<n.snps.threshold.male)-1),]
#r.female.flt1=r.female[-c(which(r.female$SNP.diff<n.snps.threshold),which(r.female$SNP.diff<n.snps.threshold)-1),]
#r.flt1=rbind(r.male.flt1,r.female.flt1)
r.double.xovers.within.X.SNPs <- r[c(which(r$SNP.diff<n.snps.threshold),which(r$SNP.diff<n.snps.threshold)-1),]

#pull out any crossovers that are duplicated within or across families
dup.pos <- r[which(duplicated(r[,c("chr","START","END")])),c("chr","START","END")]
dup.pos$xover=paste(dup.pos$chr,dup.pos$START,dup.pos$END,sep="_")
#remove these
r$xover=paste(r$chr,r$START,r$END,sep="_")
r.flt2 <- r[!r$xover %in% dup.pos$xover,]
#r.flt2 <- r[-which(r$START%in%dup.pos[,2] & r$END%in%dup.pos[,3] & r$chr%in%dup.pos[,1]),]
r.dup.pos <- r[r$xover %in% dup.pos$xover,]

#pull out any crossovers that are duplicated within families
r$father=fam[match(r$CHILD,fam[,2]),3]
r$mother=fam[match(r$CHILD,fam[,2]),4]
dup.pos.within.fam <- r[which(duplicated(r[,c("chr","START","END","father","mother")])),c("chr","START","END","father","mother")]
dup.pos.within.fam$xover=paste(dup.pos.within.fam$chr,dup.pos.within.fam$START,dup.pos.within.fam$END,dup.pos.within.fam$father,dup.pos.within.fam$mother,sep="_")
r$xover2=paste(r$chr,r$START,r$END,r$father,r$mother,sep="_")
r.flt3 <- r[!r$xover2 %in% dup.pos.within.fam$xover,]
r.dup.pos.within.fam  <- r[r$xover2 %in% dup.pos.within.fam$xover,]

#after removal of double crossovers within X SNPs, pull out any that are duplicated within or across families
dup.pos2 <- r.flt1[which(duplicated(r.flt1[,c("chr","START","END")])),c("chr","START","END")]
dup.pos2$xover=paste(dup.pos2$chr,dup.pos2$START,dup.pos2$END,sep="_")
#remove these
r.flt1$xover=paste(r.flt1$chr,r.flt1$START,r.flt1$END,sep="_")
r.flt1.2 <- r.flt1[!r.flt1$xover %in% dup.pos2$xover,]
r.dup.pos.2 <- r.flt1[r.flt1$xover %in% dup.pos2$xover,]

plot.outlier.snps = function(r.spurious,maps,filename,filter.threshold){
    n.spurious.by.chr=sapply(1:22,function(chr){
        chr.map=maps[[chr]]
        chr.r.spurious=r.spurious[r.spurious$chr==chr,]
        x=chr.map$pos
        n.spurious=sapply(chr.map$pos,function(pos){return(sum( chr.r.spurious$START < pos & chr.r.spurious$END > pos))})
    })
    all.spurious = unlist(n.spurious.by.chr)
    cutoff=quantile(all.spurious,filter.threshold)
    pdf(filename,height=12,width=18)
    par(mfrow=c(4,6))
    for(i in 1:22){
        x=maps[[i]]$pos
        n.spurious=n.spurious.by.chr[[i]]
        plot(x/1E6,n.spurious,main=paste0("chr",i),type="l",col="blue",xlab="Position (Mb)",ylab="Number of events")
            abline(h=cutoff,type=2)
    }
    dev.off()
    return(n.spurious.by.chr)
}
freq.duplicate.xovers=plot.outlier.snps(r.dup.pos,maps,paste0(argv[3],".frequency_of_duplicate_crossovers.pdf"),0.999)
freq.duplicate.xovers.after.double.xover.removal=plot.outlier.snps(r.dup.pos.2,maps,paste0(argv[3],".frequency_of_duplicate_crossovers_after_removing_double_crossovers_within_X_SNPs.pdf"),0.999)
freq.duplicate.xovers.within.fam=plot.outlier.snps(r.dup.pos.within.fam,maps,paste0(argv[3],".frequency_of_duplicate_crossovers_within_fam.pdf"),0.999)
freq.double.xovers.within.X.SNPs=plot.outlier.snps(r.double.xovers.within.X.SNPs,maps,paste0(argv[3],".frequency_of_double_crossovers_within_X_SNPs.pdf"),0.999)
freq.double.xovers.within.N.Mb = plot.outlier.snps(r.double.xovers.within.N.Mb,maps,paste0(argv[3],".frequency_of_double_crossovers_within_",thr/1E6,"Mb.pdf"),0.999)

library(GenomicRanges)

define.clusters=function(maps,n.spurious.by.chr,filter.threshold,window.size){
   all.spurious = unlist(n.spurious.by.chr)
   cutoff=quantile(all.spurious,filter.threshold)
   outlier.SNPs = lapply(n.spurious.by.chr,function(n.spurious){
       bad.SNPs = (n.spurious >= cutoff)
       return(which(bad.SNPs))
   })

   all.bad.snps=cbind(unlist(sapply(1:22,function(c){rep(c,length(outlier.SNPs[[c]]))})),unlist(outlier.SNPs))
   all.bad.snps2=GRanges(seqnames=all.bad.snps[,1],ranges=IRanges(start=all.bad.snps[,2],end=all.bad.snps[,2]+1))
   windows=as.data.frame(reduce(all.bad.snps2))
   windows$end=windows$end-1

   windows.pos=NULL
   for(chr in 1:22){
       map=maps[[chr]]
       windows.chr=windows[windows$seqnames==chr,]
       windows.pos=rbind(windows.pos,cbind(map[match(windows.chr$start,map$SNP.n),c("chr","pos","SNP.n")],map[match(windows.chr$end,map$SNP.n),c("chr","pos","SNP.n")]))
   }
   colnames(windows.pos)=c("SNP1.chr","SNP1.pos","SNP1.n","SNP2.chr","SNP2.pos","SNP2.n")
   windows=cbind(windows,windows.pos)
   windows.bp = GRanges(seqnames=windows$seqnames,ranges=IRanges(start=windows$SNP1.pos - window.size,end=windows$SNP2.pos + window.size))
   windows.final=reduce(windows.bp)
   return(windows.final)
}

identify.crossovers.in.clusters=function(r,windows){
    r.2=GRanges(seqnames=r$chr,ranges=IRanges(start=r$START,end=r$END))
    return(subsetByOverlaps(r.2,windows))
}

windows.duplicate.xovers.after.double.xover.removal = define.clusters(maps,freq.duplicate.xovers.after.double.xover.removal,0.999,5E5)
xovers.in.clusters.of.duplicate.xovers.after.double.xover.removal = as.data.frame(identify.crossovers.in.clusters(r.dup.pos.2,windows.duplicate.xovers.after.double.xover.removal))
xovers.in.clusters.of.duplicate.xovers.after.double.xover.removal$xover = paste(xovers.in.clusters.of.duplicate.xovers.after.double.xover.removal$seqnames,
    xovers.in.clusters.of.duplicate.xovers.after.double.xover.removal$start,xovers.in.clusters.of.duplicate.xovers.after.double.xover.removal$end,sep="_")
r.dup.pos.2$xover2=paste0(r.dup.pos.2$xover,r.dup.pos.2$CHILD,r.dup.pos.2$PARENT,sep="_")
r.flt1$xover2=paste0(r.flt1$xover,r.flt1$CHILD,r.flt1$PARENT,sep="_")
r.duplicate.xovers.to.remove = r.dup.pos.2[r.dup.pos.2$xover %in% xovers.in.clusters.of.duplicate.xovers.after.double.xover.removal$xover,]
r.excl.windows.of.duplicate.xovers.after.double.xovers.removal = r.flt1[!r.flt1$xover2 %in% r.duplicate.xovers.to.remove$xover2,]

#define and remove clusters of double xovers within X SNPs   
windows.double.xovers.within.X.SNPs = define.clusters(maps,freq.double.xovers.within.X.SNPs,0.999,5E5)
xovers.in.clusters.of.double.xovers.within.X.SNPs = as.data.frame(identify.crossovers.in.clusters(r.double.xovers.within.X.SNPs,windows.double.xovers.within.X.SNPs))
xovers.in.clusters.of.double.xovers.within.X.SNPs$xover=paste(xovers.in.clusters.of.double.xovers.within.X.SNPs$seqnames,xovers.in.clusters.of.double.xovers.within.X.SNPs$start,
    xovers.in.clusters.of.double.xovers.within.X.SNPs$end,sep="_")
r.double.xovers.within.X.SNPs$xover=paste(r.double.xovers.within.X.SNPs$chr,r.double.xovers.within.X.SNPs$START,r.double.xovers.within.X.SNPs$END,sep="_")
r.double.xovers.within.X.SNPs$xover2=paste(r.double.xovers.within.X.SNPs$xover,r.double.xovers.within.X.SNPs$CHILD,r.double.xovers.within.X.SNPs$PARENT,sep="_")
r$xover3=paste(r$xover,r$CHILD,r$PARENT,sep="_")
r.double.xovers.within.X.SNPs.to.remove=r.double.xovers.within.X.SNPs[r.double.xovers.within.X.SNPs$xover %in% xovers.in.clusters.of.double.xovers.within.X.SNPs$xover,]
r.excl.windows.of.double.xovers.within.X.SNPs = r[!r$xover3 %in% r.double.xovers.within.X.SNPs.to.remove$xover2,]

stopifnot(nrow(r)>0)

write.table(r,gsub("recombinations","recombinations.p_0.5.annotated",gsub("1.","all.",argv[1],fixed=T)),quote=F,sep="\t",row.names=F)
write.table(r.flt,gsub("recombinations","recombinations.p_0.5.annotated.double_xovers_within_1Mb_removed",gsub("1.","all.",argv[1],fixed=T)),quote=F,sep="\t",row.names=F)
write.table(r.flt1,gsub("recombinations","recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_removed",gsub("1.","all.",argv[1],fixed=T)),quote=F,sep="\t",row.names=F)
write.table(r.flt2,gsub("recombinations","recombinations.p_0.5.annotated.duplicate_xovers_removed",gsub("1.","all.",argv[1],fixed=T)),quote=F,sep="\t",row.names=F)
write.table(r.excl.windows.of.duplicate.xovers.after.double.xovers.removal,gsub("recombinations","recombinations.p_0.5.annotated.double_xovers_within_X_SNPs_and_clustered_duplicate_xovers_removed",
                                                                                gsub("1.","all.",argv[1],fixed=T)), quote=F,sep="\t",row.names=F)

#count # of crossovers for informative families, either for all crossovers, after removing double crossovers, or after removing duplicate crossovers
rec <- tabulateRecombination(r,"informative")
#remove crossovers within L Mb
rec.flt <- tabulateRecombination(r.flt,"informative")
#remove crossovers within X SNPs
rec.flt1 <- tabulateRecombination(r.flt1,"informative")
#no duplicate crossovers within or between families
rec.flt2 <- tabulateRecombination(r.flt2,"informative")
#no duplicate crossovers within families
rec.flt3 <- tabulateRecombination(r.flt3,"informative")
#remove crossovers in clusters of double crossovers within X SNPs
rec.flt4 <- tabulateRecombination(r.excl.windows.of.double.xovers.within.X.SNPs,"informative")
#remove crossovers in clusters of duplicate crossovers, after removing crossovers within X SNPs
rec.flt5 <- tabulateRecombination(r.excl.windows.of.duplicate.xovers.after.double.xovers.removal,"informative")


#count # of crossovers for informative families, either for all crossovers, after removing double crossovers, or after removing duplicate crossovers
rec.all <- tabulateRecombination(r,"dummy")
#remove crossovers within L Mb
rec.flt.all <- tabulateRecombination(r.flt,"dummy")
#remove crossovers within X SNPs
rec.flt1.all <- tabulateRecombination(r.flt1,"dummy")
#no duplicate crossovers within or between families
rec.flt2.all <- tabulateRecombination(r.flt2,"dummy")
#no duplicate crossovers within families
rec.flt3.all <- tabulateRecombination(r.flt3,"dummy")
#remove crossovers in clusters of double crossovers within X SNPs
rec.flt4.all <- tabulateRecombination(r.excl.windows.of.double.xovers.within.X.SNPs,"dummy")
#remove crossovers in clusters of duplicate crossovers, after removing crossovers within X SNPs
rec.flt5.all <- tabulateRecombination(r.excl.windows.of.duplicate.xovers.after.double.xovers.removal,"dummy")



write.table(rec.all,gsub("recombinations","recombination_count.p_0.5.annotated",gsub("1.","all.",argv[1],fixed=T)),quote=F,sep="\t",row.names=F)
write.table(rec.flt.all,gsub("recombinations","recombination_count.p_0.5.annotated.double_xovers_within_1Mb_removed",gsub("1.","all.",argv[1],fixed=T)),quote=F,sep="\t",row.names=F)
write.table(rec.flt1.all,gsub("recombinations","recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_removed",gsub("1.","all.",argv[1],fixed=T)),quote=F,sep="\t",row.names=F)
write.table(rec.flt2.all,gsub("recombinations","recombination_count.p_0.5.annotated.duplicate_xovers_removed",gsub("1.","all.",argv[1],fixed=T)),quote=F,sep="\t",row.names=F)
write.table(rec.flt5.all,gsub("recombinations","recombination_count.p_0.5.annotated.double_xovers_within_X_SNPs_and_clustered_duplicate_xovers_removed",
                                                                                gsub("1.","all.",argv[1],fixed=T)), quote=F,sep="\t",row.names=F)
cohort.stats=c(argv[5],"male")
if(nrow(subset(rec,sex=="Male"))>0){
    total=subset(rec,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt,sex=="Male"))>0){
    total=subset(rec.flt,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1,sex=="Male"))>0){
    total=subset(rec.flt1,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt4,sex=="Male"))>0){
    total=subset(rec.flt4,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt5,sex=="Male"))>0){
    total=subset(rec.flt5,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt2,sex=="Male"))>0){
    total=subset(rec.flt2,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt3,sex=="Male"))>0){
    total=subset(rec.flt3,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
male.cohort.stats=cohort.stats
cohort.stats=c(argv[5],"female")
if(nrow(subset(rec,sex=="Female"))>0){
    total=subset(rec,sex=="Female")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt,sex=="Female"))>0){
    total=subset(rec.flt,sex=="Female")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1,sex=="Female"))>0){
    total=subset(rec.flt1,sex=="Female")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt4,sex=="Female"))>0){
    total=subset(rec.flt4,sex=="Female")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt5,sex=="Female"))>0){
    total=subset(rec.flt5,sex=="Female")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt2,sex=="Female"))>0){
    total=subset(rec.flt2,sex=="Female")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt3,sex=="Female"))>0){
    total=subset(rec.flt3,sex=="Female")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}

stats=rbind(male.cohort.stats,cohort.stats)
colnames(stats)=c("cohort","sex",sapply(c("all","no_double_xovers_within_N_MB","no_double_xovers_within_X_SNPs","no_double_xovers_within_clusters_of_double_crossovers_within_X_SNPs",
            "no_double_xovers_within_X_SNPs_or_duplicate_crossovers_within_clusters","no_duplicate_xovers","no_duplicate_xovers_within_fams"),
            function(x){paste(x,c("n_duos","mean","min","25th_pc","median","75th_pc","max","stdev"),sep=".")}))

write.table(stats,paste0(argv[3],"_summary.different_filtering_strategies.informative_families_only.txt"),row.names=F,quote=F,sep="\t")

#### Make QQ plots for different filtering schemes
#stopifnot(nrow(rec)>0)
if(nrow(rec)>0){
print("QQ-PLOT")
pdf(paste0(argv[3],"-qq.pdf"),width=12,height=6)
#par(mar=c(2,2,1,1),mfrow=c(1,2),cex=1.6,oma=c(2,2,0,0))
par(mar=c(2,2,1,1),mfrow=c(1,2),cex=1,oma=c(2,2,0,0))
tmp <- qqr(subset(rec,sex=="Male")$nrec,malesize,cols[1],title=paste(argv[5]," paternal"),
           limits=range(c(subset(rec,sex=="Male")$nrec,subset(rec.flt,sex=="Male")$nrec,subset(rec.flt1,sex=="Male")$nrec,subset(rec.flt2,sex=="Male")$nrec)),xlowerlimit=10)
#points(jitter(qqr(subset(rec.flt,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[2])
points(jitter(qqr(subset(rec.flt1,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[2])
points(jitter(qqr(subset(rec.flt2,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[3])
points(jitter(qqr(subset(rec.flt3,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[4])
#points(jitter(qqr(subset(rec.flt4,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[5])
points(jitter(qqr(subset(rec.flt5,sex=="Male")$nrec,malesize,makeplot=0)),pch=16,col=cols[5])
#legend("topleft",c("Unfiltered",paste0("<",thr/1000,"kb filtered")," < X SNPs apart filtered","Duplicate positions removed"),col=cols[c(1:2,4,3)],pch=16,bg="white")
legend("topleft",c("Unfiltered","Filter 1: double xovers < X SNPs apart filtered","Filter 1 + clustered duplicate positions removed","Filter 3: Duplicate positions removed",
                   "Filter4: Duplicate positions within families removed"),col=cols[c(1,2,5,3,4)],pch=16,
       bg="white")
tmp <- qqr(subset(rec,sex=="Female")$nrec,femalesize,cols[1],title=paste(argv[5]," maternal"),
           limits=range(c(subset(rec,sex=="Female")$nrec,subset(rec.flt,sex=="Female")$nrec,subset(rec.flt1,sex=="Female")$nrec,subset(rec.flt2,sex=="Female")$nrec)),xlowerlimit=20)
#points(jitter(qqr(subset(rec.flt,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[2])
points(jitter(qqr(subset(rec.flt1,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[2])
points(jitter(qqr(subset(rec.flt2,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[3])
points(jitter(qqr(subset(rec.flt3,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[4])
#points(jitter(qqr(subset(rec.flt4,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[5])
points(jitter(qqr(subset(rec.flt5,sex=="Female")$nrec,femalesize,makeplot=0)),pch=16,col=cols[5])
mtext("Expected recombinations genome-wide",1,outer=T,cex=1.6,padj=1.7)
mtext("Observed recombinations genome-wide",2,outer=T,cex=1.6,padj=-1)
dev.off()
}



#rec.flt1=rec.flt5

rec.flt1.all= tabulateRecombination(r.flt1,"dummy")
rec.flt1.informative.2gen = tabulateRecombination(r.flt1,"informative.2gen")
rec.flt1.informative.2gen.2parents = tabulateRecombination(r.flt1,"informative.2gen.2parents")
rec.flt1.informative.2gen.1parent = tabulateRecombination(r.flt1,"informative.2gen.1parent")
rec.flt1.noninf.2kids.2gen.2parents = tabulateRecombination(r.flt1,"noninf.2kids.2gen.2parents")
rec.flt1.noninf.2kids.2gen.1parent = tabulateRecombination(r.flt1,"noninf.2kids.2gen.1parent")
rec.flt1.noninf.1kid.2gen.2parents = tabulateRecombination(r.flt1,"noninf.1kid.2gen.2parents")
rec.flt1.noninf.1kid.2gen.1parent = tabulateRecombination(r.flt1,"noninf.1kid.2gen.1parent")
rec.flt1.informative.3gen = tabulateRecombination(r.flt1,"informative.3gen")
rec.flt1.informative.3gen.2parents = tabulateRecombination(r.flt1,"informative.3gen.2parents")
rec.flt1.informative.3gen.1parent = tabulateRecombination(r.flt1,"informative.3gen.1parent")

save.image(paste0(argv[3],".multiple_filters.RData"))  
#if(FALSE){
cohort.stats=c(argv[5],"male")
if(nrow(subset(rec.flt1.all,sex=="Male"))>0){
    total=subset(rec.flt1.all,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1,sex=="Male"))>0){
    total=subset(rec.flt1,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.informative.2gen.2parents,sex=="Male"))>0){
    total=subset(rec.flt1.informative.2gen.2parents,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.informative.2gen.1parent,sex=="Male"))>0){
    total=subset(rec.flt1.informative.2gen.1parent,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Male"))>0){
    total=subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.noninf.2kids.2gen.1parent,sex=="Male"))>0){
    total=subset(rec.flt1.noninf.2kids.2gen.1parent,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Male"))>0){
    total=subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Male"))>0){
    total=subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.informative.3gen.2parents,sex=="Male"))>0){
    total=subset(rec.flt1.informative.3gen.2parents,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.informative.3gen.1parent,sex=="Male"))>0){
    total=subset(rec.flt1.informative.3gen.1parent,sex=="Male")$nrec
    cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
} else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}

male.cohort.stats=cohort.stats
cohort.stats=c(argv[5],"female")
if(nrow(subset(rec.flt1.all,sex=="Female"))>0){
        total=subset(rec.flt1.all,sex=="Female")$nrec
            cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
    } else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1,sex=="Female"))>0){
        total=subset(rec.flt1,sex=="Female")$nrec
            cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
    } else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.informative.2gen.2parents,sex=="Female"))>0){
        total=subset(rec.flt1.informative.2gen.2parents,sex=="Female")$nrec
            cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
    } else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.informative.2gen.1parent,sex=="Female"))>0){
        total=subset(rec.flt1.informative.2gen.1parent,sex=="Female")$nrec
            cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
    } else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Female"))>0){
        total=subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Female")$nrec
            cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
    } else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.noninf.2kids.2gen.1parent,sex=="Female"))>0){
        total=subset(rec.flt1.noninf.2kids.2gen.1parent,sex=="Female")$nrec
            cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
    } else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Female"))>0){
        total=subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Female")$nrec
            cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
    } else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Female"))>0){
        total=subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Female")$nrec
            cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
    } else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.informative.3gen.2parents,sex=="Female"))>0){
        total=subset(rec.flt1.informative.3gen.2parents,sex=="Female")$nrec
            cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
    } else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}
if(nrow(subset(rec.flt1.informative.3gen.1parent,sex=="Female"))>0){
        total=subset(rec.flt1.informative.3gen.1parent,sex=="Female")$nrec
            cohort.stats=c(cohort.stats,length(total),mean(total),quantile(total,c(0,0.25,0.5,0.75,1)),sd(total))
    } else {    cohort.stats=c(cohort.stats,0,NA,NA,NA,NA,NA,NA,NA)}

stats=rbind(male.cohort.stats,cohort.stats)
colnames(stats)=c("cohort","sex",sapply(c("all_families","informative","informative.2gen.2parents","informative.2gen.1parent","noninf.2kids.2gen.2parents","noninf.2kids.2gen.1parent","noninf.1kid.2gen.2parents",
            "noninf.1kid.2gen.1parent","informative.3gen.2parents","informative.3gen.1parent"),function(x){paste(x,c("n_duos","mean","min","25th_pc","median","75th_pc","max","stdev"),sep=".")}))

write.table(stats,paste0(argv[3],"_summary.different_family_types.X_SNPs_apart_filtered.txt"),row.names=F,quote=F,sep="\t")

mycols=c("darkgreen","green","red3","red","purple","purple3","orange","orange2")
mypch=c(21,24,21,24,21,24,21,24)

print("QQ-PLOT")
pdf(paste0(argv[3],"-qq_by_family_type.X_SNPs_apart_filtered.pdf"),width=12,height=6)
par(mar=c(2,2,1,1),mfrow=c(1,2),cex=1,oma=c(2,2,0,0))
tmp <- qqr(subset(rec.flt1.all,sex=="Male")$nrec,malesize,"white",title=paste(argv[5]," paternal"),limits=range(subset(rec.flt1.all,sex=="Male")$nrec),xlowerlimit=10)
if(nrow(subset(rec.flt1.informative.2gen.2parents,sex="Male"))>0){points(jitter(qqr(subset(rec.flt1.informative.2gen.2parents,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[1],col=mycols[1],bg=mycols[1])}
if(nrow(subset(rec.flt1.informative.2gen.1parent,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.informative.2gen.1parent,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[2],col=mycols[2],bg=mycols[2])}
if(nrow(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[3],col=mycols[3],bg=mycols[3])}
if(nrow(subset(rec.flt1.noninf.2kids.2gen.1parent,sex="Male"))>0){points(jitter(qqr(subset(rec.flt1.noninf.2kids.2gen.1parent,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[4],col=mycols[4],bg=mycols[4])}
if(nrow(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[5],col=mycols[5],bg=mycols[5])}
if(nrow(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[6],col=mycols[6],bg=mycols[6])}
if(nrow(subset(rec.flt1.informative.3gen.2parents,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.informative.3gen.2parents,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[7],col=mycols[7],bg=mycols[7])}
if(nrow(subset(rec.flt1.informative.3gen.1parent,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.informative.3gen.1parent,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[8],col=mycols[8],bg=mycols[8])}
legend("topleft",c("3 kids, both parents","3 kids, one parent","2 kids, both parents","2 kids, one parent","1 kid, both parents","1 kid, one parent","3 generations, both parents","3 generations, one parent"),
       pch=mypch,col=mycols,cex=0.8,pt.bg=mycols)
#legend("topleft",c("Unfiltered",paste0("<",thr/1000,"kb filtered")," < X SNPs apart filtered","Duplicate positions removed"),col=cols[c(1:2,4,3)],pch=16,bg="white")
tmp <- qqr(subset(rec.flt1.all,sex=="Female")$nrec,femalesize,"white",title=paste(argv[5]," maternal"),limits=range(subset(rec.flt1.all,sex=="Female")$nrec),xlowerlimit=20)
if(nrow(subset(rec.flt1.informative.2gen.2parents,sex="Female"))>0){points(jitter(qqr(subset(rec.flt1.informative.2gen.2parents,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[1],col=mycols[1],bg=mycols[1])}
if(nrow(subset(rec.flt1.informative.2gen.1parent,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.informative.2gen.1parent,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[2],col=mycols[2],bg=mycols[2])}
if(nrow(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[3],col=mycols[3],bg=mycols[3])}
if(nrow(subset(rec.flt1.noninf.2kids.2gen.1parent,sex="Female"))>0){points(jitter(qqr(subset(rec.flt1.noninf.2kids.2gen.1parent,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[4],col=mycols[4],bg=mycols[4])}
if(nrow(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[5],col=mycols[5],bg=mycols[5])}
if(nrow(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[6],col=mycols[6],bg=mycols[6])}
if(nrow(subset(rec.flt1.informative.3gen.2parents,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.informative.3gen.2parents,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[7],col=mycols[7],bg=mycols[7])}
if(nrow(subset(rec.flt1.informative.3gen.1parent,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.informative.3gen.1parent,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[8],col=mycols[8],bg=mycols[8])}
   mtext("Expected recombinations genome-wide",1,outer=T,cex=1.6,padj=1.7)
mtext("Observed recombinations genome-wide",2,outer=T,cex=1.6,padj=-1)
dev.off()
#}

####histograms
pdf(paste0(argv[3],"-histograms_by_family_type.X_SNPs_apart_filtered.pdf"),width=8,height=8)
par(mfrow=c(2,2))
if(nrow(subset(rec.flt1.informative.2gen.2parents,sex="Male"))>0){hist(subset(rec.flt1.informative.2gen.2parents,sex=="Male")$nrec,col="lightblue",xlab="Number of recombinations",main="male, >2 kids")
    abline(v=25.9,col="darkgreen",lwd=2)
    abline(v=26.6,col="blue",lwd=2)
    abline(v=mean(subset(rec.flt1.informative.2gen.2parents,sex=="Male")$nrec),col="black",lwd=2)
     legend("topleft",c("deCODE","Adam","these data"),col=c("darkgreen","blue","black"),lwd=2,lty=1)
}else {par(mfg=c(1,2))}

if(nrow(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Male"))>0){hist(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Male")$nrec,col="lightblue",xlab="Number of recombinations",main="male, 2 kids")
    abline(v=25.9,col="darkgreen",lwd=2)
    abline(v=26.6,col="blue",lwd=2)
    abline(v=mean(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Male")$nrec),col="black",lwd=2)
    abline(v=mean(subset(rec.flt1.informative.2gen.2parents,sex=="Male")$nrec),col="black",lwd=2,lty=2)
    legend("topleft",">2 kids",col="black",lwd=2,lty=2)
    }else {par(mfg=c(2,1))}

if(nrow(subset(rec.flt1.informative.2gen.2parents,sex="Female"))>0){hist(subset(rec.flt1.informative.2gen.2parents,sex=="Female")$nrec,col="pink",xlab="Number of recombinations",main="female, >2 kids")
    abline(v=42.81,col="darkgreen",lwd=2)
    abline(v=41.6,col="blue",lwd=2)
    abline(v=mean(subset(rec.flt1.informative.2gen.2parents,sex=="Female")$nrec),col="black",lwd=2)
                                     }else {par(mfg=c(2,2))}
if(nrow(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Female"))>0){hist(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Female")$nrec,col="pink",xlab="Number of recombinations",main="female, 2 kids")
    abline(v=42.81,col="darkgreen",lwd=2)
abline(v=41.6,col="blue",lwd=2)
abline(v=mean(subset(rec.flt1.noninf.2kids.2gen.2parents,sex=="Female")$nrec),col="black",lwd=2)
abline(v=mean(subset(rec.flt1.informative.2gen.2parents,sex=="Female")$nrec),col="black",lwd=2,lty=2)
    legend("topleft",">2 kids",col="black",lwd=2,lty=2)

}                                                                      
                                                                  
dev.off()

if(FALSE){
if(nrow(subset(rec.flt1.informative.2gen.1parent,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.informative.2gen.1parent,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[2],col=mycols[2],bg=mycols[2])}
if(nrow(subset(rec.flt1.noninf.2kids.2gen.1parent,sex="Male"))>0){points(jitter(qqr(subset(rec.flt1.noninf.2kids.2gen.1parent,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[4],col=mycols[4],bg=mycols[4])}
if(nrow(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[5],col=mycols[5],bg=mycols[5])}
if(nrow(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[6],col=mycols[6],bg=mycols[6])}
if(nrow(subset(rec.flt1.informative.3gen.2parents,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.informative.3gen.2parents,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[7],col=mycols[7],bg=mycols[7])}
if(nrow(subset(rec.flt1.informative.3gen.1parent,sex=="Male"))>0){points(jitter(qqr(subset(rec.flt1.informative.3gen.1parent,sex=="Male")$nrec,malesize,makeplot=0)),pch=mypch[8],col=mycols[8],bg=mycols[8])}



if(nrow(subset(rec.flt1.informative.2gen.1parent,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.informative.2gen.1parent,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[2],col=mycols[2],bg=mycols[2])}
if(nrow(subset(rec.flt1.noninf.2kids.2gen.1parent,sex="Female"))>0){points(jitter(qqr(subset(rec.flt1.noninf.2kids.2gen.1parent,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[4],col=mycols[4],bg=mycols[4])}
if(nrow(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.noninf.1kid.2gen.2parents,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[5],col=mycols[5],bg=mycols[5])}
if(nrow(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.noninf.1kid.2gen.1parent,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[6],col=mycols[6],bg=mycols[6])}
if(nrow(subset(rec.flt1.informative.3gen.2parents,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.informative.3gen.2parents,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[7],col=mycols[7],bg=mycols[7])}
if(nrow(subset(rec.flt1.informative.3gen.1parent,sex=="Female"))>0){points(jitter(qqr(subset(rec.flt1.informative.3gen.1parent,sex=="Female")$nrec,femalesize,makeplot=0)),pch=mypch[8],col=mycols[8],bg=mycols[8])}
dev.off()

#print(head(rec))
print(paste0(argv[3],".rec"))
rec <- rec[order(rec$sex,rec$nrec),]
write.table(x=rec,file=paste0(argv[3],".rec"),quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")


r <- subset(r,informative)
if(nrow(r)>0){
nmale   <- sum(r[!duplicated(r[,1:2]),"sex"]=="Male")
nfemale <- sum(r[!duplicated(r[,1:2]),"sex"]=="Female")

pdf(paste0(argv[3],"-perchromosome.pdf"),width=12,height=6)

par(mar=c(2,2,1,1),mfrow=c(1,2),cex=1.6,oma=c(2,2,0,0))
xhat <- seq(0,10,.01)
y <- table(subset(r,sex=="Male")$chr)/nmale
l <- c(0,max(male.lengths,y))
plot(male.lengths,y,ylab="Mean per chromosome",xlab="Expected per chromosome",type='n',xlim=l,ylim=l)
mtext(paste("Paternal (n =",nmale,")"),3,cex=1.6)
y1 <- xhat - 1.96*sqrt(xhat)/sqrt(nmale)
y2 <- xhat + 1.96*sqrt(xhat)/sqrt(nmale)
polygon(c(xhat,rev(xhat)),c(y1,rev(y2)),col='light grey',border=NA)
abline(0,1,col=1,lty=2)
text(male.lengths,y,labels=paste(1:22),col=cols[1])

y <- table(subset(r,sex=="Female")$chr)/nfemale
l <- c(0,max(female.lengths,y))
plot(female.lengths,y,ylab="Mean per chromosome",xlab="Expected per chromosome",type='n',xlim=l,ylim=l)
mtext(paste("Maternal (n =",nfemale,")"),3,cex=1.6)
y1 <- xhat - 1.96*sqrt(xhat)/sqrt(nfemale)
y2 <- xhat + 1.96*sqrt(xhat)/sqrt(nfemale)
polygon(c(xhat,rev(xhat)),c(y1,rev(y2)),col='light grey',border=NA)
abline(0,1,col=1,lty=2)
text(female.lengths,y,labels=paste(1:22),col=cols[1])

mtext("Expected recombinations per chromosome",1,outer=T,cex=1.6,padj=1.7)
mtext("Observed recombinations per chromosome",2,outer=T,cex=1.6,padj=-1)

#legend("topleft",c("SHAPEIT2","Merlin"),col=cols[1:2],pch=16,bg='white')
dev.off()
}

#save.image(paste0(argv[3],".multiple_filters.RData"))
                                        #save.image(paste0(argv[3],".RData"))
}
