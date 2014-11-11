options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

recomb.file=args[1]
sample.file = args[2]#sample.file="SHAPEIT_results/FVG.chr1.generr_removed.SHAPEIT_phased.sample"
p.cutoff = as.numeric(args[3])
out.prefix = args[4]#out.prefix="duoHMM_results/FVG.generr_removed.total_recombination_counts"

recomb.file="duoHMM_results/QTR_minus_MZs.370K.1.post_MERLIN.max_0.05_missing.generr_removed.recombinations.txt"
sample.file="SHAPEIT_results/QTR_minus_MZs.370K.1.post_MERLIN.max_0.05_missing.generr_removed.SHAPEIT_phased.sample"
p.cutoff=0.5
out.prefix="duoHMM_results/QTR.370K.post_MERLIN.max_0.05_missing.generr_removed.no_adjacent_crossovers"  

samples=read.delim(sample.file,header=F,sep="",stringsAsFactors=F,skip=2)
samples=samples[,-c(3,7)]
samples2=cbind(as.character(samples[,1]),as.character(samples[,2]),as.character(samples[,3]),as.character(samples[,4]),as.character(samples[,5]))
rownames(samples2) = samples2[,2]
samples2=as.data.frame(samples2,stringsAsFactors=F)
#remove parents not genotyped
samples2[! samples[,3] %in% samples[,2],3]="0"
samples2[! samples[,4] %in% samples[,2],4]="0"

samples2$total_fam_size=sapply(samples2[,1],function(x){sum(samples2[,1]==x)})
samples2$total_same_father=sapply(samples2[,2],function(x){sum(samples2[,3]==samples2[x,3])})
samples2$total_same_mother=sapply(samples2[,2],function(x){sum(samples2[,4]==samples2[x,4])})
samples2$total_same_parents=sapply(samples2[,2],function(x){sum(samples2[,4]==samples2[x,4] & samples2[,3]==samples2[x,3])})
samples2$total_same_father[samples2[,3]=="0"]=0
samples2$total_same_mother[samples2[,4]=="0"]=0
samples2$total_same_parents[samples2[,4]=="0" & samples2[,3]=="0"]=0

samples2$maternal_grandmother="0"
samples2$maternal_grandmother[samples2[,4] %in% samples2[,2]] = samples2[samples2[samples2[,4]%in% samples2[,2],4],4]
samples2$maternal_grandfather="0"
samples2$maternal_grandfather[samples2[,4] %in% samples2[,2]] = samples2[samples2[samples2[,4]%in% samples2[,2],4],3]

samples2$paternal_grandmother="0"
samples2$paternal_grandmother[samples2[,3] %in% samples2[,2]] = samples2[samples2[samples2[,3]%in% samples2[,2],3],4]
samples2$paternal_grandfather="0"
samples2$paternal_grandfather[samples2[,3] %in% samples2[,2]] = samples2[samples2[samples2[,3]%in% samples2[,2],3],3]

samples2$n_genotyped_maternal_grandparents =rowSums(samples2[,c("maternal_grandmother","maternal_grandfather")]!="0")
samples2$n_genotyped_paternal_grandparents =rowSums(samples2[,c("paternal_grandmother","paternal_grandfather")]!="0")
samples2$n_genotyped_grandparents =samples2$n_genotyped_maternal_grandparents+samples2$n_genotyped_paternal_grandparents
samples2$n_genotyped_children = sapply(samples2[,2],function(x){sum(samples[,3] == x | samples[,4] == x)})

halfsibs=samples2[(samples2$total_same_parents!=samples2$total_same_mother&samples2[,4]!="0")|(samples2$total_same_parents!=samples2$total_same_father & samples2[,3]!="0"),]
samples2$possible_halfsib=((samples2$total_same_parents!=samples2$total_same_mother&samples2[,4]!="0")|(samples2$total_same_parents!=samples2$total_same_father & samples2[,3]!="0"))
colnames(samples2)[1:5]=c("family","child","father","mother","sex")
if(p.cutoff==0){
write.table(halfsibs,paste(out.prefix,".annotated_as_halfsibs.txt",sep=""),quote=F,sep="\t",row.names=F)
write.table(samples2,paste0(out.prefix,"_sample_information.txt"),quote=F,sep="\t",row.names=F)
}


maternal=samples2[samples2[,4] != "0",c(1,2,4,3,6:ncol(samples2))]
paternal=samples2[samples2[,3] != "0",c(1,2,3,4,6:ncol(samples2))]

rownames(maternal)=maternal[,2]
rownames(paternal)=paternal[,2]
colnames(maternal)[1:4]=c("family","child","mother","father")
colnames(paternal)[1:4]=c("family","child","father","mother")

all.recomb=NULL
keep.adjacent=NULL
for(i in 1:22){
recomb.file2=gsub("1.",paste(i,".",sep=""),recomb.file,fixed=T)

recomb=read.delim(recomb.file2,header=T,sep="\t",stringsAsFactors=F)
#recomb=read.delim(recomb.file2,header=F,sep="\t",stringsAsFactors=F)
#colnames(recomb)=c("CHILD","PARENT","START","END","PROB_RECOMBINATION")
#recomb=recomb[order(recomb$CHILD,recomb$PARENT),]
recomb$chr=i
#recomb$adjacent=c(FALSE,(recomb$START[2:nrow(recomb)]==recomb$END[1:(nrow(recomb)-1)] & recomb$CHILD[2:nrow(recomb)]==recomb$CHILD[1:(nrow(recomb)-1)]))
#recomb$adjacent=c(FALSE,(recomb$START[2:nrow(recomb)]==recomb$END[1:(nrow(recomb)-1)]))
#keep.adjacent=rbind(keep.adjacent,recomb[sort(c(which(recomb$adjacent),which(recomb$adjacent)-1)),])
#print(recomb[recomb$adjacent,])

recomb2=recomb[recomb$PROB_RECOMBINATION > p.cutoff,]
#recomb2=recomb2[order(recomb2$CHILD,recomb2$PARENT),]
#recomb2$adjacent=c(FALSE,(recomb2$START[2:nrow(recomb2)]==recomb2$END[1:(nrow(recomb2)-1)] & recomb2$CHILD[2:nrow(recomb2)]==recomb2$CHILD[1:(nrow(recomb2)-1)]))
#keep.adjacent=rbind(keep.adjacent,recomb2[sort(c(which(recomb2$adjacent),which(recomb2$adjacent)-1)),])
#print(recomb2[recomb2$adjacent,])

all.recomb=rbind(all.recomb,recomb2)
#if(FALSE){
maternal=data.frame(maternal,rep(0,nrow(maternal)),stringsAsFactors=F)
colnames(maternal)[ncol(maternal)]=paste("chr",i,sep="")
maternal.counts=apply(maternal[,2:3],1,function(indivs){sum(recomb2[,1]==indivs[1] & recomb2[,2]==indivs[2])})
maternal[names(maternal.counts),ncol(maternal)]=maternal.counts

paternal=data.frame(paternal,rep(0,nrow(paternal)),stringsAsFactors=F)
colnames(paternal)[ncol(paternal)]=paste("chr",i,sep="")
paternal.counts=apply(paternal[,2:3],1,function(indivs){sum(recomb2[,1]==indivs[1] & recomb2[,2]==indivs[2])})
paternal[names(paternal.counts),ncol(paternal)]=paternal.counts
#}
}
#write.table(keep.adjacent,"old_data/GPC_old_duoHMM_data.adjacent_crossovers.p_0.txt",quote=F,sep="\t",row.names=F)
#write.table(keep.adjacent,"old_data/GPC_old_duoHMM_data.adjacent_crossovers.p_0.5.txt",quote=F,sep="\t",row.names=F)
#write.table(keep.adjacent,"duoHMM_results/QTR.610K.post_MERLIN.max_0.05_missing.generr_removed.adjacent_crossovers.txt",quote=F,sep="\t",row.names=F)

write.table(all.recomb,paste0(out.prefix,".all_crossovers.p_",p.cutoff,".txt"),quote=F,sep="\t",row.names=F)

maternal.recomb=all.recomb[all.recomb$PARENT %in% maternal$mother,]
paternal.recomb=all.recomb[all.recomb$PARENT %in% paternal$father,]

maternal.recomb.inform=all.recomb[all.recomb$PARENT %in% maternal[maternal$n_genotyped_maternal_grandparents >0 |maternal$total_same_mother>2 & maternal$total_same_father>2 & maternal$total_same_mother==maternal$total_same_father,"mother"],]
paternal.recomb.inform=all.recomb[all.recomb$PARENT %in% paternal[paternal$n_genotyped_paternal_grandparents >0 |paternal$total_same_father>2 & paternal$total_same_mother>2 & paternal$total_same_father==paternal$total_same_mother,"father"],]

write.table(maternal.recomb,paste0(out.prefix,".maternal_crossovers.p_",p.cutoff,".txt"),quote=F,sep="\t",row.names=F)
write.table(paternal.recomb,paste0(out.prefix,".paternal_crossovers.p_",p.cutoff,".txt"),quote=F,sep="\t",row.names=F)
write.table(maternal.recomb.inform,paste0(out.prefix,".maternal_crossovers_from_informative_meioses.p_",p.cutoff,".txt"),quote=F,sep="\t",row.names=F)
write.table(paternal.recomb.inform,paste0(out.prefix,".paternal_crossovers_from_informative_meioses.p_",p.cutoff,".txt"),quote=F,sep="\t",row.names=F)

maternal$total=rowSums(maternal[,grep("chr",colnames(maternal))])
paternal$total=rowSums(paternal[,grep("chr",colnames(paternal))])


write.table(maternal,paste(out.prefix,".total_recombination_counts.maternal.p_",p.cutoff,".txt",sep=""),quote=F,sep="\t",row.names=F)
write.table(paternal,paste(out.prefix,".total_recombination_counts.paternal.p_",p.cutoff,".txt",sep=""),quote=F,sep="\t",row.names=F)



