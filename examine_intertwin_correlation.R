age.files=c("NTR/NTR_v2_minus_MZs.with_parental_age.txt",
        "QTR/QTR_minus_MZs.370K.with_parental_age.txt",
        "QTR/QTR_minus_MZs.610K.with_parental_age.txt")

fams=c("NTR/duoHMM_results/NTR.more_stringent.sample_information.txt","QTR/duoHMM_results/QTR370.more_stringent.sample_information.txt",
    "QTR/duoHMM_results/QTR610.more_stringent.sample_information.txt")
    
names(age.files)=c("NTR","QTR370","QTR610")

twin.info=list()
for(i in 1:length(age.files)){
    ages=read.delim(age.files[i],header=T,colClasses=c(rep("character",4),"integer","integer","numeric","numeric"))
    ages$unique.id = paste(ages$family,ages$father,ages$mother,ages$maternal_age,ages$paternal_age,sep="_")
    counts=table(ages$unique.id)
    ages$count=counts[ages$unique.id]
   
    twins=ages[ages$count==2 & !(ages$father == "0" & ages$mother=="0"),]
#    twins=twins[order(twins$family),]
    twins=twins[order(twins$mother),] #need to order this way else we get errors for families with >1 pair of twins
    
    twin1=twins[1:nrow(twins) /2 -1:nrow(twins)%/% 2==0.5,]
    twin2=twins[1:nrow(twins) /2 -1:nrow(twins)%/% 2==0,]
    
    twins=data.frame(cbind(twin1[,c(1:5,7:8)],twin2[,c(2,5)]))
    colnames(twins)=c("family","twin1","father","mother","twin1.sex","twin.mat.age","twin.pat.age","twin2","twin2.sex")
    
    sibs=ages[ages$count==1 & !(ages$father == "0" & ages$mother=="0"),]
    nsibs=c()
    for(j in 1:nrow(twin1)){
        nsibs=c(nsibs,length(sibs[sibs$father==twin1[j,"father"] & sibs$mother==twin1[j,"mother"],"individual"]))
    }
    twins$nsibs=nsibs
    max.n.sibs=max(nsibs)
    for(c in 1:max.n.sibs){
        dummy=matrix(NA,nrow=nrow(twins),ncol=4)
        colnames(dummy)=paste0("sib",c,c("",".sex",".mat.age",".pat.age"))
        twins=cbind(twins,dummy)
    }
    for(j in 1:nrow(twin1)){
        these.sibs=    sibs[sibs$father==twin1[j,"father"] & sibs$mother==twin1[j,"mother"],c("individual","sex","maternal_age","paternal_age")]
        if(nrow(these.sibs)>0){
            for(c in 1:nrow(these.sibs)){
                twins[j,paste0("sib",c,c("",".sex",".mat.age",".pat.age"))]=these.sibs[c,]
            }
        }
    }
    twin.info[[i]] = twins
}
names(twin.info)=c("NTR","QTR370","QTR610")

twin.pairs=rbind(twin.info[[1]][,c("twin1","twin2")],twin.info[[2]][,c("twin1","twin2")],twin.info[[3]][,c("twin1","twin2")])
twin.pairs1=twin.pairs
twin.pairs2=twin.pairs[,c(2,1)]
colnames(twin.pairs2)  = colnames(twin.pairs1)
twin.pairs=rbind(twin.pairs1,twin.pairs2)




load("NFTOOLS_data_for_RSTAN.more_stringent.RData")
nftools.data1.mat=data1.mat
nftools.data1.pat=data1.pat

load("duoHMM_data_for_RSTAN.more_stringent.RData")
duohmm.data0.mat=data0.mat # inform nuc fams
duohmm.data1.mat=data1.mat # inform meioses
duohmm.data2.mat=data2.mat # all meioses from families with >1 kid

duohmm.data0.pat=data0.pat # inform nuc fams
duohmm.data1.pat=data1.pat # inform meioses
duohmm.data2.pat=data2.pat # all meioses from families with >1 kid


duohmm.data1.mat$other.twin=NA
duohmm.data1.mat$other.twin = twin.pairs[match(duohmm.data1.mat[,"CHILD"],twin.pairs[,"twin1"]),"twin2"]
duohmm.data1.mat$other.twin[!(duohmm.data1.mat$other.twin %in% duohmm.data1.mat$CHILD)] = NA


duohmm.data1.pat$other.twin=NA
duohmm.data1.pat$other.twin = twin.pairs[match(duohmm.data1.pat[,"CHILD"],twin.pairs[,"twin1"]),"twin2"]
duohmm.data1.pat$other.twin[!(duohmm.data1.pat$other.twin %in% duohmm.data1.pat$CHILD)] = NA


duohmm.data1.mat.simple=duohmm.data1.mat[,c("CHILD","PARENT","nrec","age.at.birth","cohort","other.twin")]
duohmm.data1.pat.simple=duohmm.data1.pat[,c("CHILD","PARENT","nrec","age.at.birth","cohort","other.twin")]
colnames(duohmm.data1.mat.simple)[3] = "n.crossovers"
colnames(duohmm.data1.pat.simple)[3] = "n.crossovers"

write.table(duohmm.data1.mat.simple,"maternal_crossovers.informative_meioses.with_twins_indicated.txt",quote=F,sep="\t",row.names=F)
write.table(duohmm.data1.pat.simple,"paternal_crossovers.informative_meioses.with_twins_indicated.txt",quote=F,sep="\t",row.names=F)

twin.rec.info=list()
for(i in 1:length(twin.info)){
    twins=twin.info[[i]]
    twins$inf.nuc.fam.mat=twins$twin1 %in% duohmm.data0.mat$CHILD
    twins$inf.mat=twins$twin1 %in% duohmm.data1.mat$CHILD
    twins$multikid.mat = twins$twin1 %in% duohmm.data2.mat$CHILD

    twins$inf.nuc.fam.pat=twins$twin1 %in% duohmm.data0.pat$CHILD
    twins$inf.pat=twins$twin1 %in% duohmm.data1.pat$CHILD
    twins$multikid.pat = twins$twin1 %in% duohmm.data2.pat$CHILD
    
    twins$twin1.duohmm.mat.nrec=duohmm.data2.mat[match(twins$twin1,duohmm.data2.mat$CHILD),"nrec"]
    twins$twin2.duohmm.mat.nrec=duohmm.data2.mat[match(twins$twin2,duohmm.data2.mat$CHILD),"nrec"]

    for(j in 1:max(twins$nsibs)){
        nrec = duohmm.data2.mat[match(twins[,paste0("sib",j)],duohmm.data2.mat$CHILD),"nrec"]
        twins=cbind(twins,nrec)
        colnames(twins)[ncol(twins)] = paste0("sib",j,".duohmm.mat.nrec")
    }
    twins$twin1.duohmm.pat.nrec=duohmm.data2.pat[match(twins$twin1,duohmm.data2.pat$CHILD),"nrec"]
    twins$twin2.duohmm.pat.nrec=duohmm.data2.pat[match(twins$twin2,duohmm.data2.pat$CHILD),"nrec"]
    for(j in 1:max(twins$nsibs)){
        nrec = duohmm.data2.pat[match(twins[,paste0("sib",j)],duohmm.data2.pat$CHILD),"nrec"]
        twins=cbind(twins,nrec)
        colnames(twins)[ncol(twins)] = paste0("sib",j,".duohmm.pat.nrec")
    }
    twins$twin1.nftools.mat.nrec=nftools.data1.mat[match(twins$twin1,nftools.data1.mat$child),"Freq"]
    twins$twin2.nftools.mat.nrec=nftools.data1.mat[match(twins$twin2,nftools.data1.mat$child),"Freq"]
    for(j in 1:max(twins$nsibs)){
        nrec = nftools.data1.mat[match(twins[,paste0("sib",j)],nftools.data1.mat$child),"Freq"]
        twins=cbind(twins,nrec)
        colnames(twins)[ncol(twins)] = paste0("sib",j,".nftools.mat.nrec")
    }
    twins$twin1.nftools.pat.nrec=nftools.data1.pat[match(twins$twin1,nftools.data1.pat$child),"Freq"]
    twins$twin2.nftools.pat.nrec=nftools.data1.pat[match(twins$twin2,nftools.data1.pat$child),"Freq"]
    for(j in 1:max(twins$nsibs)){
        nrec = nftools.data1.pat[match(twins[,paste0("sib",j)],nftools.data1.pat$child),"Freq"]
        twins=cbind(twins,nrec)
        colnames(twins)[ncol(twins)] = paste0("sib",j,".nftools.pat.nrec")
    }
twin.rec.info[[i]]=twins    
}
names(twin.rec.info) = c("NTR","QTR370","QTR610")


save(twin.rec.info,file="recombination_counts_for_twins_and_siblings.RData")

all.twin.info=rbind(twin.rec.info[[1]][,colnames(twin.rec.info[[1]]) %in% colnames(twin.rec.info[[2]])],twin.rec.info[[2]],twin.rec.info[[3]])

all.twin.info.inf.nuc.fam = all.twin.info[all.twin.info$inf.nuc.fam.mat,]

binom.test(sum(all.twin.info.inf.nuc.fam$twin1.duohmm.mat.nrec > all.twin.info.inf.nuc.fam$sib1.duohmm.mat.nrec,na.rm=T),nrow(all.twin.info.inf.nuc.fam[!is.na(all.twin.info.inf.nuc.fam$twin1.duohmm.mat.nrec) & !is.na(all.twin.info.inf.nuc.fam$sib1.duohmm.mat.nrec),]),p=0.5)

binom.test(sum(all.twin.info.inf.nuc.fam$twin2.duohmm.mat.nrec > all.twin.info.inf.nuc.fam$sib1.duohmm.mat.nrec,na.rm=T),nrow(all.twin.info.inf.nuc.fam[!is.na(all.twin.info.inf.nuc.fam$twin2.duohmm.mat.nrec) & !is.na(all.twin.info.inf.nuc.fam$sib1.duohmm.mat.nrec),]),
           p=0.5)



binom.test(sum(all.twin.info$twin1.duohmm.mat.nrec > all.twin.info$sib1.duohmm.mat.nrec,na.rm=T),nrow(all.twin.info[!is.na(all.twin.info$twin1.duohmm.mat.nrec) & !is.na(all.twin.info$sib1.duohmm.mat.nrec),]),p=0.5)

binom.test(sum(all.twin.info$twin2.duohmm.mat.nrec > all.twin.info$sib1.duohmm.mat.nrec,na.rm=T),nrow(all.twin.info[!is.na(all.twin.info$twin2.duohmm.mat.nrec) & !is.na(all.twin.info$sib1.duohmm.mat.nrec),]), p=0.5)

#all.twin.info=all.twin.info[all.twin.info$inf.nuc.fam.mat & all.twin.info$inf.nuc.fam.pat,]
#all.twin.info=all.twin.info[all.twin.info$inf.mat & all.twin.info$inf.pat,]

all.twin.info$twin.duohmm.mat.nrec.diff = all.twin.info$twin1.duohmm.mat.nrec- all.twin.info$twin2.duohmm.mat.nrec
all.twin.info$twin1.sib.duohmm.mat.nrec.diff = all.twin.info$twin1.duohmm.mat.nrec- all.twin.info$sib1.duohmm.mat.nrec
all.twin.info$twin2.sib.duohmm.mat.nrec.diff = all.twin.info$twin2.duohmm.mat.nrec- all.twin.info$sib1.duohmm.mat.nrec

all.twin.info$twin.duohmm.pat.nrec.diff = all.twin.info$twin1.duohmm.pat.nrec- all.twin.info$twin2.duohmm.pat.nrec
all.twin.info$twin1.sib.duohmm.pat.nrec.diff = all.twin.info$twin1.duohmm.pat.nrec- all.twin.info$sib1.duohmm.pat.nrec
all.twin.info$twin2.sib.duohmm.pat.nrec.diff = all.twin.info$twin2.duohmm.pat.nrec- all.twin.info$sib1.duohmm.pat.nrec

#95% CI for tausq from Model 1
#(55.6595,63.9898) maternal
#(16.0409,18.5702) paternal



pdf("/home/hilary/maternal_age_recombination/interwin_correlation/histograms_of_difference_in_n_crossovers_between_twins_vs_sibs.duoHMM.pdf",height=10,width=15)
par(mfrow=c(2,3))

hist(all.twin.info$twin.duohmm.mat.nrec.diff,main="Maternal, twins",xlab=("R1-R2"),breaks=20,prob=T)
legend("topleft",c(paste0("mean = ",round(mean(all.twin.info$twin.duohmm.mat.nrec.diff,na.rm=T),2)),paste0("var = ",round(var(all.twin.info$twin.duohmm.mat.nrec.diff,na.rm=T),2))))
curve(dnorm(x,0,sqrt(2*55.6595)),add=T)
curve(dnorm(x,0,sqrt(2*63.9898)),add=T,lty=2)
curve(dnorm(x,mean(all.twin.info$twin.duohmm.mat.nrec.diff,na.rm=T),sd(all.twin.info$twin.duohmm.mat.nrec.diff,na.rm=T)),add=T,lty=1,col="red")
legend("topright",c("N(0,2*tausq_2.5th)","N(0,2*tausq_97.5th)","N(emp_mean,emp_var)"),lty=c(1,2,1),col=c("black","black","red"))

hist(all.twin.info$twin1.sib.duohmm.mat.nrec.diff,main="Maternal, twin1-sib",xlab=("R1-R2"),breaks=20,prob=T)
legend("topleft",c(paste0("mean = ",round(mean(all.twin.info$twin1.sib.duohmm.mat.nrec.diff,na.rm=T),2)),paste0("var = ",round(var(all.twin.info$twin1.sib.duohmm.mat.nrec.diff,na.rm=T),2))))
curve(dnorm(x,0,sqrt(2*55.6595)),add=T)
curve(dnorm(x,0,sqrt(2*63.9898)),add=T,lty=2)
curve(dnorm(x,mean(all.twin.info$twin1.sib.duohmm.mat.nrec.diff,na.rm=T),sd(all.twin.info$twin1.sib.duohmm.mat.nrec.diff,na.rm=T)),add=T,lty=1,col="red")

hist(all.twin.info$twin2.sib.duohmm.mat.nrec.diff,main="Maternal, twin2-sib",xlab=("R1-R2"),breaks=20,prob=T)
legend("topleft",c(paste0("mean = ",round(mean(all.twin.info$twin2.sib.duohmm.mat.nrec.diff,na.rm=T),2)),paste0("var = ",round(var(all.twin.info$twin2.sib.duohmm.mat.nrec.diff,na.rm=T),2))))
curve(dnorm(x,0,sqrt(2*55.6595)),add=T)
curve(dnorm(x,0,sqrt(2*63.9898)),add=T,lty=2)
curve(dnorm(x,mean(all.twin.info$twin2.sib.duohmm.mat.nrec.diff,na.rm=T),sd(all.twin.info$twin2.sib.duohmm.mat.nrec.diff,na.rm=T)),add=T,lty=1,col="red")

hist(all.twin.info$twin.duohmm.pat.nrec.diff,main="Paternal, twins",xlab=("R1-R2"),breaks=20,prob=T)
legend("topleft",c(paste0("mean = ",round(mean(all.twin.info$twin.duohmm.pat.nrec.diff,na.rm=T),2)),paste0("var = ",round(var(all.twin.info$twin.duohmm.pat.nrec.diff,na.rm=T),2))))
curve(dnorm(x,0,sqrt(2*16.0409)),add=T)
curve(dnorm(x,0,sqrt(2*18.5702)),add=T,lty=2)
curve(dnorm(x,mean(all.twin.info$twin.duohmm.pat.nrec.diff,na.rm=T),sd(all.twin.info$twin.duohmm.pat.nrec.diff,na.rm=T)),add=T,lty=1,col="red")

hist(all.twin.info$twin1.sib.duohmm.pat.nrec.diff,main="Paternal, twin1-sib",xlab=("R1-R2"),breaks=20,prob=T)
legend("topleft",c(paste0("mean = ",round(mean(all.twin.info$twin1.sib.duohmm.pat.nrec.diff,na.rm=T),2)),paste0("var = ",round(var(all.twin.info$twin1.sib.duohmm.pat.nrec.diff,na.rm=T),2))))
curve(dnorm(x,0,sqrt(2*16.0409)),add=T)
curve(dnorm(x,0,sqrt(2*18.5702)),add=T,lty=2)
curve(dnorm(x,mean(all.twin.info$twin1.sib.duohmm.pat.nrec.diff,na.rm=T),sd(all.twin.info$twin1.sib.duohmm.pat.nrec.diff,na.rm=T)),add=T,lty=1,col="red")

hist(all.twin.info$twin2.sib.duohmm.pat.nrec.diff,main="Paternal, twin2-sib",xlab=("R1-R2"),breaks=20,prob=T)
legend("topleft",c(paste0("mean = ",round(mean(all.twin.info$twin2.sib.duohmm.pat.nrec.diff,na.rm=T),2)),paste0("var = ",round(var(all.twin.info$twin2.sib.duohmm.pat.nrec.diff,na.rm=T),2))))
curve(dnorm(x,0,sqrt(2*16.0409)),add=T)
curve(dnorm(x,0,sqrt(2*18.5702)),add=T,lty=2)
curve(dnorm(x,mean(all.twin.info$twin2.sib.duohmm.pat.nrec.diff,na.rm=T),sd(all.twin.info$twin2.sib.duohmm.pat.nrec.diff,na.rm=T)),add=T,lty=1,col="red")

dev.off()

pdf("/home/hilary/maternal_age_recombination/interwin_correlation/ecdf_of_difference_in_n_crossovers_between_twins_vs_sibs.duoHMM.pdf",height=7,width=12)
par(mfrow=c(2,3))

test1 = ks.test(all.twin.info$twin.duohmm.mat.nrec.diff,pnorm,sd=sqrt(2*55.6595))
test2 = ks.test(all.twin.info$twin.duohmm.mat.nrec.diff,pnorm,sd=sqrt(2*63.9898))
plot(ecdf(all.twin.info$twin.duohmm.mat.nrec.diff),main="Maternal, twins",xlab=("R1-R2"),verticals=T,pch=NA)
lines(ecdf(rnorm(10000,0,sqrt(2*55.6595))),verticals=T,col="blue")
lines(ecdf(rnorm(10000,0,sqrt(2*63.9898))),verticals=T,col="red")
legend("topleft",c("observed data","N(0,2*tausq_2.5th)","N(0,2*tausq_97.5th)"),lty=c(1,2,1),col=c("black","blue","red"))
legend("left",c(paste0("KS test p-value = ",round(test1$p.value,3)),paste0("KS test p-value = ",round(test2$p.value,3))),text.col=c("blue","red"))

test1 = ks.test(all.twin.info$twin1.sib.duohmm.mat.nrec.diff,pnorm,sd=sqrt(2*55.6595))
test2 = ks.test(all.twin.info$twin1.sib.duohmm.mat.nrec.diff,pnorm,sd=sqrt(2*63.9898))
plot(ecdf(all.twin.info$twin1.sib.duohmm.mat.nrec.diff),main="Maternal, twin1-sib",xlab=("R1-R2"),verticals=T,pch=NA)
lines(ecdf(rnorm(10000,0,sqrt(2*55.6595))),verticals=T,col="blue")
lines(ecdf(rnorm(10000,0,sqrt(2*63.9898))),verticals=T,col="red")
legend("left",c(paste0("KS test p-value = ",round(test1$p.value,3)),paste0("KS test p-value = ",round(test2$p.value,3))),text.col=c("blue","red"))

test1 = ks.test(all.twin.info$twin2.sib.duohmm.mat.nrec.diff,pnorm,sd=sqrt(2*55.6595))
test2 = ks.test(all.twin.info$twin2.sib.duohmm.mat.nrec.diff,pnorm,sd=sqrt(2*63.9898))
plot(ecdf(all.twin.info$twin2.sib.duohmm.mat.nrec.diff),main="Maternal, twin2-sib",xlab=("R1-R2"),verticals=T,pch=NA)
lines(ecdf(rnorm(10000,0,sqrt(2*55.6595))),verticals=T,col="blue")
lines(ecdf(rnorm(10000,0,sqrt(2*63.9898))),verticals=T,col="red")
legend("left",c(paste0("KS test p-value = ",round(test1$p.value,3)),paste0("KS test p-value = ",round(test2$p.value,3))),text.col=c("blue","red"))

test1 = ks.test(all.twin.info$twin.duohmm.pat.nrec.diff,pnorm,sd=sqrt(2*16.0409))
test2 = ks.test(all.twin.info$twin.duohmm.pat.nrec.diff,pnorm,sd=sqrt(2*18.5702))
plot(ecdf(all.twin.info$twin.duohmm.pat.nrec.diff),main="Paternal, twins",xlab=("R1-R2"),verticals=T,pch=NA)
lines(ecdf(rnorm(10000,0,sqrt(2*16.0409))),verticals=T,col="blue")
lines(ecdf(rnorm(10000,0,sqrt(2*18.5702))),verticals=T,col="red")
legend("left",c(paste0("KS test p-value = ",round(test1$p.value,3)),paste0("KS test p-value = ",round(test2$p.value,3))),text.col=c("blue","red"))

test1 = ks.test(all.twin.info$twin1.sib.duohmm.pat.nrec.diff,pnorm,sd=sqrt(2*16.0409))
test2 = ks.test(all.twin.info$twin1.sib.duohmm.pat.nrec.diff,pnorm,sd=sqrt(2*18.5702))
plot(ecdf(all.twin.info$twin1.sib.duohmm.pat.nrec.diff),main="Paternal, twin1-sib",xlab=("R1-R2"),verticals=T,pch=NA)
lines(ecdf(rnorm(10000,0,sqrt(2*16.0409))),verticals=T,col="blue")
lines(ecdf(rnorm(10000,0,sqrt(2*18.5702))),verticals=T,col="red")
legend("left",c(paste0("KS test p-value = ",round(test1$p.value,3)),paste0("KS test p-value = ",round(test2$p.value,3))),text.col=c("blue","red"))

test1 = ks.test(all.twin.info$twin2.sib.duohmm.pat.nrec.diff,pnorm,sd=sqrt(2*16.0409))
test2 = ks.test(all.twin.info$twin2.sib.duohmm.pat.nrec.diff,pnorm,sd=sqrt(2*18.5702))
plot(ecdf(all.twin.info$twin2.sib.duohmm.pat.nrec.diff),main="Paternal, twin2-sib",xlab=("R1-R2"),verticals=T,pch=NA)
lines(ecdf(rnorm(10000,0,sqrt(2*16.0409))),verticals=T,col="blue")
lines(ecdf(rnorm(10000,0,sqrt(2*18.5702))),verticals=T,col="red")
legend("left",c(paste0("KS test p-value = ",round(test1$p.value,3)),paste0("KS test p-value = ",round(test2$p.value,3))),text.col=c("blue","red"))

dev.off()


cor.diff.test = function(r1, r2, n1, n2, alternative = c("two.sided", "less", "greater")) {
    Z1 = 0.5 * log( (1+r1)/(1-r1) )
    Z2 = 0.5 * log( (1+r2)/(1-r2) )
    diff = Z1 - Z2
    SEdiff = sqrt( 1 / (n1 - 3) + 1 / (n2 - 3))
    diff.Z = diff / SEdiff
    if (alternative == "less") {
        return(pnorm(diff.Z, lower.tail=F))
    } else if (alternative == "greater") {
        return(pnorm(-diff.Z, lower.tail=F))
    } else if (alternative == "two.sided") {
        return(2 * pnorm( abs(diff.Z), lower.tail=F))
    } else {
        warning(paste("Invalid alterantive", alternative), domain=NA)
        return(NA)
    }
}

calculate.corr.diff.stat = function(data,indiv1,indiv2,indiv3,alternative="two.sided"){
    data=data[,c(indiv1,indiv2,indiv3)]
    data=data[rowSums(is.na(data))==0,]
    n1=nrow(data)
    n2=n1
    r1=cor.test(data[,indiv1],data[,indiv2])
    r2=cor.test(data[,indiv1],data[,indiv3])
    r3 = cor.test(data[,indiv2],data[,indiv3])
    Z1 = 0.5 * log((1+r1$estimate)/(1-r1$estimate))
    Z2 = 0.5 * log((1+r2$estimate)/(1-r2$estimate))
    Z3 = 0.5 * log((1+r3$estimate)/(1-r3$estimate))
    diff.Z.1 = (Z1-Z2)/sqrt( 1/(n1-3) + 1/(n2-3))
    diff.Z.2 = (Z1-Z3)/sqrt( 1/(n1-3) + 1/(n2-3))
    if (alternative == "less") {
        p=pnorm(diff.Z.1, lower.tail=F)
        p2=pnorm(diff.Z.2, lower.tail=F)
    } else if (alternative == "greater") {
        p=pnorm(-diff.Z.1, lower.tail=F)
        p2=pnorm(-diff.Z.2, lower.tail=F)
    } else if (alternative == "two.sided") {
        p=2 * pnorm( abs(diff.Z.1), lower.tail=F)
        p2=2 * pnorm( abs(diff.Z.2), lower.tail=F)
    } else {
        warning(paste("Invalid alterantive", alternative), domain=NA)
        p=NA
        p2=NA
    }
    output=c(n1,r1$estimate,r2$estimate,r3$estimate,r1$conf.int,r2$conf.int,r3$conf.int,r1$p.value,r2$p.value,r3$p.value,Z1,Z1,Z3,diff.Z.1,diff.Z.2,p,p2)
    names(output)=c("N","cor1vs2","cor1vs3","cor2vs3","L.CI.cor1vs2","U.CI.cor1vs2","L.CI.cor1vs3","U.CI.cor1vs3","UL.CI.cor2vs3","U.CI.cor2vs3","p.cor1vs2","p.cor1vs3","p.cor2vs3","Z.cor1vs2","Z.cor1vs3","Z.cor2vs3",
             "diff.Z.1.2.vs.1.3","diff.Z.1.2.vs.2.3",             "p.diff.Z.1.2.vs.1.3","p.diff.Z.1.2.vs.2.3")
return(output)
}


cor.results=list()
#correlations between twins
all.rec0.mat=NULL
all.rec0.pat=NULL
    for(i in 1:length(twin.rec.info)){
    twin.rec.info[[i]]$cohort=names(twin.rec.info)[[i]]        
    rec=twin.rec.info[[i]][,grep("nrec",colnames(twin.rec.info[[i]]))]
    rec0.mat=rec[twin.rec.info[[i]]$inf.nuc.fam.mat,]
    rec0.pat=rec[twin.rec.info[[i]]$inf.nuc.fam.pat,]
    all.rec0.mat=rbind(all.rec0.mat,twin.rec.info[[i]][twin.rec.info[[i]]$inf.nuc.fam.mat,c("cohort","family","twin1","father","mother","twin.mat.age","twin.pat.age","twin2","nsibs","sib1",
        "sib1.mat.age",   "sib1.pat.age","twin1.duohmm.mat.nrec", "twin2.duohmm.mat.nrec","sib1.duohmm.mat.nrec","twin1.duohmm.pat.nrec","twin2.duohmm.pat.nrec","sib1.duohmm.pat.nrec",
        "twin1.nftools.mat.nrec",        "twin2.nftools.mat.nrec","sib1.nftools.mat.nrec","twin1.nftools.pat.nrec","twin2.nftools.pat.nrec","sib1.nftools.pat.nrec" )])
    all.rec0.pat=rbind(all.rec0.pat,twin.rec.info[[i]][twin.rec.info[[i]]$inf.nuc.fam.pat,c("cohort","family","twin1","father","mother","twin.mat.age","twin.pat.age","twin2","nsibs","sib1",
        "sib1.mat.age","sib1.pat.age","twin1.duohmm.mat.nrec","twin2.duohmm.mat.nrec","sib1.duohmm.mat.nrec","twin1.duohmm.pat.nrec","twin2.duohmm.pat.nrec","sib1.duohmm.pat.nrec",
        "twin1.nftools.mat.nrec","twin2.nftools.mat.nrec","sib1.nftools.mat.nrec", "twin1.nftools.pat.nrec","twin2.nftools.pat.nrec","sib1.nftools.pat.nrec" )])
                                        #    rec1.mat=rec[twin.rec.info[[i]]$inf.mat,]
#    rec1.pat=rec[twin.rec.info[[i]]$inf.pat,]
#    rec2.mat=rec[twin.rec.info[[i]]$multikid.mat,]
#    rec2.pat=rec[twin.rec.info[[i]]$multikid.pat,]
   #duohmm inf nuc fams
    cor.twins.0.duohmm.mat = calculate.corr.diff.stat(rec0.mat,"twin1.duohmm.mat.nrec","twin2.duohmm.mat.nrec","sib1.duohmm.mat.nrec")
    cor.twins.0.duohmm.pat = calculate.corr.diff.stat(rec0.pat,"twin1.duohmm.pat.nrec","twin2.duohmm.pat.nrec","sib1.duohmm.pat.nrec")
#nftools inf nuc fams
    cor.twins.0.nftools.mat = calculate.corr.diff.stat(rec0.mat,"twin1.nftools.mat.nrec","twin2.nftools.mat.nrec","sib1.nftools.mat.nrec")
    cor.twins.0.nftools.pat = calculate.corr.diff.stat(rec0.pat,"twin1.nftools.pat.nrec","twin2.nftools.pat.nrec","sib1.nftools.pat.nrec")
#no point including families with < 3 kids since we need a sib for comparison anyway
        #duohmm inf 
                                        #   cor.twins.1.duohmm.mat = calculate.corr.diff.stat(rec1.mat,"twin1.duohmm.mat.nrec","twin2.duohmm.mat.nrec","sib1.duohmm.mat.nrec")
 #   cor.twins.1.duohmm.pat = calculate.corr.diff.stat(rec1.pat,"twin1.duohmm.pat.nrec","twin2.duohmm.pat.nrec","sib1.duohmm.pat.nrec")
#duohmm all families with >1 kid
 #   cor.twins.2.duohmm.mat = calculate.corr.diff.stat(rec2.mat,"twin1.duohmm.mat.nrec","twin2.duohmm.mat.nrec","sib1.duohmm.mat.nrec")
 #   cor.twins.2.duohmm.pat = calculate.corr.diff.stat(rec2.pat,"twin1.duohmm.pat.nrec","twin2.duohmm.pat.nrec","sib1.duohmm.pat.nrec")
 #   cor.results[[i]]=rbind(cor.twins.0.nftools.mat,cor.twins.0.duohmm.mat, cor.twins.1.duohmm.mat,cor.twins.2.duohmm.mat ,cor.twins.0.nftools.pat,cor.twins.0.duohmm.pat,cor.twins.1.duohmm.pat,
 #                  cor.twins.2.duohmm.pat)
 #   rownames(cor.results[[i]]) = c("mat.inf.nuc.fams.nftools","mat.inf.nuc.fams.duohmm","mat.info.duohmm","mat.multikid.duohmm",
 #               "pat.inf.nuc.fams.nftools","pat.inf.nuc.fams.duohmm","pat.info.duohmm","pat.multikid.duohmm")
    cor.results[[i]]=rbind(cor.twins.0.nftools.mat,cor.twins.0.duohmm.mat, cor.twins.0.nftools.pat,cor.twins.0.duohmm.pat)
    rownames(cor.results[[i]]) = c("mat.inf.nuc.fams.nftools","mat.inf.nuc.fams.duohmm",                "pat.inf.nuc.fams.nftools","pat.inf.nuc.fams.duohmm")
}

all.rec0.mat = all.rec0.mat[rowSums(is.na(all.rec0.mat[,c("twin1.duohmm.mat.nrec","twin2.duohmm.mat.nrec","sib1.duohmm.mat.nrec","twin1.duohmm.pat.nrec","twin2.duohmm.pat.nrec","sib1.duohmm.pat.nrec",
    "twin1.nftools.mat.nrec","twin2.nftools.mat.nrec","sib1.nftools.mat.nrec","twin1.nftools.pat.nrec","twin2.nftools.pat.nrec","sib1.nftools.pat.nrec","twin.mat.age","sib1.pat.age","twin.pat.age",
    "sib1.pat.age")]))==0,]
all.rec0.pat=all.rec0.mat
cor.twins.0.duohmm.mat = calculate.corr.diff.stat(all.rec0.mat,"twin1.duohmm.mat.nrec","twin2.duohmm.mat.nrec","sib1.duohmm.mat.nrec")
cor.twins.0.duohmm.pat = calculate.corr.diff.stat(all.rec0.pat,"twin1.duohmm.pat.nrec","twin2.duohmm.pat.nrec","sib1.duohmm.pat.nrec")

cor.twins.0.nftools.mat = calculate.corr.diff.stat(all.rec0.mat,"twin1.nftools.mat.nrec","twin2.nftools.mat.nrec","sib1.nftools.mat.nrec")
cor.twins.0.nftools.pat = calculate.corr.diff.stat(all.rec0.pat,"twin1.nftools.pat.nrec","twin2.nftools.pat.nrec","sib1.nftools.pat.nrec")

cor.results[[4]]=rbind(cor.twins.0.nftools.mat,cor.twins.0.duohmm.mat, cor.twins.0.nftools.pat,cor.twins.0.duohmm.pat)
rownames(cor.results[[4]]) = c("mat.inf.nuc.fams.nftools","mat.inf.nuc.fams.duohmm",                "pat.inf.nuc.fams.nftools","pat.inf.nuc.fams.duohmm")
if(FALSE){
write.table(cor.results[[4]],"/home/hilary/maternal_age_recombination/interwin_correlation/correlations_between_twins_and_twin_sib_pair.informative_nuclear_families.more_stringent.txt",quote=F,sep="\t")

#Simulate same number of families under a model of no age effect, or very small age effect (can the difference in correlations be explained by a small age effect, as we think exists?)
#calculate the correlation for two twins, and between twin1 and sib\\
#Calculate transformed difference in statistics
#calculate this statistic for the actual data for corr(twin1,twin2) vs corr(twin1,sib), then ask what proportion of simulated datasets had more extreme value

#what is the correct way to simualte? using posterior means? posterior means for all parameters except beta age? best draw (although parameters may not be within 95% CI? draw from posterior differentlyfor each simulation? but should presumably use same parameters for each

    simulate.counts = function(all.rec0.mat,mu.post.mean,sigmasq.post.mean,beta_Age,twin.age,sib1.age,tausq){
all.sim.data=NULL
    for(j in 1:1000){
#        if(j %/% 100 ==0){
#            print(j)
#        }
        a0=sapply(1:nrow(all.rec0.mat),function(i){rnorm(1,mean=mu.post.mean[all.rec0.mat[i,"cohort"],"mean"],sd=sqrt(sigmasq.post.mean[all.rec0.mat[i,"cohort"],"mean"]))})
        twin1.fake.y=c()
        twin2.fake.y=c()
        sib1.fake.y=c()
        for(i in 1:length(a0)){
            twin1.fake.y=c(twin1.fake.y,rnorm(1,a0[i]+beta_Age * twin.age[i],sd=sqrt(tausq)))
            twin2.fake.y=c(twin2.fake.y,rnorm(1,a0[i]+beta_Age * twin.age[i],sd=sqrt(tausq)))
            sib1.fake.y=c(sib1.fake.y,rnorm(1,a0[i]+beta_Age * sib1.age[i],sd=sqrt(tausq)))
        }
        all.fake.y=cbind(twin1.fake.y,twin2.fake.y,sib1.fake.y)
        colnames(all.fake.y)=c("twin1","twin2","sib1")
        cor.sim.data = calculate.corr.diff.stat(all.fake.y,"twin1","twin2","sib1")
        all.sim.data=rbind(all.sim.data,cor.sim.data)
    }
return(all.sim.data)
}

    
#description="model1.6.2.maternal.mu_m_N_41_100.sigmasq_m_IG_2_40"
#description="model1.6.2.paternal.mu_m_N_27_64.sigmasq_m_IG_2_15"
#mydir="RSTAN_output_on_duoHMM_more_stringent/model1.62"
#load(best.draw,file=paste0(mydir,"/",description,".best_draw.RData"))
#mat.posterior.summary=read.delim("RSTAN_output_on_duoHMM_more_stringent/Model1_summary_of_posteriors_maternal.duohmm.txt",header=T)
mat.posterior.summary=read.delim("RSTAN_output_with_NFTOOLS_more_stringent/Model1_summary_of_posteriors_maternal.nftools.txt",header=T)
#pat.posterior.summary=read.delim("RSTAN_output_on_duoHMM_more_stringent/Model1_summary_of_posteriors_paternal.duohmm.txt",header=T)
pat.posterior.summary=read.delim("RSTAN_output_with_NFTOOLS_more_stringent/Model1_summary_of_posteriors_paternal.nftools.txt",header=T)
    
beta.age.mat=c( mat.posterior.summary[1,"mean"],0,0.082)

beta.age.pat = c(pat.posterior.summary[1,"mean"],0,0.082)

names(beta.age.mat) = c("posterior.mean","no","decode")

for(b in 1:3){
    print(b)
#maternal simulations
sigmasq.post.mean=mat.posterior.summary[grep("sigmasq_m",mat.posterior.summary[,3]),]
mu.post.mean=mat.posterior.summary[grep("mu_m",mat.posterior.summary[,3]),]
rownames(sigmasq.post.mean)  = sigmasq.post.mean$Cohort
rownames(mu.post.mean)  = mu.post.mean$Cohort

mother = sort(rep(1:nrow(all.rec0.mat),3))
twin.age=all.rec0.mat$twin.mat.age
sib1.age=all.rec0.mat$sib1.mat.age

beta_Age = beta.age.mat[b]
tausq = mat.posterior.summary[2,"mean"]

mat.sims=  simulate.counts(all.rec0.mat,mu.post.mean,sigmasq.post.mean,beta_Age,twin.age,sib1.age,tausq)

#paternal
sigmasq.post.mean=pat.posterior.summary[grep("sigmasq_m",pat.posterior.summary[,3]),]
mu.post.mean=pat.posterior.summary[grep("mu_m",pat.posterior.summary[,3]),]
rownames(sigmasq.post.mean)  = sigmasq.post.mean$Cohort
rownames(mu.post.mean)  = mu.post.mean$Cohort

twin.age=all.rec0.mat$twin.pat.age
sib1.age=all.rec0.mat$sib1.pat.age

beta_Age = beta.age.pat[b]
tausq = pat.posterior.summary[2,"mean"]

pat.sims=  simulate.counts(all.rec0.mat,mu.post.mean,sigmasq.post.mean,beta_Age,twin.age,sib1.age,tausq)

#compare to actual results
true.results = cor.results[[4]]

emp.pvals.mat.nftools=sapply(1:ncol(true.results),function(i){sum(mat.sims[,i] > true.results[1,i])})/nrow(mat.sims)
emp.pvals.mat.duohmm=sapply(1:ncol(true.results),function(i){sum(mat.sims[,i] > true.results[2,i])})/nrow(mat.sims)
emp.pvals.pat.nftools=sapply(1:ncol(true.results),function(i){sum(pat.sims[,i] > true.results[3,i])})/nrow(pat.sims)
emp.pvals.pat.duohmm=sapply(1:ncol(true.results),function(i){sum(pat.sims[,i] > true.results[4,i])})/nrow(pat.sims)

emp.pvals = rbind(emp.pvals.mat.nftools,emp.pvals.mat.duohmm,emp.pvals.pat.nftools,emp.pvals.pat.duohmm)
colnames(emp.pvals) = colnames(true.results)
rownames(emp.pvals) = c("mat.nftools","mat.duohmm","pat.nftools","pat.duohmm")
emp.pvals=emp.pvals[,c("Z.cor1vs2","Z.cor1vs3","Z.cor2vs3","diff.Z.1.2.vs.1.3","diff.Z.1.2.vs.2.3" )]

#write.table(emp.pvals,paste0("/home/hilary/maternal_age_recombination/interwin_correlation/empirical_pvalues_for_difference_in_transformed_correlations_between_twins_and_twin_sib_pair.simulated_with_",names(beta.age.mat)[b],".beta.Age.txt"),quote=F,sep="\t")
    write.table(emp.pvals,paste0("/home/hilary/maternal_age_recombination/interwin_correlation/empirical_pvalues_for_difference_in_transformed_correlations_between_twins_and_twin_sib_pair.simulated_with_",names(beta.age.mat)[b],".beta.Age.nftools.txt"),quote=F,sep="\t")


}
    
 }   



#### theoretical calculations


beta.ages = c(0,0.05,0.1,0.5,1,5)

twin.ages=all.rec0.mat$twin.mat.age
sib.ages=all.rec0.mat$sib1.mat.age

mat.posterior.summary=read.delim("RSTAN_output_with_NFTOOLS_more_stringent/Model1_summary_of_posteriors_maternal.nftools.txt",header=T)
sigmasq.post.mean=mat.posterior.summary[grep("sigmasq_m",mat.posterior.summary[,3]),]
mu.post.mean=mat.posterior.summary[grep("mu_m",mat.posterior.summary[,3]),]
rownames(sigmasq.post.mean)  = sigmasq.post.mean$Cohort
rownames(mu.post.mean)  = mu.post.mean$Cohort

sigmasq_z=66 #tausq
sigmasq_alpha=8 #pick same one for all cohorts

theoretical.values.mat= NULL
for(b in 1:length(beta.ages)){
  beta_age=beta.ages[b]
  print(beta_age)
  theoretical.cor.twins = (sigmasq_alpha + beta_age^2 *var(twin.ages))/(sigmasq_alpha + beta_age^2 *var(twin.ages)+ sigmasq_z)
    theoretical.cor.sibs = (sigmasq_alpha + beta_age^2 *cov(twin.ages,sib.ages))/(sqrt((sigmasq_alpha + beta_age^2*var(twin.ages)+ sigmasq_z) * (sigmasq_alpha + beta_age^2*var(sib.ages)+ sigmasq_z)))
theoretical.values.mat=rbind(theoretical.values.mat,c(beta_age,sigmasq_alpha,sigmasq_z,var(twin.ages),var(sib.ages),cov(twin.ages,sib.ages),theoretical.cor.twins,theoretical.cor.sibs))
}
colnames(theoretical.values.mat) = c("beta.age","sigmasq_alpha","sigmasq_z","var_twin_ages","var_sib_ages","cov_twin_sib_ages","E_r_twins","E_r_sibs")

write.table(theoretical.values.mat,"/home/hilary/maternal_age_recombination/interwin_correlation/expected_correlations_for_realistic_values_of_parameters.maternal.txt",quote=F,sep="\t")
                                        #pat.posterior.summary=read.delim("RSTAN_output_on_duoHMM_more_stringent/Model1_summary_of_posteriors_paternal.duohmm.txt",header=T)

twin.ages=all.rec0.mat$twin.pat.age
sib.ages=all.rec0.mat$sib1.pat.age

mat.posterior.summary=read.delim("RSTAN_output_with_NFTOOLS_more_stringent/Model1_summary_of_posteriors_paternal.nftools.txt",header=T)

sigmasq.post.mean=mat.posterior.summary[grep("sigmasq_m",mat.posterior.summary[,3]),]
mu.post.mean=mat.posterior.summary[grep("mu_m",mat.posterior.summary[,3]),]
rownames(sigmasq.post.mean)  = sigmasq.post.mean$Cohort
rownames(mu.post.mean)  = mu.post.mean$Cohort

sigmasq_z= 18 #tausq
sigmasq_alpha=2 #pick same one for all cohorts


theoretical.values.pat= NULL
for(beta_age in beta.ages){
  theoretical.cor.twins = (sigmasq_alpha + beta_age^2 *var(twin.ages))/(sigmasq_alpha + beta_age^2 *var(twin.ages)+ sigmasq_z)
    theoretical.cor.sibs = (sigmasq_alpha + beta_age^2 *cov(twin.ages,sib.ages))/(sqrt((sigmasq_alpha + beta_age^2*var(twin.ages)+ sigmasq_z) * (sigmasq_alpha + beta_age^2*var(sib.ages)+ sigmasq_z)))
theoretical.values.pat=rbind(theoretical.values.pat,c(beta_age,sigmasq_alpha,sigmasq_z,var(twin.ages),var(sib.ages),cov(twin.ages,sib.ages),theoretical.cor.twins,theoretical.cor.sibs))
}
colnames(theoretical.values.pat) = c("beta.age","sigmasq_alpha","sigmasq_z","var_twin_ages","var_sib_ages","cov_twin_sib_ages","E_r_twins","E_r_sibs")

write.table(theoretical.values.pat,"/home/hilary/maternal_age_recombination/interwin_correlation/expected_correlations_for_realistic_values_of_parameters.paternal.txt",quote=F,sep="\t")
