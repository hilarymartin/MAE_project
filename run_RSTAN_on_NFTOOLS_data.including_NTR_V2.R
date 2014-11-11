library(rstan)
library(scales)
library(parallel)
load("NFTOOLS_data_for_RSTAN.NTR_v2.RData")
argv <- commandArgs(trailingOnly = TRUE)
print(argv)    
mycols=c("black","blue","red","green","orange","red4","purple","darkgreen")
names(mycols)=c("CARL","FC","FVG","GPC","NTR","QTR370","QTR610","VB")


if(argv[1]==1){
### Model 1.
#for informative meioses only with nkid=>2, all cohorts, fit cohort effect and family effect
#cohort-specific age effect

    model1.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.stan",data=data1.mat2,iter=10,chains=0,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0"))
#    model1.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0")))

    if(FALSE){
#MAP
   beta_Age = c(0.002663752,    -0.12945442,    0.092491647,    0.015474999,     -0.004568869,    0.020810457,    0.150696835,    0.210443124)
    beta_Cohort = c(4.415177335,    -0.570331424,    1.06200512,    3.047638636,    1.179672138,    -1.753863924,    -4.460305374)
    tausq = 68.41491232
    sigmasq_m = 0.028946522
    mu_m = 39.30664645
}
 #   cat("Running model 1, initializing estimates close to MAP estimates\n")
       cat("Running model 1, initializing estimates at MAP estimates\n")
#    model1.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0"),
#                                         init=list(list(sigmasq_m=0.05,beta_Age=rep(0,8),beta_Cohort=rep(0,7),tausq=19,mu_m=31,a0=rnorm(n=data1.mat2$I,mean=31,sd=sqrt(0.05))))))
    model1.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0"),
                                         init=list(list(sigmasq_m=0.028946522,beta_Age=    c(0.002663752,    -0.12945442,    0.092491647,    0.015474999,     -0.004568869,    0.020810457,    0.150696835,    0.210443124),
                                             beta_Cohort=   c(4.415177335,    -0.570331424,    1.06200512,    3.047638636,    1.179672138,    -1.753863924,    -4.460305374) ,tausq=68.41491232,mu_m=39.30664645,
                                             a0=rnorm(n=data1.mat2$I,mean=    39.30664645,sd=sqrt(  0.028946522))))))

    model1.mat <- sflist2stanfit(model1.mat.list)

                                print(model1.mat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
                                        #    save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model1.mat.initialized_close_to_MAP.RData")
                                save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model1.mat.initialized_at_MAP.RData")
                                        #    save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model1.mat.RData")
                                        #load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.mat.RData")
if(FALSE){
    mysim<-extract(model1.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    my.ylim=7
    pdf("RSTAN_output_with_NFTOOLS/model1.maternal.posteriors.beta_Age.pdf",height=5,width=5)
    for(i in 1:ncol(mysim$beta_Age)){
        if(i==1){
            plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
        abline(v=0,lwd=2)
        polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
        abline(v=0.067,lty=2,lwd=2)
        polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
        abline(v=0.082,lty=3,lwd=2)
        polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
        abline(v=0.19,lty=4,lwd=2)
        abline(v=-0.42,lty=4,lwd=2,col="grey")
        legend("topleft",as.character(cohort.codes[order(cohort.codes[,1]),2]),col=mycols[as.character(cohort.codes[,2])],lty=1,lwd=2,cex=0.5)
        legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
    }
    dev.off()
}
}else if(argv[1] ==2){

    model1.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.stan",data=data1.pat2,iter=10,chains=0,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0"))
    model1.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0")))
    model1.pat <- sflist2stanfit(model1.pat.list)
    print(model1.pat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
    save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model1.pat.RData")
                                        #load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.pat.RData")
    mysim<-extract(model1.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    my.ylim=15
    pdf("RSTAN_output_with_NFTOOLS/model1.paternal.posteriors.beta_Age.pdf",height=5,width=5)
    for(i in 1:ncol(mysim$beta_Age)){
        if(i==1){
            plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
        legend("topleft",as.character(cohort.codes[order(cohort.codes[,1]),2]),col=mycols[as.character(cohort.codes[,2])],lty=1,lwd=2,cex=0.5)
        abline(v=0,lwd=2)
    }
    dev.off()
    
}else if(argv[1] ==3){
#common age effect
    model1.2.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.2.stan",data=data1.mat2,iter=10,chains=0,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0"))
    model1.2.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.2.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0")))
    model1.2.mat <- sflist2stanfit(model1.2.mat.list)
    
    print(model1.2.mat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
    save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model1.2.mat.RData")
                                        #load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.2.mat.RData")
    mysim<-extract(model1.2.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    my.ylim=12
    pdf("RSTAN_output_with_NFTOOLS/model1.2.maternal.posteriors.beta_Age.pdf",height=5,width=5)
    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
    polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
    abline(v=0.067,lty=2,lwd=2)
    polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
    abline(v=0.082,lty=3,lwd=2)
    polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
    abline(v=0.19,lty=4,lwd=2)
    abline(v=-0.42,lty=4,lwd=2,col="grey")
    legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
    dev.off()
}else if(argv[1] ==4){
    model1.2.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.2.stan",data=data1.pat2,iter=10,chains=0,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0"))
    model1.2.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.2.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0")))
    model1.2.pat <- sflist2stanfit(model1.2.pat.list)
    
    print(model1.2.pat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
    save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model1.2.pat.RData")
                                        #load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.2.pat.RData")
    mysim<-extract(model1.2.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    my.ylim=12
    pdf("RSTAN_output_with_NFTOOLS/model1.2.paternal.posteriors.beta_Age.pdf",height=5,width=5)
    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col="black",lwd=2)
    abline(v=0,lwd=2)
    dev.off()
    
}else if(argv[1] == 9){
#cohort-specific age effects, but drawn from N(beta_global,sigmasq_global); with cohort effect
    model1.3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.3.stan",data=data1.mat2,iter=10,chains=0,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0","beta_global","sigmasq_global"))
    model1.3.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.3.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0","beta_global","sigmasq_global")))
    model1.3.mat <- sflist2stanfit(model1.3.mat.list)

    print(model1.3.mat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
    save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model1.3.mat.RData")
                                        #load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.3.mat.RData")
    print(model1.3.mat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
    mysim<-extract(model1.3.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    pdf("RSTAN_output_with_NFTOOLS/model1.3.maternal.posteriors.beta_Age.pdf",height=5,width=5)
    my.ylim=7
    for(i in 1:ncol(mysim$beta_Age)){
        if(i==1){
            plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
        lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
        abline(v=0,lwd=2)
        polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
        abline(v=0.067,lty=2,lwd=2)
        polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
        abline(v=0.082,lty=3,lwd=2)
        polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
        abline(v=0.19,lty=4,lwd=2)
        abline(v=-0.42,lty=4,lwd=2,col="grey")
        legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
        legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
    }
    dev.off()
    pdf("RSTAN_output_with_NFTOOLS/model1.3.maternal.posteriors.sigmasq_global.pdf",height=5,width=5)
    plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
    dev.off()
}else if(argv[1] ==10){

    model1.3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.3.stan",data=data1.pat2,iter=10,chains=0,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0","beta_global","sigmasq_global"))
    model1.3.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.3.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","a0","beta_global","sigmasq_global")))
    model1.3.pat <- sflist2stanfit(model1.3.pat.list)

    print(model1.3.pat,pars=c("beta_Age","beta_Cohort","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
    save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model1.3.pat.RData")
    
    mysim<-extract(model1.3.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    my.ylim=14
    pdf("RSTAN_output_with_NFTOOLS/model1.3.paternal.posteriors.beta_Age.pdf",height=5,width=5)
    for(i in 1:ncol(mysim$beta_Age)){
        if(i==1){
            plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
        lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
        abline(v=0,lwd=2)
        legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    }
    dev.off()
    pdf("RSTAN_output_with_NFTOOLS/model1.3.paternal.posteriors.sigmasq_global.pdf",height=5,width=5)
    plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
    dev.off()
    
    
}else if(argv[1] == 11){
#cohort-specific age effects, but drawn from N(beta_global,sigmasq_global); no cohort effect

    model1.4.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.4.stan",data=data1.mat2,iter=10,chains=0,pars=c("beta_Age","tausq","sigmasq_m","mu_m","a0","beta_global","sigmasq_global"))
    model1.4.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.4.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","sigmasq_m","mu_m","a0","beta_global","sigmasq_global")))
    model1.4.mat <- sflist2stanfit(model1.4.mat.list)

    print(model1.4.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
    save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model1.4.mat.RData")
                                        #load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.4.mat.RData")
    
    mysim<-extract(model1.4.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    my.ylim=12
    pdf("RSTAN_output_with_NFTOOLS/model1.4.maternal.posteriors.beta_Age.pdf",height=5,width=5)
    
    for(i in 1:ncol(mysim$beta_Age)){
        if(i==1){
            plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
        lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
        abline(v=0,lwd=2)
        polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
        abline(v=0.067,lty=2,lwd=2)
        polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
        abline(v=0.082,lty=3,lwd=2)
        polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
        abline(v=0.19,lty=4,lwd=2)
        abline(v=-0.42,lty=4,lwd=2,col="grey")
        legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
        legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
    }
    dev.off()
    pdf("RSTAN_output_with_NFTOOLS/model1.4.maternal.posteriors.sigmasq_global.pdf",height=5,width=5)
    plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
    dev.off()
    
}else if(argv[1] ==12){
    model1.4.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.4.stan",data=data1.pat2,iter=10,chains=0,pars=c("beta_Age","tausq","sigmasq_m","mu_m","a0","beta_global","sigmasq_global"))
    model1.4.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.4.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","sigmasq_m","mu_m","a0","beta_global","sigmasq_global")))
    model1.4.pat <- sflist2stanfit(model1.4.pat.list)
    print(model1.4.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
    save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model1.4.pat.RData")
    
#load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model1.4.pat.RData")

mysim<-extract(model1.4.pat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
my.ylim=15
pdf("RSTAN_output_with_NFTOOLS/model1.4.paternal.posteriors.beta_Age.pdf",height=5,width=5)

for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
}
dev.off()
pdf("RSTAN_output_with_NFTOOLS/model1.4.paternal.posteriors.sigmasq_global.pdf",height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
    
}else if(argv[1] ==5){
### Model 2.
#for informative meioses only, all cohorts, fit cohort effect, BUT NO FAMILY EFFECT
#cohort-specific age effect

    model2.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model2.stan",data=data1.2.mat2,iter=10,chains=0,pars=c("beta_Age","beta_Cohort","tausq","mu_m"))
    model2.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model2.mat.compiled, data = data1.2.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","mu_m")))
    model2.mat <- sflist2stanfit(model2.mat.list)

print(model2.mat,pars=c("beta_Age","beta_Cohort","tausq","mu_m"),digits=4)
save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model2.mat.RData")

#load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model2.mat.RData")
mysim<-extract(model2.mat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_2.txt",header=T)
my.ylim=8
pdf("RSTAN_output_with_NFTOOLS/model2.maternal.posteriors.beta_Age.pdf",height=5,width=5)

for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
abline(v=0,lwd=2)
polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")
legend("topleft",as.character(cohort.codes[order(cohort.codes[,1]),2]),col=mycols[as.character(cohort.codes[,2])],lty=1,lwd=2,cex=0.5)
legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
}
dev.off()

    
}else if(argv[1] ==6){

    model2.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model2.stan",data=data1.2.pat2,iter=10,chains=0,pars=c("beta_Age","beta_Cohort","tausq","mu_m"))
    model2.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model2.pat.compiled, data = data1.2.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","mu_m")))
    model2.pat <- sflist2stanfit(model2.pat.list)

    print(model2.pat,pars=c("beta_Age","beta_Cohort","tausq","mu_m"),digits=4)
    save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model2.pat.RData")
                                
#load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model2.pat.RData")
mysim<-extract(model2.pat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_2.txt",header=T)
my.ylim=15
pdf("RSTAN_output_with_NFTOOLS/model2.paternal.posteriors.beta_Age.pdf",height=5,width=5)

for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
legend("topleft",as.character(cohort.codes[order(cohort.codes[,1]),2]),col=mycols[as.character(cohort.codes[,2])],lty=1,lwd=2,cex=0.5)
abline(v=0,lwd=2)
}
dev.off()

}else if(argv[1] ==7){
#common age effect
   model2.2.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model2.2.stan",data=data1.2.mat2,iter=10,chains=0,pars=c("beta_Age","beta_Cohort","tausq","mu_m"))
    model2.2.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model2.2.mat.compiled, data = data1.2.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","mu_m")))
    model2.2.mat <- sflist2stanfit(model2.2.mat.list)
    print(model2.2.mat,pars=c("beta_Age","beta_Cohort","tausq","mu_m"),digits=4)
    save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model2.2.mat.RData")

#load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model2.2.mat.RData")
mysim<-extract(model2.2.mat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_2.txt",header=T)
pdf("RSTAN_output_with_NFTOOLS/model2.2.maternal.posteriors.beta_Age.pdf",height=5,width=5)
my.ylim=12

    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
abline(v=0,lwd=2)
polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")

legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
dev.off()

   
}else if(argv[1] ==8){
    model2.2.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model2.2.stan",data=data1.2.pat2,iter=10,chains=0,pars=c("beta_Age","beta_Cohort","tausq","mu_m"))
    model2.2.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model2.2.pat.compiled, data = data1.2.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","tausq","mu_m")))
    model2.2.pat <- sflist2stanfit(model2.2.pat.list)

    print(model2.2.pat,pars=c("beta_Age","beta_Cohort","tausq","mu_m"),digits=4)
save.image("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model2.2.pat.RData")

#load("/well/donnelly/hilary/maternal_age_and_recombination/running_RSTAN.model2.2.pat.RData")
mysim<-extract(model2.2.pat,permuted=T)
cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_2.txt",header=T)
pdf("RSTAN_output_with_NFTOOLS/model2.2.paternal.posteriors.beta_Age.pdf",height=5,width=5)
my.ylim=12

    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col="black",lwd=2)

abline(v=0,lwd=2)

dev.off()

}

#### Model 3.
if(argv[1] ==25){
    cat("Model 3 maternal, with adjusted priors\n")
#    model3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.adjusted_priors.stan",data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
#    model3.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.mat.compiled, data = data2.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
#    model3.mat <- sflist2stanfit(model3.mat.list)
#    save(model3.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.model3.adjusted_priors.all_chains.mat.including_alpha.RData"))
#    print(model3.mat,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.model3.adjusted_priors.all_chains.mat.including_alpha.RData")
    mysim<-extract(model3.mat,permuted=T)
    parent="adjusted_priors.maternal"
    cohort.codes=read.delim("RSTAN_output_with_NFTOOLS/key_for_maternal_cohorts_to_include_in_model_3.NTR_v2.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_maternal.including_NTR_v2.txt",header=T)[,1:5]

}

if(argv[1] ==26){
    cat("Model 3 paternal, with adjusted priors\n")
#    model3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.adjusted_priors.stan",data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
#    model3.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.pat.compiled, data = data2.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
#    model3.pat <- sflist2stanfit(model3.pat.list)
#    save(model3.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.model3.adjusted_priors.all_chains.pat.including_alpha.RData"))
#    print(model3.pat,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.model3.adjusted_priors.all_chains.pat.including_alpha.RData")
    mysim<-extract(model3.pat,permuted=T)
    parent="adjusted_priors.paternal"
    cohort.codes=read.delim("RSTAN_output_with_NFTOOLS/key_for_paternal_cohorts_to_include_in_model_3.NTR_v2.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_paternal.including_NTR_v2.txt",header=T)[,1:5]
}


if(argv[1]==25|argv[1]==26){
    cohort.codes=cohort.codes[order(cohort.codes[,1]),]
    cohort.codes$Pop=unlist(lapply(strsplit(cohort.codes$Cohort,".",fixed=T),function(x){return(x[[1]])}))
    cohort.codes$Fam.type=sapply(1:nrow(cohort.codes),function(x){gsub(paste0(cohort.codes$Pop[x],"."),"",cohort.codes$Cohort[x])})
    beta.parameters=beta.parameters[order(beta.parameters$Cohort.type),]
    mylty=c(1,2,4,1,2,4)
    names(mylty)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
    mytype=c("l","l","l","b","b","b")
    names(mytype)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
    mylwd=c(2,2,2,1,1,1)
    names(mylwd)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
        mypch=c(NA,NA,NA,19,19,19)
    names(mypch)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")

    pdf(paste0("RSTAN_output_with_NFTOOLS/model3.",parent,".posteriors.beta_Age.pdf"),height=5,width=5)
    my.ylim=27
    for(i in 1:ncol(mysim$beta_Age)){
        if(i==1){
            plot(density(mysim$beta_Age[,i]),xlim=c(-0.2,0.2),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],ylim=c(0,my.ylim),
                 lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])])
        }else {
            lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,"Pop"])],lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])])
        }
        lines(density(mysim$beta_global),col="black",lwd=3,lty=3)
          legend("topleft",c(names(mycols),"global","informative, 2 generations","informative, 3 generations","uninformative, 2 kids","both parents","one parent"),col=c(mycols,rep("black",6)),lty=c(rep(1,length(mycols)),3,1,2,4,1,1),
                 lwd=c(rep(2,length(mycols)),3,2,2,2,2,1),cex=0.5)
        abline(v=0,lwd=2)
        if(parent=="adjusted_priors.maternal"){
            polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
            abline(v=0.067,lty=2,lwd=2)
            polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
            abline(v=0.082,lty=3,lwd=2)
            polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
            abline(v=0.19,lty=4,lwd=2)
            abline(v=-0.42,lty=4,lwd=2,col="grey")
            legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
        }
    }
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS/model3.",parent,".posteriors.p.pdf"),height=5,width=5)
    for(i in 1:ncol(mysim$p_by_cohort)){
        mylabels=c(paste(c("informative, 2 generations","informative, 3 generations","uninformative, 2 kids"),"both parents",sep=", "),paste(c("informative, 2 generations","informative, 3 generations","uninformative, 2 kids"),"1 parent",sep=", "))
        names(mylabels)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
        sample.size= sum(data2.mat$cohort.family.type==cohort.codes[cohort.codes[,1]==i,"Cohort"])
        plot(density(mysim$p_by_cohort[,i]),xlab="p",col=mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],xlim=c(0,1),
             lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],main=paste0(cohort.codes[cohort.codes[,1]==i,"Pop"],", ",mylabels[cohort.codes[cohort.codes[,1]==i,"Fam.type"]]))
        curve(dbeta(x,beta.parameters[beta.parameters$Cohort.type==i,"alpha"],beta.parameters[beta.parameters$Cohort.type==i,"beta"]),add=T)
        abline(v=median(mysim$p_by_cohort[,i]),col=mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])])
        abline(v=beta.parameters[beta.parameters$Cohort.type==i,"mean.relative.to.all.informative"],col="black")
        legend("topleft",c("posterior","prior",paste0("n = ",sample.size)),col=c(mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],"black",NA),lty=c(mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],1,NA),
               lwd=c(mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],1,NA))
    }
dev.off()
    pdf(paste0("RSTAN_output_with_NFTOOLS/model3.",parent,".posteriors.sigmasq_global.pdf"),height=5,width=5)
      plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
    dev.off()
}


####model 3.2

if(argv[1] ==27){
cat("Model 3.2 maternal, with adjusted priors\n")
#model3.2.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.stan",data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
#model3.2.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.2.mat.compiled, data = data2.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
#model3.2.mat <- sflist2stanfit(model3.2.mat.list)
#save(model3.2.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.model3.2.adjusted_priors.all_chains.mat.including_alpha.RData"))
#print(model3.2.mat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.model3.2.adjusted_priors.all_chains.mat.including_alpha.RData")
mysim<-extract(model3.2.mat,permuted=T)
parent="adjusted_priors.maternal"
cohort.codes=read.delim("RSTAN_output_with_NFTOOLS/key_for_maternal_cohorts_to_include_in_model_3.NTR_v2.txt",header=T,stringsAsFactors=F)
beta.parameters=read.delim("parameters_for_beta_distribution.model3_maternal.including_NTR_v2.txt",header=T)[,1:5]

}

if(argv[1] ==28){
    cat("Model 3.2 paternal, with adjusted priors\n")
#    model3.2.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.stan",data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
#    model3.2.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.2.pat.compiled, data = data2.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
#    model3.2.pat <- sflist2stanfit(model3.2.pat.list)
#    save(model3.2.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.model3.2.adjusted_priors.all_chains.pat.including_alpha.RData"))
#    print(model3.2.pat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.model3.2.adjusted_priors.all_chains.pat.including_alpha.RData")
    mysim<-extract(model3.2.pat,permuted=T)
    parent="adjusted_priors.paternal"
    cohort.codes=read.delim("RSTAN_output_with_NFTOOLS/key_for_paternal_cohorts_to_include_in_model_3.NTR_v2.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_paternal.including_NTR_v2.txt",header=T)[,1:5]
}


if(argv[1] %in% c(27,28)){
    cohort.codes=cohort.codes[order(cohort.codes[,1]),]
    cohort.codes$Pop=unlist(lapply(strsplit(cohort.codes$Cohort,".",fixed=T),function(x){return(x[[1]])}))
    cohort.codes$Fam.type=sapply(1:nrow(cohort.codes),function(x){gsub(paste0(cohort.codes$Pop[x],"."),"",cohort.codes$Cohort[x])})
    beta.parameters=beta.parameters[order(beta.parameters$Cohort.type),]
    mylty=c(1,2,4,1,2,4)
    names(mylty)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
    mytype=c("l","l","l","b","b","b")
    names(mytype)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
    mylwd=c(2,2,2,1,1,1)
    names(mylwd)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
        mypch=c(NA,NA,NA,19,19,19)
    names(mypch)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")

    pdf(paste0("RSTAN_output_with_NFTOOLS/model3.2.",parent,".posteriors.beta_Age.pdf"),height=5,width=5)
    plot(density(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age")
#    legend("topleft",c(names(mycols),"global","informative, 2 generations","informative, 3 generations","uninformative, 2 kids","both parents","one parent"),col=c(mycols,rep("black",6)),lty=c(rep(1,length(mycols)),3,1,2,4,1,1),
#                 lwd=c(rep(2,length(mycols)),3,2,2,2,2,1),cex=0.5)
    abline(v=0,lwd=2)
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS/model3.2.",parent,".posteriors.p.pdf"),height=5,width=5)
    for(i in 1:ncol(mysim$p_by_cohort)){
        mylabels=c(paste(c("informative, 2 generations","informative, 3 generations","uninformative, 2 kids"),"both parents",sep=", "),paste(c("informative, 2 generations","informative, 3 generations","uninformative, 2 kids"),"1 parent",sep=", "))
        names(mylabels)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
        sample.size= sum(data2.mat$cohort.family.type==cohort.codes[cohort.codes[,1]==i,"Cohort"])
        plot(density(mysim$p_by_cohort[,i]),xlab="p",col=mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],xlim=c(0,1),
             lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],main=paste0(cohort.codes[cohort.codes[,1]==i,"Pop"],", ",mylabels[cohort.codes[cohort.codes[,1]==i,"Fam.type"]]))
        curve(dbeta(x,beta.parameters[beta.parameters$Cohort.type==i,"alpha"],beta.parameters[beta.parameters$Cohort.type==i,"beta"]),add=T)
        abline(v=median(mysim$p_by_cohort[,i]),col=mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])])
        abline(v=beta.parameters[beta.parameters$Cohort.type==i,"mean.relative.to.all.informative"],col="black")
        legend("topleft",c("posterior","prior",paste0("n = ",sample.size)),col=c(mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],"black",NA),lty=c(mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],1,NA),
               lwd=c(mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],1,NA))
    }
    dev.off()
  
 
}

if(argv[1] ==29 | argv[1]==30){
#### Model 1.6 -- beta_Age drawn from distribution; alphas drawn from N(mu_cohort,sigmasq)
if(argv[1]==29){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    mean_alpha_prior=38    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta

    if(FALSE){
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.mat2,iter=10000,chains=0)
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
 print(model1.6.mat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
 #    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.all_chains.mat.including_alpha.RData")
}
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family

    mean_alpha_prior=30    
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 15
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
if(FALSE){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.pat2,iter=10000,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    #load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.all_chains.pat.including_alpha.RData")
}
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
if(myname=="maternal"){my.ylim=9}else {my.ylim=14}
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
   lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
   abline(v=0,lwd=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
   if(myname=="maternal"){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
   }
    dev.off()

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.mu_m.pdf"),height=5,width=5)
if(myname=="maternal"){my.ylim=0.3}else {my.ylim=0.5}
for(i in 1:ncol(mysim$mu_m)){
    if(i==1){
        plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
    curve(dnorm(x,mean_alpha_prior,sqrt(5)),add=T,lwd=3,lty=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()


library(pscl)
    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.2
        my.xlim=100
    }else {
        my.ylim=0.6
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i],n=256),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i],n=256),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
if(FALSE){
    emp.variances=read.delim("variance_in_mean_duoHMM_count_by_parent.informative_meioses_only.min_2_kids.txt",header=T,stringsAsFactors=F)
    for(z in 1:nrow(cohort.codes)){
        if(argv[1]==29){
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"maternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        } else {
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"paternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        }
    }
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior","empirical variance in parent means"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),2,2),lwd=c(rep(2,nrow(cohort.codes)),3,1),cex=0.5)
}
legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
    lines(density(rigamma(10000,3,0.5)),lwd=2,lty=2)
legend("topright",c("posterior","prior"),lty=c(1,2),lwd=c(2,2))
dev.off()
pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.tausq.pdf"),height=5,width=5)
    plot(density(mysim$tausq),xlab="tausq",main="Posterior for tausq",lwd=2)
    lines(density(rigamma(10000,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

}



if(argv[1] ==31 | argv[1]==32){
#### Model 1.7 -- beta_Age drawn from distribution; with family effect and cohort effect; using t distribution for Y and prior on beta_Age
if(argv[1]==31){
    data1.mat2$df_y <- 5
    data1.mat2$df_beta_age <- 5

    model1.7.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.7.stan",data=data1.mat2,iter=10000,chains=0)
    model1.7.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.7.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","mu_m","sigmasq_m","a0")))
    model1.7.mat <- sflist2stanfit(model1.7.mat.list)
print(model1.7.mat,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.7.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.7.all_chains.mat.including_alpha.RData"))
#load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.7.all_chains.mat.including_alpha.RData")

    mysim<-extract(model1.7.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
    data1.pat2$df_y <- 5
    data1.pat2$df_beta_age <- 5

    model1.7.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.7.stan",data=data1.pat2,iter=10000,chains=0)
    model1.7.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.7.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","mu_m","sigmasq_m","a0")))
    model1.7.pat <- sflist2stanfit(model1.7.pat.list)
print(model1.7.pat,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.7.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.7.all_chains.pat.including_alpha.RData"))
    #load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.7.all_chains.pat.including_alpha.RData")

    mysim<-extract(model1.7.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output_with_NFTOOLS/model1.7.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
for(i in 1:ncol(mysim$beta_Age)){
    if(myname=="maternal"){my.ylim=7}else{my.ylim=12}
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
   if(myname=="maternal"){
       my.ylim=20
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
   }
    dev.off()
    
pdf(paste0("RSTAN_output_with_NFTOOLS/model1.7.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}

if(argv[1] ==33 | argv[1]==34){
#### Model 1.5 -- negative binomial model for informative families only
if(argv[1]==33){
    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.stan",data=data1.mat2,iter=10000,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
print(model1.5.mat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.all_chains.mat.including_alpha.RData"))
    #load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.all_chains.mat.including_alpha.RData")

    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
    model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.stan",data=data1.pat2,iter=10000,chains=0)
    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.pat <- sflist2stanfit(model1.5.pat.list)
print(model1.5.pat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.all_chains.pat.including_alpha.RData"))
    #load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.all_chains.pat.including_alpha.RData")
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
my.ylim=55
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
dev.off()
    

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}

if(argv[1] ==35 | argv[1]==36){
#### Model 1.5, with different link function    -- negative binomial model for informative families only, with different link function
    cat("#### Model 1.5, with different link function\n")

if(argv[1]==35){
    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.stan",data=data1.mat2,iter=10000,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
print(model1.5.mat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
   save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.different_link.all_chains.mat.including_alpha.RData"))
    #load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.different_link.all_chains.mat.including_alpha.RData")

    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
    model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.stan",data=data1.pat2,iter=10000,chains=0)
    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.pat <- sflist2stanfit(model1.5.pat.list)
 print(model1.5.pat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.different_link.all_chains.pat.including_alpha.RData"))
#load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.different_link.all_chains.pat.including_alpha.RData")
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}


    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.different_link.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=200}else{my.ylim=300}
    for(i in 1:ncol(mysim$beta_Age)){
        if(i==1){
            plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
    abline(v=0,lwd=2)
    dev.off()
    

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.different_link.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}

if(argv[1] ==37 | argv[1]==38){
#if(argv[1] ==37 ){
#### Model 1.5, with different link function  and weaker priors  -- negative binomial model for informative families only, with different link function
    cat("#### Model 1.5, with different link function\n")
if(argv[1]==37){
    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.weaker_priors.stan",data=data1.mat2,iter=10000,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
 print(model1.5.mat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.different_link.weaker_priors.all_chains.mat.including_alpha.RData"))
#load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.different_link.weaker_priors.all_chains.mat.including_alpha.RData")

    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
####Need to run this for longer - didn't converge
    model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.weaker_priors.stan",data=data1.pat2,iter=10000,chains=0)
    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
#    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=20000,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.pat <- sflist2stanfit(model1.5.pat.list)
#    
print(model1.5.pat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.different_link.weaker_priors.all_chains.pat.including_alpha.RData"))
#load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.different_link.weaker_priors.all_chains.pat.including_alpha.RData")
    
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.different_link.weaker_priors.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=200}else{my.ylim=200}
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
   lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
   abline(v=0,lwd=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()
    
pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.different_link.weaker_priors.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}

####################################################

if(argv[1] ==39 | argv[1]==40){
#### Model 1.3, weaker priors on age effect -- beta_Age drawn from distribution
cat("    #### Model 1.3, weaker priors on age effect -- beta_Age drawn from distribution\n")
if(argv[1]==39){
    model1.3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.3.weaker_priors_on_age_effect.stan",data=data1.mat2,iter=10000,chains=0)
    model1.3.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.3.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","mu_m","sigmasq_m","a0")))
    model1.3.mat <- sflist2stanfit(model1.3.mat.list)
    print(model1.3.mat,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.3.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.3.weaker_priors_on_age_effect.all_chains.mat.including_alpha.RData"))
#    #load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.3.weaker_priors_on_age_effect.all_chains.mat.including_alpha.RData")

    mysim<-extract(model1.3.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
    model1.3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.3.weaker_priors_on_age_effect.stan",data=data1.pat2,iter=10000,chains=0)
    model1.3.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.3.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","mu_m","sigmasq_m","a0")))
    model1.3.pat <- sflist2stanfit(model1.3.pat.list)
    print(model1.3.pat,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.3.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.3.weaker_priors_on_age_effect.all_chains.pat.including_alpha.RData"))
#    #load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.3.weaker_priors_on_age_effect.all_chains.pat.including_alpha.RData")
    mysim<-extract(model1.3.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.3.weaker_priors_on_age_effect.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
for(i in 1:ncol(mysim$beta_Age)){
if(myname=="maternal"){    my.ylim=7}else {my.ylim=12}
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
if(myname=="maternal"){
    polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")

legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
}
dev.off()

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.3.weaker_priors_on_age_effect.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()

}


if(argv[1] ==41 | argv[1]==42){
#### Model 1.3, uniform_prior_on_sigmasq_m. -- beta_Age drawn from distribution
cat(" #### Model 1.3, uniform_prior_on_sigmasq_m. -- beta_Age drawn from distribution\n")
if(argv[1]==41){
#    model1.3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.3.uniform_prior_on_sigmasq_m.stan",data=data1.mat2,iter=10000,chains=0)
#    model1.3.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.3.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","mu_m","sigmasq_m","a0")))
#    model1.3.mat <- sflist2stanfit(model1.3.mat.list)
#    print(model1.3.mat,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
#    save(model1.3.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.3.uniform_prior_on_sigmasq_m.all_chains.mat.including_alpha.RData"))
   load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.3.uniform_prior_on_sigmasq_m.all_chains.mat.including_alpha.RData")
   mysim<-extract(model1.3.mat,permuted=T)
   cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
   myname="maternal"
} else {
   # model1.3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.3.uniform_prior_on_sigmasq_m.stan",data=data1.pat2,iter=10000,chains=0)
   # model1.3.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.3.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","mu_m","sigmasq_m","a0")))
   # model1.3.pat <- sflist2stanfit(model1.3.pat.list)
   # print(model1.3.pat,pars=c("beta_Age","beta_global","sigmasq_global","beta_Cohort","tausq","sigmasq_m","mu_m"),digits=4)
   # save(model1.3.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.3.uniform_prior_on_sigmasq_m.all_chains.pat.including_alpha.RData"))
load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.3.uniform_prior_on_sigmasq_m.all_chains.pat.including_alpha.RData")
    mysim<-extract(model1.3.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.3.uniform_prior_on_sigmasq_m.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
for(i in 1:ncol(mysim$beta_Age)){
if(myname=="maternal"){    my.ylim=7} else {my.ylim=15}
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
if(myname=="maternal"){
polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")

legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
}
dev.off()

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.3.uniform_prior_on_sigmasq_m.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.3.uniform_prior_on_sigmasq_m.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_m),xlab="sigmasq_m",main="Posterior for sigmasq_m",lwd=2)
dev.off()

}

if(argv[1] ==43 | argv[1]==44){
#### Model 1.62 -- common beta_Age for all cohorts; N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently from IG for each cohort
    df.alpha=5
    if(argv[1]==43){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    data1.mat2$df_alpha = df.alpha
##    mean_alpha_prior=38
    mean_alpha_prior=41    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=10
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
    if(FALSE){   
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.stan",data=data1.mat2,iter=10000,chains=0)
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
#    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.mat.including_alpha.RData"))
#    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
 # load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.mat.including_alpha.RData")
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData")) 
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha
#    mean_alpha_prior=30
     mean_alpha_prior=27
     data1.pat2$mean_alpha_prior = mean_alpha_prior
     sigmasq_m_alpha = 2
     sigmasq_m_beta = 15
     data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
     data1.pat2$sigmasq_m_beta = sigmasq_m_beta
     sigmasq_mu_m_prior=8
     data1.pat2$sigmasq_mu_m_prior=sigmasq_mu_m_prior
     if(FALSE){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.stan",data=data1.pat2,iter=10000,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
     save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
#    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
#    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.RData")
}
#     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}


    
pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    my.ylim=100
    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
    if(myname=="maternal"){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
   }
    dev.off()
    
    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=0.6}else {my.ylim=0.7}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
#//        curve(dnorm(x,mean_alpha_prior,sqrt(5)),add=T,lwd=3,lty=2)
    curve(dnorm(x,mean_alpha_prior,sigmasq_mu_m_prior),add=T,lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()
    library(pscl)
    
    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.2
        my.xlim=100
    }else {
        my.ylim=0.7
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
if(FALSE){
    emp.variances=read.delim("variance_in_mean_duoHMM_count_by_parent.informative_meioses_only.min_2_kids.txt",header=T,stringsAsFactors=F)
    for(z in 1:nrow(cohort.codes)){
        if(argv[1]==43){
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"maternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        } else {
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"paternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        }
    }
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior","empirical variance in parent means"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),2,2),lwd=c(rep(2,nrow(cohort.codes)),3,1),cex=0.5)
}
       legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.tausq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),lwd=2,main="Posterior for tausq",xlab="tausq")
    lines(density(rigamma(10000,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)
    dev.off()
}

if(argv[1] ==45 | argv[1]==46){
#### Model 1.5 -- negative binomial model for informative families only; common beta_Age for all cohorts; alpha drawn from N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently from IG for each cohort     
if(argv[1]==45){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    mean_alpha_prior=36
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 5
    sigmasq_m_beta = 5
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
if(FALSE){
    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.stan",data=data1.mat2,iter=10000,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
print(model1.5.mat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
#    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.all_chains.mat.including_alpha.RData")
}
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    mean_alpha_prior=33
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 5
    sigmasq_m_beta = 5
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
    if(FALSE){
        model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.stan",data=data1.pat2,iter=10000,chains=0)
        model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
        model1.5.pat <- sflist2stanfit(model1.5.pat.list)
        print(model1.5.pat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
        save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
                                        #    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.all_chains.pat.including_alpha.RData")
}
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.4.",myname,".posteriors.beta_Age.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
abline(v=0,lwd=2)
dev.off()

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.4.",myname,".posteriors.mu_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
      if(myname=="maternal"){my.ylim=1.6}else {my.ylim=1.6}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sd=sqrt(6)),lty=2,lwd=3,add=T)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

library(pscl)
    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.4.",myname,".posteriors.sigmasq_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)

if(myname=="maternal"){
    my.ylim=3
    my.xlim=4
}else{
    my.ylim=5
    my.xlim=4
}
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lty=2,lwd=3)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.4.",myname,".posteriors.omega.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
dev.off()


}


if(argv[1] ==47 | argv[1]==48|argv[1]==55|argv[1]==56|argv[1]==57|argv[1]==58|argv[1]==59|argv[1]==60|argv[1]==61|argv[1]==62){
#### Model 1.62 -- common beta_Age for all cohorts; N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently; uniform priors on all parameters
    df.alpha=5
    if(argv[1]==47|argv[1]==55|argv[1]==57|argv[1]==59|argv[1]==61){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    data1.mat2$df_alpha = df.alpha
        mean_alpha_prior=38    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
#if(FALSE){
    if(argv[1]==47){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors.stan",data=data1.mat2,iter=10000,chains=0)
    }
    if(argv[1]==55){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_beta_Age.stan",data=data1.mat2,iter=10000,chains=0)
    }
    if(argv[1]==57){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_mu_m_cohort.stan",data=data1.mat2,iter=10000,chains=0)
    }
    if(argv[1]==59){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_tau.stan",data=data1.mat2,iter=10000,chains=0)
    }
    if(argv[1]==61){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_sigma_m.stan",data=data1.mat2,iter=10000,chains=0)
    }
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=20000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
#}
    if(argv[1]==47){
                myname="uniform_priors.maternal"
#    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors.all_chains.mat.including_alpha.RData"))
#        save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors.v4.all_chains.mat.including_alpha.RData"))
#     load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors.all_chains.mat.including_alpha.RData")
#    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors.v4.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==55){
        myname="uniform_priors_on_beta_Age.maternal"
        save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.mat.including_alpha.RData"))
#               load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==57){
        myname="uniform_priors_on_mu_m_cohort.maternal"
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.mat.including_alpha.RData"))
 #       load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==59){
        myname="uniform_priors_on_tau.maternal"
        save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.mat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==61){
        myname="uniform_priors_on_sigma_m.maternal"
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.mat.including_alpha.RData"))
 #       load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.mat.including_alpha.RData"))
    }
    
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)

} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha

        mean_alpha_prior=33
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 15
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
#if(FALSE){
    if(argv[1]==48){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors.stan",data=data1.pat2,iter=10000,chains=0)
    }
    if(argv[1]==56){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_beta_Age.stan",data=data1.pat2,iter=10000,chains=0)
    }
    if(argv[1]==58){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_mu_m_cohort.stan",data=data1.pat2,iter=10000,chains=0)
    }
    if(argv[1]==60){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_tau.stan",data=data1.pat2,iter=10000,chains=0)
    }
    if(argv[1]==62){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_sigma_m.stan",data=data1.pat2,iter=10000,chains=0)
    }

    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=20000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
#}
    if(argv[1]==48){
        myname="uniform_priors.paternal"
    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors.v4.all_chains.pat.including_alpha.RData"))
                                        #    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors.all_chains.pat.including_alpha.RData")
#    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors.v4.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==56){
        myname="uniform_priors_on_beta_Age.paternal"
        save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.pat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==58){
        myname="uniform_priors_on_mu_m_cohort.paternal"
        save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.pat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==60){
        save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.pat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.pat.including_alpha.RData"))
        myname="uniform_priors_on_tau.paternal"
    }
    if(argv[1]==62){
        save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.pat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.pat.including_alpha.RData"))
        myname="uniform_priors_on_sigma_m.paternal"
    }
    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)

}


#pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.uniform_priors.v4.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
    my.ylim=100
    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
    #    if(argv[1] %in% c(47,48,55,56,57,58,59,60,61,62)){
    if(argv[1] %in% c(57,58,59,60,61,62)){
        curve(dnorm(x,0,1),add=T,lwd=2,lty=2)
        legend("topleft",c("posterior","prior"),lty=c(1,2),lwd=2)
    }
    if(argv[1] %in% c(47,55,57,59,61)){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
   }
    dev.off()
    
#    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.uniform_priors.v4.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
        pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
        if(argv[1] %in% c(47,55,57,59,61)){
        my.ylim=0.6
        my.xlim=c(32,46)
    }else {
        my.ylim=0.7
        my.xlim=c(23,32)
    }

    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=my.xlim,xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    if(argv[1] %in% c(55,56,59,60,61,62)){
        curve(dnorm(x,mean_alpha_prior,sqrt(5)),add=T,lwd=2,lty=2)
    }
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),2),cex=0.5)
    dev.off()
   library(pscl)
#    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.uniform_priors.v4.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
        pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
        if(argv[1] %in% c(47,55,57,59,61)){
        my.ylim=0.2
        my.xlim=100
    }else {
        my.ylim=0.7
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    if(argv[1] %in% c(55,56,57,58,59,60)){
        lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    }
if(FALSE){
    emp.variances=read.delim("variance_in_mean_duoHMM_count_by_parent.informative_meioses_only.min_2_kids.txt",header=T,stringsAsFactors=F)
    for(z in 1:nrow(cohort.codes)){
        if(argv[1]==47|argv[1]==55|argv[1]==57|argv[1]==59|argv[1]==61){
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"maternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        } else {
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"paternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        }
    }

    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),2),cex=0.5)
}
    dev.off()

#    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.uniform_priors.v4.",myname,".posteriors.tausq.pdf"),height=5,width=5)
        pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.tausq.pdf"),height=5,width=5)
    plot(density(mysim$tausq),main="Posterior for tausq",lwd=2,xlab="tausq")
    if(argv[1] %in% c(55,56,57,58,61,62)){
        lines(density(rigamma(10000,2,70)),lwd=3,lty=2)
    }
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)

    dev.off()
    
}


if(argv[1] ==49 | argv[1]==50){
#### Model 1.5 -- negative binomial model for informative families only; common beta_Age for all cohorts; alpha drawn from N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently  for each cohort ; uniform priors
if(argv[1]==49){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
if(FALSE){
   model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.uniform_priors.stan",data=data1.mat2,iter=10000,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
print(model1.5.mat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.uniform_priors.v4.all_chains.mat.including_alpha.RData"))
}
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.uniform_priors.v4.all_chains.mat.including_alpha.RData"))
#    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.uniform_priors.all_chains.mat.including_alpha.RData")
#    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.uniform_priors.all_chains.mat.including_alpha.RData"))
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family

    data1.pat2$mean_alpha_prior = 33
    if(FALSE){
        model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.uniform_priors.stan",data=data1.pat2,iter=20000,chains=0)
        model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
        model1.5.pat <- sflist2stanfit(model1.5.pat.list)
        print(model1.5.pat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
        save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.uniform_priors.v4.all_chains.pat.including_alpha.RData"))
                                        #    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.uniform_priors.all_chains.pat.including_alpha.RData")
#    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.uniform_priors.all_chains.pat.including_alpha.RData"))
    }
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.5.4.uniform_priors.v4.all_chains.pat.including_alpha.RData"))
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.4.uniform_priors.v4.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
abline(v=0,lwd=2)
dev.off()


pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.4.uniform_priors.v4.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
if(myname=="maternal"){
    my.ylim=1.7
} else {
    my.ylim=1.7
}
                                        #    if(myname=="maternal"){my.ylim=0.3}else {my.ylim=0.5}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.4.uniform_priors.v4.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=1
        my.xlim=10
    }else {
        my.ylim=10
        my.xlim=7
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2])),col=c(mycols[as.character(cohort.codes[,2])]),lty=c(rep(1,nrow(cohort.codes))),lwd=c(rep(2,nrow(cohort.codes))),cex=0.5)
dev.off()
    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.5.4.uniform_priors.v4.",myname,".posteriors.omega.pdf"),height=5,width=5)
plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
dev.off()


}


if(argv[1] ==51 | argv[1]==52){
#### Model 1.62 -- common beta_Age for all cohorts; N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently from IG for each cohort
    df.alpha=5
    if(argv[1]==51){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    data1.mat2$df_alpha = df.alpha
    mean_alpha_prior=38    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    
if(FALSE){
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas.stan",data=data1.mat2,iter=10000,chains=0)
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".RData"))

# load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".RData"))
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha
    mean_alpha_prior=30
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 15
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
 if(FALSE){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas.stan",data=data1.pat2,iter=10000,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distributions_on_alphas.df_",df.alpha,".RData"))
}
     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distributions_on_alphas.df_",df.alpha,".RData"))
#    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.RData")
#    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.beta_Age.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".pdf"),height=5,width=5)
    my.ylim=100
    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
    if(myname=="maternal"){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
   }
    dev.off()
    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.mu_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=0.6}else {my.ylim=0.7}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sqrt(5)),add=T,lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()
    library(pscl)
    
    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.sigmasq_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.2
        my.xlim=100
    }else {
        my.ylim=0.7
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.tausq.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),xlab="tausq",main="Posterior for tausq",lwd=2)
    lines(density(rigamma(10000,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()
}

if(argv[1] ==53 | argv[1]==54){
#### Model 1.62 -- common beta_Age for all cohorts; N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently from IG for each cohort
    df.alpha=5
    if(argv[1]==53){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    data1.mat2$df_alpha = df.alpha
    mean_alpha_prior=38    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
if(FALSE){
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas_and_Y.stan",data=data1.mat2,iter=10000,chains=0)
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".RData"))
}
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".RData"))
# load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="maternal"
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha
    mean_alpha_prior=30
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 15
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
 if(FALSE){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas_and_Y.stan",data=data1.pat2,iter=10000,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distributions_on_alphas_and_Y.df_",df.alpha,".RData"))
#    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.RData")
}
     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distributions_on_alphas_and_Y.df_",df.alpha,".RData"))
#     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.modell1.6.2.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.beta_Age.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".pdf"),height=5,width=5)
    my.ylim=100
    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
    if(myname=="maternal"){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
   }
    dev.off()
    
    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.mu_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=0.6}else {my.ylim=0.7}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sqrt(5)),add=T,lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()
    library(pscl)
    
    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.sigmasq_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.2
        my.xlim=100
    }else {
        my.ylim=0.7
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()
    pdf(paste0("RSTAN_output_with_NFTOOLS/model1.6.2.",myname,".posteriors.tausq.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),xlab="tausq",main="Posterior for tausq",lwd=2)
    lines(density(rigamma(10000,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

}else if(argv[1] ==63|argv[1]==64){
###fit model to adjusted counts y_adj~N(beta_Age * age_adj + beta_nkids*nkids,tausq_cohort)
    if(argv[1]==63){
    sigmasq_beta_nkids=5
    data1.adjusted.mat2$sigmasq_beta_nkids=sigmasq_beta_nkids
    tausq_alpha_prior=2
     data1.adjusted.mat2$tausq_alpha_prior=tausq_alpha_prior
    tausq_beta_prior=70
     data1.adjusted.mat2$tausq_beta_prior=tausq_beta_prior
    myname=paste0("sigmasq_beta_nkids_",sigmasq_beta_nkids,".tausq_",tausq_alpha_prior,"_",tausq_beta_prior,".maternal")
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T) #should be model 1 since we're only including families with >1 kid
    if(FALSE){
                                        #common age effect
    model2.3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model2.3.stan",data=data1.adjusted.mat2,iter=10,chains=0,pars=c("beta_Age","beta_nkids","tausq"))
    model2.3.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model2.3.mat.compiled, data = data1.adjusted.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_nkids","tausq")))
    model2.3.mat <- sflist2stanfit(model2.3.mat.list)

    print(model2.3.mat,pars=c("beta_Age","beta_nkids","tausq"),digits=4)
    save(model2.3.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model2.3.",myname,".RData"))
}
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model2.3.",myname,".RData"))
    
    mysim<-extract(model2.3.mat,permuted=T)

}else{
    sigmasq_beta_nkids=5
    data1.adjusted.pat2$sigmasq_beta_nkids=sigmasq_beta_nkids
    tausq_alpha_prior=2
     data1.adjusted.pat2$tausq_alpha_prior=tausq_alpha_prior
    tausq_beta_prior=30
     data1.adjusted.pat2$tausq_beta_prior=tausq_beta_prior
    myname=paste0("sigmasq_beta_nkids_",sigmasq_beta_nkids,".tausq_",tausq_alpha_prior,"_",tausq_beta_prior,".paternal")
    if(FALSE){
    model2.3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model2.3.stan",data=data1.adjusted.pat2,iter=10,chains=0,pars=c("beta_Age","beta_nkids","tausq"))
    model2.3.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model2.3.pat.compiled, data = data1.adjusted.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_nkids","tausq")))
    model2.3.pat <- sflist2stanfit(model2.3.pat.list)
    print(model2.3.pat,pars=c("beta_Age","beta_nkids","tausq"),digits=4)
    save(model2.3.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model2.3.",myname,".RData"))
}
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model2.3.",myname,".RData"))
    mysim<-extract(model2.3.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
}

pdf(paste("RSTAN_output_with_NFTOOLS/model2.3.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
my.ylim=12
plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
        abline(v=0,lwd=2)
    if(argv[1]==63){

        polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
        abline(v=0.067,lty=2,lwd=2)
        polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
        abline(v=0.082,lty=3,lwd=2)
        polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
        abline(v=0.19,lty=4,lwd=2)
        abline(v=-0.42,lty=4,lwd=2,col="grey")
        legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
    }
    abline(v=0,lwd=2)
    curve(dnorm(x,0,1),lty=2,lwd=2,add=T)
    legend("topleft",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS/model2.3.",myname,".posteriors.beta_nkids.pdf"),height=5,width=5)
    plot(density(mysim$beta_nkids),xlim=range(mysim$beta_Age),xlab="beta_nkids",main="Posterior for beta_nkids",lwd=2)
    curve(dnorm(x,0,sigmasq_beta_nkids),lty=2,lwd=2,add=T)
    abline(v=0,lwd=2)
    legend("topleft",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

library(pscl)
pdf(paste0("RSTAN_output_with_NFTOOLS/model2.3.",myname,".posteriors.tausq.pdf"),height=5,width=5)
    my.ylim=0.75
    for(i in 1:ncol(mysim$tausq)){
        if(i==1){
            plot(density(mysim$tausq[,i]),xlim=range(mysim$tausq),xlab="tausq",main="Posterior for tausq",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$tausq[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(10000,tausq_alpha_prior,tausq_beta_prior)),lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()
        

}else if(argv[1] ==65|argv[1]==66){
###fit model to adjusted counts y_adj~N(beta_Age * age_adj + beta_nkids*nkids,tausq_cohort)
    if(argv[1]==65){
    sigmasq_beta_nkids=5
    data1.adjusted.mat2$sigmasq_beta_nkids=sigmasq_beta_nkids
    tausq_alpha_prior=2
     data1.adjusted.mat2$tausq_alpha_prior=tausq_alpha_prior
    tausq_beta_prior=70
     data1.adjusted.mat2$tausq_beta_prior=tausq_beta_prior
    myname=paste0("sigmasq_beta_nkids_",sigmasq_beta_nkids,".tausq_",tausq_alpha_prior,"_",tausq_beta_prior,".maternal")
if(FALSE){
                                        #common age effect
    model2.3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model2.4.stan",data=data1.adjusted.mat2,iter=10,chains=0,pars=c("beta_Age","beta_nkids","tausq"))
    model2.3.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model2.3.mat.compiled, data = data1.adjusted.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_nkids","tausq")))
    model2.4.mat <- sflist2stanfit(model2.3.mat.list)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T) #should be model 1 since we're only including families with >1 kid

    print(model2.4.mat,pars=c("beta_Age","beta_nkids","tausq"),digits=4)
    save(model2.4.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model2.4.",myname,".RData"))
}
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model2.4.",myname,".RData"))
    
    mysim<-extract(model2.4.mat,permuted=T)

}else{
    sigmasq_beta_nkids=5
    data1.adjusted.pat2$sigmasq_beta_nkids=sigmasq_beta_nkids
    tausq_alpha_prior=2
     data1.adjusted.pat2$tausq_alpha_prior=tausq_alpha_prior
    tausq_beta_prior=30
     data1.adjusted.pat2$tausq_beta_prior=tausq_beta_prior
    myname=paste0("sigmasq_beta_nkids_",sigmasq_beta_nkids,".tausq_",tausq_alpha_prior,"_",tausq_beta_prior,".paternal")
    if(FALSE){
        model2.3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model2.4.stan",data=data1.adjusted.pat2,iter=10,chains=0,pars=c("beta_Age","beta_nkids","tausq"))
        model2.3.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model2.3.pat.compiled, data = data1.adjusted.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_nkids","tausq")))
        model2.4.pat <- sflist2stanfit(model2.3.pat.list)
        print(model2.4.pat,pars=c("beta_Age","beta_nkids","tausq"),digits=4)
        save(model2.4.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model2.4.",myname,".RData"))
    }
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/running_RSTAN.model2.4.",myname,".RData"))
    mysim<-extract(model2.4.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
}

pdf(paste("RSTAN_output_with_NFTOOLS/model2.4.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
my.ylim=12
plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
if(argv[1]==65){
        abline(v=0,lwd=2)
        polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
        abline(v=0.067,lty=2,lwd=2)
        polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
        abline(v=0.082,lty=3,lwd=2)
        polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
        abline(v=0.19,lty=4,lwd=2)
        abline(v=-0.42,lty=4,lwd=2,col="grey")
        legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
    }
    abline(v=0,lwd=2)
    curve(dnorm(x,0,1),lty=2,lwd=2,add=T)
    legend("topleft",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS/model2.4.",myname,".posteriors.beta_nkids.pdf"),height=5,width=5)
    plot(density(mysim$beta_nkids),xlim=range(mysim$beta_Age),xlab="beta_nkids",main="Posterior for beta_nkids",lwd=2)
    curve(dnorm(x,0,sigmasq_beta_nkids),lty=2,lwd=2,add=T)
    abline(v=0,lwd=2)
    legend("topleft",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

library(pscl)
    pdf(paste0("RSTAN_output_with_NFTOOLS/model2.4.",myname,".posteriors.tausq.pdf"),height=5,width=5)
    plot(density(mysim$tausq),xlim=range(mysim$tausq),xlab="tausq",main="Posterior for tausq",col="black",lwd=2)
    lines(density(rigamma(10000,tausq_alpha_prior,tausq_beta_prior)),lwd=2,lty=2)
    legend("topleft",c("posterior","prior"),lty=c(1,2),lwd=2,cex=0.5)
    dev.off()
}

if(argv[1] ==69 | argv[1]==70|argv[1] ==71 | argv[1]==72){
        ######## add in 73-76 - using data0 (only informative nuclear families)
#### Model 1.82 -- common beta_Age for all cohorts; N(mu_m,sigmasq_m) i.e. like 1.6.2 but no cohort effect - so we can compare to model 1.8 which has cohort-specific age effects but needs to avoid identifiability problem
    if(as.numeric(argv[1]) %%2 !=0){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family

    mean_alpha_prior=41    
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
#    sigmasq_mu_m_prior=100
    sigmasq_mu_m_prior=4
    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
#    if(FALSE){
    if(argv[1]==69){
#if(FALSE){#for running jsut FCs alone
        data1.mat2$y=data1.mat2$y[data1.mat2$cohort==1]
        data1.mat2$Age=data1.mat2$Age[data1.mat2$cohort==1]
        data1.mat2$family=data1.mat2$family[data1.mat2$cohort==1]
        data1.mat2$J=length(data1.mat2$Age)
        data1.mat2$C=1
        data1.mat2$I=length(unique(data1.mat2$family))       
        data1.mat2$cohort=data1.mat2$cohort[data1.mat2$cohort==1]
 #   }
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.2.stan",data=data1.mat2,iter=10000,chains=0)
#        myname="model1.8.2.maternal"
        myname="model1.8.2.maternal.FCs_alone.tighter_prior_on_mu_m"
        model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
        cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    }
    if(argv[1]==71){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.stan",data=data1.mat2,iter=10000,chains=0)
#        myname="model1.8.maternal"
        myname="model1.8.maternal.uniform_prior_on_beta_Age"
#           model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0","beta_global","sigmasq_global")))
        model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
            cohort.codes=read.delim("RSTAN_output/key_for_maternal_NFTOOLS_cohorts_to_include_in_model_1.txt",header=T)
    }
     
 model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.",myname,".all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior, ".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
 
#    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.",myname,".all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior, ".sigmasq_m_IG_",
#                sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
#}
    mysim<-extract(model1.6.mat,permuted=T)
    if("beta_global" %in% names(mysim)){
        print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
    } else {
        print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    }
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
     mean_alpha_prior=27
     sigmasq_m_alpha = 2
     sigmasq_m_beta = 15
     sigmasq_mu_m_prior=8^2
     
     data1.pat2$mean_alpha_prior = mean_alpha_prior
     data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
     data1.pat2$sigmasq_m_beta = sigmasq_m_beta
     data1.pat2$sigmasq_mu_m_prior=sigmasq_mu_m_prior


                                        #     if(FALSE){
     if(argv[1] ==70){
         model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.2.stan",data=data1.pat2,iter=10000,chains=0)
         myname="model1.8.2.paternal"
         model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
         cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
     }
     if(argv[1]==72){
         model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.stan",data=data1.pat2,iter=10000,chains=0)
         myname="model1.8.paternal"
         model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
             cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.txt",header=T)
     }
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
     save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.",myname,".all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
     #load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.",myname,".all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
      print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
     mysim<-extract(model1.6.pat,permuted=T)
 }

pdf(paste0("RSTAN_output_with_NFTOOLS/",myname,".traceplots.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    traceplot(model1.6.mat,pars=c("beta_Age"),window=c(10^2,10^4))
    traceplot(model1.6.mat,pars=c("mu_m"),window=c(10^2,10^4))
    traceplot(model1.6.mat,pars=c("sigmasq_m"),window=c(10^2,10^4))
    traceplot(model1.6.mat,pars=c("tausq"),window=c(10^2,10^4))
    if("beta_global" %in% names(mysim)){
        traceplot(model1.6.mat,pars=c("beta_global"),window=c(10^2,10^4))
        traceplot(model1.6.mat,pars=c("sigmasq_global"),window=c(10^2,10^4))
    }

    dev.off()
    
    
pdf(paste0("RSTAN_output_with_NFTOOLS/",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    my.ylim=100
#    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)

    if(!"beta_global" %in% names(mysim)){
        plot(density(mysim$beta_Age),xlim=c(-0.5,0.5),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
        abline(v=0,lwd=2)
        curve(dnorm(x,0,1),lty=2,lwd=2,add=T)
        legend("topleft",c("posterior","prior"),lty=c(1,2),lwd=2)
    } else {
        for(i in 1:ncol(mysim$beta_Age)){
            if(i==1){
                plot(density(mysim$beta_Age[,i]),xlim=c(-0.5,0.5),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,10))
            }else {
                lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
            }
        }
        curve(dnorm(x,0,1),lty=2,lwd=2,add=T,col="black")
        abline(v=0,lwd=2)
        if("beta_global" %in% names(mysim)){
            lines(density(mysim$beta_global),lty=3,lwd=2,col="red")
            legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global","prior"),col=c(mycols[as.character(cohort.codes[,2])],"red","black"),lty=c(rep(1,nrow(cohort.codes)),3,2),lwd=2,cex=0.5)
        } else {
            legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=2,cex=0.5)
        }
    }

    if(as.numeric(argv[1]) %%2 !=0){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
#       legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
       legend("topright",c("Campbell","Kong","Coop","Hussin"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
   }
    dev.off()
    
    pdf(paste0("RSTAN_output_with_NFTOOLS/",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$mu_m),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",lwd=2)
    curve(dnorm(x,mean_alpha_prior,sigmasq_mu_m_prior),lwd=2,lty=2,add=T)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

    library(pscl)
        pdf(paste0("RSTAN_output_with_NFTOOLS/",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$sigmasq_m),xlab="sigmasq_m",main="Posterior for sigmasq_m",lwd=2)
    lines(density(rigamma(10000,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

    pdf(paste0("RSTAN_output_with_NFTOOLS/",myname,".posteriors.tausq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),lwd=2,main="Posterior for tausq",xlab="tausq")
    lines(density(rigamma(10000,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)
    dev.off()
}




if(FALSE){
###################### HAVEN'T FINISHED THIS SECTION OF CODE - MAYBE VARYING P_BY_COHORT IS NOT SO IMPORTANT
    
if(argv[1] ==43 | argv[1]==44){
    #Model 3, uniform_p_priors, with adjusted priors
cat("Model 3, uniform_p_priors, with adjusted priors\n")
    if(argv[1] ==43){
        model3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.uniform_p_priors.stan",data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
        model3.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.mat.compiled, data = data2.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
        model3.mat <- sflist2stanfit(model3.mat.list)
        save(model3.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.model3.uniform_p_priors.all_chains.mat.including_alpha.RData"))
        cohort.codes=read.delim("RSTAN_output_with_NFTOOLS/key_for_maternal_cohorts_to_include_in_model_3.NTR_v2.txt",header=T)
    }else {
    model3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.uniform_p_priors.stan",data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
    model3.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.pat.compiled, data = data2.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
    model3.pat <- sflist2stanfit(model3.pat.list)
    save(model3.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_with_NFTOOLS/RSTAN.model3.uniform_p_priors.all_chains.pat.including_alpha.RData"))
    cohort.codes=read.delim("RSTAN_output_with_NFTOOLS/key_for_paternal_cohorts_to_include_in_model_3.NTR_v2.txt",header=T)
}

pdf("RSTAN_output_with_NFTOOLS/model1.3.maternal.posteriors.beta_Age.pdf",height=5,width=5)
my.ylim=7
for(i in 1:ncol(mysim$beta_Age)){
if(i==1){
    plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
}else {
    lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
abline(v=0,lwd=2)
polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
abline(v=0.067,lty=2,lwd=2)
polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
abline(v=0.082,lty=3,lwd=2)
polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
abline(v=0.19,lty=4,lwd=2)
abline(v=-0.42,lty=4,lwd=2,col="grey")
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
legend("topright",c("Adam","decode","Hutterites","Julie"),lty=c(2,3,4,4),lwd=2,col=c("black","black","black","grey"),cex=0.5)
}
dev.off()
pdf("RSTAN_output_with_NFTOOLS/model1.3.maternal.posteriors.sigmasq_global.pdf",height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()

}}
### Model 4. N.B. NOT SURE SHOULD INCLUDE FAMILIES WITH ONLY 1 KID - DOESN'T REALLY LOOK BINOMIAL
#for all meioses, NO FAMILY EFFECT; fit N_called_crossovers~Binom(Y,p=f(cohort*family_type)); only for those with
#data3.mat2=list(y=data3.mat$nrec,cohort=as.factor(data3.mat$cohort),age=as.numeric(data3.mat$age.at.birth),J=nrow(data3.mat),I=length(unique(data3.mat$PARENT)))
#data3.pat2=list(y=data3.pat$nrec,cohort=as.factor(data3.pat$cohort),age=as.numeric(data3.pat$age.at.birth),J=nrow(data3.pat),I=length(unique(data3.pat$PARENT)))

