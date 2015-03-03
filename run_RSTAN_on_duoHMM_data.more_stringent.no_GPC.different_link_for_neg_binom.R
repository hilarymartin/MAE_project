library(rstan)
library(scales)
library(parallel)
library(pscl)
load("/well/donnelly/hilary/maternal_age_and_recombination/duoHMM_data_for_RSTAN.more_stringent.no_GPC.RData")
argv <- commandArgs(trailingOnly = TRUE)
print(argv)    
mycols=c("black","blue","red","green","orange","red4","purple","darkgreen","violetred2")
names(mycols)=c("CARL","FC","FVG","GPC","NTR","QTR370","QTR610","VB","ORCADES")
niterations=10000

print.only=as.logical(as.numeric(argv[2]))

if(argv[1]==1|argv[1]==2){

  mu_alpha_prior = 3.7
  sigmasq_alpha_prior = 0.2

  data2.mat2$mu_alpha_prior = mu_alpha_prior
  data2.mat2$sigmasq_alpha_prior = sigmasq_alpha_prior

  if(argv[1] ==1){
    cat("Model 3 maternal, with adjusted priors\n")
    if(!print.only){
     model3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.beta_age_by_cohort.different_link.stan",data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
     model3.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.pat.compiled, data = data2.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
     model3.mat <- sflist2stanfit(model3.pat.list)  
        save(model3.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.beta_age_by_cohort.adjusted_priors.all_chains.mat.including_alpha.RData"))
        print(model3.mat,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
    } else {
        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.beta_age_by_cohort.adjusted_priors.all_chains.mat.including_alpha.RData")
    }
    mysim<-extract(model3.mat,permuted=T)
    parent="adjusted_priors.maternal"
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.more_stringent.no_GPC.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_maternal.more_stringent.no_GPC.txt",header=T)[,1:5]
    cohort.codes2=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.more_stringent.no_GPC.original_cohorts.txt",header=T,stringsAsFactors=F)

}
if(argv[1] ==2){
  mu_alpha_prior = 3.2
  sigmasq_alpha_prior = 0.2

  data2.pat2$mu_alpha_prior = mu_alpha_prior
  data2.pat2$sigmasq_alpha_prior = sigmasq_alpha_prior
  
  
  cat("Model 3 paternal, with adjusted priors\n")
    if(!print.only){
     model3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.beta_age_by_cohort.different_link.stan",data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
     model3.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.pat.compiled, data = data2.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
     model3.pat <- sflist2stanfit(model3.pat.list)  
     save(model3.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.beta_age_by_cohort.adjusted_priors.all_chains.pat.including_alpha.RData"))
     print(model3.pat,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
    } else {
        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.beta_age_by_cohort.adjusted_priors.all_chains.pat.including_alpha.RData")
    }
    mysim<-extract(model3.pat,permuted=T)
    parent="adjusted_priors.paternal"
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.more_stringent.no_GPC.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_paternal.more_stringent.no_GPC.txt",header=T)[,1:5]
    cohort.codes2=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.more_stringent.no_GPC.original_cohorts.txt",header=T,stringsAsFactors=F)
}

cohort.codes=cohort.codes[order(cohort.codes[,1]),]
cohort.codes2=cohort.codes2[order(cohort.codes2[,1]),]
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

    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model3.beta_age_by_cohort.",parent,".posteriors.beta_Age.pdf"),height=5,width=5)
    my.ylim=45
    for(i in 1:ncol(mysim$beta_Age)){
        if(i==1){
            plot(density(mysim$beta_Age[,i]),xlim=c(-0.15,0.1),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes2[cohort.codes2[,1]==i,2])],ylim=c(0,my.ylim),
                 lwd=2,lty=1)
        }else {
            lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes2[i,2])],lwd=2,lty=1)
        }
    }
        lines(density(mysim$beta_global),col="black",lwd=3,lty=4)
        curve(dnorm(x,0,sqrt(0.05)),add=T,lwd=2,lty=2)

legend("topleft",c(names(mycols),"posterior","posterior global","prior global"),col=c(mycols,"black","black","black"),lty=c(rep(1,length(mycols)),1,4,2),lwd=c(rep(2,length(mycols)),2,2,2),cex=0.8,bg="white")
        abline(v=0,lwd=2)

dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model3.beta_age_by_cohort.",parent,".posteriors.all_parameters.pdf"),height=7*3,width=4*3)
par(mfrow=c(7,4),oma=c(0,0,1,0),mar=c(5,4,2,1))
#sigmasq_global
        plot(density(mysim$sigmasq_global),xlab=expression(sigma^2['global']),main="",lwd=2)
        lines(density(rigamma(niterations,3,0.1)),lty=2,lwd=2)
        legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
                                        #mu_m
        plot(density(mysim$mu_m),main="",lwd=2,xlab=expression(mu['c']),cex.lab=1.2)
        curve(dnorm(x,36,sqrt(6 )),lwd=2,lty=2,add=T)
        legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
#sigmasq_m
        plot(density(mysim$sigmasq_m),xlab=expression(sigma^2['c']),main="",lwd=2,cex.lab=1.2)
        lines(density(rigamma(niterations,5,5)),lty=2,lwd=2)
        legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))

#omega        
        plot(density(mysim$omega),main="",lwd=2,xlab=expression(omega),cex.lab=1.2)
        lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2,lwd=2)
        legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
#p
        mylabels=c(paste(c("2 generations, >2 children","3 generations, 2 children","2 generations, 2 children"),"both parents",sep=", "),paste(c("2 generations, >2 children","3 generations, 2 children","2 generations, 2 children"),"1 parent",sep=", "))
        names(mylabels)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
        for(i in 1:ncol(mysim$p_by_cohort)){
            sample.size= sum(data2.mat$cohort.family.type==cohort.codes[cohort.codes[,1]==i,"Cohort"])
            plot(density(mysim$p_by_cohort[,i]),xlab=expression('p'['f,c']),xlim=c(0,1),main="",lwd=2,cex.lab=1.2)
            mtext(cohort.codes[cohort.codes[,1]==i,"Pop"],3,cex=0.8,padj=0.3,at=0.5,font=2,line=1.5)
            mtext(mylabels[cohort.codes[cohort.codes[,1]==i,"Fam.type"]],3,cex=0.8,padj=0.3,at=0.5,font=2,line=0.5)
            curve(dbeta(x,beta.parameters[beta.parameters$Cohort.type==i,"alpha"],beta.parameters[beta.parameters$Cohort.type==i,"beta"]),add=T,lty=2,lwd=2)
            abline(v=beta.parameters[beta.parameters$Cohort.type==i,"mean.relative.to.all.informative"],lwd=1,lty=2,col="blue")
            legend("topleft",c("posterior","prior","p_observed",paste0("n = ",sample.size)),col=c("black","black","blue",NA),lty=c(1,2,2,NA),lwd=c(2,2,1,NA))
        }
        dev.off() 
}

#### Model 3.

if(argv[1]==25|argv[1]==26){
if(argv[1] ==25){
    cat("Model 3 maternal, with adjusted priors\n")
    if(!print.only){
 #    model3.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.adjusted_priors.stan",data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
#        model3.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.mat.compiled, data = data2.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
#        model3.mat <- sflist2stanfit(model3.mat.list)
 model3.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.adjusted_priors.stan",data=data2.mat2,iter=niterations,chains=4,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
        save(model3.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.adjusted_priors.all_chains.mat.including_alpha.RData"))
        print(model3.mat,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
    } else {
        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.adjusted_priors.all_chains.mat.including_alpha.RData")
    }
    mysim<-extract(model3.mat,permuted=T)
    parent="adjusted_priors.maternal"
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.more_stringent.no_GPC.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_maternal.more_stringent.no_GPC.txt",header=T)[,1:5]

}
if(argv[1] ==26){
    cat("Model 3 paternal, with adjusted priors\n")
if(!print.only){
#        model3.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.adjusted_priors.stan",data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
 #       model3.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.pat.compiled, data = data2.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
 #       model3.pat <- sflist2stanfit(model3.pat.list)
   model3.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.v2.adjusted_priors.stan",data=data2.pat2,iter=niterations,chains=4,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0" ))
        save(model3.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.adjusted_priors.all_chains.pat.including_alpha.RData"))
        print(model3.pat,pars=c("beta_Age","beta_global","sigmasq_global","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
    } else {
        load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.adjusted_priors.all_chains.pat.including_alpha.RData")
    }
    mysim<-extract(model3.pat,permuted=T)
    parent="adjusted_priors.paternal"
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.more_stringent.no_GPC.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_paternal.more_stringent.no_GPC.txt",header=T)[,1:5]
}

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

    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model3.",parent,".posteriors.beta_Age.pdf"),height=5,width=5)
    my.ylim=32
    for(i in 1:ncol(mysim$beta_Age)){
        if(i==1){
            plot(density(mysim$beta_Age[,i]),xlim=c(-0.2,0.2),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[cohort.codes[,1]==i,"Pop"])],ylim=c(0,my.ylim),
                 lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])])
        }else {
            lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,"Pop"])],lwd=mylwd[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])],lty=mylty[as.character(cohort.codes[cohort.codes[,1]==i,"Fam.type"])])
        }
        lines(density(mysim$beta_global),col="red",lwd=3,lty=3)
        curve(dnorm(x,0,sqrt(0.05)),add=T,lwd=2,lty=2)
        legend("topleft",c(names(mycols),"global","informative, 2 generations","informative, 3 generations","uninformative, 2 kids","both parents","one parent","prior"),col=c(mycols,"red",rep("black",5),"black"),lty=c(rep(1,length(mycols)),3,1,2,4,1,1,2),
               lwd=c(rep(2,length(mycols)),3,2,2,2,2,1,2),cex=0.5)
        abline(v=0,lwd=2)
    }
    dev.off()

    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model3.",parent,".posteriors.p.pdf"),height=5,width=5)
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
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model3.",parent,".posteriors.sigmasq_global.pdf"),height=5,width=5)
      plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
lines(rigamma(niterations,3,0.1),lty=2,lwd=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))

    dev.off()

    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model3.",parent,".posteriors.omega.pdf"),height=5,width=5)
    plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
    lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()
    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model3.",parent,".posteriors.mu_m.pdf"),height=5,width=5)
    plot(density(mysim$mu_m),main="Posterior for mu_m",lwd=2,xlab="mu_m")
    curve(dnorm(x,36,sqrt(6 )),lwd=2,lty=2,add=T)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()
  
    
}

####model 3.2

if(argv[1] %in% c(27,28)){
  
if(argv[1] ==27){
cat("Model 3.2 maternal, with adjusted priors\n")
  mu_alpha_prior = 3.7
  sigmasq_alpha_prior = 0.2
  sigmasq_m_alpha = 5
  sigmasq_m_beta = 5
  data2.mat2$sigmasq_m_alpha = sigmasq_m_alpha
  data2.mat2$sigmasq_m_beta = sigmasq_m_beta
  data2.mat2$mu_alpha_prior = mu_alpha_prior
  data2.mat2$sigmasq_alpha_prior = sigmasq_alpha_prior


if(!print.only){
model3.2.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.stan",data=data2.mat2,iter=niterations,chains=4,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
print(model3.2.mat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
save(model3.2.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.2.adjusted_priors.all_chains.mat.including_alpha.RData"))
} else {
load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.2.adjusted_priors.all_chains.mat.including_alpha.RData")

}

mysim<-extract(model3.2.mat,permuted=T)
parent="adjusted_priors.maternal"
cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.more_stringent.no_GPC.txt",header=T,stringsAsFactors=F)
beta.parameters=read.delim("parameters_for_beta_distribution.model3_maternal.more_stringent.no_GPC.txt",header=T)[,1:5]

}

if(argv[1] ==28){
    cat("Model 3.2 paternal, with adjusted priors\n")
  mu_alpha_prior = 3.2
  sigmasq_alpha_prior = 0.2
  sigmasq_m_alpha = 5
  sigmasq_m_beta = 5
  data2.pat2$sigmasq_m_alpha = sigmasq_m_alpha
  data2.pat2$sigmasq_m_beta = sigmasq_m_beta
  data2.pat2$mu_alpha_prior = mu_alpha_prior
  data2.pat2$sigmasq_alpha_prior = sigmasq_alpha_prior
if(!print.only){
    model3.2.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.stan",data=data2.pat2,iter=niterations,chains=4,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
    save(model3.2.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.2.adjusted_priors.all_chains.pat.including_alpha.RData"))

} else {
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.2.adjusted_priors.all_chains.pat.including_alpha.RData")
    print(model3.2.pat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
}
    mysim<-extract(model3.2.pat,permuted=T)
mysim.model3.pat = mysim
    sum(mysim.model3.pat$beta_Age>0)/length(mysim.model3.pat$beta_Age)
    parent="adjusted_priors.paternal"
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.more_stringent.no_GPC.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_paternal.more_stringent.no_GPC.txt",header=T)[,1:5]
}
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

 pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model3.2.",parent,".posteriors.beta_Age.pdf"),height=5,width=5)
plot(density(mysim$beta_Age),xlab=expression(beta['age']),main="",cex.lab=1.2,xlim=c(-0.03,0.04),ylim=c(0,60))
    curve(dnorm(x,0,sqrt(0.05)),lwd=2,lty=2,add=T)
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2,cex=0.8,bg="white")
    abline(v=0,lwd=2)
dev.off()

 pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model3.2.",parent,".posteriors.all_parameters.pdf"),height=8*3,width=4*3)
par(mfrow=c(8,4),oma=c(0,0,1,0),mar=c(5,4,2,1))
    plot(density(mysim$beta_Age),xlab=expression(beta['age']),main="",cex.lab=1.2)
    curve(dnorm(x,0,sqrt(0.05)),lwd=2,lty=2,add=T)
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)
    abline(v=0,lwd=2)
    plot(density(mysim$mu_m),main="",lwd=2,xlab=expression(mu['c']),cex.lab=1.2)
    curve(dnorm(x,36,sqrt(6 )),lwd=2,lty=2,add=T)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
plot(density(mysim$sigmasq_m),xlab=expression(sigma^2['c']),main="",lwd=2,cex.lab=1.2)
lines(density(rigamma(niterations,5,5)),lty=2,lwd=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    plot(density(mysim$omega),main="",lwd=2,xlab=expression(omega),cex.lab=1.2)
    lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2,lwd=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))

mylabels=c(paste(c("2 generations, >2 children","3 generations, 2 children","2 generations, 2 children"),"both parents",sep=", "),paste(c("2 generations, >2 children","3 generations, 2 children","2 generations, 2 children"),"1 parent",sep=", "))
names(mylabels)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")
for(i in 1:ncol(mysim$p_by_cohort)){
  sample.size= sum(data2.mat$cohort.family.type==cohort.codes[cohort.codes[,1]==i,"Cohort"])
        plot(density(mysim$p_by_cohort[,i]),xlab=expression('p'['f,c']),xlim=c(0,1),main="",lwd=2,cex.lab=1.2)
mtext(cohort.codes[cohort.codes[,1]==i,"Pop"],3,cex=0.8,padj=0.3,at=0.5,font=2,line=1.5)
mtext(mylabels[cohort.codes[cohort.codes[,1]==i,"Fam.type"]],3,cex=0.8,padj=0.3,at=0.5,font=2,line=0.5)
        curve(dbeta(x,beta.parameters[beta.parameters$Cohort.type==i,"alpha"],beta.parameters[beta.parameters$Cohort.type==i,"beta"]),add=T,lty=2,lwd=2)

        abline(v=beta.parameters[beta.parameters$Cohort.type==i,"mean.relative.to.all.informative"],lwd=1,lty=2,col="blue")
        legend("topleft",c("posterior","prior","p_observed",paste0("n = ",sample.size)),col=c("black","black","blue",NA),lty=c(1,2,2,NA),lwd=c(2,2,1,NA))
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

    mean_alpha_prior=41    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=100
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
    
if(!print.only){
   #model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.mat2,iter=niterations,chains=0)
#    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0")))
 #   model1.6.mat <- sflist2stanfit(model1.6.mat.list)
   model1.6.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0"))
 print(model1.6.mat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
} else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
   mean_alpha_prior=27
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 15
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=64
    data1.pat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
 
if(!print.only){
  #  model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.pat2,iter=niterations,chains=0)
#    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0")))
#    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
 model1.6.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.pat2,iter=niterations,chains=4,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0"))
    print(model1.6.pat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
if(myname=="maternal"){my.ylim=10}else {my.ylim=10}
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=c(-0.5,0.5),xlab="",main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
lines(density(mysim$beta_global),col="black",lwd=2,lty=4)
curve(dnorm(x,0,1),lwd=2,lty=2,add=T)
abline(v=0,lwd=2)

if(myname=="maternal"){
    polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
    abline(v=0.082,lty=4,lwd=2,col="red")
    polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
    abline(v=0.19,lty=4,lwd=2)
    abline(v=-0.42,lty=4,lwd=2,col="grey")
    abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
    legend("topright",c("Previous estimates","Kong","Coop","Hussin","Bleazard"),lty=c(NA,4,4,4,4),lwd=2,col=c(NA,"red","black","grey","darkgreen"),cex=0.8,bg="white")
   }
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"posterior","posterior global","prior global"),col=c(mycols[as.character(cohort.codes[,2])],"black","black","black"),lty=c(rep(1,nrow(cohort.codes)),1,4,2),lwd=c(rep(2,2,nrow(cohort.codes)),2,2),cex=0.8,bg="white")
dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.mu_m.pdf"),height=5,width=5)
if(myname=="maternal"){my.ylim=0.18}else {my.ylim=0.28}
for(i in 1:ncol(mysim$mu_m)){
    if(i==1){
        plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.8,bg="white")
dev.off()



    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.17
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
    lines(density(rigamma(100000,sigmasq_m_alpha,sigmasq_m_beta),n=512*8),lwd=3,lty=2)
##    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior","empirical variance in parent means"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),2,2),lwd=c(rep(2,nrow(cohort.codes)),3,1),cex=0.8,bg="white")
    dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
    lines(density(rigamma(niterations,3,0.5)),lwd=2,lty=2)
legend("topright",c("posterior","prior"),lty=c(1,2),lwd=c(2,2))
dev.off()
pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.tausq.pdf"),height=5,width=5)
    plot(density(mysim$tausq),xlab="tausq",main="Posterior for tausq",lwd=2)
    lines(density(rigamma(niterations,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

}

if(argv[1] ==33 | argv[1]==34){
#### Model 1.5 -- negative binomial model for informative families only
if(argv[1]==33){
  mu_alpha_prior = 3.7
  sigmasq_alpha_prior = 0.2
  sigmasq_m_alpha = 11
  sigmasq_m_beta = 2
  data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.mat2$sigmasq_m_beta = sigmasq_m_beta
  data1.mat2$mean_alpha_prior = mu_alpha_prior
  data1.mat2$sigmasq_alpha_prior = sigmasq_alpha_prior

  myname=paste0("maternal.mu_m_N_",mu_alpha_prior,"_",sigmasq_alpha_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta)
  if(!print.only ){
 model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age.different_link.stan",data=data1.mat2,iter=niterations,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)


print(model1.5.mat,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData"))
} else {
  load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData"))
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)

} else {

  mu_alpha_prior = 3.2
  sigmasq_alpha_prior = 0.2
  sigmasq_m_alpha = 11
  sigmasq_m_beta = 2
  myname=paste0("paternal.mu_m_N_",mu_alpha_prior,"_",sigmasq_alpha_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta)  
  data1.pat2$mean_alpha_prior = mu_alpha_prior
  data1.pat2$sigmasq_alpha_prior = sigmasq_alpha_prior
  data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.pat2$sigmasq_m_beta = sigmasq_m_beta
  
  
if(!print.only){
  model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age.different_link.stan",data=data1.pat2,iter=niterations,chains=0)
    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.pat <- sflist2stanfit(model1.5.pat.list)

print(model1.5.pat,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData"))
} else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData"))
}
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)

}

cohort.codes=cohort.codes[order(cohort.codes[,1]),]

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
my.ylim=20
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
       plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
#             plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
lines(density(mysim$beta_global),col="black",lwd=2,lty=4)
curve(dnorm(x,0,sqrt(0.05)),add=T,lwd=2,lty=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global","prior"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),4,2),
       lwd=c(rep(2,nrow(cohort.codes)),3,2),cex=0.5)
abline(v=0,lwd=2)
dev.off()
    

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
lines(density(rigamma(niterations,3,0.1)),lty=2,lwd=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
my.ylim=0.7
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mu_alpha_prior,sd=sqrt(sigmasq_alpha_prior)),lty=2,lwd=3,add=T)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()


    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
#if(myname=="maternal"){
  my.ylim=3
#    my.xlim=4
#}else{
#    my.ylim=5
#    my.xlim=4
#}
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
 #                    plot(density(mysim$sigmasq_m[,i]),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lty=2,lwd=3)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.omega.pdf"),height=5,width=5)
plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2,lwd=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
dev.off()


}


if(argv[1] ==7 | argv[1]==8){
#### Model 1.5 -- negative binomial model for informative families only
if(argv[1]==7){
  mu_alpha_prior = 3.7
  sigmasq_alpha_prior = 0.2
  sigmasq_m_alpha = 11
  sigmasq_m_beta = 2
  data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.mat2$sigmasq_m_beta = sigmasq_m_beta
  data1.mat2$mean_alpha_prior = mu_alpha_prior
  data1.mat2$sigmasq_alpha_prior = sigmasq_alpha_prior

  myname=paste0("cohort_specific_omega.maternal.mu_m_N_",mu_alpha_prior,"_",sigmasq_alpha_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta)
  if(!print.only ){
  model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age_and_omega.different_link.stan",data=data1.mat2,iter=niterations,chains=0)
  model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
  model1.5.mat <- sflist2stanfit(model1.5.mat.list)
  print(model1.5.mat,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
  save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData"))
} else {
  load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData"))
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)

} else {

  mu_alpha_prior = 3.2
  sigmasq_alpha_prior = 0.2
  sigmasq_m_alpha = 11
  sigmasq_m_beta = 2
  myname=paste0("cohort_specific_omega.paternal.mu_m_N_",mu_alpha_prior,"_",sigmasq_alpha_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta)  
  data1.pat2$mean_alpha_prior = mu_alpha_prior
  data1.pat2$sigmasq_alpha_prior = sigmasq_alpha_prior
  data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.pat2$sigmasq_m_beta = sigmasq_m_beta
  
  
if(!print.only){
  model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age_and_omega.different_link.stan",data=data1.pat2,iter=niterations,chains=0)
    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.pat <- sflist2stanfit(model1.5.pat.list)

  print(model1.5.pat,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
  save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData"))
} else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData"))
}
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)

}

cohort.codes=cohort.codes[order(cohort.codes[,1]),]

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.all_parameters.pdf"),height=17,width=4)

par(mfrow=c(5,1),oma=c(0,1,1,0),mar=c(5,5,0,1))
my.ylim=25*10
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
       plot(density(mysim$beta_Age[,i]),xlim=c(-0.05,0.05),xlab=expression(beta['age']),main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim),cex.lab=1.5)
   }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
lines(density(mysim$beta_global),col="black",lwd=2,lty=4)
curve(dnorm(x,0,sqrt(0.01)),add=T,lwd=2,lty=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global","prior"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),4,2),       lwd=c(rep(2,nrow(cohort.codes)),3,2),cex=0.8)
abline(v=0,lwd=2)

plot(density(mysim$sigmasq_global),xlab=expression(sigma^2['global']),main="",lwd=2,cex.lab=1.5)
lines(density(rigamma(niterations,3,0.1)),lty=2,lwd=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))

my.ylim=0.7*10
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i],adjust=2),xlim=range(mysim$mu_m),xlab=expression(mu['c']),main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim),cex.lab=1.5)
        }else {
            lines(density(mysim$mu_m[,i],adjust=2),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mu_alpha_prior,sd=sqrt(sigmasq_alpha_prior)),lty=2,lwd=3,add=T)
if(argv[1]==7){
    my.ylim=3*10
    my.xlim=3.5
}else{
    my.ylim=5*10
    my.xlim=4
}
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlab=expression(sigma^2['c']),main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim),cex.lab=1.5)
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lty=2,lwd=3)
if(argv[1]==7){
    my.ylim=4
    my.xlim=3.5
}else{
    my.ylim=80
    my.xlim=1.3
}
for(i in 1:ncol(mysim$omega)){
        if(i==1){
            plot(density(mysim$omega[,i]),xlim=c(1,my.xlim),xlab=expression(omega['c']),main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim),cex.lab=1.5)
        }else {
            lines(density(mysim$omega[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2,lwd=2)

dev.off()


}


if(argv[1] ==35 | argv[1]==36){
#### Model 1.5, with different link function    -- negative binomial model for informative families only, with different link function
cat("#### Model 1.5, with different link function\n")
if(argv[1]==35){
  mu_alpha_prior = 3.7
  sigmasq_alpha_prior = 0.2
  sigmasq_m_alpha = 11
  sigmasq_m_beta = 2
  data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.mat2$sigmasq_m_beta = sigmasq_m_beta
  data1.mat2$mean_alpha_prior = mu_alpha_prior
  data1.mat2$sigmasq_alpha_prior = sigmasq_alpha_prior

  myname=paste0("maternal.mu_m_N_",mu_alpha_prior,"_",sigmasq_alpha_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta)

  if(!print.only){
 #   model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age.different_link.stan",data=data1.mat2,iter=niterations,chains=0)
#    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
#    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
       model1.5.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age.different_link.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0"))
print(model1.5.mat,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
   save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.different_link.all_chains.",myname,".including_alpha.RData"))
} else {
  load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.different_link.all_chains.",myname,".including_alpha.RData")
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)

} else {

  mu_alpha_prior = 3.2
  sigmasq_alpha_prior = 0.2
  sigmasq_m_alpha = 11
  sigmasq_m_beta = 2
  data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.pat2$sigmasq_m_beta = sigmasq_m_beta
  data1.pat2$mean_alpha_prior = mu_alpha_prior
  data1.pat2$sigmasq_alpha_prior = sigmasq_alpha_prior
  
  myname=paste0("paternal.mu_m_N_",mu_alpha_prior,"_",sigmasq_alpha_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta)

 if(!print.only){
#   model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age.different_link.stan",data=data1.pat2,iter=niterations,chains=0)
#    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
#    model1.5.pat <- sflist2stanfit(model1.5.pat.list)
      model1.5.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age.different_link.stan",data=data1.pat2,iter=niterations,chains=4,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0"))
 print(model1.5.pat,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.different_link.all_chains.",myname,".including_alpha.RData"))
 } else {
   load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.different_link.all_chains.",myname,".including_alpha.RData")
}
  mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)

}
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.different_link.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
#    if(myname=="maternal"){my.ylim=200}else{my.ylim=300}
    for(i in 1:ncol(mysim$beta_Age)){
        if(i==1){
 #           plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
                     plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }else {
            lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
    abline(v=0,lwd=2)
    dev.off()
    

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.different_link.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
lines(rigamma(niterations,3,0.02),add=T,lty=2,lwd=2)
legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)
dev.off()


pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.different_link.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
#my.ylim=55
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
 #       plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
             plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
curve(rnorm(x,0,sqrt(0.01)),add=T,lwd=2,lty=3)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global","prior"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),2,3),lwd=c(rep(2,nrow(cohort.codes)),3,2),cex=0.5)
abline(v=0,lwd=2)
dev.off()
    

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.different_link.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
#      if(myname=="maternal"){my.ylim=1.6}else {my.ylim=1.6}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
 #           plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
                     plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mu_alpha_prior,sd=sqrt(sigmasq_alpha_prior)),lty=2,lwd=3,add=T)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()


    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.different_link.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
#if(myname=="maternal"){
#  my.ylim=3
#    my.xlim=4
#}else{
#    my.ylim=5
#    my.xlim=4
#}
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
 #           plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
          plot(density(mysim$sigmasq_m[,i]),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }else {
          lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
      }
lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lty=2,lwd=3)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.different_link.",myname,".posteriors.omega.pdf"),height=5,width=5)
plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2,lwd=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
dev.off()



}
if(FALSE){
if(argv[1] ==37 | argv[1]==38){
#### Model 1.5, with different link function  and weaker priors  -- negative binomial model for informative families only, with different link function
    cat("#### Model 1.5, with different link function\n")
if(argv[1]==37){
if(!print.only){
    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.weaker_priors.stan",data=data1.mat2,iter=niterations,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
 print(model1.5.mat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.different_link.weaker_priors.all_chains.mat.including_alpha.RData"))
} else {
load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.different_link.weaker_priors.all_chains.mat.including_alpha.RData")
}

    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
if(!print.only){
####Need to run this for longer - didn't converge
    model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.different_link.weaker_priors.stan",data=data1.pat2,iter=niterations,chains=0)
    model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.pat <- sflist2stanfit(model1.5.pat.list)
print(model1.5.pat,pars=c("beta_Age","beta_Cohort","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.different_link.weaker_priors.all_chains.pat.including_alpha.RData"))
} else {
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.different_link.weaker_priors.all_chains.pat.including_alpha.RData")
}
    
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.different_link.weaker_priors.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
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
    
pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.different_link.weaker_priors.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
dev.off()
}
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

    mean_alpha_prior=41    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=100
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
if(!print.only){
#    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.stan",data=data1.mat2,iter=niterations,chains=0)
#    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
#    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
      model1.6.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior, ".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_",
                mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha

     mean_alpha_prior=27
     data1.pat2$mean_alpha_prior = mean_alpha_prior
     sigmasq_m_alpha = 2
     sigmasq_m_beta = 15
     data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
     data1.pat2$sigmasq_m_beta = sigmasq_m_beta
     sigmasq_mu_m_prior=8^2
     data1.pat2$sigmasq_mu_m_prior=sigmasq_mu_m_prior
     if(!print.only){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.stan",data=data1.pat2,iter=niterations,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
#           model1.6.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.stan",data=data1.pat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
     save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
} else {
     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",
                 sigmasq_m_beta,".RData"))
 }
     mysim<-extract(model1.6.pat,permuted=T)
mysim.model1.pat = mysim
print(sum(mysim.model1.pat$beta_Age>0)/length(mysim.model1.pat$beta_Age))

     cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}
    
pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    my.ylim=100
#    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
        plot(density(mysim$beta_Age),xlim=c(-0.5,0.5),xlab="beta_Age",main="Posterior for beta_Age",lwd=2,ylim=c(0,16))
    abline(v=0,lwd=2)
    curve(dnorm(x,0,1),lty=2,lwd=2,add=T)
    if(myname=="maternal"){
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.082,lty=4,lwd=2,col="red")
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
       legend("topright",c("Previous estimates","Kong","Coop","Hussin","Bleazard"),lty=c(NA,4,4,4,4),lwd=2,col=c(NA,"red","black","grey","darkgreen"),cex=0.8,bg="white")
   }
    legend("topleft",c("posterior","prior"),lty=c(1,2),lwd=2,cex=0.8,bg="white")
    dev.off()
    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(argv[1]==43){my.ylim=0.35}else {my.ylim=0.55}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.8,bg="white")
    dev.off()

    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.13
        my.xlim=100
    }else {
        my.ylim=0.55
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(100000,sigmasq_m_alpha,sigmasq_m_beta),n=512*8),lwd=3,lty=2)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.8,bg="white")
    dev.off()

    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.tausq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),lwd=2,main="Posterior for tausq",xlab="tausq")
    lines(density(rigamma(niterations,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2,bg="white",cex=0.8)
    dev.off()
}

if(argv[1] ==9 | argv[1]==10){
#### Poisson model with lambda = a0+beta_Age*age
    df.alpha=5
    if(argv[1]==9){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    data1.mat2$df_alpha = df.alpha

    mean_alpha_prior=41    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=100
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
if(!print.only){
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.Poisson_model.stan",data=data1.mat2,iter=niterations,chains=0)
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","mu_m","sigmasq_m","exp_a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    print(model1.6.mat,pars=c("beta_Age","sigmasq_m","mu_m"),digits=4)
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.Poisson_model.linear_mean.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior, ".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.Poisson_model.linear_mean.all_chains.mat.including_alpha.mu_m_N_",
                mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="Poisson_model.linear_mean.maternal"
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha

     mean_alpha_prior=27
     data1.pat2$mean_alpha_prior = mean_alpha_prior
     sigmasq_m_alpha = 2
     sigmasq_m_beta = 15
     data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
     data1.pat2$sigmasq_m_beta = sigmasq_m_beta
     sigmasq_mu_m_prior=8^2
     data1.pat2$sigmasq_mu_m_prior=sigmasq_mu_m_prior
     if(!print.only){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.Poisson_model.stan",data=data1.pat2,iter=niterations,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","mu_m","sigmasq_m","exp_a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
    print(model1.6.pat,pars=c("beta_Age","sigmasq_m","mu_m"),digits=4)
     save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.Poisson_model.linear_mean.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
} else {
     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.Poisson_model.linear_mean.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",
                 sigmasq_m_beta,".RData"))
 }
     mysim<-extract(model1.6.pat,permuted=T)
     mysim.model1.pat = mysim
     print(sum(mysim.model1.pat$beta_Age>0)/length(mysim.model1.pat$beta_Age))

     cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
     myname="Poisson_model.linear_mean.paternal"
}
    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    my.ylim=100
    plot(density(mysim$beta_Age),xlim=c(-0.5,0.5),xlab="beta_Age",main="Posterior for beta_Age",lwd=2,ylim=c(0,16))
    abline(v=0,lwd=2)
    curve(dnorm(x,0,1),lty=2,lwd=2,add=T)
    if(myname=="maternal"){
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.082,lty=4,lwd=2,col="red")
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
       legend("topright",c("Previous estimates","Kong","Coop","Hussin","Bleazard"),lty=c(NA,4,4,4,4),lwd=2,col=c(NA,"red","black","grey","darkgreen"),cex=0.8,bg="white")
   }
    legend("topleft",c("posterior","prior"),lty=c(1,2),lwd=2,cex=0.8,bg="white")
    dev.off()
    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(argv[1]==9){my.ylim=0.35}else {my.ylim=0.55}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.8,bg="white")
    dev.off()

    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(argv[1]==9){
        my.ylim=0.13
        my.xlim=100
    }else {
        my.ylim=0.5
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.8,bg="white")
    dev.off()
}


if(argv[1] ==3 | argv[1]==4){
#### Model 1.62 -- common beta_Age for all cohorts; N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently from IG for each cohort
    df.alpha=5
    if(argv[1]==3){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family
    data1.mat2$df_alpha = df.alpha

    mean_alpha_prior=41    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=100
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
if(!print.only){
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.cohort_specific_tausq.stan",data=data1.mat2,iter=niterations,chains=0)
    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.mat <- sflist2stanfit(model1.6.mat.list)

    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.cohort_specific_tausq.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior, ".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.cohort_specific_tausq.all_chains.mat.including_alpha.mu_m_N_",
                mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="cohort_specific_tausq.maternal"
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha

     mean_alpha_prior=27
     data1.pat2$mean_alpha_prior = mean_alpha_prior
     sigmasq_m_alpha = 2
     sigmasq_m_beta = 15
     data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
     data1.pat2$sigmasq_m_beta = sigmasq_m_beta
     sigmasq_mu_m_prior=8^2
     data1.pat2$sigmasq_mu_m_prior=sigmasq_mu_m_prior
     if(!print.only){
    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.cohort_specific_tausq.stan",data=data1.pat2,iter=niterations,chains=0)
    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
    model1.6.pat <- sflist2stanfit(model1.6.pat.list)

    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
     save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.cohort_specific_tausq.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
} else {
     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.cohort_specific_tausq.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",
                 sigmasq_m_beta,".RData"))
 }
     mysim<-extract(model1.6.pat,permuted=T)
mysim.model1.pat = mysim
print(sum(mysim.model1.pat$beta_Age>0)/length(mysim.model1.pat$beta_Age))

     cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="cohort_specific_tausq.paternal"
}
    
pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    my.ylim=100
#    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
        plot(density(mysim$beta_Age),xlim=c(-0.5,0.5),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
    curve(dnorm(x,0,1),lty=2,lwd=2,add=T)
    legend("topleft",c("posterior","prior"),lty=c(1,2),lwd=2)
    if(argv[1]==3){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
       legend("topright",c("Campbell","Kong","Coop","Hussin","Bleazard"),lty=c(2,3,4,4,4),lwd=2,col=c("black","black","black","grey","darkgreen"),cex=0.5)
   }
    dev.off()
    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(argv[1]==3){my.ylim=0.4}else {my.ylim=0.7}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(argv[1]==3){
        my.ylim=0.15
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
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.tausq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    for(i in 1:ncol(mysim$tausq)){
        if(i==1){
            plot(density(mysim$tausq[,i]),xlim=range(mysim$tausq),xlab="tausq",main="Posterior for tausq",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,0.1))
        }else {
            lines(density(mysim$tausq[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(niterations,2,70)),lwd=2,lty=2)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
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

  mu_alpha_prior = 3.7
  sigmasq_alpha_prior = 2
  sigmasq_m_alpha = 11
  sigmasq_m_beta = 2
  data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.mat2$sigmasq_m_beta = sigmasq_m_beta
  data1.mat2$mean_alpha_prior = mu_alpha_prior
  data1.mat2$sigmasq_alpha_prior = sigmasq_alpha_prior

    if(!print.only){
  model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.different_link.stan",data=data1.mat2,iter=niterations,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
print(model1.5.mat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family

  mu_alpha_prior = 3.2
  sigmasq_alpha_prior = 2
  sigmasq_m_alpha = 11
  sigmasq_m_beta = 2

    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.pat2$sigmasq_m_beta = sigmasq_m_beta
  data1.pat2$mean_alpha_prior = mu_alpha_prior
  data1.pat2$sigmasq_alpha_prior = sigmasq_alpha_prior

    if(!print.only){
        model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.different_link.stan",data=data1.pat2,iter=niterations,chains=0)
        model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
        model1.5.pat <- sflist2stanfit(model1.5.pat.list)
        print(model1.5.pat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
        save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.5.pat,permuted=T)
mysim.model2.pat = mysim
    sum(mysim.model2.pat$beta_Age>0)/length(mysim.model2.pat$beta_Age)

    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}
mean_alpha_prior = mu_alpha_prior

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.",myname,".posteriors.beta_Age.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
curve(dnorm(x,0,0.05),lty=2,lwd=2,add=T)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
abline(v=0,lwd=2)
dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.",myname,".posteriors.mu_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
      if(myname=="maternal"){my.ylim=1.6}else {my.ylim=1.6}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sd=sqrt(sigmasq_alpha_prior)),lty=2,lwd=3,add=T)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()


    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.",myname,".posteriors.sigmasq_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)

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
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lty=2,lwd=3)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.",myname,".posteriors.omega.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2,lwd=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
dev.off()


}



if(argv[1] ==5 | argv[1]==6){
#### Model 1.5 -- negative binomial model for informative families only; common beta_Age for all cohorts; alpha drawn from N(mu_cohort,sigmasq_m_cohort); sigmasq_m_cohort drawn independently from IG for each cohort     
if(argv[1]==5){
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
if(!print.only){
    model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_omega.stan",data=data1.mat2,iter=niterations,chains=0)
    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
    model1.5.mat <- sflist2stanfit(model1.5.mat.list)

print(model1.5.mat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.cohort_specific_omega.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.cohort_specific_omega.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
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
    if(!print.only){
        model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_omega.stan",data=data1.pat2,iter=niterations,chains=0)
        model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
        model1.5.pat <- sflist2stanfit(model1.5.pat.list)

        print(model1.5.pat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
        save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.cohort_specific_omega.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.cohort_specific_omega.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.5.pat,permuted=T)
mysim.model2.pat = mysim
    sum(mysim.model2.pat$beta_Age>0)/length(mysim.model2.pat$beta_Age)

    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}


#pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.cohort_specific_omega.",myname,".posteriors.beta_Age.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.cohort_specific_omega.",myname,".posteriors.all_parameters.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=10,width=10)
par(mfrow=c(2,2),oma=c(0,0,1,0),mar=c(5,4,2,1))
plot(density(mysim$beta_Age),xlim=range(-0.04,0.04),xlab="",main="",lwd=2,ylim=c(0,40))
curve(dnorm(x,0,0.05),lty=2,lwd=2,add=T)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2),cex=1)
abline(v=0,lwd=2)
mtext(expression(beta['age']),1,line=2.5)
#dev.off()
#pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.cohort_specific_omega.",myname,".posteriors.mu_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
      if(myname=="maternal"){my.ylim=1.5}else {my.ylim=1.5}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="",main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sd=sqrt(6)),lty=2,lwd=3,add=T)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=1)
mtext(expression(mu['c']),1,line=2.5)
#dev.off()
 #   pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.cohort_specific_omega.",myname,".posteriors.sigmasq_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
if(myname=="maternal"){
    my.ylim=3
    my.xlim=4
}else{
    my.ylim=5
    my.xlim=4
}
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i]),xlim=c(0,my.xlim),xlab="",main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lty=2,lwd=3)
#    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=1)
mtext(expression(sigma^2['c']),1,line=2.5)
#dev.off()
#pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.cohort_specific_omega.",myname,".posteriors.omega.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
for(i in 1:ncol(mysim$omega)){
print(i )
if(argv[1]==5){
    my.ylim=4
    my.xlim=3.5
} else {
    my.ylim=80
    my.xlim=1.3
}
if(i==1){
            plot(density(mysim$omega[,i]),xlim=c(1,my.xlim),xlab="",main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
print(summary(mysim$omega[,i]))
print(mycols[as.character(cohort.codes[i,2])])
        }else {
            lines(density(mysim$omega[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
print(summary(mysim$omega[,i]))
            print(mycols[as.character(cohort.codes[i,2])])
        }
    }
lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2,lwd=2)
#legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=1)
mtext(expression(omega),1,line=2.5)
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
if(!print.only){
    if(argv[1]==47){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors.stan",data=data1.mat2,iter=niterations,chains=0)
    }
    if(argv[1]==55){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_beta_Age.stan",data=data1.mat2,iter=niterations,chains=0)
    }
    if(argv[1]==57){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_mu_m_cohort.stan",data=data1.mat2,iter=niterations,chains=0)
    }
    if(argv[1]==59){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_tau.stan",data=data1.mat2,iter=niterations,chains=0)
    }
    if(argv[1]==61){
        model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_sigma_m.stan",data=data1.mat2,iter=niterations,chains=0)
    }
#    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
#    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
    model1.6.mat = stan(fit = model1.6.mat.compiled,data = data1.mat2,chains = 4,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))

    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)

    if(argv[1]==47){
                myname="uniform_priors.maternal"
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==55){
        myname="uniform_priors_on_beta_Age.maternal"
        save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==57){
        myname="uniform_priors_on_mu_m_cohort.maternal"
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==59){
        myname="uniform_priors_on_tau.maternal"
        save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==61){
        myname="uniform_priors_on_sigma_m.maternal"
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.mat.including_alpha.RData"))
    }
} else{
    if(argv[1]==47){
        myname="uniform_priors.maternal"
             load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors.all_chains.mat.including_alpha.RData")
        }
    if(argv[1]==55){
        myname="uniform_priors_on_beta_Age.maternal"
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.mat.including_alpha.RData"))
    }
    if(argv[1]==57){
        myname="uniform_priors_on_mu_m_cohort.maternal"
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.mat.including_alpha.RData"))
    }
      if(argv[1]==59){
          myname="uniform_priors_on_tau.maternal"
          load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.mat.including_alpha.RData"))
      }
        if(argv[1]==61){
            myname="uniform_priors_on_sigma_m.maternal"
            load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.mat.including_alpha.RData"))
        }
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)

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
if(!print.only){
    if(argv[1]==48){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors.stan",data=data1.pat2,iter=niterations,chains=0)
    }
    if(argv[1]==56){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_beta_Age.stan",data=data1.pat2,iter=niterations,chains=0)
    }
    if(argv[1]==58){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_mu_m_cohort.stan",data=data1.pat2,iter=niterations,chains=0)
    }
    if(argv[1]==60){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_tau.stan",data=data1.pat2,iter=niterations,chains=0)
    }
    if(argv[1]==62){
        model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.uniform_priors_on_sigma_m.stan",data=data1.pat2,iter=niterations,chains=0)
    }

#    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
#    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
       model1.6.pat = stan(fit = model1.6.pat.compiled, data = data1.pat2,chains =4,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)

    if(argv[1]==48){
        myname="uniform_priors.paternal"
    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors.v4.all_chains.pat.including_alpha.RData"))
                                        #    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors.all_chains.pat.including_alpha.RData")
#    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors.v4.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==56){
        myname="uniform_priors_on_beta_Age.paternal"
        save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.pat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==58){
        myname="uniform_priors_on_mu_m_cohort.paternal"
        save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.pat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==60){
        save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.pat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.pat.including_alpha.RData"))
        myname="uniform_priors_on_tau.paternal"
    }
    if(argv[1]==62){
        save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.pat.including_alpha.RData"))
#        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.pat.including_alpha.RData"))
        myname="uniform_priors_on_sigma_m.paternal"
    }

} else {

   if(argv[1]==48){
       myname="uniform_priors.paternal"
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors.v4.all_chains.pat.including_alpha.RData")

    }
    if(argv[1]==56){
        myname="uniform_priors_on_beta_Age.paternal"
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_beta_Age.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==58){
        myname="uniform_priors_on_mu_m_cohort.paternal"
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_mu_m_cohort.all_chains.pat.including_alpha.RData"))
    }
    if(argv[1]==60){
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_tau.all_chains.pat.including_alpha.RData"))
        myname="uniform_priors_on_tau.paternal"
    }
    if(argv[1]==62){
        load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.uniform_priors_on_sigma_m.all_chains.pat.including_alpha.RData"))
        myname="uniform_priors_on_sigma_m.paternal"
    }

}
    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
}
#pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.uniform_priors.v4.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
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
    
        pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
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


        pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
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
        lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    }

    emp.variances=read.delim("variance_in_mean_duoHMM_count_by_parent.informative_meioses_only.min_2_kids.txt",header=T,stringsAsFactors=F)
    for(z in 1:nrow(cohort.codes)){
        if(argv[1]==47|argv[1]==55|argv[1]==57|argv[1]==59|argv[1]==61){
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"maternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        } else {
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"paternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        }
    }

    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),2),cex=0.5)
    dev.off()


        pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.tausq.pdf"),height=5,width=5)
    plot(density(mysim$tausq),main="Posterior for tausq",lwd=2,xlab="tausq")
    if(argv[1] %in% c(55,56,57,58,61,62)){
        lines(density(rigamma(niterations,2,70)),lwd=3,lty=2)
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
if(!print.only){
#   model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.uniform_priors.stan",data=data1.mat2,iter=niterations,chains=0)
#    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
#    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
    model1.5.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.uniform_priors.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0" ))
print(model1.5.mat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.uniform_priors.all_chains.mat.including_alpha.RData"))
} else {
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.uniform_priors.all_chains.mat.including_alpha.RData")
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family

    data1.pat2$mean_alpha_prior = 33
    if(!print.only){
#        model1.5.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.uniform_priors.stan",data=data1.pat2,iter=niterations,chains=0)
#        model1.5.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0")))
#      model1.5.pat <- sflist2stanfit(model1.5.pat.list)
      model1.5.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.uniform_priors.stan",data=data1.pat2,iter=niterations,chains=4,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0"))
        print(model1.5.pat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
        save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.uniform_priors.all_chains.pat.including_alpha.RData"))
    } else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.uniform_priors.all_chains.pat.including_alpha.RData"))
}
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}


pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.uniform_priors.v4.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
abline(v=0,lwd=2)
dev.off()


pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.uniform_priors.v4.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
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

    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.uniform_priors.v4.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
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
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.uniform_priors.v4.",myname,".posteriors.omega.pdf"),height=5,width=5)
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
    
if(!print.only){
#    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas.stan",data=data1.mat2,iter=niterations,chains=0)
#    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
 #   model1.6.mat <- sflist2stanfit(model1.6.mat.list)
   model1.6.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))             
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".RData"))
} else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.mat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
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
 if(!print.only){
#    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas.stan",data=data1.pat2,iter=niterations,chains=0)
#    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
#    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
       model1.6.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas.stan",data=data1.pat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distributions_on_alphas.df_",df.alpha,".RData"))
}else{
     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.pat.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distributions_on_alphas.df_",df.alpha,".RData"))
 }
    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.beta_Age.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".pdf"),height=5,width=5)
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
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.mu_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".pdf"),height=5,width=5)
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

    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.sigmasq_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".pdf"),height=5,width=5)
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
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.tausq.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas.df_",df.alpha,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),xlab="tausq",main="Posterior for tausq",lwd=2)
    lines(density(rigamma(niterations,2,70)),lwd=2,lty=2)
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

    ##mean_alpha_prior=38
    mean_alpha_prior=41
    sigmasq_mu_m_prior=100
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
if(!print.only){
#    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas_and_Y.stan",data=data1.mat2,iter=niterations,chains=0)
#    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
#    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
      model1.6.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas_and_Y.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".RData"))
}else {
load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".RData"))
}
    
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
    data1.pat2$df_alpha = df.alpha

    mean_alpha_prior=27
    sigmasq_mu_m_prior=64
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    data1.pat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior

     data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 15
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
 if(!print.only){
#    model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas_and_Y.stan",data=data1.pat2,iter=niterations,chains=0)
#    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
#    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
   model1.6.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.t_distribution_on_alphas_and_Y.stan",data=data1.pat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
         save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distributions_on_alphas_and_Y.df_",df.alpha,".RData"))
} else {
load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distributions_on_alphas_and_Y.df_",df.alpha,".RData"))
}
     mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".pdf"),height=5,width=5)
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
    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=0.6}else {my.ylim=0.7}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".pdf"),height=5,width=5)
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
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.tausq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".t_distribution_on_alphas_and_Y.df_",df.alpha,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),xlab="tausq",main="Posterior for tausq",lwd=2)
    lines(density(rigamma(niterations,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

}


if(argv[1] ==69 | argv[1]==70|argv[1] ==71 | argv[1]==72|argv[1] ==73 | argv[1]==74|argv[1] ==75 | argv[1]==76){
 ######## add in 73-76 - using data0 (only informative nuclear families)
#### Model 1.82 -- common beta_Age for all cohorts; N(mu_m,sigmasq_m) i.e. like 1.6.2 but no cohort effect - so we can compare to model 1.8 which has cohort-specific age effects but needs to avoid identifiability problem
    if(as.numeric(argv[1]) %%2 !=0){
#    cohort_by_family=rep(NA,data1.mat2$I)
#    for(i in 1:data1.mat2$I){
#        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
#    }
#    data1.mat2$cohort_by_family = cohort_by_family
    mean_alpha_prior=41    
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    sigmasq_mu_m_prior=100
    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
    data0.mat2$mean_alpha_prior = mean_alpha_prior
    data0.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data0.mat2$sigmasq_m_beta = sigmasq_m_beta
    data0.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
    if(argv[1]==69){
       myname="model1.8.2.maternal"
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
if(!print.only){       
#       model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.2.stan",data=data1.mat2,iter=niterations,chains=0)
         model1.6.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.2.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
#       model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
   }
    }
    if(argv[1]==71){
        myname="model1.8.maternal"
        cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
if(!print.only){
#            model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.stan",data=data1.mat2,iter=niterations,chains=0)
              model1.6.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0","beta_global","sigmasq_global"))
#          model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0","beta_global","sigmasq_global")))
        }

    }
    if(argv[1]==73){
        myname="model1.8.2.maternal.informative_nuc_fams_only"
        cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.informative_nuclear_families_only.more_stringent.no_GPC.txt",header=T)
if(!print.only){
#    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.2.stan",data=data0.mat2,iter=niterations,chains=0)
      model1.6.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.2.stan",data=data0.mat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
#        model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data0.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
        }
    }
    if(argv[1]==75){
        myname="model1.8.maternal.informative_nuc_fams_only"
        cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.informative_nuclear_families_only.more_stringent.no_GPC.txt",header=T)
if(!print.only){
#            model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.stan",data=data0.mat2,iter=niterations,chains=0)
              model1.6.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.stan",data=data0.mat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0","beta_global","sigmasq_global"))
#        model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data0.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0","beta_global","sigmasq_global")))
        }
   }
    if(!print.only){ 
#    model1.6.mat <- sflist2stanfit(model1.6.mat.list)
        save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.",myname,".all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior, ".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.6.mat,permuted=T)
    if("beta_global" %in% names(mysim)){
        print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
    } else {
        print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    }

} else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.",myname,".all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior, ".sigmasq_m_IG_",                sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
    mysim<-extract(model1.6.mat,permuted=T)
}

} else {
 #    cohort_by_family=rep(NA,data1.pat2$I)
 #   for(i in 1:data1.pat2$I){
 #       cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
 #   }
  #  data1.pat2$cohort_by_family = cohort_by_family
     mean_alpha_prior=27
     sigmasq_m_alpha = 2
     sigmasq_m_beta = 15
     sigmasq_mu_m_prior=8^2
     
     data1.pat2$mean_alpha_prior = mean_alpha_prior
     data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
     data1.pat2$sigmasq_m_beta = sigmasq_m_beta
     data1.pat2$sigmasq_mu_m_prior=sigmasq_mu_m_prior

     data0.pat2$mean_alpha_prior = mean_alpha_prior
     data0.pat2$sigmasq_m_alpha = sigmasq_m_alpha
     data0.pat2$sigmasq_m_beta = sigmasq_m_beta
     data0.pat2$sigmasq_mu_m_prior=sigmasq_mu_m_prior

     if(argv[1] ==70){
         myname="model1.8.2.paternal"
         cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
         if(!print.only){
 #            model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.2.stan",data=data1.pat2,iter=niterations,chains=0)
                       model1.6.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.2.stan",data=data1.pat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
 #            model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
         }
     }
     if(argv[1]==72){
         myname="model1.8.paternal"
         cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
         if(!print.only){
#             model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.stan",data=data1.pat2,iter=niterations,chains=0)
                        model1.6.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.stan",data=data1.pat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
#             model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
         }
     }
     if(argv[1] ==74){
         myname="model1.8.2.paternal.informative_nuc_fams_only"
         cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.informative_nuclear_families_only.more_stringent.no_GPC.txt",header=T)
         if(!print.only){
#             model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.2.stan",data=data0.pat2,iter=niterations,chains=0)
                        model1.6.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.2.stan",data=data0.pat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
#             model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data0.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
         }
     }
     if(argv[1]==76){
         myname="model1.8.paternal.informative_nuc_fams_only"
         cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.informative_nuclear_families_only.more_stringent.no_GPC.txt",header=T)
         if(!print.only){
#             model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.stan",data=data0.pat2,iter=niterations,chains=0)
                        model1.6.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.8.stan",data=data0.pat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
#             model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data0.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0")))
         }
     }

     if(!print.only){
#         model1.6.pat <- sflist2stanfit(model1.6.pat.list)
         save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.",myname,".all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
         mysim<-extract(model1.6.pat,permuted=T)
    if("beta_global" %in% names(mysim)){
        print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m","beta_global","sigmasq_global"),digits=4)
    } else {
        print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
    }

     } else {
         load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.",myname,".all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
         mysim<-extract(model1.6.pat,permuted=T)
     }

 }
if(FALSE){
pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/",myname,".traceplots.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    traceplot(model1.6.mat,pars=c("beta_Age"),window=c(10^2,10^4))
    traceplot(model1.6.mat,pars=c("mu_m"),window=c(10^2,10^4))
    traceplot(model1.6.mat,pars=c("sigmasq_m"),window=c(10^2,10^4))
    traceplot(model1.6.mat,pars=c("tausq"),window=c(10^2,10^4))
    if("beta_global" %in% names(mysim)){
        traceplot(model1.6.mat,pars=c("beta_global"),window=c(10^2,10^4))
        traceplot(model1.6.mat,pars=c("sigmasq_global"),window=c(10^2,10^4))
    }
    dev.off()
}
    
pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    my.ylim=100
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
       abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
       legend("topright",c("Campbell","Kong","Coop","Hussin","Bleazard"),lty=c(2,3,4,4,4),lwd=2,col=c("black","black","black","grey","darkgreen"),cex=0.5)
    }
    dev.off()
    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$mu_m),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",lwd=2)
    curve(dnorm(x,mean_alpha_prior,sigmasq_mu_m_prior),lwd=2,lty=2,add=T)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()


        pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$sigmasq_m),xlab="sigmasq_m",main="Posterior for sigmasq_m",lwd=2)
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/",myname,".posteriors.tausq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),lwd=2,main="Posterior for tausq",xlab="tausq")
    lines(density(rigamma(niterations,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)
    dev.off()
}


####model 3.2

if(argv[1] %in% c(77,78)){

if(argv[1] ==77){
cat("Model 3.2 maternal, with adjusted priors\n")
if(!print.only){
data2.mat2$mu_alpha_prior= 3.7
data2.mat2$sigmasq_alpha_prior = 0.2
    model3.2.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.different_link.stan",data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
model3.2.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.2.mat.compiled, data = data2.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
model3.2.mat <- sflist2stanfit(model3.2.mat.list)

save(model3.2.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.2.adjusted_priors.different_link.all_chains.mat.including_alpha.RData"))
print(model3.2.mat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
} else {
load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.2.adjusted_priors.different_link.all_chains.mat.including_alpha.RData")
print(model3.2.mat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
}

mysim<-extract(model3.2.mat,permuted=T)
parent="adjusted_priors.maternal"
cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.more_stringent.no_GPC.txt",header=T,stringsAsFactors=F)
beta.parameters=read.delim("parameters_for_beta_distribution.model3_maternal.more_stringent.no_GPC.txt",header=T)[,1:5]

}

if(argv[1] ==78){
    cat("Model 3.2 paternal, with adjusted priors\n")
if(!print.only){
data2.pat2$mu_alpha_prior= 3.2
data2.pat2$sigmasq_alpha_prior = 0.2
    model3.2.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.adjusted_priors.different_link.stan",data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
    model3.2.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.2.pat.compiled, data = data2.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
    model3.2.pat <- sflist2stanfit(model3.2.pat.list)

    save(model3.2.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.2.adjusted_priors.different_link.all_chains.pat.including_alpha.RData"))
    print(model3.2.pat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
} else {
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.2.adjusted_priors.different_link.all_chains.pat.including_alpha.RData")
}
    mysim<-extract(model3.2.pat,permuted=T)
    parent="adjusted_priors.paternal"
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.more_stringent.no_GPC.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_paternal.more_stringent.no_GPC.txt",header=T)[,1:5]
}
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
 pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model3.2.different_link.",parent,".posteriors.all_parameters.pdf"),height=8*3,width=4*3)
par(mfrow=c(8,4))

    plot(density(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age")
    curve(dnorm(x,0,sqrt(0.05)),lwd=2,lty=2,add=T)
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)
    abline(v=0,lwd=2)

plot(density(mysim$sigmasq_m),xlab="sigmasq_m",main="Posterior for sigmasq_m",lwd=2)
lines(density(rigamma(niterations,10,2)),lty=2,lwd=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))

plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
    lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2,lwd=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))

plot(density(mysim$mu_m),main="Posterior for mu_m",lwd=2,xlab="mu_m")
    curve(dnorm(x,3.7,sqrt(0.2)),lwd=2,lty=2,add=T)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))

mylabels=c(paste(c("informative, 2 generations","informative, 3 generations","uninformative, 2 kids"),"both parents",sep=", "),paste(c("informative, 2 generations","informative, 3 generations","uninformative, 2 kids"),"1 parent",sep=", "))
names(mylabels)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")

for(i in 1:ncol(mysim$p_by_cohort)){
        sample.size= sum(data2.mat$cohort.family.type==cohort.codes[cohort.codes[,1]==i,"Cohort"])
        plot(density(mysim$p_by_cohort[,i]),xlab="p",xlim=c(0,1),main=paste0(cohort.codes[cohort.codes[,1]==i,"Pop"],", ",mylabels[cohort.codes[cohort.codes[,1]==i,"Fam.type"]]),lwd=2)
        curve(dbeta(x,beta.parameters[beta.parameters$Cohort.type==i,"alpha"],beta.parameters[beta.parameters$Cohort.type==i,"beta"]),add=T,lty=2,lwd=2)
        abline(v=beta.parameters[beta.parameters$Cohort.type==i,"mean.relative.to.all.informative"],lwd=1,lty=2,col="blue")
        legend("topleft",c("posterior","prior","p_observed",paste0("n = ",sample.size)),col=c("black","black","blue",NA),lty=c(1,2,2,NA),lwd=c(2,2,1,NA))
    }
    dev.off()


}

if(argv[1] ==7){
#### Model 1.6, but with NTR split into pill, no pill and not sure
if(argv[1]==7){

  mycols=c("black","blue","red","green","orange","skyblue","tan3","red4","purple","mediumseagreen","magenta","darkgreen","violetred2")
names(mycols)=c("CARL","FC","FVG","GPC","NTR_no_ART","NTR_ART","NTR_maybe_ART","QTR370_maybe_ART","QTR610_no_ART","QTR610_ART","QTR610_maybe_ART","VB","ORCADES")

  
  load("duoHMM_data_for_RSTAN.more_stringent.no_GPC.ART_separated.RData")
  mean_alpha_prior=41    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=100
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
    
if(!print.only){
    model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data2.pat2,iter=10,chains=0)
      model1.6.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit =  model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0")))
      model1.6.mat <- sflist2stanfit( model1.6.mat.list)  

 print(model1.6.mat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.ART_separated.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
} else {
  load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.ART_separated.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.ART_separated.txt",header=T)
    myname="maternal"
}


pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.ART_separated.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.all_parameters.pdf"),height=20,width=4)
par(mfrow=c(5,1),oma=c(0,0,1,0),mar=c(5,4,2,1))

my.ylim=5.5
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=c(-0.5,0.5),xlab=expression(beta['age']),main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim),cex.lab=1.2)
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
lines(density(mysim$beta_global),col="black",lwd=2,lty=4)
curve(dnorm(x,0,1),lwd=2,lty=2,add=T)
abline(v=0,lwd=2)
if(myname=="maternal"){
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.082,lty=4,lwd=2,col="red")
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
       legend("topright",c("Previous estimates","Kong","Coop","Hussin","Bleazard"),lty=c(NA,4,4,4,4),lwd=2,col=c(NA,"red","black","grey","darkgreen"),cex=0.8,bg="white")
   }
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global","prior"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),4,2),lwd=c(rep(2,nrow(cohort.codes)),2,2),cex=0.8,bg="white")
  plot(density(mysim$sigmasq_global),xlab=expression(sigma^2['global']),main="",lwd=2,cex.lab=1.2)
    lines(density(rigamma(niterations,3,0.5)),lwd=2,lty=2)
legend("topright",c("posterior","prior"),lty=c(1,2),lwd=c(2,2),cex=0.8)
my.ylim=0.2
for(i in 1:ncol(mysim$mu_m)){
    if(i==1){
        plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab=expression(mu['c']),main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim),cex.lab=1.2)
    }else {
        lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
    if(myname=="maternal"){
        my.ylim=0.2
        my.xlim=100
    }else {
        my.ylim=0.6
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i],n=256),xlim=c(0,my.xlim),xlab=expression(sigma^2['c']),main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim),cex.lab=1.2)
        }else {
            lines(density(mysim$sigmasq_m[,i],n=256),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    plot(density(mysim$tausq),xlab=expression(tau^2),main="",lwd=2,cex.lab=1.2)
    lines(density(rigamma(niterations,2,70)),lwd=2,lty=2)
    dev.off()

}



if(argv[1] ==8){
#### Model 1.6, but with NTR split into pill, no pill and not sure
if(argv[1]==8){

  mycols=c("black","blue","red","green","orange","skyblue","tan3","red4","purple","darkgreen","violetred2")
names(mycols)=c("CARL","FC","FVG","GPC","NTR_no_pill","NTR_pill","NTR_maybe_pill","QTR370","QTR610","VB","ORCADES")

  load("duoHMM_data_for_RSTAN.more_stringent.no_GPC.NTR_pill_separated.RData")
      mean_alpha_prior=41    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=100
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
    
if(!print.only){
   model1.6.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0"))
 print(model1.6.mat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.NTR_pill_separated.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
} else {
  load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.NTR_pill_separated.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.NTR_pill_separated.txt",header=T)
    myname="maternal"
}

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.NTR_pill_separated.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.all_parameters.pdf"),height=20,width=4)
par(mfrow=c(5,1),oma=c(0,0,1,0),mar=c(5,4,2,1))

my.ylim=5.5
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=c(-0.5,0.5),xlab=expression(beta['age']),main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim),cex.lab=1.2)
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
lines(density(mysim$beta_global),col="black",lwd=2,lty=4)
curve(dnorm(x,0,1),lwd=2,lty=2,add=T)
abline(v=0,lwd=2)
if(myname=="maternal"){
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.082,lty=4,lwd=2,col="red")
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
       legend("topright",c("Previous estimates","Kong","Coop","Hussin","Bleazard"),lty=c(NA,4,4,4,4),lwd=2,col=c(NA,"red","black","grey","darkgreen"),cex=0.8,bg="white")
   }
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global","prior"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),4,2),lwd=c(rep(2,nrow(cohort.codes)),2,2),cex=0.8,bg="white")
  plot(density(mysim$sigmasq_global),xlab=expression(sigma^2['global']),main="",lwd=2,cex.lab=1.2)
    lines(density(rigamma(niterations,3,0.5)),lwd=2,lty=2)
legend("topright",c("posterior","prior"),lty=c(1,2),lwd=c(2,2),cex=0.8)
my.ylim=0.2
for(i in 1:ncol(mysim$mu_m)){
    if(i==1){
        plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab=expression(mu['c']),main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim),cex.lab=1.2)
    }else {
        lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
    if(myname=="maternal"){
        my.ylim=0.2
        my.xlim=100
    }else {
        my.ylim=0.6
        my.xlim=30
    }
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
            plot(density(mysim$sigmasq_m[,i],n=256),xlim=c(0,my.xlim),xlab=expression(sigma^2['c']),main="",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim),cex.lab=1.2)
        }else {
            lines(density(mysim$sigmasq_m[,i],n=256),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    plot(density(mysim$tausq),xlab=expression(tau^2),main="",lwd=2,cex.lab=1.2)
    lines(density(rigamma(niterations,2,70)),lwd=2,lty=2)
    dev.off()

}

if(FALSE){

if(argv[1] ==80){
#### Model 1.5 -- negative binomial model for informative families only
if(argv[1]==80){
  mycols=c("black","blue","red","green","orange","skyblue","tan3","red4","purple","darkgreen","violetred2")
names(mycols)=c("CARL","FC","FVG","GPC","NTR_no_pill","NTR_pill","NTR_maybe_pill","QTR370","QTR610","VB","ORCADES")
  load("duoHMM_data_for_RSTAN.more_stringent.no_GPC.NTR_pill_separated.RData")

  mu_alpha_prior = 37
  sigmasq_alpha_prior = 6
  sigmasq_m_alpha = 5
  sigmasq_m_beta = 5
  data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.mat2$sigmasq_m_beta = sigmasq_m_beta
  data1.mat2$mean_alpha_prior = mu_alpha_prior
  data1.mat2$sigmasq_alpha_prior = sigmasq_alpha_prior

  myname=paste0("NTR_pill_separated.maternal.mu_m_N_",mu_alpha_prior,"_",sigmasq_alpha_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta)
  if(!print.only ){
  model1.5.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0"))
print(model1.5.mat,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData"))
} else {
  load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData"))
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.NTR_pill_separated.txt",header=T)
}
cohort.codes=cohort.codes[order(cohort.codes[,1]),]

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
my.ylim=20
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
       plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
#             plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=4)
curve(dnorm(x,0,sqrt(0.05)),add=T,lwd=2,lty=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global","prior"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),2,3),lwd=c(rep(2,nrow(cohort.codes)),4,2),cex=0.5)
abline(v=0,lwd=2)
dev.off()
    
pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
lines(density(rigamma(niterations,3,0.1)),lty=2,lwd=2)
dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
my.ylim=0.7
    for(i in 1:ncol(mysim$mu_m)){
              if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mu_alpha_prior,sd=sqrt(sigmasq_alpha_prior)),lty=2,lwd=3,add=T)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){

                     plot(density(mysim$sigmasq_m[,i]),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,3))
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lty=2,lwd=3)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.omega.pdf"),height=5,width=5)
plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2,lwd=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
dev.off()
}






if(argv[1] ==81 | argv[1]==82){
#### Model 1.6 -- beta_Age drawn from distribution; alphas drawn from N(mu_cohort,sigmasq)
if(argv[1]==81){
    cohort_by_family=rep(NA,data0.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family

    mean_alpha_prior=41    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=100
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
    
if(!print.only){
   #model1.6.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.mat2,iter=niterations,chains=0)
#    model1.6.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0")))
 #   model1.6.mat <- sflist2stanfit(model1.6.mat.list)
   model1.6.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0"))
 print(model1.6.mat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
} else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.all_chains.mat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="maternal"
} else {
    cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family
   mean_alpha_prior=27
    data1.pat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 15
    data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.pat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=64
    data1.pat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
 
if(!print.only){
  #  model1.6.pat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.pat2,iter=niterations,chains=0)
#    model1.6.pat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.6.pat.compiled, data = data1.pat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0")))
#    model1.6.pat <- sflist2stanfit(model1.6.pat.list)
 model1.6.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.stan",data=data1.pat2,iter=niterations,chains=4,pars=c("beta_Age","beta_global","sigmasq_global","tausq","mu_m","sigmasq_m","a0"))
    print(model1.6.pat,pars=c("beta_Age","beta_global","sigmasq_global","tausq","sigmasq_m","mu_m"),digits=4)
    save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.all_chains.pat.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)
    myname="paternal"
}

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
if(myname=="maternal"){my.ylim=7}else {my.ylim=14}
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
        plot(density(mysim$beta_Age[,i]),xlim=c(-0.5,0.5),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
lines(density(mysim$beta_global),col="black",lwd=2,lty=4)
curve(dnorm(x,0,1),lwd=2,lty=2,add=T)
abline(v=0,lwd=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global","prior"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),4,2),lwd=c(rep(2,nrow(cohort.codes)),2,2),cex=0.5)
if(myname=="maternal"){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
       legend("topright",c("Campbell","Kong","Coop","Hussin","Bleazard"),lty=c(2,3,4,4,4),lwd=2,col=c("black","black","black","grey","darkgreen"),cex=0.5)
   }
    dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.mu_m.pdf"),height=5,width=5)
if(myname=="maternal"){my.ylim=0.25}else {my.ylim=0.5}
for(i in 1:ncol(mysim$mu_m)){
    if(i==1){
        plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
    }else {
        lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()



    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
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
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    emp.variances=read.delim("variance_in_mean_duoHMM_count_by_parent.informative_meioses_only.min_2_kids.txt",header=T,stringsAsFactors=F)
    for(z in 1:nrow(cohort.codes)){
        if(argv[1]==81){
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"maternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        } else {
            abline(v=emp.variances[as.character(cohort.codes[z,2]),"paternal"],lty=2,col=mycols[as.character(cohort.codes[z,2])])
        }
    }
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior","empirical variance in parent means"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),2,2),lwd=c(rep(2,nrow(cohort.codes)),3,1),cex=0.5)
    dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
    lines(density(rigamma(niterations,3,0.5)),lwd=2,lty=2)
legend("topright",c("posterior","prior"),lty=c(1,2),lwd=c(2,2))
dev.off()
pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".",myname,".posteriors.tausq.pdf"),height=5,width=5)
    plot(density(mysim$tausq),xlab="tausq",main="Posterior for tausq",lwd=2)
    lines(density(rigamma(niterations,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
    dev.off()

}

if(argv[1] ==33 | argv[1]==34){
#### Model 1.5 -- negative binomial model for informative families only
if(argv[1]==33){
  mu_alpha_prior = 37
  sigmasq_alpha_prior = 6
  sigmasq_m_alpha = 5
  sigmasq_m_beta = 5
  data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.mat2$sigmasq_m_beta = sigmasq_m_beta
  data1.mat2$mean_alpha_prior = mu_alpha_prior
  data1.mat2$sigmasq_alpha_prior = sigmasq_alpha_prior

  myname=paste0("maternal.mu_m_N_",mu_alpha_prior,"_",sigmasq_alpha_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta)
  if(!print.only ){
 # model1.5.mat.compiled = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age.stan",data=data1.mat2,iter=niterations,chains=0)
#    model1.5.mat.list = mclapply(1:4, mc.cores = 4,function(i) stan(fit = model1.5.mat.compiled, data = data1.mat2,chains = 1, chain_id = i,iter=niterations,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0")))
#    model1.5.mat <- sflist2stanfit(model1.5.mat.list)
  model1.5.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0"))

print(model1.5.mat,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData"))
} else {
  load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData")
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)

} else {

  mu_alpha_prior = 32
  sigmasq_alpha_prior = 6
  sigmasq_m_alpha = 5
  sigmasq_m_beta = 5
  myname=paste0("paternal.mu_m_N_",mu_alpha_prior,"_",sigmasq_alpha_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta)  
  data1.pat2$mean_alpha_prior = mu_alpha_prior
  data1.pat2$sigmasq_alpha_prior = sigmasq_alpha_prior
  data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
  data1.pat2$sigmasq_m_beta = sigmasq_m_beta
  
  
if(!print.only){
  model1.5.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.cohort_specific_beta_Age.stan",data=data1.pat2,iter=niterations,chains=4,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m","exp_a0"))
print(model1.5.pat,pars=c("beta_Age","beta_global","sigmasq_global","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData"))
} else {
    load("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.all_chains.",myname,".including_alpha.RData")
}
    mysim<-extract(model1.5.pat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.more_stringent.no_GPC.txt",header=T)

}

cohort.codes=cohort.codes[order(cohort.codes[,1]),]

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.beta_Age.pdf"),height=5,width=5)
#my.ylim=55
for(i in 1:ncol(mysim$beta_Age)){
    if(i==1){
             plot(density(mysim$beta_Age[,i]),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }else {
        lines(density(mysim$beta_Age[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
    }
}
lines(density(mysim$beta_global),col="black",lwd=3,lty=2)
curve(rnorm(x,0,sqrt(0.05)),add=T,lwd=2,lty=3)
legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"global","prior"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),2,3),lwd=c(rep(2,nrow(cohort.codes)),3,2),cex=0.5)
abline(v=0,lwd=2)
dev.off()
    

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.sigmasq_global.pdf"),height=5,width=5)
  plot(density(mysim$sigmasq_global),xlab="sigmasq_global",main="Posterior for sigmasq_global",lwd=2)
lines(rigamma(niterations,3,0.1),add=T,lty=2,lwd=2)
dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.mu_m.pdf"),height=5,width=5)
      if(myname=="maternal"){my.ylim=1.6}else {my.ylim=1.6}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mu_alpha_prior,sd=sqrt(sigmasq_alpha_prior)),lty=2,lwd=3,add=T)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()


    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.sigmasq_m.pdf"),height=5,width=5)
    for(i in 1:ncol(mysim$sigmasq_m)){
        if(i==1){
                     plot(density(mysim$sigmasq_m[,i]),xlab="sigmasq_m",main="Posterior for sigmasq_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }else {
            lines(density(mysim$sigmasq_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lty=2,lwd=3)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.",myname,".posteriors.omega.pdf"),height=5,width=5)
plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2,lwd=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
dev.off()


}
}
if(argv[1] %in% c(79:86)){

if(argv[1] %in% c(79,81,83,85)){
cat("Model 3.2 maternal, with adjusted priors\n")

  mu_alpha_prior = 3.7
  sigmasq_alpha_prior = 0.2
  data2.mat2$mu_alpha_prior = mu_alpha_prior
  data2.mat2$sigmasq_alpha_prior = sigmasq_alpha_prior


if(argv[1] == 79){
  stan.code  = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.uniform_prior_on_beta_Age.different_link.stan"
  name="uniform_prior_on_beta_Age.maternal"
}

if(argv[1] == 81){
  stan.code  = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.uniform_prior_on_sigma_m.different_link.stan"
  name="uniform_prior_on_sigma_m.maternal"
}

if(argv[1] == 83){

  stan.code  = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.uniform_prior_on_mu_m.different_link.stan"
  name="uniform_prior_on_mu_m.maternal"
}


if(argv[1] == 85){
  stan.code  = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.uniform_prior_on_p.different_link.stan"
  name="uniform_prior_on_p.maternal"
}

if(!print.only){
model3.2.mat.compiled = stan(file = stan.code,data=data2.mat2,iter=10,chains=0,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
model3.2.mat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.2.mat.compiled, data = data2.mat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
model3.2.mat <- sflist2stanfit(model3.2.mat.list)
save(model3.2.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.2.",name,".RData"))
print(model3.2.mat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
} else {
load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.2.",name,".RData"))
print(model3.2.mat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
}
mysim<-extract(model3.2.mat,permuted=T)
parent="adjusted_priors.maternal"
cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_3.more_stringent.no_GPC.txt",header=T,stringsAsFactors=F)
beta.parameters=read.delim("parameters_for_beta_distribution.model3_maternal.more_stringent.txt",header=T)[,1:5]
}else {
    cat("Model 3.2 paternal, with adjusted priors\n")
  mu_alpha_prior = 3.2
  sigmasq_alpha_prior = 0.2
  data2.pat2$mu_alpha_prior = mu_alpha_prior
  data2.pat2$sigmasq_alpha_prior = sigmasq_alpha_prior

if(argv[1] == 80){
  stan.code  = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.uniform_prior_on_beta_Age.different_link.stan"
  name="uniform_prior_on_beta_Age.paternal"

}

if(argv[1] == 82){
  stan.code  = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.uniform_prior_on_sigma_m.different_link.stan"
  name="uniform_prior_on_sigma_m.paternal"
}


if(argv[1] == 84){

  stan.code  = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.uniform_prior_on_mu_m.different_link.stan"
  name="uniform_prior_on_mu_m.paternal"
}


if(argv[1] == 86){
  stan.code  = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model3.2.uniform_prior_on_p.different_link.stan"
  name="uniform_prior_on_p.paternal"
}

if(!print.only){
    model3.2.pat.compiled = stan(file = stan.code,data=data2.pat2,iter=10,chains=0,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0"))
    model3.2.pat.list <- mclapply(1:4, mc.cores = 4,function(i) stan(fit = model3.2.pat.compiled, data = data2.pat2,chains = 1, chain_id = i,iter=10000,pars=c("beta_Age","p_by_cohort","omega","mu_m","sigmasq_m","exp_a0")))
    model3.2.pat <- sflist2stanfit(model3.2.pat.list)
    save(model3.2.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.2.",name,".RData"))
    print(model3.2.pat,pars=c("beta_Age","p_by_cohort","mu_m","sigmasq_m","omega"),digits=4)
} else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.model3.2.",name,".RData"))
}
    mysim<-extract(model3.2.pat,permuted=T)
    parent="adjusted_priors.paternal"
    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_3.more_stringent.no_GPC.txt",header=T,stringsAsFactors=F)
    beta.parameters=read.delim("parameters_for_beta_distribution.model3_paternal.more_stringent.txt",header=T)[,1:5]
}
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



 pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model3.2.",name,".posteriors.all_parameters.pdf"),height=7*3,width=4*3)
par(mfrow=c(7,4),oma=c(0,0,3,0),mar=c(5,4,2,1))

    plot(density(mysim$beta_Age),xlab=expression(beta['age']),main="Posterior for beta_Age",cex.lab=1.2)
if(!argv[1] %in% c(79:80)){
  curve(dnorm(x,0,sqrt(0.05)),lwd=2,lty=2,add=T)
}
   legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)
   abline(v=0,lwd=2)

  #    plot(density(mysim$mu_m),main="Posterior for mu_m",lwd=2,xlab="mu_m")
    plot(density(mysim$mu_m),main="Posterior for mu_m",lwd=2,xlab=expression(mu['c']),cex.lab=1.2)
      if(!argv[1] %in% c(83:84)){
        curve(dnorm(x,36,sqrt(6 )),lwd=2,lty=2,add=T)
 }
   legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))


#   plot(density(mysim$sigmasq_m),xlab="sigmasq_m",main="Posterior for sigmasq_m",lwd=2)
   plot(density(mysim$sigmasq_m),xlab=expression(sigma^2['c']),main="Posterior for sigmasq_m",lwd=2,cex.lab=1.2)
   if(!argv[1] %in% c(81:82)){
     lines(density(rigamma(10000,5,5)),lty=2,lwd=2)
   }
      legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
      
#      plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
      plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab=expression(omega),cex.lab=1.2)
      lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2,lwd=2)
      legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
      
mylabels=c(paste(c("2 generations, >2 children","3 generations, 2 children","2 generations, 2 children"),"both parents",sep=", "),paste(c("2 generations, >2 children","3 generations, 2 children","2 generations, 2 children"),"1 parent",sep=", "))
   names(mylabels)=c("infor.2gen.2parents","infor.3gen.2parents","noninfor.2kids.2gen.2parents","infor.2gen.1parent","infor.3gen.1parent","noninfor.2kids.2gen.1parent")

   
   for(i in 1:ncol(mysim$p_by_cohort)){
     sample.size= sum(data2.mat$cohort.family.type==cohort.codes[cohort.codes[,1]==i,"Cohort"])
        plot(density(mysim$p_by_cohort[,i]),xlab=expression('p'['f,c']),xlim=c(0,1),main="",lwd=2,cex.lab=1.2)
mtext(cohort.codes[cohort.codes[,1]==i,"Pop"],3,cex=0.8,padj=0.3,at=0.5,font=2,line=1.5)
mtext(mylabels[cohort.codes[cohort.codes[,1]==i,"Fam.type"]],3,cex=0.8,padj=0.3,at=0.5,font=2,line=0.5)

      if(!argv[1] %in% c(85:86)){
        curve(dbeta(x,beta.parameters[beta.parameters$Cohort.type==i,"alpha"],beta.parameters[beta.parameters$Cohort.type==i,"beta"]),add=T,lty=2,lwd=2)
}
        abline(v=beta.parameters[beta.parameters$Cohort.type==i,"mean.relative.to.all.informative"],lwd=1,lty=2,col="blue")
        legend("topleft",c("posterior","prior","p_observed",paste0("n = ",sample.size)),col=c("black","black","blue",NA),lty=c(1,2,2,NA),lwd=c(2,2,1,NA))
    }


if(argv[1] %in% 79:80){
figname="Supplementary Figure 10A"
} 
if(argv[1] %in% 81:82){
figname="Supplementary Figure 10C"
} 

if(argv[1] %in% 83:84){
figname="Supplementary Figure 10B"
} 

if(argv[1] %in% 85:86){
figname="Supplementary Figure 10D"
} 

mtext(figname,3,outer=T,cex=1.5,padj=0.2,at=0.2,font=2,line=1.5)

dev.off()


}




if(argv[1] ==87 | argv[1]==88){
#### Model 1.62 using informative nuclear families only
data1.mat2=data0.mat2
data1.pat2=data0.pat2
    if(argv[1]==87){
    cohort_by_family=rep(NA,data1.mat2$I)
    for(i in 1:data1.mat2$I){
        cohort_by_family[i] <- data1.mat2$cohort[data1.mat2$family ==i][1]
    }
    data1.mat2$cohort_by_family = cohort_by_family


    mean_alpha_prior=41    
    data1.mat2$mean_alpha_prior = mean_alpha_prior
    sigmasq_m_alpha = 2
    sigmasq_m_beta = 40
    data1.mat2$sigmasq_m_alpha = sigmasq_m_alpha
    data1.mat2$sigmasq_m_beta = sigmasq_m_beta
    sigmasq_mu_m_prior=100
    data1.mat2$sigmasq_mu_m_prior = sigmasq_mu_m_prior
if(!print.only){
      model1.6.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
    print(model1.6.mat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
       save(model1.6.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.mat.informative_nuclear_families.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior, ".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.mat.informative_nuclear_families.including_alpha.mu_m_N_",
                mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.6.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.informative_nuclear_families_only.more_stringent.no_GPC.txt",header=T)
    myname="maternal.informative_nuclear_families"
} else {
     cohort_by_family=rep(NA,data1.pat2$I)
    for(i in 1:data1.pat2$I){
        cohort_by_family[i] <- data1.pat2$cohort[data1.pat2$family ==i][1]
    }
    data1.pat2$cohort_by_family = cohort_by_family


     mean_alpha_prior=27
     data1.pat2$mean_alpha_prior = mean_alpha_prior
     sigmasq_m_alpha = 2
     sigmasq_m_beta = 15
     data1.pat2$sigmasq_m_alpha = sigmasq_m_alpha
     data1.pat2$sigmasq_m_beta = sigmasq_m_beta
     sigmasq_mu_m_prior=8^2
     data1.pat2$sigmasq_mu_m_prior=sigmasq_mu_m_prior
     if(!print.only){
           model1.6.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.6.2.stan",data=data1.pat2,iter=niterations,chains=4,pars=c("beta_Age","tausq","mu_m","sigmasq_m","a0"))
    print(model1.6.pat,pars=c("beta_Age","tausq","sigmasq_m","mu_m"),digits=4)
     save(model1.6.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.pat.informative_nuclear_families.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
} else {
     load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.6.2.all_chains.pat.informative_nuclear_families.including_alpha.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",
                 sigmasq_m_alpha,"_",                 sigmasq_m_beta,".RData"))
 }
     mysim<-extract(model1.6.pat,permuted=T)
mysim.model1.pat = mysim
print(sum(mysim.model1.pat$beta_Age>0)/length(mysim.model1.pat$beta_Age))

     cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.informative_nuclear_families_only.more_stringent.no_GPC.txt",header=T)
    myname="paternal.informative_nuclear_families"
}
    
pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.beta_Age.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    my.ylim=100
#    plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
        plot(density(mysim$beta_Age),xlim=c(-0.5,0.5),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
    abline(v=0,lwd=2)
    curve(dnorm(x,0,1),lty=2,lwd=2,add=T)
    legend("topleft",c("posterior","prior"),lty=c(1,2),lwd=2)
    if(myname=="maternal"){
       polygon(x=rep(c(0.067-0.0215,0.067+0.0215),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("grey",0.2),border=NA)
       abline(v=0.067,lty=2,lwd=2)
       polygon(x=rep(c(0.082-0.012,0.082+0.012),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("lightblue",0.2),border=NA)
       abline(v=0.082,lty=3,lwd=2)
       polygon(x=rep(c(0.19-0.092,0.19+0.092),2)[c(1,2,4,3)],y=c(0,0,my.ylim,my.ylim),col=alpha("pink",0.2),border=NA)
       abline(v=0.19,lty=4,lwd=2)
       abline(v=-0.42,lty=4,lwd=2,col="grey")
       abline(v=-0.29,lty=4,lwd=3,col="darkgreen")
       legend("topright",c("Campbell","Kong","Coop","Hussin","Bleazard"),lty=c(2,3,4,4,4),lwd=2,col=c("black","black","black","grey","darkgreen"),cex=0.5)
   }
    dev.off()
    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.mu_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(myname=="maternal"){my.ylim=0.4}else {my.ylim=0.7}
    for(i in 1:ncol(mysim$mu_m)){
        if(i==1){
            plot(density(mysim$mu_m[,i]),xlim=range(mysim$mu_m),xlab="mu_m",main="Posterior for mu_m",col=mycols[as.character(cohort.codes[i,2])],lwd=2,ylim=c(0,my.ylim))
        }else {
            lines(density(mysim$mu_m[,i]),col=mycols[as.character(cohort.codes[i,2])],lwd=2)
        }
    }
    curve(dnorm(x,mean_alpha_prior,sqrt(sigmasq_mu_m_prior)),add=T,lwd=3,lty=2)
    legend("topleft",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
    dev.off()

    
    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.sigmasq_m.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    if(myname=="maternal"){
        my.ylim=0.15
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
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lwd=3,lty=2)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior","empirical variance in parent means"),col=c(mycols[as.character(cohort.codes[,2])],"black","black"),lty=c(rep(1,nrow(cohort.codes)),2,2),lwd=c(rep(2,nrow(cohort.codes)),3,1),cex=0.5)
    dev.off()

    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.6.2.",myname,".posteriors.tausq.mu_m_N_",mean_alpha_prior,"_",sigmasq_mu_m_prior,".sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
    plot(density(mysim$tausq),lwd=2,main="Posterior for tausq",xlab="tausq")
    lines(density(rigamma(niterations,2,70)),lwd=2,lty=2)
    legend("topright",c("posterior","prior"),lty=c(1,2),lwd=2)
    dev.off()
}

    



if(argv[1] ==89 | argv[1]==90){
data1.mat2=data0.mat2
data1.pat2=data0.pat2

#### Model 1.5.4 on informative nuclear families only
if(argv[1]==89){
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
if(!print.only){
      model1.5.mat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.stan",data=data1.mat2,iter=niterations,chains=4,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0"))
print(model1.5.mat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
    save(model1.5.mat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.all_chains.mat.informative_nuclear_families.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else{
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.all_chains.mat.informative_nuclear_families.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.5.mat,permuted=T)
    cohort.codes=read.delim("RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.informative_nuclear_families_only.more_stringent.no_GPC.txt",header=T)
    myname="maternal.informative_nuclear_families"
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
    if(!print.only){
              model1.5.pat = stan(file = "/well/donnelly/hilary/maternal_age_and_recombination/bin/RSTAN.model1.5.4.stan",data=data1.pat2,iter=niterations,chains=4,pars=c("beta_Age","omega","mu_m","sigmasq_m","exp_a0"))
        print(model1.5.pat,pars=c("beta_Age","omega","mu_m","sigmasq_m"),digits=4)
        save(model1.5.pat,file=paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.all_chains.pat.informative_nuclear_families.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}else {
    load(paste0("/well/donnelly/hilary/maternal_age_and_recombination/RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/RSTAN.modell1.5.4.all_chains.pat.informative_nuclear_families.including_alpha.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".RData"))
}
    mysim<-extract(model1.5.pat,permuted=T)
mysim.model2.pat = mysim
    sum(mysim.model2.pat$beta_Age>0)/length(mysim.model2.pat$beta_Age)

    cohort.codes=read.delim("RSTAN_output/key_for_paternal_cohorts_to_include_in_model_1.informative_nuclear_families_only.more_stringent.no_GPC.txt",header=T)
    myname="paternal.informative_nuclear_families"
}


pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.",myname,".posteriors.beta_Age.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
plot(density(mysim$beta_Age),xlim=range(mysim$beta_Age),xlab="beta_Age",main="Posterior for beta_Age",lwd=2)
curve(dnorm(x,0,0.05),lty=2,lwd=2,add=T)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
abline(v=0,lwd=2)
dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.",myname,".posteriors.mu_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
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


    pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.",myname,".posteriors.sigmasq_m.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)

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
    lines(density(rigamma(niterations,sigmasq_m_alpha,sigmasq_m_beta)),lty=2,lwd=3)
    legend("topright",c(as.character(cohort.codes[order(cohort.codes[,1]),2]),"prior"),col=c(mycols[as.character(cohort.codes[,2])],"black"),lty=c(rep(1,nrow(cohort.codes)),2),lwd=c(rep(2,nrow(cohort.codes)),3),cex=0.5)
dev.off()

pdf(paste0("RSTAN_output_on_duoHMM_more_stringent/no_GPC_different_link_for_NB/model1.5.4.",myname,".posteriors.omega.sigmasq_m_IG_",sigmasq_m_alpha,"_",sigmasq_m_beta,".pdf"),height=5,width=5)
plot(density(mysim$omega),main="Posterior for omega",lwd=2,xlab="omega")
lines(1/seq(from=0,to=1,by=0.01),seq(from=0,to=1,by=0.01),type="l",lty=2,lwd=2)
legend("topright",c("posterior","prior"),lwd=2,lty=c(1,2))
dev.off()


}


