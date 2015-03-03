setwd("/well/donnelly/hilary/maternal_age_and_recombination/")
load("duoHMM_data_for_RSTAN.more_stringent.with_new_QTR.RData")
duohmm.data1.mat=data1.mat # inform meioses

duohmm.data1.mat.simple=duohmm.data1.mat[,c("CHILD","PARENT","nrec","age.at.birth","cohort")
colnames(duohmm.data1.mat.simple)[3] = "n.crossovers"

write.table(duohmm.data1.mat.simple,"maternal_crossovers.informative_meioses.with_new_QTR.txt",quote=F,sep="\t",row.names=F)



setwd("/well/donnelly/hilary/maternal_age_and_recombination/ART")
counts=read.delim("../maternal_crossovers.informative_meioses.with_new_QTR.txt",header=T,stringsAsFactors=F)

ntr.art=read.delim("pregnancy_info_NTR.txt",header=T,na.string="#NULL!")
qtr.new.art = read.delim("QTR_new_ART_data.txt",header=T,na.string="")
qtr.old.art=read.delim("QTR_original_ART_data.txt",header=T,na.string="")
qtr.final.art=read.delim("ALL_QTR_ART_scott_2015.txt",header=T,na.string="")[,1:29]
qtr.final.art$mumid = as.character(    qtr.final.art$mumid)
indiv.convert = read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/data_120215/IndivID_DeidentTable_Feb2015version_Rel7.txt",header=F,stringsAsFactors=F,sep=" ")
fam.convert = read.delim("/well/donnelly/hilary/maternal_age_and_recombination/QTR/data_120215/FamilyID_DeidentTable_Feb2015version_Rel7.txt",header=F,stringsAsFactors=F,sep=" ")
colnames(indiv.convert) = c("ident.IID","deident.IID")
colnames(fam.convert) = c("ident.FID","deident.FID")
    
qtr.final.art=cbind(qtr.final.art,indiv.convert[match(qtr.final.art$mumid,indiv.convert$ident.IID),])
qtr.final.art$any.response = rowSums(!is.na(qtr.final.art[,2:29]))>0
qtr.final.art$any.yes = rowSums(qtr.final.art[,2:29]==1,na.rm=T) >0
qtr.final.art$any.no = rowSums(qtr.final.art[,2:29]==2,na.rm=T) >0


#    1 = yes
#    2 = no

    
ntr.counts=counts[counts$cohort=="NTR",]
qtr.counts=counts[counts$cohort %in% c("QTR370","QTR610","QTRCoreExome"),]

ntr.counts$modeofpregn = ntr.art[match(ntr.counts$PARENT,ntr.art$GID),"modeofpregn"]
ntr.counts$preIVF = ntr.art[match(ntr.counts$PARENT,ntr.art$GID),"preIVF"]


ntr.counts.non.art=ntr.counts[!is.na(ntr.counts$modeofpregn) & ((ntr.counts$modeofpregn==1|(!is.na(ntr.counts$preIVF) & ntr.counts$preIVF==1))),]
ntr.counts.art=ntr.counts[!is.na(ntr.counts$modeofpregn) & (ntr.counts$modeofpregn==2),]
ntr.counts.unknown.art = ntr.counts[is.na(ntr.counts$modeofpregn)|!is.na(ntr.counts$modeofpregn) & ntr.counts$modeofpregn %in% c(3,4)|!is.na(ntr.counts$modeofpregn) & ntr.counts$modeofpregn  == 9 & is.na(ntr.counts$preIVF) ,]

                                        #NTR
#The code for the variable mode of pregn is :
#1: spontaneous
#2:non-spontaneous (hormones, ivf, etc,.)
#3:inconsistent
#4:two twin pairs, one spontaneous, one not
#9: no information or not relevant
#For the preIVF:
#    1: yes
#ever_pill: ever used birth control pills, as reported in ANTR survey data, with 0=no, 1=yes
#nreport_pill : number of reports on birth control pill use in ANTR survey data.
#age_max_report : maximum age of report on birth control pills in ANTR survey data.
#ac_pill_pr : birth control pill prior to twin pregnancy, based on data from several projects, with 0=no, 1=yes.
#time_pill : when was pill use stopped prior to twin pregnancy, based on data from several projects, with 1= stopped less than 2 months prior to pregnancy, 2= stopped 2-6 months prior to pregnancy, 3= stopped more than 6 months prior to pregnancy, 4=I do not know.
#% IVF = nrow(ntr.counts.art)/nrow(ntr.counts.non.art))


qtr.counts=cbind(qtr.counts,qtr.final.art[match(qtr.counts$PARENT,qtr.final.art$deident.IID),31:34])
qtr.counts=cbind(qtr.counts,qtr.new.art[match(qtr.counts$PARENT,qtr.new.art$MUM_ID),],qtr.old.art[match(qtr.counts$PARENT,qtr.old.art$Individual.ID),])

qtr.counts$old.any.ART=(rowSums(qtr.counts[,c("X12_IVF1","X12_IVF2","X12_IVF3","X14_IVF1", "X14_IVF2","X14_IVF3","X16_IVF1","X16_IVF2","X16_IVF3")]==1,na.rm=T)!=0)
qtr.counts$old.answered.any.ART=(rowSums(!is.na(qtr.counts[,c("X12_IVF1","X12_IVF2","X12_IVF3","X14_IVF1", "X14_IVF2","X14_IVF3","X16_IVF1","X16_IVF2","X16_IVF3")]))!=0)
qtr.counts$old.any.ART[!qtr.counts$old.answered.any.ART] = NA

qtr.counts$old.online.ART=(rowSums(qtr.counts[,c("HT","IVF","NT")]==1,na.rm=T)!=0)
qtr.counts$old.answered.online.ART=(rowSums(!is.na(qtr.counts[,c("HT","IVF","NT")]))!=0)
#if they answered other questions in the women's health section but didn't say "yes" to these, we'll assume the answer is 'no'
qtr.counts$old.answered.online.ART[qtr.counts$QSTAT==1] = TRUE
qtr.counts$old.online.ART[!qtr.counts$old.answered.online.ART] = NA

#check for inconsistencies
table(qtr.counts[,c("old.any.ART","old.online.ART","ANY_ART")])

    
qtr.counts.non.art=qtr.counts[(qtr.counts$ANY_ART=="n" & !is.na(qtr.counts$ANY_ART)) |
                                  (!qtr.counts$old.online.ART & qtr.counts$old.answered.online.ART) |
                                  (!qtr.counts$old.any.ART & qtr.counts$old.answered.any.ART) |
                                      (qtr.counts$kids.preIVF==1 & !is.na(qtr.counts$kids.preIVF))|
                                         (!qtr.counts$any.yes & qtr.counts$any.no & !is.na(qtr.counts$any.response)),]


qtr.counts.art=qtr.counts[(!qtr.counts$ANY_ART=="n" & !is.na(qtr.counts$ANY_ART)) |
                                  (qtr.counts$old.online.ART & qtr.counts$old.answered.online.ART) |
                                  (qtr.counts$old.any.ART & qtr.counts$old.answered.any.ART)|
                                     (qtr.counts$any.yes &!is.na(qtr.counts$any.response)),]



conflicting = qtr.counts.non.art[qtr.counts.non.art$old.any.ART & !qtr.counts.non.art$old.online.ART & !is.na(qtr.counts.non.art$old.any.ART ) & !is.na(qtr.counts.non.art$old.online.ART),]
conflicting2 = qtr.counts.non.art[!is.na(qtr.counts.non.art$kids.preIVF) & (qtr.counts.non.art$kids.preIVF ==1 & ((!is.na(qtr.counts.non.art$old.any.ART) & qtr.counts.non.art$old.any.ART) |(!is.na(qtr.counts.non.art$old.online.ART) &  qtr.counts.non.art$old.online.ART)
                                                                                | (!is.na(qtr.counts.non.art$ANY_ART) & qtr.counts.non.art$ANY_ART=="y" ))),]
conflicting=rbind(conflicting,conflicting2)

qtr.counts.non.art = qtr.counts.non.art[!qtr.counts.non.art$PARENT %in% conflicting$PARENT,]
qtr.counts.non.art = qtr.counts.non.art[!qtr.counts.non.art$PARENT %in% qtr.counts.art$PARENT,]

qtr.counts.unknown.art = qtr.counts[!qtr.counts$CHILD %in% c(qtr.counts.non.art$CHILD,qtr.counts.art$CHILD),]


load("../duoHMM_data_for_RSTAN.more_stringent.with_new_QTR.no_GPC.RData")

data1.mat$cohort2 = data1.mat$cohort

data1.mat$cohort2[data1.mat$PARENT %in% ntr.counts.non.art$PARENT]="NTR_no_ART"
data1.mat$cohort2[data1.mat$PARENT %in% ntr.counts.art$PARENT]="NTR_ART"
data1.mat$cohort2[data1.mat$PARENT %in% ntr.counts.unknown.art$PARENT]="NTR_maybe_ART"

data1.mat$cohort2[data1.mat$cohort=="QTR610" & data1.mat$PARENT %in% qtr.counts.non.art$PARENT]="QTR610_no_ART"
data1.mat$cohort2[data1.mat$cohort=="QTR610" & data1.mat$PARENT %in% qtr.counts.art$PARENT]="QTR610_ART"
data1.mat$cohort2[data1.mat$cohort=="QTR610" & data1.mat$PARENT %in% qtr.counts.unknown.art$PARENT]="QTR610_maybe_ART"

data1.mat$cohort2[data1.mat$cohort=="QTR370" & data1.mat$PARENT %in% qtr.counts.non.art$PARENT]="QTR370_no_ART"
data1.mat$cohort2[data1.mat$cohort=="QTR370" & data1.mat$PARENT %in% qtr.counts.art$PARENT]="QTR370_ART"
data1.mat$cohort2[data1.mat$cohort=="QTR370" & data1.mat$PARENT %in% qtr.counts.unknown.art$PARENT]="QTR370_maybe_ART"

data1.mat$cohort2[data1.mat$cohort=="QTRCoreExome" & data1.mat$PARENT %in% qtr.counts.non.art$PARENT]="QTRCoreExome_no_ART"
data1.mat$cohort2[data1.mat$cohort=="QTRCoreExome" & data1.mat$PARENT %in% qtr.counts.art$PARENT]="QTRCoreExome_ART"
data1.mat$cohort2[data1.mat$cohort=="QTRCoreExome" & data1.mat$PARENT %in% qtr.counts.unknown.art$PARENT]="QTRCoreExome_maybe_ART"


table(data1.mat$cohort2)
data1.mat=data1.mat[data1.mat$cohort2!="QTR370_ART",]
    
    
    
mat.key.model1=unique(data.frame(as.integer(as.factor(data1.mat$cohort2)),as.factor(data1.mat$cohort2)))
colnames(mat.key.model1)=c("Code","Cohort")
write.table(mat.key.model1,"../RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.with_new_QTR.no_GPC.ART_separated.txt",quote=F,sep="\t",row.names=F)


data1.mat2=list(y=data1.mat$nrec,cohort=as.integer(as.factor(data1.mat$cohort2)),family=as.integer(as.factor(data1.mat$PARENT)),Age=as.numeric(data1.mat$age.at.birth),J=nrow(data1.mat),I=length(unique(data1.mat$PARENT)),    C=length(unique(data1.mat$cohort2)))
cohort_by_family=rep(NA,data1.mat2$I)
for(i in 1:data1.mat2$I){
  cohort_by_family[i]= data1.mat2$cohort[data1.mat2$family ==i][1]
}
data1.mat2$cohort_by_family = cohort_by_family
save.image("../duoHMM_data_for_RSTAN.more_stringent.with_new_QTR.no_GPC.ART_separated.RData")

############ IVF

ntr.ivf=read.delim("NTR_pill_use_03122014.txt",header=T,stringsAsFactors=T,colClasses=c("character","integer","integer","integer","real","integer","integer"))
#need to work out whether the twins were born first

ntr.ages=read.delim("../NTR/NTR_v2_age_data.txt",header=T,colClasses=c("character","integer","integer"))
ntr.fam=read.delim("../NTR/HM2_1.fam",header=F,colClasses=c("character","character","character","character","integer","integer"),sep="")
ntr.ages$mother=ntr.fam[match(ntr.ages[,1],ntr.fam[,2]),4]

ntr.ivf$twins.last = unlist(sapply(ntr.ivf$GID_M6,function(mother){
    returned=0
    if(sum(ntr.ages$mother==mother)==0){# mother has no children in this dataset
        return(NA)
        returned=1
    } else {
    fam.ages=ntr.ages[ntr.ages$mother==mother,]
    age.counts=table(fam.ages$age_mum_birth)
    if(sum(is.na(fam.ages$age_mum_birth))>0){
        return(NA) #missing mum ages, don't trust
        returned=1
    }else if(length(age.counts)==1 & age.counts>1){# just a single pair of twins/triplets, no siblings
        return(1) #twins first
        returned=1
    } else if(length(age.counts)==1 & age.counts<2){ #only one child
        return(NA) # not clear who the twins are
        returned=1
    }else if(length(age.counts)>1){
        twin.ages = as.numeric(names(age.counts[age.counts>1]))
        sib.ages =  as.numeric(names(age.counts[age.counts==1]))
        print(mother)
#        if(length(twin.ages)==1){
#            cat(twin.ages,"\t",sib.ages,"\n")
#        }
        if(length(twin.ages)==0){
            return(NA) #only sibs, no twins
            returned=1
        } else if(length(twin.ages)==1 & twin.ages > max(sib.ages)){ #only one pair of twins/triplets, and these were born after all the sibs
            return(1)
            returned=1
        }else if(length(twin.ages)==1 & twin.ages < max(sib.ages)){ #only one pair of twins/triplets, and these were born before at least one of the sibs
            return(0)
            returned=1
        }else if(length(twin.ages)>0){ #more than one pair of twins/triplets - ignore
            return(NA)
            returned=1
        }
    }
    if(returned==0){
        print(fam.ages)
    }
}
}))

# GID_M6 gender ever_pill nreport_pill age_max_report ac_pill_pr time_pill

ntr.counts=cbind(ntr.counts,ntr.ivf[match(ntr.counts$PARENT,ntr.ivf$GID_M6),])

ntr.no.pill=ntr.counts[(!is.na(ntr.counts$ever_pill) & ntr.counts$ever_pill == 0 ) | (!is.na(ntr.counts$ac_pill_pr) & ntr.counts$ac_pill_pr == 0 & !is.na(ntr.counts$twins.last ) & ntr.counts$twins.last ==1),]

ntr.pill=ntr.counts[ (!is.na(ntr.counts$ac_pill_pr) & ntr.counts$ac_pill_pr == 1 & ntr.counts$twins.last ==1 & !is.na(ntr.counts$twins.last)),]

ntr.pill.not.sure =ntr.counts[! ntr.counts$PARENT %in% c(ntr.no.pill$PARENT,ntr.pill$PARENT),]


write.table(ntr.no.pill,"maternal_crossovers.informative_meioses.with_twins_indicated.NTR_only_non_pill_mothers.txt",quote=F,sep="\t",row.names=F)



load("../duoHMM_data_for_RSTAN.more_stringent.with_new_QTR.no_GPC.RData")

data1.mat$cohort2 = data1.mat$cohort
data1.mat$cohort2[data1.mat$PARENT %in% ntr.no.pill$PARENT]="NTR_no_pill"
data1.mat$cohort2[data1.mat$PARENT %in% ntr.pill$PARENT]="NTR_pill"
data1.mat$cohort2[data1.mat$PARENT %in% ntr.pill.not.sure$PARENT]="NTR_maybe_pill"

mat.key.model1=unique(data.frame(as.integer(as.factor(data1.mat$cohort2)),as.factor(data1.mat$cohort2)))
colnames(mat.key.model1)=c("Code","Cohort")
write.table(mat.key.model1,"../RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.with_new_QTR.no_GPC.NTR_pill_separated.txt",quote=F,sep="\t",row.names=F)


data1.mat2=list(y=data1.mat$nrec,cohort=as.integer(as.factor(data1.mat$cohort2)),family=as.integer(as.factor(data1.mat$PARENT)),Age=as.numeric(data1.mat$age.at.birth),J=nrow(data1.mat),I=length(unique(data1.mat$PARENT)),    C=length(unique(data1.mat$cohort2)))
cohort_by_family=rep(NA,data1.mat2$I)
for(i in 1:data1.mat2$I){
  cohort_by_family[i]= data1.mat2$cohort[data1.mat2$family ==i][1]
}
data1.mat2$cohort_by_family = cohort_by_family
save.image("../duoHMM_data_for_RSTAN.more_stringent.with_new_QTR.no_GPC.NTR_pill_separated.RData")

##### Now add on data for ORCADES pill

load("../duoHMM_data_for_RSTAN.more_stringent.with_new_QTR.no_GPC.NTR_pill_separated.RData")
orcades.pill=read.delim("/well/donnelly/hilary/maternal_age_and_recombination/ART/ORCADES_pill_data_Oct2014.txt",header=T,stringsAsFactors=F)

orcades.never.pill = orcades.pill[orcades.pill$PillEver ==2,]
orcades.ever.pill = orcades.pill[orcades.pill$PillEver ==1,]
orcades.unknown.pill = orcades.pill[orcades.pill$PillEver ==9,]

data1.mat$cohort2[data1.mat$PARENT %in% orcades.never.pill[,1]]="ORCADES_no_pill"
data1.mat$cohort2[data1.mat$PARENT %in% orcades.ever.pill[,1]]="ORCADES_pill"
data1.mat$cohort2[data1.mat$PARENT %in% orcades.unknown.pill[,1]]="ORCADES_maybe_pill"

mat.key.model1=unique(data.frame(as.integer(as.factor(data1.mat$cohort2)),as.factor(data1.mat$cohort2)))
colnames(mat.key.model1)=c("Code","Cohort")
write.table(mat.key.model1,"../RSTAN_output/key_for_maternal_cohorts_to_include_in_model_1.more_stringent.with_new_QTR.no_GPC.NTR_and_ORCADES_pill_separated.txt",quote=F,sep="\t",row.names=F)

data1.mat2=list(y=data1.mat$nrec,cohort=as.integer(as.factor(data1.mat$cohort2)),family=as.integer(as.factor(data1.mat$PARENT)),Age=as.numeric(data1.mat$age.at.birth),J=nrow(data1.mat),I=length(unique(data1.mat$PARENT)),    C=length(unique(data1.mat$cohort2)))
cohort_by_family=rep(NA,data1.mat2$I)
for(i in 1:data1.mat2$I){
  cohort_by_family[i]= data1.mat2$cohort[data1.mat2$family ==i][1]
}
data1.mat2$cohort_by_family = cohort_by_family
save.image("../duoHMM_data_for_RSTAN.more_stringent.with_new_QTR.no_GPC.NTR_and_ORCADES_pill_separated.RData")


twin.datasets.no.art = rbind(ntr.counts.non.art[,c("CHILD","PARENT","n.crossovers","age.at.birth","cohort")],qtr.counts.non.art[,c("CHILD","PARENT","n.crossovers","age.at.birth","cohort")])
write.table(twin.datasets.no.art,"maternal_crossovers.with_new_QTR.informative_meioses.with_twins_indicated.QTR_and_NTR_only_nonART_mothers.txt",quote=F,sep="\t",row.names=F)


#QTR old
#Q_STATDescription
#0:did not complete online Q
#1:answered other questions in women's health section
#2:completed some of questionnaire, none of women's health section

#Columns C-K
#Source 1: Clinic visit
#The field names are the same for all 3 visits (prefix 12/14/16) and correspond to the following questions.
#Field Name  Question
#IVF1                 Did any of your pregnancies (including those that may have ended in miscarriage or stillbirth), involve hormone treatment prescribed by your doctor?
#IVF2                Was IVF involved in the conception of your twins?
#IVF3                Did any of your pregnancies (including those that may have ended in miscarriage or stillbirth), involve alternative or natural fertility treatments?
#For all fields 1=yes, 2=no. Dots represent no data found.

#Columns L-N
#Source2: Online Questionnaire
#Question – If you have ever been pregnant, for any of the pregnancies did you require or experience any of the following…
#(Check boxes used by respondent to indicate positive response)
#Field Name  Response
#HT                   Hormone Treatment
#IVF                   IVF
#NT                   Natural Therapy Treatment
#For all fields 1=yes, dot represents option not chosen by respondent or no data found

#QTR new
#1.     MUM_DOB = the twins’ mother’s DOB.... whereas MOTH_DOB = the twin mother’s mother’s DOB (as does MOTH_AGE)
#2.    DAD_DOB = the twins’ father’s DOB.... whereas FATH_DOB = the twin mother’s father’s DOB (as does FATH_AGE)
#3.    Where there has been no answer given as to DOB or age of parent, both fields have a ‘?’. If data has been given for one of those fields, the other one is just left blank
#4.    y = yes, n = no
#5.    only where the mum has answered ‘y’ to ANY_ART is there further information recorded
#6.    if any information about ART has been given for the twins it is listed against twin 1’s ID (TWIN_ID)

