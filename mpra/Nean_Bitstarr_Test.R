library(data.table)
library(qqman)
library(ggplot2)
library(tidyverse)

datFolder="/rs/rs_grp_neanbit/"

filenames <- c("VitD"="Individual_Bitstarr_VitD.txt",
               "EtOH"="Individual_Bitstarr_EtOH.txt",
               "Water"="Individual_Bitstarr_Water.txt",
               "VitA"="Individual_Bitstarr_VitA.txt",
               "Selenium"="Individual_Bitstarr_Selenium.txt",
               "Caffeine"="Individual_Bitstarr_Caffeine.txt",
               "Dex"="Individual_Bitstarr_Dex.txt")


mydat <- map_dfr(names(filenames),function(ii){
  fname=paste0(datFolder,filenames[ii])
  aa= fread(fname)
  aa$Treatment=ii
  aa
})


mydat$batch = gsub("LCL.*","",mydat$.id)
table(mydat$batch,mydat$Treatment)

mydat$Tr=mydat$Treatment

mydat$Tr[mydat$Tr=="EtOH"]="Control"
mydat$Tr[mydat$Tr=="Water"]="Control"

mydat$Tr = factor(mydat$Tr) %>% relevel("Control")


mydat <- mutate(mydat, betacorr =  betas.beta.binom - qlogis(DNA_prop))

myc = mydat %>% group_by(identifier) %>% summarize(n=n(),k=length(unique(Tr)),c=sum(Tr=="Control"))

mycfilt <- myc %>% filter(n>3,k>1,c>2)

mydat2 <- mydat %>% 
  filter(identifier %in% mycfilt$identifier) %>%
  select(identifier, R, A, DNA_prop, bin, Treatment, Tr, batch) %>%
  mutate(N=R+A) %>% 
  group_by(identifier) %>%
  mutate(Nstar=mean(N))

msum <- mydat2 %>% summarise(Nstar=mean(N)) 
Nstar <- msum$Nstar


nbreaks=10
cov_breaks <- unique(c(0,quantile(Nstar,(1:nbreaks)/nbreaks)))
bin <- cut(Nstar,cov_breaks)
M <- exp((0:500)/50)
msum$bin <- bin
table(bin)

msum2 <- msum %>% group_by(bin) %>% summarize(N=mean(Nstar))


mydat2$bin <- cut(mydat2$Nstar,cov_breaks)

levels(mydat2$bin)

eps=0.001

Mvec <- sapply(levels(bin),function(mybin){
  cat("Estimating Bin",mybin,":")
  aux <- sapply(M,function(M){
    sum(logLikBetaBinomialRhoEps(mydat2$DNA_prop[mydat2$bin==mybin],
                                 eps,M,
                                 mydat2$R[mydat2$bin==mybin ],
                                 mydat2$A[mydat2$bin==mybin]))
  })
  Mmax <- M[which.max(aux)]
  cat(" ",which.max(aux)," ",Mmax,"\n")
##  Mmax
  aux/sum(aux)
})

plot(log10(M),Mvec[,1])

plot(log10(msum2$N),Mvec)

mybin=levels(bin)[10]
aux <- sapply(M,function(M){
  sum(logLikBetaBinomialRhoEps(mydat2$DNA_prop[mydat2$bin==mybin],
                               eps,M,
                               mydat2$R[mydat2$bin==mybin ],
                               mydat2$A[mydat2$bin==mybin]))
})

plot(log10(M),aux)

paux <- exp(aux-max(aux))
paux <- paux/sum(paux)
plot(log10(M),log10(paux))

which.max(aux)


head(mydat2)

dim(Mvec)



# Function to implement linear model and include DF
mylm <- function(data){
  ##  SNP_Individual=paste0(data$rsID[1],"_",data$Individual[1])
  obj <- lm(betacorr ~ batch + Tr, data = data, weights = 1/betas_se^2)
  mysum <- summary(obj)
  mycoeff <- as.data.frame(mysum$coefficients)
  colnames(mycoeff) <- c("estimate","std.error","statistic","p.value")
  mycoeff %>% rownames_to_column("term") %>%
    mutate(deg.f=obj$df.residual,
           rank=obj$qr$rank)#,
  #kappa=kappa.qr(obj$qr))
}

aux <- mydat %>% group_by(identifier) %>% nest()
res <- aux %>% mutate(res=map(data,mylm))

##res2 <- aux %>% mutate(res=map_dfr(data,mylm))

res2 <- res %>% select(identifier=identifier,data=res) %>% unnest(cols=c(data))
rescInt <- res2 %>% filter(term=="(Intercept)") %>% ungroup %>%  mutate(padj=p.adjust(p.value))
#Make qq plots and save dfs
resDex <- res2 %>% filter(term=="TrDex")
resDex$padj<-p.adjust(resDex$p.value)
resVitD <- res2 %>% filter(term=="TrVitD")
resVitD$padj<-p.adjust(resVitD$p.value)
resVitA <- res2 %>% filter(term=="TrVitA")
resVitA$padj<-p.adjust(resVitA$p.value)
resCaffeine <- res2 %>% filter(term=="TrCaffeine")
resCaffeine$padj<-p.adjust(resCaffeine$p.value)
resSelenium <- res2 %>% filter(term=="TrSelenium")
resSelenium$padj<-p.adjust(resSelenium$p.value)


png(file="/nfs/rprscratch/wwwShare/carly/ASE_Plots/dex_case_manually_calc_qq_linear_model.png")
qq(resDex$p.value)
dev.off()
png(file="/nfs/rprscratch/wwwShare/carly/ASE_Plots/vitd_case_manually_calc_qq_linear_model.png")
qq(resVitD$p.value)
dev.off()
png(file="/nfs/rprscratch/wwwShare/carly/ASE_Plots/vita_case_manually_calc_qq_linear_model.png")
qq(resVitA$p.value)
dev.off()
png(file="/nfs/rprscratch/wwwShare/carly/ASE_Plots/caf_case_manually_calc_qq_linear_model.png")
qq(resCaffeine$p.value)
dev.off()
png(file="/nfs/rprscratch/wwwShare/carly/ASE_Plots/selenium_case_manually_calc_qq_linear_model.png")
qq(resSelenium$p.value)
dev.off()
png(file="/nfs/rprscratch/wwwShare/carly/ASE_Plots/ase_manually_calc_qq_linear_model.png")
qq(rescInt$p.value)
dev.off()

write.table(resDex,file="case_linearmodel_Bitstarr_dex.txt",quote=FALSE,row.names=FALSE,sep="\t")
write.table(resSelenium,file="case_linearmodel_Bitstarr_selenium.txt",quote=FALSE,row.names=FALSE,sep="\t")
write.table(resCaffeine,file="case_linearmodel_Bitstarr_caffeine.txt",quote=FALSE,row.names=FALSE,sep="\t")
write.table(resVitD,file="case_linearmodel_Bitstarr_vitd.txt",quote=FALSE,row.names=FALSE,sep="\t")
write.table(resVitA,file="case_linearmodel_Bitstarr_vita.txt",quote=FALSE,row.names=FALSE,sep="\t")
write.table(rescInt,file="ase_linearmodel_Bitstarr.txt",quote=FALSE,row.names=FALSE,sep="\t")