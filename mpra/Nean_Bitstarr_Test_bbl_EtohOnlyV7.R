library(data.table)
library(qqman)
library(ggplot2)
library(tidyverse)
library(QuASAR2)

datFolder="/rs/rs_grp_neanbit/"

filenames <- c("VitD"="Individual_Bitstarr_VitD.txt",
               "EtOH"="Individual_Bitstarr_EtOH.txt",
               ##"Water"="Individual_Bitstarr_Water.txt",
               "VitA"="Individual_Bitstarr_VitA.txt",
               ##"Selenium"="Individual_Bitstarr_Selenium.txt",
               ##"Caffeine"="Individual_Bitstarr_Caffeine.txt",
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
##mydat$Tr[mydat$Tr=="Water"]="Control"
mydat$Tr = factor(mydat$Tr) %>% relevel("Control")
mydat$batch<- factor(mydat$batch)

##mydat <- mutate(mydat, betacorr =  betas.beta.binom - qlogis(DNA_prop))

nbreaks = 20; M = 20; eps = 0.001;

fitQuasar2 <- function(dd, nbreaks = 20, M = 20, eps = 0.001) {
  dd <- mydat2 %>%  group_by(identifier) %>%
  mutate(MN=mean(R+A))
  MN <- dd %>% summarize(MN=mean(R+A)) 
  cov_breaks <- unique(c(0,quantile(MN$MN,(1:nbreaks)/nbreaks)))
  ##bin <- cut(MN$MN,cov_breaks)
  ##table(bin)
  dd <- dd %>% mutate(bin=cut(MN,cov_breaks)) 
  binlab <- levels(dd$bin)
  
  neweps <- eps
  
  dd$M <- M
  Mvec <- rep(M,length(binlab))
  names(Mvec)<- binlab
  
  aux <- dd %>% group_by(identifier) %>% nest()
  res <- aux %>% mutate(res=map(data,fitBetaBinomialLogistic,~ batch + Tr, eps=neweps))


  res <- res %>% mutate(data=map(res, "dd"),res=map(res,"res"))

  res2 <- res %>% select(identifier,res) %>% unnest(cols=c(res))
  res3 <- res %>% select(identifier, data) %>% unnest(cols=c(data))
  
  iteration=1
  notConverged=TRUE
  while(notConverged){
    
    
  cat("Iteration: ",iteration,"\n")
  cat("Optimizing M\n")
  Mres <- sapply(binlab,function(bin){
    selv <- res3$bin==bin;
    auxLogis <- optim(par=Mvec[bin],
                    fn=logLikBetaBinomialM,
                    gr=gLogLikBetaBinomialM,
                    p=res3$p.fit[selv],
                    R=res3$R[selv],
                    A=res3$A[selv],
                    method="L-BFGS-B",#control=c(trace=6),
                    hessian=TRUE,
                    lower=0.1,
                    upper=10000)
    data.frame(M=auxLogis$par,Mse=1/(auxLogis$hessian)^.5,conv=auxLogis$convergence)
  })
  Mres <- t(Mres)
  OldMvec <- Mvec
  Mvec <- unlist(Mres[,1])
  cat("M change max:",max(abs(Mvec-OldMvec)),"\n")
  
  
  oldeps=neweps
  auxLogis <- optim(par=neweps,
                    fn=logLikBetaBinomialEps,
                    gr=gLogLikBetaBinomialEps,
                    p=res3$p.fit,
                    R=res3$R,
                    A=res3$A,
                    M=res3$M,
                    method="L-BFGS-B",#control=c(trace=6),
                    hessian=TRUE,
                    lower=1E-10,
                    upper=0.1)
  auxLogis$convergence
  auxLogis$message
  c(auxLogis$par,1/(auxLogis$hessian)^.5)
  neweps=auxLogis$par
  cat("eps change max:",max(abs(neweps-oldeps)),"\n")
  
##  res3 <- dd %>% group_by(identifier) %>% nest()

  res3$M <- Mvec[res3$bin]
  
  aux <- res3 %>% group_by(identifier) %>% nest()

  aux$res <- res$res
  
  res <- aux %>% mutate(res=map2(data,res, ~fitBetaBinomialLogistic(.x,model=~ batch + Tr, b=.y$coeff,eps=neweps))) 

  res <- res %>% mutate(data=map(res, "dd"),res=map(res,"res"))
  
  res3 <- res %>% select(identifier, data) %>% 
    unnest(cols=c(data))
  
  oldcoeff <- res2$coeff
  res2 <- res %>% select(identifier,res) %>% unnest(cols=c(res))
  
  cat("Max coeff change",max(abs(oldcoeff-res2$coeff)),"\n")

    
  if(max(abs(Mvec-OldMvec))<0.001 & max(abs(oldcoeff-res2$coeff))<0.02 & max(abs(neweps-oldeps))<0.001 & iteration>1){
    notConverged=FALSE
  }
  iteration=iteration+1;
  }
   
  res2 
}




myc = mydat %>% group_by(identifier) %>% summarize(n=n(),k=length(unique(Tr)),c=sum(Tr=="Control"))
mycfilt <- myc %>% filter(n>3,k>1,c>2)

mydat2 <- mydat %>% filter(identifier %in% mycfilt$identifier) %>% 
  mutate(offset=qlogis(DNA_prop)) %>%
  select(identifier,R=N,A=H,offset,batch,Tr)


dd <- mydat2 %>% group_by(identifier) %>% mutate(mean_RA = mean(R + A))

res2 <- fitQuasar2(dd)


mydat2 <- mydat %>% filter(identifier %in% mycfilt$identifier) %>% 
  mutate(offset=qlogis(DNA_prop),betacorr= betas.beta.binom-offset) %>%
  select(identifier,R=N,A=H,offset,batch,Tr,betacorr,betas_se) %>%
  group_by(identifier,Tr)

mymeta <- mydat2 %>% summarize(
  meta_estimate = weighted.mean(betacorr, 1 / betas_se^2),
  # Calculate the standard error of the meta-analysis estimate
  meta_se = sqrt(1 / sum(1 / betas_se^2)),
  # Count the number of observations in each group
  n = n()
) %>%
  # Remove the grouping applied by 'group_by'
  ungroup() %>% group_by(Tr)
  # Add columns for the meta-analysis z-score, p-value, and adjusted p-value
  mutate(
    meta_z = (meta_estimate) / meta_se,
    meta_p = 2 * (1 - pnorm(abs(meta_z))),
    meta_padj = p.adjust(meta_p,method="BH")
  )

myMetaCt <- filter(mymeta,Tr=="Control")
myMetaDex <- filter(mymeta,Tr=="Dex")
myMetaVitD <- filter(mymeta,Tr=="VitD")
myMetaVitA <- filter(mymeta,Tr=="VitA")


rescInt <- res2 %>% filter(term=="(Intercept)",convergence==0,se<100) %>% ungroup %>%  mutate(padj=p.adjust(pval))
sum(rescInt$padj<0.1)

resBatch <- res2 %>% filter(grepl("batch.*",term),convergence==0,se<100) %>% ungroup() %>% mutate(padj=p.adjust(pval))

sum(resBatch$padj<0.1)
qq(resBatch$pval)
max(resBatch$pval[resBatch$padj<0.1])

selBatch <- resBatch %>% filter(pval<0.05) %>% select(identifier) %>% unlist() %>% unique()

#Make qq plots and save dfs
resDex <- res2 %>% filter(term=="TrDex",convergence==0,se<100)
resDex$padj<-p.adjust(resDex$pval)
sum(resDex$padj<0.1)
qq(resDex$pval)

resVitD <- res2 %>% filter(term=="TrVitD",convergence==0,se<100)
resVitD$padj<-p.adjust(resVitD$pval)
sum(resVitD$padj<0.1)
qq(resVitD$pval)


resVitA <- res2 %>% filter(term=="TrVitA",convergence==0,se<100)
resVitA$padj<-p.adjust(resVitA$pval)
sum(resVitA$padj<0.1)
qq(resVitA$pval)



outPath <- "/nfs/rprscratch/wwwShare/carly/ASE_Plots/Roger4/"
system(paste0("mkdir -p ",outPath))


png(file=paste0(outPath,"dex_case_manually_calc_qq_linear_model.png"))
qq(resDex$pval)
dev.off()

combDex <- inner_join(resDex,myMetaCt,by="identifier") %>% inner_join(myMetaDex,by="identifier") %>%
  mutate(batchEff=(identifier %in% selBatch))



combDex %>% arrange(-padj) %>% ggplot(aes(meta_estimate.x,meta_estimate.y,color=padj<0.1)) +
  geom_point() + geom_abline(slope=1,intercept = 0) +
  theme_bw()

combDex %>% arrange(-padj) %>% ggplot(aes(meta_z.x,meta_z.y,color=padj<0.1)) +
  geom_point() + geom_abline(slope=1,intercept = 0) +
  theme_bw()

combDex %>% filter(padj<0.1) %>% arrange(-batchEff) %>% ggplot(aes(meta_estimate.x,meta_estimate.y,color=batchEff)) +
  geom_point() + geom_abline(slope=1,intercept = 0) +
  theme_bw()

combDex %>% filter(padj<0.1) %>% arrange(-batchEff) %>% ggplot(aes(meta_z.x,meta_z.y,color=batchEff)) +
  geom_point() + geom_abline(slope=1,intercept = 0) +
  theme_bw()



png(file=paste0(outPath,"vitd_case_manually_calc_qq_linear_model.png"))
qq(resVitD$pval)
dev.off()

png(file=paste0(outPath,"vita_case_manually_calc_qq_linear_model.png"))
qq(resVitA$pval)
dev.off()

png(file=paste0(outPath,"ase_manually_calc_qq_linear_model.png"))
qq(rescInt$pval)
dev.off()

p1 <- rescInt %>% ggplot(aes(coeff,fill=pval<0.05)) + geom_histogram(alpha=0.5,bins=100) + theme_bw()
ggsave(paste0(outPath,"ase_hist_bbl_model.png"),p1)


combInt <- inner_join(rescInt,myMetaCt,by="identifier") %>% inner_join(myMetaDex,by="identifier")

combInt %>% arrange(-padj) %>% ggplot(aes(meta_estimate.x,meta_estimate.y,color=padj<0.1)) +
  geom_point() + geom_abline(slope=1,intercept = 0) +
  theme_bw()

combInt %>% arrange(-padj) %>% ggplot(aes(meta_z.x,meta_z.y,color=padj<0.1)) +
  geom_point() + geom_abline(slope=1,intercept = 0) +
  theme_bw()



outPath <- "/rs/rs_grp_neanbit/RogerStuff/BetaBinomialLogisticResEtohOnlyV4/"
system(paste0("mkdir -p ",outPath))

write.table(resDex,file=paste0(outPath,"case_linearmodel_Bitstarr_dex.txt"),quote=FALSE,row.names=FALSE,sep="\t")
##write.table(resSelenium,file=paste0(outPath,"case_linearmodel_Bitstarr_selenium.txt"),quote=FALSE,row.names=FALSE,sep="\t")
##write.table(resCaffeine,file=paste0(outPath,"case_linearmodel_Bitstarr_caffeine.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(resVitD,file=paste0(outPath,"case_linearmodel_Bitstarr_vitd.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(resVitA,file=paste0(outPath,"case_linearmodel_Bitstarr_vita.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(rescInt,file=paste0(outPath,"ase_linearmodel_Bitstarr.txt"),quote=FALSE,row.names=FALSE,sep="\t")



