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

logLikBetaBinomialLogistic <- function(b,X,D,R,A,offset=0.0){
  p <- plogis(X %*% b + offset ) ## maybe add offset
  aux <- (lgamma(D) + lgamma(R+p*D) + lgamma(A + (1-p)*D) - lgamma(R+A+D) - lgamma(p*D) - lgamma((1-p)*D));
  -sum(aux)
}

gLogLikBetaBinomialLogistic <- function(b,X,D,R,A,offset=0.0){
  p <- plogis(X %*% b + offset) 
  aux <- (p * (1-p))*D*(digamma(R+p*D) - digamma(A + (1-p)*D) - digamma(p*D) + digamma((1-p)*D)) 
  -t(X) %*% aux
}



## Consider a ridge/lasso penalty? 

M=10 
## ##"chr10:10741817_rev"
## ##"chr10:111658679_rev"
##. dd <- mydat %>% filter(identifier=="chr10:111658679_fwd") ## %>%
## ###  mutate(DNAq=qlogis(DNA_prop)) %>%
##  select(N,H,DNA_prop,Tr,batch,identifier) 

fitBetaBinomialLogistic <- function(dd,M=3.3){
  X <- model.matrix(~ batch + Tr,dd)
  b <- rep(0,ncol(X))
  auxLogis <- optim(b,
                    fn=logLikBetaBinomialLogistic,
                    gr=gLogLikBetaBinomialLogistic,
                    X=X,
                    D=M,
                    R=dd$N,
                    A=dd$H,
                    offset=qlogis(dd$DNA_prop),
                    method="L-BFGS-B",
                    hessian=TRUE,
                    lower=-8,
                    upper=8)
  ##auxLogis$convergence
  ##auxLogis$message
  p.fit <- plogis(X %*% auxLogis$par + qlogis(dd$DNA_prop))
  
  res <- tibble( identifier=dd$identifier[1],term = colnames(X),coeff=auxLogis$par,
          se=1/diag(auxLogis$hessian)^.5,
         convergence=auxLogis$convergence,
         qr.rank=qr(X)$rank,df=ncol(X),n=nrow(X)) %>% 
    mutate(stat=coeff/se,pval=2*pnorm(-abs(stat)))
  list(res=res,p.fit=p.fit)
  }


myc = mydat %>% group_by(identifier) %>% summarize(n=n(),k=length(unique(Tr)),c=sum(Tr=="Control"))
mycfilt <- myc %>% filter(n>3,k>1,c>2)

mydat2 <- mydat %>% filter(identifier %in% mycfilt$identifier)

aux <- mydat2 %>% group_by(identifier) %>% nest()
##res <- aux %>% mutate(res=map(data,mylm))
res <- aux %>% mutate(res=map(data,fitBetaBinomialLogistic))

# Extracting the p.fit  and combining with the data
res <- res %>% mutate(p.fit=map(res, "p.fit"),res=map(res,"res"))

res3 <- res %>% 
  mutate(data = map2(data, p.fit, ~mutate(.x, p.fit=.y))) %>% 
  select(identifier, data) %>% 
  unnest(cols=c(data))


res2 <- res %>% select(identifier=identifier,data=res) %>% unnest(cols=c(data))


rescInt <- res2 %>% filter(term=="(Intercept)",convergence==0,se<100) %>% ungroup %>%  mutate(padj=p.adjust(pval))
sum(rescInt$padj<0.1)

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

resCaffeine <- res2 %>% filter(term=="TrCaffeine",convergence==0,se<100)
resCaffeine$padj<-p.adjust(resCaffeine$pval)
sum(resCaffeine$padj<0.1)
qq(resCaffeine$pval)

resSelenium <- res2 %>% filter(term=="TrSelenium",convergence==0,se<100)
resSelenium$padj<-p.adjust(resSelenium$pval)
sum(resSelenium$padj<0.1)
qq(resSelenium$pval)


logLikBetaBinomialRhoEps <- function(rho,eps,D,R,A){
  p <- (rho*(1-eps)+(1-rho)*eps)
  aux <- (lgamma(R+p*D) + lgamma(A + (1-p)*D) - lgamma(R+A+D) - lgamma(p*D) - lgamma((1-p)*D) + lgamma(D)) ##+ lgamma(R+A+1) - lgamma(A+1) - lgamma(R+1)
  aux
}


logLikBetaBinomialM <- function(D,p,R,A){
#  eps=0.001
#  p <- (p*(1-eps)+(1-p)*eps)
  aux <- (lgamma(D) + lgamma(R+p*D) + lgamma(A + (1-p)*D) - lgamma(R+A+D) - lgamma(p*D) - lgamma((1-p)*D) );
  -sum(aux)
}

gLogLikBetaBinomialM <- function(D,p,R,A){
#  eps=0.001
#  p <- (p*(1-eps)+(1-p)*eps)
  aux <- (digamma(D) + p * digamma(R + p * D) + (1 - p) * digamma(A + (1 - p) * D) - digamma(R+A+D) - p * digamma(p * D) - (1-p) * digamma((1 - p) * D ))
  -sum(aux)
}

logLikBetaBinomialM(50, 
                    p=res3$p.fit,
                    R=res3$N,
                    A=res3$H)

gLogLikBetaBinomialM(50, 
                    p=res3$p.fit,
                    R=res3$N,
                    A=res3$H)


mt <- res3 %>% summarize(mt=mean(N+H)) 

hist(log10(mt$mt),breaks=20)

selid <- mt$identifier[mt$mt>1 & mt$mt<6000000]

selv <- res3$identifier %in% selid

sum(selv)


auxLogis <- optim(par=10,
                  fn=logLikBetaBinomialM,
                  gr=gLogLikBetaBinomialM,
                  p=res3$p.fit[selv],
                  R=res3$N[selv],
                  A=res3$H[selv],
                  method="L-BFGS-B",control=c(trace=6),
                  hessian=TRUE,
                  lower=0.1,
                  upper=10000)


auxLogis$convergence

auxLogis$message

c(auxLogis$par,1/(auxLogis$hessian)^.5)


M <- exp((-20:50)/5)
aux <- sapply(M,function(M){
  (logLikBetaBinomialM(M,p=res3$p.fit,
                               R=res3$N,
                               A=res3$H))
})

aux2 <- sapply(M,function(M){
  (gLogLikBetaBinomialM(M,p=res3$p.fit,
                       R=res3$N,
                       A=res3$H))
})


Mmax <- M[which.min(aux)]
cat(" ",which.max(aux)," ",Mmax,"\n")

plot(log10(M),-aux/max(aux))

lines(log10(M),-aux2/max(aux2))
abline(h=0)

M <- exp((-50:50)/5)
eps=0.001
aux <- sapply(M,function(M){
  sum(logLikBetaBinomialRhoEps(mydat2$DNA_prop,
                               eps,M,
                               mydat2$N,
                               mydat2$H))
})
Mmax <- M[which.max(aux)]
cat(" ",which.max(aux)," ",Mmax,"\n")


M <- exp((-50:50)/5)
eps=0.001
aux <- sapply(M,function(M){
  sum(logLikBetaBinomialRhoEps(res3$p.fit,
                               eps,M,
                               res3$N,
                               res3$H))
})
Mmax <- M[which.max(aux)]
cat(" ",which.max(aux)," ",Mmax,"\n")



outPath <- "/nfs/rprscratch/wwwShare/carly/ASE_Plots/Roger/"
system(paste0("mkdir -p ",outPath))


png(file=paste0(outPath,"dex_case_manually_calc_qq_linear_model.png"))
qq(resDex$pval)
dev.off()

png(file=paste0(outPath,"vitd_case_manually_calc_qq_linear_model.png"))
qq(resVitD$pval)
dev.off()

png(file=paste0(outPath,"vita_case_manually_calc_qq_linear_model.png"))
qq(resVitA$pval)
dev.off()

png(file=paste0(outPath,"caf_case_manually_calc_qq_linear_model.png"))
qq(resCaffeine$pval)
dev.off()

png(file=paste0(outPath,"selenium_case_manually_calc_qq_linear_model.png"))
qq(resSelenium$pval)
dev.off()

png(file=paste0(outPath,"ase_manually_calc_qq_linear_model.png"))
qq(rescInt$pval)
dev.off()

outPath <- "/rs/rs_grp_neanbit/RogerStuff/BetaBinomialLogisticRes/"
system(paste0("mkdir -p ",outPath))

write.table(resDex,file=paste0(outPath,"case_linearmodel_Bitstarr_dex.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(resSelenium,file=paste0(outPath,"case_linearmodel_Bitstarr_selenium.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(resCaffeine,file=paste0(outPath,"case_linearmodel_Bitstarr_caffeine.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(resVitD,file=paste0(outPath,"case_linearmodel_Bitstarr_vitd.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(resVitA,file=paste0(outPath,"case_linearmodel_Bitstarr_vita.txt"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(rescInt,file=paste0(outPath,"ase_linearmodel_Bitstarr.txt"),quote=FALSE,row.names=FALSE,sep="\t")

