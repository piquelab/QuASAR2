library(data.table)
library(qqman)
library(ggplot2)
library(tidyverse)
library(QuASAR2)


mydat=fread('/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/ASE/all_allOutput.txt')

cov.file <- "/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/covariates/GMCcol_deep_covariates.txt"

cv <- read.table(cov.file, sep="\t", stringsAsFactors=F, header=T, comment="")

mydat2 <- select(mydat,chr,pos,ref,alt,rsID,SampleID=treatment,cell.line,R=ref.reads,A=alt.reads)

mydat2$R <- as.numeric(mydat2$R)
mydat2$A <- as.numeric(mydat2$A)



colnames(mydat)["treatment"] <- "SampleID"

cv2 <- select(cv,c(SampleID,Colonoid_line,Treatment,Microbiome,Urbanism,Prevotella,Bacteroides,Firmicutes))

mydat2 <- left_join(mydat2,cv2)

mydat2 <- filter(mydat2,!is.na(mydat2$R))

table(mydat2$cell.line,mydat2$Colonoid_line)

table(mydat2$Colonoid_line,mydat2$Treatment)


mydat2$Tr=mydat2$Treatment

mydat2$Tr = factor(mydat2$Tr) %>% relevel("control")

mydat2$identifier <- paste0(mydat2$rsID,":",mydat2$Colonoid_line)

##mydat <- mutate(mydat, betacorr =  betas.beta.binom - qlogis(DNA_prop))

nbreaks = 20; M = 20; eps = 0.001;


myc = mydat2 %>% group_by(identifier) %>% summarize(n=n(),k=length(unique(Tr)),c=sum(Tr=="control"))
mycfilt <- myc %>% filter(n>3,k>1,c>2)

mydat3 <- mydat2 %>% filter(identifier %in% mycfilt$identifier) %>% 
  mutate(offset=0) %>%
  select(identifier,R,A,offset,SampleID,Colonoid_line,Tr)


dd <- mydat3 %>% group_by(identifier) %>% mutate(mean_RA = mean(R + A))

###
resObj <- fitQuasar2(dd,~Tr)

res2 <- resObj$result

resObj$Mres

resObj$Mvec

table(res2$convergence,res2$qr.rank)

####
rescInt <- res2 %>% filter(term=="(Intercept)",convergence==0,se<100) %>% ungroup %>%  mutate(padj=p.adjust(pval))
sum(rescInt$padj<0.1)

#Make qq plots and save dfs
resTrM <- res2 %>% filter(term=="Trmicrobiome",convergence==0,se<100)
resTrM$padj<-p.adjust(resTrM$pval)
sum(resTrM$padj<0.1)

qq(resTrM$pval)

qq(rescInt$pval)


outPath <- "/rs/rs_grp_gmccol/QuASAR2_results/"
system(paste0("mkdir -p ",outPath))


png(file=paste0(outPath,"Microbiome_qq_v2.png"))
qq(resTrM$pval)
dev.off()


png(file=paste0(outPath,"Microbiome_ASE_int_qq_v2.png"))
qq(rescInt$pval)
dev.off()

write_tsv(rescInt,file=paste0(outPath,"Microbiome_ASE_v2_table.tsv"))

write_tsv(resTrM,file=paste0(outPath,"Microbiome_cASE_v2_table.tsv"))
