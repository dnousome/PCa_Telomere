#################################################
#
#
#
#
#
#
#################################################

##References
#https://cebp.aacrjournals.org/content/29/3/676
#https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(19)30831-X/fulltext
library(biomaRt)
library(tidyverse)


#snps=read_tsv("~/Documents/xfer/Darryl/Telomere.tsv")
#write_tsv(snps,"TSNPs.tsv")

snps=read_tsv("TSNPs.tsv") %>% separate(`Allele*`,into=c("REF","ALT"))

ensembl = useMart("ENSEMBL_MART_SNP")
#listDatasets(ensembl)
ensembl=useDataset("hsapiens_snp",mart=ensembl)
snppos=getBM(attributes=c("chr_name",'chrom_start',"chrom_end",'refsnp_id'),filters='snp_filter',values=snps$`SNP ID`,mart=ensembl)                              


annot=snppos %>% full_join(.,snps,by=c("refsnp_id"="SNP ID"))

gt=apply(annot,1,function(x){
  system(sprintf("/home/dnousome/miniconda3/envs/bio/bin/bcftools query -f '%%CHROM\t%%POS\t%%REF\t%%ALT\t[%%GT\t]\n' /home/dnousome/Documents/Glenn/DDRG/Data/TAGC/1_VCF/2_FixedPloidy/ddrg_chr%s_setGT.bcf -r chr%s:%s-%s",
          as.numeric(x[1]),as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[3])),intern = T)
})
ids=system("/home/dnousome/miniconda3/envs/bio/bin/bcftools query -l  /home/dnousome/Documents/Glenn/DDRG/Data/TAGC/1_VCF/2_FixedPloidy/ddrg_chr1_setGT.bcf",intern=T)

gt1=lapply(1:length(gt),function(x){
  s=strsplit(gt[[x]],"\t")
  
  annot_temp=annot[x,]
  
  flip=ifelse(s[[1]][4]==annot_temp$ALT,"CORRECT",
         ifelse(s[[1]][3]==annot_temp$ALT,"FLIP",NA))
 
  
  counts=ifelse(s[[1]][-1:-4]=="./.",NA,str_count(s[[1]][-1:-4],"1"))
  flipcounts=recode(counts, `0` = 2,`2`=0,`1`=1)
  if(flip=="CORRECT"){
    counts
    }else if(flip=="FLIP"){
    flipcounts
    }else {NA}
    
  })

gt2=data.frame(do.call(rbind,gt1))
names(gt2)=ids
gt3=data.frame(Telo=colSums(gt2,na.rm=T)) %>%rownames_to_column("ID")
hist(gt3$Telo)


#######Read in the covariate data

#########Covariate
covar1=read_sas("~/Documents/G/Work_Bench/Epidemiology/Projects/Petrovics 2018 DDRG prec med grant/data/final10032019wbcr.sas7bdat")

covar=read_xlsx("~/Documents/G/Work_Bench/Epidemiology/Projects/Petrovics 2018 DDRG prec med grant/data/DDRG clinical data to TAGC 05022019.xlsx") %>% 
  arrange(match(`Sample ID`,covar1$LOG_)) %>% 
  mutate(PID=as.character(covar1$PID)) %>%
  mutate(BCR=covar1$BCR) %>%
  mutate(time_RP_BCR=covar1$time_RP_BCR) %>%
  mutate(time_Dx_mets=covar1$time_Dx_mets) %>%
  mutate(BCR=as.factor(BCR))


###Read in Fx
covar=read_xlsx("~/Documents/Glenn/DDRG/Data/DDRG_Fam_Hx_CaP_to Darryl 03022020.xlsx") %>%
  mutate(PID=as.character(PID)) %>%
  left_join(covar,.,by='PID') %>%
  mutate(FAMILY_CAP=ifelse(FAMILY_CAP=="Yes",1,ifelse(FAMILY_CAP=="No",0,NA)))

taglinker=read_delim("~/Documents/Glenn/DDRG/Data/cpdr.sample.lookup.txt",col_names = F,delim=" ")
#tagex=read_delim("~/Documents/Aldrin/DDRG/Data/cpdr.exclude.sample.qa.txt",col_names = F,delim=" ")


covarf=taglinker %>% left_join(.,covar,by=c("X2"="Sample ID")) %>% 
  left_join(.,gt3,by=c(X1="ID")) %>%
  mutate(teloq=cut(Telo,breaks = 2)) %>%
  mutate(BCR=as.factor(BCR))
            

summary(glm(as.factor(mets_status)~teloq+RACE,family = binomial,data=covarf))

library(survminer)
library(survival)
fit <- survfit(Surv(time=time_RP_BCR,as.numeric(BCR)) ~ teloq,
               data = covarf)
# Visualize with survminer
ggsurvplot(fit, data = covarf %>% filter(!is.na(time_RP_BCR)&!is.na(BCR)&!is.na(teloq)), risk.table = T)

coxph(Surv(time_RP_BCR, as.numeric(BCR))~teloq+RACE+Dx_age,data=covarf)
