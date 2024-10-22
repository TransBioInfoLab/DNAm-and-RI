## Colocalisation analysis (resilience study)
## mQTLs (11 CpGs - 1110 SNPs) / GWAS results 
## Lissette Gomez, October 2024


setwd("M:/AD/analysis/projects/Sex_strat/Aim3/Meta-analysis/Lissette/bloodMeta_mQTLs/resilience")

mQTL_GWAS_overlap <- read.csv("GWASloci_resilienceBlood_overlap_bellenguez.csv")


library(coloc)
mQTL_GWAS_overlap$chrPos <- paste0(mQTL_GWAS_overlap$seqnames, ":", mQTL_GWAS_overlap$start) #hg38
SNPs <- unique(mQTL_GWAS_overlap$chrPos) #n=549


#Find the SNPs with betas for the same allele in both datasets
mQTL_GWAS_overlap_sameEff <- mQTL_GWAS_overlap[which(mQTL_GWAS_overlap$allele1 == mQTL_GWAS_overlap$effect_allele),]
length(unique(mQTL_GWAS_overlap_sameEff$chrPos)) #301

#Find SNPs with different effect allele in GWAS and mQTL results
diffEffectAllesChrPos <- SNPs[which(!SNPs%in%mQTL_GWAS_overlap_sameEff$chrPos)] #n=248

#Change Beta and alleles for SNPs with different effect allele in each dataset
mQTL_GWAS_overlap_diffEff <- mQTL_GWAS_overlap[which(mQTL_GWAS_overlap$chrPos %in% diffEffectAllesChrPos),]

mQTL_GWAS_overlap_diffEff$beta <- -mQTL_GWAS_overlap_diffEff$beta
nonEffectAlleleFix <- mQTL_GWAS_overlap_diffEff$effect_allele
EffectAlleleFix <- mQTL_GWAS_overlap_diffEff$other_allele
mQTL_GWAS_overlap_diffEff$effect_allele <- EffectAlleleFix
mQTL_GWAS_overlap_diffEff$other_allele <- nonEffectAlleleFix




#Concat SNPs with same and different effect allele after changing the beta 
mQTL_GWAS_overlap_fixEffectAllele <- rbind(mQTL_GWAS_overlap_sameEff, mQTL_GWAS_overlap_diffEff)

mQTL_GWAS_overlap_fixEffectAllele$SNP <- paste0(mQTL_GWAS_overlap_fixEffectAllele$chrPos,
                                          "_",
                                          mQTL_GWAS_overlap_fixEffectAllele$effect_allele,
                                          "_",
                                          mQTL_GWAS_overlap_fixEffectAllele$other_allele)

d2 <- unique(mQTL_GWAS_overlap_fixEffectAllele[, c("SNP", "beta", "standard_error")])
d1 <- unique(mQTL_GWAS_overlap_fixEffectAllele[, c("SNP", "beta_a1", "se", "samplesize", "freq_a1")])

##Run co-localisation test by CpG
CpG <- unique(mQTL_GWAS_overlap$cpg) #14

dataset2<-list(beta= d2$beta,
               varbeta= (d2$standard_error^2),
               type = "cc",
               snp = d2$SNP,
               s=0.21,
               N=487511)

summaryRes <- data.frame(matrix(ncol=7,nrow=0))
allRes <- data.frame(matrix(ncol=12,nrow=0))
for (i in CpG[1:length(CpG)]){
  
  mQTLcpg <- mQTL_GWAS_overlap_fixEffectAllele[which(mQTL_GWAS_overlap_fixEffectAllele$cpg == i),] 
  
  dataset1<-list(beta=mQTLcpg$beta_a1,
                 varbeta = mQTLcpg$se^2,
                 type = "quant",
                 snp = mQTLcpg$SNP,
                 N = mQTLcpg$samplesize,
                 MAF = mQTLcpg$freq_a1
  )
  
  my.res<-coloc.abf(dataset1, dataset2)
  
  allRes <- rbind(allRes, cbind(i,my.res$results))
  summaryRes <- rbind(summaryRes, cbind(i, as.data.frame(t(my.res$summary))))
  
  
}

colnames(allRes)[1] <- "CpG"
colnames(summaryRes)[1] <- "CpG"

options(scipen=999)
write.csv(allRes, "colocalisationResults_all_mQTLs_Bellenguez.csv", row.names = FALSE)
write.csv(summaryRes, "colocalisationResults_summary_mQTLs_Bellenguez.csv", row.names = FALSE)


