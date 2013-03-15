#!/bin/bash

#This script is for performing GWAS using the QIMR GWAS dosage data stored on /hpscratch . Note that PLINK can also be used to perform GWAS with dosage data.
#However many of PLINK's normal functions don't work with dosage scores and this script provides a framework to add code in an R script to do other functions.
#The original reason was for writing the script was to look at SNP x sex interaction.
#This code will write R scripts for each chromosomal segment. These scripts should be submitted using Guo-bo Chen's Rsub.R script also located in the CTGG github.


#Here you can list your home directory
export GR=~/Jack_anxiety/1000k/

cd ${GR}

for ((i=1;i<=22;i++))

#for i in 22;
do

for((j=1;j<=100;j++))
#for j in 1;
do

if test -f /hpscratch/wrayvisscher/GWAS/GWAS/GeneralRelease/Imputed/Release6_1000G_281K/infer/1000G_281K_infer_QCpass_chr"$i"."$j".dat.gz
then

#Here the program will write multiple R scripts - one for each file in the 1KG QIMR dataset that are ready to run the GWAS
echo "data <- read.table(gzfile(\"/hpscratch/wrayvisscher/GWAS/GWAS/GeneralRelease/Imputed/Release6_1000G_281K/infer/1000G_281K_infer_QCpass_chr$i.$j.ped.gz\"))

#change the name of the second column so that it can be merged later with .ped file 
names(data)[2] <- \"IID\"

require(stats)

store <- matrix(nrow=1,ncol=4)

interaction <- matrix(nrow=1,ncol=4)

setwd(\"${GR}\")

#Here you read in the name of your pedigree file generated from e.g. Scott Gordon's GWAS_ID_remapper program
e <- read.table(\"FS2_mapped.txt\",header=F,colClasses=c(\"character\",\"character\",\"numeric\"))

#change name of columns for merging with genotypes
names(e) <- c(\"FID\",\"IID\",\"FS2\")

#Here I have a separate covariates file that merges in, but no need to have covariate and phenotype information separately.
f <- read.table(\"covars_mapped.txt\",header=F,colClasses=c(\"character\",\"character\",\"numeric\",\"numeric\"))

# Now read in Principal Components file for use as covariates in regression.
PC <- read.table(\"/hpscratch/wrayvisscher/GWAS/GWAS/GeneralRelease/Release6/info/GWAS_PrincipalComponentsScores.txt\",header=T)

#Select ID columns and 1st four PCs
PC_trim <- PC[,c(1,2,6,7,8,9,10)]

#Change names for merging
names(PC_trim)[1] <- \"FID\"

names(PC_trim)[2] <- \"IID\"

names(f) <- c(\"FID\",\"IID\",\"SEX\",\"AGE\")

# Merge covariates and PCs
f_PC <- merge(f,PC_trim,by=c(\"FID\",\"IID\"),all.x=T,all.y=F)

# Merge phenotype and covariate files (not strictly necessary to keep them separate - see above)
g <- merge(e,f_PC,by=\"IID\")

# Merge p'types, covars and genotypes
pedfile <- merge(data,g,by=\"IID\")

# Test null model with only covariates
null <- lm(pedfile\$FS2~pedfile\$SEX+pedfile\$AGE+pedfile\$PC1+pedfile\$PC2+pedfile\$PC3+pedfile\$PC4,na.action=na.exclude)

LRT_null <- logLik(null)

# Now test SNPs
# The first 5 columns of data file are IDs, so only use those columns
for (k in 1:(ncol(data)-5))

# First test SNP effect without age and sex covariates.
{

snp <- lm(pedfile\$FS2~pedfile[,(k+5)]+pedfile\$PC1+pedfile\$PC2+pedfile\$PC3+pedfile\$PC4,na.action=na.exclude)

snpstore <- matrix(nrow=1,ncol=4)

# Use if statement to identify missing SNPs
if (nrow(coef(summary(snp))) > 5) {

for (l in 1:length(snpstore))

{
# extract regression statistics for the SNP which is second row of coef(summary) 
snpstore[1,l] <- coef(summary(snp))[2,l]

}

} else {snpstore[1,] <- \"NA\"}

write.table(snpstore,\"FS2_1KG_nocovars_chr${i}_${j}_results.txt\",row.names=F,col.names=F,quote=F,append=T)

#Now run while including covariates

lin <- lm(pedfile\$FS2~pedfile[,(k+5)]+pedfile\$SEX+pedfile\$AGE+pedfile\$PC1+pedfile\$PC2+pedfile\$PC3+pedfile\$PC4,na.action=na.exclude)

if (nrow(coef(summary(lin))) > 7) {

for (l in 1:ncol(store))

{

store[1,l] <- coef(summary(lin))[2,l]

}

} else {store[1,1:4] <- \"NA\"}

write.table(store,\"FS2_1KG_covars_chr${i}_${j}_results.txt\",row.names=F,col.names=F,quote=F,append=T)

# Optional section below to look at SNP by sex interaction. Note that main effects cannot be correctly interpreted
inter <- lm(pedfile\$FS2~pedfile[,(k+5)]*pedfile\$SEX+pedfile\$AGE+pedfile\$PC1+pedfile\$PC2+pedfile\$PC3+pedfile\$PC4,na.action=na.exclude)

if (nrow(coef(summary(inter))) > 8) {


for (l in 1:ncol(interaction))

{

interaction[1,l] <- coef(summary(inter))[9,l]

}

} else {interaction[1,1:4] <- \"NA\"}

write.table(interaction,\"FS2_1KG_interaction_chr${i}_${j}_results.txt\",row.names=F,col.names=F,quote=F,append=T)
}

" > Rscript_FS2_covars_chr${i}_${j}.R

fi

done

done


# Now run the generated R scripts with Guo-bo's Rsub script e.g. Rsub Rscript Rscript_FS2_covars_chr1_1.R @FS2_chr1_1_submit




