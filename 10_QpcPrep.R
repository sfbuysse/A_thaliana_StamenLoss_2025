# Sophie Buysse
# 8/3/2023
# updated 11/16/2023
# making file with 50,000 randomly selected SNPs for Qpc

# HPCC set up if looking at this and not running it as a job
#module purge
#module load GCC/11.2.0  OpenMPI/4.1.1 R/4.2.2
# if submitting, need:
#SBATCH --constraint="[intel16|intel18|amd20|amd22]"

##### load libraries ######
library(genio)

##### Genotype Info #####

# set seed
set.seed(803)

# load files
# want the file where each input is a 0, 1, 2 - NOT where it is already in letter format.
# this is reading in the allsites plink filtered file that is the same input file for Fst, PCA, and GWAS.
SNPs_Cent <- read_plink("/mnt/scratch/buysseso/GVCF/allsites_filtered_plinkTest", verbose = TRUE)
# I think X is what I want
head(SNPs_Cent$X)[1:5, 1:5]
# check for SNPs with no missing data
complete <- SNPs_Cent$X[complete.cases(SNPs_Cent$X), ]
str(complete)
# yay! 172591 of these exist!!!
# if we don't do this subsetting way, I will need to do imputation

# now transpose
Stamen_Gt <- t(complete)
head(Stamen_Gt)[1:5, 1:5]

# and calculate into frequencies
Stamen_Gt <- Stamen_Gt / 2
head(Stamen_Gt)[1:5, 1:5]

## subset 50,000 randomly chosen SNPs ##
# this does include the centromeres.

# SNPs are cols, so find n.cols
ncol(Stamen_Gt)
cols_keep <- sample(c(1:ncol(Stamen_Gt)), size = 50000, replace = FALSE, prob = NULL)
# double check no repeat numbers
length(unique(cols_keep))
# the individual ID is the rowname, and I need to make sure that stays.
Stamen_Gt_sub <- Stamen_Gt[ , cols_keep]
head(Stamen_Gt_sub)[1:5, 1:5]
length(unique(colnames(Stamen_Gt_sub)))
str(Stamen_Gt_sub)
class(Stamen_Gt_sub)

## change id structure so can merge with phenotype info later on
# get rid of the dash in some of the rownames
rownames(Stamen_Gt_sub) <- gsub(pattern = "-", replacement = "", rownames(Stamen_Gt_sub))
# make all uppercase
rownames(Stamen_Gt_sub) <- toupper(rownames(Stamen_Gt_sub))
#make it a dataframe to merge
Stamen_Gt_sub <- as.data.frame(Stamen_Gt_sub)
Stamen_Gt_sub$LineID <- rownames(Stamen_Gt_sub)
# and organize so it is the first column
Stamen_Gt_sub <- Stamen_Gt_sub[ ,c(50001, 1:50000)]
head(Stamen_Gt_sub)[1:5, 1:5]

save(Stamen_Gt_sub, file = "/mnt/research/josephslab/Sophie/Qpc/GenotypeMatrix_Cent_50k_Nov2023.ROBJ")
write.csv(Stamen_Gt_sub, file = "/mnt/research/josephslab/Sophie/Qpc/GenotypeMatrix_Cent_50k_Nov2023.csv", row.names = FALSE)

##### Phenotype Information #####
# read in FAM... how have I done this before?
# can read in no cent or cent, the fam files are identical but file names need to match geno files to run the GWAS so there are copies
pheno <- read.table("/mnt/gs21/scratch/buysseso/GWAS/NoCent.PlinkFiltering_raw.fam")
head(pheno)
# look to be in same order as geno file, but let's make the case match too just cause
# get rid of the dash in some of them
pheno$V1 <- gsub(pattern = "-", replacement = "", pheno$V1)
head(pheno)
# make all uppercase
pheno$V1 <- toupper(pheno$V1)
# look at all to see if V1 and V2 match
pheno
# looks good to me.
# let's save just the columns I need
pheno2 <- pheno[ , c("V1", "V6")]
colnames(pheno2) <- c("LineID", "Pheno")
head(pheno2)

# save them
save(pheno2, file = "/mnt/research/josephslab/Sophie/Qpc/RawPhenotypes_Aug2023.ROBJ")
write.csv(pheno2, file = "/mnt/research/josephslab/Sophie/Qpc/RawPhenotypes_Aug2023.csv", row.names = FALSE)