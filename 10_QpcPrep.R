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
# but I want it transposed
Stamen_Gt <- t(SNPs_Cent$X)
head(Stamen_Gt)[1:5, 1:5]

## subset 50,000 randomly chosen SNPs ##
# this does include the centromeres.

# SNPs are cols, so find n.cols
ncol(Stamen_Gt)
cols_keep <- runif(n = 50000, min = 1, max = ncol(Stamen_Gt))
# double check no repeat numbers
length(unique(cols_keep))
# the individual ID is the rowname, and I need to make sure that stays.
Stamen_Gt_sub <- Stamen_Gt[ , cols_keep]
head(Stamen_Gt_sub)[1:5, 1:5]

str(Stamen_Gt_sub)
class(Stamen_Gt_sub)
tmp <- Stamen_Gt_sub[!complete.cases(Stamen_Gt_sub), ]
# hmm. these NAs will be a problem later on b/c make_k would make a matrix of all NA (even after dividing by 2).
# so do I need to.... do what?
# quaint documentation says to use random imputation to replace missing genotypes. how??

# what if I subsample just the complete cases from the initial matrix? are there any?
tmp2 <- Stamen_Gt[complete.cases(Stamen_Gt), ]
# HA it is empty.
# I think I need another program to do imputation, and that would happen before I even start the subsetting in R.
# ok. well it looks like I"ll need to use BEAGLE to do this. and that requires data in a vcf format. and downloading another program and stuff
# I don't have it in me to do that now.



#save(cent_sub, file = "/mnt/research/josephslab/Sophie/Qpc/GenotypeMatrix_Cent_50k_Aug2023.ROBJ")
#write.csv(cent_sub, file = "/mnt/research/josephslab/Sophie/Qpc/GenotypeMatrix_Cent_50k_Aug2023.csv", row.names = FALSE)


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