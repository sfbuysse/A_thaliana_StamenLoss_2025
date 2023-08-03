# Sophie Buysse
# 8/3/2023
# making file with 50,000 randomly selected SNPs for Emily to do Qpc with

# HPCC set up if looking at this and not running it as a job
#module purge
#module load GCC/11.2.0  OpenMPI/4.1.1 R/4.2.2
# if submitting, need:
#SBATCH --constraint="[intel16|intel18|amd20|amd22]"

##### Genotype Info #####

# set seed
set.seed(803)

# load files
load(file = "/mnt/research/josephslab/Sophie/Qpc/GenotypeMatrix_plinkFiltered_Cent_Dec2022.ROBJ")
cent_matrix <- tmp2
rm(tmp2)
load(file = "/mnt/research/josephslab/Sophie/Qpc/GenotypeMatrix_plinkFiltered_NoCent_Dec2022.ROBJ")
nocent_matrix <- tmp2
rm(tmp2)

## subset 50,000 randomly chosen SNPs for fst matrix ##

# this does include the centromeres.
length(cent_matrix$SNP_ID)
rows_keep <- runif(n = 50000, min = 2, max = length(cent_matrix$SNP_ID))
rows_keep <- c(1,rows_keep)
# need to keep row one for the ID, set rows_keep to start at two so I never have row 1 twice
# want to keep all columns because each column is an individual, only filter rows
cent_sub <- cent_matrix[rows_keep, ]
str(cent_sub)
class(cent_sub)
length(cent_sub$SNP_ID)
#50,001

save(cent_sub, file = "/mnt/research/josephslab/Sophie/Qpc/GenotypeMatrix_Cent_50k_Aug2023.ROBJ")
write.csv(cent_sub, file = "/mnt/research/josephslab/Sophie/Qpc/GenotypeMatrix_Cent_50k_Aug2023.csv", row.names = FALSE)

# this does not include the centromeres.
length(nocent_matrix$SNP_ID)
rows_keep2 <- runif(n = 50000, min = 2, max = length(nocent_matrix$SNP_ID))
rows_keep2 <- c(1,rows_keep2)
# need to keep row one for the ID, set rows_keep to start at two so I never have row 1 twice
# want to keep all columns because each column is an individual, only filter rows
nocent_sub <- nocent_matrix[rows_keep, ]
str(nocent_sub)
class(nocent_sub)
length(nocent_sub$SNP_ID)
#50,001

save(nocent_sub, file = "/mnt/research/josephslab/Sophie/Qpc/GenotypeMatrix_NoCent_50k_Aug2023.ROBJ")
write.csv(nocent_sub, file = "/mnt/research/josephslab/Sophie/Qpc/GenotypeMatrix_NoCent_50k_Aug2023.csv", row.names = FALSE)

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