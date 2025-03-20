# *Arabidopsis thaliana* Stamen Loss
Pipeline for the manuscript "Evaluating the Roles of Drift and Selection in Trait Loss along an Elevational Gradient" - see our preprint here: https://www.biorxiv.org/content/10.1101/2024.06.12.598645v1

The number in front of each file indicates the order of analyses. The code in the main repository is cleaned up for ease of reuse and understanding the analyses in the manuscript. For completeness, all of the code is reported in the "MessyCode" folder.

The sequence information will be publicly available on the SRA: LINK. The only processing before submitting to SRA was to combine all forward and reverse files from the same lane for each genotype (as the genotypic information was split into 1 to 7 files per genotype due to size). You can find this code at the end of the MessyCode/01_Fastq_to_Bam_allmycode.txt.

Phenotypic information is in the data folder.

Due to their size, GWAS output files are on Dryad: LINK. Most intermediates files are not shared here because they are large.

# File Descriptions:
- 01_FastQ_to_BAM: Starts with fastq files. The end result is aligned to the reference genome but not cleaned at all.
- 01B_Coverage Statistics: calculating coverage statistics for the supplemental table
- 02_CallVariants: Starts with aligned reads. The end result is the raw allsites file of genotypes.
- 03_FilterVariants: Code to go from raw to filtered variant and invariant sites before doing pixy (vcf) and to make the filtered variant sites file used for Fst (later cut from analysis), GWAS, and PCA (plink format). It also removes the centromere regions post filtering.
- 04_PopulationRawMeans_LSMs: calculates LSMs and compares them to the raw means. This makes the Elev_means_Feb2022 ROBJ that is used in scripts 6 and 11.
- 05_Pixy_pi: the code for running pixy (all code also calculated Fst with pixy) on the HPCC
- 06_Pixy_pi_04082022: Multiple regression analysis to test for drift maintaining the cline. Also does the short stamen number elevation cline quadratic regression. *makes fig 2* 
- 07_Fst_diveRsity_2022: makes the Fst heatmap for cent and no cent. analysis cut from manuscript.
- 08_PCA_2022: runs PCA and makes figures for cent and no cent for first 4 PCs. Also does regressions of PCs with elevation. *makes fig 3A, B*
- 09_GWAS: makes all the GWAS manhattan and qqplots for each of the 4 types of GWAS with and without the centromere. Also does the shared SNP analysis.
- 09B_GWAS_HighlightedSnpsFigures: makes some fancy figures and pulls out SNPs near the Royer et al. 2016 QTL. *makes fig 4*
- 09C LD decay: added analysis to calculate LD decay to estimate windows of LD when looking for genes that may be in LD with GWAS hits
- 10_QpcPrep: prepped the files by taking a random subset of 50k SNPs with no missing data
- 11_Qpc: runs the Qpc and makes figures *makes fig 3C*
- PhenotypingInfo: calculates stats for methods paragraphs about sample size
- elevation_map: code that makes fig1 map. Some packages now out of date. *makes fig 1B*


last updated: 3/6/2025
