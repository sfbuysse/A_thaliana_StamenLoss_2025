# Athaliana_Stamens
Making full pipeline for short stamen loss project

The number in front of each file indicates what order in the pipeline it is. When there are doubles, the '_all my code' version will have a whole bunch of things that didn't work and the other verision will be edited down to be a more usable version that represents the final pipeline for the project.

The sequence information will be publically available. Phenotypic information is in the data folder, as are GWAS output file. Most intermediates files are not shared here because they are large, but can be shared upon request.

# File Descriptions:
- 01_FastQ_to_BAM: Starts with fastq files. The end result is aligned to the reference genome but not cleaned at all.
- 01B_Coverage Statistics: calculating coverage statistics for the supplemental table
- 02_CallVariants: The end result is the raw allsites file of genotypes.
- 03_FilterVariants: is the code to go from raw to filtered variant and invariant sites before doing pixy (vcf) and to make the filtered variant sites file used for Fst (later cut from analysis), GWAS, and PCA (plink format). It also removes the centromere regions post filtering.
- 04_PopulationRawMeans_LSMs: calculates LSMs and compares them to the raw means. This makes the Elev_means_Feb2022 ROBJ that is used in scripts 6 and 11.
- 05_Pixy_pi: the code for running pixy (all code also calculated Fst with pixy) on the HPCC
- 06_Pixy_pi_04082022: Multiple regression analysis to test for drift maintaining the cline. ALso does the short stamen numbet elevation cline quadratic regression. 
- 07_Fst_diveRsity_2022: makes the Fst heatmap for cent and no cent. analysis cut from manuscript
- 08_PCA_2022: runs PCA and makes figures for cent and no cent for first 4 PCs. Also does regressions of PCs with elevation.
- 09_GWAS: makes all the GWAS manhattan and qqplots for each of the 4 types of GWAS with and without the centromere. Also does the shared SNP analysis.
- 09B_GWAS_HighlightedSnpsFigures: makes some fancy figures and pulls out SNPs near the Royer et al. 2016 QTL.
- 10_QpcPrep: prepped the files
- 11_Qpc: runs the Qpc and makes figures
- PhenotypingInfo: calculates stats for methods paragraphs about sample size
- elevation_map: code that makes fig1 map. Some packages now out of date.
- LD decay: added analysis to calculate LD decay to estimate windows of LD when looking for genes that may be in LD with GWAS hits
