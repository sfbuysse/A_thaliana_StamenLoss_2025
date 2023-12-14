# Athaliana_Stamens
Making full pipeline for short stamen loss project

This code is a bit of a mess at the moment. The number in front of each file indicates what order in the pipeline it is. When there are doubles, the '_all my code' version will have a whole bunch of things that didn't work and the other verision will be edited down to be a more sharable version.

I also apparently have old code, so pay attention to the date espcially for pixy. some things are also .Rmd files that have been knitted, so they have accompanying .html files but I"m not sure when the last time I actually knitted it is :)

In general this project needs new organization so that it is synced with the HPCC and the code is in a similar place to all the data, but currently the code is here and the data is there and I do lots of copying and pasting... room for workflow improvement!

Most of the data is in my scratch space, but the most important things should be backed up to the josephslab directory, but some things are backed up twice and it isn't super organized. 

# Notes on Manuscript Figures:
- 01_FastQ_to_BAM: not currently any figures from this code. The end result is aligned to the reference genome but not cleaned at all.
- 02_CallVariants: not currently any figures. The end result is the raw allsites file of genotypes.
- 03_FilterVariants:not currently any figures. is the code to go from raw to filtered variant and invariant sites before doing pixy (vcf) and to make the filtered variant sites file used for Fst, GWAS, and PCA (plink format). It also removes the centromere regions post filtering.
- 04_PopulationRawMeans_LSMs: no figures. calculates LSMs and compares them to the raw means (so the correlation values in the methods section come from here). This makes the Elev_means_Feb2022 ROBJ that is used in scripts 6 and 11.  there is another file with an identical name besides the 04 - I likely made this file with another R script and then renamed it.
- 05_Pixy_pi: no figures, just the bash code for running pixy (all code also calculated Fst with pixy)
- 06_Pixy_pi_04082022: many figures! pi by elevation (cent and no cent - fig 2), ssn by elevation, SSN by pi (cent and no cent - fig 3), both residual multiple regression figures (no cent only - fig 3) - this one has an old code file as well that I no longer use.
- 07_Fst_diveRsity_2022: makes the Fst heatmap for cent and no cent. only in supplemental, not main manuscript
- 08_PCA_2022: runs PCA and makes figures for cent and no cent for first 4 PCs (fig 2). Also does regressions of PCs with elevation.
- 09_GWAS: makes all the GWAS manhattan and qqplots for each of the 4 types of GWAS with and without the centromere. Also does the shared SNP analysis.
- 09_GWAS_HighlightedSnpsFigures: makes some fancy figures and pulls out SNPs near the Royer et al. 2016 QTL.
- 10_QpcPrep: prepped the files, no figures
- 11_Qpc: runs the Qpc and makes figures
- PhenotypingInfo: no figures, calculates stats for methods paragraphs about sample size
- Elev_Means_Feb2022.ROBJ: just a table with elevation, mean SSN for all flowers, and mean SSN for sequenced lines.
- SampleSequencingInfo_052022: attempt at getting sample sequencing info together. not sure I'll keep this but it will be a good starting point for getting the sample sequencing table together. This is not used in any genotype quality info presented in table 1 as of 12/14/2023
- elevation_map: code only makes the map that is figure 1