# Code to calculate coverage statistics from what was calculated with the sorted BWA output
# 10/9/2023
# Sophie Buysse

# the statistics I want to calculate are:
# 1) average coverage across everything
# 2) median across everything
# 3) how many rows even are there? (locations)
# 4) how many rows with at least 1x (percentage)
# 5) how many rows with at least 4x (percentage)
# 6) how many rows with at least 10x (percentage)

# # let's read in the temp file I downloaded.
# 
# PIN9 <- read.delim("C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/FastQC/PIN9_CGAGGCTG-TAGATCGC_sorted_coverage.txt", header=FALSE)
# # name the columns
# colnames(PIN9) <- c("Chr", "Pos", "Coverage")
# # do the 6 things with chloroplast and mitochondria
# avg_all <- mean(PIN9$Coverage)
# med_all <- median(PIN9$Coverage)
# #how many rows?
# tot <- nrow(PIN9)
# one <- nrow(PIN9[PIN9$Coverage >= 1, ])
# four <- nrow(PIN9[PIN9$Coverage >= 4, ])
# ten <- nrow(PIN9[PIN9$Coverage >= 10, ])
# 
# # let's see the percentages
# # at least one
# print(one / tot)
# #94.4
# 
# # at least four
# print (four / tot)
# # 89.2
# 
# # at least ten
# print(ten / tot)
# # 56.7
# 
# # a final data frame I would want could look like:
# df <- data.frame(tot, one, four, ten, avg_all, med_all)
# 
# # do the 6 things without chloroplast and mitochondria
# sub_PIN9 <- PIN9[PIN9$Chr %in% c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), ]
# # then immediately remove the big one for space
# rm(PIN9)
# 
# avg_sub <- mean(sub_PIN9$Coverage)
# med_sub <- median(sub_PIN9$Coverage)
# #how many rows?
# tot_sub <- nrow(sub_PIN9)
# one_sub <- nrow(sub_PIN9[sub_PIN9$Coverage >= 1, ])
# four_sub <- nrow(sub_PIN9[sub_PIN9$Coverage >= 4, ])
# ten_sub <- nrow(sub_PIN9[sub_PIN9$Coverage >= 10, ])
# 
# # let's see the percentages
# # at least one
# print(one_sub / tot_sub)
# #94.4 -> the same
# 
# # at least four
# print (four_sub / tot_sub)
# # 89.14 -> practically the same
# 
# # at least ten
# print(ten_sub / tot_sub)
# # 56.5 -> practically the same
# 
# # a final data frame I would want could look like:
# df_sub <- data.frame(tot_sub, one_sub, four_sub, ten_sub, avg_sub, med_sub)
# 
# tot_sub / tot
# # 0.9956 so this is all basically the same so I just want to get rid of the chloroplast and mitochondria because I think that makes sense because they are then removed during the haplocaller step?
# # based on this example, practically it should not matter... so actually let's just leave it in and make the calculating code a tiny bit easier for me.

############## start here pm HPCC ###################
## load in modules
#module purge
#module load  GCC/11.2.0  OpenMPI/4.1.1 R/4.2.2
#R


## need a loop to do the calculations on each file and then make an output table at the end
# okay so first I want to read in all the files that follow this pattern:
files_all <- list.files(path = "/mnt/scratch/buysseso/BWA_bam/coverage", pattern = "*_sorted_coverage.txt", full.names = TRUE)
# finds an extra ALE10 fine, but I don't think that is a big deal.

# need to figure out how I want to name all the files as I read them in. 
# because I won't really be looking at the files, why not name them just the long string of things before _sorted_coverage?
names <- substr(files_all, 40, 47)
# then write my loop

# run a loop to read in the files to the pre-determined names.
#for (i in c(1:length(names))){
#  assign(names[i], read.delim(files_all[i]))
#}
# well this aborted my R studio sesison twice using an interactive HPCC session, so maybe don't do this.
# even just reading in one file aborted the R studio interactive HPCC session with default memory (750MB) and requesting 5GB memory. Then just went to running
# the session within MobaXterm and loading the R module

# a final data frame I would want could look like:
df <- data.frame("tot" = rep(NA, times = 63), "one" = rep(NA, times = 63), "per_1" = rep(NA, times = 63), "four" = rep(NA, times = 63), "per_4" = rep(NA, times = 63), "ten" = rep(NA, times = 63), "per_10" = rep(NA, times = 63), "avg_all" = rep(NA, times = 63), "med_all" = rep(NA, times = 63))
rownames(df) <- names

for (i in c(1:length(files_all))){
  # track progress
  print(names[i])
  # read in the file
  tmp_file <- read.delim(files_all[i], header = FALSE)
  # and rename the columns on each one to be matching
  colnames(tmp_file) <- c("Chr", "Pos", "Coverage")
  
  # do the 6 things with chloroplast and mitochondria
  avg_all <- mean(tmp_file$Coverage)
  med_all <- median(tmp_file$Coverage)
  #how many rows?
  tot <- nrow(tmp_file)
  one <- nrow(tmp_file[tmp_file$Coverage >= 1, ])
  four <- nrow(tmp_file[tmp_file$Coverage >= 4, ])
  ten <- nrow(tmp_file[tmp_file$Coverage >= 10, ])
  per_1 <- one / tot
  per_4 <- four / tot
  per_10 <- ten / tot
  vals <- c(tot, one, per_1, four, per_4, ten, per_10, avg_all, med_all)
  df[i, ] <- vals
  
  # remove file
  rm(tmp_file)
  # track progress
  print("done")
}

# each one takes about 2 minutes, so I expect the full code to take about 2 hours on my local machine
# it would probably be smart to submit this as a job so it actually finishes,
# but I"m not feeling smart today and just started it locally at 3:27pm on 10/10/23
# done at 5:30

# check out the output
head(df)
tail(df)
# I want to save this dataframe! let's do a .csv instead of an R file
write.csv(df, "/mnt/scratch/buysseso/BWA_bam/coverage/coverage_summary.csv", row.names = TRUE, quote = FALSE)
