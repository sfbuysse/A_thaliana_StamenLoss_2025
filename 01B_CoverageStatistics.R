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

# # let's read in the temp file I downloaded to try stuff out
# 
PIN9 <- read.delim("file/location/PIN9_CGAGGCTG-TAGATCGC_sorted_coverage.txt", header=FALSE)
# # name the columns
colnames(PIN9) <- c("Chr", "Pos", "Coverage")
# # do the 6 things with chloroplast and mitochondria
avg_all <- mean(PIN9$Coverage)
med_all <- median(PIN9$Coverage)
# #how many rows?
tot <- nrow(PIN9)
one <- nrow(PIN9[PIN9$Coverage >= 1, ])
four <- nrow(PIN9[PIN9$Coverage >= 4, ])
ten <- nrow(PIN9[PIN9$Coverage >= 10, ])

# let's see the percentages
# at least one
print(one / tot)
#94.4

# at least four
print (four / tot)
# 89.2

# at least ten
print(ten / tot)
# 56.7

# a final data frame I would want could look like:
df <- data.frame(tot, one, four, ten, avg_all, med_all)

# do the 6 things without chloroplast and mitochondria
sub_PIN9 <- PIN9[PIN9$Chr %in% c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"), ]
# then immediately remove the big one for space
rm(PIN9)

avg_sub <- mean(sub_PIN9$Coverage)
med_sub <- median(sub_PIN9$Coverage)
#how many rows?
tot_sub <- nrow(sub_PIN9)
one_sub <- nrow(sub_PIN9[sub_PIN9$Coverage >= 1, ])
four_sub <- nrow(sub_PIN9[sub_PIN9$Coverage >= 4, ])
ten_sub <- nrow(sub_PIN9[sub_PIN9$Coverage >= 10, ])

# let's see the percentages
# at least one
print(one_sub / tot_sub)
#94.4 -> the same

# at least four
print (four_sub / tot_sub)
# 89.14 -> practically the same

# at least ten
print(ten_sub / tot_sub)
# 56.5 -> practically the same

# a final data frame I would want could look like:
df_sub <- data.frame(tot_sub, one_sub, four_sub, ten_sub, avg_sub, med_sub)

tot_sub / tot
# 0.9956 

#basically the same, so leaving chloroplast and mitochondria in for calculating these statistics

############## do them all! ###################
# load in modules if doing on HPCC
#module purge
#module load  GCC/11.2.0  OpenMPI/4.1.1 R/4.2.2
#R

## need a loop to do the calculations on each file and then make an output table at the end
# read in all files that follow pattern
files_all <- list.files(path = "/file/location/BWA_bam/coverage", pattern = "*_sorted_coverage.txt", full.names = TRUE)

# pull out identifier to use for names
names <- substr(files_all, 40, 47)

# a final data frame I would want could look like:
df <- data.frame("tot" = rep(NA, times = 63), "one" = rep(NA, times = 63), "per_1" = rep(NA, times = 63), "four" = rep(NA, times = 63), "per_4" = rep(NA, times = 63), "ten" = rep(NA, times = 63), "per_10" = rep(NA, times = 63), "avg_all" = rep(NA, times = 63), "med_all" = rep(NA, times = 63))
rownames(df) <- names

# loop
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
  # percentages
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

# check out the output
head(df)
tail(df)
# save this dataframe!
write.csv(df, "/file/location/BWA_bam/coverage/coverage_summary.csv", row.names = TRUE, quote = FALSE)
