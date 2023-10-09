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

## need a loop to do the calculations on each file and then make an output table at the end
# okay so first I want to read in all the files that follow this pattern:
files_all <- list.files(path = "/mnt/scratch/buysseso/BWA_bam/coverage", pattern = "*_sorted_coverage.txt", full.names = TRUE)

############## start here ###################
# need to figure out how I want to name all the files as I read them in. 
# then write my loop

# make list of object names to use (doing this manually but could be automated with some substr code I bet)
names_GS1 <- c("Phenology_1", "Plant_1", "Stratify_1", "Bolting_1", "PreVern_1", "Root_1")
# run a loop to read in the files to the pre-determined names.
for (i in c(1:length(names_GS1))){
  assign(names_GS1[i], read_csv(files_GS1[i], show_col_types = FALSE))
}


PIN9 <- read.delim("C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/FastQC/PIN9_CGAGGCTG-TAGATCGC_sorted_coverage.txt", header=FALSE)
# and rename the columns on each one to be matching
colnames(PIN9) <- c("Chr", "Pos", "Coverage")

# now let's make an empty dataframe
# a final data frame I would want could look like:
df <- data.frame("tot" = NA, "one" = NA, "per_1" = NA, "four" = NA, "per_4" = NA, "ten" = NA, "per_10" = NA, "avg_all" = NA, "med_all" = NA)

# do the 6 things with chloroplast and mitochondria
avg_all <- mean(PIN9$Coverage)
med_all <- median(PIN9$Coverage)
#how many rows?
tot <- nrow(PIN9)
one <- nrow(PIN9[PIN9$Coverage >= 1, ])
four <- nrow(PIN9[PIN9$Coverage >= 4, ])
ten <- nrow(PIN9[PIN9$Coverage >= 10, ])

# let's see the percentages
# at least one
per_1 <- one / tot
#94.4

# at least four
per_4 <- four / tot
# 89.2

# at least ten
per_10 <- ten / tot
# 56.7

# a final data frame I would want could look like:
df <- data.frame(tot, one, four, ten, avg_all, med_all)