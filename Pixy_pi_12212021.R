########## Visualize the Pixy results  ##########


library(dplyr)
library(ggplot2)
library(ggpubr)
####Need to read in the results
###dat <- read.delim("C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/pixy_round2/all_filtered.50k.out_pi.txt", header=TRUE)
###head(dat)
###tail(dat)
###
########## Calculate Population global pi #######
### use dyplyr to group by population
# make a facto
dat$pop <- as.factor(dat$pop)
# calculate the aggregated value
## (window 1 count_diffs + window 2 count_diffs) / (window 1 comparisons + window 2 comparisons)
sums <- dat %>% group_by(pop) %>% summarize(Sum.count.diff = sum(count_diffs, na.rm = T), Sum.count.comp = sum(count_comparisons, na.rm = T))
sums$Mean.pi <- sums$Sum.count.diff/sums$Sum.count.comp
##### this should be right following the pixy documentation
###
#### add the other population information
#### ARB -> ARU
#### SPE -> SAL
#### I must have actually fixed this in previous code so I don't have to do it for every analysis
###elev <- data.frame("pop" = c("ALE","ARU","BAR","BIS","BOS","COC","HOR","MUR","PAL","PAN","PIN","POB","RAB","SAL",
###                            "VDM","VIE"), elev_m = c(1229,416,441,1444,715,519,413,836,1585,1706,177,665,61,303,991,1605))
#### merge the info together
###comp <- merge(sums, elev, by = "pop")
##### add in mean pop ssn as well. just using the raw means because they are highly correlated with the LSMs
###pheno.all <- read.csv("C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Ind_trait_data.csv", stringsAsFactors = FALSE)
###ave2 <- pheno.all %>% group_by(Population) %>% summarize(Mean.ssn = mean(Short_Stamens, na.rm=TRUE))
#### change ARB and SPE to ARU and SAL
###ave2$Population <- c("ALE","ARU","BAR","BIS","BOS","COC","HOR","MUR","PAL","PAN","PIN","POB","RAB","SAL",
###                    "VDM","VIE")
###head(ave2)
#### join the means with other info
###comp2 <- merge(comp, ave2, by.x = "pop", by.y = "Population")
#### save as an R object
###save(comp2, file = "C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/PopGlobalPi_allsite_12212021.ROBJ")

load("C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/PopGlobalPi_allsite_12212021.ROBJ")

##### Correlation/regressions between Variables #######
### Pi and Elevation ###
plot(Mean.pi ~ elev_m, data = comp2)
cor.test(comp2$Mean.pi, comp2$elev_m, method = "pearson")
### strong correlation. r = -0.81 p = 0.0001807
m.pi.elev <- lm(Mean.pi ~ elev_m, data = comp2)
summary(m.pi.elev)
### r2 = 0.6447. p val is significant (but from stats I recall I should not look at this)
### means that elevation does explain pi. 

### SSN and Elevation ###
plot(Mean.ssn ~ elev_m, data = comp2)
cor.test(comp2$Mean.ssn, comp2$elev_m, method = "pearson")
### r = 0.568 p = 0.02165
###plot(ShrtLSM ~ elev_m, data = comp2)
###cor.test(comp2$ShrtLSM, comp2$elev_m, method = "pearson")
### should really be a linear and quadratic regression
### which means I do the quadratic and look at both the linear and quadratic components of the model (coefficients)
### for quadratic regression
comp2$elev_m2 <- (comp2$elev_m)^2
m.ssn.elev <- lm(Mean.ssn ~ elev_m + elev_m2, data = comp2)
summary(m.ssn.elev)
### linear component is sig at 0.1 level. coefficient is positive to the e-03
### quadratic component is not significant. coefficient is negative to the e-07
### neg means it opens downward and small value means it is a pretty flat line. will want to plot this on the graph I make
### r2 = 0.4134
### p value of model is significant?

### Pi and SSN ###
plot(Mean.ssn ~ Mean.pi, data = comp2)
cor.test(comp2$Mean.pi, comp2$
           Mean.ssn, method = "pearson")
### r = -0.2257 p = 0.4
m.pi.ssn <- lm(Mean.ssn ~ Mean.pi, dat = comp2)
summary(m.pi.ssn)
### r2 = 0.051 so not a good fit. pi does not explain short stamen number. 
### this is an interesting result. does not provide evidence for the drift theory because 
### drift would show that at lower pi values (less nucleotide diverity) there would be less short stamen loss.
### the slope is negative, which is the direction we expect
### this graph is later replaced with a multiple regression later on 

### So, both mean.pi and mean.ssn are decently correlated with elevation, but not with each other! 
### interesting.probably makes sense though because pi is genome wide

##### simulate data to draw lines on the plots where elevation is the predictor #######
newdata <- data.frame(elev_m = seq(55,1750, length = 2000), elev_m2 = seq(3700, 2910450, length = 2000))
newdata2 <- data.frame(Mean.pi = seq(0.001,.0051, length = 2000))
### Now use predict to get the predictions
### can't use se.fit = T for mixed effects model? 
m.ssn.elev.pred = predict(m.ssn.elev, newdata = newdata, se.fit = TRUE)
m.pi.elev.pred = predict(m.pi.elev, newdata = newdata, se.fit = TRUE)
m.pi.ssn.pred = predict(m.pi.ssn, newdata = newdata2, se.fit = TRUE)

forplot2 <- data.frame('newd1' = newdata$elev_m,
                       'newd2' = newdata$elev_m2,
                       'newd3' = newdata2,
                       'pred.ssn.el' = m.ssn.elev.pred,
                       'pred.pi.el' = m.pi.elev.pred,
                       'pred.pi.ssn' = m.pi.ssn.pred)

##### fancy plot. #######
### pi by elevation
comp2 <- comp2[order(comp2$elev_m),]
comp2$pop <- as.factor(comp2$pop)
comp2$pop <- factor(comp2$pop, levels = comp2$pop[order(comp2$elev_m)])
str(comp2$pop)
### this looks like the right order
### key is that you need to order the whole sheet AND reorder the factor. BOTH, not just one.
comp2$labels <- paste0(comp2$pop, " - ", comp2$elev_m, "m")

# with a linear fit line of elevation explaining mean pi
ggplot(comp2)+
  geom_point(data=comp2, aes(x=elev_m, y=Mean.pi, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit), linetype = "solid",
            alpha = 0.7, size = 1.25)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit+1.96*pred.pi.el.se.fit), linetype = "solid",
            alpha = 0.3, size = 0.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit-1.96*pred.pi.el.se.fit), linetype = "solid",
            alpha = 0.5, size = 0.75)+
  labs(x= "Elevation (m)", y= "Mean Pi", title = "")+
  scale_color_manual(name = "Population",
                     labels = comp2$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4))))+
  scale_shape_manual(name = "Population",
                     labels = comp2$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 16),
    legend.spacing.y = unit(0.03, "cm"))
ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/PiPerPop_pixy_12212021.png", height = 7, width = 9)

### ssn by elevation
ggplot(comp2)+
  geom_point(data=comp2, aes(x=elev_m, y=Mean.ssn, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit), linetype = "solid",
            alpha = 0.7, size = 1.25)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit+1.96*pred.ssn.el.se.fit), linetype = "solid",
            alpha = 0.3, size = 0.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit-1.96*pred.ssn.el.se.fit), linetype = "solid",
            alpha = 0.5, size = 0.75)+
  labs(x= "Elevation (m)", y= "Mean Short Stamen Number", title = "")+
  scale_color_manual(name = "Population",
                     labels = comp2$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4)), "black"))+ # black not needed unless all is included
  scale_shape_manual(name = "Population",
                     labels = comp2$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4), 0))+ # extra 0 not needed unless all is included
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 16),
    legend.spacing.y = unit(0.03, "cm"),
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 16))
ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/MeanSSNByElev_12212021_quad.png", height = 7, width = 9)
## looks essentially linear which is not what we were expecting

# lowess line instead
ggplot(comp2)+
  geom_point(data=comp2, aes(x=elev_m, y=Mean.ssn, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_smooth(aes(x=elev_m, y=Mean.ssn),stat = "smooth",position = "identity",method = "loess",
    formula = y ~ x + (x**2),
    se = TRUE,na.rm = FALSE)+
  labs(x= "Elevation (m)", y= "Mean Short Stamen Number", title = "")+
  scale_color_manual(name = "Population",
                     labels = comp2$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4)), "black"))+ # black not needed unless all is included
  scale_shape_manual(name = "Population",
                     labels = comp2$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4), 0))+ # extra 0 not needed unless all is included
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 16),
    legend.spacing.y = unit(0.03, "cm"),
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 16))


### SSN by Pi
ggplot(comp2)+
  geom_point(data=comp2, aes(y=Mean.ssn, x=Mean.pi, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_line(dat = forplot2, aes(x = Mean.pi, y = pred.pi.ssn.fit), linetype = "solid",
            alpha = 0.7, size = 1.25)+
  geom_line(dat = forplot2, aes(x = Mean.pi, y = pred.pi.ssn.fit+1.96*pred.pi.ssn.se.fit), linetype = "solid",
            alpha = 0.3, size = 0.75)+
  geom_line(dat = forplot2, aes(x = Mean.pi, y = pred.pi.ssn.fit-1.96*pred.pi.ssn.se.fit), linetype = "solid",
            alpha = 0.5, size = 0.75)+
  labs(y="Mean Short Stamen Number" , x= "Mean Pi", title = "")+
  scale_color_manual(name = "Population",
                     labels = comp2$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4)), "black"))+ # black not needed unless all is included
  scale_shape_manual(name = "Population",
                     labels = comp2$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4), 0))+ # extra 0 not needed unless all is included
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 16),
    legend.spacing.y = unit(0.03, "cm"))
ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/SSNbyPi_12212021.png", height = 7, width = 9)

####### Population windowed pi #######
#the purpose of this is to find regions of elevated pi (kinda a selection test but I removed the centromere instead of doing any stats)
# this is informing which region of the centromere to exclude, though I didn't actually do the excluding until later.
# the pi values I calcualted previously in this code DO NOT exclude the centromere region. Should they?

pops <- c("ALE","ARU","BAR","BIS","BOS","COC","HOR","MUR","PAL","PAN","PIN","POB","RAB","SAL",
          "VDM","VIE")
pi <- dat %>%
  # compute chromosome size in bp
  group_by(chromosome) %>%
  summarise(chr_len=max(window_pos_1)) %>%
  
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  #add this info to the initial data set (so like adding new column and sorting by it)
  left_join(dat, ., by=c("chromosome"="chromosome")) %>%
  
  #add cum position of each SNP
  arrange(chromosome, window_pos_1) %>%
  mutate( psCum=window_pos_1+tot)

axisdf = pi %>% group_by(chromosome) %>% summarize(center=( max(psCum) + min(psCum) ) /2)
varName <- list()
for (i in 1:16){
  varName[[i]] <- ggplot(pi[pi$pop == pops[i],], aes(x=psCum, y= avg_pi)) +
    geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1.3)+
    scale_color_manual(values = rep(c("grey", "black"), 3)) +
    scale_x_continuous( label = axisdf$chromosome, breaks = axisdf$center ) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x="Chromosome", y= "Pi")+
    ggtitle(paste0("Pi in 50kb windows for ", pops[i]))+
    theme_bw() +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
    )
}

all_plot <- ggarrange(varName[[1]], varName[[2]], varName[[3]], varName[[4]], varName[[5]], varName[[6]], varName[[7]], varName[[8]], varName[[9]], varName[[10]], varName[[11]], varName[[12]], varName[[13]], varName[[14]], varName[[15]], varName[[16]], ncol = 4, nrow = 4)
ggsave("C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/50kb_pixy.pi_12212021_B.png", plot = all_plot, width = 20, height = 20)

## for Kieran
Coc_pi <- varName[[6]]

for (i in 1:16){
  varName[[i]] <- ggplot(pi[pi$pop == pops[i],], aes(x=psCum, y= no_sites)) +
    geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1.3)+
    scale_color_manual(values = rep(c("grey", "black"), 3)) +
    scale_x_continuous( label = axisdf$chromosome, breaks = axisdf$center ) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x="Chromosome", y= "Number Sites per Window")+
    ggtitle(paste0("Sites in 50kb windows for ", pops[i]))+
    theme_bw() +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
    )
}

all_plot <- ggarrange(varName[[1]], varName[[2]], varName[[3]], varName[[4]], varName[[5]], varName[[6]], varName[[7]], varName[[8]], varName[[9]], varName[[10]], varName[[11]], varName[[12]], varName[[13]], varName[[14]], varName[[15]], varName[[16]], ncol = 4, nrow = 4)
ggsave("C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/50kb_sites_12212021.png", plot = all_plot, width = 20, height = 20)
Coc_sites <- varName[[6]]

#### look for correlation between the number of sites and the pi value in each window
#### does a correlation make sense? If sites is number of variant sites, a positive correlation makes sense
### not totally sure what to think about this.
ggplot(dat, aes(x = no_sites, y = avg_pi))+
  geom_point(aes(col = pop))
ggplot(dat, aes(x = count_diffs, y = avg_pi))+
  geom_point(aes(col = pop))
ggplot(dat, aes(x = count_comparisons, y = avg_pi))+
  geom_point(aes(col = pop))
ggplot(dat, aes(x = count_missing, y = avg_pi))+
  geom_point(aes(col = pop))

## just COC pop for Kieran
ggplot(pi[pi$pop == 'COC', ], aes(x=psCum, y= avg_pi)) +
  geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1.3)+
  scale_color_manual(values = rep(c("grey", "black"), 3)) +
  scale_x_continuous( label = axisdf$chromosome, breaks = axisdf$center ) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.11)) +
  labs(x="Chromosome", y= "Pi")+
  ggtitle(paste0("Pi in 50kb windows for COC"))+
  theme_bw() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  )
# end of document
# last ran on 12/21/2021
# thoughts: do I need to exclude the centromere region and calculate pi and run all the analyses again? 
# I think that makes sense in case the centromere region is shifting the global pi calculation.

########## Remove Centromere ##########
require(dplyr)
require(ggplot2)
require(ggpubr)
#Need to read in the results

##### pi calculated in December, 2021
####dat <- read.delim("C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/pixy_round2/all_filtered.50k.out_pi.txt", header=TRUE)
####head(dat)
####tail(dat)
####
####dat$pop <- as.factor(dat$pop)
####
##### first, want to remove the centromere sections published by Clark et al. 2007 and see if the pi to elevation correlation changes
##### need to exclude a different region for each chromosome, so maybe split up and then join back together?
####
###### and filter based on no_sites
####min_sites <- 30000
##### of 50,000 so this means that 3/5 of the sites in the window are present. maybe a little high?
####
####Chr1 <- dat[dat$chromosome == 'Chr1', ]
####head(Chr1)
####tail(Chr1)
##### centromere is 13.7 to 15.9 MB which is position 13700000 to 15900000
####start <- 13700000 # Clark et al published num
#####start <- 12200000 #wide
####stop <- 15900000 # Clark et al published num
#####stop <- 17200000 #wide
####Chr1 <- Chr1[which(Chr1$window_pos_2 < start & Chr1$no_sites > min_sites | Chr1$window_pos_1 > stop & Chr1$no_sites > min_sites), ]
#####Chr1.test <- subset(Chr1, Chr1$window_pos_1 < start | Chr1$window_pos_2 > stop, select = c(pop:count_missing))
##### this does the same thing
####
##### centromere is 2.45 - 5.5 MB
####Chr2 <- dat[dat$chromosome == 'Chr2', ]
####start <- 2450000 # Clark et al published num
#####start <- 1450000 # wide
####stop <- 5500000 # Clark et al published num
#####stop <- 6800000 # wide
####Chr2 <- Chr2[which(Chr2$window_pos_2 < start & Chr2$no_sites > min_sites | Chr2$window_pos_1 > stop & Chr2$no_sites > min_sites), ]
####
##### centromere is 11.3 to 14.3 MB
####Chr3 <- dat[dat$chromosome == 'Chr3', ]
####start <- 11300000 # Clark et al published num
#####start <- 10300000 #wide
####stop <- 14300000 # Clark et al published num
#####stop <- 16800000 #wide
####Chr3 <- Chr3[which(Chr3$window_pos_2 < start & Chr3$no_sites > min_sites | Chr3$window_pos_1 > stop & Chr3$no_sites > min_sites), ]
####
##### centromere is 1.8-5.15MB
####Chr4 <- dat[dat$chromosome == 'Chr4', ]
####start <- 1800000 # Clark et al published num
#####start <- 1300000 #wide
####stop <- 5150000 # Clark et al published num
#####stop <- 6150000 #wide
####Chr4 <- Chr4[which(Chr4$window_pos_2 < start & Chr4$no_sites > min_sites | Chr4$window_pos_1 > stop & Chr4$no_sites > min_sites), ]
####
##### centromere is 11 to 13.35 MB
####Chr5 <- dat[dat$chromosome == 'Chr5', ]
####start <- 11000000 # Clark et al published num
#####start <- 10000000
####stop <- 13350000 # Clark et al published num
#####stop <- 15350000
####Chr5 <- Chr5[which(Chr5$window_pos_2 < start & Chr5$no_sites > min_sites | Chr5$window_pos_1 > stop & Chr5$no_sites > min_sites), ]
####
##### now rbind them all back together
####dat.clean <- rbind(Chr1, Chr2, Chr3, Chr4, Chr5)
####save(dat.clean, file = "C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/pi_noCent_50k_01192022.ROBJ")
load("C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/pi_noCent_50k_01192022.ROBJ")

######## Calculate population global pi ########
# use dyplyr to group by population
## to calculate aggregated value
## (window 1 count_diffs + window 2 count_diffs) / (window 1 comparisons + window 2 comparisons)

sums.c <- dat.clean %>% group_by(pop) %>% summarize(Sum.count.diff = sum(count_diffs, na.rm = T), Sum.count.comp = sum(count_comparisons, na.rm = T))
sums.c$Mean.pi <- sums.c$Sum.count.diff/sums.c$Sum.count.comp
# this should be right following the pixy documentation
sum.comparison <- data.frame('all' = sums$Mean.pi, 'noCent' = sums.c$Mean.pi)
ggplot(sum.comparison)+
  geom_point(aes(x = all, y = noCent))+
  geom_abline(aes(intercept = 0, slope = 1))+
  theme_classic()

# interesting. same patterns but the no cent values are all shifted lower by about the same amount from a 1:1 line.
# maybe more of a reduction in the pops with higher pi generally

elev <- data.frame("pop" = c("ALE","ARU","BAR","BIS","BOS","COC","HOR","MUR","PAL","PAN","PIN","POB","RAB","SAL",
                             "VDM","VIE"), elev_m = c(1229,416,441,1444,715,519,413,836,1585,1706,177,665,61,303,991,1605))

comp.c <- merge(sums.c, elev, by = "pop")

### add in mean pop ssn as well
#pheno.all <- read.csv("C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Ind_trait_data.csv", stringsAsFactors = FALSE)
#ave2 <- pheno.all %>% group_by(Population) %>% summarize(Mean.ssn = mean(Short_Stamens, na.rm=TRUE))
#ave2$Population <- c("ALE","ARU","BAR","BIS","BOS","COC","HOR","MUR","PAL","PAN","PIN","POB","RAB","SAL","VDM","VIE")
#head(ave2)
#
# use ave2 from the everything analysis. currently uses raw means and not ssn.

#comp2.c <- merge(comp.c, ave2, by.x = "pop", by.y = "Population")
#save(comp2.c, file = "C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/PopGlobalPi_NoCent_01192022.ROBJ")
load("C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/PopGlobalPi_NoCent_01192022.ROBJ")

##### Correlations #####
# these are correlations, not regressions like I am doing above.
plot(Mean.pi ~ elev_m, data = comp2.c)
cor.test(comp2.c$Mean.pi, comp2.c$elev_m, method = "pearson")

plot(Mean.ssn ~ elev_m, data = comp2.c)
cor.test(comp2.c$Mean.ssn, comp2.c$elev_m, method = "pearson")
# 0.56

plot(Mean.pi ~ Mean.ssn, data = comp2.c)
cor.test(comp2.c$Mean.pi, comp2.c$Mean.ssn, method = "pearson")

# fancy plot.
## pi
comp2.c <- comp2.c[order(comp2.c$elev_m),]
comp2.c$pop <- as.factor(comp2.c$pop)
comp2.c$pop <- factor(comp2.c$pop, levels = comp2.c$pop[order(comp2.c$elev_m)])
str(comp2.c$pop)
# this looks like the right order
# key is that you need to order the whole sheet AND reorder the factor. BOTH, not just one.
comp2.c$labels <- paste0(comp2.c$pop, " - ", comp2.c$elev_m, "m")

ggplot(comp2.c)+
  geom_point(data=comp2.c, aes(x=elev_m, y=Mean.pi, shape = as.factor(pop), col= as.factor(pop)), size = 4, stroke = 1.25)+
  labs(x= "Elevation (m)", y= "Mean Pi", title = "")+
  scale_color_manual(name = "Population",
                     labels = comp2.c$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4)), "black"))+ # black not needed unless all is included
  scale_shape_manual(name = "Population",
                     labels = comp2.c$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4), 0))+ # extra 0 not needed unless all is included
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 12),
    legend.spacing.y = unit(0.01, "cm"),
    axis.title = element_text(color = "black", size = 14),
    axis.text = element_text(color = "black", size = 12))
#ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/PiPerPop_NoCent_50k_pixy_01192022.png", height = 7, width = 9)

##ssn by mean.pi
ggplot(comp2.c)+
  geom_point(data=comp2.c, aes(x=Mean.ssn, y=Mean.pi, shape = as.factor(pop), col= as.factor(pop)), size = 4, stroke = 1.25)+
  labs(x= "Mean Short Stamen Number", y= "Mean Pi", title = "")+
  scale_color_manual(name = "Population",
                     labels = comp2.c$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4)), "black"))+ # black not needed unless all is included
  scale_shape_manual(name = "Population",
                     labels = comp2.c$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4), 0))+ # extra 0 not needed unless all is included
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 12),
    legend.spacing.y = unit(0.01, "cm"),
    axis.title = element_text(color = "black", size = 14),
    axis.text = element_text(color = "black", size = 12))
#ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/SSNBYPi_NoCent_50k_pixy_01192022.png", height = 7, width = 9)

# windowed to check that it looks like I am expecting it to. The values are lower here if I look at the y axis scale but the trends remain the same.
##### Population windowed pi #####
pi <- dat.clean

pops <- c("ALE","ARU","BAR","BIS","BOS","COC","HOR","MUR","PAL","PAN","PIN","POB","RAB","SAL",
          "VDM","VIE")

don <- pi %>%
  # compute chromosome size in bp
  group_by(chromosome) %>%
  summarise(chr_len=max(window_pos_1)) %>%
  
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  #add this info to the initial data set (so like adding new column and sorting by it)
  left_join(pi, ., by=c("chromosome"="chromosome")) %>%
  
  #add cum position of each SNP
  arrange(chromosome, window_pos_1) %>%
  mutate( psCum=window_pos_1+tot)

axisdf = don %>% group_by(chromosome) %>% summarize(center=( max(psCum) + min(psCum) ) /2)
varName <- list()
for (i in 1:16){
  varName[[i]] <- ggplot(don[don$pop == pops[i],], aes(x=psCum, y= avg_pi)) +
    geom_point( aes(color=as.factor(chromosome)), alpha=0.8, size=1.3)+
    scale_color_manual(values = rep(c("grey", "black"), 3)) +
    scale_x_continuous( label = axisdf$chromosome, breaks = axisdf$center ) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(don$avg_pi, na.rm = T)+0.002)) +
    labs(x="Chromosome", y= "Pi")+
    ggtitle(paste0("Pi in 50kb windows for ", pops[i]))+
    theme_bw() +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
    )
}

all_plot <- ggarrange(varName[[1]], varName[[2]], varName[[3]], varName[[4]], varName[[5]], varName[[6]], varName[[7]], varName[[8]], varName[[9]], varName[[10]], varName[[11]], varName[[12]], varName[[13]], varName[[14]], varName[[15]], varName[[16]], ncol = 4, nrow = 4)
ggsave("C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/50kb_NoCent_pixy_01192022.png", plot = all_plot, width = 20, height = 20)

