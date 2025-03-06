########## Visualize the Pixy results  ##########
library(dplyr)
library(ggplot2)
library(ggpubr)

####Need to read in the results
dat <- read.delim("~/pixy_042022/all_filtered.50k_pi_042022.out_pi.txt", header=TRUE)
head(dat)
tail(dat)

########## with Centromeres ##########
#most of this just goes in the supplement.
###### Calculate Population global pi ######
##### use dyplyr to group by population
#### make a factor
dat$pop <- as.factor(dat$pop)
#### calculate the aggregated value
#### (window 1 count_diffs + window 2 count_diffs) / (window 1 comparisons + window 2 comparisons)
sums <- dat %>% group_by(pop) %>% summarize(Sum.count.diff = sum(count_diffs, na.rm = T), Sum.count.comp = sum(count_comparisons, na.rm = T))
sums$Mean.pi <- sums$Sum.count.diff/sums$Sum.count.comp
####### this follows the pixy documentation
###### add the other population information
###### Should have labels ARU and SAL
###### this has both the full and only the sequenced means.
load(file = "~/R_script/StamenLossPipeline/Elev_Means_Feb2022.ROBJ")
###### merge the info together
comp <- merge(sums, Elev_Means, by.x = "pop", by.y = "Population")


###### save as an R object
save(comp, file = "~/R_script/PopGlobalPi_allsite_04082022.ROBJ")

load("~/R_script/PopGlobalPi_allsite_04082022.ROBJ")


###### Correlation/regressions between Variables ######
### Pi and Elevation ###
plot(Mean.pi ~ Elev_m, data = comp)
cor.test(comp$Mean.pi, comp$Elev_m, method = "pearson")
### strong correlation. r = -0.801 p = 0.0001942
m.pi.elev <- lm(Mean.pi ~ Elev_m, data = comp)
summary(m.pi.elev)
### r2 = 0.6411. p val is significant (but from stats I recall I should not look at this)
### means that elevation does explain pi. 

# adding this on later date. what about quadratic elevation predicting pi?
comp$Elev_c2 <- (comp$Elev_m - mean(comp$Elev_m))^2
m.pi.elev2 <- lm(Mean.pi ~ Elev_m + Elev_c2, data = comp)
summary(m.pi.elev2)
#slightly higher r2 than linear.

### SSN and Elevation ###
## Do this with the Full means because we want the best picture of each popualtion
plot(Full_PopFlwrMean ~ Elev_m, data = comp)
cor.test(comp$Full_PopFlwrMean, comp$Elev_m, method = "pearson")
### r = 0.618 p = 0.01074

### quadratic regression and look at both the linear and quadratic components of the model (coefficients)
### for quadratic regression
comp$Elev_m2 <- (comp$Elev_m)^2
m.ssn.elev <- lm(Full_PopFlwrMean ~ Elev_m + Elev_m2, data = comp)
summary(m.ssn.elev)
### linear component is sig at 0.05 level. coefficient is positive to the e-03
### quadratic component is not significant. coefficient is negative to the e-07
### neg means it opens downward and small value means it is a pretty flat line. will want to plot this on the graph I make
### r2 = 0.4818
### p value of model is significant?
# calculated a y max (using -b/2a formula) at 1,330.372250423011844331641285956

#ggplot code to test sanity
#ggplot(comp, aes(x = Elev_m, y = Full_PopFlwrMean)) + geom_point() +
#  stat_smooth(aes(x = Elev_m, y = Full_PopFlwrMean), method = "lm", formula = y ~ x, colour = "red") +
#  stat_smooth(aes(x = Elev_m, y = Full_PopFlwrMean), method = "lm", formula = y ~ poly(x, 2), colour = "blue")

#center quadratic values to decrease colinearity
comp$Elev_c2 <- (comp$Elev_m - mean(comp$Elev_m))^2
tmp4 <- lm(Full_PopFlwrMean ~ Elev_m + Elev_c2, data = comp)
summary(tmp4)

# try just a linear fit 
m.ssn.elev.linear <- lm(Full_PopFlwrMean ~ Elev_m, data = comp)
summary(m.ssn.elev.linear)

### Pi and SSN ###
## here only use the sequenced lines because that is what pi was calculated from
plot(Seq_PopFlwrMean ~ Mean.pi, data = comp)
cor.test(comp$Mean.pi, comp$
           Seq_PopFlwrMean, method = "pearson")
### r = -0.2323357 p = 0.3865
m.pi.ssn <- lm(Seq_PopFlwrMean ~ Mean.pi, dat = comp)
summary(m.pi.ssn)
### r2 = 0.054 so not a good fit. pi does not explain short stamen number. 
### this graph is later replaced with a multiple regression later on 

### So, both mean.pi and mean.ssn are decently correlated with elevation, but not with each other! 
### interesting.probably makes sense though because pi is genome wide

# adding question of if quadratic nucleotide diveristy would get a better model
comp$pi_c2 <- (comp$Mean.pi - mean(comp$Mean.pi))^2
m.pi.ssn2 <- lm(Seq_PopFlwrMean ~ Mean.pi + pi_c2, dat = comp)
summary(m.pi.ssn2)

### The two means ###
cor.test(comp$Full_PopFlwrMean, comp$Seq_PopFlwrMean, method = "pearson")
# r = 0.985 p < 0.001

###### simulate data to draw lines on the plots where elevation is the predictor ######
new_elev <- seq(55,1750, length = 2000)
center_val <- mean(comp$Elev_m)
newdata <- data.frame(Elev_m = new_elev, Elev_m2 = (new_elev**2), Elev_c2 = ((new_elev - center_val)**2))
newdata2 <- data.frame(Mean.pi = seq(0.001,.0051, length = 2000))
### Now use predict to get the predictions
m.ssn.elev.pred = predict(tmp4, newdata = newdata, se.fit = TRUE)
m.pi.elev.pred = predict(m.pi.elev, newdata = newdata, se.fit = TRUE)
m.pi.ssn.pred = predict(m.pi.ssn, newdata = newdata2, se.fit = TRUE)

forplot2 <- data.frame('newd1' = newdata$Elev_m,
                       'newd3' = newdata2,
                       'pred.ssn.el' = m.ssn.elev.pred,
                       'pred.pi.el' = m.pi.elev.pred,
                       'pred.pi.ssn' = m.pi.ssn.pred)

###### fancy plot. ######
### pi by elevation
comp <- comp[order(comp$Elev_m),]
comp$pop <- as.factor(comp$pop)
comp$pop <- factor(comp$pop, levels = comp$pop[order(comp$Elev_m)])
str(comp$pop)
comp$labels <- paste0(comp$pop, " - ", comp$Elev_m, "m")

# for topo color scheme
barplot(rep(1,6), col = topo.colors(6))
barplot(rep(1,10), col = topo.colors(10))
barplot(rep(1,16), col = topo.colors(16))

# with a linear fit line of elevation explaining mean pi
fig2b <- ggplot(comp)+
  geom_point(data=comp, aes(x=Elev_m, y=Mean.pi, shape = as.factor(pop), fill= Elev_m), col = "black", size = 3, stroke = 1, show.legend = FALSE)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit), linetype = "solid",
            alpha = 0.9, linewidth = 1.)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit+1.96*pred.pi.el.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.65)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit-1.96*pred.pi.el.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.65)+
  annotate(geom = "text", x = -Inf, y = -Inf, label = "~beta==-1.81e10^-6 *';' ~italic(p)~ bold(' < 0.001')", parse = TRUE,
           hjust = -.01, vjust = -1.5, size = 3, fontface = "bold")+
  annotate(geom = "text", x = -Inf, y = -Inf, label = "~italic(r)^2 == '0.64'", parse = TRUE,
           hjust = -.01, vjust = -0.5, size = 3)+
  labs(x= "Elevation (m)", y= "Nucleotide Diversity", title = "")+
  scale_fill_gradientn(name = "Elevation", colours = topo.colors(16))+
  scale_shape_manual(name = "Population",
                     labels = comp$label,
                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 10),
    axis.text = element_text(color = "black", size = 10),
    legend.spacing.y = unit(0.03, "cm"))

fig2b <- annotate_figure(fig2b,
                fig.lab = "B", fig.lab.face = "bold")

ggsave(filename ="~/Figures/ManuscriptFigs/PiPerPop_topo.png",
       height = 3, width = 3, device = "png", dpi = 700)

### ssn by elevation
fig2a <- ggplot(comp)+
  geom_point(data=comp, aes(x=Elev_m, y=Full_PopFlwrMean, shape = as.factor(pop), fill = Elev_m), col = "black", size = 3, stroke = 1, show.legend = FALSE)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit), linetype = "solid",
            alpha = 0.9, linewidth = 1)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit+1.96*pred.ssn.el.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.65)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit-1.96*pred.ssn.el.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.65)+
  annotate(geom = "text", x = -Inf, y = -Inf, label = "~beta==4.82e10^-4 *';' ~italic(p)~ bold(' = 0.004')", parse = TRUE,
           hjust = -0.6, vjust = -2.5, size = 3, fontface = "bold")+
  annotate(geom = "text", x = -Inf, y = -Inf, label = "gamma==-4.73e10^-7 *';' ~italic(p)~ ' = 0.137'", parse = TRUE,
           hjust = -0.6, vjust = -1.3, size = 3)+
  annotate(geom = "text", x = -Inf, y = -Inf, label = "italic(r)^2 == '0.48'", parse = TRUE,
           hjust = -1.95, vjust = -0.5, size = 3)+
  labs(x= "Elevation (m)", y= "Mean Short Stamen Number", title = "")+ # - Full Phenotype Set
  scale_fill_gradientn(name = "Elevation", colours = topo.colors(10))+
  scale_shape_manual(name = "Population",
                     labels = comp$label,
                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 12),
    legend.spacing.y = unit(0.03, "cm"),
    axis.title = element_text(color = "black", size = 10),
    axis.text = element_text(color = "black", size = 10))
fig2a <- annotate_figure(fig2a,
                         fig.lab = "A", fig.lab.face = "bold")
ggsave(filename ="~/Figures/ManuscriptFigs/MeanSSNByElev_topo.png",
       height = 3, width = 3, device = "png", dpi = 700)

### SSN by Pi
ggplot(comp)+
  geom_point(data=comp, aes(y=Seq_PopFlwrMean, x=Mean.pi, shape = as.factor(pop), fill = Elev_m), col = "black", size = 5, stroke = 1.75)+
  geom_line(dat = forplot2, aes(x = Mean.pi, y = pred.pi.ssn.fit), linetype = "solid",
            alpha = 0.7, size = 1.25)+
  geom_line(dat = forplot2, aes(x = Mean.pi, y = pred.pi.ssn.fit+1.96*pred.pi.ssn.se.fit), linetype = "solid",
            alpha = 0.3, size = 0.75)+
  geom_line(dat = forplot2, aes(x = Mean.pi, y = pred.pi.ssn.fit-1.96*pred.pi.ssn.se.fit), linetype = "solid",
            alpha = 0.5, size = 0.75)+
  labs(y="Mean Short Stamen Number" , x= "Nucleotide Diversity", title = "")+ # - Sequenced Only
  scale_fill_gradientn(name = "Elevation", colours = topo.colors(10))+
  scale_shape_manual(name = "Population",
                     labels = comp$label,
                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20),
    axis.title = element_text(color = "black", size = 20),
    axis.text = element_text(color = "black", size = 20),
    legend.spacing.y = unit(0.03, "cm"))
ggsave(filename ="~/Figures/ManuscriptFigs/SSNbyPi_topo.png",
       height = 7, width = 9, device = "png", dpi = 500)

# Compare the FUll and Seq Means
ggplot(comp)+
  geom_point(data=comp, aes(y=Seq_PopFlwrMean, x=Full_PopFlwrMean, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_abline(intercept = 0, slope = 1, size = 0.75)+
  labs(y="Mean Short Stamen Number - Sequenced Only" , x= "Mean Short Stamen Number - Full", title = "")+
  scale_color_manual(name = "Population",
                     labels = comp$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4)), "black"))+ # black not needed unless all is included
  scale_shape_manual(name = "Population",
                     labels = comp$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4), 0))+ # extra 0 not needed unless all is included
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 16),
    legend.spacing.y = unit(0.03, "cm"))
ggsave(filename ="~/Figures/MeanComparison_02022022.png", height = 7, width = 9)

########## Without Centromeres ##########
require(dplyr)
require(ggplot2)
require(ggpubr)
#Need to read in the results

#### pi calculated in December, 2021
###dat <- read.delim("~/pixy_042022/all_filtered.50k_pi_042022.out_pi.txt", header=TRUE)
###head(dat)
###tail(dat)
###
###dat$pop <- as.factor(dat$pop)
###
#### first, want to remove the centromere sections published by *Clark et al. 2007* and see if the pi to elevation correlation changes
#### need to exclude a different region for each chromosome, so maybe split up and then join back together?
###
##### and filter based on no_sites
###min_sites <- 30000
#### of 50,000 so this means that 3/5 of the sites in the window are present. maybe a little high?
###
###Chr1 <- dat[dat$chromosome == 'Chr1', ]
####9744 observations
###head(Chr1)
###tail(Chr1)
#### centromere is 13.7 to 15.9 MB which is position 13700000 to 15900000
###start <- 13700000 # Clark et al published num
####start <- 12200000 #wide
###stop <- 15900000 # Clark et al published num
####stop <- 17200000 #wide
###Chr1 <- Chr1[which(Chr1$window_pos_2 < start & Chr1$no_sites > min_sites | Chr1$window_pos_1 > stop & Chr1$no_sites > min_sites), ]
####Chr1.test <- subset(Chr1, Chr1$window_pos_1 < start | Chr1$window_pos_2 > stop, select = c(pop:count_missing))
#### this does the same thing
#### 8599 observations
###
#### centromere is 2.45 - 5.5 MB
###Chr2 <- dat[dat$chromosome == 'Chr2', ]
####6304 observations
###start <- 2450000 # Clark et al published num
####start <- 1450000 # wide
###stop <- 5500000 # Clark et al published num
####stop <- 6800000 # wide
###Chr2 <- Chr2[which(Chr2$window_pos_2 < start & Chr2$no_sites > min_sites | Chr2$window_pos_1 > stop & Chr2$no_sites > min_sites), ]
#### 5085 observations
###
###
#### centromere is 11.3 to 14.3 MB
###Chr3 <- dat[dat$chromosome == 'Chr3', ]
#### 7520 observations
###start <- 11300000 # Clark et al published num
####start <- 10300000 #wide
###stop <- 14300000 # Clark et al published num
####stop <- 16800000 #wide
###Chr3 <- Chr3[which(Chr3$window_pos_2 < start & Chr3$no_sites > min_sites | Chr3$window_pos_1 > stop & Chr3$no_sites > min_sites), ]
#### 6108 observations
###
#### centromere is 1.8-5.15MB
###Chr4 <- dat[dat$chromosome == 'Chr4', ]
#### 5952 observations
###start <- 1800000 # Clark et al published num
####start <- 1300000 #wide
###stop <- 5150000 # Clark et al published num
####stop <- 6150000 #wide
###Chr4 <- Chr4[which(Chr4$window_pos_2 < start & Chr4$no_sites > min_sites | Chr4$window_pos_1 > stop & Chr4$no_sites > min_sites), ]
#### 4684 observations
###
#### centromere is 11 to 13.35 MB
###Chr5 <- dat[dat$chromosome == 'Chr5', ]
####8640 observations
###start <- 11000000 # Clark et al published num
####start <- 10000000
###stop <- 13350000 # Clark et al published num
####stop <- 15350000
###Chr5 <- Chr5[which(Chr5$window_pos_2 < start & Chr5$no_sites > min_sites | Chr5$window_pos_1 > stop & Chr5$no_sites > min_sites), ]
#### 7324 observations
###
#### now rbind them all back together
###dat.clean <- rbind(Chr1, Chr2, Chr3, Chr4, Chr5)
###save(dat.clean, file = "~/R_script/pi_noCent_50k_04082022.ROBJ")
load("~/R_script/pi_noCent_50k_04082022.ROBJ")

###### Calculate population global pi for centromere trimmed data ######
# use dyplyr to group by population
## to calculate aggregated value
## (window 1 count_diffs + window 2 count_diffs) / (window 1 comparisons + window 2 comparisons)

# final object is saved and loaded below.
sums.c <- dat.clean %>% group_by(pop) %>% summarize(Sum.count.diff = sum(count_diffs, na.rm = T), Sum.count.comp = sum(count_comparisons, na.rm = T))
sums.c$Mean.pi.c <- sums.c$Sum.count.diff/sums.c$Sum.count.comp
# this should be right following the pixy documentation
sum.comparison <- data.frame('all' = sums$Mean.pi, 'noCent' = sums.c$Mean.pi.c, 'pop' = sums$pop, 'pop.c' = sums.c$pop)
cor.test(sum.comparison$all, sum.comparison$noCent, method = "pearson")
# r = 0.9981 p < 0.001

ggplot(data = sum.comparison)+
  geom_point(aes(x = all, y = noCent))+
  geom_abline(aes(intercept = 0, slope = 1))+
  theme_classic()

# this section is a little out of order b/c I need the labels from comp2.c dataframe
### something went wrong here and the coloring doesn't seem to  match up. probably because of the factor ordering I do for comp that I didn't do here.
ggplot(sum.comparison)+
  geom_point(aes(x = all, y = noCent, shape = as.factor(pop.c), col= as.factor(pop.c)), size = 4, stroke = 1.25)+
  geom_abline(intercept = 0, slope = 1)+
  labs(x= "All Windows Pi", y= "Centromere Excluded Pi", title = "")+
  scale_color_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4)), "black"))+ # black not needed unless all is included
  scale_shape_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4), 0))+ # extra 0 not needed unless all is included
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 12),
    legend.spacing.y = unit(0.01, "cm"),
    axis.title = element_text(color = "black", size = 14),
    axis.text = element_text(color = "black", size = 12))

# interesting. same patterns but the no cent values are all shifted lower by about the same amount from a 1:1 line.
# maybe more of a reduction in the pops with higher pi generally

# add in elevation and stamen means
#load(file = "~/R_script/StamenLossPipeline/Elev_Means_Feb2022.ROBJ")

#comp.c <- merge(sums.c, Elev_Means, by.x = "pop", by.y = "Population")

# check how similar the pi values are
cor(x=comp$Mean.pi, y=comp.c$Mean.pi.c)
# 0.9981597

#save(comp.c, file = "~/R_script/PopGlobalPi_NoCent_04082022.ROBJ")
load("~/R_script/PopGlobalPi_NoCent_04082022.ROBJ")

###### Correlations with no centromere values and variables ######
# these are correlations, not regressions like I am doing above.
# on 8/7/2023 I am adding in the new regressions (or some are already here.)

## Mean pi and Elevation
plot(Mean.pi.c ~ Elev_m, data = comp.c)
cor.test(comp.c$Mean.pi.c, comp.c$Elev_m, method = "pearson")
# quite strong negative correlation
# p val = 0.0001808, cor = -0.802886
m2.pi.elev <- lm(Mean.pi.c ~ Elev_m, data = comp.c)
summary(m2.pi.elev)
### r2 = 0.6446. p val is significant (0.000181) and exactly the same as the correlation
### means that elevation does explain pi. (same as with the centromere included)
### estimate is -1.672e-6 so I think that is the beta?


## mean pi and SSN
plot(Mean.pi.c ~ Seq_PopFlwrMean, data = comp.c)
cor.test(comp.c$Mean.pi.c, comp.c$Seq_PopFlwrMean, method = "pearson")
# p = 0.3872
# cor = -0.2320122

m2.pi.ssn <- lm(Seq_PopFlwrMean ~ Mean.pi.c, dat = comp.c)
summary(m2.pi.ssn)
# r squared is tiny. 0.05383. p value (0.387) again exactly matches the correlation
# estimate is -66.6471

### Simulate Data for Plots

newdata3 <- data.frame(Elev_m = seq(55,1750, length = 2000))
newdata4 <- data.frame(Mean.pi.c = seq(0.0001,.004, length = 2000))

### Now use predict to get the predictions
m2.pi.elev.pred = predict(m2.pi.elev, newdata = newdata3, se.fit = TRUE)
m2.pi.ssn.pred = predict(m2.pi.ssn, newdata = newdata4, se.fit = TRUE)
forplot3 <- data.frame('newd3' = newdata3,
                       'newd4' = newdata4,
                       'pred.pi.el' = m2.pi.elev.pred,
                       'pred.pi.ssn' = m2.pi.ssn.pred)

###### fancy plots. ######
## pi
comp.c <- comp.c[order(comp.c$Elev_m),]
comp.c$pop <- as.factor(comp.c$pop)
comp.c$pop <- factor(comp.c$pop, levels = comp.c$pop[order(comp.c$Elev_m)])
str(comp.c$pop)
comp.c$labels <- paste0(comp.c$pop, " - ", comp.c$Elev_m, "m")
# with a linear fit line of elevation explaining mean pi
# topo
ggplot(comp.c)+
  geom_point(data=comp.c, aes(x=Elev_m, y=Mean.pi.c, shape = as.factor(pop), fill= Elev_m), col = "black", size = 5, stroke = 1.75)+
  geom_line(dat = forplot3, aes(x = Elev_m, y = pred.pi.el.fit), linetype = "solid",
            alpha = 0.7, linewidth = 1.25)+
  geom_line(dat = forplot3, aes(x = Elev_m, y = pred.pi.el.fit+1.96*pred.pi.el.se.fit), linetype = "solid",
            alpha = 0.3, linewidth = 0.75)+
  geom_line(dat = forplot3, aes(x = Elev_m, y = pred.pi.el.fit-1.96*pred.pi.el.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.75)+
  labs(x= "Elevation (m)", y= "Nucleotide Diversity", title = "")+
  scale_fill_gradientn(name = "Elevation", colours = topo.colors(16))+
  scale_shape_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20),
    axis.title = element_text(color = "black", size = 20),
    axis.text = element_text(color = "black", size = 20),
    legend.spacing.y = unit(0.03, "cm"))
ggsave(filename ="~/Figures/supp/ManuscriptFigs/PiPerPop_NoCent_topo.png",
       height = 7, width = 9, device = "png", dpi = 500)

##ssn by mean.pi
ggplot(comp.c)+
  geom_point(data=comp.c, aes(y=Seq_PopFlwrMean, x=Mean.pi.c, shape = as.factor(pop), fill = Elev_m), col = "black", size = 5, stroke = 1.75)+
  geom_line(dat = forplot3, aes(x = Mean.pi.c, y = pred.pi.ssn.fit), linetype = "solid",
            alpha = 0.7, size = 1.25)+
  geom_line(dat = forplot3, aes(x = Mean.pi.c, y = pred.pi.ssn.fit+1.96*pred.pi.ssn.se.fit), linetype = "solid",
            alpha = 0.3, size = 0.75)+
  geom_line(dat = forplot3, aes(x = Mean.pi.c, y = pred.pi.ssn.fit-1.96*pred.pi.ssn.se.fit), linetype = "solid",
            alpha = 0.5, size = 0.75)+
  labs(y="Mean Short Stamen Number" , x= "Nucleotide Diversity", title = "")+ # - Sequenced Only
  scale_fill_gradientn(name = "Elevation", colours = topo.colors(10))+
  scale_shape_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20),
    axis.title = element_text(color = "black", size = 20),
    axis.text = element_text(color = "black", size = 20),
    legend.spacing.y = unit(0.03, "cm"))
ggsave(filename ="~/Figures/ManuscriptFigs/supp/SSNbyPi_NoCent_topo.png",
       height = 7, width = 9, device = "png", dpi = 500)

########## Multiple Regression ##########
### with cent ###
# this was done 10/26 with the cent included files and a centered quadratic model to reduce collinearity.
m.ssn.elev.pi.good <- lm(Seq_PopFlwrMean ~ Elev_m + Elev_c2 + Mean.pi, dat = comp)
summary(m.ssn.elev.pi.good)
# results in model results excel sheet (supplemental)
# linear elevation is the only significant one.

# adding squared pi
m.ssn.elev.pi2 <- lm(Seq_PopFlwrMean ~ Elev_m + Elev_c2 + Mean.pi + pi_c2, dat = comp)
summary(m.ssn.elev.pi2)

# quick check that standardizing doesn't change p values:
tmp.standard2 <- lm(Seq_PopFlwrMean ~ scale(Elev_m) + scale(Elev_c2) + scale(Mean.pi), dat = comp)
summary(tmp.standard2)
# only p value that changes is the intercept

# need to do a sequenced lines by elevation model to get the residuals from it. 
tmp4_seq <- lm(Seq_PopFlwrMean ~ Elev_m + Elev_c2, data = comp)
plot(residuals(tmp4_seq))
# decently random when not coloring by anything
summary(tmp4_seq)

comp$Elev_residuals <- residuals(tmp4_seq) # index order matches row name order.
plot(comp$Elev_residuals, comp$Elev_m)

# now need to use the residuals in a plot - quick plot to start
plot(comp$Mean.pi, comp$Elev_residuals)
# residuals range from -0.6 to 0.3
m.resid.pi2 <- lm(Elev_residuals ~ Mean.pi, dat = comp)
summary(m.resid.pi2)
# mean.pi estimate = 46.8278, p value = 0.397, r squared = 0.05168 

# then I want the inverse. so the residuals of ssn ~ mean.pi regressed with elevation
tmp.resid2 <- residuals(lm(Seq_PopFlwrMean ~ Mean.pi, data = comp))
tmp.resid2
plot(tmp.resid2)

comp$Pi_residuals <- tmp.resid2
plot(comp$Elev_m, comp$Pi_residuals)

m.resid.elev2 <- lm(Pi_residuals ~ Elev_m + Elev_c2, dat = comp)
summary(m.resid.elev2)
# Elev_m estimate = 3.169e-04, p value = 0.0613,
# Elev_m2 estimate = -4.543e-07, p value = 0.1955
# r squared = 0.2559

## 3/6/2025: reviewer comment: is the pattern driven by outliers in fig 2D? (which is made below, but addressing with a model here)
# the two outlier populations are BOS and SAL. remove those.
comp_sub <- comp[!(comp$pop %in% c("BOS", "SAL")), ]
# run model and look at output
summary(lm(Seq_PopFlwrMean ~ Elev_m + Elev_c2 + Mean.pi, dat = comp_sub))
# general takeaway does not change after removing outliers. Elev and Elev^2 significant, pi less significant than before.


# use newdata (elevations) and newdata2 (pi with cent) from earlier

### Now use predict to get the predictions
newdata5 <- cbind(newdata, newdata2, Elev_residuals = seq(-0.7,0.4, length = 2000), Pi_residuals = seq(-0.8, 0.5, length = 2000) )
m.resid.pi2.pred = predict(m.resid.pi2, newdata = newdata5, se.fit = TRUE)
m.resid.elev2.pred = predict(m.resid.elev2, newdata = newdata5, se.fit = TRUE)
forplot5 <- data.frame('Elev_m' = newdata$Elev_m,
                       'Elev_m2' = newdata$Elev_m2,
                       'newd2' = newdata2,
                       'newd_res_elev' = newdata5$Elev_residuals,
                       'pred.m.resid.pi' = m.resid.pi2.pred,
                       'newd_res_pi' = newdata5$Pi_residuals,
                       'pred.m.resid.elev' = m.resid.elev2.pred)

## and plot
# might not need this code if running everything in a line, but I did it a different day so needed to order again
comp <- comp[order(comp$Elev_m),]
comp$pop <- as.factor(comp$pop)
comp$pop <- factor(comp$pop, levels = comp$pop[order(comp$Elev_m)])
str(comp$pop)
comp$labels <- paste0(comp$pop, " - ", comp$Elev_m, "m")


fig2d <- ggplot(comp)+
  geom_point(data=comp, aes(x=Mean.pi, y=Elev_residuals, shape = as.factor(pop), fill= Elev_m), col = "black", size = 3, stroke = 1, show.legend = FALSE)+
  geom_line(dat = forplot5, aes(x = Mean.pi, y = pred.m.resid.pi.fit), linetype = "solid",
            alpha = 0.9, linewidth = 1)+
  geom_line(dat = forplot5, aes(x = Mean.pi, y = pred.m.resid.pi.fit+1.96*pred.m.resid.pi.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.65)+
  geom_line(dat = forplot5, aes(x = Mean.pi, y = pred.m.resid.pi.fit-1.96*pred.m.resid.pi.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.65)+
  labs(x= "Nucleotide Diversity", y= "Residual Short Stamen Number")+
  annotate(geom = "text", x = -Inf, y = -Inf, label = "~beta==138.5 *';' ~italic(p)~ ' = 0.167'", parse = TRUE,
           hjust = -0.025, vjust = -0.5, size = 3, fontface = "bold")+
  scale_fill_gradientn(name = "Elevation", colours = topo.colors(16))+
  scale_shape_manual(name = "Population",
                     labels = comp$label,
                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 10),
    axis.text = element_text(color = "black", size = 10),
    legend.spacing.y = unit(0.03, "cm"))
fig2d <- annotate_figure(fig2d,
                         fig.lab = "D", fig.lab.face = "bold")
ggsave(filename ="~/Figures/ManuscriptFigs/SSNElevQuadResid_topo.png",
       height = 3, width = 3, device = "png", dpi = 700)

# pi residuals
fig2c <- ggplot(comp)+
  geom_point(data=comp, aes(x=Elev_m, y=Pi_residuals, shape = as.factor(pop), fill= Elev_m), col = "black", size = 3, stroke = 1, show.legend = FALSE)+
  geom_line(dat = forplot5, aes(x = Elev_m, y = pred.m.resid.elev.fit), linetype = "solid",
            alpha = 0.9, linewidth = 1.)+
  geom_line(dat = forplot5, aes(x = Elev_m, y = pred.m.resid.elev.fit+1.96*pred.m.resid.elev.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.65)+
  geom_line(dat = forplot5, aes(x = Elev_m, y = pred.m.resid.elev.fit-1.96*pred.m.resid.elev.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.65)+
  labs(y= "Residual Short Stamen Number", x= "Elevation (m)")+
  annotate(geom = "text", x = -Inf, y = -Inf, label = "~beta==6.47e10^-4 *';' ~italic(p)~ bold(' = 0.009')", parse = TRUE,
           hjust = -0.6, vjust = -1.7, size = 3, fontface = "bold")+
  annotate(geom = "text", x = -Inf, y = -Inf, label = "gamma==-2.99e10^-7 *';' ~italic(p)~ ' = 0.345'", parse = TRUE,
           hjust = -0.6, vjust = -.5, size = 3)+
  scale_fill_gradientn(name = "Elevation", colours = topo.colors(16))+
  scale_shape_manual(name = "Population",
                     labels = comp$label,
                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 10),
    axis.text = element_text(color = "black", size = 10),
    legend.spacing.y = unit(0.03, "cm"))
fig2c <- annotate_figure(fig2c,
                         fig.lab = "C", fig.lab.face = "bold")
ggsave(filename ="~/Figures/ManuscriptFigs/SSNPiResid_topo.png",
       height = 3, width = 3, device = "png", dpi = 700)

### no cent ###
# not updated with changed y axis labels.

# this is done with the no cent file pi values.
comp.c$Elev_m2 <- (comp.c$Elev_m)^2
m.ssn.elev.pi <- lm(Seq_PopFlwrMean ~ Elev_m + Elev_m2 + Mean.pi.c, dat = comp.c)
summary(m.ssn.elev.pi)
# results from Seq_PopFlwrMean ~ Elev_m + Mean.pi.c
# elev estimate = 6.320e-4, p value is 0.00904
# mean.pi.c estimate is 1.771e02, p value is 0.09714
#r squared is 0.4507

# results from Seq_PopFlwrMean ~ Elev_m + (Elev_m^2) + Mean.pi.c
# elev_m : estimate = 1.66e-03, p 0.0583
# elev_m2: estimate = -3.095e-07, p = 0.3227
# mean.pi.c : estimate = 1.551e02, p = 0.1508
# rsquared = 0.4954
# r2 got a little better, so this model should be doing better than the one without elev_m2?

# but I think a potential issue here is that I haven't standardized things? so there are big differences in scale between elevation and pi
tmp.standard <- lm(Seq_PopFlwrMean ~ scale(Elev_m) + scale(Elev_m^2) + scale(Mean.pi.c), dat = comp.c)
summary(tmp.standard)
# but the p values match the non-scaled values and the r squared is also exactly the same, so I don't think this is the issue
# these are the estimates:
# elev_m = 0.6315
# elev_m2 = -0.3170
# mean.pi.c = 0.1749

m.ssn.elev.c <- lm(Seq_PopFlwrMean ~ Elev_m + Elev_m2, data = comp.c)
summary(m.ssn.elev.c)
identical(residuals(m.ssn.elev.c), residuals(m.ssn.elev))

plot(residuals(m.ssn.elev.c))

# now need to use the residuals in a plot
comp.c$Elev_residuals <- residuals(m.ssn.elev.c) # index order matches row name order.

plot(comp.c$Mean.pi.c, comp.c$Elev_residuals)
# residuals range from -0.6 to 0.4
m.resid.pi <- lm(Elev_residuals ~ Mean.pi.c, dat = comp.c)
summary(m.resid.pi)
# mean.pi.c estimate = 52.6855, p value = 0.379, r squared = 0.05573 

# then I want the inverse. so the residuals of ssn ~ mean.pi.c regressed with elevation
# making new model b/c need order to match comp.c row order
tmp.resid <- residuals(lm(Seq_PopFlwrMean ~ Mean.pi.c, data = comp.c))
tmp.resid
plot(tmp.resid)

comp.c$Pi_residuals <- tmp.resid
plot(comp.c$Elev_m, comp.c$Pi_residuals)

m.resid.elev <- lm(Pi_residuals ~ Elev_m + Elev_m2, dat = comp.c)
summary(m.resid.elev)
# Elev_m estimate = 1.051e-03, p value = 0.1196,
# Elev_m2 estimate = -4.481e-07, p value = 0.2022
# r squared = 0.2528 so this is not a great fit but is much better than the other residual plot!


# use newdata (elevations) and 4 (mean.pi values) from earlier

### Now use predict to get the predictions
newdata5 <- cbind(newdata, newdata4, Elev_residuals = seq(-0.7,.4, length = 2000), Pi_residuals = seq(-0.8, 0.5, length = 2000) )
#m.ssn.elev.pi.pred = predict(m.ssn.elev.pi, newdata = newdata5, se.fit = TRUE)
m.resid.pi.pred = predict(m.resid.pi, newdata = newdata5, se.fit = TRUE)
m.resid.elev.pred = predict(m.resid.elev, newdata = newdata5, se.fit = TRUE)
forplot4 <- data.frame('Elev_m' = newdata$Elev_m,
                       'Elev_m2' = newdata$Elev_m2,
                       'newd4' = newdata4,
                       'newd_res_elev' = newdata5$Elev_residuals,
                       'pred.m.resid.pi' = m.resid.pi.pred,
                       'newd_res_pi' = newdata5$Pi_residuals,
                       'pred.m.resid.elev' = m.resid.elev.pred)

## and plot
ggplot(comp.c)+
  geom_point(data=comp.c, aes(x=Mean.pi.c, y=Elev_residuals, shape = as.factor(pop), fill= Elev_m), col = "black", size = 5, stroke = 1.75)+
  geom_line(dat = forplot4, aes(x = Mean.pi.c, y = pred.m.resid.pi.fit), linetype = "solid",
            alpha = 0.7, linewidth = 1.25)+
  geom_line(dat = forplot4, aes(x = Mean.pi.c, y = pred.m.resid.pi.fit+1.96*pred.m.resid.pi.se.fit), linetype = "solid",
            alpha = 0.3, linewidth = 0.75)+
  geom_line(dat = forplot4, aes(x = Mean.pi.c, y = pred.m.resid.pi.fit-1.96*pred.m.resid.pi.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.75)+
  labs(x= "Nucleotide Diversity", y= "Residuals of SSN ~ Elevation + Elevation^2", title = "")+
  scale_fill_gradientn(name = "Elevation", colours = topo.colors(16))+
  scale_shape_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20),
    axis.title = element_text(color = "black", size = 20),
    axis.text = element_text(color = "black", size = 19.5),
    legend.spacing.y = unit(0.03, "cm"))
ggsave(filename ="~/Figures/ManuscriptFigs/SSNElevQuadResid_NoCent_topo.png",
       height = 7, width = 9, device = "png", dpi = 500)

# pi residuals
ggplot(comp.c)+
  geom_point(data=comp.c, aes(x=Elev_m, y=Pi_residuals, shape = as.factor(pop), fill= Elev_m), col = "black", size = 5, stroke = 1.75)+
  geom_line(dat = forplot4, aes(x = Elev_m, y = pred.m.resid.elev.fit), linetype = "solid",
            alpha = 0.7, linewidth = 1.25)+
  geom_line(dat = forplot4, aes(x = Elev_m, y = pred.m.resid.elev.fit+1.96*pred.m.resid.elev.se.fit), linetype = "solid",
            alpha = 0.3, linewidth = 0.75)+
  geom_line(dat = forplot4, aes(x = Elev_m, y = pred.m.resid.elev.fit-1.96*pred.m.resid.elev.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.75)+
  labs(y= "Residuals of SSN ~ Nucleotide Diversity", x= "Elevation (m)")+ #title = "line from resid ~ elev + elev2"
  scale_fill_gradientn(name = "Elevation", colours = topo.colors(16))+
  scale_shape_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20),
    axis.title = element_text(color = "black", size = 20),
    axis.text = element_text(color = "black", size = 20),
    legend.spacing.y = unit(0.03, "cm"))
ggsave(filename ="~/Figures/ManuscriptFigs/SSNPiResid_NoCent_topo.png",
       height = 7, width = 9, device = "png", dpi = 500)

########## Figure Legend ##########
ggplot(comp.c)+
  geom_point(data=comp.c, aes(x=rep(1, times = 16), y=rev(c(1:16)), shape = as.factor(pop), fill= Elev_m), col = "black", size = 6, stroke = 1.75, show.legend = FALSE)+
  geom_text(aes(x=rep(1, times = 16), y=rev(c(1:16)), label = comp.c$label), size = 20/.pt, hjust = -0.25)+
  labs(x= "just one", y= "rank value", title = "")+
  scale_fill_gradientn(name = "Elevation", colours = topo.colors(16))+
  scale_shape_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
  theme_classic()+
  theme(
    axis.text.x=element_blank(), #remove x axis labels
    axis.ticks.x=element_blank(), #remove x axis ticks
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank(),  #remove y axis ticks
    axis.title = element_blank()
    )
ggsave(filename ="~/Figures/ManuscriptFigs/Legend.png",
       height = 7, width = 9, device = "png", dpi = 500)

########## AIC scores ##########
AIC(m.pi.elev)
AIC(m.pi.elev2)
AIC(m.ssn.elev.linear)
AIC(tmp4)
AIC(m.pi.ssn)
AIC(m.pi.ssn2)
AIC(m.ssn.elev.pi.good)
AIC(m.ssn.elev.pi2)
AIC(tmp4_seq)
AIC(m.resid.pi2)
AIC(m.resid.elev2)
