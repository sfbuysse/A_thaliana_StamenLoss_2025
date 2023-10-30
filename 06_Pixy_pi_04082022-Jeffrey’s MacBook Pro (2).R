########## Visualize the Pixy results  ##########
library(dplyr)
library(ggplot2)
library(ggpubr)
####Need to read in the results
dat <- read.delim("C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/pixy_042022/all_filtered.50k_pi_042022.out_pi.txt", header=TRUE)
head(dat)
tail(dat)

########## with Centromeres ##########
#most of this just goes in the supplement.
###### Calculate Population global pi ######
##### use dyplyr to group by population
#### make a facto
###dat$pop <- as.factor(dat$pop)
#### calculate the aggregated value
#### (window 1 count_diffs + window 2 count_diffs) / (window 1 comparisons + window 2 comparisons)
sums <- dat %>% group_by(pop) %>% summarize(Sum.count.diff = sum(count_diffs, na.rm = T), Sum.count.comp = sum(count_comparisons, na.rm = T))
sums$Mean.pi <- sums$Sum.count.diff/sums$Sum.count.comp
####### this should be right following the pixy documentation
###### add the other population information
###### Should have labels ARU and SAL
###### this has both the full and only the sequenced means.
###load(file = "C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/StamenLossPipeline/04_Elev_Means_Feb2022.ROBJ")
###### merge the info together
###comp <- merge(sums, Elev_Means, by.x = "pop", by.y = "Population")
#### values are generally lower than what was calculated in february but not always consistently lower.

###### save as an R object
###save(comp, file = "C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/PopGlobalPi_allsite_04082022.ROBJ")

load("C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/PopGlobalPi_allsite_04082022.ROBJ")


###### Correlation/regressions between Variables ######
### Pi and Elevation ###
plot(Mean.pi ~ Elev_m, data = comp)
cor.test(comp$Mean.pi, comp$Elev_m, method = "pearson")
### strong correlation. r = -0.801 p = 0.0001942
m.pi.elev <- lm(Mean.pi ~ Elev_m, data = comp)
summary(m.pi.elev)
### r2 = 0.6411. p val is significant (but from stats I recall I should not look at this)
### means that elevation does explain pi. 

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


## model noodling to try type 3 to see if it matches Jeff's
#options(contrasts=c("contr.sum","contr.poly"))
#tmp <- lm(Full_PopFlwrMean ~ Elev_m + Elev_m2, data = comp)
#summary(tmp)
#
#tmp2 <- lm(Full_PopFlwrMean ~ Elev_m2 + Elev_m, data = comp)
#summary(tmp2)
#
#m.ssn.elev_3 <- lm(Full_PopFlwrMean ~ Elev_m + Elev_m2, data = comp, contrasts=list(Elev_m=contr.sum, Elev_m2=contr.sum))
#anova(m.ssn.elev) # this would be the default so it would be type I
#library(car)
#Anova(m.ssn.elev, type = 3)
comp$Elev_c2 <- (comp$Elev_m - mean(comp$Elev_m))^2
tmp4 <- lm(Full_PopFlwrMean ~ Elev_m + Elev_c2, data = comp)
summary(tmp4)
# yes, matches on 8/24/2023. need to determine if I want to center (and understand why) to decide which model to use

# see what output poly() gives me
# hmm. well this looks totally different in terms of estimates and the linear term p values
# the quadratic and full model p values have not changed; neither have the r squareds
# if I set raw = TRUE, then everything matches my original model (squared term not centered)
tmp5 <- lm(Full_PopFlwrMean ~ poly(Elev_m, 2, raw = TRUE), data = comp)
summary(tmp5)

# if I think about the decreasing colinearity argument, how correlated are my two terms? I would expect highly correlated...
cor(comp$Elev_m, comp$Elev_m2, method = "pearson")
# 0.9752654
cor(comp$Elev_m, comp$Elev_c2, method = "pearson")
# 0.4363318
# well it does decrease that correlation.

# try just a linear fit to see if it matches Jeff's model
m.ssn.elev.linear <- lm(Full_PopFlwrMean ~ Elev_m, data = comp)
summary(m.ssn.elev.linear)
# yep, this matches on 8/24/2023

### Pi and SSN ###
## here only use the sequenced lines because that is what pi was calculated from
plot(Seq_PopFlwrMean ~ Mean.pi, data = comp)
cor.test(comp$Mean.pi, comp$
           Seq_PopFlwrMean, method = "pearson")
### r = -0.2323357 p = 0.3865
m.pi.ssn <- lm(Seq_PopFlwrMean ~ Mean.pi, dat = comp)
summary(m.pi.ssn)
### r2 = 0.054 so not a good fit. pi does not explain short stamen number. 
### this is an interesting result. does not provide evidence for the drift theory because 
### drift would show that at lower pi values (less nucleotide diverity) there would be less short stamen loss.
### the slope is negative, which is the direction we expect
### this graph is later replaced with a multiple regression later on 

### So, both mean.pi and mean.ssn are decently correlated with elevation, but not with each other! 
### interesting.probably makes sense though because pi is genome wide

### The two means ###
cor.test(comp$Full_PopFlwrMean, comp$Seq_PopFlwrMean, method = "pearson")
# r = 0.985 p < 0.001

###### simulate data to draw lines on the plots where elevation is the predictor ######
new_elev <- seq(55,1750, length = 2000)
newdata <- data.frame(Elev_m = new_elev, Elev_m2 = new_elev**2)
newdata2 <- data.frame(Mean.pi = seq(0.001,.0051, length = 2000))
### Now use predict to get the predictions
### can't use se.fit = T for mixed effects model? 
m.ssn.elev.pred = predict(m.ssn.elev, newdata = newdata, se.fit = TRUE)
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
### this looks like the right order
### key is that you need to order the whole sheet AND reorder the factor. BOTH, not just one.
comp$labels <- paste0(comp$pop, " - ", comp$Elev_m, "m")

# with terrain, I can get as many colors as I want I think bc it is a function
# for PC where I didn't want it filled in before because there is overlap, maybe I should just change alpha?
barplot(rep(1,6), col = terrain.colors(6))
barplot(rep(1,10), col = terrain.colors(10))
barplot(rep(1,16), col = terrain.colors(16))
# for topo color scheme
barplot(rep(1,6), col = topo.colors(6))
barplot(rep(1,10), col = topo.colors(10))
barplot(rep(1,16), col = topo.colors(16))
# choosing topo as the color scheme, so commenting out all the terrain plots.
# with a linear fit line of elevation explaining mean pi
ggplot(comp)+
  geom_point(data=comp, aes(x=Elev_m, y=Mean.pi, shape = as.factor(pop), fill= Elev_m), col = "black", size = 5, stroke = 1.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit), linetype = "solid",
            alpha = 0.7, linewidth = 1.25)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit+1.96*pred.pi.el.se.fit), linetype = "solid",
            alpha = 0.3, linewidth = 0.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit-1.96*pred.pi.el.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.75)+
  labs(x= "Elevation (m)", y= "Nucleotide Diversity", title = "")+
  scale_fill_gradientn(name = "Elevation", colours = topo.colors(16))+
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
ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/supp/PiPerPop_topo.png",
       height = 7, width = 9, device = "png", dpi = 500)
# can also give plot name (useful if making multipanel plots in r) and give dimensions in pixels

#ggplot(comp)+
#  geom_point(data=comp, aes(x=Elev_m, y=Mean.pi, shape = as.factor(pop), fill= Elev_m), col = "black", size = 5, stroke = 1.75)+
#  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit), linetype = "solid",
#            alpha = 0.7, linewidth = 1.25)+
#  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit+1.96*pred.pi.el.se.fit), linetype = "solid",
#            alpha = 0.3, linewidth = 0.75)+
#  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit-1.96*pred.pi.el.se.fit), linetype = "solid",
#            alpha = 0.5, linewidth = 0.75)+
#  labs(x= "Elevation (m)", y= "Mean Pi", title = "")+
#  scale_fill_gradientn(name = "Elevation", colours = terrain.colors(10))+
#  scale_shape_manual(name = "Population",
#                     labels = comp$label,
#                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
#  theme_classic()+
#  theme(
#    legend.title = element_text(color = "black", size = 20),
#    legend.text = element_text(color = "black", size = 20),
#    axis.title = element_text(color = "black", size = 20),
#    axis.text = element_text(color = "black", size = 20),
#    legend.spacing.y = unit(0.03, "cm"))
#ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/PiPerPop_terrain.png", 
#       height = 7, width = 9,device = "png", dpi = 500)

# with a linear fit line of elevation explaining mean pi
ggplot(comp)+
  geom_point(data=comp, aes(x=Elev_m, y=Mean.pi, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit), linetype = "solid",
            alpha = 0.7, linewidth = 1.25)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit+1.96*pred.pi.el.se.fit), linetype = "solid",
            alpha = 0.3, linewidth = 0.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.pi.el.fit-1.96*pred.pi.el.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.75)+
  labs(x= "Elevation (m)", y= "Nucleotide Diversity", title = "")+
  scale_color_manual(name = "Population",
                     labels = comp$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4))))+
  scale_shape_manual(name = "Population",
                     labels = comp$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 16),
    legend.spacing.y = unit(0.03, "cm"))

ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/PiPerPop_pixy_04082022.png", height = 7, width = 9)

### ssn by elevation
ggplot(comp)+
  geom_point(data=comp, aes(x=Elev_m, y=Full_PopFlwrMean, shape = as.factor(pop), fill = Elev_m), col = "black", size = 5, stroke = 1.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit), linetype = "solid",
            alpha = 0.7, linewidth = 1.25)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit+1.96*pred.ssn.el.se.fit), linetype = "solid",
            alpha = 0.3, linewidth = 0.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit-1.96*pred.ssn.el.se.fit), linetype = "solid",
            alpha = 0.5, linewidth = 0.75)+
  labs(x= "Elevation (m)", y= "Mean Short Stamen Number", title = "")+ # - Full Phenotype Set
  scale_fill_gradientn(name = "Elevation", colours = topo.colors(10))+
  scale_shape_manual(name = "Population",
                     labels = comp$label,
                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20),
    legend.spacing.y = unit(0.03, "cm"),
    axis.title = element_text(color = "black", size = 20),
    axis.text = element_text(color = "black", size = 20))
ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/MeanSSNByElev_topo.png",
       height = 7, width = 9, device = "png", dpi = 500)

#ggplot(comp)+
#  geom_point(data=comp, aes(x=Elev_m, y=Full_PopFlwrMean, shape = as.factor(pop), fill = Elev_m), col = "black", size = 5, stroke = 1.75)+
#  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit), linetype = "solid",
#            alpha = 0.7, size = 1.25)+
#  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit+1.96*pred.ssn.el.se.fit), linetype = "solid",
#            alpha = 0.3, size = 0.75)+
#  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit-1.96*pred.ssn.el.se.fit), linetype = "solid",
#            alpha = 0.5, size = 0.75)+
#  labs(x= "Elevation (m)", y= "Mean Short Stamen Number - Full Phenotype Set", title = "")+
#  scale_fill_gradientn(name = "Elevation", colours = terrain.colors(10))+
#  scale_shape_manual(name = "Population",
#                     labels = comp$label,
#                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
#  theme_classic()+
#  theme(
#    legend.title = element_text(color = "black", size = 20),
#    legend.text = element_text(color = "black", size = 20),
#    legend.spacing.y = unit(0.03, "cm"),
#    axis.title = element_text(color = "black", size = 20),
#    axis.text = element_text(color = "black", size = 20))
#ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/MeanSSNByElev_terrain.png",
#       height = 7, width = 9, device = "png", dpi = 500)

# old plot
ggplot(comp)+
  geom_point(data=comp, aes(x=Elev_m, y=Full_PopFlwrMean, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit), linetype = "solid",
            alpha = 0.7, size = 1.25)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit+1.96*pred.ssn.el.se.fit), linetype = "solid",
            alpha = 0.3, size = 0.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit-1.96*pred.ssn.el.se.fit), linetype = "solid",
            alpha = 0.5, size = 0.75)+
  labs(x= "Elevation (m)", y= "Mean Short Stamen Number - Full Phenotype Set", title = "")+
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
    legend.spacing.y = unit(0.03, "cm"),
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 16))
ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/MeanSSNByElev_02082022_quad.png", height = 7, width = 9)
## haha!

### files for Jeff Dec 14, 2022
ggplot(comp)+
  geom_point(data=comp, aes(x=Elev_m, y=Full_PopFlwrMean, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75, show.legend = FALSE)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit), linetype = "solid",
            alpha = 0.7, size = 1.25)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit+1.96*pred.ssn.el.se.fit), linetype = "solid",
            alpha = 0.3, size = 0.75)+
  geom_line(dat = forplot2, aes(x = newd1, y = pred.ssn.el.fit-1.96*pred.ssn.el.se.fit), linetype = "solid",
            alpha = 0.5, size = 0.75)+
  labs(x= "Elevation (m)", y= "Mean Short Stamen Number", title = "")+
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
    legend.spacing.y = unit(0.03, "cm"),
    axis.title = element_text(color = "black", size = 20),
    axis.text = element_text(color = "black", size = 16))
ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/MeanSSNByElev_12142022_quad.jpg", height = 7, width = 7, dpi = 300)

ggplot(comp)+
  geom_point(data=comp, aes(x=Elev_m, y=Full_PopFlwrMean, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75, show.legend = FALSE)+
  labs(x= "Elevation (m)", y= "Mean Short Stamen Number", title = "")+
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
    legend.spacing.y = unit(0.03, "cm"),
    axis.title = element_text(color = "black", size = 24),
    axis.text = element_text(color = "black", size = 24))
ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/MeanSSNByElev_12142022.jpg", height = 7, width = 7, dpi = 300)


# lowess line instead
ggplot(comp)+
  geom_point(data=comp, aes(x=Elev_m, y=Full_PopFlwrMean, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_smooth(aes(x=Elev_m, y=Full_PopFlwrMean),stat = "smooth",position = "identity",method = "loess",
    formula = y ~ x + (x**2),
    se = TRUE, span = 0.99, na.rm = FALSE)+
  labs(x= "Elevation (m)", y= "Mean Short Stamen Number - Full", title = "")+
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
    legend.spacing.y = unit(0.03, "cm"),
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 16))


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
ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/supp/SSNbyPi_topo.png",
       height = 7, width = 9, device = "png", dpi = 500)

#ggplot(comp)+
#  geom_point(data=comp, aes(y=Seq_PopFlwrMean, x=Mean.pi, shape = as.factor(pop), fill = Elev_m), col = "black", size = 5, stroke = 1.75)+
#  geom_line(dat = forplot2, aes(x = Mean.pi, y = pred.pi.ssn.fit), linetype = "solid",
#            alpha = 0.7, size = 1.25)+
#  geom_line(dat = forplot2, aes(x = Mean.pi, y = pred.pi.ssn.fit+1.96*pred.pi.ssn.se.fit), linetype = "solid",
#            alpha = 0.3, size = 0.75)+
#  geom_line(dat = forplot2, aes(x = Mean.pi, y = pred.pi.ssn.fit-1.96*pred.pi.ssn.se.fit), linetype = "solid",
#            alpha = 0.5, size = 0.75)+
#  labs(y="Mean Short Stamen Number - Sequenced Only" , x= "Mean Pi", title = "")+
#  scale_fill_gradientn(name = "Elevation", colours = terrain.colors(10))+
#  scale_shape_manual(name = "Population",
#                     labels = comp$label,
#                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
#  theme_classic()+
#  theme(
#    legend.title = element_text(color = "black", size = 20),
#    legend.text = element_text(color = "black", size = 20),
#    axis.title = element_text(color = "black", size = 20),
#    axis.text = element_text(color = "black", size = 20),
#    legend.spacing.y = unit(0.03, "cm"))
#ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/SSNbyPi_terrain.png",
#       height = 7, width = 9, device = "png", dpi = 500)

#old plot
ggplot(comp)+
  geom_point(data=comp, aes(y=Seq_PopFlwrMean, x=Mean.pi, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_line(dat = forplot2, aes(x = Mean.pi, y = pred.pi.ssn.fit), linetype = "solid",
            alpha = 0.7, size = 1.25)+
  geom_line(dat = forplot2, aes(x = Mean.pi, y = pred.pi.ssn.fit+1.96*pred.pi.ssn.se.fit), linetype = "solid",
            alpha = 0.3, size = 0.75)+
  geom_line(dat = forplot2, aes(x = Mean.pi, y = pred.pi.ssn.fit-1.96*pred.pi.ssn.se.fit), linetype = "solid",
            alpha = 0.5, size = 0.75)+
  labs(y="Mean Short Stamen Number - Sequenced Only" , x= "Mean Pi", title = "")+
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
ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/SSNbyPi_04082022.png", height = 7, width = 9)

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
ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/MeanComparison_02022022.png", height = 7, width = 9)

###### Population windowed pi ######
#the purpose of this is to find regions of elevated pi (kinda a selection test but I removed the centromere instead of doing any stats)
# this is informing which region of the centromere to exclude, though I didn't actually do the excluding until later.
# the pi values I calculated previously in this code DO NOT exclude the centromere region. Should they? Depends on what I show below.

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
ggsave("C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/50kb_pixy.pi_04082022.png", plot = all_plot, width = 20, height = 20)

# there are definitely still peaks in pi at the centromeres.

###### for Kieran ######
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
ggsave("C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/50kb_sites_04082022.png", plot = all_plot, width = 20, height = 20)
Coc_sites <- varName[[6]]

## still dips in sites per window at the centromere as well.

#### look for correlation between the number of sites and the pi value in each window
#### does a correlation make sense? If sites is number of variant sites, a positive correlation makes sense
### not totally sure what to think about this.
ggplot(dat, aes(x = no_sites, y = avg_pi))+
  geom_point(aes(col = pop))
# negative trend
ggplot(dat, aes(x = count_diffs, y = avg_pi))+
  geom_point(aes(col = pop))
# positive but not really a clear trend. more dnesity near a cutoff line at the bottom
ggplot(dat, aes(x = count_comparisons, y = avg_pi))+
  geom_point(aes(col = pop))
# lots at zero and then kind of a hill. huh. This one maybe looks weird but idk how to interpret it or what it 'should' look like
ggplot(dat, aes(x = count_missing, y = avg_pi))+
  geom_point(aes(col = pop))
# nothing really here. where there isthe most count_missing, it is an intermediate pi value.

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

# last ran on 4/12/2022
# thoughts: do I need to exclude the centromere region and calculate pi and run all the analyses again? 
# I think that makes sense in case the centromere region is shifting the global pi calculation.

########## Without Centromeres ##########
require(dplyr)
require(ggplot2)
require(ggpubr)
#Need to read in the results

#### pi calculated in December, 2021
###dat <- read.delim("C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/pixy_042022/all_filtered.50k_pi_042022.out_pi.txt", header=TRUE)
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
###save(dat.clean, file = "C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/pi_noCent_50k_04082022.ROBJ")
load("C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/pi_noCent_50k_04082022.ROBJ")

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

ggplot()+
  geom_point()+
  geom_abline(aes(intercept = 0, slope = 1))+
  theme_classic()

# this section is a little out of order b/c I need the labels from comp2.c dataframe
### something went wrong here and the coloring doesn't seem to  match up
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
#load(file = "C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/StamenLossPipeline/04_Elev_Means_Feb2022.ROBJ")

#comp.c <- merge(sums.c, Elev_Means, by.x = "pop", by.y = "Population")

# check how similar the pi values are
cor(x=comp$Mean.pi, y=comp.c$Mean.pi.c)
# 0.9981597

#save(comp.c, file = "C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/PopGlobalPi_NoCent_04082022.ROBJ")
load("C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/R_script/PopGlobalPi_NoCent_04082022.ROBJ")

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
### r2 = 0.6446. p val is significant (0.000181) and exactly the same as the correlation(but from stats I recall I should not look at this)
### means that elevation does explain pi. (same as with the centromere included)
### estimate is -1.672e-6 so I think that is the beta?


## mean pi and SSN
plot(Mean.pi.c ~ Seq_PopFlwrMean, data = comp.c)
cor.test(comp.c$Mean.pi.c, comp.c$Seq_PopFlwrMean, method = "pearson")
# p = 0.3872
# cor = -0.2320122
# pretty weak like nothing going on here
m2.pi.ssn <- lm(Seq_PopFlwrMean ~ Mean.pi.c, dat = comp.c)
summary(m2.pi.ssn)
# r squared is tiny. 0.05383. p value (0.387) again exactly matches the correlation
# so really, no relationship between these two variables shown here
# so this does not support drift b/c genetic variation doesn't explain short stamen number
# OK, so no evidence for drift which changes what my previous narrative was.
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
# this looks like the right order
# key is that you need to order the whole sheet AND reorder the factor. BOTH, not just one.
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
ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/PiPerPop_NoCent_topo.png",
       height = 7, width = 9, device = "png", dpi = 500)
# can also give plot name (useful if making multipanel plots in r) and give dimensions in pixels

#ggplot(comp.c)+
#  geom_point(data=comp.c, aes(x=Elev_m, y=Mean.pi.c, shape = as.factor(pop), fill= Elev_m), col = "black", size = 5, stroke = 1.75)+
#  geom_line(dat = forplot3, aes(x = Elev_m, y = pred.pi.el.fit), linetype = "solid",
#            alpha = 0.7, linewidth = 1.25)+
#  geom_line(dat = forplot3, aes(x = Elev_m, y = pred.pi.el.fit+1.96*pred.pi.el.se.fit), linetype = "solid",
#            alpha = 0.3, linewidth = 0.75)+
#  geom_line(dat = forplot3, aes(x = Elev_m, y = pred.pi.el.fit-1.96*pred.pi.el.se.fit), linetype = "solid",
#            alpha = 0.5, linewidth = 0.75)+
#  labs(x= "Elevation (m)", y= "Mean Pi", title = "")+
#  scale_fill_gradientn(name = "Elevation", colours = terrain.colors(10))+
#  scale_shape_manual(name = "Population",
#                     labels = comp.c$label,
#                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
#  theme_classic()+
#  theme(
#    legend.title = element_text(color = "black", size = 20),
#    legend.text = element_text(color = "black", size = 20),
#    axis.title = element_text(color = "black", size = 20),
#    axis.text = element_text(color = "black", size = 20),
#    legend.spacing.y = unit(0.03, "cm"))
#ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/PiPerPop_NoCent_terrain.png", 
#       height = 7, width = 9,device = "png", dpi = 500)
# old plot
ggplot(comp.c)+
  geom_point(data=comp.c, aes(x=Elev_m, y=Mean.pi.c, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_line(dat = forplot3, aes(x = Elev_m, y = pred.pi.el.fit), linetype = "solid",
            alpha = 0.7, size = 1.25)+
  geom_line(dat = forplot3, aes(x = Elev_m, y = pred.pi.el.fit+1.96*pred.pi.el.se.fit), linetype = "solid",
            alpha = 0.3, size = 0.75)+
  geom_line(dat = forplot3, aes(x = Elev_m, y = pred.pi.el.fit-1.96*pred.pi.el.se.fit), linetype = "solid",
            alpha = 0.5, size = 0.75)+
  labs(x= "Elevation (m)", y= "Mean Pi", title = "")+
  scale_color_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4))))+
  scale_shape_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 16),
    legend.spacing.y = unit(0.03, "cm"))
ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/PiPerPop_NoCent_50k_pixy_04082022.png", height = 7, width = 9)

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
ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/SSNbyPi_NoCent_topo.png",
       height = 7, width = 9, device = "png", dpi = 500)

#ggplot(comp.c)+
#  geom_point(data=comp.c, aes(y=Seq_PopFlwrMean, x=Mean.pi.c, shape = as.factor(pop), fill = Elev_m), col = "black", size = 5, stroke = 1.75)+
#  geom_line(dat = forplot3, aes(x = Mean.pi.c, y = pred.pi.ssn.fit), linetype = "solid",
#            alpha = 0.7, size = 1.25)+
#  geom_line(dat = forplot3, aes(x = Mean.pi.c, y = pred.pi.ssn.fit+1.96*pred.pi.ssn.se.fit), linetype = "solid",
#            alpha = 0.3, size = 0.75)+
#  geom_line(dat = forplot3, aes(x = Mean.pi.c, y = pred.pi.ssn.fit-1.96*pred.pi.ssn.se.fit), linetype = "solid",
#            alpha = 0.5, size = 0.75)+
#  labs(y="Mean Short Stamen Number - Sequenced Only" , x= "Mean Nucleotide Diversity", title = "")+
#  scale_fill_gradientn(name = "Elevation", colours = terrain.colors(10))+
#  scale_shape_manual(name = "Population",
#                     labels = comp.c$label,
#                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
#  theme_classic()+
#  theme(
#    legend.title = element_text(color = "black", size = 20),
#    legend.text = element_text(color = "black", size = 20),
#    axis.title = element_text(color = "black", size = 20),
#    axis.text = element_text(color = "black", size = 20),
#    legend.spacing.y = unit(0.03, "cm"))
#ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/SSNbyPi_NoCent_terrain.png",
#       height = 7, width = 9, device = "png", dpi = 500)


# old plot
ggplot(comp.c)+
  geom_point(data=comp.c, aes(x=Mean.pi.c, y=Seq_PopFlwrMean, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_line(dat = forplot3, aes(x = Mean.pi.c, y = pred.pi.ssn.fit), linetype = "solid",
            alpha = 0.7, size = 1.25)+
  geom_line(dat = forplot3, aes(x = Mean.pi.c, y = pred.pi.ssn.fit+1.96*pred.pi.ssn.se.fit), linetype = "solid",
            alpha = 0.3, size = 0.75)+
  geom_line(dat = forplot3, aes(x = Mean.pi.c, y = pred.pi.ssn.fit-1.96*pred.pi.ssn.se.fit), linetype = "solid",
            alpha = 0.5, size = 0.75)+
  labs(y= "Mean Short Stamen Number for Sequenced Lines", x= "Mean Pi", title = "")+
  scale_color_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4))))+
  scale_shape_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 16),
    legend.spacing.y = unit(0.03, "cm"))

ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/SSNBYPi_NoCent_50k_pixy_04082022.png", height = 7, width = 9)

# windowed to check that it looks like I am expecting it to. The values are lower here if I look at the y axis scale but the trends remain the same.
###### Population windowed pi ######
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
ggsave("C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/50kb_NoCent_pixy_04082022.png", plot = all_plot, width = 20, height = 20)
# I think this looks as expected and is what I should use.
# would want to use same spacing as the unremoved centromere if I am actually comparing the two though.


###### more complex model! ######
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
# okay. so the estimates here have changes (because scale changed, so that makes sense)
# but the p values match the non-scaled values and the r squared is also exactly the same, so I don't think this is the issue
# just for the sake of reporting, these are the estimates:
# elev_m = 0.6315
# elev_m2 = -0.3170
# mean.pi.c = 0.1749

# to graph, I want the residuals of ssn ~ elevation and plot those against mean pi I think. not sure how I would make a line exactly, but that is what I would want intially
# should this be the residuals from the quadratic regression? I vote yes. so also changed the model above to include the quadratic term
#m.tmp <- lm(Seq_PopFlwrMean ~ Elev_m, dat = comp.c)
#plot(residuals(m.tmp))
# ok, so these seem decently random I guess

m.ssn.elev.c <- lm(Seq_PopFlwrMean ~ Elev_m + Elev_m2, data = comp.c)
summary(m.ssn.elev.c)
identical(residuals(m.ssn.elev.c), residuals(m.ssn.elev))
# why do these two not match?
# because m.ssn.elev uses Full_PopFlwrMean but I want to only use the sequenced lines here.

# in this model (m.ssn.elev.c) elev_m has ap value of 0.085 and elev_m2 has a p value of 0.21 and the r2 is 0.3964
# for reference, in this same dataset, the rsquared of just the linear term is 0.31
#summary(lm(Seq_PopFlwrMean ~ Elev_m, data = comp.c))

plot(residuals(m.ssn.elev.c))
# kinda seems like two chunks? like two arcs?

# now need to use the residuals in a plot

#comp.c$Elev_residuals <- residuals(m.tmp)
comp.c$Elev_residuals <- residuals(m.ssn.elev.c) # index order matches row name order.

plot(comp.c$Mean.pi.c, comp.c$Elev_residuals)
# there is maybe a positive relationship here? worth a model I think.
# residuals range from -0.6 to 0.4
m.resid.pi <- lm(Elev_residuals ~ Mean.pi.c, dat = comp.c)
summary(m.resid.pi)
# mean.pi.c estimate = 52.6855, p value = 0.379, r squared = 0.05573 so this is a terrible fit.

# then I want the inverse. so the residuals of ssn ~ mean.pi.c regressed with elevation
# making new model b/c need order to match comp.c row order
tmp.resid <- residuals(lm(Seq_PopFlwrMean ~ Mean.pi.c, data = comp.c))
tmp.resid
plot(tmp.resid)
# kinda two curves again?
comp.c$Pi_residuals <- tmp.resid
plot(comp.c$Elev_m, comp.c$Pi_residuals)
# not great, but doesn't include the quadratic term.
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
ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/SSNElevQuadResid_NoCent_topo.png",
       height = 7, width = 9, device = "png", dpi = 500)
#ggplot(comp.c)+
#  geom_point(data=comp.c, aes(x=Mean.pi.c, y=Elev_residuals, shape = as.factor(pop), fill= Elev_m), col = "black", size = 5, stroke = 1.75)+
#  geom_line(dat = forplot4, aes(x = Mean.pi.c, y = pred.m.resid.pi.fit), linetype = "solid",
#            alpha = 0.7, linewidth = 1.25)+
#  geom_line(dat = forplot4, aes(x = Mean.pi.c, y = pred.m.resid.pi.fit+1.96*pred.m.resid.pi.se.fit), linetype = "solid",
#            alpha = 0.3, linewidth = 0.75)+
#  geom_line(dat = forplot4, aes(x = Mean.pi.c, y = pred.m.resid.pi.fit-1.96*pred.m.resid.pi.se.fit), linetype = "solid",
#            alpha = 0.5, linewidth = 0.75)+
#  labs(x= "Mean Pi", y= "Residuals of SSN ~ Elevation + Elevation^2", title = "")+
#  scale_fill_gradientn(name = "Elevation", colours = terrain.colors(16))+
#  scale_shape_manual(name = "Population",
#                     labels = comp.c$label,
#                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
#  theme_classic()+
#  theme(
#    legend.title = element_text(color = "black", size = 20),
#    legend.text = element_text(color = "black", size = 20),
#    axis.title = element_text(color = "black", size = 20),
#    axis.text = element_text(color = "black", size = 20),
#    legend.spacing.y = unit(0.03, "cm"))
#ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/SSNElevQuadResid_NoCent_terrain.png",
#       height = 7, width = 9, device = "png", dpi = 500)

ggplot(comp.c)+
  geom_point(data=comp.c, aes(x=Mean.pi.c, y=Elev_residuals, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_line(dat = forplot4, aes(x = Mean.pi.c, y = pred.m.resid.pi.fit), linetype = "solid",
            alpha = 0.7, size = 1.25)+
  geom_line(dat = forplot4, aes(x = Mean.pi.c, y = pred.m.resid.pi.fit+1.96*pred.m.resid.pi.se.fit), linetype = "solid",
            alpha = 0.3, size = 0.75)+
  geom_line(dat = forplot4, aes(x = Mean.pi.c, y = pred.m.resid.pi.fit-1.96*pred.m.resid.pi.se.fit), linetype = "solid",
            alpha = 0.5, size = 0.75)+
  labs(y= "Residuals of SSN ~ Elevation + Elevation^2", x= "Mean Pi")+
  scale_color_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4))))+
  scale_shape_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 16),
    legend.spacing.y = unit(0.03, "cm"))
#ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/SSNElevQuadResid_ByPi_06232022.png", height = 7, width = 9)

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
ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/SSNPiResid_NoCent_topo.png",
       height = 7, width = 9, device = "png", dpi = 500)

#ggplot(comp.c)+
#  geom_point(data=comp.c, aes(x=Elev_m, y=Pi_residuals, shape = as.factor(pop), fill= Elev_m), col = "black", size = 5, stroke = 1.75)+
#  geom_line(dat = forplot4, aes(x = Elev_m, y = pred.m.resid.elev.fit), linetype = "solid",
#            alpha = 0.7, linewidth = 1.25)+
#  geom_line(dat = forplot4, aes(x = Elev_m, y = pred.m.resid.elev.fit+1.96*pred.m.resid.elev.se.fit), linetype = "solid",
#            alpha = 0.3, linewidth = 0.75)+
#  geom_line(dat = forplot4, aes(x = Elev_m, y = pred.m.resid.elev.fit-1.96*pred.m.resid.elev.se.fit), linetype = "solid",
#            alpha = 0.5, linewidth = 0.75)+
#  labs(y= "Residuals of SSN ~ Mean.pi", x= "Elevation (m)", title = "line from resid ~ elev + elev2")+
#  scale_fill_gradientn(name = "Elevation", colours = terrain.colors(16))+
#  scale_shape_manual(name = "Population",
#                     labels = comp.c$label,
#                     values = c(rep(c(22, 21, 24, 23, 25), times = 4)))+
#  theme_classic()+
#  theme(
#    legend.title = element_text(color = "black", size = 20),
#    legend.text = element_text(color = "black", size = 20),
#    axis.title = element_text(color = "black", size = 20),
#    axis.text = element_text(color = "black", size = 20),
#    legend.spacing.y = unit(0.03, "cm"))
#ggsave(filename ="C:/Users/Sophia/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/SSNPiResid_NoCent_terrain.png",
#       height = 7, width = 9, device = "png", dpi = 500)

ggplot(comp.c)+
  geom_point(data=comp.c, aes(x=Elev_m, y=Pi_residuals, shape = as.factor(pop), col= as.factor(pop)), size = 5, stroke = 1.75)+
  geom_line(dat = forplot4, aes(x = Elev_m, y = pred.m.resid.elev.fit), linetype = "solid",
            alpha = 0.7, size = 1.25)+
  geom_line(dat = forplot4, aes(x = Elev_m, y = pred.m.resid.elev.fit+1.96*pred.m.resid.pi.se.fit), linetype = "solid",
            alpha = 0.3, size = 0.75)+
  geom_line(dat = forplot4, aes(x = Elev_m, y = pred.m.resid.elev.fit-1.96*pred.m.resid.pi.se.fit), linetype = "solid",
            alpha = 0.5, size = 0.75)+
  labs(y= "Residuals of SSN ~ Mean.pi", x= "Elevation (m)", title = "line from resid ~ elev + elev2")+
  scale_color_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c("red", "orange", "green", "blue"), times = c(4,4,4,4))))+
  scale_shape_manual(name = "Population",
                     labels = comp.c$label,
                     values = c(rep(c(0, 1, 2, 5), times = 4)))+
  theme_classic()+
  theme(
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 16),
    legend.spacing.y = unit(0.03, "cm"))

### test an interaction term
m2.ssn.elev.pi <- lm(Seq_PopFlwrMean ~ Elev_m * Mean.pi.c, dat = comp.c)
summary(m2.ssn.elev.pi)
# if I add the interaction, nothing is significant
# do I want an interaction?
AIC(m.ssn.elev.pi)
# 6.71003
# this one is a better fit than the interaction model
AIC(m2.ssn.elev.pi)
# 8.528473

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
ggsave(filename ="C:/Users/Sophie/Michigan State University/Conner, Jeffrey - SophieAnalyses/Figures/ManuscriptFigs/Legend.png",
       height = 7, width = 9, device = "png", dpi = 500)