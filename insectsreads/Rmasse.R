############## ITS STATS TIME ##############################

setwd("E:/Cours/SupAgro/Stage M1/Analyses/Maya/")

#packages
library(data.table)
library(dplyr)
library(ggplot2)

rm(list=ls()) #to clear the memory if needed

##### INPUT DATA #####

bug<-read.table("newbug.txt",header=T)
var<-read.table("variables_mkb.txt", header=T)

##### merge rows/data treatment #####

aggbug<-aggregate(cbind(p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p1, p20, p21, p22, p23,
                        p24, p25, p26, p27, p28, p29, p2, p30, p31, p32, p33, p34, p35, p36, p37,
                        p38, p39, p3, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p4, p50, p51,
                        p52, p53, p54, p55, p56, p57, p58, p59, p5, p60, p61, p62, p63, p64, p65, p66,
                        p67, p68, p69, p6, p70, p71, p72, p73, p74, p75, p76, p77, p78, p79, p7, p80,
                        p81, p82, p83, p84, p85, p86, p87,p88, p89, p8, p90, p91, p92, p9) ~ spp + Guild, data = bug, sum)
#to merge the same insects together
aggbug<-aggbug[,-2] #to remove the guilds

tbug<-transpose(aggbug)
rownames(tbug)<-colnames(aggbug) #transpose titles too
colnames(tbug)<-tbug[1,] #first line is species names, we want them as titles
tbug<-tbug[-1,] #delete first line
tbug2<-data.frame(sample=row.names(tbug), tbug) #we want samples# into the table

readstable<-left_join(var,tbug2, by = "sample") # merge the tables with bugs and variables 

#We want to work on the mass of insects

entirebugs<-readstable[-1:-80,-1:-4] #first we keep samples of entire bugs only
entirebugs$mass = c(3,2,1,2,3,2,1,1,2,2,2,1) #we add a column of approximative/comparative masses (1 smallest species, 3 biggest)
#we want the values need to be read as numbers, not text
insectscol<- c(2:42)
entirebugs[,insectscol]<-lapply(entirebugs[,insectscol], as.numeric) 
#we create a new variable which is the sum of reads for each insect
sum_colbugs <- rowSums(entirebugs[, 2:41]) 
#we add it to the table
sumbugs <- cbind(entirebugs, sum_col = sum_colbugs)
#we delete all species columns 
smallbugs <- sumbugs[,-2:-41]

##### ANOVA on reads/mass

#check hypotheses first

res_aov <- aov(smallbugs$sum_col ~ smallbugs$mass,
               data = smallbugs)

#test normality residuals : 
hist(res_aov$residuals) # histogram #normality assumption not met

#non parametric test then
kruskal.test(smallbugs$sum_col ~ smallbugs$mass)
#NON
