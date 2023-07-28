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

###CREATING TABLE WITH BINARY VALUES (PRESENT/MISSING)
binarytable<-readstable #second table created based on reads one
w<-(readstable[c(6:45)]) #cut table to only work on insects var
w[w>0]<-1 #if 1 read or more, take value 1 (presence) (else, still 0, absence)
binarytable[c(6:45)]<-w #replace columns in table with binary values created

#we want the values need to be read as numbers, not text
insectscol<- c(6:45)
binarytable[,insectscol]<-lapply(binarytable[,insectscol], as.numeric) 
readstable[,insectscol]<-lapply(readstable[,insectscol], as.numeric) 

#we create a new variable which is the sum
sum_colbin <- rowSums(binarytable[, 6:45]) 
sum_colreads <- rowSums(readstable[, 6:45]) 

#we add it to the table
newbin <- cbind(binarytable, sum_col = sum_colbin)
newreads <- cbind(readstable, sum_col = sum_colreads)

#same with a new variable which is presence/absence of detection for each sample
z <- sum_colbin
z[z>0]<-1
newbin2<- cbind(newbin, colbin = z)
newreads2<- cbind(newreads, colbin = z)

#we delete all species columns and the samples of entire insects
smallbin <- newbin2[-(81:92),-(5:45)]
smallreads <- newreads2[-(81:92),-(5:45)]

#mean of sum_col depending on swipes
meanbin<-tapply(smallbin$sum_col,smallbin$swipes,mean)
meanreads<-tapply(smallreads$sum_col,smallreads$swipes,mean)

###############GRAPHS

hist(smallbin$sum_col) 
hist(smallreads$sum_col) 
#super simple hist

###nb swipes

#for species only
dev.new()                         	#open new graphic window
par(mfrow=c(1, 2))                	#split window in two
hist(smallbin$sum_col[smallbin$swipes=="40"])
hist(smallbin$sum_col[smallbin$swipes=="80"])

#both reads and species
boxplot(smallbin$sum_col~smallbin$swipes)
boxplot(smallreads$sum_col~smallreads$swipes)

######################STUDENT TESTS FOR NUMBER OF SWIPES

#On the number of species found for each sample
t.test(smallbin$sum_col[smallbin$swipes=="40"],y=smallbin$sum_col[smallbin$swipes=="80"])  

#On the number of reads for each sample
t.test(smallreads$sum_col[smallreads$swipes=="40"],y=smallreads$sum_col[smallreads$swipes=="80"])  


###########ANOVA ON SAMPLER

#check hypotheses first

res_aov <- aov(smallbin$sum_col ~ smallbin$sampler,
               data = smallbin)

#test normality residuals : 
par(mfrow = c(1, 2)) # combine plots

hist(res_aov$residuals) # histogram ##NORMALITY ASSUMPTION NOT MET

library(car)
qqPlot(res_aov$residuals, # QQ-plot
       id = FALSE) # id = FALSE to remove point identification
#not really a straight line either

#We can't do an anova ): let's try kruskal-wallis
kruskal.test(smallbin$sum_col ~ smallbin$sampler)
#no significative diff between samplers


#######CHI2 tests on presence/absence in sample

#on swipes
table <- table(smallbin$colbin, smallbin$swipes)
chisq.test(table)
#pvalue=0,48

