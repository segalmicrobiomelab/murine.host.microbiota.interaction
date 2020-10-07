### Scrits for MiSeq data for Host Imunne Response to Dysbiosis and Pathogen Challenge

#Load Phyloseq
library(phyloseq)
library(ade4)
library(vegan)
library(biomformat)
library(devtools)
library(readr)
library(readtext)
library(qiime2R)
library("ape")
library(phangorn)

install.packages('devtools')
install.packages('tidyverse')
install.packages('readr')
install.packages('readtext')
install.packages('vegan')
install.packages('ade4')
install.packages('biomformat')
install.packages('phyloseq')
install.packages("devtools") #if you don't have devtools installed
install.packages("ape")
install.packages("phangorn")
devtools::install_github("jbisanz/qiime2R")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("phyloseq")

update.packages()

#################################################################
#################################################################
#################################################################

save.image(file="mouse.RData")
load(file="mouse.RData")

#################################################################
#################################################################
#################################################################

sample_data(phy.relative.table)$Sample

# Base Lung samples 
Lung.phy.table = subset_samples(phy.relative.table, Sample %in% c('Lung'))

# Filter out only experiments of interest
Lung.timeexperiments.phy.table = subset_samples(Lung.phy.table, Experiment %in% c('Shortexpsoure.Single.Day.2', 'Platelet.single.D5', 'MOC.single.D14', 'Platelet.single.D5.repeat', 'Shortexpsoure.Single.Day.3', 'Shortexpsoure.Single.Day.2.repeat'))

# Filter out only MOC and PBS group
MOC.PBS.Lung.timeexperiments.phy.table = subset_samples(Lung.timeexperiments.phy.table, Exposure %in% c('MOC', 'PBS'))

#################################################################
#################################################################
#################################################################

library(DESeq2)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(pathfindR)
library(scales)
library(data.table)
library(fBasics)
library(forcats)
library(omu)
library(maptools)
library(phyloseq)
library(SpiecEasi)
library("vegan")
library(ade4)
library("reshape2")
library(dplyr)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Glimma")
install.packages('gplots')
install.packages('RColorBrewer')
install.packages('pheatmap')
install.packages('ggplot2')
install.packages('ggrepel')
install.packages('pathfindR')
install.packages('scales')
install.packages('data.table')
install.packages('fBasics')
install.packages('forcats')
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("KEGGREST")
install.packages('omu')
install_github("connor-reid-tiffany/Omu.git")
 *** need to re-install 
install.packages('maptools')
install.packages(phyloseq)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SpiecEasi")
library(devtools)
install_github("zdk123/SpiecEasi")
install.packages("vegan")
install.packages('ade4')
install.packages("reshape2")
install.packages('dplyr')
install.packages('ggthemes')

#How to plot PCOA with labels
bray.dist = phyloseq::distance(MOC.PBS.Lung.timeexperiments.phy.table, method = "bray")

### estimate number of axesb
bray.pco = dudi.pco(cailliez(bray.dist), scannf = FALSE, nf = 3)

### Code with all samples: 

#Create Distance Matrix with Bray (or wUniFrac depending what you are using)
vegdist   = vegdist(t(otu_table(MOC.PBS.Lung.timeexperiments.phy.table)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = sample_data(MOC.PBS.Lung.timeexperiments.phy.table), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Day.Exposure,data= newResults, mean) #Here you would use your grouping variable (e.g., days)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Day.Exposure",suffixes=c("",".centroid")) #Here you would use your grouping variable (e.g., days)

pdf("200119_Murine.exposure.day.MOC.PBS.pdf", height = 10, width = 10)
    ggplot(newResults, aes(PC1, PC2, color= Day.Exposure)) + # Graph PC1 and PC2
    geom_point(size=5) + # Set the size of the points
    xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
    #Set colors for each category, should be the same number of categories
    scale_color_manual(values=c("red", "dodgerblue", "darkred", "gold2", "navy", "darkgreen")) + 
    #plot point and lines from centroid
    geom_point(data=centroids, aes(x=PC1, y=PC2, color=Day.Exposure), size=0) +  #Here you would use your grouping variable (e.g., days)
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Day.Exposure))+  #Here you would use your grouping variable (e.g., days)
    #If you want to identify specific samples use the code bellow
    #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BAL.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
    #labels centroids should be same number of categories
    geom_label_repel(data = centroids, aes(x=PC1, y=PC2), label=c('D1.MOC', 'D5.MOC', 'D2.MOC', 'D3.MOC', 'D14.MOC', 'PBS'), size=10) +  #Here you can label the way you want
    #Use the code bellow if you want to switch the X axis around
    scale_x_reverse() +
    theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

adonis(vegdist  ~ Day.Exposure, data=data.frame(sample_data(MOC.PBS.Lung.timeexperiments.phy.table)))

Call:
adonis(formula = vegdist ~ Day.Exposure, data = data.frame(sample_data(MOC.PBS.Lung.timeexperiments.phy.table))) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Day.Exposure  5    8.7688 1.75376   7.734 0.33723  0.001 ***
Residuals    76   17.2339 0.22676         0.66277           
Total        81   26.0027                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#################################################################
#################################################################
#################################################################

save.image(file="mouse.RData")
load(file="mouse.RData")

#################################################################
#################################################################
#################################################################

#Load libraries 
library(phyloseq)
library(ade4)
library("vegan")
library("ggplot2")
library("ggthemes")
library("reshape2")

### Intergroup diversity 1, 2, 3, 5, 14 day all exposures Lung 
Nasal.phy.table = subset_samples(phy.relative.table, Sample %in% c('Nasal'))
Oral.phy.table = subset_samples(phy.relative.table, Sample %in% c('Oral'))
Cecum.phy.table = subset_samples(phy.relative.table, Sample %in% c('Cecum'))
MOC.PBS.Lung.timeexperiments.phy.table

Nasal.experiments.phy.table = subset_samples(Nasal.phy.table, Experiment %in% c('Shortexpsoure.Single.Day.2', 'Platelet.single.D5', 'MOC.single.D14', 'Platelet.single.D5.repeat', 'Shortexpsoure.Single.Day.3', 'Shortexpsoure.Single.Day.2.repeat'))
Oral.experiments.phy.table = subset_samples(Oral.phy.table, Experiment %in% c('Shortexpsoure.Single.Day.2', 'Platelet.single.D5', 'MOC.single.D14', 'Platelet.single.D5.repeat', 'Shortexpsoure.Single.Day.3', 'Shortexpsoure.Single.Day.2.repeat'))
Cecum.experiments.phy.table = subset_samples(Cecum.phy.table, Experiment %in% c('Shortexpsoure.Single.Day.2', 'Platelet.single.D5', 'MOC.single.D14', 'Platelet.single.D5.repeat', 'Shortexpsoure.Single.Day.3', 'Shortexpsoure.Single.Day.2.repeat'))

# Filter out only experiments of interest
Nasal.timeexperiments.phy.table = subset_samples(Nasal.experiments.phy.table, Experiment %in% c('Shortexpsoure.Single.Day.2', 'Platelet.single.D5', 'MOC.single.D14', 'Platelet.single.D5.repeat', 'Shortexpsoure.Single.Day.3', 'Shortexpsoure.Single.Day.2.repeat'))
Oral.timeexperiments.phy.table = subset_samples(Oral.experiments.phy.table, Experiment %in% c('Shortexpsoure.Single.Day.2', 'Platelet.single.D5', 'MOC.single.D14', 'Platelet.single.D5.repeat', 'Shortexpsoure.Single.Day.3', 'Shortexpsoure.Single.Day.2.repeat'))
Cecum.timeexperiments.phy.table = subset_samples(Cecum.experiments.phy.table, Experiment %in% c('Shortexpsoure.Single.Day.2', 'Platelet.single.D5', 'MOC.single.D14', 'Platelet.single.D5.repeat', 'Shortexpsoure.Single.Day.3', 'Shortexpsoure.Single.Day.2.repeat'))

# Filter out only MOC and PBS group
MOC.PBS.Nasal.timeexperiments.phy.table = subset_samples(Nasal.timeexperiments.phy.table, Exposure %in% c('MOC', 'PBS'))
MOC.PBS.Oral.timeexperiments.phy.table = subset_samples(Oral.timeexperiments.phy.table, Exposure %in% c('MOC', 'PBS'))
MOC.PBS.Cecum.timeexperiments.phy.table = subset_samples(Cecum.timeexperiments.phy.table, Exposure %in% c('MOC', 'PBS'))

MOC.PBS.Nasal.timeexperiments.phy.table

sample_data(MOC.PBS.Nasal.timeexperiments.phy.table)$Day.Exposure
sample_data(MOC.PBS.Oral.timeexperiments.phy.table)$Day.Exposure
sample_data(MOC.PBS.Cecum.timeexperiments.phy.table)$Day.Exposure

MOC.PBS.Oral.timeexperiments.phy.table
MOC.PBS.Cecum.timeexperiments.phy.table 
MOC.PBS.Lung.timeexperiments.phy.table

#Calculate Distance - you will use this matrix 
wunif.Lung.Time.dist = phyloseq::distance(MOC.PBS.Lung.timeexperiments.phy.table, method="bray")

#estimate number of axes
wunif.Lung.Time.pco = dudi.pco(cailliez(wunif.Lung.Time.dist), scannf = FALSE, nf = 3)

#Plot first to see colors
s.class(wunif.Lung.Time.pco $li, sample_data(MOC.PBS.Lung.timeexperiments.phy.table)$Day.Exposure, col=c("grey", "blue", "grey", "grey", "grey", "darkgreen", "orange", "grey", "grey", "red", "darkred"))

# Subset for patient's sample, removing background
# Patient.Sample.otu.relative.table = subset_samples(Lung.Time.Strep.PBS.MOC.otu.relative.table, topograph !="BKG")
# sample_data(Patient.Sample.otu.relative.table)$topograph

#############
### Statistics functions ###
#We want to create these functions in order to plot on top of our plot
# Create function to calculate the following statistics: mean, std deviation, median, min value, max value, 10%ile, 25%ile, 75%ile, 90%ile
stats.all = function(x) {
  mean <- mean(x)
  stddev <- sd(x)
  median <- median(x)
  val_min <- min(x)
  val_max <- max(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(mean = mean, sd = stddev, median = median, 
           val_min = val_min, val_max = val_max, 
           per10 = per10, per25 = per25, per75 = per75,  per90 = per90))
}

#Create function for boxplot statistics: median, 25%ile, 75%ile
stats.boxplot <- function(x) {
  m <- median(x)
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  return(c(y = m, ymin = per25, ymax = per75))
}

#Create function for whiskers statistics: median, min value, max value
stats.whiskers = function(x) {
  m <- median(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(y = m, ymin = per10, ymax = per90))
}

###################
#Set up for plotting 
# Colors
col1 <- c("green", "blue", "darkmagenta ", "cadetblue1 ", "purple", "red", "orange")
my.cols = col1
#Labels
label1 <- c("MOC.Day5", "PBS.all", "MOC.Day3", "MOC.Day1", "MOC.Day2", "MOC.Day14")
mylab = label1
#theme
theme_bw()


###################
#Begin code to plot 

#melt the distance matrix into columns with each sample site and distance
b <- melt(as.matrix(wunif.Lung.Time.dist))


#Then need to remove self distances and duplicated distances
p    <- t(apply(b[,c(1,2)],1,FUN=sort))
rmv1 <- which(p[,1] == p[,2])

p    <- paste(p[,1],p[,2],sep="|")
rmv2 <- which(duplicated(p))

#establish new data frame that removes those values 
b.df   <- b[-c(rmv1,rmv2),] 


##Now we need to replace variable rows with topo 
#set up new data frame
new.df <- b.df
#create new data frame with the variable you want to use for group comparisons. This code goes through the distance data frame and makes columns for with the group for each sample
new.df[] <- lapply(b.df, function(x) sample_data(MOC.PBS.Lung.timeexperiments.phy.table)$Day.Exposure [match(x, rownames(sample_data(MOC.PBS.Lung.timeexperiments.phy.table)))])

#create two lists of the group variable 
topo.var1 <- new.df[,1]
topo.var2 <-new.df[,2]

#Add the two columns of group variable data onto the end of the distance data frame 
b.var.df <- cbind(b.df, topo.var1, topo.var2)


##We will now need to make sure we do not have intra groups, so we will remove those 
#create new data frame 
btw.b.var.df <- b.var.df
#set row names to re-zero
rownames(btw.b.var.df) <- NULL

#establish matrix for input of indexes to be removed 
toremove<-numeric()

#select indexes to remove 
for (i in 1:nrow(btw.b.var.df)) {
	if (btw.b.var.df$topo.var1[i] == btw.b.var.df$topo.var2[i]) {
		toremove <- append(toremove, i)

	} 
}

#remove indexes we selected
btw.b.var.df <- btw.b.var.df[-toremove,]

#Now the intragroup should be removed, can confirm 
head(btw.b.var.df)

##Now we need to see that the between groups we have are not reverse permutations of each other
#Use the two group categories two create a new category of the permutation and set as a data frame
new.cat.btw = paste(btw.b.var.df$topo.var1, "to", btw.b.var.df$topo.var2)
new.cat.btw.df <- data.frame(btw.b.var.df$topo.var1, btw.b.var.df$topo.var2)

#create a list of combinations from our specific data frame and select for the unique ones for comparison
dat.sort = t(apply(new.cat.btw.df, 1, sort))
unique.new.cat.btw <- unique(new.cat.btw.df[!duplicated(dat.sort),])
colnames(unique.new.cat.btw) <- NULL
rownames(unique.new.cat.btw) <- NULL
unique.new.cat.btw <- paste(unique.new.cat.btw[,1], "to", unique.new.cat.btw[,2])


#create new data frame 
clean.btw.b.var.df <- btw.b.var.df

#reset row names
rownames(clean.btw.b.var.df) <- NULL

#this code checks if any of the reverse combinations exist in the unique list of permutations and will reverse them if so. Reversing them allows them to be plotted as one group rather than deleting any data
for (i in 1:nrow(clean.btw.b.var.df)){
	if (paste(clean.btw.b.var.df$topo.var2[i], "to", clean.btw.b.var.df$topo.var1[i]) %in% unique.new.cat.btw) {
		clean.btw.b.var.df$topo.var1[i] <- btw.b.var.df$topo.var2[i]
		clean.btw.b.var.df$topo.var2[i] <- btw.b.var.df$topo.var1[i]	
	}
}



#Use the two new categories two create a new category of the permutation without the doubles
new.cat.btw.clean = paste(clean.btw.b.var.df$topo.var1, "to", clean.btw.b.var.df$topo.var2)

#confirm permutations 
unique(new.cat.btw.clean)
  
# create all-encompassing data frame using the names, category and distance data, Var1 and Var2 are the names of the columns with the sampleIDs
comb.inter.data.braypart <- data.frame(Sample1 = clean.btw.b.var.df$Var1, Sample2 = clean.btw.b.var.df$Var2, Category = new.cat.btw.clean, Distance = as.numeric(clean.btw.b.var.df$value))

#plot 
pdf("bray_interbeta_lung.pdf", height = 10, width = 15)
ggplot(
	data= comb.inter.data.braypart,	
	#set the aesthetics based on the names from the data frame created above
	aes(x= Category, y= Distance, fill = Category, colour = Category)) + 
	#defines data as discrete x,y
	geom_jitter(size = 0.8, color = "gray34") +
	stat_summary(fun.data = stats.whiskers, geom = "errorbar", color = "black", size = 0.8, width = 0.6) +
 	 #adds boxplot statistics
 	 stat_summary(fun.data = stats.boxplot, geom = "crossbar", color = "black", size = 0.5, width = 0.5) +
	#set labels 
	xlab(NULL) + 
	ylab("Weighted UniFrac Distance") + 
	#change theme - This theme is based on what we like for these graphs. feel free to edit it as needed, but try baseline first
   theme_bw() + 
   theme(axis.text=element_text(size=8, angle=90, color = "black"), panel.border = element_blank(), axis.line =element_line(colour = "black", size = 1, linetype = 1, lineend = "butt"), axis.ticks = element_line(size=1), axis.ticks.length = unit(5, "pt"),  panel.grid=element_blank(),  axis.title.y = element_text(color="black" , size = 16), axis.text.y = element_text( color = "black", size = 14, angle=0))
dev.off()


#Save data as table
write.table(comb.inter.data.braypart, file="Inter_Beta_diversity_Lung_wuniF.txt", sep="\t")

#################################################################
#################################################################
#################################################################

save.image(file="mouse.RData")
load(file="mouse.RData")

#################################################################
#################################################################
#################################################################
#################################################################

#Calculate Distance - you will use this matrix 
wunif.Nasal.Time.dist = phyloseq::distance(MOC.PBS.Nasal.timeexperiments.phy.table, method="bray")

#estimate number of axes
wunif.Nasal.Time.pco = dudi.pco(cailliez(wunif.Nasal.Time.dist), scannf = FALSE, nf = 3)

#Plot first to see colors
s.class(wunif.Nasal.Time.pco $li, sample_data(MOC.PBS.Nasal.timeexperiments.phy.table)$Day.Exposure, col=c("grey", "blue", "grey", "grey", "grey", "darkgreen", "orange", "grey", "grey", "red", "darkred"))

table(sample_data(MOC.PBS.Nasal.timeexperiments.phy.table)$Day.Exposure)

#############
### Statistics functions ###
#We want to create these functions in order to plot on top of our plot
# Create function to calculate the following statistics: mean, std deviation, median, min value, max value, 10%ile, 25%ile, 75%ile, 90%ile
stats.all = function(x) {
  mean <- mean(x)
  stddev <- sd(x)
  median <- median(x)
  val_min <- min(x)
  val_max <- max(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(mean = mean, sd = stddev, median = median, 
           val_min = val_min, val_max = val_max, 
           per10 = per10, per25 = per25, per75 = per75,  per90 = per90))
}

#Create function for boxplot statistics: median, 25%ile, 75%ile
stats.boxplot <- function(x) {
  m <- median(x)
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  return(c(y = m, ymin = per25, ymax = per75))
}

#Create function for whiskers statistics: median, min value, max value
stats.whiskers = function(x) {
  m <- median(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(y = m, ymin = per10, ymax = per90))
}

###################
#Set up for plotting 
# Colors
col1 <- c("green", "blue", "darkmagenta ", "cadetblue1 ", "purple", "red", "orange")
my.cols = col1
#Labels
label1 <- c("MOC.Day5", "PBS.all", "MOC.Day3", "MOC.Day1", "MOC.Day2", "MOC.Day14")
mylab = label1
#theme
theme_bw()


###################
#Begin code to plot 

#melt the distance matrix into columns with each sample site and distance
b <- melt(as.matrix(wunif.Nasal.Time.dist))


#Then need to remove self distances and duplicated distances
p    <- t(apply(b[,c(1,2)],1,FUN=sort))
rmv1 <- which(p[,1] == p[,2])

p    <- paste(p[,1],p[,2],sep="|")
rmv2 <- which(duplicated(p))

#establish new data frame that removes those values 
b.df   <- b[-c(rmv1,rmv2),] 


##Now we need to replace variable rows with topo 
#set up new data frame
new.df <- b.df
#create new data frame with the variable you want to use for group comparisons. This code goes through the distance data frame and makes columns for with the group for each sample
new.df[] <- lapply(b.df, function(x) sample_data(MOC.PBS.Nasal.timeexperiments.phy.table)$Day.Exposure [match(x, rownames(sample_data(MOC.PBS.Nasal.timeexperiments.phy.table)))])

#create two lists of the group variable 
topo.var1 <- new.df[,1]
topo.var2 <-new.df[,2]

#Add the two columns of group variable data onto the end of the distance data frame 
b.var.df <- cbind(b.df, topo.var1, topo.var2)


##We will now need to make sure we do not have intra groups, so we will remove those 
#create new data frame 
btw.b.var.df <- b.var.df
#set row names to re-zero
rownames(btw.b.var.df) <- NULL

#establish matrix for input of indexes to be removed 
toremove<-numeric()

#select indexes to remove 
for (i in 1:nrow(btw.b.var.df)) {
	if (btw.b.var.df$topo.var1[i] == btw.b.var.df$topo.var2[i]) {
		toremove <- append(toremove, i)

	} 
}

#remove indexes we selected
btw.b.var.df <- btw.b.var.df[-toremove,]

#Now the intragroup should be removed, can confirm 
head(btw.b.var.df)

##Now we need to see that the between groups we have are not reverse permutations of each other
#Use the two group categories two create a new category of the permutation and set as a data frame
new.cat.btw = paste(btw.b.var.df$topo.var1, "to", btw.b.var.df$topo.var2)
new.cat.btw.df <- data.frame(btw.b.var.df$topo.var1, btw.b.var.df$topo.var2)

#create a list of combinations from our specific data frame and select for the unique ones for comparison
dat.sort = t(apply(new.cat.btw.df, 1, sort))
unique.new.cat.btw <- unique(new.cat.btw.df[!duplicated(dat.sort),])
colnames(unique.new.cat.btw) <- NULL
rownames(unique.new.cat.btw) <- NULL
unique.new.cat.btw <- paste(unique.new.cat.btw[,1], "to", unique.new.cat.btw[,2])


#create new data frame 
clean.btw.b.var.df <- btw.b.var.df

#reset row names
rownames(clean.btw.b.var.df) <- NULL

#this code checks if any of the reverse combinations exist in the unique list of permutations and will reverse them if so. Reversing them allows them to be plotted as one group rather than deleting any data
for (i in 1:nrow(clean.btw.b.var.df)){
	if (paste(clean.btw.b.var.df$topo.var2[i], "to", clean.btw.b.var.df$topo.var1[i]) %in% unique.new.cat.btw) {
		clean.btw.b.var.df$topo.var1[i] <- btw.b.var.df$topo.var2[i]
		clean.btw.b.var.df$topo.var2[i] <- btw.b.var.df$topo.var1[i]	
	}
}



#Use the two new categories two create a new category of the permutation without the doubles
new.cat.btw.clean = paste(clean.btw.b.var.df$topo.var1, "to", clean.btw.b.var.df$topo.var2)

#confirm permutations 
unique(new.cat.btw.clean)
  
# create all-encompassing data frame using the names, category and distance data, Var1 and Var2 are the names of the columns with the sampleIDs
comb.inter.data.braypart <- data.frame(Sample1 = clean.btw.b.var.df$Var1, Sample2 = clean.btw.b.var.df$Var2, Category = new.cat.btw.clean, Distance = as.numeric(clean.btw.b.var.df$value))

#plot 
pdf("bray_interbeta_nasal.pdf", height = 10, width = 15)
ggplot(
	data= comb.inter.data.braypart,	
	#set the aesthetics based on the names from the data frame created above
	aes(x= Category, y= Distance, fill = Category, colour = Category)) + 
	#defines data as discrete x,y
	geom_jitter(size = 0.8, color = "gray34") +
	stat_summary(fun.data = stats.whiskers, geom = "errorbar", color = "black", size = 0.8, width = 0.6) +
 	 #adds boxplot statistics
 	 stat_summary(fun.data = stats.boxplot, geom = "crossbar", color = "black", size = 0.5, width = 0.5) +
	#set labels 
	xlab(NULL) + 
	ylab("Weighted UniFrac Distance") + 
	#change theme - This theme is based on what we like for these graphs. feel free to edit it as needed, but try baseline first
   theme_bw() + 
   theme(axis.text=element_text(size=8, angle=90, color = "black"), panel.border = element_blank(), axis.line =element_line(colour = "black", size = 1, linetype = 1, lineend = "butt"), axis.ticks = element_line(size=1), axis.ticks.length = unit(5, "pt"),  panel.grid=element_blank(),  axis.title.y = element_text(color="black" , size = 16), axis.text.y = element_text( color = "black", size = 14, angle=0))
dev.off()


#Save data as table
write.table(comb.inter.data.braypart, file="Inter_Beta_diversity_nasal_wuniF.txt", sep="\t")

#################################################################
#################################################################
#################################################################

save.image(file="mouse.RData")
load(file="mouse.RData")

#################################################################
#################################################################
#################################################################
#################################################################

#Calculate Distance - you will use this matrix 
wunif.Oral.Time.dist = phyloseq::distance(MOC.PBS.Oral.timeexperiments.phy.table, method="bray")

#estimate number of axes
wunif.Oral.Time.pco = dudi.pco(cailliez(wunif.Oral.Time.dist), scannf = FALSE, nf = 3)

#Plot first to see colors
s.class(wunif.Oral.Time.pco $li, sample_data(MOC.PBS.Oral.timeexperiments.phy.table)$Day.Exposure, col=c("grey", "blue", "grey", "grey", "grey", "darkgreen", "orange", "grey", "grey", "red", "darkred"))

#############
### Statistics functions ###
#We want to create these functions in order to plot on top of our plot
# Create function to calculate the following statistics: mean, std deviation, median, min value, max value, 10%ile, 25%ile, 75%ile, 90%ile
stats.all = function(x) {
  mean <- mean(x)
  stddev <- sd(x)
  median <- median(x)
  val_min <- min(x)
  val_max <- max(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(mean = mean, sd = stddev, median = median, 
           val_min = val_min, val_max = val_max, 
           per10 = per10, per25 = per25, per75 = per75,  per90 = per90))
}

#Create function for boxplot statistics: median, 25%ile, 75%ile
stats.boxplot <- function(x) {
  m <- median(x)
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  return(c(y = m, ymin = per25, ymax = per75))
}

#Create function for whiskers statistics: median, min value, max value
stats.whiskers = function(x) {
  m <- median(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(y = m, ymin = per10, ymax = per90))
}

###################
#Set up for plotting 
# Colors
col1 <- c("green", "blue", "darkmagenta ", "cadetblue1 ", "purple", "red", "orange")
my.cols = col1
#Labels
label1 <- c("MOC.Day5", "PBS.all", "MOC.Day3", "MOC.Day1", "MOC.Day2", "MOC.Day14")
mylab = label1
#theme
theme_bw()


###################
#Begin code to plot 

#melt the distance matrix into columns with each sample site and distance
b <- melt(as.matrix(wunif.Oral.Time.dist))


#Then need to remove self distances and duplicated distances
p    <- t(apply(b[,c(1,2)],1,FUN=sort))
rmv1 <- which(p[,1] == p[,2])

p    <- paste(p[,1],p[,2],sep="|")
rmv2 <- which(duplicated(p))

#establish new data frame that removes those values 
b.df   <- b[-c(rmv1,rmv2),] 


##Now we need to replace variable rows with topo 
#set up new data frame
new.df <- b.df
#create new data frame with the variable you want to use for group comparisons. This code goes through the distance data frame and makes columns for with the group for each sample
new.df[] <- lapply(b.df, function(x) sample_data(MOC.PBS.Oral.timeexperiments.phy.table)$Day.Exposure [match(x, rownames(sample_data(MOC.PBS.Oral.timeexperiments.phy.table)))])

#create two lists of the group variable 
topo.var1 <- new.df[,1]
topo.var2 <-new.df[,2]

#Add the two columns of group variable data onto the end of the distance data frame 
b.var.df <- cbind(b.df, topo.var1, topo.var2)


##We will now need to make sure we do not have intra groups, so we will remove those 
#create new data frame 
btw.b.var.df <- b.var.df
#set row names to re-zero
rownames(btw.b.var.df) <- NULL

#establish matrix for input of indexes to be removed 
toremove<-numeric()

#select indexes to remove 
for (i in 1:nrow(btw.b.var.df)) {
	if (btw.b.var.df$topo.var1[i] == btw.b.var.df$topo.var2[i]) {
		toremove <- append(toremove, i)

	} 
}

#remove indexes we selected
btw.b.var.df <- btw.b.var.df[-toremove,]

#Now the intragroup should be removed, can confirm 
head(btw.b.var.df)

##Now we need to see that the between groups we have are not reverse permutations of each other
#Use the two group categories two create a new category of the permutation and set as a data frame
new.cat.btw = paste(btw.b.var.df$topo.var1, "to", btw.b.var.df$topo.var2)
new.cat.btw.df <- data.frame(btw.b.var.df$topo.var1, btw.b.var.df$topo.var2)

#create a list of combinations from our specific data frame and select for the unique ones for comparison
dat.sort = t(apply(new.cat.btw.df, 1, sort))
unique.new.cat.btw <- unique(new.cat.btw.df[!duplicated(dat.sort),])
colnames(unique.new.cat.btw) <- NULL
rownames(unique.new.cat.btw) <- NULL
unique.new.cat.btw <- paste(unique.new.cat.btw[,1], "to", unique.new.cat.btw[,2])


#create new data frame 
clean.btw.b.var.df <- btw.b.var.df

#reset row names
rownames(clean.btw.b.var.df) <- NULL

#this code checks if any of the reverse combinations exist in the unique list of permutations and will reverse them if so. Reversing them allows them to be plotted as one group rather than deleting any data
for (i in 1:nrow(clean.btw.b.var.df)){
	if (paste(clean.btw.b.var.df$topo.var2[i], "to", clean.btw.b.var.df$topo.var1[i]) %in% unique.new.cat.btw) {
		clean.btw.b.var.df$topo.var1[i] <- btw.b.var.df$topo.var2[i]
		clean.btw.b.var.df$topo.var2[i] <- btw.b.var.df$topo.var1[i]	
	}
}



#Use the two new categories two create a new category of the permutation without the doubles
new.cat.btw.clean = paste(clean.btw.b.var.df$topo.var1, "to", clean.btw.b.var.df$topo.var2)

#confirm permutations 
unique(new.cat.btw.clean)
  
# create all-encompassing data frame using the names, category and distance data, Var1 and Var2 are the names of the columns with the sampleIDs
comb.inter.data.braypart <- data.frame(Sample1 = clean.btw.b.var.df$Var1, Sample2 = clean.btw.b.var.df$Var2, Category = new.cat.btw.clean, Distance = as.numeric(clean.btw.b.var.df$value))

#plot 
pdf("bray_interbeta_noral.pdf", height = 10, width = 15)
ggplot(
	data= comb.inter.data.braypart,	
	#set the aesthetics based on the names from the data frame created above
	aes(x= Category, y= Distance, fill = Category, colour = Category)) + 
	#defines data as discrete x,y
	geom_jitter(size = 0.8, color = "gray34") +
	stat_summary(fun.data = stats.whiskers, geom = "errorbar", color = "black", size = 0.8, width = 0.6) +
 	 #adds boxplot statistics
 	 stat_summary(fun.data = stats.boxplot, geom = "crossbar", color = "black", size = 0.5, width = 0.5) +
	#set labels 
	xlab(NULL) + 
	ylab("Bray Distance") + 
	#change theme - This theme is based on what we like for these graphs. feel free to edit it as needed, but try baseline first
   theme_bw() + 
   theme(axis.text=element_text(size=8, angle=90, color = "black"), panel.border = element_blank(), axis.line =element_line(colour = "black", size = 1, linetype = 1, lineend = "butt"), axis.ticks = element_line(size=1), axis.ticks.length = unit(5, "pt"),  panel.grid=element_blank(),  axis.title.y = element_text(color="black" , size = 16), axis.text.y = element_text( color = "black", size = 14, angle=0))
dev.off()


#Save data as table
write.table(comb.inter.data.braypart, file="Inter_Beta_diversity_oral_wuniF.txt", sep="\t")

#################################################################
#################################################################
#################################################################

#Calculate Distance - you will use this matrix 
wunif.Cecum.Time.dist = phyloseq::distance(MOC.PBS.Cecum.timeexperiments.phy.table, method="bray")

#estimate number of axes
wunif.Cecum.Time.pco = dudi.pco(cailliez(wunif.Cecun.Time.dist), scannf = FALSE, nf = 3)

#Plot first to see colors
s.class(wunif.Cecum.Time.pco $li, sample_data(MOC.PBS.Cecum.timeexperiments.phy.table)$Day.Exposure, col=c("grey", "blue", "grey", "grey", "grey", "darkgreen", "orange", "grey", "grey", "red", "darkred"))

#############
### Statistics functions ###
#We want to create these functions in order to plot on top of our plot
# Create function to calculate the following statistics: mean, std deviation, median, min value, max value, 10%ile, 25%ile, 75%ile, 90%ile
stats.all = function(x) {
  mean <- mean(x)
  stddev <- sd(x)
  median <- median(x)
  val_min <- min(x)
  val_max <- max(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(mean = mean, sd = stddev, median = median, 
           val_min = val_min, val_max = val_max, 
           per10 = per10, per25 = per25, per75 = per75,  per90 = per90))
}

#Create function for boxplot statistics: median, 25%ile, 75%ile
stats.boxplot <- function(x) {
  m <- median(x)
  per25 <- as.numeric(quantile(x, prob = c(0.25)))
  per75 <- as.numeric(quantile(x, prob = c(0.75)))
  return(c(y = m, ymin = per25, ymax = per75))
}

#Create function for whiskers statistics: median, min value, max value
stats.whiskers = function(x) {
  m <- median(x)
  per10 <- as.numeric(quantile(x, prob = c(0.10)))
  per90 <- as.numeric(quantile(x, prob = c(0.90)))
  return(c(y = m, ymin = per10, ymax = per90))
}

###################
#Set up for plotting 
# Colors
col1 <- c("green", "blue", "darkmagenta ", "cadetblue1 ", "purple", "red", "orange")
my.cols = col1
#Labels
label1 <- c("MOC.Day5", "PBS.all", "MOC.Day3", "MOC.Day1", "MOC.Day2", "MOC.Day14")
mylab = label1
#theme
theme_bw()


###################
#Begin code to plot 

#melt the distance matrix into columns with each sample site and distance
b <- melt(as.matrix(wunif.Cecum.Time.dist))


#Then need to remove self distances and duplicated distances
p    <- t(apply(b[,c(1,2)],1,FUN=sort))
rmv1 <- which(p[,1] == p[,2])

p    <- paste(p[,1],p[,2],sep="|")
rmv2 <- which(duplicated(p))

#establish new data frame that removes those values 
b.df   <- b[-c(rmv1,rmv2),] 


##Now we need to replace variable rows with topo 
#set up new data frame
new.df <- b.df
#create new data frame with the variable you want to use for group comparisons. This code goes through the distance data frame and makes columns for with the group for each sample
new.df[] <- lapply(b.df, function(x) sample_data(MOC.PBS.Cecum.timeexperiments.phy.table)$Day.Exposure [match(x, rownames(sample_data(MOC.PBS.Cecum.timeexperiments.phy.table)))])

#create two lists of the group variable 
topo.var1 <- new.df[,1]
topo.var2 <-new.df[,2]

#Add the two columns of group variable data onto the end of the distance data frame 
b.var.df <- cbind(b.df, topo.var1, topo.var2)


##We will now need to make sure we do not have intra groups, so we will remove those 
#create new data frame 
btw.b.var.df <- b.var.df
#set row names to re-zero
rownames(btw.b.var.df) <- NULL

#establish matrix for input of indexes to be removed 
toremove<-numeric()

#select indexes to remove 
for (i in 1:nrow(btw.b.var.df)) {
	if (btw.b.var.df$topo.var1[i] == btw.b.var.df$topo.var2[i]) {
		toremove <- append(toremove, i)

	} 
}

#remove indexes we selected
btw.b.var.df <- btw.b.var.df[-toremove,]

#Now the intragroup should be removed, can confirm 
head(btw.b.var.df)

##Now we need to see that the between groups we have are not reverse permutations of each other
#Use the two group categories two create a new category of the permutation and set as a data frame
new.cat.btw = paste(btw.b.var.df$topo.var1, "to", btw.b.var.df$topo.var2)
new.cat.btw.df <- data.frame(btw.b.var.df$topo.var1, btw.b.var.df$topo.var2)

#create a list of combinations from our specific data frame and select for the unique ones for comparison
dat.sort = t(apply(new.cat.btw.df, 1, sort))
unique.new.cat.btw <- unique(new.cat.btw.df[!duplicated(dat.sort),])
colnames(unique.new.cat.btw) <- NULL
rownames(unique.new.cat.btw) <- NULL
unique.new.cat.btw <- paste(unique.new.cat.btw[,1], "to", unique.new.cat.btw[,2])


#create new data frame 
clean.btw.b.var.df <- btw.b.var.df

#reset row names
rownames(clean.btw.b.var.df) <- NULL

#this code checks if any of the reverse combinations exist in the unique list of permutations and will reverse them if so. Reversing them allows them to be plotted as one group rather than deleting any data
for (i in 1:nrow(clean.btw.b.var.df)){
	if (paste(clean.btw.b.var.df$topo.var2[i], "to", clean.btw.b.var.df$topo.var1[i]) %in% unique.new.cat.btw) {
		clean.btw.b.var.df$topo.var1[i] <- btw.b.var.df$topo.var2[i]
		clean.btw.b.var.df$topo.var2[i] <- btw.b.var.df$topo.var1[i]	
	}
}



#Use the two new categories two create a new category of the permutation without the doubles
new.cat.btw.clean = paste(clean.btw.b.var.df$topo.var1, "to", clean.btw.b.var.df$topo.var2)

#confirm permutations 
unique(new.cat.btw.clean)
  
# create all-encompassing data frame using the names, category and distance data, Var1 and Var2 are the names of the columns with the sampleIDs
comb.inter.data.braypart <- data.frame(Sample1 = clean.btw.b.var.df$Var1, Sample2 = clean.btw.b.var.df$Var2, Category = new.cat.btw.clean, Distance = as.numeric(clean.btw.b.var.df$value))

#plot 
pdf("bray_interbeta_cecum.pdf", height = 10, width = 15)
ggplot(
	data= comb.inter.data.braypart,	
	#set the aesthetics based on the names from the data frame created above
	aes(x= Category, y= Distance, fill = Category, colour = Category)) + 
	#defines data as discrete x,y
	geom_jitter(size = 0.8, color = "gray34") +
	stat_summary(fun.data = stats.whiskers, geom = "errorbar", color = "black", size = 0.8, width = 0.6) +
 	 #adds boxplot statistics
 	 stat_summary(fun.data = stats.boxplot, geom = "crossbar", color = "black", size = 0.5, width = 0.5) +
	#set labels 
	xlab(NULL) + 
	ylab("Bray Distance") + 
	#change theme - This theme is based on what we like for these graphs. feel free to edit it as needed, but try baseline first
   theme_bw() + 
   theme(axis.text=element_text(size=8, angle=90, color = "black"), panel.border = element_blank(), axis.line =element_line(colour = "black", size = 1, linetype = 1, lineend = "butt"), axis.ticks = element_line(size=1), axis.ticks.length = unit(5, "pt"),  panel.grid=element_blank(),  axis.title.y = element_text(color="black" , size = 16), axis.text.y = element_text( color = "black", size = 14, angle=0))
dev.off()


#Save data as table
write.table(comb.inter.data.braypart, file="Inter_Beta_diversity_cecum_wuniF.txt", sep="\t")

#################################################################
#################################################################
#################################################################

save.image(file="mouse.RData")
load(file="mouse.RData")

#################################################################
#################################################################
#################################################################
