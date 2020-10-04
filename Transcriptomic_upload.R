### RNA.Seq.mouse
### Benjamin Wu
### Adapted from Imran Sulaiman's code

file = read.table('~/RNA.Seq.Data/Mouse/Mouse_RNA_Count_Table.txt')
RNA.Seq <- file

update.packages()

#Load Packages
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
library(vegan)
library(devtools)
library("vsn")
library("hexbin")
library("SummarizedExperiment")


#install.packages('')
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("Glimma")
BiocManager::install("phyloseq")
BiocManager::install("SpiecEasi")
BiocManager::install("KEGGREST")
BiocManager::install("pathview")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("vsn")
BiocManager::install("hexbin")
BiocManager::install("SummarizedExperiment")
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
install.packages('omu')
install.packages('curl')
install.packages('maptools')
install.packages('devtools')
install_github("zdk123/SpiecEasi")
install.packages('vegan')

#Load and Save File
getwd()
setwd('~/Dropbox/Mouse')

load(file="Mouse.RNA.seq.Analysis.Ben.RData")
save.image(file="Mouse.RNA.seq.Analysis.Ben.RData")

data.frame(file('/Users/benjaminwu/Dropbox/Mouse/Map.merge.txt'))
class(dds)
merge<- read.table("/Users/benjaminwu/Dropbox/Mouse/Map.merge.txt", header=TRUE, sep="\t", row.names="Sample.ID")
class(merge)
merge.dataframe <- data.frame(merge)
class(merge.dataframe)

#Count table --> assay(dds)
#Meta Data --> colData(dds)
#normalized Data --> rld
assay(rld)
assay(dds)
dds.deseq.results <- DESeq(dds)

colData(dds)
colData(rld)

### regularlized log change (based on negative binomial distribut)
###
write.table(assay(rld),file="rld.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
### 

ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

#Subset [x.y] remember these commands 
rld.moc <- rld[,rld$Condition=="MOC"]

#Subset by day 
rld.day01 <- rld[,rld$Day=="1"]
rld.day14 <- rld[,rld$Day=="14"]
colData(rld.day01)
colData(rld.day14)

#Subset by MyD88KO exposure 
rld.day14.C57BL6  <- rld.day14[,rld.day14$Mouse.genetics=="C57BL6"]
rld.day14.MYD88KO <- rld.day14[,rld.day14$Mouse.genetics=="MyD88KO"]
colData(rld.day14.C57BL6)
colData(rld.day14.MYD88KO)

#Subset by C57BL6 by single 
rld.day14.C57BL6.single  <- rld.day14.C57BL6[,rld.day14.C57BL6$Single.Multiple.Exposure=="Single"]
colData(rld.day14.C57BL6.single)

#Datasets to compare 
colData(rld.day01)
colData(rld.day14.C57BL6.single)

# Need to subset te data for:
# PBS D1 D14 

rld.PBS.D1.D14.D14multiple <- rld[ , rld$Condition.Day %in% c("PBS" , "MOC.14", "MOC.1") ]
​plotPCA(rld.PBS.D1.D14.D14multiple, "Condition.Day")

# PBS D1 D14 
rld.day14.day1.PBS <- rld[ , rld$combined %in% c("PBS.14.C57BL6.Single", "MOC.14.C57BL6.Single" , "PBS.1.C57BL6.Single", "MOC.1.C57BL6.Single") ]
rld.day14.PBS <- rld[ , rld$combined %in% c("PBS.14.C57BL6.Single", "MOC.14.C57BL6.Single") ]
rld.day14.all.PBS <- rld[ , rld$combined %in% c("PBS.14.C57BL6.Single", "PBS.1.C57BL6.Single","MOC.14.C57BL6.Single") ]

rld.PBS.D1.D14 <- rld[ , rld$Condition.Day %in% c("PBS", "MOC.14", "MOC.1") ]
​plotPCA(rld.PBS.D1.D14, "Condition.Day")

rld.PBS.D14 <- rld[ , rld$Condition.Day %in% c("PBS" , "MOC.14") ]
​plotPCA(rld.PBS.D14 , "Condition.Day")

​plotPCA(rld.day14.day1.PBS, "Condition.Day")
​plotPCA(rld.day14.PBS, "Condition.Day")
​plotPCA(rld.day14.all.PBS, "Condition.Day")

# This is logfold DESeq2 analysis 
dds

#Subset by day 
dds.day01 <- dds[, dds $Day=="1"]
dds.day14 <- dds[, dds $Day=="14"]
colData(dds.day01)
colData(dds.day14)

#Subset by C57BL6 by multiple or single 
dds.day14.C57BL6.single  <- dds.day14.C57BL6[,dds.day14.C57BL6$Single.Multiple.Exposure=="Single"]
colData(dds.day14.C57BL6.single)

#Subset by Day 14, Day 1, and PBS in various combinations 
dds.day14.day1.PBS <- dds[ , dds$combined %in% c("PBS.14.C57BL6.Single", "MOC.14.C57BL6.Single" , "PBS.1.C57BL6.Single", "MOC.1.C57BL6.Single") ]
dds.day14.PBS <- dds[ , dds$combined %in% c("PBS.14.C57BL6.Single", "MOC.14.C57BL6.Single") ]
dds.day14.all.PBS <- dds[ , dds$combined %in% c("PBS.14.C57BL6.Single", "PBS.1.C57BL6.Single","MOC.14.C57BL6.Single") ]
dds.day1.all.PBS <- dds[ , dds$combined %in% c("PBS.14.C57BL6.Single", "PBS.1.C57BL6.Single","MOC.1.C57BL6.Single") ]

#Datasets to compare 
colData(dds.day01)
colData(dds.day14.C57BL6.single)

#Differential Analysis of PBS versus MOC for Day 1
#Drop any uwanted levels
dds.day01$Condition <- droplevels(dds.day01$Condition)

#Set the Baseline Comparator for the Differntial Analysis
dds.day01$Condition <- relevel(dds.day01$Condition, ref ="PBS")

#Differential Analysis of PBS versus MOC for Day 1
#Drop any uwanted levels
dds.day14.C57BL6.single$Condition <- droplevels(dds.day14.C57BL6.single$Condition)

#Set the Baseline Comparator for the Differntial Analysis
dds.day14.C57BL6.single$Condition <- relevel(dds.day14.C57BL6.single$Condition, ref ="PBS")

#Differential Analysis of PBS versus MOC for Day 1
#Drop any uwanted levels
dds.day14.C57BL6.multiple$Condition <- droplevels(dds.day14.C57BL6.multiple$Condition)

#Set the Baseline Comparator for the Differntial Analysis
dds.day14.C57BL6.multiple$Condition <- relevel(dds.day14.C57BL6.multiple$Condition, ref ="PBS")

#Differential Analysis of PBS versus MOC for Day 1
#Drop any uwanted levels
dds.day14.MOC.single.multiple$Condition <- droplevels(dds.day14.MOC.single.multiple$Condition)

#Set the Baseline Comparator for the Differntial Analysis
dds.day14.MOC.single.multiple$Condition <- relevel(dds.day14.MOC.single.multiple$Condition, ref ="Single")

#Differential Analysis of PBS versus MOC for Day 1
#Drop any uwanted levels
dds.day14.day28.MOC.null.condition$Condition <- as.factor(dds.day14.day28.MOC.null.condition$Condition)

dds.day14.day28.MOC.null.condition$Condition <- droplevels(dds.day14.day28.MOC.null.condition$Condition)

#Set the Baseline Comparator for the Differntial Analysis
dds.day14.day28.MOC.null.condition$Condition <- relevel(dds.day14.day28.MOC.null.condition$Condition, ref ="28")


dds.day14.all.PBS$Condition <- as.factor(dds.day14.all.PBS$Condition)
dds.day14.all.PBS$Condition <- droplevels(dds.day14.all.PBS$Condition)

#Set the Baseline Comparator for the Differntial Analysis
dds.day14.all.PBS$Condition <- relevel(dds.day14.all.PBS$Condition, ref ="PBS")

dds.day1.all.PBS$Condition <- as.factor(dds.day1.all.PBS$Condition)
dds.day1.all.PBS$Condition <- droplevels(dds.day1.all.PBS$Condition)

#Set the Baseline Comparator for the Differntial Analysis
dds.day1.all.PBS$Condition <- relevel(dds.day1.all.PBS$Condition, ref ="PBS")

#Datasets to compare 
colData(dds.day01)
colData(dds.day14.C57BL6.single)

#Differential Analysis of PBS versus MOC for Day 1
#Run the differential Analsysis, so what is upregulated will be upregulated in MOC
dds.day01.deseq <- DESeq(dds.day01)
#Output results of differential analysis into a table
dds.day01.deseq.res <- results(dds.day01.deseq)
dds.day01.deseq.res

dds.day14.C57BL6.single.deseq <- DESeq(dds.day14.C57BL6.single)
dds.day14.C57BL6.single.deseq.res <- results(dds.day14.C57BL6.single.deseq, contrast=c("Condition","PBS", "MOC"))
dds.day14.C57BL6.single.deseq.res

dds.day01.deseq.res
dds.day14.C57BL6.single.deseq.res

write.table(dds.day01.deseq.res,file="dds.day01.deseq.res.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)
write.table(dds.day14.C57BL6.single.deseq.res,file="dds.day14.C57BL6.single.deseq.res.IPA.txt", sep="\t", col.names = NA, row.names = TRUE, quote=FALSE)

?results
?DESeq

#=========================================================
/////////////////////////PCOA PLOT///////////////////////
#=========================================================

​plotPCA(rld.day14.day1.PBS, "Condition.Day")

#colData(rld.PBS.D1.D14))
​plotPCA(rld.day14.day1.PBS, "Condition.Day")

vegdist   = vegdist(t(assay(rld.day14.day1.PBS)), method="bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =6)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
require(data.table)
newResults <- merge(x = CmdScale, y = colData(rld.day14.day1.PBS), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~Condition.Day,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Condition.Day",suffixes=c("",".centroid"))

pdf("Mouse.RNA.seq_PBS.D1.D14.new.Bray.longer.pdf", height = 5, width = 7)
    ggplot(newResults, aes(PC1, PC2, color=Condition.Day)) +
    geom_point(size=5,alpha=0.5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    #coord_fixed() +
    scale_color_manual(values=c("#FFA502", "#ff8c00", "#006400")) +
    #plot ellipse
    #stat_ellipse(type = "t") +
    #plot point and lines from centroid
    geom_point(data=centroids, aes(x=PC1, y=PC2, color=Condition.Day), size=0) +
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Condition.Day))+
    #geom_text_repel(aes(label=ifelse(newResults$name%in% c("COPD.0002.BAL.L.untouched","COPD.0030.BAL.L.untouched","COPD.0035.BAL.L.untouched") , as.character(newResults$name),'')),size=3,force=25) +
    #labels centroids
    geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("MOC.D1", "MOC.D14", "PBS")), size=10) +
    #geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=lab), parse=TRUE,size=10) +
    #scale_x_reverse() +
    theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

adonis(vegdist ~ rld.day14.day1.PBS$Condition.Day)

Call:
adonis(formula = vegdist ~ rld.day14.day1.PBS$Condition.Day) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                                 Df SumsOfSqs    MeanSqs F.Model   R2 Pr(>F)   
rld.day14.day1.PBS$Condition.Day  2 0.0012607 0.00063036   6.245 0.51  0.003 **
Residuals                        12 0.0012113 0.00010094         0.49          
Total                            14 0.0024720                    1.00          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#=========================================================
////////////////////VOLCANO PLOT///////////////////////
#=========================================================

alpha=0.05
dds.day14.all.PBS.deseq.2.res 

#Compute FDR in a log scale
dds.day14.all.PBS.deseq.2.res$sig <- -log10(dds.day14.all.PBS.deseq.2.res$padj)

#See How many are now infinte
sum(is.infinite(dds.day14.all.PBS.deseq.2.res$sig <- -log10(dds.day14.all.PBS.deseq.2.res$padj)))

#Convert any Infinite values to a maximum value for adjusted pvalue for the graph (here I choose 350)
dds.day14.all.PBS.deseq.2.res[is.infinite(dds.day14.all.PBS.deseq.2.res$sig),"sig"] <- 350 

# Create variable for color gradient of the dots in the Volcano Plot
cols <- densCols(dds.day14.all.PBS.deseq.2.res$log2FoldChange, dds.day14.all.PBS.deseq.2.res$sig)
cols[dds.day14.all.PBS.deseq.2.res$pvalue ==0] <- "purple"
cols[dds.day14.all.PBS.deseq.2.res$log2FoldChange > 1 & dds.day14.all.PBS.deseq.2.res$padj < alpha ] <- "darkorange"
cols[dds.day14.all.PBS.deseq.2.res$log2FoldChange < -1 & dds.day14.all.PBS.deseq.2.res$padj < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
dds.day14.all.PBS.deseq.2.res$pch <- 1
dds.day14.all.PBS.deseq.2.res$pch[dds.day14.all.PBS.deseq.2.res$pvalue ==0] <- 6

dds.day14.all.PBS.deseq.2.res
head(dds.day14.all.PBS.deseq.2.res)
rownames(dds.day14.all.PBS.deseq.2.res)

pdf(file="D14.FDR.1.all.PBS.pdf", width=10, height=10)
    ggplot(data.frame(dds.day14.all.PBS.deseq.2.res), aes(x = log2FoldChange, y = sig,label=rownames(dds.day14.all.PBS.deseq.2.res))) +
    geom_point(color=cols, size = 2) + #Chose Colors and size for dots
    geom_text(aes(label=ifelse(dds.day14.all.PBS.deseq.2.res$log2FoldChange>=3.75 & dds.day14.all.PBS.deseq.2.res$padj < 0.05 , as.character(rownames(dds.day14.all.PBS.deseq.2.res)), ifelse(dds.day14.all.PBS.deseq.2.res$log2FoldChange<=-3.75 & dds.day14.all.PBS.deseq.2.res$padj < 0.05, as.character(rownames(dds.day14.all.PBS.deseq.2.res)),''))),size=2,force=8, nudge_x = .2, nudge_y =.2) + #Label values based on parameters, including pcal and logFC
    theme(legend.position = "none") + #Remove Legend
    geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
    xlab("Effect size: log2(fold-change)") + #label X Axis
    ylab("-log10(adjusted p-value)") + #label Y Axis
    theme #Set Theme
dev.off()

#=========================================================
////////////////////VOLCANO PLOT///////////////////////
#=========================================================

dds.day1.all.PBS.deseq.2.res

alpha = 0.05

#Compute FDR in a log scale
dds.day1.all.PBS.deseq.2.res$sig <- -log10(dds.day1.all.PBS.deseq.2.res$padj)

#See How many are now infinte
sum(is.infinite(dds.day1.all.PBS.deseq.2.res$sig <- -log10(dds.day1.all.PBS.deseq.2.res$padj)))

#Convert any Infinite values to a maximum value for adjusted pvalue for the graph (here I choose 350)
dds.day1.all.PBS.deseq.2.res[is.infinite(dds.day1.all.PBS.deseq.2.res$sig),"sig"] <- 350 ### error here  object 'res.bal' not found

# Create variable for color gradient of the dots in the Volcano Plot
cols <- densCols(dds.day1.all.PBS.deseq.2.res$log2FoldChange, dds.day1.all.PBS.deseq.2.res$sig)
cols[dds.day1.all.PBS.deseq.2.res$pvalue ==0] <- "purple"
cols[dds.day1.all.PBS.deseq.2.res$log2FoldChange > 1 & dds.day1.all.PBS.deseq.2.res$padj < alpha ] <- "darkorange"
cols[dds.day1.all.PBS.deseq.2.res$log2FoldChange < -1 & dds.day1.all.PBS.deseq.2.res$padj < alpha ] <- "darkgreen"

# Create a Variable for the size of the dots in the Volcano Plot
dds.day1.all.PBS.deseq.2.res$pch <- 1
dds.day1.all.PBS.deseq.2.res$pch[dds.day1.all.PBS.deseq.2.res$pvalue ==0] <- 6

dds.day1.all.PBS.deseq.2.res
head(dds.day1.all.PBS.deseq.2.res)
rownames(dds.day1.all.PBS.deseq.2.res)


pdf(file="D01.FDR.1.all.PBS.pdf", width=10, height=10)
    ggplot(data.frame(dds.day1.all.PBS.deseq.2.res), aes(x = log2FoldChange, y = sig,label=rownames(dds.day1.all.PBS.deseq.2.res))) +
    geom_point(color=cols, size = 2) + #Chose Colors and size for dots
    geom_text(aes(label=ifelse(dds.day1.all.PBS.deseq.2.res$log2FoldChange>=5 & dds.day1.all.PBS.deseq.2.res$padj < 0.05 , as.character(rownames(dds.day1.all.PBS.deseq.2.res)), ifelse(dds.day1.all.PBS.deseq.2.res$log2FoldChange<=-5 & dds.day1.all.PBS.deseq.2.res$padj < 0.05, as.character(rownames(dds.day1.all.PBS.deseq.2.res)),''))),size=2,force=8, nudge_x = .4, nudge_y =.4) + #Label values based on parameters, including pcal and logFC
    theme(legend.position = "none") + #Remove Legend
    geom_hline(yintercept=-log10(alpha), color="red",linetype="dashed") + #Create Reference line for FDR
    xlab("Effect size: log2(fold-change)") + #label X Axis
    ylab("-log10(adjusted p-value)") + #label Y Axis
    theme #Set Theme
dev.off()
