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

phy<-qza_to_phyloseq(features="filtered.merge.table.two.qza", tree="rooted-tree.qza", taxonomy="taxonomy.two.qza",metadata="mouse.metadata.txt")

random_tree = rtree(ntaxa(phy), rooted=TRUE, tip.label=taxa_names(phy))
plot(random_tree)

#Give a colnames to separate different taxonomic levels
colnames(tax_table(phy))=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "OTU")

rownames(sample_data(phy))
colnames(sample_data(phy))

# Remove taxa with 0 abundance
phy = subset_taxa(phy, rowSums(otu_table(phy)) != 0)

##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {
    x/sum(x)
}
phy.relative.table = transformSampleCounts(phy, normalizeSample)

rownames(sample_data(phy.relative.table))
colnames(sample_data(phy.relative.table))

Phylum.rel.table = tax_glom(phy.relative.table, taxrank = "Phylum")
Class.rel.table = tax_glom(phy.relative.table, taxrank = "Class")
Order.rel.table = tax_glom(phy.relative.table, taxrank = "Order")
Family.rel.table = tax_glom(phy.relative.table, taxrank = "Family")
Genus.rel.table = tax_glom(phy.relative.table, taxrank = "Genus")
OTU.rel.table = tax_glom(phy.relative.table, taxrank = "OTU")

#################################################################
#################################################################
#################################################################

save.image(file="mouse.upstream.RData")
load(file="mouse.upstream.RData")

#################################################################
#################################################################
#################################################################