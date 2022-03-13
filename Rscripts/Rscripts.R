# Vibrio genome analysis

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#install.packages("RColorBrewer")
#install.packages("geosphere")

suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
# Load dependencies
library("ggplot2")
library("ggdendro")
library("reshape2")
library("grid")
library(circlize)
library(vegan)
library(geosphere)

#######################################################################
##### Genotypic analysis
#######################################################################
#most basic heatmap with genus level
genusdata = read.csv("~/OneDrive - Flinders/EDgenomes/genotypic_analysis/EdgenomesOnly-ANI.csv", row.names = 1)  # load the data w row names 
genus = as.matrix(genusdata)
Heatmap(genus,row_names_gp=gpar(cex= .6),row_dend_width = unit(0, "cm"),column_dend_height = unit(4,"cm"), cluster_columns = FALSE, cluster_rows = FALSE )

#Clustering 
genus.dendro <- as.dendrogram(hclust(d = dist(x = genus)))
# Create dendro
dendro.plot <- ggdendrogram(data = genus.dendro, rotate = FALSE)
# Preview the plot
print(dendro.plot)


### comparing genotype data ANI versus functions clustering 
functions = read.csv("OneDrive - Flinders/EDgenomes/Edgenomes_PATRIC_gene_presence_absence.csv", row.names = 1)  # load the data w row names 
functions = as.matrix(t(functions))
Heatmap(functions,row_names_gp=gpar(cex= .6),row_dend_width = unit(0, "cm"),column_dend_height = unit(4,"cm"), cluster_columns = FALSE, cluster_rows = FALSE )

#Clustering 
functions.dendro <- as.dendrogram(hclust(d = dist(x = functions)))
# Create dendro
dendro.plot <- ggdendrogram(data = functions.dendro, rotate = FALSE)
# Preview the plot
print(dendro.plot)

#######################################################################
### Phenotypic analysis
#######################################################################
pheno_data = read.csv("Downloads/EDgenomes/phenotypic/phenotypic_data_filtered.csv", row.names = 1)  # load the data w row names 
phenotype = as.matrix(pheno_data)
#Clustering 
pheno.dendro <- as.dendrogram(hclust(d = dist(x = phenotype)))
dendro.plot <- ggdendrogram(data = pheno.dendro, rotate = FALSE)
print(dendro.plot)

Heatmap(phenotype,row_names_gp=gpar(cex= 1.4),column_names_gp =gpar(cex= .7),row_dend_width = unit(20, "mm"),
        heatmap_legend_param = list(title = "Growth level"),
        column_dend_height = unit(10,"cm"), cluster_columns=FALSE)

#heatmap for Carbon sources
carbon=read.csv("Downloads/EDgenomes/phenotypic/phenotypic_carbon_data.csv", row.names = 1)
matc= as.matrix(carbon)
Heatmap(matc)

#heatmap for nitrogen sources
nitrogen=read.csv("Downloads/EDgenomes/phenotypic/phenotypic_nitrogen_data.csv", row.names = 1)
matn=as.matrix(nitrogen)
Heatmap(matn)

#heatmap for phosphorous sources
phosphorous=read.csv("Downloads/EDgenomes/phenotypic/phenotypic_phophorous_data.csv", row.names = 1)
matp=as.matrix(phosphorous)
Heatmap(matp)

#######################################################################
## Corelating phenotype to genotype data
#######################################################################
dist.geno= vegdist(genus, method = "bray")
dist.pheno= vegdist(phenotype, method = "bray")

mantel_geno=mantel(dist.geno, dist.pheno, method = "spearman", permutations=9999)
mantel_geno

mantel_pheno=mantel(dist.pheno, dist.geno, method ="pearson", permutations =9999)
mantel_pheno

#######################################################################
## PERMANOVA
#######################################################################
#1 for kelp, 0 to water
#ecology=as.matrix(c("1", "1", "0", "1", "1", "0", "1", "0", "1"))
ecology=as.matrix(c("kelp", "kelp", "water", "kelp", "kelp", "water", "kelp", "water", "kelp"))
genotype=as.data.frame(genus)
adonis2(dist.geno~ecology, data=genotype, permutations=9999)
#no significance when comparing genotype data to location - water/kelp

phenotype=as.data.frame(phenotype)
adonis2(dist.pheno~ecology, data=phenotype, permutations=9999)
#no significance when comparing phenotype data to location - water/kelp
