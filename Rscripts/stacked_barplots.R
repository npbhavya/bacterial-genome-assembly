#install.packages("funrar")
library(ggplot2)
library(reshape2)
library(funrar)
library(vegan)

category=read.csv("~/OneDrive - Flinders/EDgenomes/genotypic_analysis/PATRIC_Superclass.csv", row.names = 1, )
category_t=t(category)
category_rel_mat = make_relative(as.matrix(category_t))
category_long = melt(category_rel_mat, id = c("Row.Labels"))

ggplot(category_long, aes(x = Var1, fill = Var2, y = value)) + 
  geom_bar(stat = "identity", colour = "black") + 
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 6, face = "bold", colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative Abundance (%)", fill = "category")

## stats
ecology=as.matrix(c("kelp", "kelp", "kelp", "kelp", "water", "kelp", "water", "water", "kelp"))
dist.category= vegdist(category_rel_mat, method = "bray")
adonis2(dist.category~ecology, data=category, permutations=9999)

######### subcategory
subcategory=read.csv("~/OneDrive - Flinders/EDgenomes/genotypic_analysis/PATRIC_class.csv", row.names = 1, header = TRUE)
subcategory_t=t(subcategory)
subcategory_rel_mat = make_relative(as.matrix(subcategory_t))
subcategory_long = melt(subcategory_rel_mat, id = c("Row.Labels"))
ggplot(subcategory_long, aes(x = Var1, fill = Var2, y = value)) + 
  geom_bar(stat = "identity", colour = "black") + 
  theme(axis.text.x = element_text(angle = 90, size = 4, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold", colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative Abundance (%)", fill = "genomes")

## stats
ecology=as.matrix(c("kelp", "kelp", "kelp", "kelp", "water", "kelp", "water", "water", "kelp"))
dist.subcategory= vegdist(subcategory_rel_mat, method = "bray", na.rm = TRUE)
adonis2(dist.subcategory~ecology, data=as.data.frame(subcategory_t), permutations=9999)


subsystem=read.csv("~/OneDrive - Flinders/EDgenomes/genotypic_analysis/PATRIC_subclass.csv", row.names = 1, header = TRUE)
subsystem_t=t(subsystem)
subsystem_rel_mat = make_relative(as.matrix(subsystem_t))
subsystem_long = melt(subsystem_rel_mat, id = c("Row.Labels"))
ggplot(subsystem_long, aes(x = Var1, fill = Var2, y = value)) + 
  geom_bar(stat = "identity", colour = "black") + 
  theme(axis.text.x = element_text(angle = 90, size = 4, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 16), 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold", colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "", y = "Relative Abundance (%)", fill = "genomes")+
  theme(legend.position=c(.5, .5))
