### This R script has been adapted from code written by Christina Pavloudi ###
########## Use this space to add contact information for Christina ##########



#Tell me where is the working directory (this command its like pwd in UNIX)
getwd()
#Change the working directory (like cd)
setwd("/Volumes/Drive/DADA2phyloseq_output/outputFiles_final")

# for the actual analysis
library(phyloseq); packageVersion("phyloseq")
# for the beautiful plots
library(ggplot2); packageVersion("ggplot2")
#to create beautiful colour palettes for your plots
library(RColorBrewer); packageVersion("RColorBrewer") 
# Define a default theme for ggplot graphics
theme_set(theme_bw()) 

library(dplyr)
library(ape)
library(vegan)

####################### Import and format your data ##########################
##############################################################################

#import the OTU table (or else biotic data)
biotic <- read.csv("seq_table2.csv", sep = ",", header=TRUE, row.names = 1)
#import the taxonomy table
taxonomy <- read.csv("seq_Taxonomy_silva.csv", sep = ",", header=FALSE, row.names = 1, na.strings=c("","NA"))
colnames(taxonomy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#check where there are NA values in the taxonomy data frame
colSums(is.na(taxonomy))

#fill in Phylum column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Phylum[i])){
      taxonomy$Phylum[i]=taxonomy$Kingdom[i]
    }
  if (sum(is.na(taxonomy$Phylum))==0) {
    break
  }
}

#fill in Class column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Class[i])){
      taxonomy$Class[i]=taxonomy$Phylum[i]
    }
  if (sum(is.na(taxonomy$Class))==0) {
    break
  }
}

#fill in Order column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Order[i])){
      taxonomy$Order[i]=taxonomy$Class[i]
    }
  if (sum(is.na(taxonomy$Order))==0) {
    break
  }
}

#fill in Family column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Family[i])){
      taxonomy$Family[i]=taxonomy$Order[i]
    }
  if (sum(is.na(taxonomy$Family))==0) {
    break
  }
}

#fill in Genus column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Genus[i])){
      taxonomy$Genus[i]=taxonomy$Family[i]
    }
  if (sum(is.na(taxonomy$Genus))==0) {
    break
  }
}

#fill in Species column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Species[i])){
      taxonomy$Species[i]=taxonomy$Genus[i]
    }
  if (sum(is.na(taxonomy$Species))==0) {
    break
  }
}

#check where there are NA values in the taxonomy data frame
colSums(is.na(taxonomy))

#convert the taxonomy data from data frame to matrix
taxonomy_matrix <- as.matrix(taxonomy)

# prepare the object for the phyloseq object
TAX = tax_table(taxonomy_matrix)
head(TAX)

#convert the biotic data from data frame to matrix
biotic_matrix <- as.matrix(biotic)
# prepare the objects for the phyloseq object
OTU = otu_table(biotic_matrix, taxa_are_rows = TRUE)
head(OTU)

#import the metadata of the samples
metadata_physeq <- read.csv("metadata2.csv", header=TRUE, row.names = 1) 
# prepare the objects for the phyloseq object
META = sample_data(metadata_physeq)
head(META)
#class(META)

###########################PHYLOSEQ analysis##################################
##############################################################################

# combine them all to create the phyloseq object
physeq = phyloseq(OTU, TAX, META)
#check what the phyloseq object contains
physeq

#Check if there are ASVs with no counts and how many there are
any(taxa_sums(physeq) == 0)
sum(taxa_sums(physeq) == 0)
#Remove ASVs/OTUs/taxa that are empty, i.e. that have no counts
physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
# Remove samples that are now empty, i.e. that have no counts
physeq <- prune_samples(sample_sums(physeq) > 0, physeq)

# create a bar plot, for the taxon rank you like, e.g Phylum 
# this bar chart is created with the default colour palette of phyloseq
barchart <- plot_bar(physeq, fill = "Phylum")
pdf("barchart_silva.pdf", width = 12)
print(barchart)
dev.off()

#get the data frame from the phyloseq object
pd <- psmelt(physeq)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd[,c("Phylum")])))
HowManyFamilies <- length(unique(unlist(pd[,c("Family")])))
HowManyGenera <- length(unique(unlist(pd[,c("Genus")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Set2".
# "Set2" is one of the colorblind friendly palettes
# Another example of a colorblind friendly palette is "Dark2"
# If you want, by running the command display.brewer.all(colorblindFriendly = TRUE)
# you can see all the colorblind friendly palettes of the RColorBrewer package.
## Now the "getPalette" variable, we set it in the if statement of the marker gene.
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
PhylaPalette = getPalette(HowManyPhyla)
FamilyPalette = getPalette(HowManyFamilies)
GenusPalette = getPalette(HowManyGenera)

#and do the actual plotting
barchart_palette_better <- ggplot(pd, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~WorkerType) +
  facet_wrap(~WorkerType, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5), axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16,face="bold"), 
        legend.title=element_text(size=16, face="bold"), legend.text=element_text(size=14), strip.text.x = element_text(size = 24, face = "bold"))
ggsave("barchart_palette_better_silva_custom.png", width = 22, height = 15, dpi=300)

## by colony
barchart_palette_better <- ggplot(pd, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  facet_wrap(~Colony, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), 
        legend.title=element_text(size=14), legend.text=element_text(size=12))
ggsave("barchart_palette_better_silva_colony.png", width = 22, height = 15, dpi=300)

## colony, family level
barchart_palette_better <- ggplot(pd, aes(x = Sample, y = Abundance, factor(Family), fill = factor(Family))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = FamilyPalette) +
  labs(fill = "Family") +
  facet_wrap(~Colony, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), 
        legend.title=element_text(size=14), legend.text=element_text(size=12))
ggsave("barchart_silva_colony_family.png", width = 22, height = 25, dpi=300)

barchart_palette_better <- ggplot(pd, aes(x = Sample, y = Abundance, factor(Family), fill = factor(Family))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = FamilyPalette) +
  labs(fill = "Family") +
  facet_wrap(~WorkerType, scales="free_x") +
  guides(fill=guide_legend(ncol=3)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), 
        legend.title=element_text(size=14), legend.text=element_text(size=12))
ggsave("barchart_silva_worker_family.png", width = 22, height = 15, dpi=300)

## colony, family level
barchart_palette_better <- ggplot(pd, aes(x = Sample, y = Abundance, factor(Family), fill = factor(Family))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = FamilyPalette) +
  labs(fill = "Family") +
  facet_wrap(~Colony, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), 
        legend.title=element_text(size=14), legend.text=element_text(size=12))
ggsave("barchart_silva_colony_families.png", width = 22, height = 15, dpi=300)

## worker, genus
barchart_palette_better <- ggplot(pd, aes(x = Sample, y = Abundance, factor(Genus), fill = factor(Genus))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = GenusPalette) +
  labs(fill = "Genus") +
  facet_wrap(~WorkerType, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), 
        legend.title=element_text(size=14), legend.text=element_text(size=12))
ggsave("barchart_silva_worker_genera.png", width = 22, height = 15, dpi=300)

##fecal grouping
barchart_palette_better <- ggplot(pd, aes(x = Sample, y = Abundance, factor(Family), fill = factor(Family))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = FamilyPalette) +
  labs(fill = "Family") +
  facet_wrap(~Colony, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0.5), axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), 
        legend.title=element_text(size=14), legend.text=element_text(size=12))
ggsave("barchart_silva_colony_families.png", width = 22, height = 15, dpi=300)

# plot the diversity indices with colour coding by e.g. dpw (info included in the metadata) 
richness <- plot_richness(physeq, measures=c("Observed", "Chao1", "ACE"), color="Colony")
pdf("richness_silva.pdf", width = 14)
print(richness)
dev.off()

richness1 <- plot_richness(physeq, measures=c("Shannon", "Simpson", "Chao1", "ACE"), color="Colony")
pdf("richness1_silva.pdf", width = 14)
print(richness1)
dev.off()

richness2 <- plot_richness(physeq, measures=c("Observed", "Chao1", "ACE"), color="WorkerType", shape="Colony")
pdf("richness2_silva.pdf", width = 14)
print(richness2)
dev.off()

richness3 <- plot_richness(physeq, measures=c("Observed", "Chao1", "ACE"), color="IDS")
pdf("richness3_silva.pdf", width = 14)
print(richness3)
dev.off()

richness4 <- plot_richness(physeq, measures=c("Observed", "Chao1", "ACE"), color="FoS")
pdf("richness4_silva.pdf", width = 14)
print(richness4)
dev.off()

richness5 <- plot_richness(physeq, measures=c("Observed", "Chao1", "ACE"), color="FeS")
pdf("richness5_silva.pdf", width = 14)
print(richness5)
dev.off()

richness5a <- plot_richness(physeq, measures=c("Shannon", "Simpson", "Chao1", "ACE"), color="FeS")
pdf("richness5a_silva.pdf", width = 14)
print(richness5a)
dev.off()

##Kruskal-Wallis
library(ggpubr)

richness <- estimate_richness(physeq1)
head(richness)

sample_data(physeq1)$FeS <- factor((sample_data(physeq1)$FeS), levels=c("High","Intermediate","Low"))

a_my_comparisons <- list( c("High", "Intermediate", "Low"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

test_kruskal <- plot_richness(physeq1, x="FeS", measures=c("Shannon", "Chao1", "ACE"), color = "FeS") + geom_boxplot(alpha=0.6) + theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12)) + stat_compare_means(test = "kruskal.test", comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args)

pdf("ks.pdf", width = 14)
print(test_kruskal)
dev.off()

kruskal.test(richness$Chao1 ~ sample_data(physeq1)$FeS)
kruskal.test(richness$ACE ~ sample_data(physeq1)$FeS)
kruskal.test(richness$Shannon ~ sample_data(physeq1)$FeS)

kruskal.test(richness$Chao1 ~ sample_data(physeq1)$WorkerType)
kruskal.test(richness$ACE ~ sample_data(physeq1)$FeS)
kruskal.test(richness$Shannon ~ sample_data(physeq1)$FeS)

#### NMDS plots ####

# create the nMDS plot colour coded by e.g. dpw (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq, ord.nmds.bray, color="Colony", title="Bray NMDS")
pdf("p1_silva.pdf")
print(p1)
dev.off()

ord.nmds.bray <- ordinate(physeq, method="NMDS", distance="bray")
p2 <- plot_ordination(physeq, ord.nmds.bray, color="WorkerType", title="Bray NMDS")
pdf("p2_silva.pdf")
print(p2)
dev.off()

ord.nmds.bray <- ordinate(physeq, method="NMDS", distance="bray")
p3 <- plot_ordination(physeq, ord.nmds.bray, color="IDS", title="Bray NMDS")
pdf("p3_silva.pdf")
print(p3)
dev.off()

ord.nmds.bray <- ordinate(physeq, method="NMDS", distance="bray")
p4 <- plot_ordination(physeq, ord.nmds.bray, color="FoS", shape="WorkerType", title="Bray NMDS")
pdf("p4_silva.pdf")
print(p4)
dev.off()

ord.nmds.bray <- ordinate(physeq, method="NMDS", distance="bray")
p5_ellipse <- plot_ordination(physeq, ord.nmds.bray, color="FeS", shape="WorkerType", title="Bray NMDS") + stat_ellipse(type = "norm", linetype=2)
pdf("p5_ellipse.pdf")
print(p5_ellipse)
dev.off()

##do more PCoA
ord.pcoa.bray <- ordinate(physeq, method="PCoA", distance="bray")
p5a <- plot_ordination(physeq, ord.pcoa.bray, color="FeS", shape="WorkerType", title="Bray PCoA")
pdf("p5_pcoa.pdf")
print(p5a)
dev.off()


#check if the clustering you see in the nMDS plot is 
#statistically significant by running PERMANOVA
library("vegan")
metadata_permanova <- as(sample_data(physeq), "data.frame")
permanova.colony <- adonis2(distance(physeq, method="bray") ~ Colony, data = metadata_permanova)
permanova.colony

permanova.worker <- adonis2(distance(physeq, method="bray") ~ WorkerType, data = metadata_permanova)
permanova.worker

permanova.intActivity <- adonis2(distance(physeq, method="bray") ~ IDS, data = metadata_permanova)
permanova.intActivity

permanova.extActivity <- adonis2(distance(physeq, method="bray") ~ FoS, data = metadata_permanova)
permanova.extActivity

permanova.fecal <- adonis2(distance(physeq, method="bray") ~ FeS, data = metadata_permanova)
permanova.fecal


###########################EXPORT ASV ABUN TABLE TO EXCEL#####################
##############################################################################

write.csv(pd, "abundance.csv")

#merges ASVs that have the same taxonomy at a the Phylum level
physeq_merged_Phylum <- tax_glom(physeq, "Phylum")
#transform counts to percentages
ps0 <- transform_sample_counts(physeq_merged_Phylum, function(x) x / sum(x))
plot_bar(ps0, fill="Phylum")
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "ASV")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "ASV")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "ASV")
#save the final phyto otu table with the taxonomies
write.csv(OTU_TAX_merged, "OTU_TAX_merged_phylum.csv")

############## PHYLOSEQ analysis for specific taxa ############################
###############################################################################

#### Network analysis ####

ig <- make_network(physeq, dist.fun="jaccard", max.dist=0.5)
networkPlot <- plot_network(ig, physeq, color="FeS", shape="WorkerType", line_weight=0.4, label=NULL)
pdf("network.pdf")
print(networkPlot)
dev.off()

networkPlot2 <- plot_net(physeq, maxdist = 0.5, color = "FeS", shape="WorkerType")
pdf("network2.pdf")
print(networkPlot2)
dev.off()


#### tree building and unifrac #### 

tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(tree)

physeq1 = merge_phyloseq(physeq, tree)
physeq1

# Create subset of the phyloseq object, e.g. 
physeq.Actinobacteriota = subset_taxa(physeq1, Phylum=="Actinobacteriota")
physeq.Proteobacteria = subset_taxa(physeq1, Phylum=="Proteobacteria")
physeq.Acidobacteriota = subset_taxa(physeq1, Phylum=="Acidobacteriota")


physeq.Entomoplasmataceae = subset_taxa(physeq1, Family=="Entomoplasmataceae")

plot_tree(physeq.Acidobacteriota, color="WorkerType", shape="FeS", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
plot_tree(physeq.Entomoplasmataceae, color="FeS", shape="WorkerType", label.tips="taxa_names", ladderize="left", plot.margin=0.3)

treeAcido <- plot_tree(physeq.Acidobacteriota, color="WorkerType", shape="FeS", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
pdf("treeAcido.pdf")
print(treeAcido)
dev.off()

unifrac <- UniFrac(physeq1, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)


ord.pcoa.unifrac1 <- ordinate(physeq1, method="PCoA", distance="unifrac", weighted=TRUE)
weighted <- plot_ordination(physeq1, ord.pcoa.unifrac1, color="FeS", shape="WorkerType", title="PCoA - Weighted UniFrac") + stat_ellipse(aes(color=FeS),type = "norm", linetype=2)
pdf("weighted_pcoa.pdf")
print(weighted)
dev.off()

weighted2 <- weighted +  stat_ellipse(type = "t", linetype = 2) ##+ stat_ellipse(type = "t")
pdf("weighted_pcoa.pdf")
print(weighted2)
dev.off()

ord.pcoa.unifrac2 <- ordinate(physeq1, method="PCoA", distance="unifrac", weighted=FALSE)
unweighted <- plot_ordination(physeq1, ord.pcoa.unifrac2, color="FeS", shape="WorkerType", title="PCoA - Unweighted UniFrac") + stat_ellipse(type = "norm", linetype=2)
pdf("unweighted_pcoa.pdf")
print(unweighted)
dev.off()

unweighted2 <- unweighted +  stat_ellipse(type = "norm", linetype = 2) ##+ stat_ellipse(type = "t")
pdf("weighted_pcoa.pdf")
print(unweighted2)
dev.off()
