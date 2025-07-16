### BETA DIVERSITY - Geographic Comparison of Gut Microbial Flexibility in TRES
# Microbiome Analysis
# Written by: Jenn Houtz, updated by Natalie Morris
# Modified from code by: Worsley et al. 2021, Animal Microbiome
# Taxa bar plot basics from Gabri Schiro
# Last updated: 7/16/2025

# need to remove singletons, if not already gone


# Load libraries

pacman::p_load("ggplot2", "grid", "gridExtra", "emmeans", "phyloseq",
               "vegan", "devtools", "sjPlot", "here", "ape", "philr", "pairwiseAdonis",
               "ggordiplots", "BiodiversityR", "dplyr", "tidyr")
library(vegan)
library(ggplot2)
if(!require("ggpubr")){      # an if condition, if ggpubr is not installed:
  install.packages("ggpubr") #install ggpubr
}
if(!require("reshape2")){      # an if condition, if package is not installed:
  install.packages("reshape2") #install package
}
library(reshape2)
library(phyloseq)

physeqBeta_rare <- readRDS("~/Library/CloudStorage/OneDrive-UniversityofArizona/Cornell/ht_njm/ht_njm_R/honors_thesis_njm/physeqBeta_rare_taxabarplot.rds")
# Set working directory


##### Create taxa bar plots #####


################# UPDATED ###############
#### Phyla - incubation
# Taxa bar plot for incubation
physeq_incub_rare<- subset_samples(physeqBeta_rare, Capture_Number=="1")
#physeq_incub_rare<- subset_samples(physeq_incub_rare, Individual_Treatment=="Control")


# Extract taxonomy table
TAXtab <- as.data.frame(tax_table(physeq_incub_rare))

# Extract ASV (OTU) table
ASVtab <- as.data.frame(otu_table(physeq_incub_rare))

# Extract sample metadata
metadata <- as.data.frame(sample_data(physeq_incub_rare))

sort(rowSums(ASVtab))


phylum_tab = aggregate(t(ASVtab) ~ TAXtab$Phylum, FUN = "sum") ## sum all the ASV together in the same sample
rownames(phylum_tab) = phylum_tab$`TAXtab$Phylum`
phylum_tab=phylum_tab[,-1]
total_ab_phyla = rowSums(phylum_tab)/sum(rowSums(phylum_tab))
total_ab_phyla = total_ab_phyla*100
names(total_ab_phyla) = rownames(phylum_tab)
total_ab_phyla
top10 = names(sort(total_ab_phyla, decreasing = TRUE))
top10 = top10[1:10]
top10

phylum_tab$Phylum = rownames(phylum_tab)
phylum_tab_melted = melt(phylum_tab, value = phylum_tab$`TAXtab$Phylum`) ## This melt the table

colnames(phylum_tab_melted)[1:2] = c("Phylum", "sample_names")
phylum_tab_melted$group = metadata$Location[match(phylum_tab_melted$sample_names, metadata$sample_names)]

phylum_tab_melted1 = phylum_tab_melted[!is.na(phylum_tab_melted$group),]
phylum_tab_melted1$Phylum = ifelse(phylum_tab_melted1$Phylum %in% top10,
                                   phylum_tab_melted1$Phylum, "others")


ggplot(phylum_tab_melted1,aes(x = sample_names, y = value, fill = Phylum))+
  geom_bar(position="fill", stat="identity", width = 1, color = NA)+
  facet_grid(.~group, scales = "free") + 
  scale_fill_brewer(palette = "Paired") + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "right",
        axis.ticks.x = element_blank()
  ) +
  ylab("Relative abundance") +
  labs(title = "Phyla during Incubation",
       x = "Sample ID",
       y = "Relative Abundance") +
  scale_y_continuous(expand = c(0, 0)) 
  
  
  


  # Family level
  f_tab = aggregate(t(ASVtab) ~ TAXtab$Family, FUN = "sum") ## sum all the ASV together in the same sample
  rownames(f_tab) = f_tab$'TAXtab$Family'
  f_tab=f_tab[,-1]
  total_ab_f = rowSums(f_tab)/sum(rowSums(f_tab))
  total_ab_f = total_ab_f*100
  names(total_ab_f) = rownames(f_tab)
  total_ab_f
  top10 = names(sort(total_ab_f, decreasing = TRUE))
  top10 = top10[1:10]
  top10
  
  f_tab$Family = rownames(f_tab)
  f_tab_melted = melt(f_tab, value = Family) ## This melt the table
  
  colnames(f_tab_melted)[1:2] = c("Family", "sample_name")
  metadata$sample_name = as.factor(rownames(metadata))
  class(f_tab_melted$sample_name)
  f_tab_melted$group = metadata$Location[match(f_tab_melted$sample_name, metadata$sample_names)]
  
  f_tab_melted1 = f_tab_melted[!is.na(f_tab_melted$group),]
  f_tab_melted1$Family = ifelse(f_tab_melted1$Family %in% top10,
                                f_tab_melted1$Family, "others")
  
  
  ggplot(f_tab_melted1,aes(x = sample_name, y = value, fill = Family))+
    geom_bar(position="fill", stat="identity", width = 1, color = NA)+
    facet_grid(.~group, scales = "free") + 
    scale_fill_brewer(palette = "Paired") + 
    theme_classic() +
    theme(axis.text.x = element_blank(),
          legend.position = "right",
          axis.ticks.x = element_blank()
    ) +
    ylab("Relative abundance") +
    labs(title = "Family during Incubation",
         x = "Sample ID",
         y = "Relative Abundance") +
    scale_y_continuous(expand = c(0, 0))
  
  
# Genera

    g_tab = aggregate(t(ASVtab) ~ TAXtab$Genus, FUN = "sum") ## sum all the ASV together in the same sample
    rownames(g_tab) = g_tab$`TAXtab$Genus`
    g_tab = g_tab[,-1]
    total_ab_g = rowSums(g_tab)/sum(rowSums(g_tab))
    total_ab_g = total_ab_g*100
    names(total_ab_g) = rownames(g_tab)
    total_ab_g
    top10 = names(sort(total_ab_g, decreasing = TRUE))
    top10 = top10[1:10]
    top10
    
    g_tab$Genus = rownames(g_tab)
    g_tab_melted = melt(g_tab, value = g_tab$`TAXtab$Genus`) ## This melt the table
    
    colnames(g_tab_melted)[1:2] = c("Genus", "sample_name")
    g_tab_melted$group = metadata$Location[match(g_tab_melted$sample_name, metadata$sample_name)]
    
    g_tab_melted1 = g_tab_melted[!is.na(g_tab_melted$group),]
    g_tab_melted1$Genus = ifelse(g_tab_melted1$Genus %in% top10,
                                 g_tab_melted1$Genus, "others")
    
    
    library(ggplot2)
    library(RColorBrewer)
    
    color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(g_tab_melted1$Genus)))
    
    ggplot(g_tab_melted1, aes(x = sample_name, y = value, fill = Genus)) +
       geom_bar(position="fill", stat="identity", width = 1, color = NA)+
      scale_fill_brewer(palette = "Paired") + 
      facet_grid(.~group, scales = "free") + 
      theme_classic() +
      theme(axis.text.x = element_blank(),
            legend.position = "right",
            axis.ticks.x = element_blank()
      ) +
      ylab("Relative abundance") +
      labs(title = "Genus during Incubation",
           x = "Sample ID",
           y = "Relative Abundance") +
      scale_y_continuous(expand = c(0, 0)) 
    
    

# Provisioning: physeqBeta_rareP (control only)

    # subset !
    physeqBetarareP<- subset_samples(physeqBeta_rare, Capture_Number=="3")
    physeqBeta_rareP<- subset_samples(physeqBetarareP, Individual_Treatment=="Control")
    
    
# Taxa bar plot for provisioning

# Extract taxonomy table
TAXtab <- as.data.frame(tax_table(physeqBeta_rareP))

# Extract ASV (OTU) table
ASVtab <- as.data.frame(otu_table(physeqBeta_rareP))

# Extract sample metadata
metadata <- as.data.frame(sample_data(physeqBeta_rareP))

sort(rowSums(ASVtab))

phylum_tab = aggregate(t(ASVtab) ~ TAXtab$Phylum, FUN = "sum") ## sum all the ASV together in the same sample
rownames(phylum_tab) = phylum_tab$`TAXtab$Phylum`
phylum_tab=phylum_tab[,-1]
total_ab_phyla = rowSums(phylum_tab)/sum(rowSums(phylum_tab))
total_ab_phyla = total_ab_phyla*100
names(total_ab_phyla) = rownames(phylum_tab)
total_ab_phyla
top10 = names(sort(total_ab_phyla, decreasing = TRUE))
top10 = top10[1:10]
top10

phylum_tab$Phylum = rownames(phylum_tab)
phylum_tab_melted = melt(phylum_tab, value = phylum_tab$`TAXtab$Phylum`) ## This melt the table

colnames(phylum_tab_melted)[1:2] = c("Phylum", "sample_names")
phylum_tab_melted$group = metadata$Location[match(phylum_tab_melted$sample_names, metadata$sample_names)]

phylum_tab_melted1 = phylum_tab_melted[!is.na(phylum_tab_melted$group),]
phylum_tab_melted1$Phylum = ifelse(phylum_tab_melted1$Phylum %in% top10,
                                   phylum_tab_melted1$Phylum, "others")


ggplot(phylum_tab_melted1,aes(x = sample_names, y = value, fill = Phylum))+
  geom_bar(position="fill", stat="identity", width = 1, color = NA)+
  facet_grid(.~group, scales = "free") + 
  scale_fill_brewer(palette = "Paired") + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "right",
        axis.ticks.x = element_blank()
  ) +
  ylab("Relative abundance") +
  labs(title = "Phyla during Provisioning (control birds)",
       x = "Sample ID",
       y = "Relative Abundance") +
  scale_y_continuous(expand = c(0, 0)) 


# Family level
f_tab = aggregate(t(ASVtab) ~ TAXtab$Family, FUN = "sum") ## sum all the ASV together in the same sample
rownames(f_tab) = f_tab$'TAXtab$Family'
f_tab=f_tab[,-1]
total_ab_f = rowSums(f_tab)/sum(rowSums(f_tab))
total_ab_f = total_ab_f*100
names(total_ab_f) = rownames(f_tab)
total_ab_f
top10 = names(sort(total_ab_f, decreasing = TRUE))
top10 = top10[1:10]
top10

f_tab$Family = rownames(f_tab)
f_tab_melted = melt(f_tab, value = Family) ## This melt the table

colnames(f_tab_melted)[1:2] = c("Family", "sample_name")
metadata$sample_name = as.factor(rownames(metadata))
class(f_tab_melted$sample_name)
f_tab_melted$group = metadata$Location[match(f_tab_melted$sample_name, metadata$sample_names)]

f_tab_melted1 = f_tab_melted[!is.na(f_tab_melted$group),]
f_tab_melted1$Family = ifelse(f_tab_melted1$Family %in% top10,
                              f_tab_melted1$Family, "others")


ggplot(f_tab_melted1,aes(x = sample_name, y = value, fill = Family))+
  geom_bar(position="fill", stat="identity", width = 1, color = NA)+
  facet_grid(.~group, scales = "free") + 
  scale_fill_brewer(palette = "Paired") + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "right",
        axis.ticks.x = element_blank()
  ) +
  ylab("Relative abundance") +
  labs(title = "Family during Provisioning (control birds)",
       x = "Sample ID",
       y = "Relative Abundance") +
  scale_y_continuous(expand = c(0, 0))


# Genera

g_tab = aggregate(t(ASVtab) ~ TAXtab$Genus, FUN = "sum") ## sum all the ASV together in the same sample
rownames(g_tab) = g_tab$`TAXtab$Genus`
g_tab = g_tab[,-1]
total_ab_g = rowSums(g_tab)/sum(rowSums(g_tab))
total_ab_g = total_ab_g*100
names(total_ab_g) = rownames(g_tab)
total_ab_g
top10 = names(sort(total_ab_g, decreasing = TRUE))
top10 = top10[1:10]
top10

g_tab$Genus = rownames(g_tab)
g_tab_melted = melt(g_tab, value = g_tab$`TAXtab$Genus`) ## This melt the table

colnames(g_tab_melted)[1:2] = c("Genus", "sample_name")
g_tab_melted$group = metadata$Location[match(g_tab_melted$sample_name, metadata$sample_name)]

g_tab_melted1 = g_tab_melted[!is.na(g_tab_melted$group),]
g_tab_melted1$Genus = ifelse(g_tab_melted1$Genus %in% top10,
                             g_tab_melted1$Genus, "others")

ggplot(g_tab_melted1, aes(x = sample_name, y = value, fill = Genus)) +
  geom_bar(position="fill", stat="identity", width = 1, color = NA)+
  scale_fill_brewer(palette = "Paired") + 
  facet_grid(.~group, scales = "free") + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "right",
        axis.ticks.x = element_blank()
  ) +
  ylab("Relative abundance") +
  labs(title = "Genus during Provisioning (control birds)",
       x = "Sample ID",
       y = "Relative Abundance") +
  scale_y_continuous(expand = c(0, 0)) 

########older

# Phyla

      phylum_tab = aggregate(t(ASVtab) ~ TAXtab$Phylum, FUN = "sum") ## sum all the ASV together in the same sample
      rownames(phylum_tab) = phylum_tab$`TAXtab$Phylum`
      total_ab_phyla = rowSums(phylum_tab1)/sum(rowSums(phylum_tab1))
      total_ab_phyla = total_ab_phyla*100
      names(total_ab_phyla) = phylum_tab$`TAXtab$Phylum`
      total_ab_phyla
      top10 = names(sort(total_ab_phyla, decreasing = TRUE))
      top10 = top10[1:10]
      top10
      
      
      phylum_tab_melted = melt(phylum_tab, value = phylum_tab$`TAXtab$Phylum`) ## This melt the table
      
      colnames(phylum_tab_melted)[1:2] = c("Phylum", "sample_name")
      phylum_tab_melted$group = metadata$Location[match(phylum_tab_melted$sample_name, metadata$sample_name)]
      
      phylum_tab_melted1 = phylum_tab_melted[!is.na(phylum_tab_melted$group),]
      phylum_tab_melted1$Phylum = ifelse(phylum_tab_melted1$Phylum %in% top10,
                                         phylum_tab_melted1$Phylum, "others")
      
      
      ggplot(phylum_tab_melted1,aes(x = group, y = value, fill = Phylum))+
        geom_bar(position="fill", stat="identity")+
        scale_fill_brewer(palette = "Paired") + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ylab("Relative abundance") +
        theme_classic() + labs(title = "Phylum at provisioning",
                               x = "Location") 

# Genera
    
      g_tab = aggregate(t(ASVtab) ~ TAXtab$Genus, FUN = "sum") ## sum all the ASV together in the same sample
      rownames(g_tab) = g_tab$`TAXtab$Genus`
      total_ab_g = rowSums(g_tab1)/sum(rowSums(g_tab1))
      total_ab_g = total_ab_g*100
      names(total_ab_g) = g_tab$`TAXtab$Genus`
      total_ab_g
      top15 = names(sort(total_ab_g, decreasing = TRUE))
      top15 = top15[1:15]
      top15
      
      
      g_tab_melted = melt(g_tab, value = g_tab$`TAXtab$Genus`) ## This melt the table
      
      colnames(g_tab_melted)[1:2] = c("Genus", "sample_name")
      g_tab_melted$group = metadata$Location[match(g_tab_melted$sample_name, metadata$sample_name)]
      
      g_tab_melted1 = g_tab_melted[!is.na(g_tab_melted$group),]
      g_tab_melted1$Genus = ifelse(g_tab_melted1$Genus %in% top15,
                                   g_tab_melted1$Genus, "others")
      
      
      library(ggplot2)
      library(RColorBrewer)
      
      color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(g_tab_melted1$Genus)))
      
      ggplot(g_tab_melted1, aes(x = group, y = value, fill = Genus)) +
        geom_bar(position = "fill", stat = "identity") +
        scale_fill_manual(values = color_palette) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ylab("Relative abundance") +
        theme_classic() +
        labs(title = "Genera at Provisioning", x = "Location")
    
    
    library(ggplot2)
    library(RColorBrewer)
    
    # Generate 15 distinct colors (you can change "Set3" to other palettes)
    color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(g_tab_melted1$Genus)))
    
    ggplot(g_tab_melted1, aes(x = group, y = value, fill = Genus)) +
      geom_bar(position = "fill", stat = "identity") +
      scale_fill_manual(values = color_palette) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylab("Relative abundance") +
      theme_classic() +
      labs(title = "Taxa - Provisioning - Genus", x = "Location")
  


#Upload ASV, taxonomy, tree and metadata files to phyloseq:

#ASV table
unrarefied_ASV_table <- read.csv ("ASV_counts.csv", row.names=1) #read in asv table with feature names as rownames
str (unrarefied_ASV_table) #10776 obs. of  534 variables
head (unrarefied_ASV_table)
unrarefied_ASV_table <- as.matrix (unrarefied_ASV_table) #make into a matrix

#taxonomy file
taxonomy <- read.csv ("taxonomy.csv", row.names=1)
taxonomy <- separate(taxonomy, Taxon, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep="; ")
str (taxonomy) #11997 obs. of  7 variables
taxonomy <- as.matrix (taxonomy)

#read in tree as a phyloseq object
phy_tree <- read_tree ("tree.nwk")

#load metadata file
#before importing, make the following changes:
# CANNOT BE ANY ZEROS!!! or NA!!!
# add quant_reading column, cross reference plate sample IDs and DNA quantification readings from sequencing facility
# for any 0s in quant_reading (PCR blanks), change to arbitrary number = 0.0001 for decontam step to work
# add sample_or_control column, for all blanks = negative_control, for all actual samples = true_sample
# sort samples into alphabetical order, ensure that all blanks have names identical to that of the ASV table
metadatafull <- read.csv("metadata_geo_njm.csv") 
str(metadatafull)  # contains 445 samples - (should match ASV table observations)
row.names(metadatafull)<-metadatafull$microbiome_sample_nb

#import all as phyloseq objects
ASV <- otu_table(unrarefied_ASV_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
META <- sample_data(metadatafull)
head(META)


#check that the ASV and sample names are consistent across objects
# i.e. make sure blanks are labeled the same in metadata and ASV table
# check formatting of sample and taxa names, taxa names do not need to be in the same order
# taxa and phy_tree should have the same number of features, ASV should have a different, but similar number
str(taxa_names(TAX)) 
str(taxa_names(ASV))
str(taxa_names(phy_tree))

str(sample_names(ASV))
str(sample_names(META))

#### MERGE INTO PHYLOSEQ OBJECT ####
physeq <- phyloseq(ASV, TAX, META, phy_tree) #ignore error "found more than one..."
physeq #10776 taxa and 445 samples #### 


#check the number of reads per sample
sample_reads<-data.frame(reads=sample_sums(physeq))
head(sample_reads)

##### Decontam #####

library(decontam)

# Remove contaminants with the combined function of prevalence and frequency

sample_data(physeq)$is.neg <- sample_data(physeq)$sample_or_control == "negative_control" # metadata column, blanks as negative_control, samples as true_sample
contamdf.prev <- isContaminant(physeq, method="combined",conc="quant_reading",neg="is.neg") # blanks must have number greater than 0, set to 0.0001
table(contamdf.prev$contaminant) #table--true column = number of contaminants

# Creating the decontaminated phyloseq object (physeq_nc)

physeq_nc<-prune_taxa(!contamdf.prev$contaminant,physeq) 

physeq_nc #identified 62 contaminants and removed them leaving 10714 taxa

# Remove negative controls (n=17)

physeq2 <- subset_samples(physeq_nc, sample_or_control == "true_sample") # ignore error about phylo
physeq2 #428 samples, 10714 taxa


# Read in phyloseq

prevdf = apply(X = otu_table(geo_phylo),
               MARGIN = ifelse(taxa_are_rows(geo_phylo), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts for each phylum to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(geo_phylo),
                    tax_table(geo_phylo))
head(prevdf)
str(prevdf)

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) #average prevelance of features within each phylum and the sum of feature prevalence within each phylum

# Plot the unique phyla: each dot will be a feature- total abundance of that feature across samples on x axis and the prevalance (the fraction of all samples it occurs in on the y axis).
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(geo_phylo, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(geo_phylo),color=Phylum)) +
  # Include filtering threshold here; 2%, total abundance across all samples set to 50
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_vline(xintercept = 50, alpha = 0.5, linetype = 2)+  geom_point(size = 1, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")


#Export feature/OTU table as a biom table for the post-filtering sequence and ASV metrics
######DO NOT TRANSFORM OR THE FEATURES AND SAMPLES WILL BE SWITCHED

library(biomformat);packageVersion("biomformat")

otu<-as(otu_table(geo_phylo),"matrix")
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"otu_biom_postfiltering_2025.biom")

##### RAREFY READS TO MIN SAMPLING DEPTH ######

# Look at lowest sample depth

sample_sums(physeq) # lowest M17A0023 = 58 reads 

#rarefy to 1000 and set seed before rarefying (28367), so that results are reproducible  
physeqRare<-rarefy_even_depth(physeq, 5000, rngseed = 28367)
# 15 samples removed

# OTUs removed after subsampling- leaves 5550 taxa and 413 samples
physeqRare
sample_sums(physeqRare)

saveRDS(physeqRare, "physeqR.rds")
##### BETA DIVERSITY ######

##### PERMANOVA - Bray Curtis location #####

# Use the rarefied phyloseq object physeqRare

##########################
### RAREFIED PRUNE RARE TAXA ######
##########################

# Compute prevalence of each feature (total number of samples in which a taxon appears at least once), store as data.frame
prevdf_rare = apply(X = otu_table(physeqRare),
                    MARGIN = ifelse(taxa_are_rows(physeqRare), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts for each phylum to this data.frame
prevdf_rare = data.frame(Prevalence = prevdf_rare,
                         TotalAbundance = taxa_sums(physeqRare),
                         tax_table(physeqRare))
head(prevdf_rare)
str(prevdf_rare)

plyr::ddply(prevdf_rare, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) #average prevelance of features within each phylum and the sum of feature prevalence within each phylum

# Plot the unique phyla: each dot will be a feature- total abundance of that feature across samples on x axis and the prevalance (the fraction of all samples it occurs in on the y axis).
prevdf1_rare = subset(prevdf_rare, Phylum %in% get_taxa_unique(physeqRare, "Phylum"))
ggplot(prevdf1_rare, aes(TotalAbundance, Prevalence / nsamples(physeqRare),color=Phylum)) +
  # Include filtering threshold here; 2%, total abundance across all samples set to 50
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_vline(xintercept = 50, alpha = 0.5, linetype = 2)+  geom_point(size = 1, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# Define abundance as 10 total reads across samples

abundanceThreshold<-10

# Execute the prevalence filter, using `prune_taxa()` function

head(prevdf1_rare)

KeepTaxa1_rare<- rownames(prevdf_rare)[(prevdf_rare$TotalAbundance >= abundanceThreshold)]
str(KeepTaxa1_rare) #3322 taxa

physeqBeta_rare<- prune_taxa(KeepTaxa1_rare,physeqRare)
physeqBeta_rare # 3322 taxa, 393 samples

physeqBeta_rare1 <- subset_samples(physeqBeta_rare, Location!="Wyoming")
physeqBeta_rare2 <- subset_samples(physeqBeta_rare, Capture_Number!="2" )
physeqBeta_rareC <- subset_samples(physeqBeta_rare, Individual_Treatment =="Control")
physeqBeta_rareC



physeqBeta_rareI <- subset_samples(physeqBeta_rare2, Capture_Number == "1")

otu<-as(otu_table(physeqBeta_rareI),"matrix")
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"otu_biom_rarefied_incubation.biom")

physeqBeta_rareP <- subset_samples(physeqBeta_rareC, Capture_Number == "3")

otu<-as(otu_table(physeqBeta_rareP),"matrix")
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"otu_biom_rarefied_provisioning.biom")

###########################################################################################

##### Export rarefied ASV table and taxonomy into qiime2 for taxabarplot #####

#Export taxonomy table as "tax.txt"

tax<-as(tax_table(physeqBeta_rare),"matrix")
tax_cols <- colnames(tax)
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "tax_2025.txt", quote=FALSE, col.names=FALSE, sep="\t")

#Export feature/OTU table as a biom table
######DO NOT TRANSFORM OR THE FEATURES AND SAMPLES WILL BE SWITCHED

library(biomformat);packageVersion("biomformat")

otu<-as(otu_table(physeqBeta_rareC),"matrix")
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"otu_biom_rarefiedC.biom")

# Subset to Incubation only

physeq_incub_rare<- subset_samples(physeqBeta_rare, Capture_Number=="1")
physeq_incub_rare #238 samples

#Extract metadata

incub_beta_data<- data.frame(sample_data(physeq_incub_rare))

# Look at structure and make sure they are correct data types

str(incub_beta_data)
incub_beta_data$Location <- factor(incub_beta_data$Location)

##### PCoA Incubation location Bray #####
location_colors <- c("#4682B4", "#B31B1B", "#f77f00", "darkgreen")

library(ggforce)
ordu = ordinate(physeq_incub_rare, "PCoA", "bray", weighted=FALSE)
P <- plot_ordination(physeq_incub_rare, ordu, color="Location")
pca_incub_bray <- P + geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_bw() +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_mark_ellipse(aes(color = Location,
                        label=Location),
                    expand = unit(0.5,"mm")) +
  scale_color_manual(values = location_colors)

pca_incub_bray


ggsave("pca_incb_bray.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches

# Calculate Bray-Curtis dissimiliary
locationI_bray <- phyloseq::distance(physeq_incub_rare, method = "bray")

perm <- how(nperm = 9999)
set.seed(5498)
permanova_incub_bray<- adonis2(locationI_bray ~ Location, data=incub_beta_data, 
                               permutations = perm, method = "bray", by = "margin")
permanova_incub_bray

#Export results

str(permanova_incub_bray)
perm_results_incub_bray <- data.frame(permanova_incub_bray[c(1,2,3,4,5)])
perm_results_incub_bray
write.csv(perm_results_incub_bray, "permanova_results_incub_bray.csv")

############################  
##### Pairwise permanova - location bray-curtis ######
############################ 

library(devtools)

# run this code before installing packages from Github AHHHHH
# Sys.setenv(R_REMOTES_STANDALONE="true", build = F)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

perm <- how(nperm = 9999)
set.seed(549838)
pairwise_perm_incub_bray <- pairwise.adonis2(locationI_bray ~ Location,
                                             data = incub_beta_data, method="euclidean", by="margin", nperm=perm) 
pairwise_perm_incub_bray  

# P value adjustment for multiple comparisons (change these to your p values)
p.adjusted_pairwise_perm_incub_bray <- p.adjust(c(0.001,0.001,0.001,0.005,0.001, 0.001),method="BH")
p.adjusted_pairwise_perm_incub_bray

# Export results

write.csv(pairwise_perm_incub_bray, "pairwise_perm_incub_bray2.csv")

write.csv(p.adjusted_pairwise_perm_incub_bray, "p.adjusted_pairwise_perm_incub_bray2.csv")


###### Subset to Provisioning only #####

physeq_prov_rare<- subset_samples(physeqBeta_rare, Capture_Number=="3")
physeq_prov_rare #12 samples

#Extract metadata

prov_beta_data<- data.frame(sample_data(physeq_prov_rare))

# Look at structure and make sure they are correct data types

str(prov_beta_data)
prov_beta_data$Location <- factor(prov_beta_data$Location)

##### PCoA Provisioning location Bray #####

ordu = ordinate(physeq_prov_rare, "PCoA", "bray", weighted=FALSE)
P <- plot_ordination(physeq_prov_rare, ordu, color="Location")
pca_prov_bray <- P + geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_bw() +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_mark_ellipse(aes(color = Location,
                        label=Location),
                    expand = unit(0.5,"mm"))
pca_prov_bray


ggsave("pca_prov_bray.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches

# Calculate Bray-Curtis dissimiliary
locationP_bray <- phyloseq::distance(physeq_prov_rare, method = "bray")

perm <- how(nperm = 9999)
set.seed(5498)
permanova_prov_bray<- adonis2(locationP_bray ~ Location, data=prov_beta_data, 
                              permutations = perm, method = "bray", by = "margin")
permanova_prov_bray

#Export results

str(permanova_prov_bray)
perm_results_prov_bray <- data.frame(permanova_prov_bray[c(1,2,3,4,5)])
perm_results_prov_bray
write.csv(perm_results_prov_bray, "permanova_results_prov_bray.csv")

############################  
##### Pairwise permanova - location bray-curtis ######
############################ 

library(devtools)

# run this code before installing packages from Github AHHHHH
# Sys.setenv(R_REMOTES_STANDALONE="true", build = F)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

perm <- how(nperm = 9999)
set.seed(549838)
pairwise_perm_prov_bray <- pairwise.adonis2(locationP_bray ~ Location,
                                            data = prov_beta_data, method="euclidean", by="margin", nperm=perm) 
pairwise_perm_prov_bray  

# P value adjustment for multiple comparisons (change these to your p values)!!!!!!!!!!!!
p.adjusted_pairwise_perm_prov_bray <- p.adjust(c(0.001, 0.001, 0.001),method="BH")
p.adjusted_pairwise_perm_prov_bray

# Export results

write.csv(pairwise_perm_prov_bray, "pairwise_perm_prov_bray2.csv")

write.csv(p.adjusted_pairwise_perm_prov_bray, "p.adjusted_pairwise_perm_prov_bray2.csv")

##### Subset to Filtered, Incubation, Control only #####

physeq_incub_rare<- subset_samples(physeqBeta_rare, Capture_Number=="1")
physeq_incubC_rare<- subset_samples(physeq_incub_rare, Individual_Treatment == "Control")
physeq_incubC_rare<- subset_samples(physeq_incubC_rare, Location != "Wyoming")
physeq_incubC_rare #60 samples

#Extract metadata

incubC_beta_data<- data.frame(sample_data(physeq_incubC_rare))
table(incubC_beta_data$Location)
# Look at structure and make sure they are correct data types

str(incubC_beta_data)
incubC_beta_data$Location <- factor(incubC_beta_data$Location)

##### PCoA IncubationC location Bray #####
location_colors <- c("#4682B4", "#B31B1B", "#f77f00")

library(ggforce)
ordu = ordinate(physeq_incubC_rare, "PCoA", "bray", weighted=FALSE)
PI <- plot_ordination(physeq_incubC_rare, ordu, color = "Location")
pca_incubC_bray <- PI + geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_bw() +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_mark_ellipse(aes(color = Location,
                        label=Location),
                    expand = unit(0.5,"mm"))+
  scale_color_manual(values = location_colors)

pca_incubC_bray


ggsave("pca_incubC_bray.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches

# Calculate Bray-Curtis dissimiliary
locationIC_bray <- phyloseq::distance(physeq_incubC_rare, method = "bray")

perm <- how(nperm = 9999)
set.seed(5498)
permanova_incubC_bray<- adonis2(locationIC_bray ~ Location, data=incubC_beta_data, 
                                permutations = perm, method = "bray", by = "margin")
permanova_incubC_bray

#Export results

str(permanova_incubC_bray)
perm_results_incubC_bray <- data.frame(permanova_incubC_bray[c(1,2,3,4,5)])
perm_results_incubC_bray
write.csv(perm_results_incubC_bray, "permanova_results_incubC_bray.csv")

############################  
##### Pairwise permanova - location bray-curtis ######
############################ 

library(devtools)

# run this code before installing packages from Github AHHHHH
# Sys.setenv(R_REMOTES_STANDALONE="true", build = F)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

perm <- how(nperm = 9999)
set.seed(549838)
pairwise_perm_incubC_bray <- pairwise.adonis2(locationIC_bray ~ Location,
                                              data = incubC_beta_data, method="euclidean", by="margin", nperm=perm) 
pairwise_perm_incubC_bray  

# P value adjustment for multiple comparisons (change these to your p values)
p.adjusted_pairwise_perm_incubC_bray <- p.adjust(c(0.001,0.001,0.131),method="BH")
p.adjusted_pairwise_perm_incubC_bray

# Export results

write.csv(pairwise_perm_incubC_bray, "pairwise_perm_incubC_bray.csv")

write.csv(p.adjusted_pairwise_perm_incubC_bray, "p.adjusted_pairwise_perm_incubC_bray.csv")


##### Subset to Filtered, Provisioning, Control only #####

physeq_prov_rare<- subset_samples(physeqBeta_rare, Capture_Number=="3")
physeq_provC_rare<- subset_samples(physeq_prov_rare, Individual_Treatment == "Control")
physeq_provC_rare #53 samples


#Extract metadata

provC_beta_data<- data.frame(sample_data(physeq_provC_rare))
table(provC_beta_data$Location)
# Look at structure and make sure they are correct data types

str(provC_beta_data)
provC_beta_data$Location <- factor(provC_beta_data$Location)

##### PCoA Provisioning C location Bray #####

ordu = ordinate(physeq_provC_rare, "PCoA", "bray", weighted=FALSE)
PP <- plot_ordination(physeq_provC_rare, ordu, color="Location")
pca_provC_bray <- PP + geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_bw() +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_mark_ellipse(aes(color = Location,
                        label=Location),
                    expand = unit(0.5,"mm")) +
  scale_color_manual(values = location_colors)
pca_provC_bray


ggsave("pca_provC_bray.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches

# Calculate Bray-Curtis dissimiliary
locationPC_bray <- phyloseq::distance(physeq_provC_rare, method = "bray")

perm <- how(nperm = 9999)
set.seed(5498)
permanova_provC_bray<- adonis2(locationPC_bray ~ Location, data=provC_beta_data, 
                               permutations = perm, method = "bray", by = "margin")
permanova_provC_bray

#Export results

str(permanova_provC_bray)
perm_results_provC_bray <- data.frame(permanova_provC_bray[c(1,2,3,4,5)])
perm_results_provC_bray
write.csv(perm_results_provC_bray, "permanova_results_provC_bray.csv")

############################  
##### Pairwise permanova - location bray-curtis ######
############################ 

library(devtools)
library(pairwiseAdonis)

perm <- how(nperm = 9999)
set.seed(549838)
pairwise_perm_provC_bray <- pairwise.adonis2(locationPC_bray ~ Location,
                                             data = provC_beta_data, method="euclidean", by="margin", nperm=perm) 
pairwise_perm_provC_bray  

# P value adjustment for multiple comparisons (change these to your p values)!!!!!!!!!!!!
# This is the value to report in draft
p.adjusted_pairwise_perm_provC_bray <- p.adjust(c(0.008, 0.008, 0.0034),method="BH")
p.adjusted_pairwise_perm_provC_bray

# Export results

write.csv(pairwise_perm_provC_bray, "pairwise_perm_provC_bray.csv")

write.csv(p.adjusted_pairwise_perm_provC_bray, "p.adjusted_pairwise_perm_provC_bray.csv")


########

combined_plot <- grid.arrange(pca_provC_bray, pca_incubC_bray, ncol = 2)
# Combine both plots on the same graph
combined_plot <- pca_incubC_bray + pca_provC_bray

# Display the combined plot
print(combined_plot)

combined_plot <- pca_provC_bray + geom_point(data = ggplot_build(pca_incubC_bray)$data[[1]],
                                             aes(x = x, y = y), alpha = 0.5)

# Display the combined plot
print(combined_plot)


########
citation('vegan')

##### SIMPER
beta_data<- data.frame(sample_data(physeqBeta_rareC))

physeqBeta_rareC

C_otus = as(otu_table(physeqBeta_rareC), "matrix")

if(taxa_are_rows(physeqBeta_rareC)){C_otus <- t(C_otus)}

# Coerce the object to a data.frame
C_OTUs_scaled = as.data.frame(C_otus)

# running the simper analysis on the dataframe and the variable of interest "Location"
C_simper <- simper(C_OTUs_scaled, beta_data$Capture_Number, permutations = 100)

# printing the top OTUs
print(C_simper)
summary(C_simper)

##### SIMPER - Incubation

incubC_otus = as(otu_table(physeq_incubC_rare), "matrix")

if(taxa_are_rows(physeq_incubC_rare)){incubC_otus <- t(incubC_otus)}

# Coerce the object to a data.frame
incubC_OTUs_scaled = as.data.frame(incubC_otus)

# running the simper analysis on the dataframe and the variable of interest "Location"
IC_simper <- simper(incubC_OTUs_scaled, incubC_beta_data$Location, permutations = 999)

# printing the top OTUs
print(IC_simper)
summary(IC_simper, digits = max(3,getOption("digits") - 3))

summary(IC_simper, digits = max(3,getOption("digits") - 3))$Alaska_Tennessee


#### SIMPER - Provisioning

provC_otus = as(otu_table(physeq_provC_rare), "matrix")

if(taxa_are_rows(physeq_provC_rare)){provC_otus <- t(provC_otus)}

# Coerce the object to a data.frame
provC_OTUs_scaled = as.data.frame(provC_otus)

# running the simper analysis on the dataframe and the variable of interest "Location"
PC_simper <- simper(provC_OTUs_scaled, provC_beta_data$Location, permutations = 100)

# printing the top OTUs
print(PC_simper)


#### taxa bar plot - older version

##### old version

# Extract taxonomy table
TAXtab <- as.data.frame(tax_table(physeqBeta_rare))

# Extract ASV (OTU) table
ASVtab <- as.data.frame(otu_table(physeqBeta_rare))

# Extract sample metadata
metadata <- as.data.frame(sample_data(physeqBeta_rare))

#saveRDS(physeqBeta_rare, "physeqBeta_rare_taxabarplot.rds")


sort(rowSums(ASVtab))

phylum_tab = aggregate(t(ASVtab) ~ TAXtab$Phylum, FUN = "sum") ## sum all the ASV together in the same sample
rownames(phylum_tab) = phylum_tab$'TAXtab$Phylum'
phylum_tab=phylum_tab[,-1]
total_ab_phyla = rowSums(phylum_tab)/sum(rowSums(phylum_tab))
total_ab_phyla = total_ab_phyla*100
names(total_ab_phyla) = rownames(phylum_tab)
total_ab_phyla
top10 = names(sort(total_ab_phyla, decreasing = TRUE))
top10 = top10[1:10]
top10

phylum_tab$Phylum = rownames(phylum_tab)
phylum_tab_melted = melt(phylum_tab, value = Phylum) ## This melt the table

colnames(phylum_tab_melted)[1:2] = c("Phylum", "sample_name")
metadata$sample_name = as.factor(rownames(metadata))
class(phylum_tab_melted$sample_name)
phylum_tab_melted$group = metadata$Location[match(phylum_tab_melted$sample_name, metadata$sample_name)]

phylum_tab_melted1 = phylum_tab_melted[!is.na(phylum_tab_melted$group),]
phylum_tab_melted1$Phylum = ifelse(phylum_tab_melted1$Phylum %in% top10,
                                   phylum_tab_melted1$Phylum, "others")


ggplot(phylum_tab_melted1,aes(x = group, y = value, fill = Phylum))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_brewer(palette = "Paired") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Relative abundance") +
  theme_classic() + labs(title = "Taxa - Incubation",
                         x = "Location") 




phylum_tab = aggregate(t(ASVtab) ~ TAXtab$Phylum, FUN = "sum") ## sum all the ASV together in the same sample
rownames(phylum_tab) = phylum_tab$`TAXtab$Phylum`
phylum_tab=phylum_tab[,-1]
total_ab_phyla = rowSums(phylum_tab)/sum(rowSums(phylum_tab))
total_ab_phyla = total_ab_phyla*100
names(total_ab_phyla) = phylum_tab$`TAXtab$Phylum`
total_ab_phyla
top10 = names(sort(total_ab_phyla, decreasing = TRUE))
top10 = top10[1:10]
top10


phylum_tab_melted = melt(phylum_tab, value = phylum_tab$`TAXtab$Phylum`) ## This melt the table

colnames(phylum_tab_melted)[1:2] = c("Phylum", "sample_name")
phylum_tab_melted$group = metadata$Location[match(phylum_tab_melted$sample_name, metadata$sample_name)]

phylum_tab_melted1 = phylum_tab_melted[!is.na(phylum_tab_melted$group),]
phylum_tab_melted1$Phylum = ifelse(phylum_tab_melted1$Phylum %in% top10,
                                   phylum_tab_melted1$Phylum, "others")


ggplot(phylum_tab_melted1,aes(x = sample_name, y = value, fill = Phylum))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_brewer(palette = "Paired") + 
  facet_wrap(~group, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Relative abundance")


geom_bar(position="fill", stat="identity")+
  facet_wrap(~group, scales = "free") +
  scale_fill_brewer(palette = "Paired") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Relative abundance") +
  theme_classic() + labs(title = "Phylum at incubation",
                         x = "Location") 
