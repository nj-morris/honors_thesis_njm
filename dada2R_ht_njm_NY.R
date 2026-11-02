# Natalie J Morris Honors Thesis #
# Geographic Comparison of Tree Swallow Cloacal Microbiome #
# DNA extraction: Qiagen DNeasy PowerSoil Pro #
# PCR: V4 of 16S rRNA - Earth Microbiome Protocol - 515F and 806R #
# Sequencing: Cornell Biotechnology Resource Center - Illumina MiSeq paired-end 2 x 250 bp #
# Code adapted from Callahan et al., 2017 #


# Best Practices # 
# Each sequencing run = Different dada2 run, then merge feature tables #
# Control for different PCR plates in models #

# Set up
library("knitr")
library("ggplot2")
library("gridExtra")
library("phyloseq")
library("DECIPHER")
library("dada2")
library("phangorn")


set.seed(100)

miseq_path <- "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs" # CHANGE to the directory containing the fastq files after unzipping.
list.files(miseq_path)


# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_R1.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2.fastq"))



# Extract sample names

# Function to extract sample names
extract_sample_name <- function(fnFs) {
  sapply(fnFs, function(x) {
    # Match sample names, prioritizing MB### over C#
    match <- regmatches(x, regexpr("_(MB\\d+|Blank\\d*|BLANK|NOCA22N\\d{3}|(?<!Plate\\d{1,2}_)C\\d+)_", x, perl = TRUE))
    if (length(match) > 0) {
      gsub("_", "", match)  # Remove underscores from the matched string
    } else {
      NA  # Return NA if no match is found
    }
  })
}



# Apply the function
sample_names <- extract_sample_name(fnFs)

# Verify results
print(sample_names)


# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)

fnFs[1:3] #check
fnRs[1:3]

# Inspect read quality profiles   
plotQualityProfile(fnFs[15:18])
plotQualityProfile(fnRs[1:2])

# Filter and trim
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample_names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(180,180), trimLeft = c(19,20), # remove primers #
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, verbose = TRUE) # On Windows set multithread=FALSE
  saveRDS(out, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/out.rds")

head(out)

# Error rates

errF <- learnErrors(filtFs, multithread=FALSE)
saveRDS(errF, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/errF.rds")

errR <- learnErrors(filtRs, multithread=FALSE)
saveRDS(errR, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/errR.rds")


plotErrors(errF, nominalQ=TRUE)

# Dereplication

derepFs <- derepFastq(filtFs, verbose=TRUE)
saveRDS(derepFs, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/derepFs.rds")

derepRs <- derepFastq(filtRs, verbose=TRUE)
saveRDS(derepRs, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/derepRs.rds")


# Name the derep-class objects by the sample names
names(derepFs) <- sample_names
names(derepRs) <- sample_names

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
saveRDS(dadaFs, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/dadaFs.rds")

dadaRs <- dada(derepRs, err=errR, multithread=FALSE)
saveRDS(dadaRs, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/dadaRs.rds")

dadaFs[[1]]

# Merge paired reads

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
saveRDS(mergers, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/mergers.rds")

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
saveRDS(seqtab, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/seqtab.rds")

# Read in other seq table
seqtab_geo <- readRDS("C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/seqtab.rds")

# Merge seq tables - 
mergetab <- mergeSequenceTables(seqtab, seqtab_geo)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(mergetab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
  saveRDS(seqtab.nochim, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/seqtab_nochim_merged.rds")
  seqtab.nochim <- readRDS("C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/seqtab_nochim_merged.rds")

# Remove NOCA samples
seqtab.nochim <- seqtab_nochim[!grepl("^NOCA", rownames(seqtab_nochim)), ]

# save
saveRDS(seqtab.nochim, file = "seqtab.nochim.noNOCA.RDS")

# Construct phylogenetic tree
  seqs <- getSequences(seqtab.nochim)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
  
  phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phangAlign)
    saveRDS(dm, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/seqtab_nochim_merged.rds")
  
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  fit = pml(treeNJ, data=phangAlign)
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
  saveRDS(fitGTR, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/fitGTR.rds")
  detach("package:phangorn", unload=TRUE)
  
# Track reads through pipeline - for each dada2 run
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
head(track)

# Export as a csv
write.csv(track, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/dada2_sample_reads_NY_geo.csv")

# Assign taxonomy - download and move to tax folder
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/tax/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=FALSE)
    saveRDS(taxa, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/taxa_geo_merged.rds")
    taxa <- readRDS("C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/taxa_geo_merged.rds")

taxa <- addSpecies(taxa, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/tax/silva_v138.2_assignSpecies.fa.gz")
    saveRDS(taxa, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/raw_tres_seqs/taxa_geo_specie_merged.rds")
    taxa <- readRDS("path/to/taxa.rds")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Metadata
# Read in metadata file
metadatafull <- read.csv("C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/metadata_ht_njm_2025.csv") 
metadatafull <- as.data.frame(metadata_ht_njm_2025)
row.names(metadatafull)<-metadatafull$sample_names
META <- sample_data(metadatafull)
head(META)

# Construct phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               META, 
               tax_table(taxa))
              # phy_tree(fitGTR$tree)#
ps

saveRDS(ps, file = "phyloseq_geo_njm_OG.RDS")

phyloseq_geo_njm_OG
sum(sample_sums(phyloseq_geo_njm_OG))

ps <- phyloseq_geo_njm_OG

# Show available ranks in the dataset
rank_names(ps)
# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)

### From F1000 ####
# Remove NA phyla
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
# Calculate 
library(dplyr)
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})


######

# Retain only bacteria
library(dplyr)

ps2 <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "Mitochondria" &             #filter out mitochondria
      Class   != "Chloroplast"                #filter out chloroplasts
  )
ps2

# look at phyla 
unique(ps2@tax_table@.Data[,"Phylum"])

#check the number of reads per sample
sample_readsB<-data.frame(reads=sample_sums(ps2))
head(sample_readsB)

saveRDS(ps2, "phyloseq_geo_njm_predecontam.RDS")
##### Decontaminate - from Jenn Houtz #####

library(decontam)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")
library(decontam)

# Remove contaminants with the combined function of prevalence and frequency

sample_data(ps2)$is.neg <- sample_data(ps2)$sample_or_control == "negative_control" # metadata column, blanks as negative_control, samples as true_sample
contamdf.prev <- isContaminant(ps2, method="combined",conc="quant_reading",neg="is.neg") # blanks must have number greater than 0, set to 0.0001
table(contamdf.prev$contaminant) #table--true column = number of contaminants

# look at what contaminants are:
contaminant_taxa <- rownames(contamdf.prev)[contamdf.prev$contaminant == TRUE]
contaminant_tax_table <- tax_table(ps2)[contaminant_taxa, ]
View(contaminant_tax_table)
otu <- otu_table(ps2)

if (!taxa_are_rows(otu)) {
  otu <- t(otu)
}
otu_contaminants <- otu[contaminant_taxa, ]
contaminant_counts <- rowSums(otu_contaminants > 0)
contaminant_tax_table <- as.data.frame(tax_table(ps2)[contaminant_taxa, ])
contaminant_tax_table$NumSamplesPresent <- contaminant_counts
sample_meta <- sample_data(ps2)
true_samples <- rownames(sample_meta)[sample_meta$sample_or_control == "true_sample"]
neg_controls <- rownames(sample_meta)[sample_meta$sample_or_control == "negative_control"]
# Subset OTU table columns
otu_true <- otu_contaminants[, true_samples]
otu_neg <- otu_contaminants[, neg_controls]
otu_neg <- as.data.frame(otu_neg)
otu_true <- as.data.frame(otu_true)

write_csv(otu_neg, "otu_neg_decontam.csv")
# Count number of samples where each contaminant is present
present_in_true <- rowSums(otu_true > 0)
present_in_neg <- rowSums(otu_neg > 0)
contam_summary <- data.frame(
  Taxon = rownames(otu_contaminants),
  PresentInTrueSamples = present_in_true,
  PresentInNegControls = present_in_neg
)

# Add taxonomy info (optional)
tax_info <- as.data.frame(tax_table(ps2)[contaminant_taxa, ])
otutrue <- cbind(otu_true, tax_info)

# View it
View(otu_neg)
View(contam_summary)
write_csv(contam_summary, "decontam_contaminants_ht_njm.csv")

# Creating the decontaminated phyloseq object (physeq_nc)

physeq_nc<-prune_taxa(!contamdf.prev$contaminant,ps2) 
physeq_nc #identified 56 contaminants and removed them leaving Y taxa

# Remove negative controls (n=17)

physeq2 <- subset_samples(physeq_nc, sample_or_control == "true_sample") 
physeq2 #

unique(physeq2@tax_table@.Data[,"Phylum"])
unique(physeq2@tax_table@.Data[,"Species"])
View(tax_table(physeq2))

saveRDS(physeq2, "phyloseq_ht_njm_decontam_true.rds")


#Export feature/OTU table as a biom table for the post-filtering sequence and ASV metrics
######DO NOT TRANSFORM OR THE FEATURES AND SAMPLES WILL BE SWITCHED

library(biomformat);packageVersion("biomformat")

otu<-as(otu_table(phyloseq_ht_njm_USE),"matrix")
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"otu_biom_htnjm.biom")


# Picrust
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

library("stringr")
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta,  "ASVs.fa")
#################### or ?????
library(dada2)
library(ShortRead)
writeFastq()
writeFasta(as.character(ASVseq), "rep-seqs.fna")
#############
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

OTU_TAX_table <- merge(asv_tab, asv_tax, by=0)

write.table(OTU_TAX_table, "OTU_TAX_table.tsv", sep = "\t", quote=F, col.names=NA)

############ 

otu  = read.delim("OTU_TAX_table.tsv", row.names = 2)
otu= otu[,-1]
otu = t(otu)
write.table(otu, "OTU_TAX_table_transposed_coltrue.tsv", sep = "\t", quote=F,
            col.names = NA)
otu <- read.table("OTU_TAX_table_transposed_coltrue.tsv", 
                  sep = "\t", header = TRUE, row.names = 1)
asv <- read.table("ASVs_counts.tsv")
asv = t(asv)

write.table(asv, "ASVs_counts.tsv", sep = "\t", quote=F,
            col.names = NA)

############## look at clustering 

#Extract metadata

ps2data<- data.frame(sample_data(ps2))

# Look at structure and make sure they are correct data types

str(incub_beta_data)
ps2data$sample_or_control <- factor(ps2data$sample_or_control)

##### check sample ordination

library(ggforce)
ordu = ordinate(ps2, "PCoA", "bray", weighted=FALSE)
P <- plot_ordination(ps2, ordu, color="sample_or_control")
pca <- P + geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_bw() +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

geom_mark_ellipse(aes(color = Location,
                        label=Location),
                    expand = unit(0.5,"mm")) +
  scale_color_manual(values = location_colors)

pca


############# figure out reads and ASVs
dim(seqtab.nochim.noNOCA)
sum(sample_reads)
