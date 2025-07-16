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

miseq_path <- "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/geo_seqs/geo_seqs" # CHANGE to the directory containing the fastq files after unzipping.
list.files(miseq_path)


# Sort ensures forward/reverse reads are in same order
  fnFs <- sort(list.files(miseq_path, pattern="_R1.fastq"))
  fnRs <- sort(list.files(miseq_path, pattern="_R2.fastq"))
  
  
  
# Extract sample names

  # Function to extract sample names
  extract_sample_name <- function(fnFs) {
    sapply(fnFs, function(x) {
      # Match sample name patterns including NOCA22NXXX
      match <- regmatches(x, regexpr("_(M\\d+[A-Z]*\\d+|Blank\\d+|BLANK|NOCA22N\\d{3})_", x))
      if (length(match) > 0) {
        gsub("_", "", match)  # Remove underscores
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
  head(out)
  
# Error rates
  
  errF <- learnErrors(filtFs, multithread=FALSE)
   saveRDS(errF, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/errF.rds")
  
  errR <- learnErrors(filtRs, multithread=FALSE)
   saveRDS(errR, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/errR.rds")
  
  
  plotErrors(errF, nominalQ=TRUE)
  
# Dereplication
  
  derepFs <- derepFastq(filtFs, verbose=TRUE)
    saveRDS(derepFs, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/derepFs.rds")
  
  derepRs <- derepFastq(filtRs, verbose=TRUE)
    saveRDS(derepRs, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/derepRs.rds")
  
  
  # Name the derep-class objects by the sample names
  names(derepFs) <- sample_names
  names(derepRs) <- sample_names
  
# Sample inference
  dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
    saveRDS(dadaFs, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/dadaFs.rds")
  
  dadaRs <- dada(derepRs, err=errR, multithread=FALSE)
    saveRDS(dadaRs, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/dadaRs.rds")
  
  dadaFs[[1]]
  
# Merge paired reads
  
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  # Inspect the merger data.frame from the first sample
  head(mergers[[1]])
  saveRDS(mergers, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/mergers.rds")
  
# Construct sequence table
  seqtab <- makeSequenceTable(mergers)
  dim(seqtab)
  # Inspect distribution of sequence lengths
  table(nchar(getSequences(seqtab)))
  saveRDS(seqtab, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/seqtab.rds")
  
# Remove chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
  dim(seqtab.nochim)
  sum(seqtab.nochim)/sum(seqtab)
    saveRDS(seqtab.nochim, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/seqtab_nochim.rds")
  
  
# Track reads through pipeline
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample_names
  head(track)
  
  # Export as a csv
  write.csv(track, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/dada2_sample_reads_NOCA_geo.csv")

# Assign taxonomy - download and move to folder
  taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/tax/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=FALSE)
  saveRDS(taxa, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/taxa_geo.rds")
    taxa <- readRDS("C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/taxa_geo.rds")
  taxa <- addSpecies(taxa, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/tax/silva_v138.2_assignSpecies.fa.gz")
  saveRDS(taxa, "C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/taxa_geo_species.rds")
    # taxa <- readRDS("path/to/taxa.rds")
  taxa.print <- taxa # Removing sequence rownames for display only
  rownames(taxa.print) <- NULL
  head(taxa.print)

# Merge taxa tables from different dada2  
  
  
# Metadata
  # Read in metadata file
  metadatafull <- read.csv("C:/Users/natal/OneDrive - University of Arizona/Cornell/ht_njm/ht_njm_R/metadata_ht_njm.csv") 
  
  row.names(metadatafull)<-metadatafull$sample_names
  META <- sample_data(metadatafull)
  head(META)
  
# Construct phyloseq object
  ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                 META, 
                 tax_table(taxa))
  ps
  
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
  sample_reads<-data.frame(reads=sample_sums(ps))
  head(sample_reads)
  
##### Decontaminate #####
  
  library(decontam)
  
  # Remove contaminants with the combined function of prevalence and frequency
  
  sample_data(ps)$is.neg <- sample_data(ps)$sample_or_control == "negative_control" # metadata column, blanks as negative_control, samples as true_sample
  contamdf.prev <- isContaminant(ps, method="combined",conc="quant_reading",neg="is.neg") # blanks must have number greater than 0, set to 0.0001
  table(contamdf.prev$contaminant) #table--true column = number of contaminants
  
  # Creating the decontaminated phyloseq object (physeq_nc)
  
  physeq_nc<-prune_taxa(!contamdf.prev$contaminant,ps) 
  
  physeq_nc #identified X contaminants and removed them leaving Y taxa
  
  # Remove negative controls (n=17)
  
  physeq2 <- subset_samples(physeq_nc, sample_or_control == "true_sample") 
  physeq2 #
  
  unique(physeq2@tax_table@.Data[,"Phylum"])
  unique(physeq2@tax_table@.Data[,"Species"])
  View(tax_table(physeq2))
  