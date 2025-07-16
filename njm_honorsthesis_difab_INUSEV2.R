#### Dif Abundance #####
# https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html

library(tidyverse)
library(DT)
library(dplyr)
library(tidyr)
library(phyloseq)
library(rcompanion)
library(microViz)
library(microshades)
library(cowplot)
library(vegan)
library(patchwork)
library(speedyseq)
#library(NBZIMM)
library(magrittr)
library(lubridate)
library(lme4)
library(lmerTest)
library(visreg)
library(ggpubr)
library(tibble)
library(permute)
library(forcats)

library(MicEco)
library(glmmTMB)

library(devtools)
install_github("Russel88/MicEco")

pseq <- readRDS("phyloseq_ht_njm_USE.rds")
sample_data(pseq)$sample_reads <- sample_sums(pseq)

pseq <- subset_samples(pseq, Capture_Number %in% c("1", "3"))

inc <- subset_samples(pseq, Capture_Number == "1")
prov <- subset_samples(pseq, Capture_Number == "3")
pseq <- subset_samples(pseq, !Location == "Wyoming")
pseq <- subset_samples(pseq, Location == "New_York")

pseq <- subset_samples(pseq, !(Individual_Treatment %in% c("Dulled", "Low_Tape", "Predator") & Capture_Number == "3"))



sample_data(pseq)$sample_reads <- sample_sums(pseq)

# To compare
# locations and incubation/provisioning
# mass loss and composition change 


###percentages of taxa
phy_plot <-
  pseq %>%
  tax_glom("Phylum") %>%
  transform_sample_counts(function(x)100* x / sum(x)) %>%
  psmelt() %>%
  as_tibble()

phy_abd<-phy_plot %>% #abundances
  group_by(Phylum) %>%
  summarise(Abundance = mean(Abundance)) %>%
  arrange(desc(Abundance))

fam_plot <-
  pseq %>%
  tax_glom("Family") %>%
  transform_sample_counts(function(x)100* x / sum(x)) %>%
  psmelt() %>%
  as_tibble()

fam_abd<-fam_plot %>% #abundances
  group_by(Family) %>%
  summarise(Abundance = mean(Abundance)) %>%
  arrange(desc(Abundance))

gen_plot <-
  pseq %>%
  tax_glom("Genus") %>%
  transform_sample_counts(function(x)100* x / sum(x)) %>%
  psmelt() %>%
  as_tibble()

gen_abd<-gen_plot %>% #abundances
  group_by(Genus) %>%
  summarise(Abundance = mean(Abundance)) %>%
  arrange(desc(Abundance))

top_fams <- unique(fam_plot$Family)
fam_abd<-subset(fam_abd, Family %in% top_fams)
values_to_remove <- c("Other", "uncultured")
# Remove values
top_fams <- setdiff(top_fams, values_to_remove)

# Print the vector of top 10 Families
fam_plot1 <- subset(fam_plot, Family %in% top_fams)

bar1<-ggplot(data=fam_plot1, aes(x=fct_reorder(Family, Abundance, .fun = median, .desc = TRUE), 
                                 y=Abundance, fill=Family)) + geom_boxplot() + theme_bw() + 
  theme(axis.text.y = element_text(size=18), axis.title.x = element_text(size=18),
        axis.text.x = element_text(size=18), plot.tag = element_text(size=22),
        legend.position = "none") + coord_flip() + 
  labs(y="Relative abundance", x="", tag="B")
bar1

##########
###########DA MODELS ### FAM LEVEL WITH GLMMTMB
#################################### #################################### #################################### #################################### 
library(glmmTMB)
pseq #15879 taxa and 376 samples

### subset to inc or prov
pseq <- subset_samples(pseq, Capture_Number == "1")


#### Phylum Level
htp = tax_glom(pseq, taxrank = "Phylum") #glom , now 35
htp #35 
# Convert counts to relative abundance
phy_rel <- transform_sample_counts(htp, function(x) x / sum(x))
# Prune taxa with mean relative abundance >= 0.005 (i.e., 0.5%)
abundant_taxa <- taxa_names(filter_taxa(phy_rel, function(x) mean(x) >= 0.01, TRUE))
# Filter by prevalence (e.g., present in >= 2 samples with at least 10 reads)
phy_pruned <- ps_prune(htp, min.samples = 2, min.reads = 10)
prevalent_taxa <- taxa_names(phy_pruned)
# Keep only taxa that are both abundant and prevalent
keep_taxa <- intersect(abundant_taxa, prevalent_taxa)
# Final subsetted phyloseq object
pruned <- prune_taxa(keep_taxa, htp)
pruned # 14 - 6 at 1%

htp_df <- psmelt(pruned) 
names(htp_df)

htp_df1 <- htp_df[,c("Abundance", "sample_names", "Individual_Band", "Individual_Treatment", "Location",
                     "Site", "Capture_Number", "Age", "Mass", "Bill.Head", "Flat_Wing", "Fledged", "Fledged_Num", "CIdate", 
                     "sample_reads", "Phylum")]
names(htp_df1) #keep columns needed
dim(htp_df1)
htp_df2=reshape(data=htp_df1, idvar = c("sample_names", "Individual_Band", "Individual_Treatment", "Location",
                                        "Site", "Capture_Number", "Age", "Mass", "Bill.Head", "Flat_Wing", "Fledged", "Fledged_Num", "CIdate", 
                                        "sample_reads", "Phylum"), timevar = "Phylum", direction = "wide")
names(htp_df2) #which cols are taxa - cols 14:70

#clean up column names
colnames(htp_df2)[15:20] <- gsub("^Abundance\\.", "", colnames(htp_df2)[15:20])
names(htp_df2) #check, good

#now standardize/scale continuous variables for modelling in the loop (like density)
# stdize=function(x) {(x - mean(x))/sd(x)}
# mat_fam_df2$local_density_50_scaled<-stdize(terr_fam_df2$local_density_50) #for comparison of estimates

# IF FAILS:
#remove some taxa that don't converge by name and run these later in another way
View(htp_df2)
mat_fam_df3 <- mat_fam_df2  %>% dplyr::select(-Cyanobiaceae, -Marine_Group_II, -Clade_I, -Actinomarinaceae, -SAR86_clade,-`UCG-010`, -Saprospiraceae, -`Marinimicrobia_(SAR406_clade)`, -`AEGEAN-169_marine_group`) # -Cryomorphaceae,Nitrosopumilaceae,Sanguibacteraceae)
names(mat_fam_df3) #check, good

htp_df3 <- htp_df2


####
ht_loop <- NULL #create NULL df to fill
for (i in 15:20){
  htp_df3$taxa<-as.numeric(htp_df3[,i]) #for i in which subset of taxa (columns) to iterate through
  mmb<-glmmTMB(taxa ~ Location + Mass + offset(log(sample_reads))  + (1|Site) , 
               data=htp_df3,
               #ziformula = ~1, #zero - inflated, depends on level aglomm
               control=glmmTMBControl(optCtrl=list(iter.max=1000)),
               family = nbinom1) #negative binomial, nbinom1 runs better than nbinom2
 # Shapiro-Wilk test on residuals
  shap <- as.numeric(shapiro.test(resid(mmb))$p.value)
  # Extract coefficients & p-values dynamically
  rep_rows <- grep("Location", rownames(summary(mmb)$coefficients$cond))
  coef_rep <- round(summary(mmb)$coefficients$cond[rep_rows, 1], 4)  # Estimates
  pval_rep <- round(summary(mmb)$coefficients$cond[rep_rows, 4], 4)  # P-values
  # Create a dataframe for storage
  df <- data.frame(Column = as.character(i), 
                   Phylum = colnames(htp_df3)[i], 
                   Shapiro = shap, 
                   est_rep = coef_rep, 
                   p_rep = pval_rep)
  # Append results
  ht_loop <- rbind(ht_loop, df)
}

ht_loop<-data.frame(ht_loop)
View(ht_loop) 

#convert class
ht_loop$Shapiro<-as.numeric(as.character(ht_loop$Shapiro))
ht_loop$est_rep<-as.numeric(as.character(ht_loop$est_rep))
ht_loop$p_rep<-as.numeric(as.character(ht_loop$p_rep))

#Adjust pvalues for multiple testing
ht_loop$p_rep_fdr<-round(p.adjust(ht_loop$p_rep, method="fdr"),8)

#Subset to significant taxa only
ht_loop_sig=subset(ht_loop, p_rep_fdr <= 0.05) 
dim(ht_loop_sig)#
View(ht_loop_sig) #
write.csv(ht_loop,"ht_loc_both_phylum_DA.csv") ### SAVE

# no FDR significant differences across capture number - for all locations clumped
# yes FDR sig across locations for incubation ! but is alaska reference ??? 


############# family level <3

htf = tax_glom(pseq, taxrank = "Family") #glom , now 35
htf # 408
# Convert counts to relative abundance
f_rel <- transform_sample_counts(htf, function(x) x / sum(x))
# Prune taxa with mean relative abundance >= 0.005 (i.e., 0.5%)
abundant_taxa <- taxa_names(filter_taxa(f_rel, function(x) mean(x) >= 0.001, TRUE))
# Filter by prevalence (e.g., present in >= 2 samples with at least 10 reads)
f_pruned <- ps_prune(htf, min.samples = 2, min.reads = 10)
prevalent_taxa <- taxa_names(f_pruned)
# Keep only taxa that are both abundant and prevalent
keep_taxa <- intersect(abundant_taxa, prevalent_taxa)
# Final subsetted phyloseq object
pruned <- prune_taxa(keep_taxa, htf)
pruned # 59

htf_df <- psmelt(pruned) 
names(htf_df)

htf_df1 <- htf_df[,c("Abundance", "sample_names", "Individual_Band", "Individual_Treatment", "Location",
                     "Site", "Capture_Number", "Age", "Mass", "Bill.Head", "Flat_Wing", "Fledged", "Fledged_Num", "CIdate", 
                     "sample_reads", "Family")]
names(htf_df1) #keep columns needed
dim(htf_df1)
htf_df2=reshape(data=htf_df1, idvar = c("sample_names", "Individual_Band", "Individual_Treatment", "Location",
                                        "Site", "Capture_Number", "Age", "Mass", "Bill.Head", "Flat_Wing", "Fledged", "Fledged_Num", "CIdate", 
                                        "sample_reads", "Family"), timevar = "Family", direction = "wide")
names(htf_df2) #which cols are taxa - cols 14:70

#clean up column names
colnames(htf_df2)[15:65] <- gsub("^Abundance\\.", "", colnames(htf_df2)[15:65])
names(htf_df2) #check, good


htf_df3 <- htf_df2


####
ht_loop <- NULL #create NULL df to fill
for (i in 15:65){
  htf_df3$taxa<-as.numeric(htf_df3[,i]) #for i in which subset of taxa (columns) to iterate through
  mmb<-glmmTMB(taxa ~ Capture_Number + Location  + Mass + offset(log(sample_reads)) + (1|Individual_Band) + (1|Site), 
               data=htf_df3,
               #ziformula = ~1, #zero - inflated, depends on level aglomm
               control=glmmTMBControl(optCtrl=list(iter.max=1000)),
               family = nbinom1) #negative binomial, nbinom1 runs better than nbinom2
  # Shapiro-Wilk test on residuals
  shap <- as.numeric(shapiro.test(resid(mmb))$p.value)
  # Extract coefficients & p-values dynamically
  rep_rows <- grep("Capture_Number", rownames(summary(mmb)$coefficients$cond))
  coef_rep <- round(summary(mmb)$coefficients$cond[rep_rows, 1], 4)  # Estimates
  pval_rep <- round(summary(mmb)$coefficients$cond[rep_rows, 4], 4)  # P-values
  # Create a dataframe for storage
  df <- data.frame(Column = as.character(i), 
                   Family = colnames(htf_df3)[i], ####### change!!!!!!!!!!!!
                   Shapiro = shap, 
                   est_rep = coef_rep, 
                   p_rep = pval_rep)
  # Append results
  ht_loop <- rbind(ht_loop, df)
}

ht_loop<-data.frame(ht_loop)
View(ht_loop) 

#convert class
ht_loop$Shapiro<-as.numeric(as.character(ht_loop$Shapiro))
ht_loop$est_rep<-as.numeric(as.character(ht_loop$est_rep))
ht_loop$p_rep<-as.numeric(as.character(ht_loop$p_rep))

#Adjust pvalues for multiple testing
ht_loop$p_rep_fdr<-round(p.adjust(ht_loop$p_rep, method="fdr"),5)

#Subset to significant taxa only
ht_loop_sig=subset(ht_loop, p_rep_fdr <= 0.05) 
dim(ht_loop_sig)#
View(ht_loop_sig) #
write.csv(ht_loop_sig,"ht_inctoprov_family_DA_noWY.csv") ### SAVE
######################### genus level


#### genus Level
htg = tax_glom(pseq, taxrank = "Genus") #glom , now 35
htg # 1399
# Convert counts to relative abundance
g_rel <- transform_sample_counts(htg, function(x) x / sum(x))
# Prune taxa with mean relative abundance 
abundant_taxa <- taxa_names(filter_taxa(g_rel, function(x) mean(x) >= 0.01, TRUE))
# Filter by prevalence (e.g., present in >= 2 samples with at least 10 reads)
g_pruned <- ps_prune(htg, min.samples = 2, min.reads = 10)
prevalent_taxa <- taxa_names(g_pruned)
# Keep only taxa that are both abundant and prevalent
keep_taxa <- intersect(abundant_taxa, prevalent_taxa)
# Final subsetted phyloseq object
pruned <- prune_taxa(keep_taxa, htg)
pruned # 129

htg_df <- psmelt(pruned) 
names(htg_df)

htg_df1 <- htg_df[,c("Abundance", "sample_names", "Individual_Band", "Individual_Treatment", "Location",
                     "Site", "Capture_Number", "Age", "Mass", "Bill.Head", "Flat_Wing", "Fledged", "Fledged_Num", "CIdate", 
                     "sample_reads", "Genus")]
names(htg_df1) #keep columns needed
dim(htg_df1)
htg_df2=reshape(data=htg_df1, idvar = c("sample_names", "Individual_Band", "Individual_Treatment", "Location",
                                        "Site", "Capture_Number", "Age", "Mass", "Bill.Head", "Flat_Wing", "Fledged", "Fledged_Num", "CIdate", 
                                        "sample_reads", "Genus"), timevar = "Genus", direction = "wide")
names(htg_df2) #which cols are taxa - cols 14:70

#clean up column names
colnames(htg_df2)[15:30] <- gsub("^Abundance\\.", "", colnames(htg_df2)[15:30])
names(htg_df2) #check, good

#now standardize/scale continuous variables for modelling in the loop (like density)
# stdize=function(x) {(x - mean(x))/sd(x)}
# mat_fam_df2$local_density_50_scaled<-stdize(terr_fam_df2$local_density_50) #for comparison of estimates

# IF FAILS:
#remove some taxa that don't converge by name and run these later in another way
View(htg_df2)
mat_fam_df3 <- mat_fam_df2  %>% dplyr::select(-Cyanobiaceae, -Marine_Group_II, -Clade_I, -Actinomarinaceae, -SAR86_clade,-`UCG-010`, -Saprospiraceae, -`Marinimicrobia_(SAR406_clade)`, -`AEGEAN-169_marine_group`) # -Cryomorphaceae,Nitrosopumilaceae,Sanguibacteraceae)
names(mat_fam_df3) #check, good

htg_df3 <- htg_df2


####
ht_loop <- NULL #create NULL df to fill
for (i in 15:30){
  htg_df3$taxa<-as.numeric(htg_df3[,i]) #for i in which subset of taxa (columns) to iterate through
  mmb<-glmmTMB(taxa ~ Capture_Number + Location + Mass + offset(log(sample_reads)) + (1|Individual_Band) + (1|Site), 
               data=htg_df3,
               #ziformula = ~1, #zero - inflated, depends on level aglomm
               control=glmmTMBControl(optCtrl=list(iter.max=1000)),
               family = nbinom1) #negative binomial, nbinom1 runs better than nbinom2
  # Shapiro-Wilk test on residuals
  shap <- as.numeric(shapiro.test(resid(mmb))$p.value)
  # Extract coefficients & p-values dynamically
  rep_rows <- grep("Capture_Number", rownames(summary(mmb)$coefficients$cond))
  coef_rep <- round(summary(mmb)$coefficients$cond[rep_rows, 1], 4)  # Estimates
  pval_rep <- round(summary(mmb)$coefficients$cond[rep_rows, 4], 4)  # P-values
  # Create a dataframe for storage
  df <- data.frame(Column = as.character(i), 
                   Genus = colnames(htg_df3)[i], ####### change!!!!!!!!!!!!
                   Shapiro = shap, 
                   est_rep = coef_rep, 
                   p_rep = pval_rep)
  # Append results
  ht_loop <- rbind(ht_loop, df)
}
mmb<-glmmTMB(Lactobacillus ~ Location  + Mass + offset(log(sample_reads)) + (1|Individual_Band) + (1|Site), 
               data=htg_df3,
               #ziformula = ~1, #zero - inflated, depends on level aglomm
               control=glmmTMBControl(optCtrl=list(iter.max=1000)),
               family = nbinom2)
summary(mmb)

ht_loop<-data.frame(ht_loop)
View(ht_loop) 

#convert class
ht_loop$Shapiro<-as.numeric(as.character(ht_loop$Shapiro))
ht_loop$est_rep<-as.numeric(as.character(ht_loop$est_rep))
ht_loop$p_rep<-as.numeric(as.character(ht_loop$p_rep))

#Adjust pvalues for multiple testing
ht_loop$p_rep_fdr<-round(p.adjust(ht_loop$p_rep, method="fdr"),5)

#Subset to significant taxa only
ht_loop_sig=subset(ht_loop, p_rep_fdr <= 0.05) 
dim(ht_loop_sig)#
View(ht_loop_sig) #
write.csv(ht_loop_sig,"ht_inctoprov_genus_DA_noWY.csv") ### SAVE


####### plot ##########

# phylaaa
taxa_cols <- htp_df2[, 15:20]  
long_taxa <- melt(taxa_cols)

ggplot(long_taxa, aes(x = value)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free") +  # Separate histograms per taxa
  theme_minimal() +
  labs(title = "Phyla:Histograms of Taxa Abundances", x = "Taxa Count", y = "Frequency")

# family
taxa_cols <- mat_fam_df2[, 17:34]  
long_taxa <- melt(taxa_cols)

ggplot(long_taxa, aes(x = value)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free") +  # Separate histograms per taxa
  theme_minimal() +
  labs(title = "Family: Histograms of Taxa Abundances", x = "Taxa Count", y = "Frequency")

# genus

taxa_cols <- htg_df3[, 41:90]  
long_taxa <- melt(taxa_cols)

ggplot(long_taxa, aes(x = value)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free") +  # Separate histograms per taxa
  theme_minimal() +
  labs(title = "Histograms of Taxa Abundances", x = "Taxa Count", y = "Frequency")





















m_fam = tax_glom(pseq, taxrank = "Family") #glom by family, now 265 families
m_fam #408 families
fam_abd_1 <- subset(fam_abd, Abundance >= 0.005) #from way above / plot
keep_tax <- fam_abd_1$Family
pruned <- subset_taxa(m_fam, Family %in% keep_tax)
pruned # 70 families with greater than avg 0.05% rel abundance 

################################## try taxa ~ rep status first before pheno ~ taxa
########### try days since part 
##### maybe try binomial for zero inflated NAN
mat_fam_df <- psmelt(pruned) 
names(mat_fam_df)

mat_fam_df1 <- mat_fam_df[,c("Abundance", "gr.x", "poop_id", "sex.x", "year.x", "readDepth.x", "squirrel_id.x", "season", "run.y",
                             "spr_density", "Family", "mast.x", "rep_status2", "days_since_part2", "litter_size")]
names(mat_fam_df1) #keep columns needed
dim(mat_fam_df1)
mat_fam_df2=reshape(data=mat_fam_df1, idvar = c("gr.x", "poop_id", "sex.x", "year.x", "readDepth.x", "squirrel_id.x", "season", "run.y",
                                                "spr_density", "Family", "mast.x", "rep_status2", "days_since_part2", "litter_size"), timevar = "Family", direction = "wide")
names(mat_fam_df2) #which cols are taxa - cols 14:70

#clean up column names
colnames(mat_fam_df2)[14:70] <- gsub("^Abundance\\.", "", colnames(mat_fam_df2)[14:70])
names(mat_fam_df2) #check, good

#now standardize/scale continuous variables for modelling in the loop (like density)
# stdize=function(x) {(x - mean(x))/sd(x)}
# mat_fam_df2$local_density_50_scaled<-stdize(terr_fam_df2$local_density_50) #for comparison of estimates

# IF FAILS:
#remove some taxa that don't converge by name and run these later in another way
View(mat_fam_df2)
mat_fam_df3 <- mat_fam_df2  %>% dplyr::select(-Cyanobiaceae, -Marine_Group_II, -Clade_I, -Actinomarinaceae, -SAR86_clade,-`UCG-010`, -Saprospiraceae, -`Marinimicrobia_(SAR406_clade)`, -`AEGEAN-169_marine_group`) # -Cryomorphaceae,Nitrosopumilaceae,Sanguibacteraceae)
names(mat_fam_df3) #check, good

mat_fam_df3 <- mat_fam_df2


####
# mat_loop <- NULL #create NULL df to fill
for (i in 14:68){
  mat_fam_df3$taxa<-as.numeric(mat_fam_df3[,i]) #for i in which subset of taxa (columns) to iterate through
  mmb<-glmmTMB(taxa ~ rep_status + mast.x + season + run.y + gr.x + offset(log(readDepth.x)) + (1|squirrel_id.x) + (1|year.x), 
               data=mat_fam_df3,
               #ziformula = ~1, #zero - inflated, depends on level aglomm
               control=glmmTMBControl(optCtrl=list(iter.max=1000)),
               family = nbinom1) #negative binomial, nbinom1 runs better than nbinom2
  shap<-as.numeric(shapiro.test(resid(mmb))[2]) #shapiro test
  coef_rep<-round(summary(mmb)$coefficients$cond[2],4) # !!!!corresponds to summary table !!!!!!!!! change!!!!!
  pval_rep<-round(summary(mmb)$coefficients$cond[26],4)
  df<-c(as.character(i),colnames(mat_fam_df3)[i],shap,coef_rep, pval_rep)
  mat_loop<-rbind(mat_loop,df)
  colnames(mat_loop)<-c("Column","Family","Shapiro","est_rep", "p_rep")
}
#summary(mmb)
# some may fail!!!! if super abundant or binomial ----- take taxa out! deal with after

# TAXA ~ REP STATUS - FAM 
mat_loop <- NULL #create NULL df to fill
for (i in 14:70) {  
  mat_fam_df3$taxa <- as.numeric(mat_fam_df3[, i])  # Convert taxa column to numeric
  # Fit the GLMM
  mmbF <- glmmTMB(taxa ~ rep_status + mast.x + season + run.y + gr.x + 
                    offset(log(readDepth.x)) + (1|squirrel_id.x) + (1|year.x), 
                  data = mat_fam_df3,
                  control = glmmTMBControl(optCtrl = list(iter.max = 1000)),
                  family = nbinom2)  # Negative binomial
  # Shapiro-Wilk test on residuals
  shap <- as.numeric(shapiro.test(resid(mmbF))$p.value)
  # Extract coefficients & p-values dynamically
  rep_rows <- grep("rep_status", rownames(summary(mmbF)$coefficients$cond))
  coef_rep <- round(summary(mmbF)$coefficients$cond[rep_rows, 1], 4)  # Estimates
  pval_rep <- round(summary(mmbF)$coefficients$cond[rep_rows, 4], 4)  # P-values
  # Create a dataframe for storage
  df <- data.frame(Column = as.character(i), 
                   Genus = colnames(mat_fam_df3)[i], 
                   Shapiro = shap, 
                   est_rep = coef_rep, 
                   p_rep = pval_rep)
  # Append results
  mat_loop <- rbind(mat_loop, df)
}

mat_loop<-data.frame(mat_loop)
View(mat_loop) 


#convert to DF
mat_loop<-data.frame(mat_loop) #X NaNs if running zero-inflated, X NaNs without zero-inflation. better.
dim(mat_loop)  #242 !! 
View(mat_loop)

#convert class
mat_loop$Shapiro<-as.numeric(as.character(mat_loop$Shapiro))
mat_loop$est_rep<-as.numeric(as.character(mat_loop$est_rep))
mat_loop$p_rep<-as.numeric(as.character(mat_loop$p_rep))

#Adjust pvalues for multiple testing
mat_loop$p_rep_fdr<-round(p.adjust(mat_loop$p_rep, method="fdr"),8)

#Subset to significant taxa only
mat_loop_sig=subset(mat_loop, p_rep_fdr <= 0.05) 
dim(mat_loop_sig)#
View(mat_loop_sig) #
##write.csv(mat_loop_sig,"mat_rep_family_DA.csv") ### SAVE



###NOW GENUS LEVEL
rs_phylo# 420 samples, 4563 ASVs (good, this is befor aegglomerateing)
mat_gen = tax_glom(rs_phylo, taxrank = "Genus") #glom by genus, 
mat_gen #453 genera 
gen_abd_1 <- subset(gen_abd, Abundance >= 0.005) #from way above / plot
keep_tax <- gen_abd_1$Genus
pruned <- subset_taxa(mat_gen, Genus %in% keep_tax)
pruned #158 genera with greater than avg 0.05% rel abundance, 147 at 0.05% abd

mat_gen_df <- psmelt(pruned) 
names(mat_gen_df)

mat_gen_df1 <- mat_gen_df[,c("Abundance", "gr.x", "poop_id", "sex.x", "year.x", "readDepth.x", "squirrel_id.x", "season", "run.y",
                             "spr_density", "Genus", "mast.x", "rep_status2", "days_since_part2", "litter_size")]
names(mat_gen_df1) #keep columns needed
dim(mat_gen_df1)
mat_gen_df2=reshape(data=mat_gen_df1, idvar = c("gr.x", "poop_id", "sex.x", "year.x", "readDepth.x", "squirrel_id.x", "season", "run.y",
                                                "spr_density", "mast.x", "rep_status2", "days_since_part2", "litter_size"), timevar = "Genus", direction = "wide")
dim(mat_gen_df2) #113 genera? probably remove uncultured, and unkonwn family
names(mat_gen_df2) #which cols are taxa - cols 14-113

#clean up column names
colnames(mat_gen_df2)[14:113] <- gsub("^Abundance\\.", "", colnames(mat_gen_df2)[14:113])
names(mat_gen_df2) #check, good

#now standardize/scale continuous variables for modelling in the loop
#stdize=function(x) {(x - mean(x))/sd(x)}
#terr_fam_df2$local_density_50_scaled<-stdize(terr_fam_df2$local_density_50) #for comparison of estimates

#remove taxa that don't converge by name and run these later in another way
mat_gen_df3 <- mat_gen_df2 %>% dplyr::select(-Synechococcus_CC9902, -`Cyanobium_PCC-6307`, -Marine_Group_II, -Candidatus_Actinomarina, -SAR86_clade, -NS5_marine_group, - Clade_Ia, -Aureispira, -`Marinimicrobia_(SAR406_clade)`, -`AEGEAN-169_marine_group`)
names(mat_gen_df3) #check, good

####)
#library(glmmTMB)
matg_loop <- NULL #create NULL df to fill
for (i in 14:107){ # adjust as removing taxa!! # 
  mat_gen_df3$taxa<-as.numeric(mat_gen_df3[,i]) #for i in which subset of taxa (columns) to iterate through
  mmbG<-glmmTMB(taxa ~ rep_status + mast.x + season + run.y + gr.x + offset(log(readDepth.x)) + (1|squirrel_id.x) + (1|year.x), 
                data=mat_gen_df3,
                #ziformula = ~1, #zero - inflation?
                control=glmmTMBControl(optCtrl=list(iter.max=1000)),
                family = nbinom1) #negative binomial, nbinom1 runs better than nbinom2
  shap<-as.numeric(shapiro.test(resid(mmbG))[2]) #shapiro test
  coef_rep<-round(summary(mmbG)$coefficients$cond[2],4)
  pval_rep<-round(summary(mmbG)$coefficients$cond[26],4)
  df<-c(as.character(i),colnames(mat_gen_df3)[i],shap,coef_rep, pval_rep)
  matg_loop<-rbind(matg_loop,df)
  colnames(matg_loop)<-c("Column","Genus","Shapiro","est_rep", "p_rep")
}

#convert to DF
matg_loop<-data.frame(matg_loop)
View(matg_loop) # NaNs here, run as binomials?

matg_loop <- NULL # Create NULL df to fill
# Taxa ~ rep_status: loop for DA
for (i in 14:103) {  
  mat_gen_df3$taxa <- as.numeric(mat_gen_df3[, i])  # Convert taxa column to numeric
  # Fit the GLMM
  mmbG <- glmmTMB(taxa ~ rep_status2 + mast.x + season + run.y + gr.x + 
                    offset(log(readDepth.x)) + (1|squirrel_id.x) + (1|year.x), 
                  data = mat_gen_df3,
                  control = glmmTMBControl(optCtrl = list(iter.max = 1000)),
                  family = nbinom2)  # Negative binomial
  # Shapiro-Wilk test on residuals
  shap <- as.numeric(shapiro.test(resid(mmbG))$p.value)
  # Extract coefficients & p-values dynamically
  rep_rows <- grep("rep_status2", rownames(summary(mmbG)$coefficients$cond))
  coef_rep <- round(summary(mmbG)$coefficients$cond[rep_rows, 1], 4)  # Estimates
  pval_rep <- round(summary(mmbG)$coefficients$cond[rep_rows, 4], 4)  # P-values
  # Create a dataframe for storage
  df <- data.frame(Column = as.character(i), 
                   Genus = colnames(mat_gen_df3)[i], 
                   Shapiro = shap, 
                   est_rep = coef_rep, 
                   p_rep = pval_rep)
  # Append results
  matg_loop <- rbind(matg_loop, df)
}

matg_loop<-data.frame(matg_loop)
View(matg_loop) 

# Ensure column names are set correctly
colnames(matg_loop) <- c("Column", "Genus", "Shapiro", "est_rep", "p_rep")

#convert class
matg_loop$Shapiro<-as.numeric(as.character(matg_loop$Shapiro))
matg_loop$est_rep<-as.numeric(as.character(matg_loop$est_rep))
matg_loop$p_rep<-as.numeric(as.character(matg_loop$p_rep))

#Adjust pvalues for multiple testing
matg_loop$p_rep_fdr<-round(p.adjust(matg_loop$p_rep, method="fdr"),8)

#Subset to significant taxa only
matg_loop_sig1=subset(matg_loop, p_rep_fdr <= 0.05) 
dim(matg_loop_sig1)# 
View(matg_loop_sig1) #
#write.csv(fly_loop_sig_age,"flying_family_DA.csv") ### SAVE

#### DA: Dealing with NAs ####
library(lme4) 
library(lmerTest)
#NaNs produced for things that are probably not zero-inflated? which
# need to figure out how to change for rep status
nans <- matg_loop %>%
  filter(is.nan(p_rep_fdr))# %>%
# distinct(Genus) 

#now rerun these without zero-inflation?
mat_gen_df4 <- mat_gen_df3 %>%
  dplyr::select(rep_status2, mast.x, season, run.y, gr.x, readDepth.x, squirrel_id.x, year.x, all_of(nans$Genus))
names(mat_gen_df4)
#check for distribution. these are all really binomial histograms
library(ggplot2)

# Select only the taxa columns
taxa_cols <- mat_gen_df4[, 9:30]  

# Convert to long format for easy ggplot handling
library(reshape2)
long_taxa <- melt(taxa_cols)

# Plot histograms in facets
ggplot(long_taxa, aes(x = value)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free") +  # Separate histograms per taxa
  theme_minimal() +
  labs(title = "Histograms of Taxa Abundances", x = "Taxa Count", y = "Frequency")


#Convert Abundance to binary variable
mat_gen_df4[, 9:28] <- lapply(mat_gen_df4[, 9:28], function(x) ifelse(x > 0, 1, 0))
mat_gen_df4[, 9:28] <- lapply(mat_gen_df4[, 9:28], as.factor)

names(mat_gen_df4)
matgNAN_loop <- NULL #create NULL df to fill
for (i in 9:28){ 
  mat_gen_df4$taxa<-mat_gen_df4[,i] #for i in which subset of taxa (columns) to iterate through
  matgNAN<-glmer(taxa ~ rep_status + mast.x + season + run.y + gr.x + 
                   offset(log(readDepth.x)) + (1|squirrel_id.x) + (1|year.x), 
                 data=mat_gen_df4,family=binomial, glmerControl(optimizer = "bobyqa",  
                                                                optCtrl = list(maxfun = 10000)))
  shap<-as.numeric(shapiro.test(resid(matgNAN))[2]) #shapiro test
  # Extract coefficients & p-values dynamically
  rep_rows <- grep("rep_status", rownames(summary(matgNAN)$coefficients$cond))
  coef_rep <- round(summary(matgNAN)$coefficients$cond[rep_rows, 1], 4)  # Estimates
  pval_rep <- round(summary(matgNAN)$coefficients$cond[rep_rows, 4], 4)  
  df<-c(as.character(i),colnames(mat_gen_df4)[i],shap,coef_rep, pval_rep)
  matgNAN_loop<-rbind(matgNAN_loop,df)
  colnames(matgNAN_loop)<-c("Column","Genus","Shapiro","est_rep", "p_rep")
}

fly<-glmer(Clade_Ia ~ local_density_50_scaled + mast.x + season + run + gr.x + 
             offset(log(readDepth)) + (1|squirrel_id.x) + (1|year.x), data=terr_fam_df4,family=binomial)
summary(fly)
#convert to DF
matgNAN_loop<-data.frame(matgNAN_loop)
dim(matgNAN_loop)  #2?
View(matgNAN_loop)

#convert class
terr_loop$Shapiro<-as.numeric(as.character(terr_loop$Shapiro))
terr_loop$est_dens<-as.numeric(as.character(terr_loop$est_dens))
terr_loop$p_dens<-as.numeric(as.character(terr_loop$p_dens))

#Adjust pvalues for multiple testing
terr_loop$p_dens_fdr<-round(p.adjust(terr_loop$p_dens, method="fdr"),8)

#Subset to significant taxa only
terr_loop_sig2=subset(terr_loop, p_dens_fdr <= 0.05) 
dim(terr_loop_sig2)#
#View(terr_loop_sig)

#combine both DFs
sigloop<-rbind(terr_loop_sig1, terr_loop_sig2)
names(sigloop)
#sigloop<-terr_loop_sig1
#write.csv(sigloop, "DA_localdens.csv") #heres the data

#read in data
sigloop<-read.csv("DA_localdens.csv")

#Significant estimates with density
sigloop$color<-ifelse(sigloop$est_dens<0,"greater at low densities","greater at high densities")
sigloop$est_dens<-as.numeric(sigloop$est_dens)
sigloop$Genus<-factor(reorder(sigloop$Genus, -sigloop$est_rep))

#add family names to the weird genera
levels(sigloop$Genus)[levels(sigloop$Genus) == "UCG-005"] <- "Oscillospiraceae_UCG-005"
levels(sigloop$Genus)[levels(sigloop$Genus) == "NK4A214_group"] <- "Oscillospiraceae_NK4A214_group"
levels(sigloop$Genus)[levels(sigloop$Genus) == "[Eubacterium]_nodatum_group"] <- "Eubacterium_nodatum_group"

#########classifiction phenotypic, 
#Prevotellaceae_Ga6A1_group --- ob anaerobic, NS
#Prevotellaceae_UCG-001 ---- ob anaerobic, NS
#UCG-005 -- (oscillospiraceae) ob anaerobic, unknown sporulation
#NK4A214_group ---- oscillospiraceae, ob anaerobic, unkonwn sporulation
#Gastranaerophilales --- obligate anaerobic, unknown sporulation
#Erysipelatoclostridium --- ob anaerobic, NS
#Frisingicoccus ------ob anaerobic, unkonwn S
#Eubacterium_nodatum_group ---- ob anaerobic, NS
#Herbinix ------ ob anaerobic, NS
#Sanguibacter ----- fac anaerobic, NS

#
un_spore<-c("Frisingicoccus", "Gastranaerophilales", "Oscillospiraceae_UCG-005_group", "Prevotellaceae_UCG-001") #nonsporulating
obs<-c("Herbinix", "Eubacterium_nodatum_group", "Gastranaerophilales", "Oscillospiraceae_NK4A214_group", 
       "Oscillospiraceae_UCG-005_group", "Prevotellaceae_UCG-001", "Prevotellaceae_Ga6A1_group", "Erysipelatoclostridium") #obligates

#Significant estimates with density
sigloop$aero<-as.factor(ifelse(sigloop$Genus %in% obs, "obligate anaerobe", "facultative anaerobe"))
sigloop$spore<-as.factor(ifelse(sigloop$Genus %in% un_spore, "unknown sporulation", "non-sporulating"))



library(scales)
library(ggpattern) #use this  fill=aero, pattern=spore

ht_loop_sig <- ht_loop_sig %>%
  mutate(Genus = fct_reorder(Genus, est_rep, .desc = FALSE))

DA_genus<-ggplot(ht_loop_sig, aes(x=Genus, y=est_rep)) + 
  #geom_bar(stat = "identity", width = .7) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025) +
  theme_bw() + coord_flip() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=10, color="black"),
        axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=12, color="black"),
        panel.grid = element_blank(),
        plot.tag = element_text(size=20)) +
  scale_fill_manual(values=c("aquamarine4","white")) +
  scale_pattern_manual(values = c("stripe", "none")) +
  labs(y="Effect size (Incubation to Provisioning)", x = "Genus", tag = "A",
       title = "Genus: Inc to Prov") +
  guides(pattern = guide_legend(override.aes = list(fill = "white"),nrow=2),
         fill = guide_legend(override.aes = list(pattern = "none"),nrow=2)) #+ 
#theme(aspect.ratio = 4/3.2)
DA_genus 


ht_loop_sig <- ht_loop_sig %>%
  mutate(Family = fct_reorder(Family, est_rep, .desc = FALSE))

DA_fam<-ggplot(ht_loop_sig, aes(x=Family, y=est_rep)) + 
  #geom_bar(stat = "identity", width = .7) + 
  geom_bar_pattern(stat = "identity", position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025) +
  theme_bw() + coord_flip() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=10, color="black"),
        axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=12, color="black"),
        panel.grid = element_blank(),
        plot.tag = element_text(size=20)) +
  scale_fill_manual(values=c("aquamarine4","white")) +
  scale_pattern_manual(values = c("stripe", "none")) +
  labs(y="Effect size (Incubation to Provisioning)", x = "Family", tag = "A", 
       title = "Family: Incubation to Provisioning") +
  guides(pattern = guide_legend(override.aes = list(fill = "white"),nrow=2),
         fill = guide_legend(override.aes = list(pattern = "none"),nrow=2)) #+ 
#theme(aspect.ratio = 4/3.2)
DA_fam



########### ANCOM BC ################
out = ancombc(data = inc, tax_level = "Phylum", 
              formula = "Location", 
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
              group = "Location", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
              n_cl = 1, verbose = TRUE)

res = out$res
res_global = out$res_global

# LFC
tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "New York - Alaska", "Tennessee - Alaska", "Wyoming - Alaska" 
)
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

samp_frac = out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out$feature_table + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
# Show the first 6 samples
round(log_corr_abn[, 1:6], 2) %>% 
  datatable(caption = "Bias-corrected log observed abundances")


df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

df_fig_age = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, age, ageSE) %>%
  dplyr::filter(age != 0) %>% 
  dplyr::arrange(desc(age)) %>%
  dplyr::mutate(direct = ifelse(age > 0, "Positive LFC", "Negative LFC"))
df_fig_age$taxon_id = factor(df_fig_age$taxon_id, levels = df_fig_age$taxon_id)
df_fig_age$direct = factor(df_fig_age$direct, 
                           levels = c("Positive LFC", "Negative LFC"))

p_age = ggplot(data = df_fig_age, 
               aes(x = taxon_id, y = age, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = age - ageSE, ymax = age + ageSE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of age") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_age