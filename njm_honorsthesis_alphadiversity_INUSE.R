#### Geographic Comparison of Gut Microbial Flexibility in TRES
# Microbiome Analysis
# By Natalie Morris
# Modified from code by: Worsley et al. 2021, Animal Microbiome

# to do list
# change lm -> lmer + add random effects and covariates
# get versions of alpha data straight that shit confusing
# check decontam --- make sure not taking out major taxa!!!
  # contamination is often from the samples themselves
library(phyloseq)
library(dplyr)
library(ggplot2)
library(kableExtra)

# read in phyloseq object
geo_phylo <- readRDS("~/Library/CloudStorage/OneDrive-UniversityofArizona/Cornell/ht_njm/ht_njm_R/honors_thesis_njm/phyloseq_ht_njm_USE.rds")
# view sample data
peek <- sample_data(geo_phylo, errorIfNULL=TRUE)

##### RAREFY READS TO MIN SAMPLING DEPTH ######

# Look at lowest sample depth

sample_sums(geo_phylo) # lowest M17A0023 = 58 reads 

#rarefy to 1000 and set seed before rarefying (28367), so that results are reproducible  
physeqRare<-rarefy_even_depth(geo_phylo, 5000, rngseed = 28367)
# 35 samples removed

# OTUs removed after subsampling- leaves 6803 taxa and 393 samples
physeqRare
sample_sums(physeqRare)
sample_data(physeqRare)$read_depth <- sample_sums(physeqRare)

##### CALCULATE ALPHA DIVERSITY METRICS ####

# calculate chao1 (species richness), shannon using estimate_richness()
set.seed(787)
richnessEst<-estimate_richness(physeqRare, split=TRUE, measures= c("Chao1", "Shannon"))
head(richnessEst)
str(richnessEst)


#add alpha diversity metrics to metadata
physeqMeta <- as.data.frame(sample_data(physeqRare))
head(physeqMeta)
physeqMeta$Chao1 <- richnessEst$Chao1
physeqMeta$Shannon <- richnessEst$Shannon
head(physeqMeta)
str(physeqMeta) #420 samples

alphaData <- read.csv("alphaData_geo.csv")

# Independent variables should be vectors = factors (age, treatment, unit, box) #
alphaData$Location <- as.factor(alphaData$Location)
alphaData$Site <- as.factor(alphaData$Site)
alphaData$Nest <- as.factor(alphaData$Nest)
alphaData$Individual_Treatment <- as.factor(alphaData$Individual_Treatment)
alphaData$Individual_Band <- as.factor(alphaData$Individual_Band)
alphaData$Status <- as.factor(alphaData$Status)
alphaData$Capture_Number <- as.factor(alphaData$Capture_Number)
alphaData$Capture_Date <- as.factor(alphaData$Capture_Date)
alphaData$Manipulated_after_1st_capture <- as.factor(alphaData$Manipulated_after_1st_capture)

#Dependent variables should be numeric (any alpha diversity metric, cort_base, cort_stress, mass, headbill. etc) #
str(alphaData) # Check factors vs numeric

library(dplyr)
alphaData <- alphaData %>% 
  filter(Capture_Number !="2")
write.csv(alphaData, "alphaData_geo_incprov.csv")
### Filter to isolate individuals with samples at both time points, 
# to be able to calculate differences in alpha diversity and mass #

F_alphaData <- alphaData %>%
  group_by(Individual_Band) %>%
  filter(all(c(1, 3) %in% Capture_Number)) %>%
  ungroup()
table(F_alphaData$Location) # 58 Alaska samples, 60 Tennessee samples, 82 New York samples

## Filter to isolate individuals with samples at both times, and control treatment only for provisioning
# to only have samples from unmanipulated birds
FC_alphaData <- F_alphaData %>%
  filter(!(Individual_Treatment %in% c("Dulled", "Low_Tape", "Predator") & Capture_Number == "3"))

# Subset to only control birds for differences
FCD_alphaData <-subset(FC_alphaData, Individual_Treatment == "Control")
table(FCD_alphaData$Location)
# Add column with difference in alpha diversity and mass between captures

FCD_alphaData <- FCD_alphaData %>%
  group_by(Individual_Band) %>%
  mutate(Chao1_Difference = diff(Chao1[Capture_Number %in% c(1, 3)]))
FCD_alphaData <- FCD_alphaData %>%
  group_by(Individual_Band) %>%  
  mutate(Shannon_Difference = diff(Shannon[Capture_Number %in% c(1, 3)]))
FCD_alphaData <- FCD_alphaData %>%
  group_by(Individual_Band) %>%
  mutate(Mass_Difference = diff(Mass[Capture_Number %in% c(1, 3)]))

# Add column Capture with Incubation or Provisioning, based on Capture Number
alphaData$Capture <- ifelse(alphaData$Capture_Number == 1, "Incubation", "Provisioning")

#Reorder so locations in alphabetical order (only need to do once and it will apply to all plots)

alphaData$Location <- factor(alphaData$Location, levels=c("Alaska", "New_York", "Tennessee", "Wyoming"))

#Reorder so captures are in order (only need to do once and it will apply to all plots)

alphaData$Capture_Number <- factor(alphaData$Capture_Number, levels=c("1", "3"))

# Filter to isolate individuals with samples at both time, control treatment during provisioning, at time point of INCUBATION
FCI_alphaData <- FCD_alphaData %>%
  filter(!(Capture_Number == "3"))

# Split data into incubation and provisioning

incubation <-subset(alphaData, Capture_Number == "1")
table(incubation$Capture_Number) #224 samples

provisioning <-subset(alphaData, Capture_Number == "3")
provisioning <-subset(provisioning, Individual_Treatment == "Control")
table(provisioning$Capture_Number) #53
samples

# Split data into Location

# Alaska

table(alphaData_Alaska$Location) #72 samples
Alaska <- subset(FCI_alphaData, Location == "Alaska") # Alaska, filtered, control, incubation


# Tennessee

Tennessee <- subset(FCI_alphaData, Location == "Tennessee") # TN, filtered, control, incubation
table(Tennessee$Location) #15 samples

# Wyoming

alphaData_Wyoming <-subset(alphaData,Location == "Wyoming")
table(alphaData_Wyoming$Location) #73 samples---ONLY INCUBATION!!!!
Wyoming <- subset(alphaData, Location == "Wyoming")

# New York

alphaData_New_York <- subset(alphaData, Location == "New_York")
table(alphaData_New_York$Location) # 159 samples
New_York <- subset(FCI_alphaData, Location == "New_York") # NY, filtered, control, incubation

################ Check Data, Transform as Needed ##############################
# check data spread
hist(alphaData$Shannon) # negatively skewed ~ tail to the left
hist(alphaData$Chao1) # slightly negatively skewed? 
hist(FCI_alphaData$Shannon_Difference)
hist(FCI_alphaData$Chao1_Difference)


# transform with Tukey
library(rcompanion)
alphaData$ShannonT <- transformTukey(alphaData$Shannon)
alphaData$Chao1T <- transformTukey(alphaData$Chao1)


location_colors <- c("#B31B1B", "#B2B4B2", "blue", "green")
ggplot(alphaData, aes(x = Location, y = ShannonT, fill = factor(Capture))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(title = "ShannonT by Location during Incubation and Provisioning",
       x = "Location", y = "ShannonT", fill = "Capture") +
  theme_bw() +
  scale_fill_manual(values = location_colors)

ggplot(alphaData, aes(x = Location, y = Chao1T, fill = factor(Capture))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(title = "ChaoT by Location during Incubation and Provisioning",
       x = "Location", y = "ChaoT", fill = "Capture") +
  theme_bw() +
  scale_fill_manual(values = location_colors)

##### LM Log Alpha Microbial Diversity x Population #####
library(lme4)
library(lmerTest)
library(emmeans)
library(gtsummary)

# Shannon Tukey at Inc x Location

##### LM Alpha Microbial Diversity x Population #####

# ShannonT at Inc x Location 
mod_STincub_location <- lmer(ShannonT ~ Location + Age + (1|Site), data = incubation) # 
summary(mod_STincub_location, ddf = "Kenward-Roger")
tab_model(mod_STincub_location, p.val = "kr", show.df = TRUE, file = "mod_results_ST_incub_location.html")

# Posthoc pairwise tests

emmeans_stI_L <- emmeans(mod_STincub_location, list(pairwise~Location), adjust = "tukey") # 
emmeans_stI_L

ggplot(incubation, aes(x = Location, y = ShannonT, color = Location)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Shannon (Tukey) at Incubation",
       x = "Shannon", y = "Location", fill = "Location") +
  scale_color_manual(values = location_colors)



library(ggplot2)
library(emmeans)
library(ggpubr)

# Extract pairwise comparisons
pairwise_df <- as.data.frame(emmeans_stI_L$`pairwise differences of Location`)

# Prepare comparisons for plotting
pairwise_df <- pairwise_df %>%
  mutate(
    p_label = ifelse(p.value < 0.001, "***",
                     ifelse(p.value < 0.01, "**",
                            ifelse(p.value < 0.05, "*", "ns"))),  # Star significance notation
    group1 = sub(" - .*", "", contrast),  # Extract first location
    group2 = sub(".*- ", "", contrast),   # Extract second location
    y.position = seq(max(incubation$ShannonT) + 0.1, by = 0.1, length.out = n()))

# Boxplot with significance brackets
ggplot(incubation, aes(x = Location, y = ShannonT, color = Location)) +
  geom_boxplot(alpha = 0.5) +  
  stat_pvalue_manual(pairwise_df, label = "p_label", 
                     xmin = "group1", xmax = "group2", y.position = max(incubation$ShannonT) + 0.1) +  
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Shannon (Tukey) at Incubation",
       x = "Location", y = "Shannon Diversity Index", fill = "Location") +
  scale_color_manual(values = location_colors)


colnames(pairwise_df)
colnames(pairwise_df)[1] <- "contrast"


library(ggplot2)
library(ggpubr)
library(dplyr)

# Step 1: Ensure the contrast column is named correctly
colnames(pairwise_df)[1] <- "contrast"

# Step 2: Format comparisons correctly
pairwise_df <- pairwise_df %>%
  mutate(
    contrast = as.character(contrast),  # Ensure it's character type
    p_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"  # Not significant
    ),
    group1 = sub(" - .*", "", contrast),  # Extract first location
    group2 = sub(".*- ", "", contrast)    # Extract second location
  ) %>%
  filter(p_label != "ns") %>%  # 🚀 Remove non-significant comparisons
  mutate(
    y.position = seq(max(incubation$ShannonT) + 0.2,  # Start higher
                     by = 0.15,  # Increase spacing between brackets
                     length.out = n())  
  )


# Check column names in 'incubation'
colnames(incubation)

# If 'Location' is not a factor, convert it
incubation$Location <- as.factor(incubation$Location)

# Step 3: Plot with proper brackets
ggplot(incubation, aes(x = Location, y = ShannonT, color = Location)) +  # Use 'fill' to color the boxes
  geom_boxplot(position = position_dodge(width = 0.8)) +
  stat_pvalue_manual(pairwise_df, label = "p_label", 
                     xmin = "group1", xmax = "group2", 
                     y.position = "y.position") +  
  scale_color_manual(values = location_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Shannon (Tukey) at Incubation",
       x = "Location", y = "Shannon Diversity Index (Tukey)", fill = "Location") 




# ChaoT at Inc x Location 
mod_CTincub_location <- lmer(Chao1T ~ Location + Age + (1|Site), data = incubation) # 
summary(mod_CTincub_location, ddf = "Kenward-Roger")
tab_model(mod_STincub_location, p.val = "kr", show.df = TRUE, file = "mod_results_ST_incub_location.html")

# Posthoc pairwise tests

emmeans_ctI_L <- emmeans(mod_CTincub_location, list(pairwise~Location), adjust = "tukey") # LOG FAITH
emmeans_ctI_L

# Extract pairwise comparisons
pairwise_df <- as.data.frame(emmeans_ctI_L$`pairwise differences of Location`)

colnames(pairwise_df)[1] <- "contrast"

# Step 2: Format comparisons correctly
pairwise_df <- pairwise_df %>%
  mutate(
    contrast = as.character(contrast),  # Ensure it's character type
    p_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"  # Not significant
    ),
    group1 = sub(" - .*", "", contrast),  # Extract first location
    group2 = sub(".*- ", "", contrast)    # Extract second location
  ) %>%
  filter(p_label != "ns") %>%  # 🚀 Remove non-significant comparisons
  mutate(
    y.position = seq(max(incubation$Chao1T) + 0.2,  # Start higher
                     by = 0.15,  # Increase spacing between brackets
                     length.out = n())  
  )


# Step 3: Plot with proper brackets
ggplot(incubation, aes(x = Location, y = Chao1T, color = Location)) +  # Use 'fill' to color the boxes
  geom_boxplot(position = position_dodge(width = 0.8)) +
  stat_pvalue_manual(pairwise_df, label = "p_label", 
                     xmin = "group1", xmax = "group2", 
                     y.position = "y.position") +  
  scale_color_manual(values = location_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Chao (Tukey) at Incubation",
       x = "Location", y = "Chao Diversity Index (Tukey)", fill = "Location") 


# Shannon Difference

mod_s_diff_location <- lmer(Shannon_Difference ~ Location + ShannonT + (1|Site), data = FCD_alphaData)
summary(mod_s_diff_location, ddf = "Kenward-Roger")
tab_model(mod_LOGfaithpd_diff_location, p.val = "kr", show.df = TRUE, file = "mod_results_LOGfaithpd_diff_location.html")

# Posthoc pairwise tests

emmeans_sd_L <- emmeans(mod_s_diff_location, list(pairwise~Location), adjust = "tukey")
emmeans_sd_L


ggplot(FCD_alphaData, aes(x = Location, y = Shannon_Difference, color = Location)) +  # Use 'fill' to color the boxes
  geom_boxplot(position = position_dodge(width = 0.8)) +
  
  scale_color_manual(values = location_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Shannon Difference Across Locations",
       x = "Location", y = "Difference in Shannon", fill = "Location") 


### Chao1 Difference


mod_c_diff_location <- lmer(Chao1_Difference ~ Location + Chao1T + (1|Site), data = FCD_alphaData)
summary(mod_c_diff_location, ddf = "Kenward-Roger")
tab_model(mod_LOGfaithpd_diff_location, p.val = "kr", show.df = TRUE, file = "mod_results_LOGfaithpd_diff_location.html")

# Posthoc pairwise tests

emmeans_cd_L <- emmeans(mod_c_diff_location, list(pairwise~Location), adjust = "tukey")
emmeans_cd_L

ggplot(FCD_alphaData, aes(x = Location, y = Chao1_Difference, color = Location)) +  # Use 'fill' to color the boxes
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_color_manual(values = location_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Chao1 Difference Across Locations",
       x = "Location", y = "Difference in Chao1", fill = "Location") 




##### LM Phenotypic Flexibility: Mass Loss x Alpha Diversity (initial and change) ##### 
### Filter to isolate individuals with samples at both time points

F_alphaData <- alphaData %>%
  group_by(Individual_Band) %>%
  filter(all(c(1, 3) %in% Capture_Number)) %>%
  ungroup()
table(F_alphaData$Location) # 58 Alaska samples, 60 Tennessee samples, 82 New York samples

# check data spread
hist(FCD_alphaData$Shannon_Difference)
# negatively skewed ~ tail to the left
hist(FCD_alphaData$Chao1_Difference)

FCIM_alphaData <- FCI_alphaData %>%
  filter(!is.na(Mass_Difference))
table(FCIM_alphaData$Location) # removes 1 NY indv

### Mass by capture
mod_mass <- lmer(Mass ~ Capture_Number*Location + Age + (1|Site), data = FC_alphaData)
summary(mod_mass, ddf = "Kenward-Roger")

emm <- emmeans(mod_mass, list(pairwise~Capture_Number*Location), adjust = "tukey")

contrast_df <- as.data.frame(emm$contrasts)
contrast_df %>%
  kable(caption = "emmeans: Mass vs Cap*Loc") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))

# Mass at incubation by Location

mod_mass_inc_location <- lmer(Mass ~ Location + Age + (1|Site), data = incubation)
summary(mod_mass_inc_location, ddf = "Kenward-Roger")
tab_model(mod_mass_inc_location, p.val = "kr", show.df = TRUE, file = "mod_results_mass_inc_location.html")

emmeans(mod_mass_inc_location, list(pairwise~Location), adjust = "tukey")
# NY - AK insignificant, TN-AK sig, TN-NY sig

# Mass at provisioning by Location

mod_mass_prov_location <- lm(Mass ~ Location, data = provisioning)
summary(mod_mass_prov_location, ddf = "Kenward-Roger")
tab_model(mod_mass_prov_location, p.val = "kr", show.df = TRUE, file = "mod_results_mass_prov_location.html")

emmeans(mod_mass_prov_location, list(pairwise~Location), adjust = "tukey")
# all insignificant

# Mass difference x Location

mod_mass_diff_location <- lmer(Mass_Difference ~ Location + Age + (1|Site), data = FCIM_alphaData)
summary(mod_mass_diff_location, ddf = "Kenward-Roger")
tab_model(mod_mass_diff_location, p.val = "kr", show.df = TRUE, file = "mod_results_mass_diff_location.html")

# Posthoc pairwise tests

emmeans(mod_mass_diff_location, list(pairwise~Location), adjust = "tukey")

# Mass Difference x ShannoT at Inc

mod_MD_stInc <- lmer(Mass_Difference ~ scale(ShannonT)*Location + Age + (1|Site), data = FCIM_alphaData)
summary(mod_MD_stInc, ddf = "Kenward-Roger")
tab_model(mod_MD_FPDInc, p.val = "kr", show.df = TRUE, file = "mod_MD_FPDInc.html")

# Posthoc pairwise tests

emmeans(mod_MD_stInc, list(pairwise~Location), adjust = "tukey")


# Mass Difference x Chao1 at Inc

mod_MD_ctInc <- lmer(Mass_Difference ~ scale(Chao1T)*Location + Age + (1|Site), data = FCIM_alphaData)
summary(mod_MD_ctInc, ddf = "Kenward-Roger")
tab_model(mod_MD_FPDInc, p.val = "kr", show.df = TRUE, file = "mod_MD_FPDInc.html")

# Posthoc pairwise tests

emmeans(mod_MD_ctInc, list(pairwise~Location), adjust = "tukey")



# Mass Difference x Faith PD at Inc - TN 

mod_MD_SIncTN <- lm(Mass_Difference ~ scale(ShannonT), data = Tennessee)
summary(mod_MD_SIncTN, ddf = "Kenward-Roger")
tab_model(mod_MD_FPDIncTN, p.val = "kr", show.df = TRUE, file = "mod_MD_FPDIncTN.html")



# Mass Difference x Faith PD at Inc - NY #################### - Significant

mod_MD_SIncNY <- lm(Mass_Difference ~ scale(ShannonT), data = New_York)
summary(mod_MD_SIncNY, ddf = "Kenward-Roger")
tab_model(mod_MD_FPDIncNY, p.val = "kr", show.df = TRUE, file = "mod_MD_FPDIncNY.html")



# Mass Difference x Faith PD at Inc - AK # - Not significant

mod_MD_STIncAK <- lm(Mass_Difference ~ scale(ShannonT)
                     , data = Alaska)
summary(mod_MD_STIncAK, ddf = "Kenward-Roger") # p = 0.887
tab_model(mod_MD_FPDIncAK, p.val = "kr", show.df = TRUE, file = "mod_MD_FPDIncAK.html")



#### Mass Difference x Faith PD Difference - Not significant ####
View(FCIM_alphaData)
################## or mass 2 ~ shannon 2 + shannon1 + mass 1
# All locations
mod_MD_SD_L<- lmer(Mass_Difference ~ Shannon_Difference + Location + Age + (1|Site), data = FCIM_alphaData)
summary(mod_MD_SD_L, ddf = "Kenward-Roger")
tab_model(mod_MD_FPDD_L, p.val = "kr", show.df = TRUE, file = "mod_MD_FPDD.html")

# Posthoc pairwise tests
emmeans(mod_MD_SD_L, list(pairwise~Location), adjust = "tukey")


# All locations
mod_MD_CD_L <- lmer(Mass_Difference ~ Chao1_Difference +Location + Age + (1|Site), data = FCIM_alphaData)
summary(mod_MD_CD_L, ddf = "Kenward-Roger")
tab_model(mod_MD_FPDD_L, p.val = "kr", show.df = TRUE, file = "mod_MD_FPDD.html")

# Posthoc pairwise tests
emmeans(mod_MD_CD_L, list(pairwise~Location), adjust = "tukey")

# Mass Difference x Faith PD Difference - TN

# TN
mod_MD_FPDD_TN <- lm(Mass_Difference ~ Shannon_Difference, data = Tennessee)
summary(mod_MD_FPDD_TN, ddf = "Kenward-Roger") # p - 0,48
tab_model(mod_MD_FPDD_TN, p.val = "kr", show.df = TRUE, file = "mod_MD_FPDD_TN.html")

# Mass Difference x Faith PD Difference - NY

# NY
mod_MD_FPDD_NY <- lm(Mass_Difference ~ Shannon_Difference, data = New_York)
summary(mod_MD_FPDD_NY, ddf = "Kenward-Roger") # p - 0.543
tab_model(mod_MD_FPDD_NY, p.val = "kr", show.df = TRUE, file = "mod_MD_FPDD_NY.html")


# Mass Difference x Faith PD Difference - AK

# AK
mod_MD_FPDD_AK <- lm(Mass_Difference ~ Shannon_Difference, data = Alaska)
summary(mod_MD_FPDD_AK, ddf = "Kenward-Roger") # p - 0.6318
tab_model(mod_MD_FPDD_AK, p.val = "kr", show.df = TRUE, file = "mod_MD_FPDD_AK.html")



p <- ggplot(FCI_alphaData, aes(x=Individual_Treatment, y=Mass_Difference, fill = Location)) +
  geom_boxplot()
p

mod_IncMass_FPD <- lm(FaithPD ~ Mass, data = Alaska)
summary(mod_IncMass_FPD, ddf = "Kenward-Roger")
# NY - 0.00489 **

ggplot(FCI_alphaData, aes(x = Mass, y = ShannonT, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Mass at Incubation vs Shannon at Incubation",
       x = "Mass (g) at Incubation", y = "Faith's PD at Incubation", fill = "Location") +
  scale_color_manual(values = location_colors)

mod_IncMass_Fledge <- lm(Fledged_Num ~ Mass, data = FCI_alphaData)
summary(mod_IncMass_Fledge, ddf = "Kenward-Roger")

ggplot(FCIF_alphaData, aes(x = Location, y = Fledged_Num, fill = Location)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(title = "Number Fledged by Location",
       x = "Location", y = "Number Fledged", fill = "Location") +
  theme_bw() +
  scale_fill_manual(values = location_colors)

mod_Fledge_Loc <- lm(Fledged_Num ~ Location, data = FCIF_alphaData)
summary(mod_Fledge_Loc, ddf = "Kenward-Roger") # not significant


ggplot(FCIF_alphaData, aes(x = Mass, y = Fledged_Num, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Mass at Incubation vs Fledged Num",
       x = "Mass (g) at Incubation", y = "Fledged Num", fill = "Location") +
  scale_color_manual(values = location_colors)

mod_MassxFledge <- lm(Fledged_Num ~ Mass*Location, data = FCIF_alphaData)
summary(mod_MassxFledge, ddf = "Kenward-Roger") # AK significant

mod_FPDxFledge <- lm(Fledged_Num ~ logfaith*Location, data = FCIF_alphaData)
summary(mod_FPDxFledge, ddf = "Kenward-Roger") # AK significant
# if error: dev.off()

avg_mass <- FC_alphaData %>%
  group_by(Capture_Number, Location) %>%
  summarise(AvgMass = mean(Mass, na.rm = TRUE))
View(avg_mass)

table(FCIM_alphaData$Location)

######## Plots
library(ggplot2)
location_colors <- c("#f77f00", "#B31B1B", "#4682B4")

# Mass Loss
avg_mass <- FC_alphaData %>%
  group_by(Capture, Location) %>%
  summarise(AvgMass = mean(Mass, na.rm = TRUE))
View(avg_mass)

# Plot avg mass at Incubation to Provisioning
ggplot(avg_mass, aes(x= Capture_Number, y = AvgMass, group = Location, color = Location)) +
  geom_line(size = .75) +
  geom_point(size = 3) +
  labs( y = "Body mass (g)") +
  scale_color_manual(values = c("New_York" = "#B31B1B",  # Cornell red for Location A
                                "Tennessee" = "#f77f00",  # Cornell red for Location B
                                "Alaska" = "#4682B4"))  +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 12, color = "black")) +
  coord_cartesian(ylim = c(18.75, 22.8)) +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  theme(axis.line = element_line(color = "black", size = 0.25))
ggsave("geo_bodymass_locations.jpeg", plot = last_plot(), device = "jpeg")


ggplot(FCI_alphaData, aes(x = Location, y = ShannonT, fill = Location)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(title = "Log FaithPD by Location",
       x = "Location", y = "Log Faith's PD at Incubation", fill = "Location") +
  theme_bw() +
  scale_fill_manual(values = location_colors)

ggplot(FCI_alphaData, aes(x = Location, y = LogFaithPD_Difference, fill = Location)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(title = "Log FaithPD Difference by Location",
       x = "Location", y = "Log Faith's PD Difference", fill = "Location") +
  theme_bw() +
  scale_fill_manual(values = location_colors)

ggplot(FC_alphaData, aes(x = Capture_Number, y = FaithPD, fill = Location)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(title = "FaithPD by Location and Capture",
       x = "Location", y = "Faith's PD", fill = "Location") +
  theme_bw() +
  scale_fill_manual(values = location_colors)

ggplot(FCIM_alphaData, aes(x = Location, y = Mass_Difference, fill = Location)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(title = "Mass Difference by Location",
       x = "Location", y = "Body Mass Difference (g)", fill = "Location") +
  theme_bw() +
  scale_fill_manual(values = location_colors)


ggplot(FCIM_alphaData, aes(x = ShannonT, y = Mass_Difference, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Mass Loss vs ShannonT At Incubation",
       x = "ShannonT at Incubation", y = "Change in Mass (g)", fill = "Location") +
  scale_color_manual(values = location_colors)

# TN - Inc x Mass Loss
ggplot(Tennessee, aes(x = FaithPD, y = Mass_Difference, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Mass Loss vs Faith's PD At Incubation - TN",
       x = "Faith's PD at Incubation", y = "Change in Mass (g)", fill = "Location") +
  xlim(2.5, 13.5) + ylim(-5.5,0) +
  coord_cartesian(ylim = c(-5.5, 0)) +
  coord_cartesian(xlim = c(2.5, 13.5)) +
  scale_x_continuous(n.breaks=6) +
  scale_y_continuous(n.breaks=5) +
  scale_color_manual(values = location_colors)

# NY - Inc x Mass Loss
ggplot(New_York, aes(x = FaithPD, y = Mass_Difference, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Mass Loss vs Faith's PD At Incubation - NY",
       x = "Faith's PD at Incubation", y = "Change in Mass (g)", fill = "Location") +
  xlim(2.5, 13.5) + ylim(-5.5,0) +
  coord_cartesian(ylim = c(-5.5, 0)) +
  coord_cartesian(xlim = c(2.5, 13.5)) +
  scale_x_continuous(n.breaks=6) +
  scale_y_continuous(n.breaks=6) +
  scale_color_manual(values = "#B31B1B")

# AK
ggplot(Alaska, aes(x = FaithPD, y = Mass_Difference, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Mass Loss vs Faith's PD At Incubation - AK",
       x = "Faith's PD at Incubation", y = "Change in Mass (g)", fill = "Location") +
  xlim(2.5, 13.5) + ylim(-5.5,0) +
  coord_cartesian(ylim = c(-5.5,0)) +
  coord_cartesian(xlim = c(2.5, 13.5)) +
  scale_x_continuous(n.breaks=6) +
  scale_y_continuous(n.breaks=6) +
  scale_color_manual(values = "#4682B4")



# TN - Inc x Mass Loss
ggplot(Tennessee, aes(x = ShannonT, y = Mass_Difference, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  labs(title = "Mass Loss vs Faith's PD At Incubation - TN",
       x = "ShannonT at Incubation", y = "Change in Mass (g)", fill = "Location") +
  xlim(3, 7) + ylim(-5, 0) +
  coord_cartesian(ylim = c(-5,0)) +
  coord_cartesian(xlim = c(3, 7)) +
  scale_color_manual(values = location_colors) +
  scale_x_continuous(breaks = seq(3, 7, by = 1)) +
  scale_y_continuous(breaks = seq(-6, 0, by = 1)) 


# NY - Inc x Mass Loss
ggplot(New_York, aes(x = FaithPD, y = Mass_Difference, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  labs(title = "Mass Loss vs Faith's PD At Incubation - NY",
       x = "Faith's PD at Incubation", y = "Change in Mass (g)", fill = "Location") +
  xlim(2.5, 13.5) + ylim(-5.5, 0) +
  coord_cartesian(ylim = c(-5.5,0)) +
  coord_cartesian(xlim = c(2.5, 13.5)) +
  scale_x_continuous(breaks = seq(2, 13, by = 2)) +
  scale_y_continuous(breaks = seq(-5, 0, by = 1)) +
  scale_color_manual(values = "#B31B1B")

# AK
ggplot(Alaska, aes(x = FaithPD, y = Mass_Difference, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  labs(title = "Mass Loss vs Faith's PD At Incubation - AK",
       x = "Faith's PD at Incubation", y = "Change in Mass (g)", fill = "Location") +
  xlim(2, 5.5) + ylim(-5.5, 0) +  
  coord_cartesian(ylim = c(-5.5,0)) +
  coord_cartesian(xlim = c(2, 5.5)) +
  scale_x_continuous(breaks = seq(2, 5, by = 1)) +
  scale_y_continuous(breaks = seq(-5, 0, by = 1)) +
  scale_color_manual(values = "#4682B4")

# Faith PD Difference x Mass Loss
ggplot(FCIM_alphaData, aes(x = Shannon_Difference, y = Mass_Difference, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  labs(title = "Mass Loss vs Difference in Shannon",
       x = "Difference in Shannon (Prov-Inc)", y = "Change in Mass (g)", fill = "Location") +
  scale_color_manual(values = location_colors)

ggplot(Tennessee, aes(x = FaithPD_Difference, y = Mass_Difference, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  labs(title = "Mass Loss vs Difference in Faith's PD - TN",
       x = "Difference in Faith's PD (Prov-Inc)", y = "Change in Mass (g)", fill = "Location") +
  coord_cartesian(ylim = c(-5.5, 0)) +
  scale_color_manual(values = location_colors)

ggplot(New_York, aes(x = FaithPD_Difference, y = Mass_Difference, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  labs(title = "Mass Loss vs Difference in Faith's PD - NY",
       x = "Difference in Faith's PD (Prov-Inc)", y = "Change in Mass (g)", fill = "Location") +
  coord_cartesian(ylim = c(-5.5, 0)) +
  scale_color_manual(values = "#B31B1B")

ggplot(Alaska, aes(x = FaithPD_Difference, y = Mass_Difference, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  labs(title = "Mass Loss vs Difference in Faith's PD - AK",
       x = "Difference in Faith's PD (Prov-Inc)", y = "Change in Mass (g)", fill = "Location") +
  coord_cartesian(ylim = c(-5.5, 0)) +
  scale_color_manual(values = "#4682B4")
####

citation("pairwiseAdonis")
citation("stats")


######## Plot: 

location_colors <- c("#f77f00", "#B31B1B", "#4682B4")

F_alphaData %>%
  ggplot(aes(x = Capture_Number, y = Chao1, group = Individual_Band)) +
  geom_point(aes(color = Individual_Band), alpha = 0.7, size = 2, 
             position = position_dodge(width = 0.4)) +
  geom_line(aes(color = Individual_Band), alpha = 0.75, 
            position = position_dodge(width = 0.4)) +
  geom_boxplot(aes(group = Capture_Number), 
               alpha = 0.1, fill = "grey", color = "black", 
               outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.4)) +
  
  facet_wrap(~Location) +  # Show separate panels by Location
  scale_color_viridis_d() +
  labs(x = "Capture Number", 
       y = "Chao1 Gut Microbial Diversity", 
       title = "Chao1 Diversity from Incubation to Provisioning across Locations") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  
        axis.title = element_text(size = 14),                  
        axis.text = element_text(size = 12),                   
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5))


############## MAASS PLOTTTT


data_TN <- FC_alphaData %>% filter(Location == "Tennessee")
data_NY <- FC_alphaData %>% filter(Location == "New_York")
data_NY <- data_NY %>% filter(!is.na(Mass))
data_AK <- FC_alphaData %>% filter(Location == "Alaska")

# Calculate mean mass and standard error for each time point at each location
summary_TN <- data_TN %>%
  group_by(Capture_Number) %>%
  summarise(mean_mass = mean(Mass),
            std_error = sd(Mass) / sqrt(n()),
            sd_mass = sd(Mass))


summary_NY <- data_NY %>%
  group_by(Capture_Number) %>%
  summarise(mean_mass = mean(Mass),
            std_error = sd(Mass) / sqrt(n()),
            sd_mass = sd(Mass))


summary_AK <- data_AK %>%
  group_by(Capture_Number) %>%
  summarise(mean_mass = mean(Mass),
            std_error = sd(Mass) / sqrt(n()),
            sd_mass = sd(Mass))


# Display summaries
print("TN:")
print(summary_TN)
print("NY:")
print(summary_NY)
print("AK:")
print(summary_AK)

summary_combined <- rbind(summary_TN, summary_NY, summary_AK)
summary_combined$Location <- factor(rep(c("Tennessee", "New_York", "Alaska"), each = 2))

# Function to calculate the slope of mass loss from capture 1 to capture 3

slope <- lm(Mass ~ Capture_Number, data = data_TN)
slope

slope <- lm(Mass ~ Capture_Number, data = data_NY)
slope

slope <- lm(Mass ~ Capture_Number, data = data_AK)
slope

# Plot
ggplot(summary_combined, aes(x = Capture_Number, y = mean_mass, group = Location, color = Location)) +
  geom_errorbar(aes(ymin = mean_mass - std_error, ymax = mean_mass + std_error), width = 0.1) +
  geom_line() +
  geom_point(size = 3) +
  labs(title = "",
       x = "Capture",
       y = "Body Mass (g)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_x_discrete(breaks = c(1, 3), labels = c("Incubation", "Provisioning")) +
  scale_color_manual(values = c("New_York" = "#B31B1B",  # Cornell red for Location A
                                "Tennessee" = "#f77f00",  # Cornell red for Location B
                                "Alaska" = "#4682B4"))
