#### Geographic Comparison of Gut Microbial Flexibility in TRES
# Microbiome Analysis
# By Natalie Morris
# Modified from code by: Worsley et al. 2021, Animal Microbiome and Jenn Houtz

# to do list
# change lm -> lmer + add random effects and covariates
# get versions of alpha data straight that shit confusing
# check decontam --- make sure not taking out major taxa!!!
  # contamination is often from the samples themselves
library(phyloseq)
library(dplyr)
library(ggplot2)
library(kableExtra)
library(ggpubr)
library(rcompanion)
library(lme4)
library(lmerTest)
library(emmeans)
library(gtsummary)
library(microbiome)
library(MiscMetabar)


# read in phyloseq object
geo_phylo <- readRDS("~/Library/CloudStorage/OneDrive-UniversityofArizona/Cornell/ht_njm/ht_njm_R/honors_thesis_njm/phyloseq_ht_njm_USE.rds")
# view sample data
peek <- sample_data(geo_phylo, errorIfNULL=TRUE)

# remove samples not in study
geo_phylo<- subset_samples(geo_phylo, !Capture_Number == "2")
geo_phylo<- subset_samples(geo_phylo, !(Individual_Treatment %in% c("Dulled", "Low_Tape", "Predator") & Capture_Number == "3"))

# add run as a covariate
samdf <- data.frame(sample_data(geo_phylo))
# add "run" variable: 1 for NY, 2 for all others
samdf$run <- ifelse(samdf$Location == "New_York", 1, 2)
# push updated sample data back into phyloseq object
sample_data(geo_phylo) <- samdf


##### RAREFY READS TO MIN SAMPLING DEPTH ######

# Look at lowest sample depth

sort(sample_sums(geo_phylo)) # lowest M17A0023 = 58 reads 
sample_data(geo_phylo)$read_depth <- sample_sums(geo_phylo)

#rarefy to 5000 and set seed before rarefying (28367), so that results are reproducible  
physeqRare<-rarefy_even_depth(geo_phylo, 5000, rngseed = 28367)
# 35 samples removed

# OTUs removed after subsampling- leaves 6803 taxa and 393 samples
physeqRare
sample_sums(physeqRare)

saveRDS(physeqRare, file = "geo_phyloseq_rare_July2025.RDS")
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
str(physeqMeta) # 277 samples

alphaData <- meta(physeqMeta) ## use this to convert to dataframe to save! :) (package microbiome)
alphaData <- read.csv("alphaData_geo.csv")

######### non - rare version .....
richnessEst<-estimate_richness(geo_phylo, split=TRUE, measures= c("Chao1", "Shannon"))
#add alpha diversity metrics to metadata
physeqMeta <- as.data.frame(sample_data(geo_phylo))
head(physeqMeta)
physeqMeta$Chao1 <- richnessEst$Chao1
physeqMeta$Shannon <- richnessEst$Shannon
head(physeqMeta)
str(physeqMeta) # 277 samples

alphaData <- meta(physeqMeta) ## use this to convert to dataframe to save! :) 
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

write.csv(alphaData, file = "alphaData_geo_incprov_notrarefied.CSV", row.names = FALSE)
table(alphaData$Location, by = alphaData$Capture_Number)
### Filter to isolate individuals with samples at both time points, 
# to be able to calculate differences in alpha diversity and mass #

F_alphaData <- alphaData %>%
  group_by(Individual_Band) %>%
  filter(all(c(1, 3) %in% Capture_Number)) %>%
  ungroup()
table(F_alphaData$Location) # X Alaska samples, X Tennessee samples, X New York samples

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

FCD_alphaData <- FCD_alphaData %>%
  group_by(Individual_Band) %>%
  mutate(Chao1_DifferenceT = diff(Chao1T[Capture_Number %in% c(1, 3)]))
FCD_alphaData <- FCD_alphaData %>%
  group_by(Individual_Band) %>%  
  mutate(Shannon_DifferenceT = diff(ShannonT[Capture_Number %in% c(1, 3)]))

# Add column Capture with Incubation or Provisioning, based on Capture Number
alphaData$Capture <- ifelse(alphaData$Capture_Number == 1, "Incubation", "Provisioning")

#Reorder so locations in alphabetical order (only need to do once and it will apply to all plots)

alphaData$Location <- factor(alphaData$Location, levels=c("Alaska", "New_York", "Tennessee", "Wyoming"))
alphaData$Location <- factor(alphaData$Location, levels=c("New_York", "Tennessee", "Wyoming", "Alaska"))

#Reorder so captures are in order (only need to do once and it will apply to all plots)

alphaData$Capture_Number <- factor(alphaData$Capture_Number, levels=c("1", "3"))

# Filter to isolate individuals with samples at both time, control treatment during provisioning, at time point of INCUBATION
FCI_alphaData <- FCD_alphaData %>%
  filter(!(Capture_Number == "3"))

FCP_alphaData <- FCD_alphaData %>%
  filter(!(Capture_Number == "1"))

# Split data into incubation and provisioning

incubation <-subset(alphaData, Capture_Number == "1")
table(incubation$Location) #224 samples

provisioning <-subset(alphaData, Capture_Number == "3")
provisioning <-subset(provisioning, Individual_Treatment == "Control")
table(provisioning$Location) #53
samples

# Split data into Location

# Alaska

table(Alaska$Location) #72 samples
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

# remove Wyoming 
alphaData2 <- filter(alphaData, !Location == "Wyoming")
################ Check Data, Transform as Needed ##############################
# check data spread
hist(alphaData$Shannon) # negatively skewed ~ tail to the left
hist(alphaData$Chao1) # slightly negatively skewed? 
hist(FCI_alphaData$Shannon_Difference)
hist(FCI_alphaData$Chao1_Difference)


# transform with Tukey
alphaData$ShannonT <- transformTukey(alphaData$Shannon)
alphaData$Chao1T <- transformTukey(alphaData$Chao1)


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


# Shannon Tukey at Inc x Location

##### LM Alpha Microbial Diversity x Population #####

# ShannonT at Inc x Location 
mod_STincub_location <- lmer(ShannonT ~ Location + Age + run + (1|Site), data = incubation) # 
summary(mod_STincub_location, ddf = "Kenward-Roger")


mod_STincub_location %>% tbl_regression(exponentiate = FALSE,
  pvalue_fun = label_style_pvalue(digits = 3),) |>
     bold_p(t = 0.05) |>
    bold_labels() |>
    italicize_levels() %>% 
    as_gt() %>% gt::tab_caption("Shannon (Tukey) across Locations during Incubation ")

# Posthoc pairwise tests

emmeans_stI_L <- emmeans(mod_STincub_location, pairwise~Location, adjust = "tukey") # 
emmeans_stI_L

contrast_df <- as.data.frame(emmeans_stI_L$contrasts)
contrast_df %>%
  kable(caption = "emmeans: ShannonT across locations during incubation") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))

ggplot(incubation, aes(x = Location, y = ShannonT, fill = Location)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Shannon (Tukey) at Incubation",
       x = "Shannon", y = "Location", fill = "Location") +
  scale_fill_manual(values = location_colors)




location_colors <- c("#4682B4", "#B31B1B","#f77f00", "darkgreen"  )


# Extract pairwise comparisons
pairwise_df <- as.data.frame(emmeans_stI_L$contrasts)
colnames(pairwise_df)[1] <- "contrast"

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



incubation$Location <- as.factor(incubation$Location)

ggplot(incubation, aes(x = Location, y = ShannonT, color = Location, fill = Location)) +  # Use 'fill' to color the boxes
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.75) +
  stat_pvalue_manual(pairwise_df, label = "p_label", 
                     xmin = "group1", xmax = "group2", 
                     y.position = "y.position") +  
  scale_color_manual(values = location_colors) +
  scale_fill_manual(values = location_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Shannon (Tukey) at Incubation",
       x = "Location", y = "Shannon Diversity Index (Tukey)", fill = "Location") 




# ChaoT at Inc x Location 
mod_CTincub_location <- lmer(Chao1T ~ Location + Age + (1|Site), data = incubation) # 
summary(mod_CTincub_location, ddf = "Kenward-Roger")
qqnorm(resid(mod_CTincub_location))


mod_CTincub_location %>% tbl_regression(exponentiate = FALSE,
                                        pvalue_fun = label_style_pvalue(digits = 3)) |>
  bold_p(t = 0.05) |>
  bold_labels() |>
  italicize_levels() %>% 
  as_gt() %>% gt::tab_caption("Chao1T (Tukey) across Locations during Incubation ")

library(gtsummary)
library(gt)

# Posthoc pairwise tests

emmeans_ctI_L <- emmeans(mod_CTincub_location, pairwise ~ Location, adjust = "tukey") # 
emmeans_ctI_L

contrast_df <- as.data.frame(emmeans_ctI_L$contrasts)
contrast_df %>%
  kable(caption = "emmeans: Chao1T across locations during incubation") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))

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
ggplot(incubation, aes(x = Location, y = Chao1T, fill = Location, color = Location)) +  # Use 'fill' to color the boxes
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.75) +
  stat_pvalue_manual(pairwise_df, label = "p_label", 
                     xmin = "group1", xmax = "group2", 
                     y.position = "y.position") +  
  scale_fill_manual(values = location_colors) +
  scale_color_manual(values = location_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Chao (Tukey) at Incubation",
       x = "Location", y = "Chao Diversity Index (Tukey)", fill = "Location") 


# Shannon Difference

mod_s_diff_location <- lmer(Shannon_DifferenceT ~ Location + ShannonT + Age + (1|Site), data = FCI_alphaData)
summary(mod_s_diff_location, ddf = "Kenward-Roger")

mod_s_diff_location <- lmer(ShannonT ~ Location*Capture + Age + (1|Site) +(1|Individual_Band), data = F_alphaData)
summary(mod_s_diff_location, ddf = "Kenward-Roger")

mod_s_diff_location %>% tbl_regression(exponentiate = FALSE, 
                                     pvalue_fun = label_style_pvalue(digits = 3),) |>
  bold_p(t = 0.05) |>
  bold_labels() |>
  italicize_levels() %>% 
  as_gt() %>% gt::tab_caption("Shannon (Tukey) across Locations and Captures ")

# Posthoc pairwise tests

emmeans_sd_L <- emmeans(mod_s_diff_location, list(pairwise~Location), adjust = "tukey")
emmeans_sd_L


# Extract pairwise comparisons
pairwise_df <- as.data.frame(emmeans_sd_L$`pairwise differences of Location, Capture`)
colnames(pairwise_df)[1] <- "contrast"

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
    y.position = seq(max(F_alphaData$ShannonT) + 0.2,  # Start higher
                     by = 0.15,  # Increase spacing between brackets
                     length.out = n())  
  )



ggplot(FCI_alphaData, aes(x = Location, y = Shannon_DifferenceT, color = Location, fill = Location)) +  # Use 'fill' to color the boxes
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.75) +
  scale_color_manual(values = location_colors) +
  scale_fill_manual(values = location_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Shannon (Tukey) between Incubation and Provisioning ",
       x = "Capture", y = "Shannon Diversity Index (Tukey)", fill = "Location") 

contrast_df <- as.data.frame(emmeans_sd_L$`pairwise differences of Location, Capture`)
contrast_df %>%
  kable(caption = "emmeans: ShannonT between captures and locations") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), font_size = 10)

### Chao1 Difference


mod_c_diff_location <- lmer(Chao1_DifferenceT ~ Location + Chao1T + Age + (1|Site), data = FCI_alphaData)
summary(mod_c_diff_location, ddf = "Kenward-Roger")

mod_c_diff_location <- lmer(Chao1T ~ Location*Capture + Age + (1|Site) + (1|Individual_Band), data = F_alphaData)
summary(mod_c_diff_location, ddf = "Kenward-Roger")

mod_c_diff_location %>% tbl_regression(exponentiate = FALSE, 
                                       pvalue_fun = label_style_pvalue(digits = 3),) |>
  bold_p(t = 0.05) |>
  bold_labels() |>
  italicize_levels() %>% 
  as_gt() %>% gt::tab_caption("Chao1 (Tukey) across Locations and Captures ")

# Posthoc pairwise tests

emmeans_cd_L <- emmeans(mod_c_diff_location, list(pairwise~Location), adjust = "tukey")
emmeans_cd_L

emmeans_cd_L <- emmeans(mod_c_diff_location, list(pairwise~Location*Capture), adjust = "tukey")
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

ggplot(F_alphaData, aes(x = Capture, y = Chao1T, color = Location, fill = Location)) +  # Use 'fill' to color the boxes
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.75) +
  facet_wrap(~Location) +
  scale_color_manual(values = location_colors) +
  scale_fill_manual(values = location_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Chao1 (Tukey) between Incubation and Provisioning ",
       x = "Capture", y = "Chao1 Diversity Index (Tukey)", fill = "Location") 


contrast_df <- as.data.frame(emmeans_cd_L$`pairwise differences of Location, Capture`)
contrast_df %>%
  kable(caption = "emmeans: Chao1T between captures and locations") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), font_size = 10)
########## Provisioning ############

# ShannonT at Prov x Location 
mod_STp_location <- lmer(ShannonT ~ Location + Age + (1|Site), data = provisioning) # 
summary(mod_STp_location, ddf = "Kenward-Roger")
tab_model(mod_STincub_location, p.val = "kr", show.df = TRUE, file = "mod_results_ST_incub_location.html")

# Posthoc pairwise tests

emmeans_stp_L <- emmeans(mod_STp_location, list(pairwise~Location), adjust = "tukey") # 
emmeans_stp_L


location_colors <- c("#4682B4", "#B31B1B","#f77f00", "darkgreen"  )


# Extract pairwise comparisons
pairwise_df <- as.data.frame(emmeans_stp_L$`pairwise differences of Location`)
colnames(pairwise_df)[1] <- "contrast"

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
 # filter(p_label != "ns") %>%  # 🚀 Remove non-significant comparisons
  mutate(
    y.position = seq(max(provisioning$ShannonT) + 0.2,  # Start higher
                     by = 0.15,  # Increase spacing between brackets
                     length.out = n())  
  )



provisioning$Location <- as.factor(provisioning$Location)

ggplot(provisioning, aes(x = Location, y = ShannonT, color = Location, fill = Location)) +  # Use 'fill' to color the boxes
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.5) +
  stat_pvalue_manual(pairwise_df, label = "p_label", 
                     xmin = "group1", xmax = "group2", 
                     y.position = "y.position") +  
  scale_color_manual(values = location_colors) +
  scale_fill_manual(values = location_colors) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Shannon (Tukey) at Provisioning",
       x = "Location", y = "Shannon Diversity Index (Tukey)", fill = "Location")


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

mod_mass_diff_location %>% tbl_regression() %>% 
  as_gt() %>% gt::tab_caption("Mass Difference across Locations")

# Posthoc pairwise tests

emm <- emmeans(mod_mass_diff_location, list(pairwise~Location), adjust = "tukey")

contrast_df <- as.data.frame(emm$`pairwise differences of Location`)
contrast_df %>%
  kable(caption = "emmeans: Mass Difference across Locations") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))

# Mass Difference x ShannoT at Inc

mod <- lmer(Mass_Difference ~ scale(ShannonT) + Location + Age +run + (1|Site), data = FCIM_alphaData)
summary(mod, ddf = "Kenward-Roger")
tab_model(mod_MD_FPDInc, p.val = "kr", show.df = TRUE, file = "mod_MD_FPDInc.html")

mod %>%
  tbl_regression(exponentiate = FALSE,
                 pvalue_fun = label_style_pvalue(digits = 3)) |>
  bold_p(t = 0.05) |>
  bold_labels() |>
  italicize_levels() %>%
  as_gt() %>%
  gt::tab_caption("Chao1T at Incubation vs Mass Loss") %>%
  gt::tab_source_note(
    source_note = paste("n =", nrow(FCIM_alphaData))
  )
# Posthoc pairwise tests

emmeans(mod_MD_stInc, list(pairwise~Location), adjust = "tukey")


# Mass Difference x Chao1 at Inc

mod <- lmer(Mass_Difference ~ scale(Chao1T) + Location + run+ Age + (1|Site), data = FCIM_alphaData)
summary(mod, ddf = "Kenward-Roger")
tab_model(mod_MD_FPDInc, p.val = "kr", show.df = TRUE, file = "mod_MD_FPDInc.html")

# Posthoc pairwise tests

emmeans(mod_MD_ctInc, list(pairwise~Location), adjust = "tukey")



# Mass Difference x GM at Inc - TN 

mod_MD_SIncTN <- lm(Mass_Difference ~ scale(ShannonT), data = Tennessee)
summary(mod_MD_SIncTN, ddf = "Kenward-Roger")

mod_MD_SIncTN %>% tbl_regression(exponentiate = FALSE,
                                        pvalue_fun = label_style_pvalue(digits = 3),) |>  bold_p(t = 0.05) |>
  bold_labels() |>
  italicize_levels() %>% 
  as_gt() %>% gt::tab_caption("TN: Shannon (Tukey) Incubation vs Mass Loss")

mod_MD_SIncTN <- lm(Mass_Difference ~ scale(Chao1T), data = Tennessee)
summary(mod_MD_SIncTN, ddf = "Kenward-Roger")
qqnorm(resid(mod_MD_SIncTN))

mod_MD_SIncTN %>% tbl_regression(exponentiate = FALSE,
                                 pvalue_fun = label_style_pvalue(digits = 3),) |>  bold_p(t = 0.05) |>
  bold_labels() |>
  italicize_levels() %>% 
  as_gt() %>% gt::tab_caption("TN: Chao (Tukey) Incubation vs Mass Loss")

# Mass Difference x GM at Inc - NY ##

mod_MD_SIncNY <- lm(Mass_Difference ~ scale(ShannonT), data = New_York)
summary(mod_MD_SIncNY, ddf = "Kenward-Roger")

mod_MD_SIncNY %>% tbl_regression(exponentiate = FALSE,
                                 pvalue_fun = label_style_pvalue(digits = 3),) |>  bold_p(t = 0.05) |>
  bold_labels() |>
  italicize_levels() %>% 
  as_gt() %>% gt::tab_caption("NY: Shannon (Tukey) Incubation vs Mass Loss")


# Mass Difference x Faith PD at Inc - AK # - Not significant

mod_MD_STIncAK <- lm(Mass_Difference ~ scale(ShannonT), data = Alaska)
summary(mod_MD_STIncAK, ddf = "Kenward-Roger") # p = 0.887

mod_MD_STIncAK %>%
  tbl_regression(exponentiate = FALSE,
                 pvalue_fun = label_style_pvalue(digits = 3)) |>
  bold_p(t = 0.05) |>
  bold_labels() |>
  italicize_levels() %>%
  as_gt() %>%
  gt::tab_caption("AK: Shannon (Tukey) Incubation vs Mass Loss") %>%
  gt::tab_source_note(
    source_note = paste("n =", nrow(Alaska))
  )

mod_MD_STIncAK <- lm(Mass_Difference ~ scale(Chao1T), data = Alaska)
summary(mod_MD_STIncAK, ddf = "Kenward-Roger") # p = 0.887

mod_MD_STIncAK %>%
  tbl_regression(exponentiate = FALSE,
                 pvalue_fun = label_style_pvalue(digits = 3)) |>
  bold_p(t = 0.05) |>
  bold_labels() |>
  italicize_levels() %>%
  as_gt() %>%
  gt::tab_caption("AK: Chao1 (Tukey) Incubation vs Mass Loss") %>%
  gt::tab_source_note(
    source_note = paste("n =", nrow(Alaska))
  )

####### mass 

mod <- lmer(Mass ~ scale(ShannonT)*Capture*Location + Age + run+  (1|Individual_Band) 
            + (1|Site), data = alphaData)
summary(mod, ddf = "Kenward-Roger") # 

visreg(mod, "ShannonT", by = "Capture",
       type = "conditional",  # Controls for random effects
       ylab = "Mass", 
       xlab = "Shannon (Tukey)",
       gg = TRUE) 

ggplot(alphaData, aes(x = Chao1T, y = Mass, color = Location, fill = Location)) + geom_point(alpha = 0.4) +
  scale_color_manual(values = location_colors) +
 facet_wrap(~Capture) + 
  geom_smooth(method = "lm") +
 scale_fill_manual(values = location_colors) +
  theme_bw() +
  labs(
    x = "Chao1 Diversity Index (Tukey)",
    y = "Mass (g)",
    title = "Chao1 diversity as a predictor of mass across the breeding season",
  )
  
m <- emmeans(mod, list(pairwise~scale(ShannonT)*Location*Capture), adjust = "tukey")

contrast_df <- as.data.frame(m$`pairwise differences of ShannonT, Location, Capture`)
contrast_df %>%
  kable(caption = "emmeans: ShannonT between captures and locations") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), font_size = 10)


mod <- lmer(Mass ~ Chao1T*Capture_Number + Location + Age + (1|Individual_Band) 
            + (1|Site), data = alphaData2)
summary(mod, ddf = "Kenward-Roger")

mod_MD_S_L<- lmer(Mass_Difference ~ scale(ShannonT)*Location + Age + Mass + (1|Site), data = FCIM_alphaData)
summary(mod_MD_S_L, ddf = "Kenward-Roger")

emmeans(mod_MD_S_L, list(pairwise~Location), adjust = "tukey")

mod_MD_C_L<- lmer(Mass_Difference ~ scale(Chao1T)*Location + Age + Mass + (1|Site), data = FCIM_alphaData)
summary(mod_MD_C_L, ddf = "Kenward-Roger")

emmeans(mod_MD_S_L, list(pairwise~Location), adjust = "tukey")

ggplot(FCIM_alphaData, aes(x = Shannon_DifferenceT, y = Mass_Difference, color = Location, fill = Location)) + geom_point(alpha = 0.5) +
  scale_color_manual(values = location_colors) +
  #  facet_wrap(~Capture) + 
  geom_smooth(method = "lm") +
  scale_fill_manual(values = location_colors) +
  theme_bw() +
  labs(
    x = "Difference in Shannon Diversity Index (Tukey)",
    y = "Mass Loss (g)",
    title = "Difference in Shannon diversity as a predictor of mass loss",
  )
##### prov
mod_MD_S_L<- lmer(Mass_Difference ~ scale(ShannonT)*Location + Age + (1|Site), data = FCP_alphaData)
summary(mod_MD_S_L, ddf = "Kenward-Roger")

emmeans(mod_MD_S_L, list(pairwise~Location), adjust = "tukey")

mod_MD_C_L<- lmer(Mass_Difference ~ scale(Chao1T)*Location + Age + (1|Site), data = FCP_alphaData)
summary(mod_MD_C_L, ddf = "Kenward-Roger")

emmeans(mod_MD_C_L, list(pairwise~Location), adjust = "tukey")

###### mass difference x GM 

mod <- lmer(Mass_Difference ~ scale(ShannonT)*Capture*Location + Age + (1|Individual_Band) 
            + (1|Site), data = FCD_alphaData)
summary(mod, ddf = "Kenward-Roger") # 

visreg(mod, "ShannonT", by = "Capture",
       type = "conditional",  # Controls for random effects
       ylab = "Mass", 
       xlab = "Shannon (Tukey)",
       gg = TRUE) 

ggplot(FCD_alphaData, aes(x = ShannonT, y = Mass_Difference, color = Location, fill = Location)) + geom_point(alpha = 0.4) +
  scale_color_manual(values = location_colors) +
  facet_wrap(~Capture) + 
  geom_smooth(method = "lm") +
  scale_fill_manual(values = location_colors) +
  theme_bw() +
  labs(
    x = "Shannon Diversity Index (Tukey)",
    y = "Mass Loss (g)",
    title = "Shannon diversity as a predictor of mass loss across the breeding season",
  )



#### Mass Difference x GM Difference - Not significant ####
View(FCIM_alphaData)
################## or mass 2 ~ shannon 2 + shannon1 + mass 1
# All locations
mod_MD_SD_L<- lmer(Mass_Difference ~ scale(Shannon_DifferenceT) + Location + Age + (1|Site), data = FCIM_alphaData)
summary(mod_MD_SD_L, ddf = "Kenward-Roger")

# Posthoc pairwise tests
emmeans(mod_MD_SD_L, list(pairwise~Location), adjust = "tukey")


# All locations
mod_MD_CD_L <- lmer(Mass_Difference ~ scale(Chao1_DifferenceT) + Location + Age + (1|Site), data = FCIM_alphaData)
summary(mod_MD_CD_L, ddf = "Kenward-Roger")

# Posthoc pairwise tests
emmeans(mod_MD_CD_L, list(pairwise~Location), adjust = "tukey")

# Mass Difference x Faith PD Difference - TN

# TN
mod <- lm(Mass_Difference ~ scale(Shannon_DifferenceT), data = Tennessee)
summary(mod, ddf = "Kenward-Roger") # p - 0,48

mod %>%
  tbl_regression(exponentiate = FALSE,
                 pvalue_fun = label_style_pvalue(digits = 3)) |>
  bold_p(t = 0.05) |>
  bold_labels() |>
  italicize_levels() %>%
  as_gt() %>%
  gt::tab_caption("TN: Difference in Chao1T vs Mass Loss") %>%
  gt::tab_source_note(
    source_note = paste("n =", nrow(Tennessee))
  )

mod <- lm(Mass_Difference ~ scale(Chao1_DifferenceT), data = Tennessee)
summary(mod_MD_FPDD_TN, ddf = "Kenward-Roger") 
# Mass Difference x Faith PD Difference - NY

# NY
mod <- lm(Mass_Difference ~ scale(Shannon_DifferenceT), data = New_York)
summary(mod_MD_FPDD_NY, ddf = "Kenward-Roger") # p - 0.543

mod %>%
  tbl_regression(exponentiate = FALSE,
                 pvalue_fun = label_style_pvalue(digits = 3)) |>
  bold_p(t = 0.05) |>
  bold_labels() |>
  italicize_levels() %>%
  as_gt() %>%
  gt::tab_caption("NY: Difference in Chao1T vs Mass Loss") %>%
  gt::tab_source_note(
    source_note = paste("n =", nrow(New_York))
  )

mod <- lm(Mass_Difference ~ scale(Chao1_DifferenceT), data = New_York)
summary(mod_MD_FPDD_NY, ddf = "Kenward-Roger")


# Mass Difference x Faith PD Difference - AK

# AK
mod <- lm(Mass_Difference ~ scale(Shannon_DifferenceT), data = Alaska)
summary(mod_MD_FPDD_AK, ddf = "Kenward-Roger")

mod %>%
  tbl_regression(exponentiate = FALSE,
                 pvalue_fun = label_style_pvalue(digits = 3)) |>
  bold_p(t = 0.05) |>
  bold_labels() |>
  italicize_levels() %>%
  as_gt() %>%
  gt::tab_caption("AK: Difference in Chao1T vs Mass Loss") %>%
  gt::tab_source_note(
    source_note = paste("n =", nrow(Alaska))
  )
mod <- lm(Mass_Difference ~ scale(Chao1_DifferenceT), data = Alaska)
summary(mod_MD_FPDD_AK, ddf = "Kenward-Roger") # p - 0.6318



ggplot(FCI_alphaData, aes(x = Mass, y = ShannonT, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Mass at Incubation vs Shannon at Incubation",
       x = "Mass (g) at Incubation", y = "Shannon (Tukey) at Incubation", fill = "Location") +
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
location_colors <- c("#f77f00", "#B31B1B", "#4682B4")

# Mass Loss
avg_mass <- alphaData %>%
  group_by(Capture, Location) %>%
  summarise(AvgMass = mean(Mass, na.rm = TRUE))
View(avg_mass)

# Plot avg mass at Incubation to Provisioning
ggplot(alphaData, aes(x = Capture, y = Mass, group = Location, color = Location, fill = Location)) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.1) +
  geom_smooth(method = lm, se = FALSE) +
  labs(y = "Body mass (g)") +
  scale_color_manual(values = location_colors) +
  scale_fill_manual(values = location_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = (c(0.8, 0.8)),
        axis.line = element_line(color = "black", size = 0.25))


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

# 
ggplot(FCIM_alphaData, aes(x = Chao1T, y = Mass_Difference, color = Location)) +
geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Mass Loss vs Chao1T At Incubation",
       x = "Chao1T at Incubation", y = "Change in Mass (g)", fill = "Location") +
  scale_color_manual(values = location_colors)



# MB Difference x Mass Loss
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

### Chao Difference 

ggplot(FCIM_alphaData, aes(x = Chao1_Difference, y = Mass_Difference, color = Location)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  labs(title = "Mass Loss vs Difference in Chao1",
       x = "Difference in Chao1 (Prov-Inc)", y = "Change in Mass (g)", fill = "Location") +
  scale_color_manual(values = location_colors)

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
data_WY <- alphaData %>% filter(Location == "Wyoming")

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


summary_WY <- data_WY %>%
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
  scale_color_manual(values = c("New_York" = "#B31B1B",  # Cornell red
                                "Tennessee" = "#f77f00",  
                                "Alaska" = "#4682B4"))
