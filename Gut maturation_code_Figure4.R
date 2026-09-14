########################### Gut microbial maturation_Result 4 ################################################

#### Figure 4. Identification of bacterial species with long-term and short-term significance on microbial maturation ####

# Load required packages
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(picante)
library(ggpubr)
library(betapart)
library(openxlsx)
library(tidyr)
library(SpiecEasi)
library(igraph)
library(ComplexUpset)
library(reshape2)
library(broom)
library(ggrepel)
library(introdataviz)
library(writexl)

# Load phyloseq object
load("phy_feces_clean.RData")
phy_feces_clean#


#####  Figure 4A. Differentially abundant microbes between groups across ages #####
phy_feces_clean#

# Species level aggregation
phy_species <- tax_glom(phy_feces_clean, taxrank = "species")
phy_species <- prune_taxa(taxa_sums(phy_species) > 0, phy_species)
#save(phy_species, file = "phy_species.RData")


###### a. Differential abundance analysis across all time points ######
phy_species#

## 1 week
# prepare data filter
phy_species_1w <- phy_species%>%subset_samples(Time_new=="1w"&!is.na(Gut_maturation_Week_two))
phy_species_1w <- prune_taxa(taxa_sums(phy_species_1w)>0,phy_species_1w)
phy_species_1w

# Wilcoxon test
res_wil_species_1w <- DA.wil(phy_species_1w, predictor = "Gut_maturation_Week_two")
res_wil_species_1w$Significant <- ifelse(res_wil_species_1w$pval.adj<0.05,"Yes","No") %>%
  factor(levels = c("Yes","No"))
table(res_wil_species_1w$ordering,res_wil_species_1w$Significant)

# Extract significant species names
res_wil_species_1w_High <- subset(res_wil_species_1w, ordering=="High>Low"&pval.adj<0.05)
res_wil_species_1w_Low <- subset(res_wil_species_1w, ordering=="Low>High"&pval.adj<0.05)
res_wil_species_1w_High_name <- res_wil_species_1w_High$species
res_wil_species_1w_Low_name <- res_wil_species_1w_Low$species


## 1 month
# prepare data filter
phy_species_1m <- phy_species%>%subset_samples(Time_new=="1m"&!is.na(Gut_maturation_Week_two))
phy_species_1m <- prune_taxa(taxa_sums(phy_species_1m)>0,phy_species_1m)
phy_species_1m

# Wilcoxon test
res_wil_species_1m <- DA.wil(phy_species_1m, predictor = "Gut_maturation_Week_two")
res_wil_species_1m$Significant <- ifelse(res_wil_species_1m$pval.adj<0.05,"Yes","No") %>%
  factor(levels = c("Yes","No"))
table(res_wil_species_1m$ordering,res_wil_species_1m$Significant)

# Extract significant species names
res_wil_species_1m_High <- subset(res_wil_species_1m, ordering=="High>Low"&pval.adj<0.05)
res_wil_species_1m_Low <- subset(res_wil_species_1m, ordering=="Low>High"&pval.adj<0.05)
res_wil_species_1m_High_name <- res_wil_species_1m_High$species
res_wil_species_1m_Low_name <- res_wil_species_1m_Low$species


## 1 year
# prepare data filter
phy_species_1y <- phy_species%>%subset_samples(Time_new=="1y"&!is.na(Gut_maturation_Week_two))
phy_species_1y <- prune_taxa(taxa_sums(phy_species_1y)>0,phy_species_1y)
phy_species_1y

# Wilcoxon test
res_wil_species_1y <- DA.wil(phy_species_1y, predictor = "Gut_maturation_Week_two")
res_wil_species_1y$Significant <- ifelse(res_wil_species_1y$pval.adj<0.05,"Yes","No") %>%
  factor(levels = c("Yes","No"))
table(res_wil_species_1y$ordering,res_wil_species_1y$Significant)

# Extract significant species names
res_wil_species_1y_High <- subset(res_wil_species_1y, ordering=="High>Low"&pval.adj<0.05)
res_wil_species_1y_Low <- subset(res_wil_species_1y, ordering=="Low>High"&pval.adj<0.05)
res_wil_species_1y_High_name <- res_wil_species_1y_High$species
res_wil_species_1y_Low_name <- res_wil_species_1y_Low$species


## 4 year
# prepare data filter
phy_species_4y <- phy_species%>%subset_samples(Time_new=="4y"&!is.na(Gut_maturation_Week_two))
phy_species_4y <- prune_taxa(taxa_sums(phy_species_4y)>0,phy_species_4y)
phy_species_4y

# Wilcoxon test
res_wil_species_4y <- DA.wil(phy_species_4y, predictor = "Gut_maturation_Week_two")
res_wil_species_4y$Significant <- ifelse(res_wil_species_4y$pval.adj<0.05,"Yes","No") %>%
  factor(levels = c("Yes","No"))
table(res_wil_species_4y$ordering,res_wil_species_4y$Significant)

# Extract significant species names
res_wil_species_4y_High <- subset(res_wil_species_4y, ordering=="High>Low"&pval.adj<0.05)
res_wil_species_4y_Low <- subset(res_wil_species_4y, ordering=="Low>High"&pval.adj<0.05)
res_wil_species_4y_High_name <- res_wil_species_4y_High$species
res_wil_species_4y_Low_name <- res_wil_species_4y_Low$species


## 5 year
# prepare data filter
phy_species_5y <- phy_species%>%subset_samples(Time_new=="5y"&!is.na(Gut_maturation_Week_two))
phy_species_5y <- prune_taxa(taxa_sums(phy_species_5y)>0,phy_species_5y)
phy_species_5y

# Wilcoxon test
res_wil_species_5y <- DA.wil(phy_species_5y, predictor = "Gut_maturation_Week_two")
res_wil_species_5y$Significant <- ifelse(res_wil_species_5y$pval.adj<0.05,"Yes","No") %>%
  factor(levels = c("Yes","No"))
table(res_wil_species_5y$ordering,res_wil_species_5y$Significant)

# Extract significant species names
res_wil_species_5y_High <- subset(res_wil_species_5y, ordering=="High>Low"&pval.adj<0.05)
res_wil_species_5y_Low <- subset(res_wil_species_5y, ordering=="Low>High"&pval.adj<0.05)
res_wil_species_5y_High_name <- res_wil_species_5y_High$species
res_wil_species_5y_Low_name <- res_wil_species_5y_Low$species


## 6 year
# prepare data filter
phy_species_6y <- phy_species%>%subset_samples(Time_new=="6y"&!is.na(Gut_maturation_Week_two))
phy_species_6y <-prune_taxa(taxa_sums(phy_species_6y)>0,phy_species_6y)
phy_species_6y

# Wilcoxon test
res_wil_species_6y <- DA.wil(phy_species_6y, predictor = "Gut_maturation_Week_two")
res_wil_species_6y$Significant <- ifelse(res_wil_species_6y$pval.adj<0.05,"Yes","No") %>%
  factor(levels = c("Yes","No"))
table(res_wil_species_6y$ordering,res_wil_species_6y$Significant)

# Extract significant species names
res_wil_species_6y_High <- subset(res_wil_species_6y, ordering=="High>Low"&pval.adj<0.05)
res_wil_species_6y_Low <- subset(res_wil_species_6y, ordering=="Low>High"&pval.adj<0.05)
res_wil_species_6y_High_name <- res_wil_species_6y_High$species
res_wil_species_6y_Low_name <- res_wil_species_6y_Low$species


## Combine all significant results into one table
res_wil_species_1w; res_wil_species_1m; res_wil_species_1y;
res_wil_species_4y; res_wil_species_5y; res_wil_species_6y;

res_wil_species_1w_sig <- res_wil_species_1w %>% subset(pval.adj<0.05) %>% 
  mutate(Time="1 week") %>% arrange(ordering)
res_wil_species_1m_sig <- res_wil_species_1m %>% subset(pval.adj<0.05) %>% 
  mutate(Time="1 month")%>% arrange(ordering)
res_wil_species_1y_sig <- res_wil_species_1y %>% subset(pval.adj<0.05) %>% 
  mutate(Time="1 year")%>% arrange(ordering)
res_wil_species_4y_sig <- res_wil_species_4y %>% subset(pval.adj<0.05) %>% 
  mutate(Time="4 year")%>% arrange(ordering)
res_wil_species_5y_sig <- res_wil_species_5y %>% subset(pval.adj<0.05) %>% 
  mutate(Time="5 year")%>% arrange(ordering)
res_wil_species_6y_sig <- res_wil_species_6y %>% subset(pval.adj<0.05) %>% 
  mutate(Time="6 year")%>% arrange(ordering)

res_wil_species_alltime_sig <- rbind(res_wil_species_1w_sig, res_wil_species_1m_sig,
                                     res_wil_species_1y_sig, res_wil_species_4y_sig,
                                     res_wil_species_5y_sig, res_wil_species_6y_sig)
write_xlsx(res_wil_species_alltime_sig, "differential_abundance_species_alltime_sig.xlsx")




###### try 2 ##########

library(Maaslin2)

# phy_species <- tax_glom(phy_feces_clean, taxrank = "species")
phy_species <- prune_taxa(taxa_sums(phy_species) > 0, phy_species)

### Prepare species-level table

# Species names
species_name_raw <- as.character(tax_table(phy_species)[, "species"])

# Replace missing / empty species names by ASV/taxa IDs
species_name_raw[is.na(species_name_raw) | species_name_raw == ""] <-
  taxa_names(phy_species)[is.na(species_name_raw) | species_name_raw == ""]

# Keep original biological names for later annotation
species_name_original <- make.unique(species_name_raw)

# MaAsLin2-safe feature IDs
species_feature_id <- make.names(species_name_original, unique = TRUE)

# Save mapping table
feature_mapping <- data.frame(feature = species_feature_id,
                              species = species_name_original)

# Rename taxa for MaAsLin2
taxa_names(phy_species) <- species_feature_id

# Feature table: samples x species
otu_species <- otu_table(phy_species)

if (taxa_are_rows(otu_species)) {
  otu_species <- t(otu_species)
}

otu_species_df <- as.data.frame(otu_species)

# Metadata
meta_species <- data.frame(sample_data(phy_species))

# Match metadata to feature table
meta_species <- meta_species[rownames(otu_species_df), ]

# Format variables
meta_species <- meta_species %>%
  mutate(Gut_maturation_Week_two = factor(Gut_maturation_Week_two, levels = c("Low", "High")),
    Time_new = factor(Time_new, levels = c("1w", "1m", "1y", "4y", "5y", "6y")),
    child_ID = factor(child_ID),
    oldchild01 = factor(oldchild01),
    catanddog = factor(catanddog),
    apartment_house = factor(apartment_house),
    abbirth_mother = factor(abbirth_mother),
    sex = factor(sex),
    delivery01 = factor(delivery01),
    education = factor(education),
    birthseason = factor(birthseason))

# Check matching
stopifnot(all(rownames(otu_species_df) == rownames(meta_species)))


### MaAsLin2 by time point

run_maaslin2_by_time <- function(tp) {
  
  message("Running MaAsLin2 for time point: ", tp)
  
  # Subset metadata
  meta_tp <- meta_species %>%
    filter(Time_new == tp, !is.na(Gut_maturation_Week_two)) %>%
    droplevels()
  
  # Subset feature table
  otu_tp <- otu_species_df[rownames(meta_tp), , drop = FALSE]
  
  # Remove taxa absent in this time point
  otu_tp <- otu_tp[, colSums(otu_tp, na.rm = TRUE) > 0, drop = FALSE]
  
  # Candidate fixed effects
  fixed_effects_use <- c(
    "Gut_maturation_Week_two",
    "oldchild01",
    "catanddog",
    "apartment_house",
    "abbirth_mother",
    "sex"#,
    #"delivery01",
    #"education",
    #"birthseason"
  )
  
  # Keep only variables with enough variation
  fixed_effects_keep <- fixed_effects_use[
    fixed_effects_use %in% colnames(meta_tp)
  ]
  
  fixed_effects_keep <- fixed_effects_keep[
    sapply(fixed_effects_keep, function(v) {
      x <- meta_tp[[v]]
      if (is.numeric(x)) {
        return(sum(!is.na(x)) > 10 && sd(x, na.rm = TRUE) > 0)
      } else {
        return(length(unique(na.omit(x))) >= 2)
      }
    })
  ]
  
  # Make sure maturity group is retained
  if (!"Gut_maturation_Week_two" %in% fixed_effects_keep) {
    warning("Gut_maturation_Week_two has fewer than 2 levels at time point: ", tp)
    return(data.frame())
  }
  
  # Keep complete cases for selected metadata
  keep_samples <- complete.cases(meta_tp[, fixed_effects_keep, drop = FALSE])
  
  meta_tp2 <- meta_tp[keep_samples, , drop = FALSE] %>%
    droplevels()
  
  otu_tp2 <- otu_tp[rownames(meta_tp2), , drop = FALSE]
  
  # Remove zero-only taxa again after sample filtering
  otu_tp2 <- otu_tp2[, colSums(otu_tp2, na.rm = TRUE) > 0, drop = FALSE]
  
  # Optional abundance/prevalence filtering
  # Nonzero values are not required in both High and Low groups, as MaAsLin2 accepts zeros as part of the abundance data.
  min_total_nonzero <- 10
  min_prevalence <- 0.01
  
  keep_features <- sapply(colnames(otu_tp2), function(feature) {
    x <- otu_tp2[[feature]]
    present <- x > 0
    sum(present, na.rm = TRUE) >= min_total_nonzero &&
      mean(present, na.rm = TRUE) >= min_prevalence
  })
  
  otu_tp2 <- otu_tp2[, keep_features, drop = FALSE]
  
  message("Samples used: ", nrow(meta_tp2))
  message("Features used after filtering: ", ncol(otu_tp2))
  message("Fixed effects: ", paste(fixed_effects_keep, collapse = ", "))
  
  if (ncol(otu_tp2) == 0) {
    warning("No features left after filtering at time point: ", tp)
    return(data.frame())
  }
  
  outdir <- paste0("MaAsLin2_species_", tp)
  
  fit <- Maaslin2(
    input_data = otu_tp2,
    input_metadata = meta_tp2,
    output = outdir,
    fixed_effects = fixed_effects_keep,
    reference = "Gut_maturation_Week_two,Low",
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    correction = "BH",
    standardize = FALSE,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    cores = 4
  )
  
  res <- read.delim(file.path(outdir, "all_results.tsv"), check.names = FALSE)
  
  # MaAsLin2 output sometimes stores categorical variable as:
  # metadata == "Gut_maturation_Week_two" and value == "High"
  # or metadata == "Gut_maturation_Week_twoHigh"
  res_maturity <- res %>%
    filter(
      metadata == "Gut_maturation_Week_two" |
        metadata == "Gut_maturation_Week_twoHigh" |
        grepl("^Gut_maturation_Week_two", metadata)
    )
  
  if ("value" %in% colnames(res_maturity)) {
    res_maturity <- res_maturity %>%
      filter(is.na(value) | value == "" | value == "High")
  }
  
  res_maturity <- res_maturity %>%
    mutate(
      Time_new = tp,
      feature = as.character(feature),
      ordering = ifelse(coef > 0, "High>Low", "Low>High"),
      pval.adj = qval
    ) %>%
    left_join(feature_mapping, by = "feature")
  
  return(res_maturity)
}  




### Run all time points
time_points <- c("1w", "1m", "1y", "4y", "5y", "6y")

maaslin2_time_results <- lapply(time_points, run_maaslin2_by_time) %>%
  bind_rows()

table(maaslin2_time_results$Time_new)
View(maaslin2_time_results)

# Significant abundance-associated species
maaslin2_time_sig <- maaslin2_time_results %>%
  filter(qval < 0.05)

View(maaslin2_time_sig)

# Save results
# write.xlsx(maaslin2_time_results,
#  "MaAsLin2_covariate_adjusted_species_by_time_all_results.xlsx",rowNames = FALSE)

# write.xlsx( maaslin2_time_sig,
#  "MaAsLin2_covariate_adjusted_species_by_time_significant.xlsx",rowNames = FALSE)

#maaslin2_time_sig_8factor <- maaslin2_time_sig  # 151 
#maaslin2_time_sig_4factor <- maaslin2_time_sig # 166

write_xlsx(maaslin2_time_sig, "maaslin2_alltime_sig.xlsx")


## Extract adjusted significant species lists

get_maaslin2_names <- function(tp, direction) {
  maaslin2_time_sig %>%
    filter(Time_new == tp, ordering == direction) %>%
    pull(species) %>%
    unique() %>%
    na.omit()
}

res_ml2_species_1w_High_name <- get_maaslin2_names("1w", "High>Low")
res_ml2_species_1w_Low_name  <- get_maaslin2_names("1w", "Low>High")

res_ml2_species_1m_High_name <- get_maaslin2_names("1m", "High>Low")
res_ml2_species_1m_Low_name  <- get_maaslin2_names("1m", "Low>High")

res_ml2_species_1y_High_name <- get_maaslin2_names("1y", "High>Low")
res_ml2_species_1y_Low_name  <- get_maaslin2_names("1y", "Low>High")

res_ml2_species_4y_High_name <- get_maaslin2_names("4y", "High>Low")
res_ml2_species_4y_Low_name  <- get_maaslin2_names("4y", "Low>High")

res_ml2_species_5y_High_name <- get_maaslin2_names("5y", "High>Low")
res_ml2_species_5y_Low_name  <- get_maaslin2_names("5y", "Low>High")

res_ml2_species_6y_High_name <- get_maaslin2_names("6y", "High>Low")
res_ml2_species_6y_Low_name  <- get_maaslin2_names("6y", "Low>High")

# Quick check
sapply(list(High_1w = res_ml2_species_1w_High_name,
            High_1m = res_ml2_species_1m_High_name,
            High_1y = res_ml2_species_1y_High_name,
            High_4y = res_ml2_species_4y_High_name,
            High_5y = res_ml2_species_5y_High_name,
            High_6y = res_ml2_species_6y_High_name,
            Low_1w = res_ml2_species_1w_Low_name,
            Low_1m = res_ml2_species_1m_Low_name,
            Low_1y = res_ml2_species_1y_Low_name,
            Low_4y = res_ml2_species_4y_Low_name,
            Low_5y = res_ml2_species_5y_Low_name,
            Low_6y = res_ml2_species_6y_Low_name ),length)







###### b. Upset plot ######
res_wil_species_1w_High_name; res_wil_species_1w_Low_name #
res_wil_species_1m_High_name; res_wil_species_1m_Low_name #
res_wil_species_1y_High_name; res_wil_species_1y_Low_name #
res_wil_species_4y_High_name; res_wil_species_4y_Low_name #
res_wil_species_5y_High_name; res_wil_species_5y_Low_name #
res_wil_species_6y_High_name; res_wil_species_6y_Low_name #

### Prepare data for upset plot
## High group species
High_species <- list(High_1w = na.omit(res_ml2_species_1w_High_name),
                     High_1m = na.omit(res_ml2_species_1m_High_name),
                     High_1y = na.omit(res_ml2_species_1y_High_name),
                     High_4y = na.omit(res_ml2_species_4y_High_name),
                     High_5y = na.omit(res_ml2_species_5y_High_name),
                     High_6y = na.omit(res_ml2_species_6y_High_name))


# Helper function to convert list to data frame
from_list <- function(list_data) {
  members = unique(unlist(list_data))
  data.frame(
    lapply(list_data, function(set) members %in% set),
    row.names=members,
    check.names=FALSE
  )
}

# Prepare High group data
High_species2 <- from_list(High_species)
High_species2$species_name <- rownames(High_species2)

# Filter species that only appear in single time points
names_High=c("High_1w","High_1m","High_1y","High_4y","High_5y","High_6y")
High_species2$species_name[rowSums(High_species2[, names_High]) == 1] <- NA

# Create upset plot for High group
High_composition_venn <- ComplexUpset::upset(High_species2, rev(names_High), sort_intersections=FALSE,
                                             intersections=list(
                                               "High_1w","High_1m","High_1y","High_4y","High_5y","High_6y",
                                               
                                               c("High_1w","High_1m"),
                                               c("High_1w","High_1y"),#
                                               c("High_1w","High_1m","High_1y"),
                                               c("High_1m","High_1y"),
                                               
                                               #c("High_1w","High_1m","High_1y","High_4y"),
                                               #c("High_1w","High_1m","High_1y","High_4y","High_5y"),#
                                               #c("High_1w","High_1m","High_1y","High_4y","High_5y","High_6y"),#
                                               #c("High_1w","High_1m","High_1y","High_4y","High_6y"),#
                                               #c("High_1w","High_4y"),#
                                               #c("High_1w","High_1y","High_4y"),#
                                               c("High_1w","High_1m","High_4y"),
                                               
                                               c("High_1m","High_4y"),#
                                               c("High_1m","High_1y","High_4y"),
                                               #c("High_1m","High_1y","High_4y","High_5y"),#
                                               #c("High_1m","High_1y","High_4y","High_5y","High_6y"), #
                                               #c("High_1m","High_1y","High_5y","High_6y"),#
                                               #c("High_1m","High_1y","High_4y","High_6y"), 
                                               
                                               #c("High_1y","High_4y"),
                                               #c("High_1y","High_5y"),#
                                               #c("High_1y","High_6y"),#
                                               c("High_1y","High_4y","High_5y"),#
                                               c("High_1y","High_4y","High_5y","High_6y"),
                                               #c("High_1y","High_4y","High_6y"),
                                               #c("High_1y","High_5y","High_6y"),#
                                               
                                               c("High_4y","High_6y"),
                                               c("High_4y","High_5y","High_6y")#,
                                               #c("High_5y","High_6y")#
                                             ),
                                             
                                             height_ratio = 0.5,
                                             width_ratio = 0.2, 
                                             sort_sets=FALSE,
                                             stripes=c("#c6dbef", "#c6dbef","#ffffb3","#ffffb3","#ccebc5","#ccebc5"),
                                             wrap = T,
                                             
                                             base_annotations=list(
                                               'Frequency'=intersection_size(counts=TRUE,text = list(vjust=-1),
                                                                             mapping=aes(fill="bars_color"))+
                                                 scale_fill_manual(values = c('bars_color'="grey"),guide='none')+
                                                 theme(panel.grid = element_blank(),axis.line = element_line(),
                                                       axis.ticks.y=element_line())+
                                                 scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2),
                                                                    expand=expansion(mult=c(0,.1)))+
                                                 geom_text(mapping=aes(label=species_name),
                                                           position=position_stack(),
                                                           na.rm=TRUE,
                                                           vjust=-20,size=3)
                                             ),
                                             
                                             set_sizes = upset_set_size(position="left")+
                                               theme(panel.grid = element_blank(),
                                                     axis.ticks.x=element_line())+
                                               ylab("Total number of species")+
                                               geom_text(aes(label=..count..),hjust=1.1,stat="count"))

## Low group species
Low_species <- list(Low_1w = na.omit(res_ml2_species_1w_Low_name),
                    Low_1m = na.omit(res_ml2_species_1m_Low_name),
                    Low_1y = na.omit(res_ml2_species_1y_Low_name),
                    Low_4y = na.omit(res_ml2_species_4y_Low_name),
                    Low_5y = na.omit(res_ml2_species_5y_Low_name),
                    Low_6y = na.omit(res_ml2_species_6y_Low_name))

# Prepare Low group data
Low_species2 <- from_list(Low_species)
Low_species2$species_name <- rownames(Low_species2)

# Filter species that only appear in single time points
names_Low=c("Low_1w","Low_1m","Low_1y","Low_4y","Low_5y","Low_6y")
Low_species2$species_name[rowSums(Low_species2[, names_Low]) == 1] <- NA

# Create upset plot for Low group
Low_composition_venn <- ComplexUpset::upset(Low_species2, rev(names_Low), sort_intersections=FALSE,
                                            intersections=list(
                                              "Low_1w","Low_1m","Low_1y","Low_4y","Low_5y","Low_6y",
                                              
                                              c("Low_1w","Low_1m"),
                                              c("Low_1w","Low_1y"),#
                                              c("Low_1w","Low_1m","Low_1y"),
                                              c("Low_1m","Low_1y"),
                                              
                                              #c("Low_1w","Low_1m","Low_1y","Low_4y"),
                                              #c("Low_1w","Low_1m","Low_1y","Low_4y","Low_5y"),#
                                              #c("Low_1w","Low_1m","Low_1y","Low_4y","Low_5y","Low_6y"),#
                                              #c("Low_1w","Low_1m","Low_1y","Low_4y","Low_6y"),#
                                              #c("Low_1w","Low_4y"),#
                                              #c("Low_1w","Low_1y","Low_4y"),#
                                              c("Low_1w","Low_1m","Low_4y"),
                                              
                                              c("Low_1m","Low_4y"),#
                                              c("Low_1m","Low_1y","Low_4y"),
                                              #c("Low_1m","Low_1y","Low_4y","Low_5y"),#
                                              #c("Low_1m","Low_1y","Low_4y","Low_5y","Low_6y"), #
                                              #c("Low_1m","Low_1y","Low_5y","Low_6y"),#
                                              #c("Low_1m","Low_1y","Low_4y","Low_6y"), 
                                              
                                              #c("Low_1y","Low_4y"),
                                              #c("Low_1y","Low_5y"),#
                                              #c("Low_1y","Low_6y"),#
                                              c("Low_1y","Low_4y","Low_5y"),#
                                              c("Low_1y","Low_4y","Low_5y","Low_6y"),
                                              #c("Low_1y","Low_4y","Low_6y"),
                                              #c("Low_1y","Low_5y","Low_6y"),#
                                              
                                              c("Low_4y","Low_6y"),
                                              c("Low_4y","Low_5y","Low_6y")#,
                                              #c("Low_5y","Low_6y")#
                                            ),
                                            
                                            height_ratio = 0.5,
                                            width_ratio = 0.2, 
                                            sort_sets=FALSE,
                                            stripes=c("#c6dbef", "#c6dbef","#ffffb3","#ffffb3","#ccebc5","#ccebc5"),
                                            wrap = T,
                                            
                                            base_annotations=list(
                                              'Frequency'=intersection_size(counts=TRUE,text = list(vjust=-1),
                                                                            mapping=aes(fill="bars_color"))+
                                                scale_fill_manual(values = c('bars_color'="grey"),guide='none')+
                                                theme(panel.grid = element_blank(),axis.line = element_line(),
                                                      axis.ticks.y=element_line())+
                                                scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2),
                                                                   expand=expansion(mult=c(0,.1)))+
                                                geom_text(mapping=aes(label=species_name),
                                                          position=position_stack(),
                                                          na.rm=TRUE,
                                                          vjust=-20,size=3)),
                                            
                                            set_sizes = upset_set_size(position="left")+
                                              theme(panel.grid = element_blank(),
                                                    axis.ticks.x=element_line())+
                                              ylab("Total number of species")+
                                              geom_text(aes(label=..count..),hjust=1.1,stat="count"))

# Display both plots
High_composition_venn
Low_composition_venn
# save as 8.2 * 4.2




#####  Figure 4B. Trajectories of representative species from each group #####
library(emmeans)
library(lmerTest)
library(lme4)

phy_species#

### relative abundance_representative figure 

# Get bacteria from upset plot
res_Ma2_1w1msig_name <- c("Streptococcus mitis","Gemella sp.")
res_Ma2_1w1ysig_name <- c("Actinomyces sp.")

res_Ma2_1w1m1ysig_name <- c("Alistipes onderdonkii","Blautia wexlerae")
res_Ma2_1m1ysig_name <- c("Lachnospiraceae sp.","Collinsella aerofaciens","Faecalibacterium prausnitzii")
res_Ma2_1w1m4ysig_name <- c("Bacteroides sp.","Staphylococcus aureus", "Streptococcus salivarius")
res_Ma2_1m4ysig_name <- c("Bacteroides uniformis","Bacteroides vulgatus","Bilophila wadsworthia")
res_Ma2_1m1y4ysig_name <- c("Parabacteroides merdae")

res_Ma2_1y4y5ysig_name <- c("Odoribacter splanchnicus")
res_Ma2_1y4y5y6ysig_name <- c("Alistipes putredinis")
res_Ma2_4y6ysig_name <- c("Christensenellaceae R-7 group sp.","Clostridia vadinBB60 group sp.","UCG-010 sp.")
res_Ma2_4y5y6ysig_name <- c("UCG-002 sp.")


select_sig_taxa <- c(res_Ma2_1w1msig_name,res_Ma2_1w1ysig_name,
                     res_Ma2_1w1m1ysig_name,res_Ma2_1m1ysig_name,res_Ma2_1w1m4ysig_name,
                     res_Ma2_1m4ysig_name, res_Ma2_1m1y4ysig_name,
                     res_Ma2_1y4y5ysig_name,res_Ma2_1y4y5y6ysig_name,
                     res_Ma2_4y6ysig_name,res_Ma2_4y5y6ysig_name )
length(select_sig_taxa)

# Relative abundance table
time_levels <- c("1w","1m","1y","4y","5y","6y")
sample_data(phy_species)$Time_new <- factor(sample_data(phy_species)$Time_new, levels = time_levels)
phy_species_rela <- transform_sample_counts(phy_species,function(x) x/sum(x))


# Sample metadata
phy_species_rela_sample_df <- data.frame(sample_data(phy_species_rela))
phy_species_rela_sample_df$SampleID <- rownames(phy_species_rela_sample_df)

phy_species_rela_sample_df <- phy_species_rela_sample_df %>%
  mutate(
    Time_new = factor(Time_new, levels = time_levels),
    Gut_maturation_Week_two = factor(Gut_maturation_Week_two, levels = c("Low", "High")),
    child_ID = factor(child_ID),
    oldchild01 = factor(oldchild01),
    catanddog = factor(catanddog),
    apartment_house = factor(apartment_house),
    abbirth_mother = factor(abbirth_mother),
    sex = factor(sex)
  )

# OTU table: samples x taxa
otu_mat <- as(otu_table(phy_species_rela), "matrix")

if (taxa_are_rows(phy_species_rela)) {
  otu_mat <- t(otu_mat)
}

otu_df <- as.data.frame(otu_mat, check.names = FALSE)
otu_df$SampleID <- rownames(otu_df)

# Taxonomy table
tax_df <- as.data.frame(tax_table(phy_species_rela))
tax_df$taxa_id <- taxa_names(phy_species_rela)
tax_df$species_name <- as.character(tax_df$species)

# Keep selected taxa
tax_keep <- tax_df %>%
  filter(species_name %in% select_sig_taxa) %>%
  select(taxa_id, species_name)

# Check missing bacteria
missing_taxa <- setdiff(select_sig_taxa, tax_keep$species_name)
missing_taxa

# Long abundance table, collapse duplicated species if needed
abun_sel_long <- otu_df %>%
  select(SampleID, all_of(tax_keep$taxa_id)) %>%
  pivot_longer(cols = -SampleID,
    names_to = "taxa_id",
    values_to = "value") %>%
  left_join(tax_keep, by = "taxa_id") %>%
  group_by(SampleID, species = species_name) %>%
  summarise(value = sum(value, na.rm = TRUE),
    .groups = "drop")

# Wide table for your old coding style
abun_sel_df <- abun_sel_long %>%
  pivot_wider(names_from = species,
    values_from = value,
    values_fill = 0) %>%
  left_join( phy_species_rela_sample_df %>%
      select(SampleID,
        child_ID,
        Gut_maturation_Week_two,
        Time_new,
        AGE,
        oldchild01,
        catanddog,
        apartment_house,
        abbirth_mother,
        sex),by = "SampleID") %>%
  as.data.frame()
rownames(abun_sel_df) <- abun_sel_df$SampleID
View(abun_sel_df)
table(abun_sel_df$Time_new, useNA = "ifany")


##### One-bacterium trajectory plot + LMM using categorical Time_new #####
col_maturity <- c("High" = "#6a2675", "Low" = "#56b4e9")

covariates_use <- c("oldchild01",
                    "catanddog",
                    "apartment_house",
                    "abbirth_mother",
                    "sex")

run_one_bacterium_LMM_plot <- function(bac_name, time_keep=c("1w","1m","1y","4y","5y","6y"), 
                                       label_x = NULL, label_y = -0.7) {
  
  message("Running: ", bac_name)
  message("Using time points: ", paste(time_keep, collapse = ", "))
  
  # Select one bacterium
  bac_df <- abun_sel_df %>%
    select(all_of(c("SampleID",
                    "child_ID",
                    "Time_new",
                    "Gut_maturation_Week_two",
                    covariates_use,
                    bac_name))) %>%
    filter(Time_new %in% time_keep,
           !is.na(Gut_maturation_Week_two),
           !is.na(child_ID)) %>%
    mutate(variable = bac_name,
           value = .data[[bac_name]],
           log_value = log10(value + 1e-5),
           Time_new = factor(Time_new, levels = time_keep),
           Time_num = as.numeric(Time_new),
           Gut_maturation_Week_two = factor(Gut_maturation_Week_two, levels = c("Low", "High")),
           child_ID = factor(child_ID))
  
  # Keep covariates with variation
  cov_keep <- covariates_use[
    sapply(covariates_use, function(v) {
      x <- bac_df[[v]]
      length(unique(na.omit(x))) >= 2
    })
  ]
  
  # Complete cases for model
  model_vars <- c("log_value",
                  "Time_new",
                  "Gut_maturation_Week_two",
                  "child_ID",
                  cov_keep)
  
  bac_df2 <- bac_df %>%
    filter(complete.cases(across(all_of(model_vars)))) %>%
    droplevels()
  
  message("Samples used in LMM: ", nrow(bac_df2))
  message("Children used in LMM: ", length(unique(bac_df2$child_ID)))
  message("Covariates used: ", paste(cov_keep, collapse = ", "))
  
  if (length(unique(bac_df2$Time_new)) < 2) {
    stop("Less than 2 time points after filtering.")
  }
  
  if (length(unique(bac_df2$Gut_maturation_Week_two)) < 2) {
    stop("Less than 2 maturity groups after filtering.")
  }
  
  # Formula
  cov_text <- if (length(cov_keep) > 0) {
    paste0(" + ", paste(cov_keep, collapse = " + "))
  } else {
    ""
  }
  
  # Full model: maturity-specific pattern across categorical time points
  formula_full <- as.formula(
    paste0("log_value ~ Time_new * Gut_maturation_Week_two",
      cov_text,
      " + (1 | child_ID)"))
  
  # Null model 1: same time pattern, no maturity group
  formula_null_no_group <- as.formula(
    paste0("log_value ~ Time_new",
           cov_text,
           " + (1 | child_ID)"))
  
  # Null model 2: group effect but no interaction
  formula_null_no_interaction <- as.formula(
    paste0( "log_value ~ Time_new + Gut_maturation_Week_two",
            cov_text," + (1 | child_ID)"))
  
  # Fit models
  fit_full <- lmerTest::lmer(formula_full,
                             data = bac_df2,
                             REML = FALSE,
                             control = lme4::lmerControl(optimizer = "bobyqa"))
  
  fit_null_no_group <- lmerTest::lmer(formula_null_no_group,
                                      data = bac_df2,
                                      REML = FALSE,
                                      control = lme4::lmerControl(optimizer = "bobyqa"))
  
  fit_null_no_interaction <- lmerTest::lmer(formula_null_no_interaction,
                                            data = bac_df2,
                                            REML = FALSE,
                                            control = lme4::lmerControl(optimizer = "bobyqa"))
  
  # Type III ANOVA
  anova_full <- anova(fit_full, type = 3)
  
  # Model comparison 1:
  # Does adding maturity group and interaction improve the model?
  test_group_profile <- anova(fit_null_no_group, fit_full)
  
  # Model comparison 2:
  # Does interaction improve the model beyond main group effect?
  test_interaction <- anova(fit_null_no_interaction, fit_full)
  
  p_group_profile <- test_group_profile$`Pr(>Chisq)`[2]
  p_interaction <- test_interaction$`Pr(>Chisq)`[2]
  
  # High vs Low at each time point using emmeans
  emm_group_by_time <- emmeans::emmeans(fit_full,
    pairwise ~ Gut_maturation_Week_two | Time_new,
    adjust = "BH")
  
  contrast_by_time <- as.data.frame(emm_group_by_time$contrasts)
  
  print(summary(fit_full))
  print(anova_full)
  print(test_group_profile)
  print(test_interaction)
  print(contrast_by_time)
  
  label_text <- paste0(#"Group-time pattern: P = ", signif(p_group_profile, 3),
    "P = ", 
    signif(p_interaction, 3))
  
  if (is.null(label_x)) {
    label_x <- mean(range(bac_df2$Time_num, na.rm = TRUE))
  }
  
  # Summary observed data for plotting
  summary_df <- bac_df2 %>%
    group_by(variable, Time_new, Time_num, Gut_maturation_Week_two) %>%
    summarise(median_value = median(log_value, na.rm = TRUE),
              lower = quantile(log_value, 0.25, na.rm = TRUE),
              upper = quantile(log_value, 0.75, na.rm = TRUE),
              n = n(),
              .groups = "drop")
  
  # Plot: loess curve based on individual samples
  p <- bac_df2 %>%
    ggplot(aes(x = Time_num,y = log_value,
               fill = Gut_maturation_Week_two)) +
    theme_bw() +
    geom_smooth(aes(group = Gut_maturation_Week_two,
                    fill = Gut_maturation_Week_two),
                method = "loess",se = TRUE,linewidth = 1,
                color = "black",alpha = 0.25) +
    annotate("text",x = label_x,y = label_y,
             label = label_text,size = 3.6) +
    scale_x_continuous( breaks = seq_along(time_keep),
                        labels = time_keep,
                        expand = c(0, 0) ) +
    coord_cartesian(xlim = c(0.85, length(time_keep) + 0.15),
                    ylim = c(-5, 0)) +
    scale_fill_manual(values = col_maturity) +
    facet_wrap(~ variable) +
    theme(legend.position = "none",
          text = element_text(size = 9),
          strip.text = element_text(face = "italic")) +
    ylab("Log10 (Relative abundance)") +
    xlab("Age")
  
  return(
    list(
      data = bac_df2,
      summary = summary_df,
      fit_full = fit_full,
      fit_null_no_group = fit_null_no_group,
      fit_null_no_interaction = fit_null_no_interaction,
      anova_full = anova_full,
      test_group_profile = test_group_profile,
      test_interaction = test_interaction,
      contrast_by_time = contrast_by_time,
      p_group_profile = p_group_profile,
      p_interaction = p_interaction,
      plot = p
    )
  )
}



###### a. Early-transient taxa ######

######__ 1w1m __######
### Streptococcus mitis # p=1.94e-06
Streptococcus_mitis_res <- run_one_bacterium_LMM_plot(
  bac_name = "Streptococcus mitis",
  time_keep = c("1w", "1m", "1y"),
  label_x = 1.5,
  label_y = -0.7)
Streptococcus_mitis_res$plot
Streptococcus_mitis_res$p_group_profile
Streptococcus_mitis_res$p_interaction
Streptococcus_mitis_res$contrast_by_time

# Gemella sp. # p=27e-05
Gemella_sp_res <- run_one_bacterium_LMM_plot(
  bac_name = "Gemella sp.",
  time_keep = c("1w", "1m", "1y"),
  label_x = 1.5,
  label_y = -0.7)
Gemella_sp_res$plot
Gemella_sp_res$contrast_by_time

######__ 1w1y __######
# Actinomyces sp.
Actinomyces_sp_res <- run_one_bacterium_LMM_plot(
  bac_name = "Actinomyces sp.",
  time_keep = c("1w", "1m", "1y"),
  label_x = 1.5,
  label_y = -0.7)
Actinomyces_sp_res$plot
Actinomyces_sp_res$contrast_by_time


###### b. Early-persistent taxa ######

######__ 1w1m1y __######
# Alistipes onderdonkii
Alistipes_onderdonkii_res <- run_one_bacterium_LMM_plot(
  bac_name = "Alistipes onderdonkii",
  time_keep = c("1w", "1m", "1y"),
  label_x = 2,
  label_y = -0.7)
Alistipes_onderdonkii_res$plot
Alistipes_onderdonkii_res$contrast_by_time

# Blautia wexlerae # P =0.000125
Blautia_wexlerae_res <- run_one_bacterium_LMM_plot(
  bac_name = "Blautia wexlerae",
  time_keep = c("1w", "1m", "1y"),
  label_x = 2,
  label_y = -0.7)
Blautia_wexlerae_res$plot
Blautia_wexlerae_res$contrast_by_time

######__ 1m1y __######
# Lachnospiraceae sp. # P=0.00864
Lachnospiraceae_sp._res <- run_one_bacterium_LMM_plot(
  bac_name = "Lachnospiraceae sp.",
  time_keep = c("1w", "1m", "1y"),
  label_x = 2,
  label_y = -0.7)
Lachnospiraceae_sp._res$plot
Lachnospiraceae_sp._res$contrast_by_time

# Collinsella aerofaciens # P=0.0246
Collinsella_aerofaciens_res <- run_one_bacterium_LMM_plot(
  bac_name = "Collinsella aerofaciens",
  time_keep = c("1w", "1m", "1y"),
  label_x = 2,
  label_y = -0.7)
Collinsella_aerofaciens_res$plot
Collinsella_aerofaciens_res$contrast_by_time

# Faecalibacterium prausnitzii # P=0.000564
Faecalibacterium_prausnitzii_res <- run_one_bacterium_LMM_plot(
  bac_name = "Faecalibacterium prausnitzii",
  time_keep = c("1w", "1m", "1y","4y"),
  label_x = 2,
  label_y = -0.7)
Faecalibacterium_prausnitzii_res$plot
Faecalibacterium_prausnitzii_res$contrast_by_time


######__ 1w1m4y __######

# Bacteroides sp.# P =0.0315
Bacteroides_sp._res <- run_one_bacterium_LMM_plot(
  bac_name = "Bacteroides sp.",
  time_keep = c("1w", "1m", "1y","4y"),
  label_x = 2,
  label_y = -0.7)
Bacteroides_sp._res$plot
Bacteroides_sp._res$contrast_by_time

# Staphylococcus aureus # P = 1.58e-42
Staphylococcus_aureus_res <- run_one_bacterium_LMM_plot(
  bac_name = "Staphylococcus aureus",
  time_keep = c("1w", "1m", "1y","4y"),
  label_x = 2,
  label_y = -0.7)
Staphylococcus_aureus_res$plot
Staphylococcus_aureus_res$contrast_by_time

# Streptococcus salivarius # P=8.18e-06
Streptococcus_salivarius_res <- run_one_bacterium_LMM_plot(
  bac_name = "Streptococcus salivarius",
  time_keep = c("1w", "1m", "1y","4y"),
  label_x = 2,
  label_y = -0.7)
Streptococcus_salivarius_res$plot
Streptococcus_salivarius_res$contrast_by_time


######__ 1m4y __######

# Bacteroides uniformis
Bacteroides_uniformis_res <- run_one_bacterium_LMM_plot(
  bac_name = "Bacteroides uniformis",
  time_keep = c( "1m", "1y","4y"),
  label_x = 2,
  label_y = -0.7)
Bacteroides_uniformis_res$plot
Bacteroides_uniformis_res$contrast_by_time

# Bacteroides vulgatus
Bacteroides_vulgatus_res <- run_one_bacterium_LMM_plot(
  bac_name = "Bacteroides vulgatus",
  time_keep = c("1m", "1y","4y"),
  label_x = 2,
  label_y = -0.7)
Bacteroides_vulgatus_res$plot
Bacteroides_vulgatus_res$contrast_by_time

# Bilophila wadsworthia
Bilophila_wadsworthia_res <- run_one_bacterium_LMM_plot(
  bac_name = "Bilophila wadsworthia",
  time_keep = c( "1m", "1y","4y"),
  label_x = 2,
  label_y = -0.7)
Bilophila_wadsworthia_res$plot
Bilophila_wadsworthia_res$contrast_by_time

######__ 1m1y4y __######
# Parabacteroides merdae
Parabacteroides_merdae_res <- run_one_bacterium_LMM_plot(
  bac_name = "Parabacteroides merdae",
  time_keep = c( "1m", "1y","4y"),
  label_x = 2,
  label_y = -0.7)
Parabacteroides_merdae_res$plot
Parabacteroides_merdae_res$contrast_by_time


###### c. Late-persistent taxa ######
######__ 1y4y5y __######
# Odoribacter splanchnicus # P=2.08e-06
Odoribacter_splanchnicus_res <- run_one_bacterium_LMM_plot(
  bac_name = "Odoribacter splanchnicus",
  time_keep = c("1w", "1m","1y","4y","5y","6y"),
  label_x = 2,
  label_y = -0.7)
Odoribacter_splanchnicus_res$plot
Odoribacter_splanchnicus_res$contrast_by_time


######__ 1y4y5y6y __######
# Alistipes putredinis # P=4.15e-11
Alistipes_putredinis_res <- run_one_bacterium_LMM_plot(
  bac_name = "Alistipes putredinis",
  time_keep = c("1w", "1m", "1y","4y","5y","6y"),
  label_x = 2,
  label_y = -0.7)
Alistipes_putredinis_res$plot
Alistipes_putredinis_res$contrast_by_time

######__ 4y6y __######
# Christensenellaceae R-7 group sp. #P=5.07e-14
Christensenellaceae_R_7_group_sp._res <- run_one_bacterium_LMM_plot(
  bac_name = "Christensenellaceae R-7 group sp.",
  time_keep = c("1w", "1m", "1y","4y","5y","6y"),
  label_x = 2,
  label_y = -0.7)
Christensenellaceae_R_7_group_sp._res$plot
Christensenellaceae_R_7_group_sp._res$contrast_by_time

# Clostridia vadinBB60 group sp. # P=2.15e-12
Clostridia_vadinBB60_group_sp._res <- run_one_bacterium_LMM_plot(
  bac_name = "Clostridia vadinBB60 group sp.",
  time_keep = c("1w", "1m", "1y","4y","5y","6y"),
  label_x = 2,
  label_y = -0.7)
Clostridia_vadinBB60_group_sp._res$plot
Clostridia_vadinBB60_group_sp._res$contrast_by_time

# UCG-010 sp. # P=1.5e-15
UCG_010_sp._res <- run_one_bacterium_LMM_plot(
  bac_name = "UCG-010 sp.",
  time_keep = c("1w", "1m", "1y","4y","5y","6y"),
  label_x = 2,
  label_y = -0.7)
UCG_010_sp._res$plot
UCG_010_sp._res$contrast_by_time

######__ 4y5y6y __######
# UCG-002 sp. # P=3.5e-24
UCG_002_sp._res <- run_one_bacterium_LMM_plot(
  bac_name = "UCG-002 sp.",
  time_keep = c("1w", "1m", "1y","4y","5y","6y"),
  label_x = 2,
  label_y = -0.7)
UCG_002_sp._res$plot
UCG_002_sp._res$contrast_by_time



###### d. Save plots ######

### Mian figure: arrange plots together in AI

# Early-transient taxa
Streptococcus mitis # p=1.94e-06
Gemella sp. # p=27e-05
# Late-persistent taxa
Blautia wexlerae # P =0.000125
Lachnospiraceae sp. # P=0.00864
Faecalibacterium prausnitzii # P=0.000564
Bacteroides sp.# P =0.0315
Staphylococcus aureus # P = 1.58e-42
Streptococcus salivarius # P=8.18e-06
# Late-persistent taxa
Odoribacter splanchnicus # P=2.08e-06
Alistipes putredinis # P=4.15e-11

# Early-transient taxa
pdf("file_path/p_Streptococcus.mitis_1w1m1y.pdf",  height=2.8,width = 0.3+0.6*3)
Streptococcus_mitis_res$plot
dev.off()
pdf("file_path/p_Gemella.sp._1w1m1y.pdf",  height=2.8,width = 0.3+0.6*3)
Gemella_sp_res$plot 
dev.off()
# Early-persistent taxa
pdf("file_path/p_Blautia.wexlerae_1w1m1y.pdf",  height=2.8,width = 0.3+0.6*3)
Blautia_wexlerae_res$plot
dev.off()
pdf("file_path/p_Lachnospiraceae.sp._1w1m1y.pdf",  height=2.8,width = 0.3+0.6*3)
Lachnospiraceae_sp._res$plot
dev.off()
pdf("file_path/p_Faecalibacterium.prausnitzii_1w1m1y4y.pdf",  height=2.8,width = 0.3+0.6*4)
Faecalibacterium_prausnitzii_res$plot
dev.off()
pdf("file_path/p_Bacteroides.sp._1w1m1y4y.pdf",  height=2.8,width = 0.3+0.6*4)
Bacteroides_sp._res$plot
dev.off()
pdf("file_path/p_Staphylococcus.aureus_1w1m1y4y.pdf",  height=2.8,width = 0.3+0.6*4)
Staphylococcus_aureus_res$plot
dev.off()
pdf("file_path/p_Streptococcus.salivarius_1w1m1y4y.pdf",  height=2.8,width = 0.3+0.6*4)
Streptococcus_salivarius_res$plot
dev.off()
# Late-persistent taxa
pdf("file_path/p_Odoribacter.splanchnicus_1w1m1y4y5y6y.pdf",  height=2.8,width = 0.3+0.6*6)
Odoribacter_splanchnicus_res$plot
dev.off()
pdf("file_path/p_Alistipes.putredinis_1w1m1y4y5y6y.pdf",  height=2.8,width = 0.3+0.6*6)
Alistipes_putredinis_res$plot
dev.off()


### Supplementary figure: arrange plots together in AI

# Early-transient taxa
Actinomyces sp.
# Late-persistent taxa
Alistipes onderdonkii
Collinsella aerofaciens # P=0.0246
Bacteroides uniformis
Bacteroides vulgatus
Bilophila wadsworthia
Parabacteroides merdae
# Late-persistent taxa
Christensenellaceae R-7 group sp. #P=5.07e-14
Clostridia vadinBB60 group sp. # P=2.15e-12
UCG-010 sp. # P=1.5e-15
UCG-002 sp. # P=3.5e-24





#####  Figure 4B. Trajectories of representative species from each group #####
phy_species#

### relative abundance_representative figure 

# Get bacteria from upset plot
res_wil_1w1msig_name <- c("Staphylococcus aureus")
res_wil_1m1ysig_name <- c("Blautia faecis","Collinsella aerofaciens","Dialister invisus",
                          "Intestinibacter bartlettii","Lachnospiraceae sp.",
                          "Ruminococcus bromii")
res_wil_1w1m1y4ysig_name <- c("Faecalibacterium prausnitzii")
res_wil_1m1y4ysig_name <- c("Bacteroides uniformis","Parabacteroides merdae")
res_wil_1m1y4y6ysig_name <- c("Alistipes onderdonkii")
res_wil_1y4ysig_name <- c("Agathobacter rectalis","Alistipes finegoldii",
                          "Alistipes ihumii","Blautia stercoris","Coprococcus eutactus","Cuneatibacter caecimuris",
                          "Neglecta timonensis","Porphyromonas asaccharolytica","Roseburia faecis","Roseburia inulinivorans",
                          "Ruminococcus callidus","[Eubacterium] eligens group sp.","[Eubacterium] ruminantium group sp.",
                          "[Eubacterium] siraeum group sp.","Lachnospiraceae NK4A136 group sp.","Christensenella sp.",
                          "Colidextribacter sp.","Faecalibacterium sp.","Frisingicoccus sp.","Kineothrix sp.",
                          "Lachnospira sp.","Oscillibacter sp.","Ruminococcaceae sp.","Ruminococcus sp.",
                          "Sporobacter sp.")
res_wil_1y4y5y6ysig_name <- c("Alistipes putredinis","Odoribacter splanchnicus",
                              "UCG-002 sp.","UCG-010 sp.")
res_wil_1y4y6ysig_name <- c("[Eubacterium] coprostanoligenes group sp.",
                            "Alistipes indistinctus","Alistipes obesi","Alistipes shahii",
                            "Barnesiella intestinihominis","Christensenellaceae R-7 group sp.",
                            "Clostridia vadinBB60 group sp.","[Ruminococcus] gnavus group sp.")
res_wil_4y6ysig_name <- c("Murimonas intestini","Sutterella wadsworthensis","Parabacteroides distasonis",
                          "[Eubacterium] hallii group sp.")
res_wil_4y5y6ysig_name <- c("Anaerostipes hadrus")

select_sig_taxa <- c(res_wil_1w1msig_name,res_wil_1m1ysig_name,res_wil_1w1m1y4ysig_name,
                     res_wil_1m1y4ysig_name,res_wil_1m1y4y6ysig_name,res_wil_1y4ysig_name,
                     res_wil_1y4y5y6ysig_name,res_wil_1y4y6ysig_name,res_wil_4y6ysig_name,
                     res_wil_4y5y6ysig_name)
length(select_sig_taxa)

# Relative abundance table
sample_data(phy_species)$Time_new <- factor(sample_data(phy_species)$Time_new, levels = c("1w","1m","1y","4y","5y","6y"))
phy_species_rela <- transform_sample_counts(phy_species,function(x) x/sum(x))

# Sample data frame
phy_species_rela_sample_df <- data.frame(
  SampleID=get_variable(phy_species_rela,"SampleID"),
  Gut_maturation_Week_two=get_variable(phy_species_rela,"Gut_maturation_Week_two")%>%as.character,
  Time_new=get_variable(phy_species_rela,"Time_new"),
  AGE=get_variable(phy_species_rela,"AGE")
)
rownames(phy_species_rela_sample_df) <- phy_species_rela_sample_df$SampleID

# OTU data frame
phy_species_rela_select_otu_df <- phy_species_rela%>% 
  subset_taxa(species %in% select_sig_taxa) 
taxa_names(phy_species_rela_select_otu_df) <- tax_table(phy_species_rela_select_otu_df)[,7]

phy_species_rela_select_otu_df <- phy_species_rela_select_otu_df %>%
  otu_table %>%
  as("matrix")%>%
  t%>%
  as.data.frame

# Make training set (OTU_table + age)
abun_sel_df <- cbind(phy_species_rela_select_otu_df,
                     Gut_maturation_Week_two=phy_species_rela_sample_df[match(rownames(phy_species_rela_otu_df),
                                                                              rownames(phy_species_rela_sample_df)),"Gut_maturation_Week_two"],
                     Time_new=phy_species_rela_sample_df[match(rownames(phy_species_rela_otu_df),rownames(phy_species_rela_sample_df)),"Time_new"],
                     AGE=phy_species_rela_sample_df[match(rownames(phy_species_rela_otu_df),rownames(phy_species_rela_sample_df)),"AGE"])%>%
  as.data.frame()


###### a. Early-transient taxa ######
abun_sel_df #
abun_sel_df$Time_new <- factor(abun_sel_df$Time_new, levels = c("1w","1m","1y","4y","5y","6y"))

######__ 1w1m __######
bac_1w1m <- abun_sel_df %>%
  select(all_of(c("Time_new","Gut_maturation_Week_two","Staphylococcus aureus")))

### Staphylococcus aureus
bac_df <- bac_1w1m %>% 
  subset(Time_new %in% c("1w","1m","1y")&!is.na(Gut_maturation_Week_two))%>%
  melt %>% 
  as.data.frame()

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1w` = 1,
                                `1m` = 2,
                                `1y` = 3) %>% as.numeric()

# Plot
col_maturity <- c("High"="#6a2675", "Low"="#56b4e9")

Staphylococcus.aureus_1w1m <- bac_df %>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", 
                     method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c("1w", "1m", "1y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 3.15)) +
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Staphylococcus.aureus_1w1m

## Add linear model
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Staphylococcus aureus")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Staphylococcus.aureus_1w1m_final <- Staphylococcus.aureus_1w1m +
  annotate("text", x = 2, y = -2, label = label_text, size = 4)
Staphylococcus.aureus_1w1m_final


######__ 1m1y __######
bac_1m1y <- abun_sel_df %>%
  select(all_of(c("Time_new","Gut_maturation_Week_two",
                  "Blautia faecis","Collinsella aerofaciens","Dialister invisus",
                  "Intestinibacter bartlettii","Lachnospiraceae sp.","Ruminococcus bromii")))

### Intestinibacter bartlettii
bac_df <- bac_1m1y %>% 
  subset(Time_new %in% c("1m","1y")&!is.na(Gut_maturation_Week_two))%>%
  melt %>% 
  as.data.frame()

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1m` = 1,
                                `1y` = 2) %>% as.numeric()

# Plot
Intestinibacter.bartlettii_1m1y <- bac_df %>%
  subset(variable=="Intestinibacter bartlettii")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", 
                     method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c( "1m", "1y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 2.15)) + 
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Intestinibacter.bartlettii_1m1y

## Add linear model
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Intestinibacter bartlettii")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Intestinibacter.bartlettii_1m1y_final <- Intestinibacter.bartlettii_1m1y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
Intestinibacter.bartlettii_1m1y_final


### Collinsella aerofaciens
bac_df <- bac_1m1y %>% 
  subset(Time_new %in% c("1m","1y")&!is.na(Gut_maturation_Week_two))%>%
  melt %>% 
  as.data.frame()

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                 `1m` = 1,
                                 `1y` = 2) %>% as.numeric()

# Plot
Collinsella.aerofaciens_1m1y <- bac_df %>%
  subset(variable=="Collinsella aerofaciens")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", 
                     method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c( "1m", "1y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 2.15)) +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Collinsella.aerofaciens_1m1y

## Add linear model
bac_df <- bac_df %>%
  subset(variable=="Collinsella aerofaciens")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Collinsella.aerofaciens_1m1y_final <-  Collinsella.aerofaciens_1m1y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
Collinsella.aerofaciens_1m1y_final


### Dialister invisus
bac_df <- bac_1m1y %>% 
  subset(Time_new %in% c("1m","1y")&!is.na(Gut_maturation_Week_two))%>%
  melt %>% 
  as.data.frame()

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                 `1m` = 1,
                                 `1y` = 2) %>% as.numeric()
Dialister.invisus_1m1y <- bac_df %>%
  subset(variable=="Dialister invisus")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", 
                     method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("1m", "1y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 2.15)) +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Dialister.invisus_1m1y

## Add linear model
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Dialister invisus")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Dialister.invisus_1m1y_final <- Dialister.invisus_1m1y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
Dialister.invisus_1m1y_final


### Lachnospiraceae sp.
bac_df <- bac_1m1y %>% 
  subset(Time_new %in% c("1m","1y")&!is.na(Gut_maturation_Week_two))%>%
  melt %>% 
  as.data.frame()

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1m` = 1,
                                `1y` = 2) %>% as.numeric()

Lachnospiraceae.sp._1m1y <- bac_df %>%
  subset(variable=="Lachnospiraceae sp.")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2),
                     labels = c( "1m", "1y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 2.15)) + 
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Lachnospiraceae.sp._1m1y

## Add linear model
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Lachnospiraceae sp.")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Lachnospiraceae.sp._1m1y_final <- Lachnospiraceae.sp._1m1y +
  annotate("text", x = 1.8, y = -2, label = label_text, size = 4)
Lachnospiraceae.sp._1m1y_final




###### b. Early-persistent taxa ######

######__ 1w1m1y4y __######
bac_1w1m1y4y <- abun_sel_df %>%
  select(all_of(c("Time_new","Gut_maturation_Week_two",
                  "Faecalibacterium prausnitzii")))

### Faecalibacterium prausnitzii
bac_df <- bac_1w1m1y4y %>% 
  subset(Time_new %in% c("1w","1m","1y","4y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>% 
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                 `1w` = 1,
                                 `1m` = 2,
                                 `1y` = 3,
                                 `4y` = 4) %>% as.numeric()

Faecalibacterium.prausnitzii_1w1m1y4y <- bac_df %>%
  subset(variable=="Faecalibacterium prausnitzii")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two), method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3,4),
                     labels = c("1w", "1m", "1y","4y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 4.15)) +
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Faecalibacterium.prausnitzii_1w1m1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Faecalibacterium prausnitzii")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Faecalibacterium.prausnitzii_1w1m1y4y_final <- Faecalibacterium.prausnitzii_1w1m1y4y +
  annotate("text", x = 2, y = -2, label = label_text, size = 4)
Faecalibacterium.prausnitzii_1w1m1y4y_final 


######__ 1m1y4y __######

### Bacteroides uniformis
bac_df <- bac_1m1y4y %>% 
  subset(Time_new %in% c("1m","1y","4y")&!is.na(Gut_maturation_Week_two))%>%
  melt %>% 
  as.data.frame()

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                 `1m` = 1,
                                 `1y` = 2,
                                 `4y` = 3) %>% as.numeric()

Bacteroides.uniformis_1m1y4y <- bac_df %>%
  subset(variable=="Bacteroides uniformis")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two), method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c( "1m", "1y","4y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 3.15)) +
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Bacteroides.uniformis_1m1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Bacteroides uniformis")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Bacteroides.uniformis_1m1y4y_final <- Bacteroides.uniformis_1m1y4y +
  annotate("text", x = 1.8, y = -2, label = label_text, size = 4)
Bacteroides.uniformis_1m1y4y_final


### Parabacteroides merdae
bac_df <- bac_1m1y4y %>% 
  subset(Time_new %in% c("1m","1y","4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt %>% 
  as.data.frame()

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1m` = 1,
                                `1y` = 2,
                                `4y` = 3,
                                `5y` = 4,
                                `6y` = 5) %>% as.numeric()

Parabacteroides.merdae_1m1y4y <- bac_df %>%
  subset(variable=="Parabacteroides merdae")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3,4,5),
                     labels = c("1m", "1y", "4y","5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 5.15)) + 
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Parabacteroides.merdae_1m1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Parabacteroides merdae")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Parabacteroides.merdae_1m1y4y_final <- Parabacteroides.merdae_1m1y4y +
  annotate("text", x = 2.5, y = -2, label = label_text, size = 4)
Parabacteroides.merdae_1m1y4y_final


######__ 1m1y4y6y __######
bac_1m1y4y6y <- abun_sel_df %>%select(all_of(c("Time_new","Gut_maturation_Week_two",
                                                 "Alistipes onderdonkii")))

### Alistipes onderdonkii
bac_df <- bac_1m1y4y6y %>% 
  subset(Time_new %in% c("1m","1y","4y","6y","5y","1w")&!is.na(Gut_maturation_Week_two))%>%
  melt()%>% as.data.frame()

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1w` = 1,
                                `1m` = 2,
                                `1y` = 3,
                                `4y` = 4,
                                `5y` = 5,
                                `6y` = 6) %>% as.numeric()

Alistipes.onderdonkii_1m1y4y6y <- bac_df %>%
  subset(variable=="Alistipes onderdonkii")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", 
                     method = "wilcox.test", color="red" ) +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  scale_x_continuous(breaks = c(1, 2, 3,4,5,6),
                     labels = c("1w", "1m", "1y","4y","5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 6.15)) + 
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Alistipes.onderdonkii_1m1y4y6y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Alistipes onderdonkii")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Alistipes.onderdonkii_1m1y4y6y_final <- Alistipes.onderdonkii_1m1y4y6y +
  annotate("text", x = 3, y = -2, label = label_text, size = 4)
Alistipes.onderdonkii_1m1y4y6y_final


###### c. Late-persistent taxa ######
######__ 1y4y __######
bac_1y4y <- abun_sel_df %>%
  select(all_of(c("Time_new","Gut_maturation_Week_two",
                  "Agathobacter rectalis","Alistipes finegoldii","Alistipes ihumii",
                  "Blautia stercoris","Coprococcus eutactus","Cuneatibacter caecimuris",
                  "Neglecta timonensis","Porphyromonas asaccharolytica","Roseburia faecis",
                  "Roseburia inulinivorans","Ruminococcus callidus",
                  "[Eubacterium] eligens group sp.","[Eubacterium] ruminantium group sp.",
                  "[Eubacterium] siraeum group sp.","Lachnospiraceae NK4A136 group sp.",
                  "Christensenella sp.","Colidextribacter sp.","Faecalibacterium sp.",
                  "Frisingicoccus sp.","Kineothrix sp.","Lachnospira sp.","Oscillibacter sp.",
                  "Ruminococcaceae sp.","Ruminococcus sp.","Sporobacter sp.")))

### Agathobacter rectalis: 1y4y5y6y
bac_df <- bac_1y4y %>% 
  subset(Time_new %in% c("1y","4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1y` = 1,
                                `4y` = 2,
                                `5y` = 3,
                                `6y` = 4) %>% as.numeric()

Agathobacter.rectalis_1y4y <- bac_df %>%
  subset(variable=="Agathobacter rectalis")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3,4),
                     labels = c("1y", "4y", "5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 4.15)) +
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Agathobacter.rectalis_1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Agathobacter rectalis")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Agathobacter.rectalis_1y4y_final <- Agathobacter.rectalis_1y4y +
  annotate("text", x = 2.5, y = -2, label = label_text, size = 4)
Agathobacter.rectalis_1y4y_final



### Alistipes finegoldii: 1y4y5y6y
bac_df <- bac_1y4y %>% 
  subset(Time_new %in% c("1y","4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                 `1y` = 1,
                                 `4y` = 2,
                                 `5y` = 3,
                                 `6y` = 4) %>% as.numeric()

Alistipes.finegoldii_1y4y <- bac_df %>%
  subset(variable=="Alistipes finegoldii")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", 
                     method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  scale_x_continuous(breaks = c(1, 2, 3,4),
                     labels = c("1y", "4y", "5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 4.15)) + 
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Alistipes.finegoldii_1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Alistipes finegoldii")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Alistipes.finegoldii_1y4y_final <- Alistipes.finegoldii_1y4y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
Alistipes.finegoldii_1y4y_final


### Roseburia inulinivorans: 1y4y5y6y
bac_df <- bac_1y4y %>% 
  subset(Time_new %in% c("1y","4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1y` = 1,
                                `4y` = 2,
                                `5y` = 3,
                                `6y` = 4) %>% as.numeric()

Roseburia.inulinivorans_1y4y <- bac_df %>%
  subset(variable=="Roseburia inulinivorans")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous( breaks = c(1, 2, 3,4),
                      labels = c("1y", "4y", "5y","6y"),
                      expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 4.15)) +
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Roseburia.inulinivorans_1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Roseburia inulinivorans")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Roseburia.inulinivorans_1y4y_final <- Roseburia.inulinivorans_1y4y +
  annotate("text", x = 1.8, y = -2, label = label_text, size = 4)
Roseburia.inulinivorans_1y4y_final


### Ruminococcus callidus: 1y4y5y6y
bac_df <- bac_1y4y %>% 
  subset(Time_new %in% c("1y","4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1y` = 1,
                                `4y` = 2,
                                `5y` = 3,
                                `6y` = 4) %>% as.numeric()

Ruminococcus.callidus_1y4y <- bac_df %>%
  subset(variable=="Ruminococcus callidus")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", 
                     method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two), method = "loess") +
  scale_x_continuous( breaks = c(1, 2, 3,4),
                      labels = c("1y", "4y", "5y","6y"),
                      expand = c(0, 0) )+
  coord_cartesian(xlim = c(0.85, 4.15)) + 
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Ruminococcus.callidus_1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Ruminococcus callidus")%>%
  mutate(log_value = log10(value + 1e-5))

# Plot
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Ruminococcus.callidus_1y4y_final <- Ruminococcus.callidus_1y4y +
  annotate("text", x = 2.5, y = -2, label = label_text, size = 4)
Ruminococcus.callidus_1y4y_final


### Sporobacter sp.: 1y4y5y6y
bac_df <- bac_1y4y %>% 
  subset(Time_new %in% c("1y","4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1y` = 1,
                                `4y` = 2,
                                `5y` = 3,
                                `6y` = 4) %>% as.numeric()

Sporobacter.sp._1y4y <- bac_df %>%
  subset(variable=="Sporobacter sp.")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", 
                     method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  scale_x_continuous( breaks = c(1, 2, 3,4),
                      labels = c("1y", "4y", "5y","6y"),
                      expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 4.15)) +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Sporobacter.sp._1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Sporobacter sp.")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Sporobacter.sp._1y4y_final <- Sporobacter.sp._1y4y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
Sporobacter.sp._1y4y_final


### Coprococcus eutactus: 1y4y5y
bac_df <- bac_1y4y %>% 
  subset(Time_new %in% c("1y","4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1y` = 1,
                                `4y` = 2,
                                `5y` = 3,
                                `6y` = 4) %>% as.numeric()

Coprococcus.eutactus_1y4y <- bac_df %>%
  subset(variable=="Coprococcus eutactus")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3,4),
                     labels = c("1y", "4y", "5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 4.15)) + 
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Coprococcus.eutactus_1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Coprococcus eutactus")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Coprococcus.eutactus_1y4y_final <- Coprococcus.eutactus_1y4y +
  annotate("text", x = 1.8, y = -2, label = label_text, size = 4)
Coprococcus.eutactus_1y4y_final


### Neglecta timonensis: 1y4y
bac_df <- bac_1y4y %>% 
  subset(Time_new %in% c("1y","4y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                 `1y` = 1,
                                 `4y` = 2) %>% as.numeric()

Neglecta.timonensis_1y4y <- bac_df %>%
  subset(variable=="Neglecta timonensis")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", 
                     method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("1y", "4y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 2.15)) + 
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Neglecta.timonensis_1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Neglecta timonensis")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Neglecta.timonensis_1y4y_final <- Neglecta.timonensis_1y4y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
Neglecta.timonensis_1y4y_final


### Faecalibacterium sp.: 1y4y
bac_df <- bac_1y4y %>% 
  subset(Time_new %in% c("1y","4y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                 `1y` = 1,
                                 `4y` = 2) %>% as.numeric()

Faecalibacterium.sp._1y4y <- bac_df %>%
  subset(variable=="Faecalibacterium sp.")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", 
                     method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("1y", "4y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 2.15)) + 
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Faecalibacterium.sp._1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Faecalibacterium sp.")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Faecalibacterium.sp._1y4y_final <- Faecalibacterium.sp._1y4y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
Faecalibacterium.sp._1y4y_final


### Oscillibacter sp.: 1y4y
bac_df <- bac_1y4y %>% 
  subset(Time_new %in% c("1y","4y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1y` = 1,
                                `4y` = 2) %>% as.numeric()

Oscillibacter.sp._1y4y <- bac_df %>%
  subset(variable=="Oscillibacter sp.")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", 
                     method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("1y", "4y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 2.15)) + 
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Oscillibacter.sp._1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Oscillibacter sp.")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Oscillibacter.sp._1y4y_final <- Oscillibacter.sp._1y4y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
Oscillibacter.sp._1y4y_final


### Ruminococcaceae sp.: 1y4y 
bac_df <- bac_1y4y %>% 
  subset(Time_new %in% c("1y","4y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1y` = 1,
                                `4y` = 2) %>% as.numeric()

Ruminococcaceae.sp._1y4y <- bac_df %>%
  subset(variable=="Ruminococcaceae sp.")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", 
                     method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("1y", "4y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 2.15)) +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Ruminococcaceae.sp._1y4y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Ruminococcaceae sp.")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Ruminococcaceae.sp._1y4y_final <- Ruminococcaceae.sp._1y4y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
Ruminococcaceae.sp._1y4y_final




######__ 1y4y5y6y __######
bac_1y4y5y6y <- abun_sel_df %>%
  select(all_of(c("Time_new","Gut_maturation_Week_two",
                  "Alistipes putredinis","Odoribacter splanchnicus",
                  "UCG-002 sp.","UCG-010 sp.")))

### Alistipes putredinis: 1y4y5y6y
bac_df <- bac_1y4y5y6y %>% 
  subset(Time_new %in% c("1y","4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1y` = 1,
                                `4y` = 2,
                                `5y` = 3,
                                `6y` = 4) %>% as.numeric()

Alistipes.putredinis_1y4y5y6y <- bac_df %>%
  subset(variable=="Alistipes putredinis")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3,4),
                     labels = c("1y", "4y", "5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 4.15)) +
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Alistipes.putredinis_1y4y5y6y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Alistipes putredinis")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Alistipes.putredinis_1y4y5y6y_final <- Alistipes.putredinis_1y4y5y6y +
  annotate("text", x = 1.8, y = -2, label = label_text, size = 4)
Alistipes.putredinis_1y4y5y6y_final


### Odoribacter splanchnicus: 1y4y5y6y
bac_df <- bac_1y4y5y6y %>% 
  subset(Time_new %in% c("1y","4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                 `1y` = 1,
                                 `4y` = 2,
                                 `5y` = 3,
                                 `6y` = 4) %>% as.numeric()

Odoribacter.splanchnicus_1y4y5y6y <- bac_df %>%
  subset(variable=="Odoribacter splanchnicus")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3,4),
                     labels = c("1y", "4y", "5y","6y"),
                     expand = c(0, 0) )+
  coord_cartesian(xlim = c(0.85, 4.15)) + 
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Odoribacter.splanchnicus_1y4y5y6y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Odoribacter splanchnicus")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Odoribacter.splanchnicus_1y4y5y6y_final <- Odoribacter.splanchnicus_1y4y5y6y +
  annotate("text", x = 2.5, y = -2, label = label_text, size = 4)
Odoribacter.splanchnicus_1y4y5y6y_final


### UCG_002_sp.: 1y4y5y6y
tax_table(phy_species)[tax_table(phy_species)[,7]=="UCG-002 sp.",]
bac_1y4y5y6y <- abun_sel_df %>%select(all_of(c("Time_new","Gut_maturation_Week_two",
                                                 "Alistipes putredinis","Odoribacter splanchnicus",
                                                 "UCG-002 sp.","UCG-010 sp.")))
bac_df <- bac_1y4y5y6y %>% 
  subset(Time_new %in% c("1y","4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1y` = 1,
                                `4y` = 2,
                                `5y` = 3,
                                `6y` = 4) %>% as.numeric()

UCG_002_sp._1y4y5y6y <- bac_df %>%
  subset(variable=="UCG-002 sp.")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3,4),
                     labels = c("1y", "4y", "5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 4.15)) + 
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
UCG_002_sp._1y4y5y6y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="UCG-002 sp.")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
UCG_002_sp._1y4y5y6y_final <- UCG_002_sp._1y4y5y6y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
UCG_002_sp._1y4y5y6y_final


### UCG_010_sp.: 1y4y5y6y
bac_df <- bac_1y4y5y6y %>% 
  subset(Time_new %in% c("1y","4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                 `1y` = 1,
                                 `4y` = 2,
                                 `5y` = 3,
                                 `6y` = 4) %>% as.numeric()

UCG_010_sp._1y4y5y6y <- bac_df %>%
  subset(variable=="UCG-010 sp.")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3,4),
                     labels = c("1y", "4y", "5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 4.15)) + 
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
UCG_010_sp._1y4y5y6y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="UCG-010 sp.")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
UCG_010_sp._1y4y5y6y_final <- UCG_010_sp._1y4y5y6y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
UCG_010_sp._1y4y5y6y_final


######__  1y4y6y __######
bac_1y4y6y <- abun_sel_df %>%
  select(all_of(c("Time_new","Gut_maturation_Week_two",
                  "[Eubacterium] coprostanoligenes group sp.","Alistipes indistinctus",
                  "Alistipes obesi","Alistipes shahii","Barnesiella intestinihominis",
                  "Christensenellaceae R-7 group sp.","Clostridia vadinBB60 group sp.",
                  "[Ruminococcus] gnavus group sp.")))

### Christensenellaceae_R_7_group_sp.: 1y4y6y
bac_df <- bac_1y4y6y %>% 
  subset(Time_new %in% c("1y","4y","6y","5y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `1y` = 1,
                                `4y` = 2,
                                `5y` = 3,
                                `6y` = 4) %>% as.numeric()

Christensenellaceae_R_7_group_sp._1y4y6y <- bac_df %>%
  subset(variable=="Christensenellaceae R-7 group sp.")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous( breaks = c(1, 2, 3,4),
                      labels = c("1y", "4y", "5y","6y"),
                      expand = c(0, 0) )+
  coord_cartesian(xlim = c(0.85, 4.15)) + 
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Christensenellaceae_R_7_group_sp._1y4y6y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Christensenellaceae R-7 group sp.")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Christensenellaceae_R_7_group_sp._1y4y6y_final <- Christensenellaceae_R_7_group_sp._1y4y6y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
Christensenellaceae_R_7_group_sp._1y4y6y_final


### Ruminococcus_gnavus_group_sp.: 1y4y6y
bac_df <- bac_1y4y6y %>% 
  subset(Time_new %in% c("1y","4y","6y","5y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                 `1y` = 1,
                                 `4y` = 2,
                                 `5y` = 3,
                                 `6y` = 4) %>% as.numeric()

Ruminococcus_gnavus_group_sp._1y4y6y <- bac_df %>%
  subset(variable=="[Ruminococcus] gnavus group sp.")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3,4),
                     labels = c("1y", "4y", "5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 4.15)) +
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Ruminococcus_gnavus_group_sp._1y4y6y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="[Ruminococcus] gnavus group sp.")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Ruminococcus_gnavus_group_sp._1y4y6y_final <- Ruminococcus_gnavus_group_sp._1y4y6y +
  annotate("text", x = 2.5, y = -2, label = label_text, size = 4)
Ruminococcus_gnavus_group_sp._1y4y6y_final


### Eubacterium_coprostanoligenes_group_sp.: 1y4y6y
bac_df <- bac_1y4y6y %>% 
  subset(Time_new %in% c("1y","4y","6y","5y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                 `1y` = 1,
                                 `4y` = 2,
                                 `5y` = 3,
                                 `6y` = 4) %>% as.numeric()

Eubacterium_coprostanoligenes_group_sp._1y4y6y <- bac_df %>%
  subset(variable=="[Eubacterium] coprostanoligenes group sp.")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two), method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous( breaks = c(1, 2, 3,4),
                      labels = c("1y", "4y", "5y","6y"),
                      expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 4.15)) + 
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Eubacterium_coprostanoligenes_group_sp._1y4y6y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="[Eubacterium] coprostanoligenes group sp.")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Eubacterium_coprostanoligenes_group_sp._1y4y6y_final <- Eubacterium_coprostanoligenes_group_sp._1y4y6y +
  annotate("text", x = 1.5, y = -2, label = label_text, size = 4)
Eubacterium_coprostanoligenes_group_sp._1y4y6y_final


### Alistipes shahii: 1y4y6y
bac_df <- bac_1y4y6y %>% 
  subset(Time_new %in% c("1y","4y","6y","5y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                 `1y` = 1,
                                 `4y` = 2,
                                 `5y` = 3,
                                 `6y` = 4) %>% as.numeric()

Alistipes.shahii_1y4y6y <- bac_df %>%
  subset(variable=="Alistipes shahii")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3,4),
                     labels = c("1y", "4y", "5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 4.15)) + 
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Alistipes.shahii_1y4y6y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Alistipes shahii")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Alistipes.shahii_1y4y6y_final <-Alistipes.shahii_1y4y6y +
  annotate("text", x = 2.5, y = -2, label = label_text, size = 4)
Alistipes.shahii_1y4y6y_final



######__ 4y6y __######
bac_4y6y <- abun_sel_df %>%
  select(all_of(c("Time_new","Gut_maturation_Week_two",
                  "Murimonas intestini","Sutterella wadsworthensis",
                  "Parabacteroides distasonis","[Eubacterium] hallii group sp.")))

### Eubacterium_hallii_group_sp.: 4y6y
bac_df <- bac_4y6y %>% 
  subset(Time_new %in% c("4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `4y` = 1,
                                `5y` = 2,
                                `6y` = 3) %>% as.numeric()

Eubacterium_hallii_group_sp._4y6y <- bac_df %>%
  subset(variable=="[Eubacterium] hallii group sp.")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c("4y", "5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 3.15)) +
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Eubacterium_hallii_group_sp._4y6y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="[Eubacterium] hallii group sp.")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Eubacterium_hallii_group_sp._4y6y_final <-Eubacterium_hallii_group_sp._4y6y+
  annotate("text", x = 1.8, y = -2, label = label_text, size = 4)
Eubacterium_hallii_group_sp._4y6y_final


### Parabacteroides distasonis: 4y6y
bac_df <- bac_4y6y %>% 
  subset(Time_new %in% c("4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `4y` = 1,
                                `5y` = 2,
                                `6y` = 3) %>% as.numeric()

Parabacteroides.distasonis_4y6y <- bac_df %>%
  subset(variable=="Parabacteroides distasonis")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c("4y", "5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 3.15)) +
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Parabacteroides.distasonis_4y6y

## Add linear model # Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Parabacteroides distasonis")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Parabacteroides.distasonis_4y6y_final <-Parabacteroides.distasonis_4y6y+
  annotate("text", x = 1.8, y = -2, label = label_text, size = 4)
Parabacteroides.distasonis_4y6y_final




######__ 4y5y6y __######
bac_4y5y6y <- abun_sel_df %>%
  select(all_of(c("Time_new","Gut_maturation_Week_two", "Anaerostipes hadrus")))

### Anaerostipes hadrus_4y5y6y
bac_df <- bac_4y5y6y %>% 
  subset(Time_new %in% c("4y","5y","6y")&!is.na(Gut_maturation_Week_two))%>%
  melt%>%
  as.data.frame

# Add time numerical variable (e.g., 1y=1, 4y=4, 5y=5, 6y=6)
bac_df$Time_num <- recode_factor(bac_df$Time_new,
                                `4y` = 1,
                                `5y` = 2,
                                `6y` = 3) %>% as.numeric()

Anaerostipes.hadrus_4y5y6y <- bac_df %>%
  subset(variable=="Anaerostipes hadrus")%>%
  ggplot(aes(x=Time_num,y=log10(value+0.00001), fill=Gut_maturation_Week_two))+
  theme_bw()+
  stat_compare_means(aes(group = Gut_maturation_Week_two), label = "p.signif", method = "wilcox.test", color="red") +
  geom_smooth(aes(x = Time_num, y = log10(value + 0.00001),
                  fill=Gut_maturation_Week_two,
                  color=Gut_maturation_Week_two),method = "loess") +
  theme(legend.position = "none")+
  ylab("Log10 (Relative abundance)")+
  xlab("Age")+
  ylim(-5, 0)+
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c( "4y", "5y","6y"),
                     expand = c(0, 0))+
  coord_cartesian(xlim = c(0.85, 3.15)) +
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  facet_wrap(~variable)+
  theme(text = element_text(size=11))
Anaerostipes.hadrus_4y5y6y

## Add linear model 
# Add log transformation column
bac_df <- bac_df %>%
  subset(variable=="Anaerostipes hadrus")%>%
  mutate(log_value = log10(value + 1e-5))

# Model
fit <- lm(log_value ~ Time_new * Gut_maturation_Week_two, data = bac_df)
summary(fit)
anova(fit)

# Extract interaction term P-value
interaction_p <- anova(fit)[["Pr(>F)"]][3]
label_text <- paste0("Interaction effect: P = ", signif(interaction_p, 3))

# Plot
Anaerostipes.hadrus_4y5y6y_final <-Anaerostipes.hadrus_4y5y6y +
  annotate("text", x = 1.8, y = -2, label = label_text, size = 4)
Anaerostipes.hadrus_4y5y6y_final



###### d. Save plots ######

### Mian figure: arrange plots together in AI
# Early-transient taxa
Staphylococcus.aureus_1w1m_final
Lachnospiraceae.sp._1m1y_final #can ferment dietary fiber and other complex polysaccharides
# Early-persistent taxa
Faecalibacterium.prausnitzii_1w1m1y4y_final #probiotics
Bacteroides.uniformis_1m1y4y_final #probiotics #lm not significant
# Late-persistent taxa
Agathobacter.rectalis_1y4y_final #probiotics
Roseburia.inulinivorans_1y4y_final #probiotics
Odoribacter.splanchnicus_1y4y5y6y_final #lm not significant
Alistipes.shahii_1y4y6y_final #probiotics
Anaerostipes.hadrus_4y5y6y_final #lm not significant


### Supplementary figure: arrange plots together in AI
# Early-transient taxa
Collinsella.aerofaciens_1m1y_final #probiotics
Intestinibacter.bartlettii_1m1y_final #common bacteria
# Early-persistent taxa
Parabacteroides.merdae_1m1y4y_final #common bacteria
# Late-persistent taxa
Faecalibacterium.sp._1y4y_final #probiotics
Oscillibacter.sp._1y4y_final #probiotics
Ruminococcaceae.sp._1y4y_final #probiotics
Alistipes.finegoldii_1y4y_final # unclear, potentially beneficial
Coprococcus.eutactus_1y4y_final # unclear, potentially beneficial
##### Ruminococcus.callidus_1y4y_final # unclear, potentially beneficial
Neglecta.timonensis_1y4y_final # rare species
Sporobacter.sp._1y4y_final # little research


### Other candidate bacteria (do not show in main/supplementary figure)
## Late-persistent taxa
Alistipes.putredinis_1y4y5y6y_final #lm not significant
Coprococcus.eutactus_1y4y_final #probiotics #lm not significant 
Parabacteroides.distasonis_4y6y_final #probiotics
UCG_002_sp._1y4y5y6y_final
Eubacterium_hallii_group_sp._4y6y_final#probiotics
Eubacterium_coprostanoligenes_group_sp._1y4y6y_final # This bacteria has expression levels but not plotted due to unknown early or late persistent taxa group
Christensenellaceae_R_7_group_sp._1y4y6y_final
Ruminococcus_gnavus_group_sp._1y4y6y_final



