########################### Gut microbial maturation - Result 2 ################################################

#### Figure 2. Development of gut microbiota and environmental influences in early childhood ####

# Load required packages
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(randomForest)
library(gtools)
library(rstatix)
library(tidyr)
library(tibble)
library(rcompanion)
library(scales)
library(grid)
library(reshape2)
library(writexl)
library(geepack)
library(broom)
library(ggalluvial)
library(forestplot)

# Load phyloseq object
load("phy_feces_clean.RData")
phy_feces_clean#




#####  Figure 2A. Correlation between predicted microbiota age and chronological age #####

###### a. Random forest model - all ASVs ######
phy_all6 <- prune_taxa(taxa_sums(phy_feces_clean) > 0, phy_feces_clean)
all6_sample_df <- data.frame(
  SampleID = get_variable(phy_all6, "SampleID"),
  age = get_variable(phy_all6, "AGE") %>% as.character %>% as.numeric
)
rownames(all6_sample_df) <- all6_sample_df$SampleID

# Extract OTU table and preprocess
all6_otu_df <- transform_sample_counts(phy_all6, function(x) log(x/sum(x) + 0.0001)) %>%
  otu_table() %>%
  as("matrix") %>%
  t() %>%
  as.data.frame() %>%
  scale()# scale is generic function whose default method centers and/or scales the columns of a numeric matrix.

# Format column names for random forest compatibility
colnames(all6_otu_df) <- paste0("x", colnames(all6_otu_df))

# Prepare training data
train_all6 <- cbind(all6_otu_df,
                    age = all6_sample_df[match(rownames(all6_otu_df),
                                               rownames(all6_sample_df)), "age"])
train_all6 <- train_all6[!is.na(train_all6[, "age"]), ]
train_all6 <- as.data.frame(train_all6)
save(train_all6, file = "train_all6.RData")

## Train random forest model
set.seed(33)
rf_trained_all6_ASV_10000 <- randomForest(
  age ~ ., 
  data = train_all6,
  ntree = 10000,
  keep.forest = TRUE,
  importance = TRUE
)
#save(rf_trained_all6_ASV_10000, file = "rf_trained_all6_ASV_10000.RData")

print(rf_trained_all6_ASV_10000)
# Result note: % var explained_91.94%

# Predict microbiota age for all samples
colnames(all6_otu_df) <- make.names(colnames(all6_otu_df))
pred_all_all6 <- predict(rf_trained_all6_ASV_10000, all6_otu_df)


###### b. MAZ score calculation ######
# Convert age to weeks
sample_data(phy_feces_clean)$Week <- paste0(as.integer(sample_data(phy_feces_clean)$AGE/7), "w")
sample_data(phy_feces_clean)$Week <- factor(
  sample_data(phy_feces_clean)$Week,
  levels = unique(mixedsort(sample_data(phy_feces_clean)$Week))
) # here MAZ is calculated for each of the 3 timepoints as factor variables “Time” (in my setup, it was 1week,1 month and 1 year)

# Create MAZ dataframe
maz_all6 <- data.frame(
  SampleID = get_variable(phy_feces_clean, "SampleID") %>% as.character,
  AGE = get_variable(phy_feces_clean, "AGE"),
  Time_new = get_variable(phy_feces_clean, "Time_new"),
  Week = get_variable(phy_feces_clean, "Week"),
  child_ID = get_variable(phy_feces_clean, "child_ID"),
  oldchild01 = get_variable(phy_feces_clean, "oldchild01"),
  catanddog = get_variable(phy_feces_clean, "catanddog"),
  apartment_house = get_variable(phy_feces_clean, "apartment_house"),
  birthseason = get_variable(phy_feces_clean, "birthseason"),
  abbirth_mother = get_variable(phy_feces_clean, "abbirth_mother"),
  sex = get_variable(phy_feces_clean, "sex"),
  education = get_variable(phy_feces_clean, "education"),
  delivery01 = get_variable(phy_feces_clean, "delivery01"),
  Microbiota_age = pred_all_all6
) %>%
  group_by(Week) %>%
  mutate(
    maz_Week = (Microbiota_age - median(Microbiota_age)) / sd(Microbiota_age),
    maz_Week2 = ifelse(is.na(maz_Week), 0, maz_Week)
  )
#save(maz_all6, file = "maz_all6.RData")

# Create binary maturity classification
maz_all6_df <- maz_all6 %>% 
  group_by(Week) %>%
  mutate(
    Microbiota_age = as.numeric(Microbiota_age),
    Week_level2 = ntile(Microbiota_age, 2)
  ) %>%
  mutate(
    Gut_maturation_Week_two = case_when(
      Week_level2 == 1 ~ "Low",
      Week_level2 == 2 ~ "High"
    )
  )
rownames(maz_all6_df) <- maz_all6_df$SampleID
# Result: now set gut maturation to two level: Low, High
#save(maz_all6_df, file = "maz_all6_df.RData")


###### c. Chronological vs Microbiota age relationship ######
maz_all6#
maz_all6_df#

formula <- y ~ poly(x, 2, raw = TRUE)
Chronologic_Microbiota_day <- ggplot(maz_all6_df, aes(x = AGE, y = Microbiota_age)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_point(color = "azure4", alpha = 0.4) +
  stat_poly_line(formula = formula, color = "red") +
  stat_poly_eq(formula = formula, use_label(c("adj.R2", "p")), size = 3) +
  xlab("Chronologic age (days)") +
  ylab("Microbiota age (days)") +
  scale_x_continuous(breaks = seq(0, 2500, 500)) +
  guides(fill = guide_legend(title = "Time"), color = guide_legend(title = "Time")) +
  theme(text = element_text(size = 11))
Chronologic_Microbiota_day




#####  Figure 2C. Chronological age distributions in high- and low-maturity groups #####

###### a. Age influence on maturity classification ######
col_maturity <- c("High" = "#ff4040", "Low" = "#56b4e9")
age_influence_Week_1456y <- maz_all6_df %>%
  subset(!is.na(Gut_maturation_Week_two) & Time_new %in% c("1y", "4y", "5y", "6y")) %>%
  ggplot(aes(x = Time_new, y = AGE, color = Gut_maturation_Week_two)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_boxplot(aes(fill = Gut_maturation_Week_two), alpha = 0.7,
               outlier.colour = "#bcbbbb", outlier.alpha = 0.5) +
  xlab("Time") +
  ylab("Chronologic age (days)") +
  guides(fill = guide_legend(title = "Maturity"), color = guide_legend(title = "Maturity")) +
  scale_fill_manual(values = c("High" = "#ff4040", "Low" = "#56b4e9")) +
  scale_color_manual(values = c("High" = "#ff4040", "Low" = "#56b4e9")) +
  stat_compare_means(paired = FALSE, label.y = c(1000, 1900, 2300, 2600), size = 3) +
  theme(text = element_text(size = 11))
age_influence_Week_1456y


###### b. Microbiota age across timepoints ######
my_comparisons <- list(c("1w", "1m"), c("1m", "1y"), c("1y", "4y"), c("4y", "5y"), c("5y", "6y"))
col_time <- c("1w" = "#FF8C94", "1m" = "#49d261", "1y" = "#EE4266", "4y" = "#33B5BF",
              "5y" = "#9B59B6", "6y" = "#FFD159")
maz_all6_df$Time_new <- factor(maz_all6_df$Time_new, levels = c("1w", "1m", "1y", "4y", "5y", "6y"))

microbiota_age_difference_6time <- maz_all6_df %>% 
  ggplot(aes(x = Time_new, y = Microbiota_age)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_jitter(aes(fill = Time_new, color = Time_new), alpha = 0.4) +
  geom_boxplot(aes(fill = Time_new), alpha = 0.5,
               outlier.colour = "#bcbbbb", outlier.alpha = 0.5) +
  xlab("Time") +
  ylab("Microbiota age (days)") +
  scale_color_manual(values = col_time) +
  scale_fill_manual(values = col_time) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons,
                     paired = FALSE, label.y = c(500, 1100, 2000, 2200, 2400), size = 3) +
  theme(text = element_text(size = 11))
microbiota_age_difference_6time


###### c. Sankey diagram of the developmental trajectory of gut microbial maturation from one week to four years ######
maz_all6_df# 

# Select subjects with complete 4-timepoint data
maz_all6_first4timepoints_paired <- maz_all6_df %>%
  subset(Time_new %in% c("1w", "1m", "1y", "4y")) %>%
  group_by(child_ID) %>%
  filter(n() == 4) %>%
  ungroup()

# Convert to wide format
maz_all6_first4timepoints_paired_two <- maz_all6_first4timepoints_paired %>%
  pivot_wider(
    id_cols = c("child_ID"),
    names_from = "Time_new",
    values_from = "Gut_maturation_Week_two",
    values_fn = list
  ) %>%
  as.data.frame()

maz_all6_first4timepoints_paired_two <- do.call(cbind, lapply(maz_all6_first4timepoints_paired_two, unlist)) %>%
  data.frame()

# Calculate transition frequencies
maz_all6_first4timepoints_paired_two2 <- maz_all6_first4timepoints_paired_two %>% 
  group_by(X1w, X1m, X1y, X4y) %>% 
  dplyr::summarize(freq = n(), .groups = 'drop') %>% 
  ungroup() %>%
  data.frame()

# Create Sankey plot
maz_all6_first4timepoints_paired_two2$mode <- seq_along(maz_all6_first4timepoints_paired_two2[, 1]) %>% factor()
maz_all6_first4timepoints_paired_two2$X1w <- factor(maz_all6_first4timepoints_paired_two2$X1w,
                                                    levels = c("High", "Low"))

sankey_plot <- ggplot(data = maz_all6_first4timepoints_paired_two2,
                      aes(axis1 = X1w, axis2 = X1m, axis3 = X1y, axis4 = X4y, y = freq)) +
  geom_alluvium(aes(fill = mode)) +
  geom_stratum(fill = "grey93") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.title.y.left = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_continuous(breaks = 1:4, labels = c("1w", "1m", "1y", "4y")) + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
  scale_fill_viridis_d() +
  theme(text = element_text(size = 24)) +
  guides(fill = guide_legend(title = "Trajectory"), color = guide_legend(title = "Trajectory"))
sankey_plot


###### d. Number difference among gut maturation levels (mode) ######
maz_all6_first4timepoints_paired_two #

#### Transitions between Week 1 and Month 1
# Calculate frequency of transition patterns
maz_all6_first4timepoints_1w1m_paired_two2 <- maz_all6_first4timepoints_paired_two %>% 
  group_by(X1w,X1m) %>% 
  dplyr::summarize(freq = n(),.groups = 'drop')%>%
  ungroup()%>%
  data.frame()

# Define transition modes
maz_all6_first4timepoints_1w1m_paired_two2$mode=ifelse(
  maz_all6_first4timepoints_1w1m_paired_two2$X1w=="High" & 
    maz_all6_first4timepoints_1w1m_paired_two2$X1m=="High","mode1",
  ifelse(maz_all6_first4timepoints_1w1m_paired_two2$X1w=="High" & 
           maz_all6_first4timepoints_1w1m_paired_two2$X1m=="Low","mode2",
         ifelse(maz_all6_first4timepoints_1w1m_paired_two2$X1w=="Low" & 
                  maz_all6_first4timepoints_1w1m_paired_two2$X1m=="High","mode3","mode4")))

# # Prepare data for Fisher's exact test (pairwise_fisher_test)
maz_all6_first4timepoints_1w1m_paired_two2$freq_no<-sum(maz_all6_first4timepoints_1w1m_paired_two2$freq)-maz_all6_first4timepoints_1w1m_paired_two2$freq
df_1w1m_paired <- data.frame(maz_all6_first4timepoints_1w1m_paired_two2$mode,maz_all6_first4timepoints_1w1m_paired_two2$freq,maz_all6_first4timepoints_1w1m_paired_two2$freq_no)
df_1w1m_paired <- data.frame(t(df_1w1m_paired[,-1]))

# Set dimension names for contingency table
n <- ncol(df_1w1m_paired)
test_name <- character(n)
for (i in 1:n) {
  test_name[i] <- paste0("mode",i)
}
dimnames(df_1w1m_paired) <- list(maz_all6_first4timepoints_1w1m_paired_two2=c("Yes","No"),
                                 Groups=test_name)

# Perform pairwise Fisher's exact tests
fisher_1w1m <- pairwise_fisher_test(df_1w1m_paired,detailed = TRUE)
fisher_1w1m_0.05 <- fisher_1w1m[fisher_1w1m$p<0.05,]
print(fisher_1w1m_0.05)# No significant differences found



### Transitions between Month 1 and Year 1
# Calculate frequency of transition patterns
maz_all6_first4timepoints_1m1y_paired_two2 <- maz_all6_first4timepoints_paired_two %>% 
  group_by(X1m,X1y) %>% 
  dplyr::summarize(freq = n(),.groups = 'drop')%>%
  ungroup()%>%
  data.frame()

# Define transition modes
maz_all6_first4timepoints_1m1y_paired_two2$mode=ifelse(
  maz_all6_first4timepoints_1m1y_paired_two2$X1m=="High" & 
    maz_all6_first4timepoints_1m1y_paired_two2$X1y=="High","mode1",
  ifelse(maz_all6_first4timepoints_1m1y_paired_two2$X1m=="High" & 
           maz_all6_first4timepoints_1m1y_paired_two2$X1y=="Low","mode2",
         ifelse(maz_all6_first4timepoints_1m1y_paired_two2$X1m=="Low" & 
                  maz_all6_first4timepoints_1m1y_paired_two2$X1y=="High","mode3","mode4")))

# Prepare data for Fisher's exact test
maz_all6_first4timepoints_1m1y_paired_two2$freq_no<-sum(maz_all6_first4timepoints_1m1y_paired_two2$freq)-maz_all6_first4timepoints_1m1y_paired_two2$freq
df_1m1y_paired <- data.frame(maz_all6_first4timepoints_1m1y_paired_two2$mode,maz_all6_first4timepoints_1m1y_paired_two2$freq,maz_all6_first4timepoints_1m1y_paired_two2$freq_no)
df_1m1y_paired <- data.frame(t(df_1m1y_paired[,-1]))

# Set dimension names for contingency table
n <- ncol(df_1m1y_paired)
test_name <- character(n)
for (i in 1:n) {
  test_name[i] <- paste0("mode",i)
}
dimnames(df_1m1y_paired) <- list(maz_all6_first4timepoints_1m1y_paired_two2=c("Yes","No"),
                                 Groups=test_name)

# Perform pairwise Fisher's exact tests
fisher_1m1y <- pairwise_fisher_test(df_1m1y_paired,detailed = TRUE)
fisher_1m1y_0.05 <- fisher_1m1y[fisher_1m1y$p<0.05,]
print(fisher_1m1y_0.05) # No significant differences found


## Transitions between Year 1 and Year 4
# Calculate frequency of transition patterns
maz_all6_first4timepoints_1y4y_paired_two2 <- maz_all6_first4timepoints_paired_two %>% 
  group_by(X1y,X4y) %>% 
  dplyr::summarize(freq = n(),.groups = 'drop')%>%
  ungroup()%>%
  data.frame()

# Define transition modes
maz_all6_first4timepoints_1y4y_paired_two2$mode=ifelse(
  maz_all6_first4timepoints_1y4y_paired_two2$X1y=="High" & 
    maz_all6_first4timepoints_1y4y_paired_two2$X4y=="High","mode1",
  ifelse(maz_all6_first4timepoints_1y4y_paired_two2$X1y=="High" & 
           maz_all6_first4timepoints_1y4y_paired_two2$X4y=="Low","mode2",
         ifelse(maz_all6_first4timepoints_1y4y_paired_two2$X1y=="Low" & 
                  maz_all6_first4timepoints_1y4y_paired_two2$X4y=="High","mode3","mode4")))

# Prepare data for Fisher's exact test
maz_all6_first4timepoints_1y4y_paired_two2$freq_no<-sum(maz_all6_first4timepoints_1y4y_paired_two2$freq)-maz_all6_first4timepoints_1y4y_paired_two2$freq
df_1y4y_paired <- data.frame(maz_all6_first4timepoints_1y4y_paired_two2$mode,maz_all6_first4timepoints_1y4y_paired_two2$freq,maz_all6_first4timepoints_1y4y_paired_two2$freq_no)
df_1y4y_paired <- data.frame(t(df_1y4y_paired[,-1]))

# Set dimension names for contingency table
n <- ncol(df_1y4y_paired)
test_name <- character(n)
for (i in 1:n) {
  test_name[i] <- paste0("mode",i)
}
dimnames(df_1y4y_paired) <- list(maz_all6_first4timepoints_1y4y_paired_two2=c("Yes","No"),
                                 Groups=test_name)

# Perform pairwise Fisher's exact tests
fisher_1y4y <- pairwise_fisher_test(df_1y4y_paired,detailed = TRUE)
fisher_1y4y_0.05 <- fisher_1y4y[fisher_1y4y$p<0.05,]
print(fisher_1y4y_0.05)# mode3 -> mode4, p.adj = 0.0223

### Combine all statistical results
# Add time interval information to each result set
fisher_1w1m; fisher_1m1y; fisher_1y4y#

fisher_1w1m <- data.frame(Time="1w_1m",fisher_1w1m)
fisher_1m1y <- data.frame(Time="1m_1y",fisher_1m1y)
fisher_1y4y <- data.frame(Time="1y_4y",fisher_1y4y)

# Combine all results
fisher_first4timepoints_numdif <- rbind(fisher_1w1m,fisher_1m1y,fisher_1y4y)

# Export results to Excel
write_xlsx(fisher_first4timepoints_numdif,"fisher_first4timepoints_numdif.xlsx")




#####  Supplementary Figure S6. Impact of eight environmental factors on trajectory of gut microbial maturation #####
maz_all6_df #
phy_feces_clean #

# Add maturity data to phyloseq object
sample_data(phy_feces_clean) <- data.frame(sample_data(phy_feces_clean)) %>%
  mutate(
    Microbiota_age = ifelse(SampleID == maz_all6_df$SampleID, maz_all6_df$Microbiota_age, NA),
    maz_Week = ifelse(SampleID == maz_all6_df$SampleID, maz_all6_df$maz_Week, NA),
    Gut_maturation_Week_two = ifelse(SampleID == maz_all6_df$SampleID,
                                     maz_all6_df$Gut_maturation_Week_two, NA)
  )
sample_df <- data.frame(sample_data(phy_feces_clean))

# Define environmental factors to test
env_factors <- c("sex", "delivery01", "birthseason", "abbirth_mother", "education",
                 "oldchild01", "apartment_house", "catanddog")

# Function to perform Fisher tests and Cramer's V calculation
analyze_environmental_factors <- function(timepoint_data, timepoint_name) {
  # Fisher's exact tests
  fisher_results <- data.frame(Env_factor = character(0), p_value = numeric(0))
  cramer_results <- data.frame(Env_factor = character(0), cramerV = numeric(0))
  
  for (factor in env_factors) {
    # Fisher test
    contingency_table <- table(timepoint_data[[factor]], timepoint_data$Gut_maturation_Week_two)
    fisher_result <- fisher.test(contingency_table)
    fisher_results <- rbind(fisher_results, 
                            data.frame(Env_factor = factor, p_value = fisher_result$p.value))
    
    # Cramer's V
    v <- cramerV(contingency_table)
    cramer_results <- rbind(cramer_results,
                            data.frame(Env_factor = factor, cramerV = round(v, 3)))
  }
  
  # FDR adjustment
  fisher_results$p_adj <- p.adjust(fisher_results$p_value, method = "BH")
  
  return(list(fisher = fisher_results, cramer = cramer_results))
}

# Analyze each timepoint
timepoints <- c("1w", "1m", "1y", "4y", "5y", "6y")
fisher_list <- list()
cramer_list <- list()

for (tp in timepoints) {
  tp_data <- subset(sample_df, Time_new == tp)
  results <- analyze_environmental_factors(tp_data, tp)
  fisher_list[[tp]] <- results$fisher
  cramer_list[[tp]] <- results$cramer
}

# Combine results
fisher_env_time <- bind_cols(
  fisher_list[["1w"]] %>% select(p_1w = "p_value", p_adj_1w = "p_adj"),
  fisher_list[["1m"]] %>% select(p_1m = "p_value", p_adj_1m = "p_adj"),
  fisher_list[["1y"]] %>% select(p_1y = "p_value", p_adj_1y = "p_adj"),
  fisher_list[["4y"]] %>% select(p_4y = "p_value", p_adj_4y = "p_adj"),
  fisher_list[["5y"]] %>% select(p_5y = "p_value", p_adj_5y = "p_adj"),
  fisher_list[["6y"]] %>% select(p_6y = "p_value", p_adj_6y = "p_adj")
) %>% data.frame() %>% round(digits = 4)
fisher_env_time$Env_factor <- env_factors
fisher_env_time

cramer_env_time <- bind_cols(
  cramer_list[["1w"]] %>% select(cramerV_1w = "cramerV"),
  cramer_list[["1m"]] %>% select(cramerV_1m = "cramerV"),
  cramer_list[["1y"]] %>% select(cramerV_1y = "cramerV"),
  cramer_list[["4y"]] %>% select(cramerV_4y = "cramerV"),
  cramer_list[["5y"]] %>% select(cramerV_5y = "cramerV"),
  cramer_list[["6y"]] %>% select(cramerV_6y = "cramerV")
) %>% data.frame() %>% round(digits = 3)
cramer_env_time$Env_factor <- env_factors
cramer_env_time

# Export results
write_xlsx(fisher_env_time, "fisher_env_time.xlsx")
write_xlsx(cramer_env_time, "cramer_env_time.xlsx")

# Create combined visualization
fisher_cramer_env1 <- reshape2::melt(fisher_env_time) %>% 
  subset(variable %in% c("p_adj_1w", "p_adj_1m", "p_adj_1y", "p_adj_4y", "p_adj_5y", "p_adj_6y"))
fisher_cramer_env1$variable <- ifelse(fisher_cramer_env1$variable == "p_adj_1w", "1w",
                                      ifelse(fisher_cramer_env1$variable == "p_adj_1m", "1m",
                                             ifelse(fisher_cramer_env1$variable == "p_adj_1y", "1y",
                                                    ifelse(fisher_cramer_env1$variable == "p_adj_4y", "4y",
                                                           ifelse(fisher_cramer_env1$variable == "p_adj_5y", "5y",
                                                                  ifelse(fisher_cramer_env1$variable == "p_adj_6y", "6y", NA))))))

fisher_cramer_env2 <- reshape2::melt(cramer_env_time)
fisher_cramer_env2$variable <- ifelse(fisher_cramer_env2$variable == "cramerV_1w", "1w",
                                      ifelse(fisher_cramer_env2$variable == "cramerV_1m", "1m",
                                             ifelse(fisher_cramer_env2$variable == "cramerV_1y", "1y",
                                                    ifelse(fisher_cramer_env2$variable == "cramerV_4y", "4y",
                                                           ifelse(fisher_cramer_env2$variable == "cramerV_5y", "5y",
                                                                  ifelse(fisher_cramer_env2$variable == "cramerV_6y", "6y", NA))))))

# Set factor levels
factor_levels <- c("oldchild01", "catanddog", "apartment_house", "abbirth_mother",
                   "sex", "delivery01", "birthseason", "education")
time_levels <- c("1w", "1m", "1y", "4y", "5y", "6y")

fisher_cramer_env1$Env_factor <- factor(fisher_cramer_env1$Env_factor, levels = factor_levels)
fisher_cramer_env1$variable <- factor(fisher_cramer_env1$variable, levels = time_levels)
fisher_cramer_env2$Env_factor <- factor(fisher_cramer_env2$Env_factor, levels = factor_levels)
fisher_cramer_env2$variable <- factor(fisher_cramer_env2$variable, levels = time_levels)

# Color scheme
fisher_cramer_col <- c("oldchild01" = "#E69F00", "abbirth_mother" = "#56B4E9", 
                       "catanddog" = "purple", "sex" = "#69b3a2",
                       "apartment_house" = "#00798c", "delivery01" = "#ffeda0", 
                       "birthseason" = "#fdbf6f", "education" = "#bf812d")

# Create dual-axis plot
max_logp <- max(-log10(fisher_cramer_env1$value), na.rm = TRUE)
max_v <- max(fisher_cramer_env2$value, na.rm = TRUE)

p_fisher_cramer <- ggplot() +
  geom_bar(data = fisher_cramer_env1, 
           aes(x = variable, y = -log10(value), fill = Env_factor),
           stat = "identity", position = "dodge", alpha = 0.8) +
  geom_bar(data = fisher_cramer_env2, 
           aes(x = variable, y = -value * max_logp / max_v, fill = Env_factor),
           stat = "identity", position = "dodge", alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#e31a1c") +
  geom_hline(yintercept = 0, color = "black") +
  scale_y_continuous(
    name = "-log10(p value)",
    breaks = seq(0, max_logp, by = 0.8),
    sec.axis = sec_axis(~ -. * max_v / max_logp, name = "Cramer's V",
                        breaks = seq(0, max_v, by = 0.04))
  ) +
  scale_color_manual(values = fisher_cramer_col) +
  scale_fill_manual(values = fisher_cramer_col) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1),
        legend.position = "bottom",
        axis.text.y.left = element_text(color = "#e31a1c"),
        axis.ticks.y.left = element_line(color = "#e31a1c"),
        axis.title.y.left = element_text(color = "#e31a1c", hjust = 0.85),
        axis.line.y.left = element_line(color = "#e31a1c"),
        axis.text.y.right = element_text(color = "blue"),
        axis.ticks.y.right = element_line(color = "blue"),
        axis.title.y.right = element_text(color = "blue", hjust = 0.8),
        axis.line.y.right = element_line(color = "blue")
  )
p_fisher_cramer




##### Figure 2D. Differences in microbiota-by-age z-score (MAZ) at one year old according to environmental exposures #####
phy_feces_clean#
maz_all6_df#

### a. oldchild01_1y
sample_df <- maz_all6_df

sample_df$oldchild01 <- ifelse(sample_df$oldchild01=="0","No",
                               ifelse(sample_df$oldchild01=="1","Yes",sample_df$oldchild01)) %>% factor(levels = c("No","Yes"))

# Microbiota_age difference on oldchild01 group_all data at 1y
maz_Week_all6_1y_oldchild01 <- sample_df%>%
  subset(!is.na(oldchild01)&Time_new %in%c("1y"))%>%
  group_by(oldchild01)%>%
  ggplot(aes(x=oldchild01,y=maz_Week2,fill=oldchild01))+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_jitter(aes(color=oldchild01),alpha=0.2,size=3)+
  geom_boxplot(alpha=0.8,outlier.colour="#f0f0f0",outlier.alpha = 0.2)+
  xlab("Siblings")+
  ylab("Microbiota by age z-scores (MAZ)")+
  theme(legend.position = "none")+
  stat_compare_means(aes(label=paste0("italic(P)==",after_stat(p.format))),
                     parse=TRUE,size = 3)+# defalt: method="wilcox.test",paired=FALSE
  scale_color_manual(values=c("grey","#E69F00"))+
  scale_fill_manual(values=c("grey","#E69F00"))+
  theme(text = element_text(size=11))
maz_Week_all6_1y_oldchild01

### b. abbirth_mother_1y
sample_df$abbirth_mother <- ifelse(sample_df$abbirth_mother=="0","No",
                                   ifelse(sample_df$abbirth_mother=="1","Yes",sample_df$abbirth_mother)) %>%factor(levels = c("No","Yes"))

maz_Week_all6_1y_abbirth_mother <- sample_df%>%
  subset(!is.na(abbirth_mother)&Time_new =="1y")%>%
  group_by(abbirth_mother)%>%
  ggplot(aes(x=abbirth_mother,y=maz_Week2,fill=abbirth_mother))+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_jitter(aes(color=abbirth_mother),alpha=0.2,size=3)+
  geom_boxplot(alpha=0.8,outlier.colour="#f0f0f0",outlier.alpha = 0.2)+
  xlab("Antibiotics (mother)")+
  ylab("Microbiota by age z-scores (MAZ)")+
  theme(legend.position = "none")+
  stat_compare_means(aes(label=paste0("italic(P)==",after_stat(p.format))),
                     parse=TRUE,size = 3)+
  scale_color_manual(values=c("grey","#56B4E9"))+
  scale_fill_manual(values=c("grey","#56B4E9"))+
  theme(text = element_text(size=11))
maz_Week_all6_1y_abbirth_mother

### c. cat and dog_1y
sample_df$catanddog <- ifelse(sample_df$catanddog=="0","No",
                              ifelse(sample_df$catanddog=="1","Yes",sample_df$catanddog)) %>%factor(levels = c("No","Yes"))

maz_Week_all6_1y_catanddog <- sample_df%>%
  subset(!is.na(catanddog)&Time_new%in%c("1y"))%>%
  group_by(catanddog)%>%
  ggplot(aes(x=catanddog,y=maz_Week2,fill=catanddog))+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_jitter(aes(color=catanddog),alpha=0.2,size=3)+
  geom_boxplot(alpha=0.8,outlier.colour="#f0f0f0",outlier.alpha = 0.2)+
  xlab("Cat and dog")+
  ylab("Microbiota by age z-scores (MAZ)")+
  theme(legend.position = "none")+
  stat_compare_means(aes(label=paste0("italic(P)==",after_stat(p.format))),
                     parse=TRUE,size = 3)+
  theme(text = element_text(size=11))+
  scale_color_manual(values=c("grey", "purple"))+
  scale_fill_manual(values=c("grey", "purple"))
maz_Week_all6_1y_catanddog

### d. apartment_house_1y
maz_Week_all6_1y_apartment_house <- sample_df%>%
  subset(!is.na(apartment_house)&Time_new=="1y")%>%
  group_by(apartment_house)%>%
  ggplot(aes(x=apartment_house,y=maz_Week2,fill=apartment_house))+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_jitter(aes(color=apartment_house),alpha=0.2,size=3)+
  geom_boxplot(alpha=0.8,outlier.colour="#f0f0f0",outlier.alpha = 0.2)+
  xlab("House type")+
  ylab("Microbiota by age z-scores (MAZ)")+
  theme(legend.position = "none")+
  stat_compare_means(aes(label=paste0("italic(P)==",after_stat(p.format))),
                     parse=TRUE,size = 3)+
  theme(text = element_text(size=11))+
  scale_color_manual(values=c("grey","#00798c"))+
  scale_fill_manual(values=c("grey","#00798c"))
maz_Week_all6_1y_apartment_house

### e. sex_1y
maz_Week_all6_1y_sex <- sample_df%>%
  subset(!is.na(sex)&Time_new%in%c("1y"))%>%
  group_by(sex)%>%
  ggplot(aes(x=sex,y=maz_Week2,fill=sex))+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_jitter(aes(color=sex),alpha=0.2,size=3)+
  geom_boxplot(alpha=0.8,outlier.colour="#f0f0f0",outlier.alpha = 0.2)+
  xlab("Sex")+
  ylab("Microbiota by age z-scores (MAZ)")+
  theme(legend.position = "none")+
  stat_compare_means(aes(label=paste0("italic(P)==",after_stat(p.format))),
                     parse=TRUE,size = 3)+
  theme(text = element_text(size=11))+
  scale_color_manual(values=c("grey","#69b3a2"))+
  scale_fill_manual(values=c("grey","#69b3a2"))
maz_Week_all6_1y_sex

## arrange plot together
maz_Week_all6_1y_envfactor <- ggarrange(maz_Week_all6_1y_oldchild01,
                                        maz_Week_all6_1y_catanddog,
                                        maz_Week_all6_1y_apartment_house,
                                        maz_Week_all6_1y_abbirth_mother,
                                        # maz_Week_all6_1y_sex,
                                        ncol=4, nrow=1)
maz_Week_all6_1y_envfactor




##### Figure 2E. Environmental influences on the rate of microbiota maturation #####
maz_all6_df #

### Effect on the slope of gut microbiota age development

# Subset: individuals with paired samples at 1w, 1m, 1y, and 4y
maz_all6_first4timepoints_paired <- maz_all6_df%>%
  subset(Time_new=="1w"|Time_new=="1m"|Time_new=="1y"|Time_new=="4y")%>%
  group_by(child_ID)%>%
  filter(n()==4)%>%
  ungroup() 
# 278 people at 4 time points

# Convert from long to wide format
maz_all6_first4timepoints_paired_two <- maz_all6_first4timepoints_paired%>%
  pivot_wider( id_cols = c("child_ID"),
               names_from = "Time_new",
               values_from="Gut_maturation_Week_two",
               values_fn = list)%>%
  as.data.frame()

# Convert list-columns to simple data.frame
maz_all6_first4timepoints_paired_two <- do.call(cbind,lapply(maz_all6_first4timepoints_paired_two,unlist))%>%
  data.frame()

# Define longitudinal trajectory patterns
test1 <- maz_all6_first4timepoints_paired_two

# Bring in sample metadata
test2 <- data.frame(sample_data(phy_feces_clean))
test2$child_ID <- as.character(test2$child_ID)

# Select necessary columns
test3 <- test2 %>% dplyr::select(child_ID,AGE, Microbiota_age)

# Compute individual slopes of microbiota-age increase across time
slopes <- test3 %>%
  group_by(child_ID) %>%
  filter(n() == 4) %>%
  do(slope_data = coef(lm(Microbiota_age ~ AGE, data = .))[2]) %>%
  unnest(cols = c(slope_data)) %>%
  rename(slope = slope_data)

# Merge slopes with environmental factors
test4 <- left_join(slopes, test2, by = "child_ID")

# Linear regression model: environmental predictors of microbiota-age slope
lm_model <- lm(slope ~ oldchild01 + catanddog + apartment_house+abbirth_mother+sex+delivery01+education+birthseason, 
               data = test4)

summary(lm_model)
anova(lm_model)
# Result interpretation:
# Each coefficient (estimate) indicates the linear influence of a given environmental factor 
# on the rate of gut microbiota maturation:
#   - Positive estimate: accelerates maturation
#   - Negative estimate: slows down maturation
#   - p-value < 0.05: statistically significant effect.

# Forest-style visualization of model effects
# Extract coefficients and confidence intervals
lm_summary <- tidy(lm_model) %>%
  filter(term != "(Intercept)") %>%
  select(term2 = "term", estimate, std.error, p_summary="p.value")

lm_summary$term <- ifelse(lm_summary$term2=="sexMale","sex",
                          ifelse(lm_summary$term2=="apartment_houseHouse","apartment_house",
                                 ifelse(lm_summary$term2%in%c("education2Medium","education3High"),"education",
                                        ifelse(lm_summary$term2%in%c("birthseasonspring","birthseasonsummer","birthseasonwinter"),
                                               "birthseason",lm_summary$term2))))

# Extract ANOVA results (sum of squares and p-values)
anova_summary <- anova(lm_model) %>%
  tibble::rownames_to_column("term")%>% 
  dplyr::select(term, sumsq = `Sum Sq`, p_anova = `Pr(>F)`)%>% 
  as.data.frame()%>% 
  subset(term != "Residuals")

anova_summary$col <- ifelse(anova_summary$p_anova<0.05, "Significant","No_significant")

# Merge coefficient and ANOVA summaries
merged_df <- left_join(lm_summary, anova_summary, by = "term")

merged_df$term <- factor(merged_df$term, 
                         levels = c("education","birthseason","delivery01","sex",
                                    "abbirth_mother","apartment_house","catanddog","oldchild01"))

# Plot 1: coefficient estimates and 95% CI
p_summary_lm_slope <-ggplot(merged_df, aes(x =reorder(term,sumsq), y = estimate)) +
  geom_point(aes(#size = sumsq, 
    color = p_summary)) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +
  coord_flip() +
  scale_color_gradient(low = "red", high = "black", name = "p-value",
                       breaks=c(0.05,0.25,0.50,0.75)
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Effect Size (Estimate)",
    title = "Effect on the slope of gut microbiota age development"
  )
p_summary_lm_slope

# Plot 2: explained variance (Sum of Squares) per factor
p_anova_lm_slope <-ggplot(merged_df, aes(x =reorder(term,sumsq), y = 0)) +
  geom_point(aes(size = sumsq, 
                 color = p_anova)) +
  coord_flip() +
  scale_color_gradient(low = "red", high = "black", name = "p-value",
                       breaks=c(0.05,0.25,0.50,0.75)
  ) +
  scale_size(range = c(2, 8), name = "Sum Sq") +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Explained Variance (Sum of Squares)",
    title = "Effect on the slope of gut microbiota age development"
  )
p_anova_lm_slope




##### Supplementary Figure S8. GEE model examining associations between early-life environmental exposures and gut microbiota maturity across four early time points #####
maz_all6_df #

###### Process paired samples from first 4 timepoints
maz_all6_first4timepoints_paired <- maz_all6_df %>%
  subset(Time_new %in% c("1w", "1m", "1y", "4y")) %>%
  group_by(child_ID) %>%
  filter(n() == 4) %>%  # Keep only subjects with all 4 timepoints
  ungroup()  # 278 subjects with complete 4-timepoint data

# Convert from long to wide format
maz_all6_first4timepoints_paired_wide <- maz_all6_first4timepoints_paired %>%
  pivot_wider(id_cols = c("child_ID"),
              names_from = "Time_new",
              values_from = "Gut_maturation_Week_two",
              values_fn = list) %>%
  as.data.frame()

# Convert nested lists to data frame
maz_all6_first4timepoints_paired_wide <- do.call(cbind, 
                                                 lapply(maz_all6_first4timepoints_paired_wide, unlist)) %>% 
  data.frame()

# Prepare data for GEE analysis
analysis_data <- data.frame(sample_data(phy_feces_clean))
analysis_data$child_ID <- as.character(analysis_data$child_ID)

# Subset to subjects with complete 4-timepoint data
analysis_data <- subset(analysis_data, 
                        Time_new %in% c("1w", "1m", "1y", "4y") & 
                          child_ID %in% maz_all6_first4timepoints_paired_wide$child_ID)

analysis_data$child_ID <- as.numeric(as.character(analysis_data$child_ID))

# Convert maturity to binary variable
analysis_data$Gut_maturation_Week_two_bin <- ifelse(analysis_data$Gut_maturation_Week_two == "High", 1, 
                                                    ifelse(analysis_data$Gut_maturation_Week_two == "Low", 0, NA))

# Select variables and remove missing values
analysis_data_clean <- analysis_data %>%
  dplyr::select(Gut_maturation_Week_two_bin, oldchild01, catanddog, apartment_house,
                abbirth_mother, sex, delivery01, education, birthseason, child_ID, Time_new) %>%
  na.omit()

# GEE Model - accounting for longitudinal correlation
analysis_data_clean$Time_new <- factor(analysis_data_clean$Time_new, 
                                       levels = c("1w", "1m", "1y", "4y"))

gee_model <- geeglm(Gut_maturation_Week_two_bin ~ Time_new +
                      oldchild01 + catanddog + apartment_house + abbirth_mother + 
                      sex + delivery01 + education + birthseason + 
                      oldchild01:Time_new + catanddog:Time_new + apartment_house:Time_new + 
                      abbirth_mother:Time_new + sex:Time_new + delivery01:Time_new + 
                      education:Time_new + birthseason:Time_new,
                    id = child_ID,
                    family = binomial,
                    corstr = "exchangeable",
                    data = analysis_data_clean)
summary(gee_model)

# Extract OR and 95% CI
OR <- exp(coef(gee_model))
CI <- exp(confint.default(gee_model))

results_df <- data.frame(
  Variable = names(OR),
  OR = round(OR, 3),
  CI_lower = round(CI[, 1], 3),
  CI_upper = round(CI[, 2], 3),
  p_value = round(summary(gee_model)$coefficients[, 4], 4))
print(results_df)

# Create forest plot using ggplot2
results_df$Variable <- factor(results_df$Variable, levels = rev(results_df$Variable))

p1_GEE_model <- ggplot(results_df, aes(x = OR, y = Variable)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_log10() +
  theme_minimal(base_size = 14) +
  xlab("Odds Ratio (log scale)") +
  ylab(NULL) +
  theme(panel_grid_minor = element_blank()) +
  labs(title = "GEE Model: Gut Maturity vs. Environmental Factors")
p1_GEE_model

# Create medical-style forest plot
tabletext <- cbind(
  Variable = results_df$Variable,
  OR = sprintf("%.2f (%.2f - %.2f)", results_df$OR, results_df$CI_lower, results_df$CI_upper),
  p = sprintf("%.3f", results_df$p_value)
)

plot_data <- cbind(results_df$OR, results_df$CI_lower, results_df$CI_upper)

p2_GEE_model <- forestplot(labeltext = tabletext,
                           mean = plot_data[, 1],
                           lower = plot_data[, 2],
                           upper = plot_data[, 3],
                           zero = 1,
                           xlog = TRUE,
                           col = fpColors(box = "black", lines = "gray30", zero = "gray50"),
                           boxsize = 0.3,
                           line.margin = 0.1,
                           txt_gp = fpTxtGp(label = gpar(fontsize = 12), 
                                            ticks = gpar(fontsize = 10)),
                           title = "Odds Ratios from GEE Model",
                           xlab = "Odds Ratio (log scale)")
p2_GEE_model




##### Figure 2 F. Longitudinal MAZ trajectories stratified by birth season, based on children with samples at all four time points #####

# Prepare birth season analysis data
birthseason_data <- maz_all6_df %>%
  subset(Time_new %in% c("1w", "1m", "1y", "4y")) %>%
  group_by(child_ID) %>%
  filter(n() == 4) %>%
  ungroup() %>%
  group_by(Time_new, birthseason) %>%
  summarise(
    median_MAZ = median(maz_Week2, na.rm = TRUE),
    se = sd(maz_Week2, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Define colors
season_colors <- c("spring" = "#fc8d59", "summer" = "#2171b5", 
                   "autumn" = "#9ecae1", "winter" = "#d73027")

# Create trajectory plot
birthseason_data$birthseason <- factor(birthseason_data$birthseason, levels = c("spring","summer","autumn","winter"))
p_birthseason_trajectory <- ggplot(birthseason_data, 
                                   aes(x = Time_new, y = median_MAZ, 
                                       color = birthseason, group = birthseason)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = median_MAZ - se, ymax = median_MAZ + se), width = 0.2) +
  scale_color_manual(values = season_colors) +
  labs(y = "MAZ score", x = "Age", color = "Birth season") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8))
p_birthseason_trajectory

# Statistical comparisons between birth seasons
birthseason_comparisons <- list(
  c("spring", "summer"), c("spring", "autumn"), c("spring", "winter"),
  c("summer", "autumn"), c("summer", "winter"), c("autumn", "winter")
)

# Perform comparisons for each timepoint
### 1w comparisons
Wil_test_adj_season_MAZ_1w <- birthseason_data%>% 
  subset(Time_new%in%c("1w"))%>%
  compare_means(maz_Week2 ~ birthseason, 
                comparisons = birthseason_com, 
                p.adjust.method = "fdr")
Wil_test_adj_season_MAZ_1w$Time <- "1w"

### 1m comparisons
Wil_test_adj_season_MAZ_1m <- birthseason_data%>% 
  subset(Time_new%in%c("1m"))%>%
  compare_means(maz_Week2 ~ birthseason, 
                comparisons = birthseason_com, 
                p.adjust.method = "fdr")
Wil_test_adj_season_MAZ_1m$Time <- "1m"

### 1y comparisons
Wil_test_adj_season_MAZ_1y <- birthseason_data%>% 
  subset(Time_new%in%c("1y"))%>%
  compare_means(maz_Week2 ~ birthseason, 
                comparisons = birthseason_com, 
                p.adjust.method = "fdr")
Wil_test_adj_season_MAZ_1y$Time <- "1y"

### 4y comparisons
Wil_test_adj_season_MAZ_4y <- birthseason_data%>% 
  subset(Time_new%in%c("4y"))%>%
  compare_means(maz_Week2 ~ birthseason, 
                comparisons = birthseason_com, 
                p.adjust.method = "fdr")
Wil_test_adj_season_MAZ_4y$Time <- "4y"

# Combine all statistical test results
Wil_test_adj_season_MAZ <- rbind(Wil_test_adj_season_MAZ_1w,Wil_test_adj_season_MAZ_1m,
                                 Wil_test_adj_season_MAZ_1y,Wil_test_adj_season_MAZ_4y)
# Format p-values

Wil_test_adj_season_MAZ <- Wil_test_adj_season_MAZ %>% 
  mutate(p.adj = round(p.adj, 3),
         p.format = formatC(p.adj, format = "f", digits = 3))

# Export results
write_xlsx(Wil_test_adj_season_MAZ,"Wil_test_adj_season_MAZ.xlsx")

