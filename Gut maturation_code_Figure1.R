########################### Gut microbial maturation_Result 1 ######################

#### Figure 1. The developmental trajectory of microbial composition in the gut ####

# Load required packages
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(rstatix)
library(gratia)
library(mgcv)
library(writexl)

# Load phyloseq object
load("phy_feces_clean.RData")
phy_feces_clean#


colnames(sample_data(phy_feces_clean))
### Prepare data
alpha_phy <-estimate_richness(phy_feces_clean, 
                              measures = c("Observed","Shannon","Chao1","InvSimpson"))
alpha_phy_plot <- merge(alpha_phy, sample_data(phy_feces_clean),by="row.names")
alpha_phy_plot$Time_new <- factor(as.character(alpha_phy_plot$Time_new),
                                  levels = c("1w","1m","1y","4y","5y","6y"))

# Define time comparisons and color palette
my_comparisons <- list(c("1w","1m"),c("1m","1y"),c("1y","4y"),c("4y","5y"),c("5y","6y"))
col_time <- c("1w"="#FF8C94","1m"="#49d261","1y"="#EE4266","4y"="#33B5BF","5y"="#9B59B6","6y"="#FFD159")





#####  Figure 1A. Alpha diversity across all time points #####


### Faith's PD
library(picante)
alpha_fpd <- pd(t(otu_table(phy_feces_clean)), phy_tree(phy_feces_clean))
alpha_fpd <- merge(data.frame(sample_data(phy_feces_clean)),data.frame(alpha_fpd),by="row.names")
alpha_fpd$Time_new <- factor(alpha_fpd$Time_new,levels = c("1w","1m","1y","4y","5y","6y"))
col_time <- c("1w"="#FF8C94","1m"="#49d261","1y"="#EE4266","4y"="#33B5BF","5y"="#9B59B6",
              "6y"="#FFD159")
my_comparisons <- list(c("1w","1m"),c("1m","1y"),c("1y","4y"),c("4y","5y"),c("5y","6y"))

Time_alpha_fpd <- ggplot(data=alpha_fpd,aes(x=Time_new,y=PD)) +
  theme_bw()+theme(panel.grid = element_blank())+
  geom_point(aes(fill=Time_new,color=Time_new),alpha=0.1)+
  geom_line(aes(group = child_ID,color=Time_new), alpha = 0.1)+
  geom_boxplot(aes(fill=Time_new),alpha=0.5)+
  
  xlab("Time")+
  ylab("Phylogenetic diversity")+
  scale_color_manual(values=col_time)+
  scale_fill_manual(values=col_time)+
  theme(legend.position ="none")+
  guides(fill=guide_legend(title="Time"),color=guide_legend(title="Time"))+
  stat_compare_means(comparisons = my_comparisons,label.y = c(23,25,28,30,32), paried=FALSE)+ # wilcox.test
  theme(text = element_text(size=11))
Time_alpha_fpd


model_Time_fpd <- aov(PD ~ Time_new,data=alpha_fpd); summary(model_Time_fpd)
TukeyHSD(model_Time_fpd)

# test
library(ggpubr)
Wil_test_adj_Time_PD <- compare_means(PD ~ Time_new, 
                                      comparisons = my_comparisons, 
                                      p.adjust.method = "fdr",
                                      data=alpha_fpd)
Wil_test_adj_Time_PD
View(Wil_test_adj_Time_PD)

# LMM for Shannon as sensitivity
library(emmeans)
colnames(alpha_fpd)
m_fpd <- lmer(PD ~ Time_new +
                    oldchild01+catanddog+apartment_house+abbirth_mother+sex+delivery01+education+birthseason+
                    (1 | child_ID),
                  data = alpha_fpd,
                  REML = FALSE)

emm_fpd <- emmeans(m_fpd, ~ Time_new)

fpd_pairwise <- contrast(emm_fpd,method = "consec", adjust = "BH") %>%
  as.data.frame()
fpd_pairwise






### Shannon plot
Time_Shannon <- ggplot(alpha_phy_plot,aes(x=Time_new,y=Shannon))+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_point(aes(fill=Time_new,color=Time_new),alpha=0.1)+
  geom_line(aes(group = child_ID,color=Time_new), alpha = 0.1)+
  geom_boxplot(aes(fill=Time_new),alpha=0.5)+
  xlab("Time")+
  ylab("Shannon")+
  scale_color_manual(values=col_time)+
  scale_fill_manual(values=col_time)+
  theme(legend.position ="none")+
  guides(fill=guide_legend(title="Time"),color=guide_legend(title="Time"))+
  stat_compare_means(comparisons = my_comparisons,label.y = c(3,3.5,4, 4.5, 5), paried=FALSE)+
  theme(text = element_text(size=11))
print(Time_Shannon)

# Wilcoxon test with FDR correction
Wil_test_adj_Time_Shannon <- compare_means(Shannon ~ Time_new, 
                                           comparisons = my_comparisons, 
                                           p.adjust.method = "fdr",
                                           data=alpha_phy_plot)
print(Wil_test_adj_Time_Shannon)
View(Wil_test_adj_Time_Shannon)
write_xlsx(Wil_test_adj_Time_Shannon, path = "Wil_test_adj_Time_Shannon.xlsx")

table(alpha_phy_plot$Time_new)

# LMM for Shannon as sensitivity
m_Shannon <- lmer(Shannon ~ Time_new +
    oldchild01+catanddog+apartment_house+abbirth_mother+sex+delivery01+education+birthseason+
    (1 | child_ID),
  data = alpha_phy_plot,
  REML = FALSE)

emm_Shannon <- emmeans(m_Shannon, ~ Time_new)

Shannon_pairwise <- contrast(emm_Shannon,method = "consec", adjust = "BH") %>%
  as.data.frame()
Shannon_pairwise


### Shannon increase speed 
alpha_phy_plot#

# Subset by time points
df_1y <- subset(alpha_phy_plot, Time_new%in%c("1y"))
df_4y <- subset(alpha_phy_plot, Time_new%in%c("4y"))
df_5y <- subset(alpha_phy_plot, Time_new%in%c("5y"))

# Mean Shannon at each time (NA-safe)
mean_1y <- mean(df_1y$Shannon, na.rm = TRUE)
mean_4y <- mean(df_4y$Shannon, na.rm = TRUE)
mean_5y <- mean(df_5y$Shannon, na.rm = TRUE)

# Calculate the ratio of consecutive annual gains
delta_1yto4y <- (mean_4y - mean_1y)/3
delta_4yto5y <- mean_5y - mean_4y

# Per-year growth from 1y to 4y vs. the single-year growth 4y to 5y
if (delta_4yto5y != 0) {
  ratio <- delta_1yto4y / delta_4yto5y
} else {
  ratio <- NA
  warning("Denominator is zero; ratio set to NA.")
}

# Print 
print(paste("Ratio =", round(ratio, 5)))
# "Relative ratio =  2.57589"



### Observed plot
Time_Observed <- ggplot(alpha_phy_plot,aes(x=Time_new,y=Observed))+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_point(aes(fill=Time_new,color=Time_new),alpha=0.1)+
  geom_line(aes(group = child_ID,color=Time_new), alpha = 0.1)+
  geom_boxplot(aes(fill=Time_new),alpha=0.5)+
  xlab("Time")+
  ylab("Observed")+
  scale_color_manual(values=col_time)+
  scale_fill_manual(values=col_time)+
  theme(legend.position ="none")+
  guides(fill=guide_legend(title="Time"),color=guide_legend(title="Time"))+
  stat_compare_means(comparisons = my_comparisons,label.y = c(200,220, 240, 260,280), paried=FALSE)+
  theme(text = element_text(size=11))
print(Time_Observed)

# Wilcoxon test with FDR correction
Wil_test_adj_Time_Observed <- compare_means(Observed ~ Time_new, 
                                            comparisons = my_comparisons, 
                                            p.adjust.method = "fdr",
                                            data=alpha_phy_plot)
print(Wil_test_adj_Time_Observed)
write_xlsx(Wil_test_adj_Time_Observed, path = "Wil_test_adj_Time_Observed.xlsx")

table(alpha_phy_plot$Time_new)

# LMM for Observed as sensitivity
m_Observed <- lmer(Observed ~ Time_new +
                     oldchild01+catanddog+apartment_house+abbirth_mother+sex+delivery01+education+birthseason+
                     (1 | child_ID),
                   data = alpha_phy_plot,
                   REML = FALSE)

emm_Observed <- emmeans(m_Observed, ~ Time_new)

Observed_pairwise <- contrast(emm_Observed,method = "consec", adjust = "BH") %>%
  as.data.frame()
Observed_pairwise










##### Figure 1B. Principal coordinates analysis (PCoA) for Bray-Curtis dissimilarity by sampling time ####
## Bray Curtis distance matrix
phy_feces_clean#

# Calculate bray curtis distance
bray_phy <- distance(phy_feces_clean,method = "bray",type = "samples")

# # PERMANOVA (adonis test)
set.seed(33);
adonis_phy_bray <- adonis2(bray_phy~Time_new+
                             oldchild01+catanddog+apartment_house+abbirth_mother+sex+delivery01+education+birthseason,
                           data = data.frame(sample_data(phy_feces_clean)),
                           na.action = na.omit,permutations = 999)
adonis_phy_bray
# P=0.001. R2=0.11361
# P=0.001. R2=0.12547


# PCoA ordination
ord_bray_phy <- ordinate(phy_feces_clean,"PCoA","bray")

# Transform multidimensional scaling (MDS) to a data matrix
pcoa_bray_phy <- cmdscale(bray_phy,eig = TRUE) #Classical multidimensional scaling (MDS) of a data matrix. 
# Collect in data.frame
pcoa_bray_phy_df <- data.frame(PC1_bray=c(pcoa_bray_phy$points[,1]),
                               PC2_bray=c(pcoa_bray_phy$points[,2]),
                               Sample=rownames(pcoa_bray_phy$points))
sample_data(phy_feces_clean)$PC1_bray <- pcoa_bray_phy_df$PC1_bray
sample_data(phy_feces_clean)$PC2_bray <- pcoa_bray_phy_df$PC2_bray


# Bray_curtis PCoA scatter plot
sample_data(phy_feces_clean)$Time_new <- factor(sample_data(phy_feces_clean)$Time_new,
                                                levels = c("1w","1m","1y","4y","5y","6y"))
plot_Time_bray_curtis <- plot_ordination(phy_feces_clean,ord_bray_phy,type = "sample",color = "Time_new")+ # plot_bray_phy
  guides(fill=guide_legend(title="Time"),color=guide_legend(title="Time"))+
  xlab(paste0("PCoA 1 [10.9%]"))+
  ylab(paste0("PCoA 2 [6.9%]"))+
  scale_color_manual(values=col_time)+
  scale_fill_manual(values=col_time)+
  theme_bw()+theme(axis.line = element_line(),
                   panel.border = element_blank())+
  theme(text = element_text(size = 11))+
  annotate(geom="text", x=-0.3, y=0.42, size=3, color="black",
           label="italic(P)==0.001",parse=TRUE)
plot_Time_bray_curtis


# Bray_curtis_PC1 plot
Time_bary_PC1 <- ggplot(sample_data(phy_feces_clean),
                        aes(x=Time_new,y=PC1_bray,fill=Time_new))+
  stat_boxplot(geom="errorbar",aes(ymin=..ymax..),width=0.4)+
  stat_boxplot(geom="errorbar",aes(ymax=..ymin..),width=0.4)+
  geom_boxplot(varwidth = TRUE,notch=TRUE,notchwidth = 0.5, 
               outlier.shape=NA,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..),
               varwidth = TRUE,notch=TRUE,notchwidth = 0.5,outlier.shape=NA)+
  scale_color_manual(values=col_time)+
  scale_fill_manual(values=col_time)+
  coord_flip()+ 
  theme_void()+ #Empty theme, useful for plots with non-standard coordinates or for drawings
  guides(fill=guide_legend(title="Time"))+
  theme(text = element_text(size = 11))+
  stat_compare_means(comparisons = my_comparisons, paried=FALSE,
                     label.y = c(0.56,0.5,0.44,0.38,0.32))
Time_bary_PC1

# Wilcoxon test with FDR correction
Wil_test_adj_Time_bary_PC1 <- compare_means(PC1_bray ~ Time_new, 
                                            comparisons = my_comparisons, 
                                            p.adjust.method = "fdr",
                                            data=data.frame(sample_data(phy_feces_clean)))
Wil_test_adj_Time_bary_PC1

# LMM for Shannon as sensitivity
m_PC1_bray <- lmer(PC1_bray ~ Time_new +
                    oldchild01+catanddog+apartment_house+abbirth_mother+sex+delivery01+education+birthseason+
                    (1 | child_ID),
                  data = data.frame(sample_data(phy_feces_clean)),
                  REML = FALSE)

emm_PC1_bray <- emmeans(m_PC1_bray, ~ Time_new)

PC1_bray_pairwise <- contrast(emm_PC1_bray,method = "consec", adjust = "BH") %>%
  as.data.frame()
PC1_bray_pairwise



# Bray_curtis_PC2 plot
Time_bary_PC2 <- ggplot(data.frame(sample_data(phy_feces_clean)),
                        aes(x=Time_new,y=PC2_bray,fill=Time_new))+
  stat_boxplot(geom="errorbar",aes(ymin=..ymax..),width=0.4)+
  stat_boxplot(geom="errorbar",aes(ymax=..ymin..),width=0.4)+
  geom_boxplot(varwidth = TRUE,notch=TRUE,notchwidth = 0.5,
               outlier.shape=NA,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..),
               varwidth = TRUE,notch=TRUE,notchwidth = 0.5,outlier.shape=NA)+
  scale_color_manual(values=col_time)+
  scale_fill_manual(values=col_time)+
  theme_void()+#Empty theme, useful for plots with non-standard coordinates or for drawings
  guides(fill=guide_legend(title="Time"))+
  theme(text = element_text(size = 11))+
  stat_compare_means(comparisons = my_comparisons,label.y = c(0.54,0.48,0.42,0.36,0.3), paried=FALSE)
Time_bary_PC2

# Wilcoxon test with FDR correction
Wil_test_adj_Time_bary_PC2 <- compare_means(PC2_bray ~ Time_new, 
                                            comparisons = my_comparisons, 
                                            p.adjust.method = "fdr",
                                            data=data.frame(sample_data(phy_feces_clean)))
Wil_test_adj_Time_bary_PC2

# LMM for Shannon as sensitivity
m_PC2_bray <- lmer(PC2_bray ~ Time_new +
                     oldchild01+catanddog+apartment_house+abbirth_mother+sex+delivery01+education+birthseason+
                     (1 | child_ID),
                   data = data.frame(sample_data(phy_feces_clean)),
                   REML = FALSE)

emm_PC2_bray <- emmeans(m_PC2_bray, ~ Time_new)

PC2_bray_pairwise <- contrast(emm_PC2_bray,method = "consec", adjust = "BH") %>%
  as.data.frame()
PC2_bray_pairwise


## Arrange plots together
bray_PC1_PC2_Time_combine <- ggarrange(Time_bary_PC1 ,NULL,plot_Time_bray_curtis,Time_bary_PC2,
                                       ncol = 2,nrow=2,align = "hv",
                                       widths = c(2,1),heights = c(1,2),
                                       common.legend = TRUE,legend="right")
bray_PC1_PC2_Time_combine





##### Figure 1C. Beta-dispersion (distances to centroid) for each time point ####
bray_phy # bray curtis distance
ord_bray_phy # ordinate for bray curtis distance

# Create a data frame with time point information
Time_new_df <- data.frame(sample_data(phy_feces_clean)) %>%
  mutate(
    Time_new = factor(Time_new, levels = c("1w", "1m", "1y", "4y", "5y", "6y")),
    child_ID = factor(child_ID),
    oldchild01 = factor(oldchild01),
    catanddog = factor(catanddog),
    apartment_house = factor(apartment_house),
    abbirth_mother = factor(abbirth_mother),
    sex = factor(sex),
    delivery01 = factor(delivery01),
    education = factor(education),
    birthseason = factor(birthseason)
  )

# Compute group centroids using betadisper()
distance_centroid_Time <- betadisper(bray_phy,Time_new_df$Time_new)

# View the group centroids and anova test
distance_centroid_Time
# Add the centroid distances to the time_df data frame
Time_new_df$centroid_dist <- distance_centroid_Time$distances

# Beta diversity group final plot
bray_distance_centroid_Time_p <-ggplot(data = Time_new_df, 
                                       aes(x=Time_new, 
                                           y=centroid_dist,
                                           fill=Time_new))+
  stat_boxplot(geom="errorbar",aes(ymin=..ymax..),width=0.4)+
  stat_boxplot(geom="errorbar",aes(ymax=..ymin..),width=0.4)+
  geom_boxplot(varwidth = TRUE,notch=TRUE,notchwidth = 0.5, 
               outlier.shape=NA,linetype="dashed")+
  stat_boxplot(aes(ymin=..lower..,ymax=..upper..),
               varwidth = TRUE,notch=TRUE,notchwidth = 0.5,outlier.shape=NA)+
  scale_color_manual(values=col_time)+
  scale_fill_manual(values=col_time)+
  xlab("Time")+
  ylab("Distance to centroid")+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(legend.position = "none")+
  theme(text = element_text(size = 11))+
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(0.8,0.78,0.75,0.72,0.69), 
                     paried=FALSE) # wilcox.test
bray_distance_centroid_Time_p

# Statistical tests and p.adj test
Wil_test_adj_bray_distance_centroid_Time <- compare_means(
  centroid_dist ~ Time_new, 
  comparisons = my_comparisons, 
  p.adjust.method = "fdr",
  data=Time_new_df)
Wil_test_adj_bray_distance_centroid_Time


## 6. LMM sensitivity analysis for beta-dispersion
m_centroid <- lmer(centroid_dist ~ Time_new +
    oldchild01 + catanddog + apartment_house + abbirth_mother +
    sex + delivery01 + education + birthseason +
    (1 | child_ID),
  data = Time_new_df,
  REML = FALSE)

## Overall adjusted time effect
anova_centroid <- anova(m_centroid, type = 3)
anova_centroid

## Adjusted marginal means by time
emm_centroid <- emmeans(m_centroid, ~ Time_new)
emm_centroid

## Adjacent time-point comparisons after covariate adjustment
centroid_pairwise <- contrast(emm_centroid,method = "consec",adjust = "BH") %>%
  as.data.frame()
centroid_pairwise






######## Pairwise beta diversity by time (Bray–Curtis)
phy_feces_clean#

# Time pairs to compare
time_pairs <- list(c("1w", "1m"),c("1m", "1y"),c("1y", "4y"),c("4y", "5y"),c("5y", "6y"))

# Container for results
adonis_age_results <- list()

### Adjusted without 8 covariates
# Loop over pairs and run PERMANOVA (adonis2)
for (pair in time_pairs) {
  time1 <- pair[1]
  time2 <- pair[2]
  pair_label <- paste(time1, "vs", time2)
  
  # Subset samples for the current pair
  phy_sub <- subset_samples(phy_feces_clean, Time_new %in% c(time1, time2))
  phy_sub <- prune_taxa(taxa_sums(phy_sub) > 0, phy_sub)
  
  # Extract metadata
  meta_df <- data.frame(sample_data(phy_sub))
 
  # IMPORTANT: prune phyloseq object to match metadata after NA filtering
  phy_sub <- prune_samples(rownames(meta_df), phy_sub)
  phy_sub <- prune_taxa(taxa_sums(phy_sub) > 0, phy_sub)
  
  # Compute Bray–Curtis distance
  bray_dist <- phyloseq::distance(phy_sub, method = "bray", type = "samples")
  
  # PERMANOVA (adonis test)
  set.seed(33)
  adonis_res_age <- adonis2(bray_dist ~ Time_new, data = meta_df, permutations = 999)
 
  # Save R² and p-value
  adonis_age_results[[pair_label]] <- data.frame(
    Comparison = pair_label,
    R2 = adonis_res_age$R2[1],
    P_value = adonis_res_age$`Pr(>F)`[1]
  )
}

# Combine all results
adonis_age_table <- bind_rows(adonis_age_results)

# Add FDR-adjusted p-values (Benjamini–Hochberg)
adonis_age_table$FDR <- p.adjust(adonis_age_table$P_value, method = "fdr")

# Save as Excel
write_xlsx(adonis_age_table, path = "adonis_pairwise_age_results.xlsx")

# Print results
print(adonis_age_table)



### Adjusted with 8 covariates
# Loop over pairs and run PERMANOVA (adonis2)
for (pair in time_pairs) {
  time1 <- pair[1]
  time2 <- pair[2]
  pair_label <- paste(time1, "vs", time2)
  
  # Subset samples for the current pair
  phy_sub <- subset_samples(phy_feces_clean, Time_new %in% c(time1, time2))
  phy_sub <- prune_taxa(taxa_sums(phy_sub) > 0, phy_sub)
  
  
  # Extract metadata
  meta_df <- data.frame(sample_data(phy_sub))
  vars_needed <- c("Time_new", "oldchild01","catanddog","apartment_house",
                   "abbirth_mother","sex", "delivery01","education","birthseason")
  # Remove samples with missing values
  meta_df <- meta_df[complete.cases(meta_df[, vars_needed]), ]
  
  # IMPORTANT: prune phyloseq object to match metadata after NA filtering
  phy_sub <- prune_samples(rownames(meta_df), phy_sub)
  phy_sub <- prune_taxa(taxa_sums(phy_sub) > 0, phy_sub)
  
  # Re-extract metadata to ensure sample order matches phy_sub
  meta_df <- data.frame(sample_data(phy_sub))
  meta_df$Time_new <- factor(meta_df$Time_new, levels = c(time1, time2))
  
  # Compute Bray–Curtis distance
  bray_dist <- phyloseq::distance(phy_sub, method = "bray", type = "samples")
  
  # PERMANOVA (adonis test)
  set.seed(33)
 # adonis_res_age <- adonis2(bray_dist ~ Time_new, data = meta_df, permutations = 999)
  
  adonis_res_age <- adonis2(bray_dist ~ Time_new+
                              oldchild01 + catanddog + apartment_house + abbirth_mother +
                              sex + delivery01 + education + birthseason, 
                            data = meta_df, permutations = 999)
 
  # Save R² and p-value
  adonis_age_results[[pair_label]] <- data.frame(
    Comparison = pair_label,
    R2 = adonis_res_age$R2[1],
    P_value = adonis_res_age$`Pr(>F)`[1]
  )
}

# Combine all results
adonis_age_table <- bind_rows(adonis_age_results)

# Add FDR-adjusted p-values (Benjamini–Hochberg)
adonis_age_table$FDR <- p.adjust(adonis_age_table$P_value, method = "fdr")

# Save as Excel
write_xlsx(adonis_age_table, path = "adonis_pairwise_age_results.xlsx")

# Print results
print(adonis_age_table)






##### Supplementary Figure S1 ####

### verify and report LMM assumptions
models <- list(
  "Phylogenetic diversity" = m_fpd,
  "Shannon diversity" = m_Shannon,
  "PCoA1 Bray-Curtis" = m_PC1_bray,
  "PCoA2 Bray-Curtis" = m_PC2_bray,
  "Beta-dispersion" = m_centroid
)

pdf("LMM_diagnostic_plots.pdf", width = 10, height = 14)

old_par <- par(no.readonly = TRUE)

par(
  mfrow = c(length(models), 2),  # 5 rows × 2 columns
  mar = c(3.5, 3.5, 2.5, 1),
  cex = 0.75,
  cex.main = 0.9
)

for (model_name in names(models)) {
  
  m <- models[[model_name]]
  
  # Q-Q plot: residual normality
  qqnorm(resid(m),
         main = paste(model_name, "\nQ-Q plot"),
         xlab = "Theoretical quantiles",
         ylab = "Sample quantiles")
  qqline(resid(m), col = "red")
  
  # Residuals vs fitted: homoscedasticity
  plot(fitted(m),
       resid(m),
       xlab = "Fitted values",
       ylab = "Residuals",
       main = paste(model_name, "\nResiduals vs fitted"))
  abline(h = 0, col = "red")
}

par(old_par)

dev.off()






##### Figure 1D. Rate of change in Shannon diversity and identification of developmental turning points ####
alpha_phy_plot#

# Fit GAM (Shannon ~ smooth(age))
gam_model <- mgcv::gam(Shannon ~ s(AGE, bs = "cs"), data = alpha_phy_plot)

# First derivative of the smooth with 95% CI
derivs <- gratia::derivatives(gam_model, select = "s(AGE)", interval = "confidence")

# Compatibility shim: handle possible column name differences across gratia versions
nm <- names(derivs)
if ("derivative" %in% nm)  names(derivs)[names(derivs) == "derivative"] <- ".derivative"
if ("lower" %in% nm)       names(derivs)[names(derivs) == "lower"]      <- ".lower_ci"
if ("upper" %in% nm)       names(derivs)[names(derivs) == "upper"]      <- ".upper_ci"

# Keep the derivative dataframe tidy
deriv_df <- derivs %>%dplyr::select(AGE, .derivative, .lower_ci, .upper_ci)

# Turning points:
# 1) Maximum positive rate (fastest increase)
peak1 <- deriv_df %>% slice_max(order_by = .derivative, n = 1)

# 2) Maximum negative rate (fastest decrease)
peak3 <- deriv_df %>% slice_min(order_by = .derivative, n = 1)

# 3) Stabilization point: derivative ~ 0 and CI includes 0; choose the smallest |derivative|
deriv_df <- deriv_df %>%
  mutate(abs_deriv = abs(.derivative))

stable_zone <- deriv_df %>%
  filter(.lower_ci < 0 & .upper_ci > 0) %>% 
  slice_min(abs_deriv, n = 1)

# Combine turning points
turning_points <- bind_rows(
  peak1 %>% mutate(type = "max increase"),
  stable_zone %>% mutate(type = "stable"),
  peak3 %>% mutate(type = "max decrease")
)
print(turning_points)


# Plot first derivative with shaded 95% CI and vertical lines at turning points
p_GAM_Shannon_turningpoint_prediction <- ggplot(deriv_df, aes(x = AGE, y = .derivative)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "blue", alpha = 0.2) +
  geom_vline(data = turning_points, aes(xintercept = AGE, color = type),
             linetype = "dashed", size = 1) +
  scale_color_manual(values = c("max increase" = "red",
                                "stable" = "orange",
                                "max decrease" = "green")) +
  scale_x_continuous(breaks = seq(0, max(alpha_phy_plot$AGE, na.rm = TRUE), by = 365),
                     labels = c("0", "1y", "2y", "3y", "4y", "5y", "6y")) +
  labs(title = "First derivative of Shannon diversity: three turning points",
       x = "Age (days)", 
       y = "Rate of change in shannon diversity(d(Shannon)/d(Age))", 
       color = "Turning Point") +
  # d(Shannon)/d(Age) quantifies how much Shannon changes per infinitesimal increase in age
  theme_minimal(base_size = 14)
p_GAM_Shannon_turningpoint_prediction







##### Supplementary Figure S2. Rate of change in phylogenetic diversity ####

library(mgcv)
library(gratia)
library(dplyr)
library(ggplot2)

# Make sure AGE and PD are available
alpha_fpd2 <- alpha_fpd %>%
  filter(!is.na(AGE), !is.na(PD))

# Fit GAM: phylogenetic diversity ~ smooth(age)
gam_model_PD <- mgcv::gam(
  PD ~ s(AGE, bs = "cs"),
  data = alpha_fpd2,
  method = "REML"
)

# First derivative of the smooth with 95% CI
derivs_PD <- gratia::derivatives(
  gam_model_PD,
  select = "s(AGE)",
  interval = "confidence"
)

# Compatibility shim for different gratia versions
nm <- names(derivs_PD)
if ("derivative" %in% nm) names(derivs_PD)[names(derivs_PD) == "derivative"] <- ".derivative"
if ("lower" %in% nm)      names(derivs_PD)[names(derivs_PD) == "lower"]      <- ".lower_ci"
if ("upper" %in% nm)      names(derivs_PD)[names(derivs_PD) == "upper"]      <- ".upper_ci"

# Keep derivative dataframe tidy
deriv_PD_df <- derivs_PD %>%
  dplyr::select(AGE, .derivative, .lower_ci, .upper_ci) %>%
  mutate(abs_deriv = abs(.derivative))

# Turning points
peak_PD_increase <- deriv_PD_df %>%
  slice_max(order_by = .derivative, n = 1)

peak_PD_decrease <- deriv_PD_df %>%
  slice_min(order_by = .derivative, n = 1)

# Stabilization point:
# choose the point after the maximum increase where derivative is closest to zero and CI includes zero
stable_PD <- deriv_PD_df %>%
  filter(AGE > peak_PD_increase$AGE,
         .lower_ci < 0,
         .upper_ci > 0) %>%
  slice_min(abs_deriv, n = 1)

# Combine turning points
turning_points_PD <- bind_rows(peak_PD_increase %>% mutate(type = "max increase"),
                               stable_PD %>% mutate(type = "stable"),
                               peak_PD_decrease %>% mutate(type = "max decrease"))

print(turning_points_PD)

# Plot PD derivative
p_GAM_PD_turningpoint_prediction <- ggplot(deriv_PD_df,aes(x = AGE, y = .derivative)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci),
              fill = "blue",alpha = 0.2) +
  geom_vline(data = turning_points_PD,
             aes(xintercept = AGE, color = type),
             linetype = "dashed",linewidth = 1) +
  scale_color_manual(values = c("max increase" = "red",
                                "stable" = "orange",
                                "max decrease" = "green")) +
  scale_x_continuous(breaks = seq(0, max(alpha_phy_plot$AGE, na.rm = TRUE), by = 365),
                     labels = c("0", "1y", "2y", "3y", "4y", "5y", "6y")) +
  labs(title = "First derivative of phylogenetic diversity",
       x = "Age (days)",
       y = "Rate of change in phylogenetic diversity",
       color = "Turning point") +
  theme_minimal(base_size = 11)

p_GAM_PD_turningpoint_prediction






##### Supplementary Figure S3. Cohort-stratified benchmarking of gut microbiota maturation in COPSAC2010 and Roswall et al. ####
library(phyloseq)
library(biomformat)
library(Biostrings)
library(readxl)
library(ape)

required_files <- c("feature-table.biom", "taxonomy.tsv", "PRJEB38986_child_metadata_merged.xlsx",
                    "dna-sequences.fasta", "mapped_human.txt", "tree.nwk")

if (!all(file.exists(required_files))) {
  stop("Missing files in the current working directory:",
       paste(required_files[!file.exists(required_files)], collapse = ", "))
}

output_dir <- "phyloseq_results"
dir.create(output_dir, showWarnings = FALSE)


# 2. Read the ASV abundance table
# Rows: ASVs; columns: samples

counts <- as.matrix(biom_data(read_biom("feature-table.biom")))
storage.mode(counts) <- "numeric"

stopifnot(!is.null(rownames(counts)), !is.null(colnames(counts)), !anyDuplicated(rownames(counts)),
          !anyDuplicated(colnames(counts)), all(is.finite(counts)), all(counts >= 0),
          all(counts == round(counts)))

cat("Number of input ASVs:", nrow(counts), "\n")
cat("Number of input samples:", ncol(counts), "\n")


# 3. Read taxonomy and align by ASV ID
taxonomy <- read.delim("taxonomy.tsv", sep = "\t", quote = "", comment.char = "", check.names = FALSE,
                       stringsAsFactors = FALSE)

stopifnot(all(c("Feature ID", "Taxon") %in% names(taxonomy)), !anyDuplicated(taxonomy[["Feature ID"]]))

if (!all(rownames(counts) %in% taxonomy[["Feature ID"]])) {
  stop("Taxonomy is missing for some ASVs.")
}

taxonomy <- taxonomy[match(rownames(counts), taxonomy[["Feature ID"]]), , drop = FALSE]


# 4. Remove ASVs on the human list, mitochondria, and chloroplasts
# Human ASVs have already been removed from the current BIOM; rechecking does not remove them twice
human_ids <- unique(trimws(readLines("mapped_human.txt", warn = FALSE)))
human_ids <- human_ids[nzchar(human_ids)]

taxon_text <- taxonomy$Taxon
taxon_text[is.na(taxon_text)] <- ""

is_human <- rownames(counts) %in% human_ids

is_organelle <- grepl("mitochondria|chloroplast", taxon_text, ignore.case = TRUE)

remove <- is_human | is_organelle

write.csv(data.frame(ASV = rownames(counts)[remove], Taxon = taxon_text[remove],
                     Human_list = is_human[remove], Mitochondria_or_chloroplast = is_organelle[remove],
                     Total_count = rowSums(counts)[remove]), file.path(output_dir, "removed_ASVs.csv"), row.names = FALSE)

cat("ASVs on the human list in the current table:", sum(is_human), "\n")
cat("Mitochondrial/chloroplast ASVs:", sum(is_organelle), "\n")

counts <- counts[!remove, , drop = FALSE]
taxonomy <- taxonomy[!remove, , drop = FALSE]

keep <- rowSums(counts) > 0
counts <- counts[keep, , drop = FALSE]
taxonomy <- taxonomy[keep, , drop = FALSE]

if (nrow(counts) == 0) {
  stop("No ASVs remain after filtering.")
}

# 5. Split SILVA taxonomy into seven ranks
library(stringr)

# Confirm that the ASV order matches
stopifnot(identical(rownames(counts), taxonomy[["Feature ID"]]))

# Original taxonomy strings
tax_text <- taxonomy$Taxon

# Extract the seven ranks separately; missing ranks remain NA
Kingdom <- str_match(tax_text, "(?:^|;)\\s*D_0__([^;]*)")[, 2]
Phylum  <- str_match(tax_text, "(?:^|;)\\s*D_1__([^;]*)")[, 2]
Class   <- str_match(tax_text, "(?:^|;)\\s*D_2__([^;]*)")[, 2]
Order   <- str_match(tax_text, "(?:^|;)\\s*D_3__([^;]*)")[, 2]
Family  <- str_match(tax_text, "(?:^|;)\\s*D_4__([^;]*)")[, 2]
Genus   <- str_match(tax_text, "(?:^|;)\\s*D_5__([^;]*)")[, 2]
Species <- str_match(tax_text, "(?:^|;)\\s*D_6__([^;]*)")[, 2]

# Combine
tax_matrix <- cbind(Kingdom, Phylum, Class, Order, Family, Genus, Species)

# Trim whitespace while preserving the matrix structure
tax_matrix[] <- trimws(tax_matrix)

# Set empty labels and ambiguous taxonomy to NA
tax_matrix[tax_matrix == "" & !is.na(tax_matrix)] <- NA_character_
tax_matrix[tax_matrix == "Ambiguous_taxa" &
             !is.na(tax_matrix)] <- NA_character_

# Add ASV names
rownames(tax_matrix) <- taxonomy[["Feature ID"]]

# Check
head(tax_matrix)
dim(tax_matrix)
colSums(!is.na(tax_matrix))

stopifnot(identical(rownames(counts), rownames(tax_matrix)))


# 6. Read metadata from Excel
# Sample ID column: run_accession
metadata <- as.data.frame(read_excel("PRJEB38986_child_metadata_merged.xlsx", sheet = 1,
                                     col_types = "text", .name_repair = "minimal"), stringsAsFactors = FALSE)

stopifnot(all(c("run_accession", "child_id", "age_month", "Time_new", "sample_type"
) %in% names(metadata)))

metadata$run_accession <- trimws(metadata$run_accession)

if (anyNA(metadata$run_accession) ||
    any(!nzchar(metadata$run_accession)) ||
    anyDuplicated(metadata$run_accession)) {
  stop("The metadata run_accession column contains missing, empty, or duplicate values.")
}

metadata$child_id <- as.character(metadata$child_id)

age_original <- trimws(metadata$age_month)
metadata$age_month <- suppressWarnings(as.numeric(age_original))

bad_age <- !is.na(age_original) &
  nzchar(age_original) &
  is.na(metadata$age_month)

if (any(bad_age)) {
  stop("Some age_month values cannot be converted to numeric; please check.")
}

# Record unmatched samples
writeLines(setdiff(colnames(counts), metadata$run_accession),
           file.path(output_dir, "BIOM_samples_without_metadata.txt"))

writeLines(setdiff(metadata$run_accession, colnames(counts)),
           file.path(output_dir, "metadata_samples_without_BIOM.txt"))

# Retain all sequenced samples and fill unmatched metadata with NA
idx <- match(colnames(counts), metadata$run_accession)
metadata_all <- metadata[idx, , drop = FALSE]

metadata_all$run_accession <- colnames(counts)
metadata_all$metadata_matched <- !is.na(idx)
rownames(metadata_all) <- colnames(counts)

cat("Samples with matching metadata:", sum(!is.na(idx)), "\n")
cat("Samples without matching metadata:", sum(is.na(idx)), "\n")


# 7. Read representative sequences and the phylogenetic tree
seqs <- readDNAStringSet("dna-sequences.fasta")
names(seqs) <- sub("\\s.*$", "", names(seqs))

tree <- read.tree("tree.nwk")

stopifnot(!anyDuplicated(names(seqs)), !anyDuplicated(tree$tip.label))

final_ids <- rownames(counts)

if (!all(final_ids %in% names(seqs))) {
  stop("Some retained ASVs are missing from the FASTA file.")
}

if (!all(final_ids %in% tree$tip.label)) {
  stop("Some retained ASVs are missing from the phylogenetic tree.")
}

seqs <- seqs[final_ids]
tree <- keep.tip(tree, final_ids)


# 8. Build a phyloseq object containing all samples
ps_all <- phyloseq(otu_table(counts, taxa_are_rows = TRUE), tax_table(tax_matrix),
                   sample_data(metadata_all), phy_tree(tree), seqs)

# Check for unintended loss of samples or ASVs during merging
stopifnot(nsamples(ps_all) == ncol(counts), ntaxa(ps_all) == nrow(counts))

# Record and exclude zero-count samples
zero_samples <- sample_names(ps_all)[sample_sums(ps_all) == 0]

writeLines(zero_samples, file.path(output_dir, "zero_count_samples.txt"))

if (length(zero_samples) > 0) {
  ps_all <- prune_samples(sample_sums(ps_all) > 0, ps_all)
}

ps_all <- prune_taxa(taxa_sums(ps_all) > 0, ps_all)


# 9. Build the analysis object for samples with metadata
matched_ids <- rownames(metadata_all)[metadata_all$metadata_matched]
matched_ids <- intersect(matched_ids, sample_names(ps_all))

if (length(matched_ids) == 0) {
  stop("No samples match the metadata.")
}

ps_matched <- prune_samples(matched_ids, ps_all)
ps_matched <- prune_taxa(taxa_sums(ps_matched) > 0, ps_matched)

print(ps_all) # 2033 children+mother
print(ps_matched) # 1676 children

md <- as(sample_data(ps_matched), "data.frame")

table(md$sample_type, useNA = "ifany")
table(md$Time_new, useNA = "ifany")
summary(md$age_month)
summary(sample_sums(ps_matched))

ps_matched #

colnames(sample_data(ps_matched))
View(otu_table(ps_matched))
View(tax_table(ps_matched))

# Extract the taxonomy matrix
tax <- as(tax_table(ps_matched), "matrix")

# Replace missing values and empty strings
tax[is.na(tax)] <- "Unknown"
tax[trimws(tax) == ""] <- "Unknown"

# Assign the taxonomy back to phyloseq
tax_table(ps_matched) <- tax_table(tax)

# Check
head(tax_table(ps_matched))
sum(is.na(tax_table(ps_matched)))



#### Figure S3A. Study design timeline
timeline_df <- dplyr::bind_rows(data.frame(Cohort = "COPSAC2010",
                                           Age_year = c(7/365, 30/365, 1, 4, 5, 6), Time_label = c("1w", "1m", "1y", "4y", "5y", "6y")),
                                data.frame(Cohort = "Roswall et al.", Age_year = c(0, 4/12, 1, 3, 5),
                                           Time_label = c("0m", "4m", "1y", "3y", "5y")))

timeline_df$Cohort <- factor(timeline_df$Cohort, levels = c("Roswall et al.","COPSAC2010"))

timeline_segment <- timeline_df %>%
  group_by(Cohort) %>%
  summarise(start = min(Age_year), end = max(Age_year), .groups = "drop")

p_timeline <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = timeline_segment,
                        ggplot2::aes(x = start, xend = end, y = Cohort, yend = Cohort,colour = Cohort), linewidth = 1 ) +
  ggplot2::geom_point(data = timeline_df, ggplot2::aes(x = Age_year, y = Cohort,colour = Cohort),
                      size = 3) +
  ggplot2::geom_text(data = timeline_df, ggplot2::aes(x = Age_year, y = Cohort, label = Time_label),
                     vjust = -0.8, size = 3) +
  ggplot2::scale_colour_manual(values = c("#b15928","#0072B2")) +
  ggplot2::scale_fill_manual(values = c("#b15928","#0072B2")) +
  ggplot2::scale_x_continuous(breaks = 0:6, labels = c("0", "1y", "2y", "3y", "4y", "5y", "6y"),
                              limits = c(-0.1, 6.2)) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_blank()) +
  ggplot2::theme(legend.position = "none")+
  ggplot2::xlab("Age") +
  ggplot2::ylab(NULL) +
  ggplot2::ggtitle("Study design and sampling windows")
p_timeline


#### Figure S3B. Annualized Faith’s PD gain 
calc_faith_pd <- function(ps, cohort_name, id_col, time_levels, age_map) {
  
  otu_mat <- as(otu_table(ps), "matrix")
  
  ## picante::pd requires samples as rows and taxa as columns
  if (taxa_are_rows(ps)) {
    otu_sample_taxa <- t(otu_mat)
  } else {
    otu_sample_taxa <- otu_mat
  }
  
  pd_df <- picante::pd(otu_sample_taxa, phy_tree(ps), include.root = FALSE)
  
  ## Extract metadata safely
  meta_df <- data.frame(sample_data(ps), check.names = FALSE)
  
  ## Avoid duplicated SampleID column
  ## If metadata already has SampleID, keep it;
  ## otherwise create SampleID from rownames.
  if (!"SampleID" %in% colnames(meta_df)) {
    meta_df$SampleID <- rownames(meta_df)
  }
  
  ## Make sure SampleID is character
  meta_df$SampleID <- as.character(meta_df$SampleID)
  
  ## Prepare PD table
  pd_out <- data.frame(SampleID = rownames(pd_df), PD = pd_df$PD, stringsAsFactors = FALSE)
  
  ## Merge PD with metadata
  out <- pd_out %>%
    left_join(meta_df, by = "SampleID") %>%
    mutate(Cohort = cohort_name, Time_new = factor(as.character(Time_new), levels = time_levels),
           Age_year = unname(age_map[as.character(Time_new)]), child_id_for_model = factor(.data[[id_col]]))
  
  return(out)
}

age_map_cop <- c("1w" = 7 / 365.25, "1m" = 1 / 12, "1y" = 1, "4y" = 4, "5y" = 5, "6y" = 6)

age_map_ros <- c("0m" = 0, "4m" = 4 / 12, "1y" = 1, "3y" = 3, "5y" = 5)

pd_cop <- calc_faith_pd(ps = ps_cop, cohort_name = "COPSAC2010", id_col = "child_ID",
                        time_levels = c("1w", "1m", "1y", "4y", "5y", "6y"), age_map = age_map_cop)

pd_ros <- calc_faith_pd(ps = ps_ros, cohort_name = "Roswall et al.", id_col = "child_id",
                        time_levels = c("0m", "4m", "1y", "3y", "5y"), age_map = age_map_ros)


table(pd_cop$Time_new, useNA = "ifany")
table(pd_ros$Time_new, useNA = "ifany")
head(pd_cop[, c("SampleID", "PD", "Time_new", "Age_year", "child_ID")])
head(pd_ros[, c("SampleID", "PD", "Time_new", "Age_year", "child_id")])


### B2. Calculate annualized PD gains between adjacent time points
## Function: annualized PD gain
get_annualized_pd_gain <- function(pd_df, cohort_name, age_map) {
  
  pd_df <- pd_df %>%
    filter(!is.na(PD), !is.na(Time_new), !is.na(child_id_for_model))
  
  m <- lmer(PD ~ Time_new + (1 | child_id_for_model), data = pd_df, REML = FALSE)
  
  emm <- emmeans(m, ~ Time_new)
  
  pair_df <- contrast(emm, method = "consec", adjust = "BH") %>%
    as.data.frame()
  
  ## contrast labels are usually like "1m - 1w"
  pair_df2 <- pair_df %>%
    tidyr::separate(contrast, into = c("later", "earlier"), sep = " - ", remove = FALSE) %>%
    mutate(Cohort = cohort_name, age_later = unname(age_map[later]),
           age_earlier = unname(age_map[earlier]), interval_years = age_later - age_earlier,
           interval_mid = (age_later + age_earlier) / 2, annualized_estimate = estimate / interval_years,
           annualized_SE = SE / interval_years,
           annualized_lower = annualized_estimate - 1.96 * annualized_SE,
           annualized_upper = annualized_estimate + 1.96 * annualized_SE,
           interval_label = paste0(earlier, "→", later), p_adj = p.value)
  
  return(list(model = m, emmeans = emm, pairwise = pair_df2))
}

pd_gain_cop_res <- get_annualized_pd_gain(pd_df = pd_cop, cohort_name = "COPSAC2010",
                                          age_map = age_map_cop)

pd_gain_ros_res <- get_annualized_pd_gain( pd_df = pd_ros, cohort_name = "Roswall et al.",
                                           age_map = age_map_ros)

pd_gain_all <- bind_rows(pd_gain_cop_res$pairwise, pd_gain_ros_res$pairwise)

pd_gain_all

#write_xlsx(list(Annualized_PD_gain = pd_gain_all,
COPSAC_PD_pairwise = pd_gain_cop_res$pairwise, Roswall_PD_pairwise = pd_gain_ros_res$pairwise),
path = file.path(OUT_COMPARE, "B_annualized_PD_gain_results.xlsx"))


# Generate significance labels
pd_gain_all <- pd_gain_all %>%
  dplyr::mutate(significance = dplyr::case_when(is.na(p_adj)  ~ NA_character_, p_adj < 0.001 ~ "***",
                                                p_adj < 0.01  ~ "**", p_adj < 0.05  ~ "*",
                                                TRUE         ~ "ns"),plot_label = paste0(interval_label, "\n", significance))

# Plot
p_pd_gain <- ggplot2::ggplot(pd_gain_all, ggplot2::aes(x = interval_mid, y = annualized_estimate,
                                                       color = Cohort, group = Cohort)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = annualized_lower, ymax = annualized_upper), width = 0.08,
                         linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(y = annualized_upper, label = plot_label), vjust = -0.4, size = 3,
                     show.legend = FALSE, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = 0:6, labels = c("0", "1y", "2y", "3y", "4y", "5y", "6y"),
                              limits = c(0, 6)) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.20)) ) +
  ggplot2::scale_colour_manual(values = c("COPSAC2010" = "#0072B2", "Roswall et al." = "#b15928") ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "top") +
  ggplot2::labs( x = "Midpoint age of interval", y = "Change in Faith's PD", title = "Faith's PD gain"#,
                 #caption = "BH-adjusted P: *** < 0.001; ** < 0.01; * < 0.05; ns ≥ 0.05"
  )

p_pd_gain


#### Figure S3C. Annualized Shannon gain 
library(phyloseq)
library(dplyr)
library(lmerTest)
library(emmeans)

# 1. Map time points to age in years
# If age_map_cop and age_map_ros are already defined, retain their existing values
age_map_cop <- c("1w" = 7 / 365.25, "1m" = 1 / 12, "1y" = 1, "4y" = 4, "5y" = 5, "6y" = 6)

age_map_ros <- c("0m" = 0, "4m" = 4 / 12, "1y" = 1, "3y" = 3, "5y" = 5)

# 2. COPSAC2010 Shannon
shannon_cop <- phyloseq::estimate_richness(ps_cop, measures = "Shannon")

meta_cop <- data.frame(sample_data(ps_cop), check.names = FALSE)

shannon_cop <- shannon_cop[rownames(meta_cop), , drop = FALSE]

shannon_cop_df <- meta_cop
shannon_cop_df$SampleID <- rownames(meta_cop)
shannon_cop_df$Shannon <- shannon_cop$Shannon
shannon_cop_df$Cohort <- "COPSAC2010"

shannon_cop_df$Time_new <- factor(as.character(meta_cop$Time_new), levels = names(age_map_cop))

shannon_cop_df$child_id_for_model <- factor(meta_cop$child_ID)


# 3. Roswall / H2GS Shannon
shannon_ros <- phyloseq::estimate_richness(ps_ros, measures = "Shannon")

meta_ros <- data.frame(sample_data(ps_ros), check.names = FALSE)

shannon_ros <- shannon_ros[rownames(meta_ros), , drop = FALSE]

shannon_ros_df <- meta_ros
shannon_ros_df$SampleID <- rownames(meta_ros)
shannon_ros_df$Shannon <- shannon_ros$Shannon
shannon_ros_df$Cohort <- "H2GS"

# Standardize NB to 0m
ros_time <- as.character(meta_ros$Time_new)
ros_time[ros_time == "NB"] <- "0m"

shannon_ros_df$Time_new <- factor(ros_time, levels = names(age_map_ros))

shannon_ros_df$child_id_for_model <- factor(meta_ros$child_id)


# Check time points and Shannon diversity
table(shannon_cop_df$Time_new, useNA = "ifany")
table(shannon_ros_df$Time_new, useNA = "ifany")

summary(shannon_cop_df$Shannon)
summary(shannon_ros_df$Shannon)


### Fit models for both cohorts:
# 4. Retain samples with complete model variables
cop_model_data <- shannon_cop_df %>%
  filter(is.finite(Shannon), !is.na(Time_new), !is.na(child_id_for_model))

ros_model_data <- shannon_ros_df %>%
  filter(is.finite(Shannon), !is.na(Time_new), !is.na(child_id_for_model))

# Confirm that every planned time point has samples
stopifnot(all(table(cop_model_data$Time_new) > 0), all(table(ros_model_data$Time_new) > 0))

# Use child ID as a random intercept
m_shannon_cop <- lmerTest::lmer(Shannon ~ Time_new + (1 | child_id_for_model), data = cop_model_data,
                                REML = FALSE)

m_shannon_ros <- lmerTest::lmer(Shannon ~ Time_new + (1 | child_id_for_model), data = ros_model_data,
                                REML = FALSE)


# 5. Estimated marginal means at each time point
emm_shannon_cop <- emmeans::emmeans(m_shannon_cop, ~ Time_new, lmer.df = "satterthwaite")

emm_shannon_ros <- emmeans::emmeans(m_shannon_ros, ~ Time_new, lmer.df = "satterthwaite")

# Adjacent time points: later minus earlier
contrast_cop <- emmeans::contrast(emm_shannon_cop, method = "consec")

contrast_ros <- emmeans::contrast(emm_shannon_ros, method = "consec")

# Extract unadjusted 95% confidence intervals and P values
gain_cop <- as.data.frame(summary(contrast_cop, infer = c(TRUE, TRUE), adjust = "none"))

gain_ros <- as.data.frame(summary(contrast_ros, infer = c(TRUE, TRUE), adjust = "none"))

# Apply BH correction separately within each cohort
gain_cop$p_adj <- p.adjust(gain_cop$p.value, method = "BH")
gain_ros$p_adj <- p.adjust(gain_ros$p.value, method = "BH")


### Divide by interval length to obtain annualized changes:
# 6. Adjacent intervals in COPSAC2010
gain_cop$Cohort <- "COPSAC2010"
gain_cop$earlier <- head(names(age_map_cop), -1)
gain_cop$later <- tail(names(age_map_cop), -1)

gain_cop$age_earlier <- unname(age_map_cop[gain_cop$earlier])
gain_cop$age_later <- unname(age_map_cop[gain_cop$later])


# 7. Adjacent intervals in Roswall et al.
gain_ros$Cohort <- "Roswall et al."
gain_ros$earlier <- head(names(age_map_ros), -1)
gain_ros$later <- tail(names(age_map_ros), -1)

gain_ros$age_earlier <- unname(age_map_ros[gain_ros$earlier])
gain_ros$age_later <- unname(age_map_ros[gain_ros$later])


# 8. Combine results and calculate annualized changes
shannon_gain_all <- bind_rows(gain_cop, gain_ros) %>%
  mutate(interval_years = age_later - age_earlier, interval_mid = (age_later + age_earlier) / 2,
         
         annualized_estimate = estimate / interval_years, annualized_SE = SE / interval_years,
         
         # Directly scale the model-derived 95% confidence intervals
         annualized_lower = lower.CL / interval_years, annualized_upper = upper.CL / interval_years,
         
         interval_label = paste0(earlier, "→", later),
         
         significance = case_when(is.na(p_adj) ~ NA_character_, p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                                  p_adj < 0.05 ~ "*", TRUE ~ "ns"),
         
         plot_label = paste0(interval_label, "\n", significance))

shannon_gain_all


### Plot
# 9. Annualized Shannon change trajectory
p_shannon_gain <- ggplot2::ggplot(shannon_gain_all, ggplot2::aes(x = interval_mid,
                                                                 y = annualized_estimate, colour = Cohort, group = Cohort)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = annualized_lower, ymax = annualized_upper), width = 0.08,
                         linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(y = annualized_upper, label = plot_label), vjust = -0.4, size = 3,
                     show.legend = FALSE, na.rm = TRUE ) +
  ggplot2::scale_x_continuous(breaks = 0:6, labels = c("0", "1y", "2y", "3y", "4y", "5y", "6y"),
                              limits = c(0, 6)) +
  ggplot2::scale_y_continuous( expand = ggplot2::expansion(mult = c(0.05, 0.20))) +
  ggplot2::scale_colour_manual(values = c( "COPSAC2010" = "#0072B2", "Roswall et al." = "#b15928"  ) ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "top") +
  ggplot2::labs(x = "Midpoint age of interval", y = "Change in Shannon diversity",
                title = "Shannon diversity change",
                #caption = paste0("Error bars: pointwise 95% CI. ",
                #                 "Within-cohort BH-adjusted P:\n",
                #                 "*** < 0.001; ** < 0.01; * < 0.05; ns ≥ 0.05")
  )
p_shannon_gain


# Save to the current directory
write.csv(shannon_gain_all, "B_annualized_Shannon_gain_results.csv", row.names = FALSE)


#### Figure S3D. Bray-Curtis / beta-dispersion comparison 

### D1. Calculate beta-dispersion and PERMANOVA

get_beta_dispersion <- function(ps, cohort_name, id_col, time_levels, age_map) {
  
  ps2 <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  # Read metadata only once
  meta_df <- data.frame(sample_data(ps2), check.names = FALSE)
  
  # Use sample row names directly without adding a duplicate column
  meta_df$SampleID <- rownames(meta_df)
  
  meta_df <- meta_df %>%
    mutate(Cohort = cohort_name, Time_new = factor(as.character(Time_new), levels = time_levels),
           Age_year = unname(age_map[as.character(Time_new)]), child_id_for_model = factor(.data[[id_col]])
    ) %>%
    filter(!is.na(Time_new), !is.na(Age_year), !is.na(child_id_for_model))
  
  # Subset phyloseq to the same samples
  ps2 <- prune_samples(meta_df$SampleID, ps2)
  ps2 <- prune_taxa(taxa_sums(ps2) > 0, ps2)
  
  # Align metadata with the phyloseq sample order
  meta_df <- meta_df[match(sample_names(ps2), meta_df$SampleID), , drop = FALSE]
  
  rownames(meta_df) <- meta_df$SampleID
  
  stopifnot(identical(meta_df$SampleID, sample_names(ps2)))
  
  # Resume the second half of the original code here
  bray <- phyloseq::distance(ps2, method = "bray", type = "samples")
  
  # Retain the original PERMANOVA, betadisper,
  # mixed-model, summary_df, and return code below,
  # and retain the closing } at the end of the function
  
  bray <- phyloseq::distance(ps2, method = "bray", type = "samples")
  
  ## PERMANOVA: within-cohort age structure
  set.seed(33)
  adonis_res <- vegan::adonis2(bray ~ Time_new, data = meta_df, permutations = 999,
                               strata = meta_df$child_id_for_model)
  
  disp <- vegan::betadisper(bray, meta_df$Time_new)
  meta_df$centroid_dist <- disp$distances
  
  m_disp <- lmer( centroid_dist ~ Time_new + (1 | child_id_for_model), data = meta_df, REML = FALSE)
  
  emm_disp <- emmeans(m_disp, ~ Time_new)
  
  pair_disp <- contrast(emm_disp, method = "consec", adjust = "BH") %>%
    as.data.frame() %>%
    tidyr::separate(contrast, into = c("later", "earlier"), sep = " - ", remove = FALSE) %>%
    mutate(Cohort = cohort_name, age_later = unname(age_map[later]),
           age_earlier = unname(age_map[earlier]), interval_years = age_later - age_earlier,
           interval_mid = (age_later + age_earlier) / 2, annualized_estimate = estimate / interval_years,
           annualized_SE = SE / interval_years,
           annualized_lower = annualized_estimate - 1.96 * annualized_SE,
           annualized_upper = annualized_estimate + 1.96 * annualized_SE,
           interval_label = paste0(earlier, "→", later), p_adj = p.value)
  
  summary_df <- meta_df %>%
    group_by(Cohort, Time_new, Age_year) %>%
    summarise(n = n(), median_centroid = median(centroid_dist, na.rm = TRUE),
              mean_centroid = mean(centroid_dist, na.rm = TRUE),
              se_centroid = sd(centroid_dist, na.rm = TRUE) / sqrt(n()), .groups = "drop")
  
  return(list(ps = ps2, bray = bray, meta = meta_df, adonis = adonis_res, dispersion_model = m_disp,
              dispersion_pairwise = pair_disp, dispersion_summary = summary_df))
}

beta_cop <- get_beta_dispersion(ps = ps_cop, cohort_name = "COPSAC2010", id_col = "child_ID",
                                time_levels = c("1w", "1m", "1y", "4y", "5y", "6y"),
                                age_map = age_map_cop)

beta_ros <- get_beta_dispersion(ps = ps_ros, cohort_name = "Roswall et al.", id_col = "child_id",
                                time_levels = c("0m", "4m", "1y", "3y", "5y"), age_map = age_map_ros)

beta_disp_all <- bind_rows(beta_cop$dispersion_pairwise, beta_ros$dispersion_pairwise)

beta_summary_all <- bind_rows(beta_cop$dispersion_summary, beta_ros$dispersion_summary)

beta_disp_all

write_xlsx(list(Beta_dispersion_annualized_change = beta_disp_all,
                Beta_dispersion_summary = beta_summary_all),
           path = file.path(OUT_COMPARE, "C_beta_dispersion_results.xlsx"))

beta_cop$adonis
beta_ros$adonis


### D2. Plot the original beta-dispersion trend
# This plot shows age-related changes in distance-to-centroid within each cohort. Focus on trends rather than absolute differences between cohorts.
p_beta_disp_trend <- ggplot2::ggplot(beta_summary_all,
                                     ggplot2::aes(x = Age_year, y = mean_centroid, color = Cohort, group = Cohort)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_centroid - 1.96 * se_centroid,
                                      ymax = mean_centroid + 1.96 * se_centroid), width = 0.08,
                         linewidth = 0.5) +
  ggplot2::scale_x_continuous(breaks = 0:6, labels = c("0", "1y", "2y", "3y", "4y", "5y", "6y"),
                              limits = c(0, 6)) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::scale_colour_manual(values = c("COPSAC2010" = "#0072B2","Roswall et al." = "#b15928") ) +
  ggplot2::xlab("Age") +
  ggplot2::ylab("Mean distance to centroid") +
  ggplot2::ggtitle("Age-related changes in beta-dispersion")
p_beta_disp_trend


### D3. Plot annualized changes in beta-dispersion
# This plot shows whether dispersion increases or decreases.

# Generate P-value labels
beta_disp_all <- beta_disp_all %>%
  dplyr::mutate(p_label = dplyr::case_when(is.na(p_adj) ~ "", p_adj < 0.001 ~ "BH P < 0.001",
                                           TRUE ~ paste0("BH P = ", sprintf("%.3f", p_adj))),
                plot_label = paste0(interval_label, "\n", p_label) )

# Plot
p_beta_disp_change <- ggplot2::ggplot(beta_disp_all, ggplot2::aes(x = interval_mid,
                                                                  y = annualized_estimate, colour = Cohort, group = Cohort )) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = annualized_lower, ymax = annualized_upper), width = 0.08,
                         linewidth = 0.5 ) +
  ggplot2::geom_text(ggplot2::aes(y = annualized_upper, label = plot_label), vjust = -0.4, size = 3,
                     show.legend = FALSE) +
  ggplot2::scale_x_continuous(breaks = 0:6, labels = c("0", "1y", "2y", "3y", "4y", "5y", "6y"),
                              limits = c(0, 6)) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.25))) +
  ggplot2::scale_colour_manual(values = c("COPSAC2010" = "#0072B2","Roswall et al." = "#b15928")) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme( panel.grid.minor = ggplot2::element_blank(), legend.position = "top") +
  ggplot2::labs(x = "Midpoint age of interval", y = "Change in beta-dispersion",
                title = "Change in beta-dispersion"#,
                #caption = "P values: adjacent age-group contrasts, BH-adjusted within each cohort."
  )
p_beta_disp_change

# Save
ggplot2::ggsave(filename = file.path(OUT_COMPARE, "C_annualized_beta_dispersion_change.pdf"),
                plot = p_beta_disp_change, width = 9, height = 5)


#### Figure S3E. Genus-level developmental trajectories 
### E1. Extract relative abundances of target genera
## Function: target genus trajectory

clean_genus_name <- function(x) {
  x <- as.character(x)
  x <- gsub("^g__", "", x)
  x <- gsub("^D_[0-9]__", "", x)
  x <- trimws(x)
  x[x == "" | is.na(x)] <- "Unclassified"
  return(x)
}

## Agglomerate to genus
rank_names(ps_cop)
rank_names(ps_ros)

head(as(tax_table(ps_cop), "matrix"), 3)
head(as(tax_table(ps_ros), "matrix"), 3)

## Standardize the COPSAC column name from lowercase genus to Genus:
tax_cop <- as(tax_table(ps_cop), "matrix")

# Trim whitespace from column names
colnames(tax_cop) <- trimws(colnames(tax_cop))

# Standardize genus, GENUS, and other case variants to Genus
colnames(tax_cop)[tolower(colnames(tax_cop)) == "genus"] <- "Genus"

tax_table(ps_cop) <- tax_table(tax_cop)

# Check
rank_names(ps_cop)


##
make_target_genus_long <- function(ps, cohort_name, id_col, time_levels, age_map, target_genera) {
  
  # Agglomerate to genus level, then calculate relative abundance
  ps_genus <- phyloseq::tax_glom(ps, taxrank = "Genus", NArm = FALSE)
  
  # Relative abundance cannot be calculated for zero-count samples
  ps_genus <- phyloseq::prune_samples(phyloseq::sample_sums(ps_genus) > 0, ps_genus)
  
  ps_genus_rel <- phyloseq::transform_sample_counts(ps_genus, function(x) x / sum(x) )
  
  # Abundance matrix: samples as rows
  otu_mat <- as(phyloseq::otu_table(ps_genus_rel), "matrix")
  
  if (phyloseq::taxa_are_rows(ps_genus_rel)) {
    otu_mat <- t(otu_mat)
  }
  
  otu_df <- as.data.frame(otu_mat, check.names = FALSE)
  otu_df$SampleID <- rownames(otu_df)
  
  # Taxonomy information
  tax_df <- as.data.frame(as(phyloseq::tax_table(ps_genus_rel), "matrix"), stringsAsFactors = FALSE)
  
  tax_df$FeatureID <- rownames(tax_df)
  tax_df$Genus_clean <- clean_genus_name(tax_df$Genus)
  
  # Sample metadata: overwrite SampleID directly without adding a duplicate column
  meta_df <- data.frame(phyloseq::sample_data(ps_genus_rel), check.names = FALSE)
  
  meta_df$SampleID <- rownames(meta_df)
  
  meta_df <- meta_df %>%
    dplyr::mutate(Cohort = cohort_name, Time_new = factor(as.character(Time_new), levels = time_levels),
                  Age_year = unname(age_map[as.character(Time_new)]),
                  child_id_for_model = factor(.data[[id_col]]))
  
  # Convert to long format and extract target genera
  genus_long <- otu_df %>%
    tidyr::pivot_longer(cols = -SampleID, names_to = "FeatureID", values_to = "Relative_abundance") %>%
    dplyr::left_join(tax_df[, c("FeatureID", "Genus_clean")], by = "FeatureID") %>%
    dplyr::filter(Genus_clean %in% target_genera) %>%
    dplyr::group_by(SampleID, Genus_clean) %>%
    dplyr::summarise(Relative_abundance = sum(Relative_abundance), .groups = "drop")
  
  # Complete missing sample-genus combinations with zero abundance
  complete_grid <- expand.grid(SampleID = phyloseq::sample_names(ps_genus_rel),
                               Genus_clean = target_genera, stringsAsFactors = FALSE)
  
  genus_long2 <- complete_grid %>%
    dplyr::left_join(genus_long, by = c("SampleID", "Genus_clean")) %>%
    dplyr::mutate(Relative_abundance = dplyr::coalesce(Relative_abundance, 0)) %>%
    dplyr::left_join(meta_df, by = "SampleID") %>%
    dplyr::mutate(Genus_clean = factor(Genus_clean, levels = target_genera))
  
  return(genus_long2)
}

target_genera <- c("Bifidobacterium", "Bacteroides", "Faecalibacterium", "Blautia", "Alistipes")

genus_cop <- make_target_genus_long(ps = ps_cop, cohort_name = "COPSAC2010", id_col = "child_ID",
                                    time_levels = c("1w", "1m", "1y", "4y", "5y", "6y"),
                                    age_map = age_map_cop, target_genera = target_genera)

genus_ros <- make_target_genus_long(ps = ps_ros, cohort_name = "Roswall et al.", id_col = "child_id",
                                    time_levels = c("0m", "4m", "1y", "3y", "5y"), age_map = age_map_ros,
                                    target_genera = target_genera)

genus_all <- bind_rows(genus_cop, genus_ros)

table(genus_all$Cohort, genus_all$Time_new)
table(genus_all$Genus_clean, genus_all$Cohort)


### E2. Summarize temporal trends for each genus
genus_summary_all <- genus_all %>%
  group_by(Cohort, Time_new, Age_year, Genus_clean) %>%
  summarise( n = n(), median_abun = median(Relative_abundance, na.rm = TRUE),
             q25 = quantile(Relative_abundance, 0.25, na.rm = TRUE),
             q75 = quantile(Relative_abundance, 0.75, na.rm = TRUE),
             mean_abun = mean(Relative_abundance, na.rm = TRUE),
             se_abun = sd(Relative_abundance, na.rm = TRUE) / sqrt(n()), .groups = "drop")

write_xlsx(list(Genus_sample_level = genus_all, Genus_summary = genus_summary_all),
           path = file.path(OUT_COMPARE, "D_target_genus_trajectory_results.xlsx"))


# Perform tests using sample-level data
genus_pairwise <- genus_all %>%
  dplyr::filter(!is.na(Time_new), !is.na(Age_year), !is.na(child_id_for_model),
                is.finite(Relative_abundance)) %>%
  dplyr::group_by(Cohort, Genus_clean) %>%
  dplyr::group_modify(~ {
    
    dat <- .x
    
    # Order time points by age
    time_order <- dat %>%
      dplyr::distinct(Time_new, Age_year) %>%
      dplyr::arrange(Age_year)
    
    dat$Time_new <- factor(as.character(dat$Time_new), levels = as.character(time_order$Time_new))
    
    dat$child_id_for_model <- droplevels(dat$child_id_for_model)
    dat$log_abun <- log10(dat$Relative_abundance + 1e-5)
    
    # Require at least two time points and variation in abundance
    if (nrow(time_order) < 2 || length(unique(dat$log_abun)) < 2) {
      return(tibble::tibble(interval_label = "Not testable", p.value = NA_real_))
    }
    
    model <- lmerTest::lmer(log_abun ~ Time_new + (1 | child_id_for_model), data = dat, REML = TRUE)
    
    emm <- emmeans::emmeans(model, ~ Time_new, lmer.df = "satterthwaite")
    
    result <- as.data.frame(emmeans::contrast(emm, method = "consec", adjust = "none"))
    
    result$interval_label <- paste(head(as.character(time_order$Time_new), -1),
                                   tail(as.character(time_order$Time_new), -1), sep = " → ")
    
    result
  }) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Cohort) %>%
  dplyr::mutate(
    # Apply BH correction across all adjacent contrasts and target genera within each cohort
    p_adj = p.adjust(p.value, method = "BH")) %>%
  dplyr::ungroup()

genus_pairwise %>%
  dplyr::select(Cohort, Genus_clean, interval_label, p_adj)

## Add P values to the existing p_genus_traj
# Prepare labels
genus_p_labels <- genus_pairwise %>%
  dplyr::mutate(p_text = dplyr::case_when(is.na(p_adj) ~ "P = NA", p_adj < 0.001 ~ "P < 0.001",
                                          TRUE ~ paste0("P = ", sprintf("%.3f", p_adj)) ),
                label = paste0(interval_label, ": ", p_text)) %>%
  dplyr::group_by(Cohort, Genus_clean) %>%
  dplyr::summarise(label = paste(c(as.character(first(Cohort)), label), collapse = "\n"),
                   .groups = "drop")


### E3. Plot genus-level developmental trajectories

# Set label heights for each panel according to the genus abundance range
genus_label_height <- genus_summary_all %>%
  dplyr::group_by(Genus_clean) %>%
  dplyr::summarise(ymax = max(c(q75, 0.001), na.rm = TRUE), .groups = "drop")

genus_p_labels <- genus_p_labels %>%
  dplyr::left_join(genus_label_height, by = "Genus_clean") %>%
  dplyr::mutate(label_y = ifelse(Cohort == "COPSAC2010", ymax * 2.5, ymax * 1.75))


# Reference height for each genus panel
genus_label_height <- genus_summary_all %>%
  dplyr::group_by(Genus_clean) %>%
  dplyr::summarise(ymax = max(c(q75, 0.001), na.rm = TRUE), .groups = "drop" )

# Remove existing position columns to allow rerunning
genus_p_labels <- genus_p_labels %>%
  dplyr::select(-dplyr::any_of(c("ymax", "label_x", "label_y"))) %>%
  dplyr::left_join(genus_label_height, by = "Genus_clean") %>%
  dplyr::mutate(label_x = ifelse(Cohort == "COPSAC2010", 0.1, 3.2), label_y = ymax * 1.05 )


# Plot
p_genus_traj <- ggplot2::ggplot(genus_summary_all,
                                ggplot2::aes(x = Age_year, y = median_abun, color = Cohort, group = Cohort)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = q25, ymax = q75), width = 0.08, linewidth = 0.4) +
  ggplot2::facet_wrap(~ Genus_clean, scales = "free_y", nrow = 1) +
  ggplot2::scale_x_continuous(breaks = 0:6, labels = c("0", "1y", "2y", "3y", "4y", "5y", "6y"),
                              limits = c(0, 6)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 strip.background = ggplot2::element_rect(fill = "grey90"),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1) ) +
  ggplot2::scale_colour_manual(values = c("COPSAC2010" = "#0072B2", "Roswall et al." = "#b15928") ) +
  
  ggplot2::geom_text(data = genus_p_labels, ggplot2::aes(x = label_x, y = label_y, label = label,
                                                         colour = Cohort), inherit.aes = FALSE, hjust = 0, vjust = 1,
                     size = 2.5, lineheight = 1, show.legend = FALSE) +
  #ggplot2::labs(caption = paste("Points and bars: median and interquartile range.",
  #                              "P: adjacent-age contrasts from mixed models of log10(abundance + 1e-5).",
  #                              "BH-adjusted across target genera and contrasts within each cohort.",sep = "\n"))+
  
  ggplot2::xlab("Age") +
  ggplot2::ylab("Relative abundance") +
  ggplot2::ggtitle("Genus-level developmental trajectories")+
  ggplot2::theme(legend.position="none")
p_genus_traj

