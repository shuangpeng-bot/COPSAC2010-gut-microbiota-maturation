########################### Gut microbial maturation_Result 5 ################################################

#### Figure 5. Predicting and functional partitioning of gut microbiome driving maturation from age 1 to 4 years ####

# Load required packages
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)
library(writexl)
library(tidyverse)
library(betapart)
library(openxlsx)
library(randomForest)
library(patchwork)
library(ggbreak)
library(pROC)
library(caret)
library(rfUtilities)
library(rfPermute)
library(utils)
library(selbal) # selbal installation:https://github.com/malucalle/selbal
library(Hmisc)
library(readxl)
library(microbiome)
library(SpiecEasi)
library(igraph)
library(zCompositions)
library(compositions)
library(biomformat)
library(Biostrings)
library(readr)
library(ggpicrust2)
library(tibble)
library(ggprism)
library(pheatmap)

# Load phyloseq object
load("phy_feces_clean.RData")
load("phy_species.RData")
phy_feces_clean#
phy_species#


### add maturity (1y&4y) information to sample_data information
phy_feces_clean#
phy_species#

### Add 4y maturity information to sample_data, and carry it to 1y samples (paired design)
df <- data.frame(sample_data(phy_feces_clean))

# Extract 4y maturity per child (one 4y sample per child)
test4y <- data.frame(child_ID = get_variable(df, "child_ID"),
                     Time_new=get_variable(df, "Time_new"),
                     Gut_maturation_4yWeektwo=get_variable(df, "Gut_maturation_Week_two")) %>%
  subset(Time_new=="4y")%>%
  group_by(child_ID) %>%
  filter(n()==1) %>%
  ungroup()
test4y <- test4y[,-2]
dim(test4y);length(unique(test4y$child_ID)) # 4y children (duplicates removed upstream)

# Merge back to full sample metadata by child
df4y <- merge(df,test4y,by="child_ID",all.x = TRUE)

# Create a column that maps 4y maturity onto 1y samples (for downstream pairing)
df4y <- df4y %>%
  mutate(Maturity_4y= ifelse(Time_new=="1y"&!is.na(Gut_maturation_4yWeektwo),Gut_maturation_4yWeektwo, NA))
rownames(df4y)<- df4y$SampleID
df4y$Gut_maturation_4yWeektwo <- NULL
#save(df4y,file="df4y.RData")

# Update sample_data in both phyloseq objects
sample_data(phy_feces_clean)<- df4y
sample_data(phy_species)<- df4y




#####  Figure 5A. Divergence of gut microbiota from age 1 to 4 years, stratified by maturity group #####

# Keep taxa present
phy_all6_species <- prune_taxa(taxa_sums(phy_species) > 0, phy_species) # 2063 taxa, 2689 samples

# Keep 1y samples that have a mapped 4y maturity label
phy_all6_species_1y_4yWeek <- phy_all6_species%>%
  subset_samples(Time_new %in%c("1y")& !is.na(Maturity_4y))# 371 samples
phy_all6_species #

### Build paired 1y/4y dataset, labeling each child by 4y maturity

# Children with known 4y maturity (from 1y subset)
t <- data.frame(sample_data(phy_all6_species_1y_4yWeek)) %>%
  select(child_ID,Maturity_4y)

# Pull 1y and 4y samples for those children
t2 <- as.data.frame(sample_data(phy_species))%>%
  subset( Time_new %in% c("1y","4y")&child_ID %in% t$child_ID)%>% data.frame

# Attach 4y maturity group to both 1y and 4y samples of the same child
t2$Maturity_1y4y <- ifelse(!is.na(match(t2$child_ID, t$child_ID)),t$Maturity_4y[match(t2$child_ID, t$child_ID)],NA)
table(table(t2$child_ID, t2$Maturity_1y4y)==2)#

# Write the label back to sample_data
t3 <- t2 %>% as.data.frame %>% select(SampleID, Maturity_1y4y)
match_index <- match(sample_data(phy_species)$SampleID, t3$SampleID)
sample_data(phy_species)$Maturity_1y4y <- ifelse(!is.na(match_index), t3$Maturity_1y4y[match_index], NA)


# Prepare the paired 1y/4y phyloseq object
phy_1y4y_4y_all6_Week_two <- subset_samples(phy_species,!is.na(Maturity_1y4y)&Time_new%in%c("1y","4y"))
phy_1y4y_4y_all6_Week_two <- prune_taxa(taxa_sums(phy_1y4y_4y_all6_Week_two) > 0, 
                                        phy_1y4y_4y_all6_Week_two)
sample_data(phy_1y4y_4y_all6_Week_two)$Time_new <- as.character(sample_data(phy_1y4y_4y_all6_Week_two)$Time_new )
table(sample_data(phy_1y4y_4y_all6_Week_two)$Maturity_1y4y,sample_data(phy_1y4y_4y_all6_Week_two)$Maturity_4y)

# Relative abundance
ps <- subset_samples(phy_1y4y_4y_all6_Week_two, !is.na(Maturity_1y4y) & Time_new %in% c("1y","4y"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps_rel <- transform_sample_counts(ps, function(x) if (sum(x) == 0) x else x / sum(x))

# Extract OTU/ASV matrix and metadata
otu <- as(otu_table(ps_rel), "matrix")
if (!taxa_are_rows(ps_rel)) otu <- t(otu)

meta <- as(phyloseq::sample_data(ps_rel), "data.frame")
if (!"SampleID" %in% colnames(meta)) {
  meta$SampleID <- rownames(meta)
}
stopifnot(all(c("SampleID", "child_ID", "Time_new", "Maturity_1y4y") %in% colnames(meta)))


# If any child has multiple samples at the same time point, average their profiles (per child-time)
meta_key <- meta %>% transmute(SampleID, child_ID, Time_new, Maturity_1y4y)
dup_tbl <- meta_key %>% count(child_ID, Time_new) %>% filter(n > 1)

if (nrow(dup_tbl) > 0) {
  meta_key <- meta_key %>%
    mutate(merged_id = paste0(child_ID, "_", Time_new))
  split_idx <- split(colnames(otu), meta_key$merged_id[match(colnames(otu), meta_key$SampleID)])
  otu_merged <- sapply(split_idx, function(cols) {
    if (length(cols) == 1) otu[, cols, drop = FALSE] else rowMeans(otu[, cols, drop = FALSE, ])
  })
  if (is.null(dim(otu_merged))) otu_merged <- matrix(otu_merged, nrow = nrow(otu))
  rownames(otu_merged) <- rownames(otu)
  colnames(otu_merged) <- names(split_idx)
  
  meta_merged <- meta_key %>%
    distinct(child_ID, Time_new, Maturity_1y4y, merged_id) %>%
    rename(SampleID = merged_id) %>%
    as.data.frame()
  
  otu <- otu_merged
  meta <- meta_merged
} else {
  stopifnot(identical(colnames(otu), rownames(meta)))
  meta$SampleID <- rownames(meta)
}

# Keep only paired children (have both 1y and 4y)
pair_ids <- names(which(table(meta$child_ID, meta$Time_new)[, "1y"] > 0 &
                          table(meta$child_ID, meta$Time_new)[, "4y"] > 0))
keep <- meta$child_ID %in% pair_ids & meta$Time_new %in% c("1y","4y")
otu <- otu[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]

# Compute within-child Bray–Curtis distance (1y vs 4y)
calc_bc <- function(v1, v2) as.numeric(vegdist(rbind(v1, v2), method = "bray")[1])

bc_df <- meta %>%
  dplyr::select(SampleID, child_ID, Time_new, Maturity_1y4y) %>%
  dplyr::group_by(child_ID) %>%
  dplyr::summarize(
    s1y = SampleID[match("1y", Time_new)],
    s4y = SampleID[match("4y", Time_new)],
    Maturity_1y4y = first(Maturity_1y4y),
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(s1y), !is.na(s4y)) %>%
  dplyr::mutate(BC_1y_to_4y = mapply(function(a, b) calc_bc(otu[, a], otu[, b]), s1y, s4y)) %>%
  dplyr::arrange(Maturity_1y4y)

# Plot: within-subject beta change (1y → 4y) by 4y maturity group
col_maturity <- c("High"="#E69F00","Low"="#56B4E9")
bray_dissimilarity_1y4y_p2 <- ggplot(bc_df, aes(x = Maturity_1y4y, y = BC_1y_to_4y)) +
  geom_jitter(aes(fill=Maturity_1y4y),width = 0.3, alpha = 0.55, size = 1.8) +
  geom_violin(aes(fill=Maturity_1y4y),alpha=0.8,outlier.shape = NA) +
  geom_boxplot(aes(fill=Maturity_1y4y),width=0.1,alpha=0.1)+
  xlab("Maturity at 4y")+
  ylab("Distance from 1y to 4y")+
  theme_bw(base_size = 14)+
  stat_compare_means()+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)
print(bray_dissimilarity_1y4y_p2)



  
#####  Figure 5B. Relative contribution of nestedness to overall beta diversity in each maturity group #####

### Beta diversity partitioning/ecological partitioning of beta diversity
phy_all6_species#

# Prepare data
phy_species3  <- phy_all6_species %>% subset_samples(Time_new %in%c("1y","4y"))
phy_species3  <- phy_species3%>%
  subset_samples(!is.na(Maturity_4y)|!is.na(Maturity_1y4y))# 371 *2 samples
phy_species3 <- prune_taxa(taxa_sums(phy_species3) > 0, phy_species3) # 1403 taxa

# Create a unified maturity variable
sample_data(phy_species3)$Maturity_1y4y_2 <- ifelse(
  !is.na(sample_data(phy_species3)$Maturity_4y),sample_data(phy_species3)$Maturity_4y,
  ifelse(!is.na(sample_data(phy_species3)$Maturity_1y4y),sample_data(phy_species3)$Maturity_1y4y,NA))

# Define maturity groups
groups <- c("High", "Low")
results_list <- list()

for (gp in groups) {
  
  cat("Processing Group:", gp, "\n")
  
  # Subset samples by group (High / Low) and timepoint (1y, 4y)
  phy_sub <- subset_samples(phy_species3,Maturity_1y4y_2 == gp & Time_new %in% c("1y","4y"))
  phy_sub <- prune_taxa(taxa_sums(phy_sub) > 0, phy_sub)
  
  # Skip if the group does not contain both 1y and 4y
  meta_sub <- as.data.frame(sample_data(phy_sub))
  if (!all(c("1y","4y") %in% meta_sub$Time_new)) {
    warning("Group ", gp, " does not have both 1y and 4y; skipping.")
    next
  }
  
  # Extract OTU table (samples × taxa)
  otu <- as(otu_table(phy_sub), "matrix")
  if (taxa_are_rows(phy_sub)) otu <- t(otu)
  
  # Convert to presence/absence
  otu_pa <- (otu > 0) * 1
  
  # Collapse samples by timepoint (sum → presence/absence)
  idx_1y <- rownames(otu_pa)[meta_sub$Time_new == "1y"]
  idx_4y <- rownames(otu_pa)[meta_sub$Time_new == "4y"]
  
  pa_1y <- as.numeric(colSums(otu_pa[idx_1y, , drop = FALSE]) > 0)
  pa_4y <- as.numeric(colSums(otu_pa[idx_4y, , drop = FALSE]) > 0)
  
  # Combine into a 2 × species matrix (rows: 1y, 4y)
  X_pa <- rbind(`1y` = pa_1y, `4y` = pa_4y)
  
  # Partition Jaccard dissimilarity into turnover and nestedness components
  bp <- betapart::beta.pair(X_pa, index.family = "jaccard")
  Turnover   <- as.numeric(bp$beta.jtu) # Bacterial turnover
  Nestedness <- as.numeric(bp$beta.jne) # Bacterial nesedness
  Total      <- as.numeric(bp$beta.jac) # Total
  
  results_list[[gp]] <- data.frame(
    Group = gp,
    Time  = "1y_vs_4y",
    Turnover = Turnover,
    Nestedness = Nestedness,
    Total = Total,
    stringsAsFactors = FALSE
  )
}

# Combine results and compute nestedness contribution (%)
beta_all_df <- bind_rows(results_list) %>%
  mutate(Nestedness_pct = Nestedness / (Turnover + Nestedness) * 100)
beta_all_df
# Nestedness_pct represents the percentage contribution of nestedness to the overall community dissimilarity.

# Visualization
beta_all_df2 <- melt(beta_all_df)
col_maturity <- c("High"="#E69F00","Low"="#56B4E9")

p_Nestedness_pct_jaccard_partitioning <- beta_all_df2 %>%
  subset(variable %in%c("Nestedness_pct"))%>%
  ggplot(aes(x=Group, y=value, fill=Group))+
  geom_bar(stat = "identity", position = "dodge", alpha=0.75) +
  geom_text(aes(label = paste0(round(value, 1), "%")),
            position = position_dodge(width = 0.9),
            vjust = 1.5, size = 5) +
  scale_fill_manual("Maturity_4y",values = col_maturity) +
  labs(y = "Nestedness Contribution (%)",
       x = "1y vs. 4y",
       title = "Relative Contribution of Nestedness to Beta Diversity") +
  theme_minimal(base_size = 14)
p_Nestedness_pct_jaccard_partitioning




#####  Figure S11. Ten-fold cross-validation error as a function of the number of input species used to differentiate between high- and low- gut microbial maturation in order of variable importance #####
phy_species #

###  used to differentiate high vs low gut microbial maturation (ordered by variable importance)

###### a. Prepare data ######
phy_all6_species <- prune_taxa(taxa_sums(phy_species) > 0, phy_species)

# Use 1y samples that have 4y maturity mapping
phy_all6_species_1y_4yWeek <- phy_all6_species%>%
  subset_samples(Time_new %in%c("1y")& !is.na(Maturity_1y4y))

phy_all6_species_1y_4yWeek <- prune_taxa(taxa_sums(phy_all6_species_1y_4yWeek) > 0, 
                                         phy_all6_species_1y_4yWeek)
phy_all6_species_1y_4yWeek

# Sample metadata
all6_species_1y_4yWeek_sample_df <- data.frame(
  SampleID=get_variable(phy_all6_species_1y_4yWeek,"SampleID"),
  Maturity_1y4y=get_variable(phy_all6_species_1y_4yWeek,"Maturity_1y4y")%>%as.character,
  Microbiota_age=get_variable(phy_all6_species_1y_4yWeek,"Microbiota_age"),
  AGE=get_variable(phy_all6_species_1y_4yWeek,"AGE"),
  Time_new=get_variable(phy_all6_species_1y_4yWeek,"Time_new"))
rownames(all6_species_1y_4yWeek_sample_df) <- all6_species_1y_4yWeek_sample_df$SampleID

# Normalize OTU table to relative abundance (rows = samples)
all6_species_1y_4yWeek_otu_df <- transform_sample_counts(phy_all6_species_1y_4yWeek,
                                                         function(x) x/sum(x))%>%
  otu_table%>%
  as("matrix")%>%
  t%>%
  as.data.frame

# Make column names unique/safe
colnames(all6_species_1y_4yWeek_otu_df) <- paste0("x",colnames(all6_species_1y_4yWeek_otu_df))

# Construct RF input: OTU matrix + maturity label
all6_species_1y_4yWeek <- cbind(
  all6_species_1y_4yWeek_otu_df,
  Maturity_1y=all6_species_1y_4yWeek_sample_df[
    match(rownames(all6_species_1y_4yWeek_otu_df),
          rownames(all6_species_1y_4yWeek_sample_df)),
    "Maturity_1y4y"])%>%
  as.data.frame()

# Ensure grouping is a factor for classification
all6_species_1y_4yWeek$Maturity_1y4y <- factor(all6_species_1y_4yWeek$Maturity_1y4y)
# If the outcome were numeric, that would be regression; here it's classification.

# Train/test split
table(sample_data(phy_all6_species_1y_4yWeek)$Maturity_1y4y) # High 173, Low 198

set.seed(77777)
split_species_1y_4yWeek <- sample(2,nrow(all6_species_1y_4yWeek),replace = TRUE,prob = c(0.7,0.3))
train_all6_species_1y_4yWeek <- all6_species_1y_4yWeek[split_species_1y_4yWeek==1,]
test_all6_species_1y_4yWeek <- all6_species_1y_4yWeek[split_species_1y_4yWeek==2,]

table(train_all6_species_1y_4yWeek$Maturity_1y4y)# 122 High, Low 132
table(test_all6_species_1y_4yWeek$Maturity_1y4y)# 51 High, Low 66

#save(train_all6_species_1y_4yWeek, file="train_all6_species_1y_4yWeek.RData")
#save(test_all6_species_1y_4yWeek, file="test_all6_species_1y_4yWeek.RData")



###### b. Random Forest cross-validation for feature selection ######
train_all6_species_1y_4yWeek #
dim(train_all6_species_1y_4yWeek) # 254 * 790
anyNA(train_all6_species_1y_4yWeek)

# rfcv: 10-fold CV; step = -5 removes features in steps of 5 (scale=FALSE is required with negative step) 
set.seed(12345);
rfcv_r10_all6_species_1y_4yWeek_10taxa <- replicate(10,rfcv(trainx = train_all6_species_1y_4yWeek[,-790],
                                                            trainy = train_all6_species_1y_4yWeek[,790],
                                                            cv.fold = 10,scale=FALSE,step = -10),
                                                    simplify = FALSE)
#save(rfcv_r10_all6_species_1y_4yWeek_10taxa,file="rfcv_r10_all6_species_1y_4yWeek_10taxa.RData")
#load("rfcv_r10_all6_species_1y_4yWeek_10taxa.RData")

# Summarize CV error across replicates and compute the mean
rfcv_r10_all6_species_1y_4yWeek_10taxa.cv <- data.frame(
  sapply(rfcv_r10_all6_species_1y_4yWeek_10taxa,"[[","error.cv"),
  mean=rowMeans(sapply(rfcv_r10_all6_species_1y_4yWeek_10taxa, "[[","error.cv")))

# Row with the minimum mean CV error
rfcv_r10_all6_species_1y_4yWeek_10taxa.cv[which.min(rfcv_r10_all6_species_1y_4yWeek_10taxa.cv[,11]),]

# Order by mean CV error and compute SD/SE for error bars
rfcv_10taxa_mean <- rfcv_r10_all6_species_1y_4yWeek_10taxa.cv[order(rfcv_r10_all6_species_1y_4yWeek_10taxa.cv$mean),]
rfcv_10taxa_mean$sd <- apply(rfcv_10taxa_mean[,1:10],1,sd,na.rm=T) # 1 means by row, 2 means by column
rfcv_10taxa_mean$se <- rfcv_10taxa_mean$sd/sqrt(10)
rfcv_10taxa_mean$taxa_No <- rownames(rfcv_10taxa_mean) %>% as.numeric() # a simple index for x-axis

# Plot method: ggplot with mean ± SE (index on x)
rfcv_r10_all6_species_1y_4yWeek_10taxa <- ggplot(rfcv_10taxa_mean, aes(y=mean, x= taxa_No))+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=2 )+
  scale_x_continuous(breaks = seq(9,799,by=18))+
  ggbreak::scale_x_break(c(110,174,204,765))+ # optional visual break
  ylim(0.39,0.47)+
  geom_line(color="blue")+
  geom_point(color="black",size=2)+
  xlab("Number of identity species")+
  ylab("Increase in mean squared error")+
  ggtitle("Cross-validation error")+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_vline(xintercept = 99, linetype="dashed", color = "black", size=0.2)+
  geom_hline(yintercept=0.3956693, linetype="dashed", color="black", size=0.2)+
  
  theme( strip.background = element_blank(),
         text=element_text(size=11))
rfcv_r10_all6_species_1y_4yWeek_10taxa + theme(plot.title = element_text(hjust = 0.5))




#####  Figure 5C. ROC curve and AUC values for test samples, using the top 99 bacterial species #####
phy_species #
train_all6_species_1y_4yWeek #
test_all6_species_1y_4yWeek #


###### a. Random Forest model ######
train_all6_species_1y_4yWeek$Maturity_4y<- factor(train_all6_species_1y_4yWeek$Maturity_4y)
test_all6_species_1y_4yWeek$Maturity_4y<- factor(test_all6_species_1y_4yWeek$Maturity_4y)

set.seed(333);
Model1 <- randomForest(
  Maturity_4y~.,
  data=train_all6_species_1y_4yWeek,
  xtest=train_all6_species_1y_4yWeek[,-which(names(train_all6_species_1y_4yWeek)=="Maturity_4y")],
  ytest=train_all6_species_1y_4yWeek[,"Maturity_4y"],
  ntree=10000,keep.forest=TRUE,#4200
  importance=TRUE,
  proximity=TRUE)

Model1 # OOB 40.16%
varImpPlot(Model1)
plot(Model1)
legend(x=4500,y=0.51, col=1:5,cex=0.85,fill=1:5,colnames(Model1$err.rate))
rfPermute::confusionMatrix(Model1)

# Classification performance (Sensitivity and specificity)
pred_train <- predict(Model1, train_all6_species_1y_4yWeek, type = "class")
caret::confusionMatrix(pred_train ,train_all6_species_1y_4yWeek [,"Maturity_4y"])# Kappa : 1
pred_test <- predict(Model1, test_all6_species_1y_4yWeek, type = "class")
caret::confusionMatrix(pred_test,test_all6_species_1y_4yWeek [,"Maturity_4y"])#  Kappa : 0.1933, Accuracy : 0.6068

# ROC plot (test-data)
ROC_Model1_testdatda <- pROC::plot.roc(test_all6_species_1y_4yWeek$Maturity_4y,
                                       predict(Model1, newdata = test_all6_species_1y_4yWeek, type="prob")[,2],
                                       ylim=c(0,1),xlim=c(1,0),
                                       smooth=FALSE,ci=TRUE, col="red",
                                       lwd=2,
                                       legacy.axes=T)
ci.lower <- round(as.numeric(ROC_Model1_testdatda$ci[1]),3) 
ci.upper <- round(as.numeric(ROC_Model1_testdatda$ci[3]),3) 
legend.name <- paste("AUC:", round(as.numeric(ROC_Model1_testdatda$auc),3), paste(ci.lower,ci.upper,sep="-"))
legend("bottomright",
       legend=legend.name,
       col="red",
       lwd=2,
       bty="n")



###### b. Select top 99 important taxa by variable importance ######
Model1#
train_all6_species_1y_4yWeek#
test_all6_species_1y_4yWeek#

# Make a data frame with predictor names and their importance
imp_df1  <- data.frame(importance(Model1,scale=TRUE),check.names = FALSE)
imp_df1 <- data.frame(predictors=rownames(imp_df1),imp_df1)

# Order the predictor levels by importance
imp_df2<- arrange(imp_df1,desc(MeanDecreaseGini))
imp_df2$predictors <- factor(imp_df2$predictors,levels=imp_df2$predictors)

# Map taxa names back to species
imp_df3 <- imp_df2[1:99,]
taxa_names_99_1 <- imp_df3$predictors
imp_df <- substring(taxa_names_99_1, first = 2)
phy_all6_species_1y_99 <- prune_taxa(taxa_names(phy_all6_species_1y_4yWeek) %in% imp_df,phy_all6_species_1y_4yWeek)
taxa_names_99_speciesname <- as.data.frame(tax_table(phy_all6_species_1y_99)[,7])

taxa_names_99_1#
phy_all6_species_1y_99#
taxa_names_99_speciesname #
write_xlsx(taxa_names_99_speciesname,  "taxa_names_99_speciesname.xlsx")

# Subset training and testing data to top 99 taxa
train_all6_species_1y_4yWeek_99taxa <- dplyr::select(train_all6_species_1y_4yWeek,"Maturity_4y",all_of(taxa_names_99_1))
test_all6_species_1y_4yWeek_99taxa <- dplyr::select(test_all6_species_1y_4yWeek,"Maturity_4y",all_of(taxa_names_99_1))




###### c. Hyperparameter tuning ######
#### Hyperparameter tuning: explore ntree=8600




###### d. Random Forest model (Model2: top 99 taxa) ######
train_all6_species_1y_4yWeek_99taxa #
test_all6_species_1y_4yWeek_99taxa #

### Model2
set.seed(333);
Model2 <- randomForest(
  Maturity_4y~.,
  data=train_all6_species_1y_4yWeek_99taxa,
  xtest=train_all6_species_1y_4yWeek_99taxa[,-which(names(train_all6_species_1y_4yWeek_99taxa)=="Maturity_4y")],
  ytest=train_all6_species_1y_4yWeek_99taxa[,"Maturity_4y"],
  ntree=8600,keep.forest=TRUE,
  importance=TRUE,
  proximity=TRUE)

Model2# OOB 38.98%
varImpPlot(Model2)
plot(Model2)
legend(x=350,y=0.52, col=1:5,cex=0.85,fill=1:5,colnames(Model2$err.rate))

# Overall confusion matrix from rfPermute (model-level summary)
rfPermute::confusionMatrix(Model2) 

# Confusion matrices (Sensitivity and specificity)
pred_train <- predict(Model2, train_all6_species_1y_4yWeek_99taxa, type = "class")
caret::confusionMatrix(pred_train,train_all6_species_1y_4yWeek_99taxa [,"Maturity_4y"]) # Kappa : 1
pred_test <- predict(Model2, test_all6_species_1y_4yWeek, type = "class")
caret::confusionMatrix(pred_test,test_all6_species_1y_4yWeek [,"Maturity_4y"])#  Kappa : 0.2284, Accuracy : 0.6239

### Variable importance (top 20)
Model2 #
Model2$importance #
taxa_names_99_speciesname #

# Make sure species name table has a predictor key like "x<taxon>"
taxa_names_99_speciesname$predictors <- paste0("x",rownames(taxa_names_99_speciesname))
"[Clostridium] innocuum group sp."%in% taxa_names_99_speciesname$species

# Plot top-importance variables
MDG_MDA_predict_Model2 <- varImpPlot(Model2,n.var=min(20, nrow(x$importance)))
MDG_MDA_predict_Model2
# MeanDecreaseAccuracy: drop in RF accuracy when the variable is permuted (larger = more important)
# MeanDecreaseGini: total decrease in node impurity contributed by the variable (larger = more important)


### ROC plot of test data in Model2
# ROC plot
ROC_Model2 <- pROC::plot.roc(test_all6_species_1y_4yWeek_99taxa$Maturity_4y,
                             predict(Model2, newdata = test_all6_species_1y_4yWeek_99taxa, type="prob")[,2],
                             ylim=c(0,1),xlim=c(1,0),
                             smooth=FALSE,ci=TRUE, col="red",
                             lwd=2,
                             print.auc=T,
                             legacy.axes=T)
ci.lower <- round(as.numeric(ROC_Model2_testdatda$ci[1]),3) 
ci.upper <- round(as.numeric(ROC_Model2_testdatda$ci[3]),3) 
legend.name <- paste("AUC:", round(as.numeric(ROC_Model2_testdatda$auc),3), paste(ci.lower,ci.upper,sep="-"))
legend("bottomright",
       legend=legend.name,
       col="red",
       lwd=2,
       bty="n")

# Overlay Model1 and Model2
pROC::plot.roc(test_all6_species_1y_4yWeek$Maturity_4y,
               predict(Model1, newdata = test_all6_species_1y_4yWeek, type="prob")[,2],
               ylim=c(0,1),xlim=c(1,0),
               smooth=FALSE,ci=TRUE, col="blue",
               lwd=2,
               print.auc=T,
               legacy.axes=T)

pROC::plot.roc(test_all6_species_1y_4yWeek$Maturity_4y,
               predict(Model2, newdata = test_all6_species_1y_4yWeek, type="prob")[,2],
               ylim=c(0,1),xlim=c(1,0),
               smooth=FALSE,ci=TRUE, col="red",
               lwd=2,
               #lty = 2,
               print.auc=T,add=T,print.auc.adj = c(0,3),
               legacy.axes=T)
legend('bottom',legend = c('Model1','Model2'),col=c('blue','red'),lwd=2)




#####  Figure 5D. Top 20 bacterial species in the 1-year gut microbiota identified by random forest classification #####
Model2 #
train_all6_species_1y_4yWeek_99taxa #
test_all6_species_1y_4yWeek_99taxa #


###### a. Top 20 bacterial species in the 1-year gut microbiota identified by random forest classification #####
## top 20 ~ significant
# Make a data frame with predictor names and their importance
test1  <- data.frame(importance(Model2,scale=TRUE),check.names = FALSE)
test1 <- data.frame(predictors=rownames(test1),test1)

# merge the species name information
test2 <- merge(test1, taxa_names_99_speciesname, by="predictors",all.x = T)

# Order the predictor levels by importance and select top 99
test3 <- arrange(test2,desc(MeanDecreaseAccuracy))
test4 <- test3[1:99,]

# How many of the top features overlap with prior significant species lists
table(test4$species %in% Sig_species_1y)
table(test4$species %in% RF_Sig_60taxa$.)
test4$species
test4$species[test4$species %in% RF_Sig_60taxa$.==T]

RF_importance_top20 <- test4$species

## Phylum annotations for the selected species
tax_df <- as.data.frame(tax_table(phy_species), stringsAsFactors = FALSE)
colnames(tax_df) <- tolower(colnames(tax_df))
sp <- unique(na.omit(test4$species))
test5 <- merge(data.frame(species = sp), tax_df, by = "species", all.x = TRUE, sort = FALSE)

# If another phyloseq object 'phy_species2' is used to filter species, keep the original order
species_keep <- unique(na.omit(as.data.frame(tax_table(phy_species2))$species))
test5 <- test5 %>%subset(test4$species %in% species_keep)
test5$species


###### b. Permutation test (significance of classification) #####
train_all6_species_1y_4yWeek_99taxa#
test_all6_species_1y_4yWeek_99taxa#
Model2 #
Model2$importance #
taxa_names_99_speciesname #

# Split labels and features; align factor levels
y_train <- train_all6_species_1y_4yWeek_99taxa[, 1]
X_train <- train_all6_species_1y_4yWeek_99taxa[, 2:100]
y_test  <- test_all6_species_1y_4yWeek_99taxa[, 1]
X_test  <- test_all6_species_1y_4yWeek_99taxa[, 2:100]

# Ensure classification labels are factors and levels match between train/test
y_train <- factor(y_train)                       
y_test  <- factor(y_test, levels = levels(y_train))  

# Observed accuracy on the held-out test set
set.seed(333)
rf_real <- randomForest(
  x = X_train, y = y_train,
  xtest = X_test, ytest = y_test,
  ntree = 8600, importance = TRUE, keep.forest = FALSE
)
acc_real <- mean(rf_real$test$predicted == y_test)
acc_real# 0.6153846 with ntree =8600


## Permutation test (right-tailed): shuffle training labels only
# Ho: features and labels are independent
n_perm <- 1000
acc_perm <- numeric(n_perm)
set.seed(333)

for (i in 1:n_perm) {
  # Shuffle training labels; keep features and the test set unchanged
  y_perm <- sample(y_train)                                   
  y_perm <- factor(y_perm, levels = levels(y_train))          
  
  rf_perm <- randomForest(
    x = X_train, y = y_perm,
    xtest = X_test, ytest = y_test,
    ntree = 8600, keep.forest = FALSE                         
  )
  
  # Test accuracy under the null (permuted labels)
  acc_perm[i] <- mean(rf_perm$test$predicted == y_test)
 }

# Right-tailed permutation p-value: P(null accuracy >= observed accuracy)
p_perm <- mean(acc_perm >= acc_real)
cat(sprintf("Observed accuracy = %.4f; permutation p = %.4g\n", acc_real, p_perm))
# Observed accuracy = 0.6154; permutation p = 0.033
# Result note: With 1000 permutations, observed accuracy = 0.6154; permutation p = 0.033




#####  Figure S12. Microbiome compositional balance discriminated gut microbial maturity groups by age four #####
phy_feces_clean #
train_all6_species_1y_4yWeek #
test_all6_species_1y_4yWeek #

###### a. reparation (align columns & factors levels) ######

set.seed(333)

# Keep only ASV columns shared by train and test (example: starting with "x")
asv_cols <- intersect(
  grep("^x", names(train_all6_species_1y_4yWeek), value = TRUE),
  grep("^x", names(test_all6_species_1y_4yWeek),  value = TRUE))

# Align factor levels: low/High (treat "High" as the positive class for AUC)
train <- train_all6_species_1y_4yWeek %>%
  dplyr::select(all_of(asv_cols), Maturity_4y) %>%
  dplyr::mutate(Maturity_4y = factor(Maturity_4y, levels = c("Low", "High")))

test  <- test_all6_species_1y_4yWeek %>%
  dplyr::select(all_of(asv_cols), Maturity_4y) %>%
  dplyr::mutate(Maturity_4y = factor(Maturity_4y, levels = c("Low", "High")))

X_tr <- as.matrix(train[, asv_cols, drop = FALSE])
y_tr <- droplevels(train$Maturity_4y)
sum(!is.finite(X_tr)); sum(!is.finite(X_te)) # Check NA/Inf (should be 0)
# Check for NA/Inf in the training matrix (both counts should be 0)

X_te <- as.matrix(test[,  asv_cols, drop = FALSE])
y_te <- droplevels(test$Maturity_4y)


###### b. Utility functions (zero-replacement, stratified K-fold, selbal wrappers, balance computation) ######

# Row-wise half-min zero replacement (if a row is all zeros, fill with tiny constant)
zrep_halfmin <- function(M){
  M <- as.matrix(M)
  rows <- lapply(seq_len(nrow(M)), function(i){
    r  <- M[i, ]
    # make NA as 0 to deal with
    is_pos <- is.finite(r) & r > 0
    nz <- r[is_pos]
    if (length(nz) == 0L) {
      rep(1e-12, length(r))
    } else {
      halfmin <- 0.5 * min(nz)
      r[!is.finite(r) | r == 0] <- halfmin
      r
    }
  })
  M2 <- do.call(rbind, rows)
  dimnames(M2) <- dimnames(M)
  storage.mode(M2) <- "double"
  M2
}

# Simple stratified K-fold splitter
make_folds <- function(y, K = 5, seed = 333){
  set.seed(seed)
  y <- factor(y)
  idx_by_class <- split(seq_along(y), y)
  fold_id <- integer(length(y))
  for (cl in names(idx_by_class)){
    ids <- sample(idx_by_class[[cl]])
    fold_id[ids] <- rep(seq_len(K), length.out = length(ids))
  }
  fold_id
}

# Extract numerator/denominator taxa from a selbal object (version-robust)
get_pos_neg <- function(sb_obj){
  # 
  cand_pos <- c("pos","numerator","Numerator","Numerador","positive","Positivo")
  cand_neg <- c("neg","denominator","Denominator","Denominador","negative","Negativo")
  dig <- function(x, keys){
    # 
    if (is.list(x)) {
      nx <- names(x)
      if (!is.null(nx)) {
        hit <- which(tolower(nx) %in% tolower(keys))
        for (h in hit) if (is.character(x[[h]]) && length(x[[h]])>0) return(as.character(x[[h]]))
      }
      for (i in seq_along(x)) {
        res <- dig(x[[i]], keys)
        if (!is.null(res)) return(res)
      }
    }
    NULL
  }
  # 
  pick_path <- function(obj, nm1, nm2){
    if (!is.null(obj[[nm1]]) && is.list(obj[[nm1]])) {
      for (k in cand_pos) if (!is.null(obj[[nm1]][[k]])) pos <<- as.character(obj[[nm1]][[k]])
      for (k in cand_neg) if (!is.null(obj[[nm1]][[k]])) neg <<- as.character(obj[[nm1]][[k]])
    }
  }
  pos <- neg <- NULL
  pick_path(sb_obj, "global.balance", NULL)
  if (is.null(pos)) pos <- dig(sb_obj, cand_pos)
  if (is.null(neg)) neg <- dig(sb_obj, cand_neg)
  
  # Fallback: parse print output
  if (is.null(pos) || is.null(neg) || length(pos)==0 || length(neg)==0){
    cap <- utils::capture.output(print(sb_obj))
    # Find Numerator/Denominator rows
    lineN <- grep("(?i)\\b(numerator|numerador|pos(itive)?)\\b", cap, perl=TRUE, value=TRUE)
    lineD <- grep("(?i)\\b(denominator|denominador|neg(ative)?)\\b", cap, perl=TRUE, value=TRUE)
    parse_vec <- function(s){
      s <- gsub(".*?:", "", s)            
      s <- gsub("[\\|/\\+]", ",", s)     
      x <- unique(trimws(unlist(strsplit(s, "\\s*,\\s*"))))
      x[nchar(x) > 0]
    }
    if (length(lineN)) pos <- parse_vec(paste(lineN, collapse=","))
    if (length(lineD)) neg <- parse_vec(paste(lineD, collapse=","))
    # 
    if ((is.null(pos) || !length(pos)) || (is.null(neg) || !length(neg))) {
      expr_line <- cap[grep("\\(", cap)]
      if (length(expr_line)) {
        s <- paste(expr_line, collapse=" ")
        m <- regmatches(s, gregexpr("\\(([^\\)]*)\\)", s))
        if (length(m) && length(m[[1]]) >= 2) {
          pos <- unique(trimws(unlist(strsplit(gsub("[()]", "", m[[1]][1]), "\\+"))))
          neg <- unique(trimws(unlist(strsplit(gsub("[()]", "", m[[1]][2]), "\\+"))))
          pos <- pos[nchar(pos)>0]; neg <- neg[nchar(neg)>0]
        }
      }
    }
  }
  
  if (is.null(pos) || is.null(neg) || length(pos)==0 || length(neg)==0) {
    stop("Unable to identify the numerator/denominator from the selbal results despite attempts using version-compatible extraction and text parsing.")
  }
  list(pos = as.character(pos), neg = as.character(neg))
}

# selbal() wrapper: apply half-min zero replacement before fitting
fit_selbal_once <- function(X, y, maxV = 20){
  X2 <- zrep_halfmin(X)
  X2 <- as.data.frame(X2, check.names = FALSE)  
  # 
  out <- tryCatch(selbal::selbal(x = X2, y = y, maxV = maxV),
                  error = function(e) e)
  if (inherits(out, "error")) {
    # 
    out2 <- tryCatch(selbal::selbal(x = X2, y = y, maxV = maxV, zero.rep = "bayes"),
                     error = function(e) e)
    if (inherits(out2, "error")) {
      stop(paste0("selbal() failed：", conditionMessage(out2)))
      #。
    } else {
      return(get_pos_neg(out2))
    }
  } else {
    return(get_pos_neg(out))
  }
}

# Compute balance score given numerator/denominator taxa
compute_balance <- function(X, pos, neg){
  X <- as.matrix(X)
  # Only keep the bacterial taxa name in current matrix
  pos2 <- intersect(pos, colnames(X))
  neg2 <- intersect(neg, colnames(X))
  # 
  if (length(pos2) < 1L || length(neg2) < 1L) {
    return(rep(NA_real_, nrow(X)))
  }
  Xsub <- X[, unique(c(pos2, neg2)), drop = FALSE]
  Xsub <- zrep_halfmin(Xsub)
  Xsub <- as.data.frame(Xsub, check.names = FALSE)  
  
  # 
  lp <- rowMeans(log(Xsub[, pos2, drop = FALSE]))
  ln <- rowMeans(log(Xsub[, neg2, drop = FALSE]))
  as.numeric(lp - ln)
}


###### c.  CV for hyperparameter selection + multi-balance construction ######
cv_multibalance <- function(X, y, K = 5,
                            grid_maxV = c(5, 10, 20, 30),#maxV：maximum bacterial in one selbal balance
                            grid_nB = 1:3,#1:3 how many selbal balances in order
                            seed = 333){
  y <- factor(y, levels = c("Low","High"))#fix the factor level (Low/High)
  folds <- make_folds(y, K = K, seed = seed)
  res <- expand.grid(maxV = grid_maxV, nB = grid_nB)# create the matrix（maxV × nB）
  res$AUC <- NA_real_#record the AUC。
  
  ## —— LOG: initialization
  cv_log <<- NULL
  .cv_log_rows <- list()
  .log_add <- function(...){
    .cv_log_rows[[length(.cv_log_rows)+1]] <<- data.frame(..., check.names = FALSE, stringsAsFactors = FALSE)
  }
  #
  on.exit({
    cv_log <<- tryCatch(
      if (length(.cv_log_rows)) dplyr::bind_rows(.cv_log_rows) else data.frame(),
      error = function(e) { message("log merge failed: ", conditionMessage(e)); data.frame() }
    )
  }, add = TRUE)
  
  for (gi in seq_len(nrow(res))){
    maxV <- res$maxV[gi]; nB <- res$nB[gi]
    aucs <- c()
    # aucs to collect AUC in each K fold
    
    for (fold in seq_len(K)){
      tr_idx <- which(folds != fold)
      va_idx <- which(folds == fold)
      Xtr <- X[tr_idx, , drop = FALSE]
      ytr <- y[tr_idx]
      Xva <- X[va_idx, , drop = FALSE]
      yva <- y[va_idx]
      # the fold fold：devide the train/test dataset
      
      if (nlevels(droplevels(ytr)) < 2L) {
        .log_add(step = "skip_fold_one_class_tr", ok = FALSE,
                 maxV = maxV, nB = nB, fold = fold)
        next
      }
      if (nlevels(droplevels(yva)) < 2L) {
        .log_add(step = "skip_fold_one_class_va", ok = FALSE,
                 maxV = maxV, nB = nB, fold = fold)
        next
      }
      # 
      
      #
      .log_add(step = "start_fold", ok = TRUE,
               maxV = maxV, nB = nB, fold = fold,
               n_tr = length(tr_idx), n_va = length(va_idx))
      
      # 
      remain <- colnames(Xtr)# remain the avaliable ASV table
      bal_tr_list <- list(); bal_va_list <- list()
      
      for (b in seq_len(nB)){
        Xtr_b <- Xtr[, remain, drop = FALSE]
        # filter the rare species
        keep <- colMeans(Xtr_b > 0) >= 0.05
        ## —— LOG: filtering results —— ##
        .log_add(step = paste0("b", b, "_filter"), 
                 ok = as.logical(sum(keep) >= 2),
                 maxV = maxV, nB = nB, fold = fold,
                 n_before = ncol(Xtr_b), n_keep = sum(keep))
        
        if (!any(keep)) break
        Xtr_b <- Xtr_b[, keep, drop = FALSE]
        if (ncol(Xtr_b) < 2) break
        # the b balance：subset the matrix form remaining data
        sb <- tryCatch(
          fit_selbal_once(Xtr_b, ytr, maxV = maxV),
          error = function(e){
            ## —— LOG: selbal failed —— ##
            .log_add(step = paste0("b", b, "_selbal"), ok = FALSE,
                     maxV = maxV, nB = nB, fold = fold,
                     err = substr(conditionMessage(e), 1, 200))
            NULL
          }
        )
        if (is.null(sb)) break
        pos <- sb$pos; neg <- sb$neg
        
        # the balance score on training/test data
        tr_bal <- compute_balance(Xtr[, unique(c(pos,neg)), drop = FALSE], pos, neg)
        va_bal <- compute_balance(Xva[, unique(c(pos,neg)), drop = FALSE], pos, neg)
        bal_tr_list[[length(bal_tr_list)+1]] <- tr_bal
        bal_va_list[[length(bal_va_list)+1]] <- va_bal
        ## —— LOG: balance calculation —— ##
        .log_add(step = paste0("b", b, "_balance"),
                 ok = !(all(is.na(tr_bal)) || all(is.na(va_bal))),
                 maxV = maxV, nB = nB, fold = fold,
                 pos_len = length(pos), neg_len = length(neg),
                 anyNA_tr = any(!is.finite(tr_bal)),
                 anyNA_va = any(!is.finite(va_bal)))
        
        # 
        remain <- setdiff(remain, unique(c(pos,neg)))
        if (length(remain) < 2) break
      }
      
      if (length(bal_tr_list) == 0L) {
        ## —— LOG: no useful balance in this part —— ##
        .log_add(step = "no_balance", ok = FALSE,
                 maxV = maxV, nB = nB, fold = fold)
        next  
      }
      
      bal_tr <- as.data.frame(bal_tr_list, optional = TRUE, check.names = FALSE)
      bal_va <- as.data.frame(bal_va_list, optional = TRUE, check.names = FALSE)
      colnames(bal_tr) <- paste0("B", seq_len(ncol(bal_tr)))
      colnames(bal_va) <- paste0("B", seq_len(ncol(bal_va)))
      
      # Drop constant/all-NA balance columns using training fold
      keep_num <- vapply(bal_tr, function(x) {
        any(is.finite(x)) && sd(x, na.rm = TRUE) > 0
      }, logical(1))
      bal_tr <- bal_tr[, keep_num, drop = FALSE]
      bal_va <- bal_va[, keep_num, drop = FALSE]
      if (ncol(bal_tr) == 0L) {
        .log_add(step = "drop_all_constant", ok = FALSE,
                 maxV = maxV, nB = nB, fold = fold)
        next
      }
      
      # Standardize by training fold moments
      bal_tr_sc <- scale(bal_tr)
      mu <- attr(bal_tr_sc, "scaled:center"); mu[!is.finite(mu)] <- 0
      sdv <- attr(bal_tr_sc, "scaled:scale");  sdv[!is.finite(sdv) | sdv == 0] <- 1
      # 
      bal_va_sc <- sweep(sweep(bal_va, 2, mu, "-"), 2, sdv, "/")
     
      ## —— LOG: standardize information  —— ##
      .log_add(step = "scale", ok = TRUE,
               maxV = maxV, nB = nB, fold = fold,
               B_selected = ncol(bal_tr),
               anyNA_mu = any(is.na(mu)),
               anyNA_sdv = any(is.na(sdv)))
      
      fit <- try(glm(ytr ~ ., data = data.frame(ytr, bal_tr_sc), family = binomial), silent = TRUE)
      if (inherits(fit, "try-error")) {
        .log_add(step = "glm_error", ok = FALSE,
                 maxV = maxV, nB = nB, fold = fold)
        next
      }
      pr <- try(as.numeric(predict(fit, newdata = as.data.frame(bal_va_sc), type = "response")), silent = TRUE)
      if (inherits(pr, "try-error")) {
        .log_add(step = "predict_error", ok = FALSE,
                 maxV = maxV, nB = nB, fold = fold)
        next
      }#
      
      ## —— LOG: check for prediction (NA) —— ##
      .log_add(step = "predict", 
               ok = isTRUE(length(unique(round(pr, 6))) >= 2) && all(is.finite(pr)),
               maxV = maxV, nB = nB, fold = fold,
               pred_const = length(unique(round(pr, 6))) < 2,
               anyNA_pr = any(!is.finite(pr)))
      
      r   <- try(pROC::roc(response = yva, predictor = pr, levels = c("Low","High")), silent = TRUE)
      if (!inherits(r, "try-error")) {
        A <- as.numeric(pROC::auc(r))
        aucs <- c(aucs, A)
        ## —— LOG: AUC —— ##
        .log_add(step = "auc", ok = TRUE,
                 maxV = maxV, nB = nB, fold = fold,
                 auc = A)
       #
        } else {
        ## —— LOG: ROC calculation failed —— ##
        .log_add(step = "roc_error", ok = FALSE,
                 maxV = maxV, nB = nB, fold = fold)
      }
    } # folds
    
    res$AUC[gi] <- if (length(aucs)) mean(aucs, na.rm = TRUE) else NA_real_
  } # grid
  
  if (all(is.na(res$AUC))) stop("CV 未得到任何有效 AUC，请检查数据或参数网格。")
  out <- res[order(-res$AUC), , drop = FALSE]
 
  ## —— LOG: the end —— ##
  cv_log <<- tryCatch(
    if (length(.cv_log_rows)) dplyr::bind_rows(.cv_log_rows) else data.frame(),
    error = function(e) { message("log merge failed: ", conditionMessage(e)); data.frame() }
  )
  out
}



## Run CV to select best hyperparameters
cv_grid <- cv_multibalance(X_tr, y_tr, K = 5,
                           grid_maxV = c(10,20,30),
                           grid_nB  = 1:3,
                           seed = 333)
print(cv_grid)
#save(cv_grid, file="cv_grid.RData")
#load("cv_grid.RData")

best <- cv_grid[which.max(cv_grid$AUC), ]
cat("\nBest params from CV: nB =", best$nB, ", maxV =", best$maxV,
    ", AUC =", round(best$AUC, 4), "\n")

## Quick CV diagnostics
subset(cv_log, grepl("_selbal$", step) & ok == FALSE)
subset(cv_log, step == "predict" & (pred_const | anyNA_pr))
table(cv_log$step, cv_log$ok)
subset(cv_log, step == "roc_error")
cv_grid %>% arrange(desc(AUC))




###### d.  Fit with selected hyperparameters on full training set; evaluate on test ######

# Choose number of balances; for stability you may prefer 2 even if 3 is slightly higher AUC
nB <- 2
maxV <- best$maxV

# repeat the steps
remain <- colnames(X_tr)
pos_list <- list(); neg_list <- list()

for (b in seq_len(nB)){
  Xtr_b <- X_tr[, remain, drop = FALSE]
  keep <- colMeans(Xtr_b > 0) >= 0.05
  if (!any(keep)) break
  Xtr_b <- Xtr_b[, keep, drop = FALSE]
  if (ncol(Xtr_b) < 2) break
 
  sb <- try(fit_selbal_once(Xtr_b, y_tr, maxV = maxV), silent = TRUE)
  if (inherits(sb, "try-error")) break
  
  pos_list[[length(pos_list)+1]] <- sb$pos
  neg_list[[length(neg_list)+1]] <- sb$neg
  
  remain <- setdiff(remain, unique(c(sb$pos, sb$neg)))
  if (length(remain) < 2) break
}

valid_idx <- which(lengths(pos_list) > 0 & lengths(neg_list) > 0)
pos_list  <- pos_list[valid_idx]
neg_list  <- neg_list[valid_idx]
if (length(pos_list) == 0L) stop("No balances were selected. Please adjust the grid or thresholds.")

get_bal_df <- function(X, pos_list, neg_list){
  out <- lapply(seq_along(pos_list), function(i){
    pos <- pos_list[[i]]; neg <- neg_list[[i]]
    compute_balance(X[, unique(c(pos,neg)), drop = FALSE], pos, neg)
  })
  if (length(out) == 0) stop("No balances were selected.")
  out <- do.call(cbind, out)
  colnames(out) <- paste0("B", seq_len(ncol(out)))
  as.data.frame(out)
}

bal_tr <- get_bal_df(X_tr, pos_list, neg_list)
bal_te <- get_bal_df(X_te, pos_list, neg_list)

# Standardize using training set
bal_tr_sc <- scale(bal_tr)
mu <- attr(bal_tr_sc, "scaled:center")
sdv <- attr(bal_tr_sc, "scaled:scale")
sdv[sdv == 0 | is.na(sdv)] <- 1# set the sdv of 0/NA as 1，for avoiding devide the 0
bal_te_sc <- sweep(sweep(bal_te, 2, mu, "-"), 2, sdv, "/")

# Final logistic model (covariates can be added here if needed)
fit_final <- glm(y_tr ~ ., data = data.frame(y_tr, bal_tr_sc), family = binomial)
cat("\nFinal model summary:\n")
print(summary(fit_final))
anova(fit_final)
summary(fit_final)$coefficients

# Overall model significance (LRT vs null)
fit_null <- glm(y_tr ~ 1, data = data.frame(y_tr, bal_tr_sc), family = binomial)
anova_result <- anova(fit_null, fit_final, test = "LRT")
print(anova_result)

# Get tye P value of this model
model_p_value <- anova_result$`Pr(>Chi)`[2]
cat("Overall model P-value:", model_p_value, "\n")

# Test-set evaluation
prob_te <- as.numeric(predict(fit_final, newdata = as.data.frame(bal_te_sc), type = "response"))
roc_te  <- pROC::roc(response = y_te, predictor = prob_te, levels = c("Low","High"))
cat("\nTest AUC:", round(pROC::auc(roc_te), 4), "\n")
# AUC: 0.6307
# B1+B2, AUC：0.6205

opt <- pROC::coords(roc_te, "best", best.method = "youden", transpose = TRUE)
thr <- opt["threshold"]
pred_te <- ifelse(prob_te >= thr, "High", "Low")
pred_te <- factor(pred_te, levels = c("Low","High"))

cat("\nConfusion matrix (test):\n")
print(table(Pred = pred_te, Truth = y_te))

cat("\nSelected balances (numerator/denominator):\n")
for (i in seq_along(pos_list)) {
  cat(paste0("B", i, ":\n  Numerator:   "), paste(pos_list[[i]], collapse = ", "), "\n",
      "  Denominator: ", paste(neg_list[[i]], collapse = ", "), "\n\n", sep = "")
}




###### e. ROC plot ######
# selectt the balance B1 and balance B2
cat("Test AUC:", round(pROC::auc(roc_te), 4), "\n")
getS3method("print", "roc", asNamespace("pROC"))(roc_te)
# Or (less recommended)
pROC:::print.roc(roc_te)
p_roc <- pROC::plot.roc(roc_te, print.auc = TRUE, col = "#2c7fb8")
stopifnot(length(prob_te) == length(y_te), !anyNA(prob_te), !anyNA(y_te))

# Base ROC (ggplot)
roc_df <- data.frame(
  fpr = 1 - roc_te$specificities,
  tpr = roc_te$sensitivities
)
auc_val <- as.numeric(pROC::auc(roc_te))

p_roc_gg <- ggplot(roc_df, aes(fpr, tpr)) +
  geom_path(linewidth = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 3) +
  annotate("text", x = 0.75, y = 0.15,
           label = sprintf("AUC-ROC\n%.3f", auc_val), size = 3.2) +
  labs(title = "ROC curve", x = "FPR", y = "TPR") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(hjust = 0.5, size = 9))
p_roc_gg

# Combine B1+B2 taxa IDs (for top annotation if needed)
num_ids2 <- unique(unlist(pos_list[1:2]))
den_ids2 <- unique(unlist(neg_list[1:2]))

# Use final model linear predictor as score (train/test)
lp_tr <- as.numeric(predict(fit_final, newdata = as.data.frame(bal_tr_sc), type = "link"))
lp_te <- as.numeric(predict(fit_final, newdata = as.data.frame(bal_te_sc), type = "link"))
df_all_lp <- rbind(
  data.frame(Balance = lp_tr, Group = y_tr, Set = "Train"),
  data.frame(Balance = lp_te, Group = y_te, Set = "Test")
)

# Boxplot (test set)
df_all_lp$Group <- factor(df_all_lp$Group, levels = c("High","Low"))
p_box <- df_all_lp%>%
  subset(Set=="Test")%>%
  ggplot( aes(x = Group, y = Balance, fill = Group)) +
  geom_boxplot(width = 0.65, alpha = 0.75, outlier.size = 0.9) +
  scale_fill_manual(values = c("High"="#E69F00","Low"="#56B4E9")) +
  labs(x = "Factor", y = "Balance") +
  theme_bw(base_size = 13) + theme(legend.position = "none")+
  stat_compare_means()
p_box

# Density (test set), horizontal for paneling symmetry
p_dens <- df_all_lp%>%
  subset(Set=="Test")%>%
  ggplot(aes(x = Balance, fill = Group)) +
  geom_density(alpha = .35) +
  scale_fill_manual(values = c("High"="#E69F00","Low"="#56B4E9")) +
  coord_flip() +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(x = NULL, y = NULL)
p_dens 

# Optional helper to wrap IDs for top annotation panels
wrap_ids <- function(x, n = 6) {
  if (length(x) == 0) return("")
  paste(sapply(split(x, ceiling(seq_along(x)/n)), paste, collapse = " "),
        collapse = "\n")
}

#
txt_b1 <- if (length(pos_list) >= 1) sprintf(
  "B1 DENOMINATOR\n%s\n\nB1 NUMERATOR\n%s",
  wrap_ids(neg_list[[1]], n = 6),
  wrap_ids(pos_list[[1]], n = 6)
) else ""

txt_b2 <- if (length(pos_list) >= 2) sprintf(
  "B2 DENOMINATOR\n%s\n\nB2 NUMERATOR\n%s",
  wrap_ids(neg_list[[2]], n = 6),
  wrap_ids(pos_list[[2]], n = 6)
) else ""

# Map balance IDs to taxonomy (for reporting tables)
balance_taxa_df <- function(phy, ids){
  
  ids_clean <- sub("^x", "", ids)
  
  tt <- as.data.frame(phyloseq::tax_table(phy), stringsAsFactors = FALSE)
  tt$ASV <- rownames(tt)
  
  ids_keep <- ids_clean[ids_clean %in% tt$ASV]
  if (length(ids_keep) == 0) {
    return(data.frame(ASV=character(0)))
  }
  out <- tt[match(ids_keep, tt$ASV), , drop = FALSE]
  
  rank_names <- colnames(tt)
  has_species <- any(grepl("species", tolower(rank_names)))
  has_genus   <- any(grepl("genus",   tolower(rank_names)))
  sp_col <- if (has_species) rank_names[grep("species", tolower(rank_names))[1]] else NA
  ge_col <- if (has_genus)   rank_names[grep("genus",   tolower(rank_names))[1]] else NA
  
  best_name <- if (!is.na(sp_col)) {
    paste(out[[ge_col]], out[[sp_col]])
  } else if (!is.na(ge_col)) {
    out[[ge_col]]
  } else {
    out$ASV
  }
  out$BestName <- best_name
  
  out <- out[, c("ASV", "BestName", setdiff(colnames(out), c("ASV", "BestName")))]
  rownames(out) <- NULL
  out
}

# Example taxonomy tables for B1/B2 (print or export as needed)
b1_num <- balance_taxa_df(phy_all6_species_1y_4yWeek, pos_list[[1]])
b1_num$species
b1_den <- balance_taxa_df(phy_all6_species_1y_4yWeek, neg_list[[1]])
b1_den$species
b2_num <- if (length(pos_list) >= 2) balance_taxa_df(phy_all6_species_1y_4yWeek, pos_list[[2]]) else NULL
b2_num$species
b2_den <- if (length(neg_list) >= 2) balance_taxa_df(phy_all6_species_1y_4yWeek, neg_list[[2]]) else NULL
b2_den$species

# Check
cat("\n=== B1 NUMERATOR ===\n"); print(b1_num)
cat("\n=== B1 DENOMINATOR ===\n"); print(b1_den)
if (!is.null(b2_num)) { cat("\n=== B2 NUMERATOR ===\n"); print(b2_num) }
if (!is.null(b2_den)) { cat("\n=== B2 DENOMINATOR ===\n"); print(b2_den) }

# Final panel (boxplot + ROC). Save as 9 x 3.2 inches if needed.
p_roc_gg#
p_box#

final_plot <- ggarrange(p_box,p_roc_gg+rremove("ylab"),
                        ncol=2,nrow=1, widths = c(2.8,1.2),
                        labels=NULL,
                        common.legend = TRUE, legend="bottom")
final_plot




#####  Figure 5E.  Co-occurrence network of the top 99 bacterial species #####
phy_all6_species #

### prepare data
phy_all6_species_1y_4yWeek <- phy_all6_species%>%
  subset_samples(Time_new =="1y"& !is.na(Maturity_4y))# 371 samples

phy_all6_species_1y_4yWeek  <- prune_taxa(taxa_sums(phy_all6_species_1y_4yWeek) > 0, 
                                          phy_all6_species_1y_4yWeek)

# Extract OTU_table to matrix and add pseudocount and normalize
all6_species_1y_4yWeek_sample_df <- data.frame(
  SampleID=get_variable(phy_all6_species_1y_4yWeek,"SampleID"),
  Maturity_4y=get_variable(phy_all6_species_1y_4yWeek,"Maturity_4y")%>%as.character)
rownames(all6_species_1y_4yWeek_sample_df) <- all6_species_1y_4yWeek_sample_df$SampleID

# Relative abundance OTU matrix (rows: taxa, columns: samples)
all6_species_1y_4yWeek <- transform_sample_counts(phy_all6_species_1y_4yWeek,function(x) x/sum(x))%>%
  otu_table%>%
  as("matrix")%>%
  as.data.frame


###### a. 99 genera related network ####

### High/Low_spearman correlation
all6_species_1y_4yWeek # 371 samples * 789 taxa
taxa_names_99_speciesname #

# Filter taxa that are present in at least one sample
otu_filter <- all6_species_1y_4yWeek[rowSums(all6_species_1y_4yWeek > 0) > 0, ]

# Keep only the top 99 RF-selected taxa
otu_filter_99 <- otu_filter[
  rownames(otu_filter) %in% rownames(taxa_names_99_speciesname),
  ,
  drop = FALSE
]

# Check
dim(otu_filter_99)
length(intersect(rownames(otu_filter), rownames(taxa_names_99_speciesname)))

# Convert to samples x taxa
otu_filter_99_sample <- t(otu_filter_99)

# Remove samples with zero total abundance, if any
otu_filter_99_sample <- otu_filter_99_sample[
  rowSums(otu_filter_99_sample, na.rm = TRUE) > 0,
  ,
  drop = FALSE
]

# Re-normalize after subsetting to the 99-taxon composition
otu_filter_99_sample <- sweep(
  otu_filter_99_sample,
  1,
  rowSums(otu_filter_99_sample, na.rm = TRUE),
  "/"
)

# Zero replacement before CLR transformation
otu_filter_99_repl <- zCompositions::cmultRepl(
  as.matrix(otu_filter_99_sample),
  label = 0,
  method = "CZM",
  output = "prop",
  z.delete = FALSE
)

# CLR transformation
otu_filter_99_clr <- as.data.frame(compositions::clr(otu_filter_99_repl))
colnames(otu_filter_99_clr) <- colnames(otu_filter_99_sample)
rownames(otu_filter_99_clr) <- rownames(otu_filter_99_sample)

# Calculate Spearman correlation on CLR-transformed abundance profiles
res <- Hmisc::rcorr(as.matrix(otu_filter_99_clr), type = "spearman")

# Function for converting matrix to long format
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
cor <- flattenCorrMatrix(res$r,res$P)

# Adjust p-values using Benjamini-Hochberg correction
cor$p.adj <- p.adjust(cor$p,method = "BH")
range(cor$p.adj);range(cor$r)
names(cor) <- c("Source","Target","r","p","p.adj")
cor$linktype <- ifelse(cor$r>0,"1","-1")
cor$abs <- abs(cor$r)

# Subset correlations involving the 99 key taxa
cor_only99 <- subset(cor,Source %in% rownames(taxa_names_99_speciesname)&Target %in% 
                       rownames(taxa_names_99_speciesname))
dim(cor_only99) # 4851 * 7
length(unique(union(unique(cor_only99$Target),unique(cor_only99$Source))))

# Determine the threshold for significant interspecies interactions # abs(x) computes the absolute value of x
r_only99_sig <- subset(cor_only99,p.adj<0.05&r>0.40)
dim(r_only99_sig) # 161 * 7, # 345 * 7 # 171 * 7
length(unique(union(unique(r_only99_sig$Target),unique(r_only99_sig$Source)))) # 65

# Add taxa information to table
taxa_info<- data.frame(tax_table(phy_all6_species_1y_4yWeek))
taxa_info$OTU_ID <- rownames(taxa_info)

# Match OTU IDs to species names for the 99 taxa network
matching_indices1 <- match(r_only99_sig$Target, taxa_info$OTU_ID)
matching_indices2 <- match(r_only99_sig$Source, taxa_info$OTU_ID)

r_only99_sig$Target2[!is.na(matching_indices1)] <- taxa_info$species[matching_indices1[!is.na(matching_indices1)]]
r_only99_sig$Source2[!is.na(matching_indices2)] <- taxa_info$species[matching_indices2[!is.na(matching_indices2)]]

# Save result for visualization in Gephi
write_xlsx(r_only99_sig,"D:/OneDrive - University of Copenhagen/01Project/COPSAC/COPSAC2010_feces/copsac2010_feces_species/network_3/cooccurrence_only99species_p0.05_r0.4.xlsx")


###### b. Operations in Gephi ####
# Result note: 4 consortium, R = 4
# Then get a file of bacterial cooccurrence moduarity. file name: cooccurrence_only99species_p0.05_r0.3_Moduarity.xlsx



###### c. Correlation network among bacterial clusters based on Spearman correlation analysis ####
all6_species_1y_4yWeek_otu_df #
all6_species_1y_4yWeek_sample_df #
taxa_names_99_speciesname #

library(readxl)
library(zCompositions)
library(compositions)

MCODE_bac <- read_excel("D:/OneDrive - University of Copenhagen/01Project/COPSAC/COPSAC2010_feces/copsac2010_feces_species/network_3/cooccurrence_only99species_p0.05_r0.4_Moduarity.xlsx")
MCODE_cluster1 <- MCODE_bac$Taxa[MCODE_bac$Modularity=="C1"]
MCODE_cluster2 <- MCODE_bac$Taxa[MCODE_bac$Modularity=="C2"]
MCODE_cluster3 <- MCODE_bac$Taxa[MCODE_bac$Modularity=="C3"]
MCODE_cluster4 <- MCODE_bac$Taxa[MCODE_bac$Modularity=="C4"]

# Subset OTU table to include only the 99 predictor species
test <- all6_species_1y_4yWeek_otu_df[,colnames(all6_species_1y_4yWeek_otu_df) %in% 
                                        #taxa_names_99_speciesname$predictors2
                                        taxa_names_99_speciesname$predictors
]

# Replace column names in test with corresponding species names
colnames(test) <- taxa_names_99_speciesname[match(colnames(test),
                                                  taxa_names_99_speciesname$predictors),"species"]

# Sum raw relative abundance per module, then CLR transformation
otu_cluster1 <- rowSums(test %>% dplyr::select(all_of(MCODE_cluster1)), na.rm = TRUE) %>% as.data.frame
colnames(otu_cluster1) <- "C1"
otu_cluster2 <- rowSums(test %>% dplyr::select(all_of(MCODE_cluster2)), na.rm = TRUE) %>% as.data.frame
colnames(otu_cluster2) <- "C2"
otu_cluster3 <- rowSums(test %>% dplyr::select(all_of(MCODE_cluster3)), na.rm = TRUE) %>% as.data.frame
colnames(otu_cluster3) <- "C3"
otu_cluster4 <- rowSums(test %>% dplyr::select(all_of(MCODE_cluster4)), na.rm = TRUE) %>% as.data.frame
colnames(otu_cluster4) <- "C4"


# Merge all cluster data frames by row names
merge.all <- function(x, ..., by = "row.names") {
  L <- list(...)
  for (i in seq_along(L)) {
    x <- merge(x, L[[i]], by = by)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
  }
  return(x)
}

otu_clusters_mid <- merge.all(otu_cluster1, otu_cluster2, otu_cluster3, otu_cluster4)

# Remove samples with zero total module abundance, if any
otu_clusters_mid <- otu_clusters_mid[
  rowSums(otu_clusters_mid, na.rm = TRUE) > 0,
  ,
  drop = FALSE
]

# Re-normalize C1-C4 into compositional data
otu_clusters_prop <- sweep(
  otu_clusters_mid,
  1,
  rowSums(otu_clusters_mid, na.rm = TRUE),
  "/"
)

# Zero replacement before CLR
otu_clusters_repl <- as.matrix(otu_clusters_prop)

# CLR transformation
otu_clusters_clr <- as.data.frame(compositions::clr(otu_clusters_repl))

colnames(otu_clusters_clr) <- colnames(otu_clusters_prop)
rownames(otu_clusters_clr) <- rownames(otu_clusters_prop)

# Keep your original downstream structure
otu_clusters_mid2 <- data.frame(t(otu_clusters_clr))

# Spearman correlation between modules
library(Hmisc)
res_community <- rcorr(t(otu_clusters_mid2), type = "spearman")

# Function to convert correlation matrix to long format
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

cor_community <- flattenCorrMatrix(res_community$r,res_community$P)

# Adjust p-values using Benjamini-Hochberg method
cor_community$p.adj <- p.adjust(cor_community$p,method = "BH")
range(cor_community$p.adj);range(cor_community$r)

# Rename columns and add correlation metadata
names(cor_community) <- c("Source","Target","r","p","p.adj")
cor_community$linktype <- ifelse(cor_community$r>0,"1","-1")
cor_community$abs <- abs(cor_community$r)

# Filter significant correlations (adjusted p < 0.05)
r_community_sig <- subset(cor_community,p.adj<0.05&abs>0.4)

# Save results
write_xlsx(r_community_sig,"D:/OneDrive - University of Copenhagen/01Project/COPSAC/COPSAC2010_feces/copsac2010_feces_species/network_3/cooccurrence_community_allspecies_p0.05.xlsx")

### Operations in Gephi, and represent the plot as Figure S13




#####  Figure 5F.  Co-occurrence network of the top 99 bacterial species #####
otu_clusters_clr # CLR-transformed module abundance, used for plotting
otu_clusters_mid # actual module relative abundance, used for statistics
phy_all6_species_1y_4yWeek #

# Prepare metadata for sample and maturity status
all6_species_1y_4yWeek_sample_df <- data.frame(
  SampleID = get_variable(phy_all6_species_1y_4yWeek, "SampleID"),
  Maturity_4y = get_variable(phy_all6_species_1y_4yWeek, "Maturity_4y") %>% as.character())

rownames(all6_species_1y_4yWeek_sample_df) <- all6_species_1y_4yWeek_sample_df$SampleID

# Keep the same samples in CLR data, actual relative abundance data, and metadata
common_ids <- Reduce(intersect,list(rownames(otu_clusters_clr),
                                    rownames(otu_clusters_mid),
                                    rownames(all6_species_1y_4yWeek_sample_df) ))

#  CLR data for plotting
otu_clusters_plot <- cbind(otu_clusters_clr[common_ids, c("C1", "C2", "C3", "C4"), drop = FALSE],
                           Maturity_4y = all6_species_1y_4yWeek_sample_df[common_ids, "Maturity_4y"]) %>%
  as.data.frame()

otu_clusters_plot <- subset(otu_clusters_plot, !is.na(Maturity_4y))

otu_clusters_plot$Maturity_4y <- factor(otu_clusters_plot$Maturity_4y,levels = c("High","Low"))

#  Actual relative abundance data for statistics
otu_clusters_RA <- cbind(otu_clusters_mid[common_ids, c("C1", "C2", "C3", "C4"), drop = FALSE],
  Maturity_4y = all6_species_1y_4yWeek_sample_df[common_ids, "Maturity_4y"]) %>%
  as.data.frame()

otu_clusters_RA <- subset(otu_clusters_RA, !is.na(Maturity_4y))

otu_clusters_RA$Maturity_4y <- factor(otu_clusters_RA$Maturity_4y,levels = c("High","Low"))

# Check sample numbers
table(otu_clusters_plot$Maturity_4y, useNA = "ifany")
table(otu_clusters_RA$Maturity_4y, useNA = "ifany")

# Statistical comparison using actual relative abundance
otu_clusters_RA_long <- otu_clusters_RA %>%
  pivot_longer(cols = C1:C4,
               names_to = "cluster",
               values_to = "total_RA") %>%
  dplyr::select(cluster, total_RA, Maturity_4y)

cluster_stats_all_C_RA <- otu_clusters_RA_long %>%
  group_by(cluster) %>%
  summarise( p_value = wilcox.test(total_RA ~ Maturity_4y, exact = FALSE)$p.value,
             mean_high = mean(total_RA[Maturity_4y == "High"], na.rm = TRUE),
             mean_low = mean(total_RA[Maturity_4y == "Low"], na.rm = TRUE),
             median_high = median(total_RA[Maturity_4y == "High"], na.rm = TRUE),
             median_low = median(total_RA[Maturity_4y == "Low"], na.rm = TRUE),
             fold_change = ifelse(mean_low == 0, NA, mean_high / mean_low),
             .groups = "drop" ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))
print(cluster_stats_all_C_RA)

# Optional: save statistics table
writexl::write_xlsx(cluster_stats_all_C_RA,
  "D:/OneDrive - University of Copenhagen/01Project/COPSAC/COPSAC2010_feces/copsac2010_feces_species/network_3/Figure5F_module_RA_statistics.xlsx")

#  Plotting using CLR-transformed module abundance
col_maturity <- c("High" = "#E69F00", "Low" = "#56B4E9")

format_p <- function(p) {
  ifelse(is.na(p), "NA",
    ifelse(p < 0.001, "< 0.001", as.character(signif(p, 3))))
}

plot_cluster <- function(cluster_name) {
  
  p_value_use <- cluster_stats_all_C_RA %>%
    dplyr::filter(cluster == cluster_name) %>%
    dplyr::pull(p_value)
  
  p_adj_use <- cluster_stats_all_C_RA %>%
    dplyr::filter(cluster == cluster_name) %>%
    dplyr::pull(p_adj)
  
  y_range <- range(otu_clusters_plot[[cluster_name]], na.rm = TRUE)
  y_pos <- y_range[2] + 0.10 * diff(y_range)
  
  ggplot( otu_clusters_plot,
          aes(x = Maturity_4y, y = .data[[cluster_name]]) ) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_jitter(aes(color = Maturity_4y),
                size = 3,
                alpha = 0.4,
                width = 0.2 ) +
    geom_boxplot(aes(fill = Maturity_4y),
                 alpha = 0.1,
                 outlier.shape = NA) +
    scale_fill_manual(values = col_maturity) +
    scale_color_manual(values = col_maturity) +
    xlab(cluster_name) +
    ylab("CLR-transformed module abundance") +
    annotate("text",
             x = 1.5,
             y = y_pos,
             label = paste0("P = ", format_p(p_adj_use)),
             size = 3.5) +
    theme(text = element_text(size = 11),
          legend.position = "bottom")
}

C1 <- plot_cluster("C1")
C2 <- plot_cluster("C2")
C3 <- plot_cluster("C3")
C4 <- plot_cluster("C4")

# Arrange plots
C_sig <- ggarrange(C1,C2 + rremove("ylab"), C4 + rremove("ylab"),
  ncol = 3,nrow = 1,labels = NULL,
  common.legend = TRUE,legend = "bottom")
C_sig
# save as 6 * 3.2

C_nosig <- ggarrange(C3,
  ncol = 1,nrow = 1,labels = NULL,
  common.legend = TRUE, legend = "bottom")
C_nosig
# save as 3 * 3.6














##### SparCC validation of associations among the top 99 predictive species #####

library(SpiecEasi)
library(pheatmap)

### 1. Prepare top 99 taxa count table: samples x taxa

# Extract count table from phyloseq object
otu_count <- as(otu_table(phy_all6_species_1y_4yWeek), "matrix")

# Make sure taxa are rows
if (!taxa_are_rows(phy_all6_species_1y_4yWeek)) {
  otu_count <- t(otu_count)
}

# Keep only top 99 RF-selected taxa
taxa_99_ids <- rownames(taxa_names_99_speciesname)

otu_count_99 <- otu_count[
  rownames(otu_count) %in% taxa_99_ids,
  ,
  drop = FALSE
]

# Convert to samples x taxa
otu_count_99_sample <- t(otu_count_99)

# Remove all-zero samples or taxa, if any
otu_count_99_sample <- otu_count_99_sample[
  rowSums(otu_count_99_sample, na.rm = TRUE) > 0,
  ,
  drop = FALSE
]

otu_count_99_sample <- otu_count_99_sample[
  ,
  colSums(otu_count_99_sample, na.rm = TRUE) > 0,
  drop = FALSE
]

# Check dimensions
dim(otu_count_99_sample)
length(intersect(colnames(otu_count_99_sample), taxa_99_ids))



### 2. Rename taxa IDs to species names for plotting

tax_info_99 <- data.frame(
  taxa_id = rownames(taxa_names_99_speciesname),
  species = taxa_names_99_speciesname$species,
  stringsAsFactors = FALSE)

tax_info_99$species_plot <- make.unique(tax_info_99$species)

match_id <- match(colnames(otu_count_99_sample), tax_info_99$taxa_id)

colnames(otu_count_99_sample) <- tax_info_99$species_plot[match_id]

# Check
head(colnames(otu_count_99_sample))


### 3. Run SparCC

set.seed(333)

sparcc_res <- SpiecEasi::sparcc(
  otu_count_99_sample,
  iter = 20,
  inner_iter = 10,
  th = 0.1
)

sparcc_cor <- sparcc_res$Cor

# Make sure row/column names are retained
rownames(sparcc_cor) <- colnames(otu_count_99_sample)
colnames(sparcc_cor) <- colnames(otu_count_99_sample)

# Quick check
range(sparcc_cor, na.rm = TRUE)


### 4. Add C1-C4 module annotation to full SparCC heatmap
### Sort heatmap by Module

# Keep all SparCC correlations
sparcc_cor_all <- sparcc_cor
diag(sparcc_cor_all) <- 1

module_df <- data.frame(
  species = c(MCODE_cluster1, MCODE_cluster2, MCODE_cluster3, MCODE_cluster4),
  Module = c(rep("C1", length(MCODE_cluster1)),
             rep("C2", length(MCODE_cluster2)),
             rep("C3", length(MCODE_cluster3)),
             rep("C4", length(MCODE_cluster4))),stringsAsFactors = FALSE)

# Because heatmap names may be make.unique(species), remove duplicated suffix if needed
module_df$species_base <- module_df$species

annotation_row <- data.frame(
  Module = module_df$Module[
    match( sub("\\.\\d+$", "", rownames(sparcc_cor_all)),
      module_df$species_base)], stringsAsFactors = FALSE)

rownames(annotation_row) <- rownames(sparcc_cor_all)

# Keep only taxa with module annotation
common_heat_taxa <- rownames(annotation_row)[!is.na(annotation_row$Module)]

sparcc_cor_all_sub <- sparcc_cor_all[common_heat_taxa,common_heat_taxa,drop = FALSE]

annotation_row_all <- annotation_row[common_heat_taxa,,drop = FALSE]

# Sort taxa by Module
annotation_row_all$Module <- factor(annotation_row_all$Module,levels = c("C1", "C2", "C3", "C4"))

taxa_order <- rownames(annotation_row_all)[
  order(annotation_row_all$Module)]

sparcc_cor_all_sub <- sparcc_cor_all_sub[taxa_order,taxa_order,drop = FALSE]

annotation_row_all <- annotation_row_all[ taxa_order,,drop = FALSE]

# Add gaps between modules
module_sizes <- table(annotation_row_all$Module)
gaps_module <- cumsum(module_sizes)
gaps_module <- gaps_module[-length(gaps_module)]

# Optional: module colors
module_colors <- list( Module = c("C1" = "#BABB26","C2" = "#11398A",
                                  "C3" = "#71BEC0","C4" = "#D26A1B"))


p_SparCC_top99_all_association_heatmap_with_modules <- pheatmap::pheatmap(
  sparcc_cor_all_sub,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  annotation_row = annotation_row_all,
  annotation_col = annotation_row_all,
  annotation_colors = module_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gaps_module,
  gaps_col = gaps_module,
  border_color = NA,
  fontsize_row = 6,
  fontsize_col = 6,
  main = "All SparCC associations among top 99 predictive species" )
p_SparCC_top99_all_association_heatmap_with_modules
# save as 9 * 8




#####  Figure 5G.  Module-specific metabolic pathways by functional category #####
phy_all6_species_1y_4yWeek #
path_contrib #
MCODE_cluster1#
MCODE_cluster2#
MCODE_cluster3#
MCODE_cluster4#

# Export ASV table and sequences for PICRUSt2
ps_asv <- phy_all6_species_1y_4yWeek
otu <- otu_table(ps_asv)
if (!taxa_are_rows(ps_asv)) otu <- t(otu) # PICRUSt2 requires the features of being arranged in rows
bm  <- biomformat::make_biom(data = as.matrix(otu))
biomformat::write_biom(bm, "asv.biom")

seqs <- Biostrings::DNAStringSet(refseq(ps_asv))
names(seqs) <- taxa_names(ps_asv)
Biostrings::writeXStringSet(seqs, "asv.fna")

# Quick check: if the row names and sequence names correspond one by one?
stopifnot(identical(rownames(as(otu, "matrix")), names(seqs)))


### On Linux server
### Run PICRUSt2 externally (command-line)
# Code: picrust2_pipeline.py -s asv.fna -i asv.biom -o picrust2_out -p 8


## Load pathway contribution data
path_contrib <- read.table("./picrust2_out/pathways_out/path_abun_contrib.tsv.gz",
                           header=TRUE, sep="\t", check.names=FALSE, comment.char="")

# Merge species annotation
taxon_table <- data.frame(tax_table( ps_asv))%>% dplyr::select(species)
taxon_table$taxon <- rownames(taxon_table)
path_contrib <- merge(path_contrib, taxon_table, by=c("taxon"))

# Merge sample metadata
sample_info <- data.frame(sample_data(phy_all6_species_1y_4yWeek))%>% 
  dplyr::select(Maturity_4y,SampleID)
path_contrib<- dplyr::rename(path_contrib, SampleID=sample)

path_contrib <- merge(path_contrib, sample_info, by=c("SampleID"))

# Check how many species from each MCODE cluster are present
table(MCODE_cluster1%in%path_contrib$species) # TRUE 16
table(MCODE_cluster2%in%path_contrib$species)# TRUE 9
table(MCODE_cluster3%in%path_contrib$species)# TRUE 3
table(MCODE_cluster4%in%path_contrib$species)# TRUE 9

path_contrib$cluster <- ifelse (
  path_contrib$species %in% MCODE_cluster1,"C1",
  ifelse(path_contrib$species %in% MCODE_cluster2,"C2",
         ifelse(path_contrib$species %in% MCODE_cluster3,"C3",
                ifelse(path_contrib$species %in% MCODE_cluster4,"C4",NA))))
table(path_contrib$cluster)

# Subset data to cluster-assigned species only
path_contrib_community <- path_contrib%>%
  subset(cluster %in%c("C1","C2","C3","C4"))

# Rename column for clarity
path_contrib_community <- rename(path_contrib_community, "function_name"="function")

# Calculate total relative abundance of each function per sample and cluster
path_contrib_community2 <- path_contrib_community %>%
  group_by(SampleID, cluster) %>%
  dplyr::mutate(total_taxon_rel_function_abun = sum(taxon_rel_function_abun, na.rm = TRUE))%>%
  ungroup()

# Identify cluster-specific functions
cluster_specific_functions <- path_contrib_community2 %>%
  group_by(cluster, function_name) %>%
  summarise(mean_abundance = mean(total_taxon_rel_function_abun, na.rm = TRUE),
            prevalence = n_distinct(SampleID),
            .groups = "drop" ) %>%
  group_by(function_name) %>%
  mutate(total_abundance = sum(mean_abundance),
         cluster_specificity = mean_abundance / total_abundance * 100) %>%
  ungroup()
cluster_specific_functions

# Save results
write_xlsx(cluster_specific_functions, "cluster_specific_functions.xlsx")

# Extract top 5 most specific functions per cluster
top_cluster_functions <- cluster_specific_functions %>%
  group_by(cluster) %>%
  arrange(cluster, desc(cluster_specificity)) %>%
  slice_head(n = 10) %>%  # Select the top 5 most unique features from each cluster
  ungroup()


## Load pathway description file

### Retrieve pathway descriptions using the KEGG API (if corresponding KEGG IDs are available)
pathway_descrip <- read.table("./picrust2_out/pathways_out/path_abun_unstrat_descrip.tsv.gz",
                              header=TRUE, sep="\t", check.names=FALSE, comment.char="")

# Initialize and update description
table(top_cluster_functions$function_name %in% pathway_descrip$pathway)
# "PWY-5747" "PWY-6107" "PWY-6654" "PWY-6944" "PWY-5823"
table("PWY1G-6944" %in% pathway_descrip$pathway)

top_cluster_functions$description <- ifelse(top_cluster_functions$function_name %in% pathway_descrip$pathway,
                                            pathway_descrip$description,top_cluster_functions$function_name)

# Step 1: Check the data
print("=== Data Check ===")
print(paste("Number of rows in top_cluster_functions:", nrow(top_cluster_functions)))
print(paste("Number of rows in pathway_descrip:", nrow(pathway_descrip)))
print("Unique pathways in top_cluster_functions:")
print(unique(top_cluster_functions$function_name))

# Step 2: Initialize the description column
top_cluster_functions$description <- top_cluster_functions$function_name

# Step 3: Exact matching
matched_count <- 0
for(i in 1:nrow(top_cluster_functions)) {
  pathway_id <- top_cluster_functions$function_name[i]
  if(pathway_id %in% pathway_descrip$pathway) {
    desc <- pathway_descrip$description[pathway_descrip$pathway == pathway_id][1]
    top_cluster_functions$description[i] <- desc
    matched_count <- matched_count + 1
  }
}
print(paste("Exact matches:", matched_count, "pathways"))

# Step 4: Handle unmatched pathways
unmatched <- top_cluster_functions$function_name[!top_cluster_functions$function_name %in% pathway_descrip$pathway]
print(paste("Number of unmatched pathways", length(unique(unmatched))))
unmatched

# Manually verify unmatched pathways
pathway_verification <- data.frame(
  pathway_id = c("PWY-5747", "PWY-6107","PWY-6654","PWY-6944",  "PWY-5823"
                 ),
  official_url = c(
    "https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=PWY-5747",
    "https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=PWY-6107", 
    "https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=PWY-6654",
    "https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=PWY-6944",
    "https://metacyc.org/META/NEW-IMAGE?type=PATHWAY&object=PWY-5823"
  ),
  stringsAsFactors = FALSE
)
print("Please verify the pathway descriptions at the following URLs:")
for(i in 1:nrow(pathway_verification)) {
  cat(paste0(pathway_verification$pathway_id[i], ": ", pathway_verification$official_url[i], "\n"))
}

# Step 5: Manually map common MetaCyc pathways
manual_descriptions <- list(
  "PWY-5747" = "2-methylcitrate cycle II",
  "PWY-6107" = "chlorosalicylate degradation",
  "PWY-6654" = "phosphopantothenate biosynthesis III (archaea)",
  "PWY-6944" = "androstenedione degradation I (aerobic)",
  "PWY-5823" = "superpathway of CDP-glucose-derived O-antigen building blocks biosynthesis")

# Apply manual mappings
for(pathway_id in names(manual_descriptions)) {
  if(pathway_id %in% top_cluster_functions$function_name) {
    top_cluster_functions$description[top_cluster_functions$function_name == pathway_id] <- manual_descriptions[[pathway_id]]
  }
}

# Step 6: Final verification
print("=== Final Result ===")
print("Pathway descriptions:")
result_summary <- top_cluster_functions %>%
  dplyr::select(function_name, description) %>%
  distinct() %>%
  arrange(function_name)
print(result_summary)

# Define pathway description fix for HTML character codes
top_cluster_functions$description <-
  ifelse(top_cluster_functions$description=="3-phenylpropanoate and 3-(3-hydroxyphenyl)propanoate degradation to 2-oxopent-4-enoate","Phenylpropanoate/hydroxyphenylpropanoate degradation",
         ifelse(top_cluster_functions$description=="catechol degradation to &beta;-ketoadipate","catechol degradation to β-ketoadipate",
                ifelse(top_cluster_functions$description=="aromatic compounds degradation via &beta;-ketoadipate","Aromatic compound degradation",
                top_cluster_functions$description)))

# Recreate matrix
pathway_matrix <- top_cluster_functions %>%
  dplyr::select(cluster, description, cluster_specificity) %>%
  pivot_wider(
    names_from = cluster, 
    values_from = cluster_specificity, 
    values_fill = 0
  ) %>%
  column_to_rownames("description") %>%
  as.matrix()

print("Heatmap matrix dimensions:")
print(dim(pathway_matrix))

# Create full classification
# Accurate category mapping from MetaCyc and KEGG
# MetaCyc: https://metacyc.org/
# KEGG: https://www.kegg.jp/kegg/pathway.html

# top 10 * 4
pathway_classification <- data.frame(
  Pathway = c(
    "L-arginine degradation II (AST pathway)",
    "catechol degradation to β-ketoadipate",
    "glucose degradation (oxidative)",
    "Phenylpropanoate/hydroxyphenylpropanoate degradation",
    "benzoyl-CoA degradation I (aerobic)",
    
    "vitamin E biosynthesis (tocopherols)",
    "catechol degradation III (ortho-cleavage pathway)",
    "Aromatic compound degradation",
    "2-methylcitrate cycle II",
    "chlorosalicylate degradation",
    
    "1,4-dihydroxy-6-naphthoate biosynthesis II",
    "superpathway of demethylmenaquinol-6 biosynthesis II",
    "1,4-dihydroxy-6-naphthoate biosynthesis I",
    "methanogenesis from H2 and CO2",
    "CDP-archaeol biosynthesis",
    "archaetidylinositol biosynthesis",
    "phosphopantothenate biosynthesis III (archaea)",
    
    "7-(3-amino-3-carboxypropyl)-wyosine biosynthesis",
    "sucrose biosynthesis III",
    "sucrose biosynthesis I (from photosynthesis)",
    "sucrose degradation II (sucrose synthase)",
    "androstenedione degradation I (aerobic)",
    
    "coenzyme B biosynthesis",
    "L-valine degradation I",
    "meta cleavage pathway of aromatic compounds",
    "L-methionine salvage cycle III",
    "protocatechuate degradation I (meta-cleavage pathway)",
    
    "mycolyl-arabinogalactan-peptidoglycan complex biosynthesis",
    "superpathway of hexuronide and hexuronate degradation",
    "D-galacturonate degradation I",
    "superpathway of sulfur oxidation (Acidianus ambivalens)",
    "CMP-pseudaminate biosynthesis",
    
    "coenzyme M biosynthesis I",
    "coenzyme B biosynthesis",
    "tetrahydromethanopterin biosynthesis",
    "mevalonate pathway II (archaea)",
    "meta cleavage pathway of aromatic compounds",
    
    "ergothioneine biosynthesis I (bacteria)",
    "superpathway of CDP-glucose-derived O-antigen building blocks biosynthesis",
    "nitrate reduction I (denitrification)"
  ),Category = c(
    "Amino acid metabolism",
    "Xenobiotics biodegradation and metabolism",
    "Carbohydrate metabolism",
    "Xenobiotics biodegradation and metabolism",
    "Xenobiotics biodegradation and metabolism",
    
    "Metabolism of cofactors and vitamins",
    "Xenobiotics biodegradation and metabolism",
    "Xenobiotics biodegradation and metabolism",
    "Carbohydrate metabolism",
    "Xenobiotics biodegradation and metabolism",
    
    "Metabolism of cofactors and vitamins",
    "Metabolism of cofactors and vitamins",
    "Metabolism of cofactors and vitamins",
    "Energy metabolism",
    "Lipid metabolism",
    "Lipid metabolism",
    "Metabolism of cofactors and vitamins",
    
    "Nucleotide metabolism / RNA modification",
    "Carbohydrate metabolism",
    "Carbohydrate metabolism",
    "Carbohydrate metabolism",
    "Xenobiotics biodegradation and metabolism",
    
    "Metabolism of cofactors and vitamins",
    "Amino acid metabolism",
    "Xenobiotics biodegradation and metabolism",
    "Amino acid metabolism",
    "Xenobiotics biodegradation and metabolism",
    
    "Glycan biosynthesis and metabolism",
    "Carbohydrate metabolism",
    "Carbohydrate metabolism",
    "Energy metabolism",
    "Carbohydrate metabolism",
    
    "Metabolism of cofactors and vitamins",
    "Metabolism of cofactors and vitamins",
    "Metabolism of cofactors and vitamins",
    "Metabolism of terpenoids and polyketides",
    "Xenobiotics biodegradation and metabolism",
    
    "Amino acid metabolism",
    "Glycan biosynthesis and metabolism",
    "Energy metabolism"), stringsAsFactors = FALSE)

# Check
stopifnot(nrow(pathway_classification) == 40)
stopifnot(length(pathway_classification$Pathway) == length(pathway_classification$Category))
table(pathway_classification$Category)

# Convert detailed categories to plotting categories matching category_colors
pathway_classification$Category[
  pathway_classification$Pathway == "7-(3-amino-3-carboxypropyl)-wyosine biosynthesis"
] <- "Metabolism of cofactors and vitamins"

pathway_classification$Category[
  pathway_classification$Pathway == "mevalonate pathway II (archaea)"
] <- "Lipid metabolism"


# Create category color mapping
# Define colors for each category
category_colors <- c(
  "Amino acid metabolism" = "#FF6B6B",
  #"Antimicrobial resistance"="#4ECDC4",
  "Carbohydrate metabolism" = "#377eb8",
  "Energy metabolism" = "#ff7f00",
  "Glycan biosynthesis and metabolism" = "#ffff33",
  "Lipid metabolism" = "#984ea3",
  "Metabolism of cofactors and vitamins" = "#33a02c",
  "Xenobiotics biodegradation and metabolism" = "#a65628"
)

matrix_rownames <- rownames(pathway_matrix)
print("Row names of matrix:")
print(matrix_rownames)

# Check pathway names in classification data
print("Pathways in classification data:")
print(pathway_classification$Pathway)

# Create row annotation from categories
row_categories <- pathway_classification$Category
names(row_categories) <- pathway_classification$Pathway

# Check for unmatched pathway names
clean_rownames <- rownames(pathway_matrix)
mismatched_pathways <- setdiff(clean_rownames, names(row_categories))
print("Unmatched pathway names:")
print(mismatched_pathways)

# Check for encoding differences
print("Encoding mismatch examples:")
for(pathway in mismatched_pathways) {
  # Find possible matches
  possible_matches <- names(row_categories)[grepl(gsub("&alpha;", "α", gsub("&beta;", "β", pathway)), names(row_categories))]
  if(length(possible_matches) > 0) {
    print(paste("Possible match:", pathway, "->", possible_matches[1]))
  }
}

# Remap with cleaned row names
row_categories_clean <- row_categories[clean_rownames]
row_categories_clean# Check for NA values

# Create annotation dataframe
row_annotation <- data.frame(Category = row_categories_clean,
                             row.names = names(row_categories_clean))

# Sort by category
row_order <- order(row_annotation$Category)
pathway_matrix_ordered_p <- pathway_matrix[row_order, ]
row_annotation_ordered <- row_annotation[row_order, , drop = FALSE]

# Draw heatmap
headmap_pathway_p <- pheatmap(pathway_matrix_ordered_p,
                              color = colorRampPalette(c("#ffffe5", "red"))(50),
                              main = "Cluster-Specific Metabolic Pathways by Functional Category",
                              cluster_rows = TRUE,
                              cluster_cols = TRUE,
                              fontsize_row = 6,
                              fontsize_col = 10,
                              display_numbers = FALSE,
                              annotation_row = row_annotation_ordered,
                              annotation_colors = list(Category = category_colors),
                              gaps_row = cumsum(rle(row_annotation_ordered$Category)$lengths)[-length(unique(row_annotation_ordered$Category))],
                              cellwidth = 15,
                              cellheight = 8)
headmap_pathway_p
# save as PDF 8 * 3.6
# save as PDF 6 * 8




#####  Figure 5H.  Relationships between cluster-level total pathway abundance and gut microbial maturity at age 4 #####
phy_all6_species_1y_4yWeek #
path_contrib_community #

### Calculate total pathway abundance per bacterial module and maturity group
cluster_maturity_analysis <- path_contrib_community2 %>%
  group_by(SampleID, cluster, Maturity_4y) %>%
  summarise(total_function = sum(total_taxon_rel_function_abun, na.rm = TRUE),
            .groups = "drop")

# Filter significant clusters (FDR < 0.05)
significant_clusters <- cluster_stats %>% filter(p_adj < 0.05)
col_maturity <- c("High"="#E69F00","Low"="#56B4E9")

# Define plotting function
plot_pathway_abundance <- function(cluster_id) {
  cluster_data <- cluster_maturity_analysis %>% filter(cluster == cluster_id)
  ggplot(cluster_data, 
         aes(x = cluster, y = total_function, group = Maturity_4y, fill = Maturity_4y)) +
    geom_point(aes(color = Maturity_4y),
               position = position_jitterdodge(jitter.width = 0.6, dodge.width = 0.8),
               size = 3, alpha = 0.4) +
    geom_boxplot(alpha = 0.1) +
       scale_y_log10()+
    labs(x = "Bacterial Cluster", y = "Log10 (relative pathway abundance)") +
    stat_compare_means() +
    scale_fill_manual(values = col_maturity) +
    scale_color_manual(values = col_maturity) +
    theme_bw()
}

# # Plot each cluster
C1_pathway_abundance_p <- plot_pathway_abundance("C1")
C2_pathway_abundance_p <- plot_pathway_abundance("C2")
C3_pathway_abundance_p <- plot_pathway_abundance("C3")
C4_pathway_abundance_p <- plot_pathway_abundance("C4")

# Statistical comparison between maturity groups
cluster_stats <- cluster_maturity_analysis %>%
  group_by(cluster) %>%
  summarise(
    p_value = wilcox.test(total_function ~ Maturity_4y)$p.value,
    mean_high = mean(total_function[Maturity_4y == "High"], na.rm = TRUE),
    mean_low = mean(total_function[Maturity_4y == "Low"], na.rm = TRUE),
    fold_change = mean_high / mean_low) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))
cluster_stats

### arrange 9 plots together
C1_pathway_abundance_p# p= 0.0206
C2_pathway_abundance_p# p= 0.0206
C3_pathway_abundance_p# p= 0.0206
C4_pathway_abundance_p# p= 0.0177

# Arrange significant clusters (C6, C7, C8)
pathway_abundance_sig <- ggarrange(C1_pathway_abundance_p,
                                   C2_pathway_abundance_p+rremove("ylab"),
                                   C3_pathway_abundance_p+rremove("ylab"),
                                   C4_pathway_abundance_p+rremove("ylab"),
                                   ncol=4,nrow=1,labels=NULL,
                                   common.legend = TRUE, legend="bottom")
pathway_abundance_sig
# save as 8.2 * 3.2


