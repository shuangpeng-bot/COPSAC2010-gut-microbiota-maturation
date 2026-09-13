########################### Gut microbial maturation_Result 3 ################################################

#### Figure 3. Characterization of the gut microbiota in the high and low maturity groups across ages ####

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


# Load phyloseq object
load("phy_feces_clean.RData")
phy_feces_clean#


#####  Figure 3A. Phylogenetic diversity of the gut microbiota in each group over time #####
phy_feces_clean#

### Calculate Faith's Phylogenetic Diversity
alpha_fpd <- pd(t(otu_table(phy_feces_clean)), phy_tree(phy_feces_clean))
alpha_fpd <- merge(data.frame(sample_data(phy_feces_clean)),data.frame(alpha_fpd),by="row.names")
alpha_fpd$Time_new <- factor(alpha_fpd$Time_new,levels = c("1w","1m","1y","4y","5y","6y"))

# Define comparisons and colors
my_comparisons <- list(c("1w","1m"),c("1m","1y"),c("1y","4y"),c("4y","5y"),c("5y","6y"))
col_maturity <- c("High"="#ff4040", "Low"="#56b4e9")

# Create phylogenetic diversity plot
Gut_maturation_all6_two_Week_phd <- alpha_fpd%>%
  subset(!is.na(Gut_maturation_Week_two))%>%
  group_by(Gut_maturation_Week_two)%>%
  
  ggplot(aes(x=Gut_maturation_Week_two,y=PD))+
  theme_bw()+theme(panel.grid = element_blank())+
  geom_jitter(aes(color=Gut_maturation_Week_two),
              alpha=0.3,size=3)+
  geom_boxplot(aes(fill=Gut_maturation_Week_two),alpha=0.7,
               outlier.colour="#bcbbbb",outlier.alpha = 0.2)+
  xlab("Microbial maturity")+
  ylab("Phylogenetic diversity")+
  scale_fill_manual(values=col_maturity)+
  scale_color_manual(values=col_maturity)+
  theme(text = element_text(size=11))+
  guides(fill=guide_legend(title="Maturity"),color=guide_legend(title="Maturity"))+
  facet_wrap(~Time_new)+
  stat_compare_means(label.y = 24.5
  )+# defalt: method="wilcox.test",paired=FALSE
  theme(legend.position="top")
Gut_maturation_all6_two_Week_phd




#####  Figure 3B. NMDS plot based on Bray-Curtis distance, showing compositional separation between groups #####
phy_feces_clean#

#### Bray-Curtis distance analysis for each timepoint

## Week 1 analysis
col_maturity <- c("High"="#ff4040", "Low"="#56b4e9")
phy_1w_all6_Week_twoall <- subset_samples(phy_feces_clean,!is.na(Gut_maturation_Week_two)&Time_new=="1w")

# PERMANOVA test (adonis)
bray_phy_1w_all6_Week_twoall <- distance(phy_1w_all6_Week_twoall,method = "bray",type="samples")
set.seed(33);
adonis_bray_phy_1w_all6_Week_twoall <- adonis2(bray_phy_1w_all6_Week_twoall~Gut_maturation_Week_two,
                                               data = data.frame(sample_data(phy_1w_all6_Week_twoall)),
                                               na.action = na.omit,permutations = 999)
adonis_bray_phy_1w_all6_Week_twoall# pr=0.001,R2=0.00667

# NMDS plot
ord_bray_phy_1w_all6_Week_twoall <- ordinate(phy_1w_all6_Week_twoall,"NMDS","bray")
plot_bray_1w_all6_Week_twoall <- plot_ordination(phy_1w_all6_Week_twoall,ord_bray_phy_1w_all6_Week_twoall,
                                                 type = "sample", color = "Gut_maturation_Week_two")+ # plot_bray_phy
  stat_ellipse(size=0.8)+ 
  guides(fill=guide_legend(title="Maturity"),color=guide_legend(title="Maturity"))+
  scale_fill_manual(values =col_maturity)+
  scale_color_manual(values =col_maturity)+
  facet_wrap(~Time_new)+
  theme_bw()+theme(panel.grid = element_blank())+# general ordination plotter based on ggplot2
  theme(text = element_text(size = 11))

# Add statistical annotation
plot_bray_1w_all6_Week_twoall <- plot_bray_1w_all6_Week_twoall +
  annotate(geom="text", x=-0.2, y=0.4,size=3, color="black",
           label="italic(P)==0.001~R^2==0.0067", parse=TRUE) # R²
plot_bray_1w_all6_Week_twoall


## Month 1 analysis
phy_1m_all6_Week_twoall <- subset_samples(phy_feces_clean,!is.na(Gut_maturation_Week_two)&Time_new=="1m")

# PERMANOVA test (adonis test)
bray_phy_1m_all6_Week_twoall <- distance(phy_1m_all6_Week_twoall,method = "bray",type="samples")
set.seed(33);
adonis_bray_phy_1m_all6_Week_twoall <- adonis2(bray_phy_1m_all6_Week_twoall~Gut_maturation_Week_two,
                                               data = data.frame(sample_data(phy_1m_all6_Week_twoall)),
                                               na.action = na.omit,permutations = 999)
adonis_bray_phy_1m_all6_Week_twoall# pr=0.001,R2=0.00656

# NMDS plot
ord_bray_phy_1m_all6_Week_twoall <- ordinate(phy_1m_all6_Week_twoall,"NMDS","bray")
plot_bray_1m_all6_Week_twoall <- plot_ordination(phy_1m_all6_Week_twoall,ord_bray_phy_1m_all6_Week_twoall,
                                                 type = "sample",color = "Gut_maturation_Week_two")+ # plot_bray_phy
  stat_ellipse(size=0.8)+ 
  guides(fill=guide_legend(title="Maturity"),color=guide_legend(title="Maturity"))+
  scale_fill_manual(values =col_maturity)+
  scale_color_manual(values =col_maturity)+
  facet_wrap(~Time_new)+
  theme_bw()+theme(panel.grid = element_blank())+# general ordination plotter based on ggplot2
  theme(text = element_text(size = 11))

# Add statistical annotation
plot_bray_1m_all6_Week_twoall <- plot_bray_1m_all6_Week_twoall + 
  annotate(geom="text", x=-0.152, y=0.29,size=3, color="black",
           label="italic(P)==0.001~R^2==0.0066", parse=TRUE) # R²
plot_bray_1m_all6_Week_twoall


## Year 1 analysis
phy_1y_all6_Week_twoall <- subset_samples(phy_feces_clean,!is.na(Gut_maturation_Week_two)&Time_new=="1y")

# PERMANOVA test (adonis test)
bray_phy_1y_all6_Week_twoall <- distance(phy_1y_all6_Week_twoall,method = "bray",type="samples")
set.seed(33);
adonis_bray_phy_1y_all6_Week_twoall <- adonis2(bray_phy_1y_all6_Week_twoall~Gut_maturation_Week_two,
                                               data = data.frame(sample_data(phy_1y_all6_Week_twoall)),
                                               na.action = na.omit,permutations = 999)
adonis_bray_phy_1y_all6_Week_twoall# pr=0.001,R2=0.01269

# NMDS plot
ord_bray_phy_1y_all6_Week_twoall <- ordinate(phy_1y_all6_Week_twoall,"NMDS","bray")
plot_bray_1y_all6_Week_twoall <- plot_ordination(phy_1y_all6_Week_twoall,ord_bray_phy_1y_all6_Week_twoall,
                                                 type = "sample",color = "Gut_maturation_Week_two")+ # plot_bray_phy
  stat_ellipse(size=0.8)+ 
  guides(fill=guide_legend(title="Maturity"),color=guide_legend(title="Maturity"))+
  scale_fill_manual(values =col_maturity)+
  scale_color_manual(values =col_maturity)+
  facet_wrap(~Time_new)+
  theme_bw()+theme(panel.grid = element_blank())+# general ordination plotter based on ggplot2
  theme(text = element_text(size = 11))

# Add statistical annotation
plot_bray_1y_all6_Week_twoall <- plot_bray_1y_all6_Week_twoall + 
  annotate(geom="text", x=-0.076, y=0.24,size=3, color="black",
           label="italic(P)==0.001~R^2==0.0127", parse=TRUE) # R²
plot_bray_1y_all6_Week_twoall


## Year 4 analysis
phy_4y_all6_Week_twoall <- subset_samples(phy_feces_clean,!is.na(Gut_maturation_Week_two)&Time_new=="4y")

# PERMANOVA test (adonis test)
bray_phy_4y_all6_Week_twoall <- distance(phy_4y_all6_Week_twoall,method = "bray",type="samples")
set.seed(33);
adonis_bray_phy_4y_all6_Week_twoall <- adonis2(bray_phy_4y_all6_Week_twoall~Gut_maturation_Week_two,
                                               data = data.frame(sample_data(phy_4y_all6_Week_twoall)),
                                               na.action = na.omit,permutations = 999)
adonis_bray_phy_4y_all6_Week_twoall# pr=0.001,R2=0.02286

# NMDS plot
ord_bray_phy_4y_all6_Week_twoall <- ordinate(phy_4y_all6_Week_twoall,"NMDS","bray")
plot_bray_4y_all6_Week_twoall <- plot_ordination(phy_4y_all6_Week_twoall,ord_bray_phy_4y_all6_Week_twoall,
                                                 type = "sample",color = "Gut_maturation_Week_two")+ # plot_bray_phy
  stat_ellipse(size=0.8)+ 
  guides(fill=guide_legend(title="Maturity"),color=guide_legend(title="Maturity"))+
  scale_fill_manual(values =col_maturity)+
  scale_color_manual(values =col_maturity)+
  facet_wrap(~Time_new)+
  theme_bw()+theme(panel.grid = element_blank())+# general ordination plotter based on ggplot2
  theme(text = element_text(size = 11))

# Add statistical annotation
plot_bray_4y_all6_Week_twoall <- plot_bray_4y_all6_Week_twoall + 
  annotate(geom="text", x=-0.05, y=1.25,size=3, color="black",
           label="italic(P)==0.001~R^2==0.0229", parse=TRUE) # R²
plot_bray_4y_all6_Week_twoall


## Year 5 analysis
phy_5y_all6_Week_twoall <- subset_samples(phy_feces_clean,!is.na(Gut_maturation_Week_two)&Time_new=="5y")

# PERMANOVA test (adonis test)
bray_phy_5y_all6_Week_twoall <- distance(phy_5y_all6_Week_twoall,method = "bray",type="samples")
set.seed(33);
adonis_bray_phy_5y_all6_Week_twoall <- adonis2(bray_phy_5y_all6_Week_twoall~Gut_maturation_Week_two,
                                               data = data.frame(sample_data(phy_5y_all6_Week_twoall)),
                                               na.action = na.omit,permutations = 999)
adonis_bray_phy_5y_all6_Week_twoall# pr=0.001,R2=0.01412

# NMDS plot
ord_bray_phy_5y_all6_Week_twoall <- ordinate(phy_5y_all6_Week_twoall,"NMDS","bray")
plot_bray_5y_all6_Week_twoall <- plot_ordination(phy_5y_all6_Week_twoall,ord_bray_phy_5y_all6_Week_twoall,
                                                 type = "sample",color = "Gut_maturation_Week_two")+ # plot_bray_phy
  stat_ellipse(size=0.8)+
  guides(fill=guide_legend(title="Maturity"),color=guide_legend(title="Maturity"))+
  scale_fill_manual(values =col_maturity)+
  scale_color_manual(values =col_maturity)+
  facet_wrap(~Time_new)+
  theme_bw()+theme(panel.grid = element_blank())+# general ordination plotter based on ggplot2
  theme(text = element_text(size = 11))

# Add statistical annotation
plot_bray_5y_all6_Week_twoall <- plot_bray_5y_all6_Week_twoall + 
  annotate(geom="text", x=-0.264, y=0.75,size=3, color="black",
           label="italic(P)==0.001~R^2==0.0141", parse=TRUE) # R²
plot_bray_5y_all6_Week_twoall


## Year 6 analysis
phy_6y_all6_Week_twoall <- subset_samples(phy_feces_clean,!is.na(Gut_maturation_Week_two)&Time_new=="6y")

# PERMANOVA test (adonis test)
bray_phy_6y_all6_Week_twoall <- distance(phy_6y_all6_Week_twoall,method = "bray",type="samples")
set.seed(33);
adonis_bray_phy_6y_all6_Week_twoall <- adonis2(bray_phy_6y_all6_Week_twoall~Gut_maturation_Week_two,
                                               data = data.frame(sample_data(phy_6y_all6_Week_twoall)),
                                               na.action = na.omit,permutations = 999)
adonis_bray_phy_6y_all6_Week_twoall# pr=0.001,R2=0.01614

# NMDS plot
ord_bray_phy_6y_all6_Week_twoall <- ordinate(phy_6y_all6_Week_twoall,"NMDS","bray")
plot_bray_6y_all6_Week_twoall <- plot_ordination(phy_6y_all6_Week_twoall,ord_bray_phy_6y_all6_Week_twoall,
                                                 type = "sample",color = "Gut_maturation_Week_two")+ # plot_bray_phy
  stat_ellipse(size=0.8)+ 
  guides(fill=guide_legend(title="Maturity"),color=guide_legend(title="Maturity"))+
  scale_fill_manual(values =col_maturity)+
  scale_color_manual(values =col_maturity)+
  facet_wrap(~Time_new)+
  theme_bw()+theme(panel.grid = element_blank())+# general ordination plotter based on ggplot2
  theme(text = element_text(size = 11))

# Add statistical annotation
plot_bray_6y_all6_Week_twoall <- plot_bray_6y_all6_Week_twoall + 
  annotate(geom="text", x=-0.948, y=0.675,size=3, color="black",
           label="italic(P)==0.001~R^2==0.0161", parse=TRUE) # R²
plot_bray_6y_all6_Week_twoall


#### Combine all NMDS plots into one figure
plot_bray_1w_all6_Week_twoall; plot_bray_1m_all6_Week_twoall#
plot_bray_1y_all6_Week_twoall; plot_bray_4y_all6_Week_twoall#
plot_bray_5y_all6_Week_twoall; plot_bray_6y_all6_Week_twoall#

bray_all6_Week_two_alltime <- ggarrange(plot_bray_1w_all6_Week_twoall + rremove("xlab"),
                                        plot_bray_1m_all6_Week_twoall+rremove("ylab") + rremove("xlab"),
                                        plot_bray_1y_all6_Week_twoall+rremove("ylab") + rremove("xlab"),
                                        plot_bray_4y_all6_Week_twoall,
                                        plot_bray_5y_all6_Week_twoall+rremove("ylab") ,
                                        plot_bray_6y_all6_Week_twoall+rremove("ylab") ,
                                        ncol=3,nrow=2,labels=NULL,
                                        common.legend = TRUE, legend="top")
bray_all6_Week_two_alltime




#####  Figure 3C. Beta-dispersion (distances to the group centroid) across time points #####

###### a.  Beta-dispersion analysis across all timepoints (1w-6y) ######

### Week 1 analysis
phy_1w_all6_Week_twoall #

# Extract grouping information
group_1w <- as.factor(sample_data(phy_1w_all6_Week_twoall)$Gut_maturation_Week_two)
Time_1w <- as.factor(sample_data(phy_1w_all6_Week_twoall)$Time_new)

# Calculate beta-dispersion
disp_1w <- betadisper(bray_phy_1w_all6_Week_twoall, group_1w)

# Prepare data for visualization
plot_df <- data.frame(
  DistanceToCentroid = disp_1w$distances,
  Group = group_1w,
  Time_new=Time_1w
)

# Statistical test
stat_test <- compare_means(
  DistanceToCentroid ~ Group, 
  data = plot_df, 
  method = "wilcox.test")

# Calculate medians for each group
medians_df <- plot_df %>%
  group_by(Time_new, Group) %>%
  summarise(median_val = median(DistanceToCentroid, na.rm = TRUE), .groups = "drop")

# Create density plot
disp_maturationtwo_1w <- ggplot(plot_df, aes(x = DistanceToCentroid, fill = Group)) +
  geom_density(alpha = 0.4) +
  geom_vline(data = medians_df, aes(xintercept = median_val, color = Group),
             linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("High" = "#70287c", "Low" = "#56b4e9")) +
  scale_color_manual(values = c("High" = "#70287c", "Low" = "#56b4e9"), guide = "none") +
  labs(y = "Density",
       x = "Distance to Group Centroid") +
  facet_wrap(~Time_new) +
  annotate("text", x = 0.7, y = 8, 
           label = paste0("P = ", signif(stat_test$p, 3))) +
  theme_bw(base_size = 14) +
  scale_x_continuous(limits = c(0.35, 0.8), breaks = seq(0.2, 0.8, by = 0.1)) +
  scale_y_continuous(limits = c(0, 8.2), breaks = seq(0, 8, by = 2))
disp_maturationtwo_1w



### Month 1 analysis
phy_1m_all6_Week_twoall #

# Extract grouping information
group_1m <- as.factor(sample_data(phy_1m_all6_Week_twoall)$Gut_maturation_Week_two)
Time_1m <- as.factor(sample_data(phy_1m_all6_Week_twoall)$Time_new)

# Calculate beta-dispersion
disp_1m <- betadisper(bray_phy_1m_all6_Week_twoall, group_1m)

# Prepare data for visualization
plot_df <- data.frame(
  DistanceToCentroid = disp_1m$distances,
  Group = group_1m,
  Time_new=Time_1m
)

# Statistical test
stat_test <- compare_means(
  DistanceToCentroid ~ Group, 
  data = plot_df, 
  method = "wilcox.test")

# Calculate medians
medians_df <- plot_df %>%
  group_by(Time_new, Group) %>%
  summarise(median_val = median(DistanceToCentroid, na.rm = TRUE), .groups = "drop")

# Create density plot
disp_maturationtwo_1m <- ggplot(plot_df, aes(x = DistanceToCentroid, fill = Group)) +
  geom_density(alpha = 0.4) +
  geom_vline(data = medians_df, aes(xintercept = median_val, color = Group),
             linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("High" = "#70287c", "Low" = "#56b4e9")) +
  scale_color_manual(values = c("High" = "#70287c", "Low" = "#56b4e9"), guide = "none") +
  labs(y = "Density",
       x = "Distance to Group Centroid") +
  facet_wrap(~Time_new) +
  annotate("text", x = 0.7, y = 8, 
           label = paste0("P = ", signif(stat_test$p, 3))) +
  theme_bw(base_size = 14) +
  scale_x_continuous(limits = c(0.35, 0.8), breaks = seq(0.2, 0.8, by = 0.1)) +
  scale_y_continuous(limits = c(0, 8.2), breaks = seq(0, 8, by = 2))
disp_maturationtwo_1m


### Year 1 analysis
phy_1y_all6_Week_twoall #

# Extract grouping information
group_1y <- as.factor(sample_data(phy_1y_all6_Week_twoall)$Gut_maturation_Week_two)

# Calculate beta-dispersion
disp_1y <- betadisper(bray_phy_1y_all6_Week_twoall, group_1y)
Time_1y <- as.factor(sample_data(phy_1y_all6_Week_twoall)$Time_new)

# Prepare data for visualization
plot_df <- data.frame(
  DistanceToCentroid = disp_1y$distances,
  Group = group_1y,
  Time_new=Time_1y
)

# Statistical test
stat_test <- compare_means(
  DistanceToCentroid ~ Group, 
  data = plot_df, 
  method = "wilcox.test")

# Calculate medians
medians_df <- plot_df %>%
  group_by(Time_new, Group) %>%
  summarise(median_val = median(DistanceToCentroid, na.rm = TRUE), .groups = "drop")

# Create density plot
disp_maturationtwo_1y <- ggplot(plot_df, aes(x = DistanceToCentroid, fill = Group)) +
  geom_density(alpha = 0.4) +
  geom_vline(data = medians_df, aes(xintercept = median_val, color = Group),
             linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("High" = "#70287c", "Low" = "#56b4e9")) +
  scale_color_manual(values = c("High" = "#70287c", "Low" = "#56b4e9"), guide = "none") +
  labs(y = "Density",
       x = "Distance to Group Centroid") +
  facet_wrap(~Time_new) +
  annotate("text", x = 0.7, y = 8, 
           label = paste0("P = ", signif(stat_test$p, 3))) +
  theme_bw(base_size = 14) +
  scale_x_continuous(limits = c(0.35, 0.8), breaks = seq(0.2, 0.8, by = 0.1)) +
  scale_y_continuous(limits = c(0, 8.2), breaks = seq(0, 8, by = 2))
disp_maturationtwo_1y


### Year 4 analysis
phy_4y_all6_Week_twoall #

# Extract grouping information
group_4y <- as.factor(sample_data(phy_4y_all6_Week_twoall)$Gut_maturation_Week_two)
Time_4y <- as.factor(sample_data(phy_4y_all6_Week_twoall)$Time_new)

# Calculate beta-dispersion
disp_4y <- betadisper(bray_phy_4y_all6_Week_twoall, group_4y)

# Prepare data for visualization
plot_df <- data.frame(
  DistanceToCentroid = disp_4y$distances,
  Group = group_4y,
  Time_new=Time_4y
)

# Statistical test
stat_test <- compare_means(
  DistanceToCentroid ~ Group, 
  data = plot_df, 
  method = "wilcox.test")

# Calculate medians
medians_df <- plot_df %>%
  group_by(Time_new, Group) %>%
  summarise(median_val = median(DistanceToCentroid, na.rm = TRUE), .groups = "drop")

# Create density plot
disp_maturationtwo_4y <- ggplot(plot_df, aes(x = DistanceToCentroid, fill = Group)) +
  geom_density(alpha = 0.4) +
  geom_vline(data = medians_df, aes(xintercept = median_val, color = Group),
             linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("High" = "#70287c", "Low" = "#56b4e9")) +
  scale_color_manual(values = c("High" = "#70287c", "Low" = "#56b4e9"), guide = "none") +
  labs(y = "Density",
       x = "Distance to Group Centroid") +
  facet_wrap(~Time_new) +
  annotate("text", x = 0.7, y = 8, 
           label = paste0("P = ", signif(stat_test$p, 3))) +
  theme_bw(base_size = 14) +
  scale_x_continuous(limits = c(0.35, 0.8), breaks = seq(0.2, 0.8, by = 0.1)) +
  scale_y_continuous(limits = c(0, 8.2), breaks = seq(0, 8, by = 2))
disp_maturationtwo_4y


### Year 5 analysis
phy_5y_all6_Week_twoall #

# Extract grouping information
group_5y <- as.factor(sample_data(phy_5y_all6_Week_twoall)$Gut_maturation_Week_two)
Time_5y <- as.factor(sample_data(phy_5y_all6_Week_twoall)$Time_new)

# Calculate beta-dispersion
disp_5y <- betadisper(bray_phy_5y_all6_Week_twoall, group_5y)

# Prepare data for visualization
plot_df <- data.frame(
  DistanceToCentroid = disp_5y$distances,
  Group = group_5y,
  Time_new=Time_5y
)

# Statistical test
stat_test <- compare_means(
  DistanceToCentroid ~ Group, 
  data = plot_df, 
  method = "wilcox.test")

# Calculate medians
medians_df <- plot_df %>%
  group_by(Time_new, Group) %>%
  summarise(median_val = median(DistanceToCentroid, na.rm = TRUE), .groups = "drop")

# Create density plot
disp_maturationtwo_5y <- ggplot(plot_df, aes(x = DistanceToCentroid, fill = Group)) +
  geom_density(alpha = 0.4) +
  geom_vline(data = medians_df, aes(xintercept = median_val, color = Group),
             linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("High" = "#70287c", "Low" = "#56b4e9")) +
  scale_color_manual(values = c("High" = "#70287c", "Low" = "#56b4e9"), guide = "none") +
  labs(y = "Density",
       x = "Distance to Group Centroid") +
  facet_wrap(~Time_new) +
  annotate("text", x = 0.7, y = 8, 
           label = paste0("P = ", signif(stat_test$p, 3))) +
  theme_bw(base_size = 14) +
  scale_x_continuous(limits = c(0.35, 0.8), breaks = seq(0.2, 0.8, by = 0.1)) +
  scale_y_continuous(limits = c(0, 8.2), breaks = seq(0, 8, by = 2))
disp_maturationtwo_5y


### Year 6 analysis
phy_6y_all6_Week_twoall #

# Extract grouping information
group_6y <- as.factor(sample_data(phy_6y_all6_Week_twoall)$Gut_maturation_Week_two)
Time_6y <- as.factor(sample_data(phy_6y_all6_Week_twoall)$Time_new)

# Calculate beta-dispersion
disp_6y <- betadisper(bray_phy_6y_all6_Week_twoall, group_6y)

# Prepare data for visualization
plot_df <- data.frame(
  DistanceToCentroid = disp_6y$distances,
  Group = group_6y,
  Time_new=Time_6y
)

# Statistical test
stat_test <- compare_means(
  DistanceToCentroid ~ Group, 
  data = plot_df, 
  method = "wilcox.test")

# Calculate medians
medians_df <- plot_df %>%
  group_by(Time_new, Group) %>%
  summarise(median_val = median(DistanceToCentroid, na.rm = TRUE), .groups = "drop")

# Create density plot
disp_maturationtwo_6y <- ggplot(plot_df, aes(x = DistanceToCentroid, fill = Group)) +
  geom_density(alpha = 0.4) +
  geom_vline(data = medians_df, aes(xintercept = median_val, color = Group),
             linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("High" = "#70287c", "Low" = "#56b4e9")) +
  scale_color_manual(values = c("High" = "#70287c", "Low" = "#56b4e9"), guide = "none") +
  labs(y = "Density",
       x = "Distance to Group Centroid") +
  facet_wrap(~Time_new) +
  annotate("text", x = 0.7, y = 8, 
           label = paste0("P = ", signif(stat_test$p, 3))) +
  theme_bw(base_size = 14) +
  scale_x_continuous(limits = c(0.35, 0.8), breaks = seq(0.2, 0.8, by = 0.1)) +
  scale_y_continuous(limits = c(0, 8.2), breaks = seq(0, 8, by = 2))
disp_maturationtwo_6y


### Combine all beta-dispersion plots
disp_maturationtwo_1w; disp_maturationtwo_1m; disp_maturationtwo_1y; 
disp_maturationtwo_4y; disp_maturationtwo_5y; disp_maturationtwo_6y; 

bdisp_maturationtwo_all6_alltime <- ggarrange(disp_maturationtwo_1w ,
                                              disp_maturationtwo_1m ,
                                              disp_maturationtwo_1y,
                                              disp_maturationtwo_4y,
                                              disp_maturationtwo_5y,
                                              disp_maturationtwo_6y,
                                              ncol=6,nrow=1,labels=NULL,
                                              common.legend = TRUE, legend="top")
bdisp_maturationtwo_all6_alltime




###### b. beta diversity partitioning /Ecological partitioning of beta diversity ######

# Define timepoints and initialize results list
time_points <- c("1w", "1m", "1y", "4y", "5y", "6y")
results_list <- list()

# Calculate beta diversity partitioning for each timepoint
for (tp in time_points) {
  
  cat("Processing:", tp, "\n")
  
  # Subset data for current timepoint
  phy_sub <- subset_samples(phy_feces_clean, Time_new == tp)
  phy_sub <- prune_taxa(taxa_sums(phy_sub) > 0, phy_sub)
  
  # Convert to presence-absence matrix
  otu_pa <- otu_table(phy_sub)
  if (taxa_are_rows(otu_pa)) { otu_pa <- t(otu_pa) }
  otu_pa[otu_pa > 0] <- 1
  otu_pa_df <- as.data.frame(otu_pa)
  
  # Separate by maturity groups
  group <- sample_data(phy_sub)$Gut_maturation_Week_two
  otu_high <- otu_pa_df[group == "High", ]
  otu_low  <- otu_pa_df[group == "Low", ]
  
  # Calculate beta diversity components
  betahigh <- beta.multi(otu_high, index.family = "jaccard")
  betalow  <- beta.multi(otu_low,  index.family = "jaccard")
  
  # Store results
  results_list[[tp]] <- data.frame(
    Time = tp,
    Group = c("High", "Low"),
    Turnover = c(betahigh$beta.JTU, betalow$beta.JTU),
    Nestedness = c(betahigh$beta.JNE, betalow$beta.JNE)
  )
}

# Combine all results and calculate nestedness percentage
beta_all_df <- bind_rows(results_list)
beta_all_df <- beta_all_df %>%
  mutate(Nestedness_pct = Nestedness / (Turnover + Nestedness) * 100)
beta_all_df

# Export results to Excel
write.xlsx(beta_all_df, file = "beta_nestedness_contribution_all_time.xlsx", rowNames = FALSE)





#####  Figure 3D. Microbial association networks inferred by SPIEC-EASI for both groups at six time points #####
phy_feces_clean#

# Define timepoints for analysis
time_points <- c("1w", "1m", "1y", "4y", "5y", "6y")

# Initialize list to store network statistics
network_stats_all <- list()

# Define data filtering and normalization function
filter_fun <- function(phy) {
  phy <- prune_taxa(taxa_sums(phy) > 10, phy)
  phy <- transform_sample_counts(phy, function(x) x / sum(x))
  return(phy)
}

# Define network metrics calculation function
network_stats <- function(net){
  list(
    nodes = vcount(net),
    edges = ecount(net),
    avg_degree = mean(degree(net)),
    clustering = transitivity(net, type="average"),
    modularity = modularity(cluster_fast_greedy(net))
  )
}

# Main analysis loop for each timepoint
for (tp in time_points) {
  message("Processing time point: ", tp)
  
  # Subset data by timepoint and maturity groups
  phy_sub <- subset_samples(phy_feces_clean, Time_new == tp & !is.na(Gut_maturation_Week_two))
  phy_high <- subset_samples(phy_sub, Gut_maturation_Week_two == "High")
  phy_low  <- subset_samples(phy_sub, Gut_maturation_Week_two == "Low")
  
  # Apply filtering and normalization
  phy_high <- filter_fun(phy_high)
  phy_low  <- filter_fun(phy_low)
  
  # Construct SPIEC-EASI networks
  spiec_high <- spiec.easi(phy_high, method='mb', lambda.min.ratio=1e-2, nlambda=20,
                           pulsar.params = list(thresh=0.05))
  spiec_low  <- spiec.easi(phy_low, method='mb', lambda.min.ratio=1e-2, nlambda=20,
                           pulsar.params = list(thresh=0.05))
  
  # Save SPIEC-EASI results
  #save(spiec_high, file = paste0("spiec_high_", tp, ".RData"))
  #save(spiec_low,  file = paste0("spiec_low_", tp,  ".RData"))
  
  # Convert to igraph networks
  net_high <- graph_from_adjacency_matrix(as(getRefit(spiec_high), "dgCMatrix"),
                                          mode = "undirected", diag = FALSE)
  net_low  <- graph_from_adjacency_matrix(as(getRefit(spiec_low), "dgCMatrix"),
                                          mode = "undirected", diag = FALSE)
  
  # Export networks for Cytoscape and Gephi
  write_graph(net_high, file = paste0("net_high_", tp, ".gml"), format = "gml")
  write_graph(net_low,  file = paste0("net_low_", tp,  ".gml"), format = "gml")
  
  write_graph(net_high, file = paste0("net_high_", tp, ".graphml"), format = "graphml")
  write_graph(net_low,  file = paste0("net_low_", tp,  ".graphml"), format = "graphml")
  
  # Calculate and store network statistics
  network_stats_all[[tp]] <- list(
    High = network_stats(net_high),
    Low  = network_stats(net_low)
  )
}

# Export network statistics to CSV
stats_df <- do.call(rbind, lapply(names(network_stats_all), function(tp) {
  data.frame(
    Time = tp,
    Group = c("High", "Low"),
    Nodes = c(network_stats_all[[tp]]$High$nodes, network_stats_all[[tp]]$Low$nodes),
    Edges = c(network_stats_all[[tp]]$High$edges, network_stats_all[[tp]]$Low$edges),
    AvgDegree = c(network_stats_all[[tp]]$High$avg_degree, network_stats_all[[tp]]$Low$avg_degree),
    Clustering = c(network_stats_all[[tp]]$High$clustering, network_stats_all[[tp]]$Low$clustering),
    Modularity = c(network_stats_all[[tp]]$High$modularity, network_stats_all[[tp]]$Low$modularity)
  )
}))
stats_df
write.csv(stats_df, "network_stats_summary.csv", row.names = FALSE)




#####  Figure 3E. Node degree distribution of the SPIEC-EASI networks across ages #####

# Load network files
net_high_1w <- read_graph("net_high_1w.gml", format = "gml")
net_high_1m <- read_graph("net_high_1m.gml", format = "gml")
net_high_1y <- read_graph("net_high_1y.gml", format = "gml")
net_high_4y <- read_graph("net_high_4y.gml", format = "gml")
net_high_5y <- read_graph("net_high_5y.gml", format = "gml")
net_high_6y <- read_graph("net_high_6y.gml", format = "gml")

net_low_1w <- read_graph("net_low_1w.gml", format = "gml")
net_low_1m <- read_graph("net_low_1m.gml", format = "gml")
net_low_1y <- read_graph("net_low_1y.gml", format = "gml")
net_low_4y <- read_graph("net_low_4y.gml", format = "gml")
net_low_5y <- read_graph("net_low_5y.gml", format = "gml")
net_low_6y <- read_graph("net_low_6y.gml", format = "gml")

# Organize networks in a list for easier management
net_list <- list(
  "1w_High" = net_high_1w,
  "1m_High" = net_high_1m,
  "1y_High" = net_high_1y,
  "4y_High" = net_high_4y,
  "5y_High" = net_high_5y,
  "6y_High" = net_high_6y,
  "1w_Low"  = net_low_1w,
  "1m_Low"  = net_low_1m,
  "1y_Low"  = net_low_1y,
  "4y_Low"  = net_low_4y,
  "5y_Low"  = net_low_5y,
  "6y_Low"  = net_low_6y
)

# Calculate degree distributions and create consolidated data frame
degree_df <- do.call(rbind, lapply(names(net_list), function(name) {
  deg <- degree(net_list[[name]])
  data.frame(
    Time = gsub("_.*", "", name),
    Group = gsub(".*_", "", name),
    Degree = deg
  )
}))

# Calculate degree frequency distribution
deg_freq_df <- degree_df %>%
  group_by(Time, Group, Degree) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  group_by(Time, Group) %>%
  mutate(Frequency = Freq / sum(Freq))

#Enhanced degree distribution plot with median annotations
  # Prepare data for median annotations
  label_df <- deg_freq_df %>%
    group_by(Time, Group) %>%
    summarise(
      median_deg = median(Degree, na.rm = TRUE),
      max_y = max(Frequency, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      x = median_deg / 10,
      y = max_y + 0.01,  # 稍微上移一点，防止重叠
      label = paste0("Median: ", round(median_deg, 1))
    )
  
  # Create enhanced plot with median annotations
  deg_freq_df$Time <- factor(deg_freq_df$Time, levels = c("1w","1m","1y","4y","5y","6y"))
  
  p_deg_freq_SPIEC_EASI_1w1m1y4y5y6y <- deg_freq_df %>%
    ggplot(aes(x = Degree, y = Frequency, color = Group)) +
    geom_line(size = 1) +
    geom_text(data = label_df, aes(x = x, y = y, label = label, color = Group),
              size = 3.2, show.legend = FALSE, inherit.aes = FALSE) +
    facet_wrap(~ Time, ncol = 6) +
    theme_minimal(base_size = 14) +
    labs(x = "Node degree of networks",
         y = "Frequency") +
    scale_color_manual(values = c("High" = "#70287c", "Low" = "#56b4e9")) +
    scale_x_continuous(limits = c(0, 120)) +
    theme(legend.position = "top") +
    guides(fill = guide_legend(title = "Maturity"), color = guide_legend(title = "Maturity"))
  p_deg_freq_SPIEC_EASI_1w1m1y4y5y6y


  
  
  
  
  ##### Supplementary Fig. S9 sensitivity analysis: equal-size subsampled SPIEC-EASI networks #####
  
  colnames(sample_data(phy_feces_clean))
  #save(phy_feces_clean,file="phy_feces_clean2.RData")
  
  
  ##### Server parallel version
  
  #####  Set working directory
  setwd("/mnt/data/ShuangPeng/Project/COPSAC2010/high_low_network/")
  
  ##### 1. Load required packages
  
  library(phyloseq)
  library(SpiecEasi)
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(parallel)
  
  ##### 2. Avoid nested over-threading
  Sys.setenv(OMP_NUM_THREADS = "1",
             OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1",
             VECLIB_MAXIMUM_THREADS = "1")
  
  ##### 3. Load phyloseq object
  load("phy_feces_clean2.RData")
  phy_feces_clean
  
  colnames(sample_data(phy_feces_clean))
  
  ##### 4. Set parameters
  time_points <- c("1w", "1m", "1y", "4y", "5y", "6y")
  n_iter <- 100
  n_cores <- 120 # 128
  
  
  ##### 5. Create output folders
  dir.create("SPIEC_EASI_equal_size_results", showWarnings = FALSE)
  dir.create("SPIEC_EASI_equal_size_results/job_rds", showWarnings = FALSE)
  dir.create("SPIEC_EASI_equal_size_results/tables", showWarnings = FALSE)
  dir.create("SPIEC_EASI_equal_size_results/plots", showWarnings = FALSE)
  
  ##### 6. Define filtering and normalization function
  filter_fun <- function(phy) {
    phy <- prune_samples(sample_sums(phy) > 0, phy)
    phy <- prune_taxa(taxa_sums(phy) > 10, phy)
    phy <- transform_sample_counts(phy, function(x) x / sum(x))
    return(phy)
  }
  
  ##### 7. Define network metrics function
  network_stats <- function(net) {
    data.frame(
      Nodes = vcount(net),
      Edges = ecount(net),
      AvgDegree = mean(degree(net)),
      MedianDegree = median(degree(net)),
      Clustering = ifelse(ecount(net) > 0, transitivity(net, type = "average"), NA),
      Modularity = ifelse(ecount(net) > 0, modularity(cluster_fast_greedy(net)), NA)
    )
  }
  
  ##### 8. Function to construct SPIEC-EASI network
  run_spiec_network <- function(phy_obj) {
    
    phy_obj <- filter_fun(phy_obj)
    
    if (nsamples(phy_obj) < 10 | ntaxa(phy_obj) < 5) {
      return(NULL)
    }
    
    spiec_obj <- tryCatch(
      spiec.easi(
        phy_obj,
        method = "mb",
        lambda.min.ratio = 1e-2,
        nlambda = 20,
        pulsar.params = list(
          thresh = 0.05,
          ncores = 1
        )
      ),
      error = function(e) {
        message("SPIEC-EASI failed: ", e$message)
        return(NULL)
      }
    )
    
    if (is.null(spiec_obj)) {
      return(NULL)
    }
    
    net <- graph_from_adjacency_matrix(
      as(getRefit(spiec_obj), "dgCMatrix"),
      mode = "undirected",
      diag = FALSE
    )
    
    return(net)
  }
  
  ##### 9. Prepare data for each time point
  
  phy_by_time <- list()
  
  for (tp in time_points) {
    
    message("Preparing time point: ", tp)
    
    phy_sub <- subset_samples(
      phy_feces_clean,
      Time_new == tp & !is.na(Gut_maturation_Week_two)
    )
    
    phy_high_full <- subset_samples(phy_sub, Gut_maturation_Week_two == "High")
    phy_low_full  <- subset_samples(phy_sub, Gut_maturation_Week_two == "Low")
    
    n_high <- nsamples(phy_high_full)
    n_low  <- nsamples(phy_low_full)
    n_equal <- min(n_high, n_low)
    
    message("High N = ", n_high,
            "; Low N = ", n_low,
            "; equal subsampled N = ", n_equal)
    
    phy_by_time[[tp]] <- list(
      High = phy_high_full,
      Low = phy_low_full,
      n_high = n_high,
      n_low = n_low,
      n_equal = n_equal
    )
  }
  # Preparing time point: 1w  High N = 276; Low N = 277; equal subsampled N = 276
  # Preparing time point: 1m  High N = 305; Low N = 309; equal subsampled N = 305
  # Preparing time point: 1y  High N = 303; Low N = 324; equal subsampled N = 303
  # Preparing time point: 4y  High N = 196; Low N = 219; equal subsampled N = 196
  # Preparing time point: 5y  High N = 101; Low N = 124; equal subsampled N = 101
  # Preparing time point: 6y  High N = 119; Low N = 136; equal subsampled N = 119
  
  
  
  
  ##### 10. Create all parallel jobs
  job_df <- expand.grid(
    Time = time_points,
    Iteration = seq_len(n_iter),
    stringsAsFactors = FALSE
  )
  
  job_df$JobID <- seq_len(nrow(job_df))
  
  head(job_df)
  dim(job_df)
  
  # n_iter = 100 ，in total 6 × 100 = 600  jobs, each job will create High and Low networks
  # 1200  SPIEC-EASI networks in total
  
  ##### 11. Define one-job function
  
  run_one_job <- function(job_id) {
    
    job <- job_df[job_id, ]
    
    tp <- job$Time
    i  <- job$Iteration
    
    message("Running job ", job_id, "/", nrow(job_df),
            " | Time = ", tp,
            " | Iteration = ", i)
    
    out_rds <- file.path(
      "SPIEC_EASI_equal_size_results/job_rds",
      paste0("job_", tp, "_iter_", i, ".rds")
    )
    
    if (file.exists(out_rds)) {
      return(readRDS(out_rds))
    }
    
    result <- tryCatch({
      
      phy_high_full <- phy_by_time[[tp]]$High
      phy_low_full  <- phy_by_time[[tp]]$Low
      
      n_high  <- phy_by_time[[tp]]$n_high
      n_low   <- phy_by_time[[tp]]$n_low
      n_equal <- phy_by_time[[tp]]$n_equal
      
      # set seed
      set.seed(100000 + match(tp, time_points) * 1000 + i)
      
      # Equal-size random subsampling
      high_ids <- sample(sample_names(phy_high_full), n_equal, replace = FALSE)
      low_ids  <- sample(sample_names(phy_low_full),  n_equal, replace = FALSE)
      
      phy_high_sub <- prune_samples(high_ids, phy_high_full)
      phy_low_sub  <- prune_samples(low_ids,  phy_low_full)
      
      # Construct SPIEC-EASI networks
      net_high <- run_spiec_network(phy_high_sub)
      net_low  <- run_spiec_network(phy_low_sub)
      
      ##### High group results:
      
      if (!is.null(net_high)) {
        
        high_stats <- network_stats(net_high) %>%
          mutate(
            Time = tp,
            Iteration = i,
            Group = "High",
            Full_N_High = n_high,
            Full_N_Low = n_low,
            Subsampled_N = n_equal,
            Error = NA_character_
          )
        
        high_degree <- data.frame(
          Time = tp,
          Iteration = i,
          Group = "High",
          Degree = degree(net_high),
          Subsampled_N = n_equal,
          Error = NA_character_
        )
        
      } else {
        
        high_stats <- data.frame(
          Nodes = NA,
          Edges = NA,
          AvgDegree = NA,
          MedianDegree = NA,
          Clustering = NA,
          Modularity = NA,
          Time = tp,
          Iteration = i,
          Group = "High",
          Full_N_High = n_high,
          Full_N_Low = n_low,
          Subsampled_N = n_equal,
          Error = "SPIEC-EASI returned NULL"
        )
        
        high_degree <- data.frame()
      }
      
      ##### Low group results:
      
      if (!is.null(net_low)) {
        
        low_stats <- network_stats(net_low) %>%
          mutate(
            Time = tp,
            Iteration = i,
            Group = "Low",
            Full_N_High = n_high,
            Full_N_Low = n_low,
            Subsampled_N = n_equal,
            Error = NA_character_
          )
        
        low_degree <- data.frame(
          Time = tp,
          Iteration = i,
          Group = "Low",
          Degree = degree(net_low),
          Subsampled_N = n_equal,
          Error = NA_character_
        )
        
      } else {
        
        low_stats <- data.frame(
          Nodes = NA,
          Edges = NA,
          AvgDegree = NA,
          MedianDegree = NA,
          Clustering = NA,
          Modularity = NA,
          Time = tp,
          Iteration = i,
          Group = "Low",
          Full_N_High = n_high,
          Full_N_Low = n_low,
          Subsampled_N = n_equal,
          Error = "SPIEC-EASI returned NULL"
        )
        
        low_degree <- data.frame()
      }
      
      list(
        metrics = bind_rows(high_stats, low_stats),
        degree = bind_rows(high_degree, low_degree)
      )
      
    }, error = function(e) {
      
      error_metrics <- data.frame(
        Nodes = NA,
        Edges = NA,
        AvgDegree = NA,
        MedianDegree = NA,
        Clustering = NA,
        Modularity = NA,
        Time = tp,
        Iteration = i,
        Group = NA_character_,
        Full_N_High = NA,
        Full_N_Low = NA,
        Subsampled_N = NA,
        Error = e$message
      )
      
      list(
        metrics = error_metrics,
        degree = data.frame()
      )
    })
    
    saveRDS(result, out_rds)
    
    return(result)
  }
  
  ##### 12. Run jobs in parallel
  
  message("Start parallel SPIEC-EASI subsampling analysis")
  message("Number of jobs: ", nrow(job_df))
  message("Number of cores: ", n_cores)
  
  res_list <- parallel::mclapply(
    X = seq_len(nrow(job_df)),
    FUN = run_one_job,
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
  
  message("Parallel SPIEC-EASI subsampling analysis finished")
  
  
  
  ### Check error jobs
  error_summary <- subsample_network_df %>%
    group_by(Time, Group, Error) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(Time, Group, desc(n))
  
  error_summary
  
  
  ##### 13. Combine all results
  
  ## On R windows
  load("Linux_network_high_low.RData")
  
  subsample_network_df <- bind_rows(lapply(res_list, function(x) x$metrics))
  subsample_degree_df  <- bind_rows(lapply(res_list, function(x) x$degree))
  
  subsample_network_df$Time <- factor(
    subsample_network_df$Time,
    levels = c("1w", "1m", "1y", "4y", "5y", "6y")
  )
  
  if (nrow(subsample_degree_df) > 0) {
    subsample_degree_df$Time <- factor(
      subsample_degree_df$Time,
      levels = c("1w", "1m", "1y", "4y", "5y", "6y")
    )
  }
  
  subsample_network_df
  subsample_degree_df
  
  View(subsample_network_df)
  View(subsample_degree_df)
  
  ##### 14. Save raw results
  
  # write.xlsx(subsample_network_df,file.path("SPIEC_EASI_equal_size_results/tables/SPIEC_EASI_equal_size_subsampling_network_metrics_all_iterations.xlsx"),rowNames = FALSE)
  # write.xlsx(subsample_degree_df, file.path("SPIEC_EASI_equal_size_results/tables/SPIEC_EASI_equal_size_subsampling_degree_all_iterations.xlsx"),rowNames = FALSE)

  
  ##### 15. Summary table
  
  subsample_summary <- subsample_network_df %>%
    filter(is.na(Error)) %>%
    group_by(Time, Group) %>%
    summarise(
      N_iterations = n(),
      Subsampled_N = unique(Subsampled_N)[1],
      
      Nodes_mean = mean(Nodes, na.rm = TRUE),
      Nodes_sd = sd(Nodes, na.rm = TRUE),
      
      Edges_mean = mean(Edges, na.rm = TRUE),
      Edges_sd = sd(Edges, na.rm = TRUE),
      
      AvgDegree_mean = mean(AvgDegree, na.rm = TRUE),
      AvgDegree_sd = sd(AvgDegree, na.rm = TRUE),
      
      MedianDegree_mean = mean(MedianDegree, na.rm = TRUE),
      MedianDegree_sd = sd(MedianDegree, na.rm = TRUE),
      
      Clustering_mean = mean(Clustering, na.rm = TRUE),
      Clustering_sd = sd(Clustering, na.rm = TRUE),
      
      Modularity_mean = mean(Modularity, na.rm = TRUE),
      Modularity_sd = sd(Modularity, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  subsample_summary
  View(subsample_summary)
  
  write.xlsx(subsample_summary,
             file.path("D:/OneDrive - University of Copenhagen/01Project/COPSAC/COPSAC2010_feces/copsac2010_feces_species/SPIEC_EASI_equal_size_subsampling_summary_by_group.xlsx"),
             rowNames = FALSE)
  
  

  ##### 16. Paired High-vs-Low comparison
  
  wilcox_p <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2 || length(unique(x)) < 2) {
      return(NA_real_)
    } else {
      return(wilcox.test(x, mu = 0)$p.value)
    }
  }
  
  subsample_wide <- subsample_network_df %>%
    filter(is.na(Error)) %>%
    select(Time, Iteration, Group, Nodes, Edges, AvgDegree, MedianDegree, Clustering, Modularity) %>%
    pivot_wider(
      names_from = Group,
      values_from = c(Nodes, Edges, AvgDegree, MedianDegree, Clustering, Modularity)
    ) %>%
    mutate(
      Nodes_diff = Nodes_High - Nodes_Low,
      Edges_diff = Edges_High - Edges_Low,
      AvgDegree_diff = AvgDegree_High - AvgDegree_Low,
      MedianDegree_diff = MedianDegree_High - MedianDegree_Low,
      Clustering_diff = Clustering_High - Clustering_Low,
      Modularity_diff = Modularity_High - Modularity_Low
    )
  
  subsample_comparison <- subsample_wide %>%
    group_by(Time) %>%
    summarise(
      N_iterations = n(),
      
      Nodes_diff_mean = mean(Nodes_diff, na.rm = TRUE),
      Nodes_High_greater_prop = mean(Nodes_diff > 0, na.rm = TRUE),
      Nodes_P = wilcox_p(Nodes_diff),
      
      Edges_diff_mean = mean(Edges_diff, na.rm = TRUE),
      Edges_High_greater_prop = mean(Edges_diff > 0, na.rm = TRUE),
      Edges_P = wilcox_p(Edges_diff),
      
      AvgDegree_diff_mean = mean(AvgDegree_diff, na.rm = TRUE),
      AvgDegree_High_greater_prop = mean(AvgDegree_diff > 0, na.rm = TRUE),
      AvgDegree_P = wilcox_p(AvgDegree_diff),
      
      MedianDegree_diff_mean = mean(MedianDegree_diff, na.rm = TRUE),
      MedianDegree_High_greater_prop = mean(MedianDegree_diff > 0, na.rm = TRUE),
      MedianDegree_P = wilcox_p(MedianDegree_diff),
      
      Clustering_diff_mean = mean(Clustering_diff, na.rm = TRUE),
      Clustering_High_greater_prop = mean(Clustering_diff > 0, na.rm = TRUE),
      Clustering_P = wilcox_p(Clustering_diff),
      
      Modularity_diff_mean = mean(Modularity_diff, na.rm = TRUE),
      Modularity_High_greater_prop = mean(Modularity_diff > 0, na.rm = TRUE),
      Modularity_P = wilcox_p(Modularity_diff),
      
      .groups = "drop"
    ) %>%
    mutate(
      Nodes_FDR = p.adjust(Nodes_P, method = "BH"),
      Edges_FDR = p.adjust(Edges_P, method = "BH"),
      AvgDegree_FDR = p.adjust(AvgDegree_P, method = "BH"),
      MedianDegree_FDR = p.adjust(MedianDegree_P, method = "BH"),
      Clustering_FDR = p.adjust(Clustering_P, method = "BH"),
      Modularity_FDR = p.adjust(Modularity_P, method = "BH")
    )
  
  subsample_comparison
  View(subsample_comparison)
  
  #write.xlsx(subsample_comparison,
  #           file.path("D:/OneDrive - University of Copenhagen/01Project/COPSAC/COPSAC2010_feces/copsac2010_feces_species/SPIEC_EASI_equal_size_subsampling_High_vs_Low_paired.xlsx"),rowNames = FALSE)
  
  
 
  
  ##### 18. Visualization: paired High-Low differences
  
  subsample_diff_plot_df <- subsample_wide %>%
    select(Time, Iteration,
           Nodes_diff, Edges_diff, AvgDegree_diff,
           MedianDegree_diff, Clustering_diff, Modularity_diff) %>%
    pivot_longer(cols = ends_with("_diff"),
                 names_to = "Metric",
                 values_to = "High_minus_Low") %>%
    mutate( Metric = gsub("_diff", "", Metric),
            Metric = factor(Metric,
                            levels = c("Nodes", "Edges", "AvgDegree", "MedianDegree", "Clustering", "Modularity")))
  
  col_time <- c("1w"="#FF8C94","1m"="#49d261","1y"="#EE4266","4y"="#33B5BF","5y"="#9B59B6",
                "6y"="#FFD159")
  
  p_subsample_network_diff <- ggplot(subsample_diff_plot_df,
                                     aes(x = Time, y = High_minus_Low, fill=Time,color=Time)) +
    
    geom_boxplot(outlier.shape = NA, alpha=0.8, color="black") +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    facet_wrap(~ Metric, scales = "free_y", ncol = 6) +
    scale_color_manual(values=col_time)+
    scale_fill_manual(values=col_time)+
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank()) +
    labs(x = "Time",
         y = "High minus Low")
  
  p_subsample_network_diff
  

  ##### 20. Export Excel file
  
  openxlsx::write.xlsx(
    list("All_network_metrics" = subsample_network_df,
      "Summary_by_group" = subsample_summary,
      "High_vs_Low_paired" = subsample_comparison,
      "Degree_all_iterations" = subsample_degree_df ),
    file = "D:/OneDrive - University of Copenhagen/01Project/COPSAC/COPSAC2010_feces/copsac2010_feces_species/SPIEC_EASI_equal_size_subsampling_sensitivity.xlsx",
    rowNames = FALSE)
  
  message("All done!")
  message("Results are saved in: SPIEC_EASI_equal_size_results/")
  
  

  