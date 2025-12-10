#############################################Alpha-diversity analysis#############################################
# Remove all stored objects & unload non-base packages
rm(list = ls())

detachAllPackages <- function() {
  base_pkgs <- c("stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")
  loaded_pkgs <- search()[grep("package:", search())]
  pkgs_to_detach <- setdiff(loaded_pkgs, paste("package:", base_pkgs, sep = ""))
  
  if (length(pkgs_to_detach) > 0) {
    lapply(pkgs_to_detach, function(pkg) detach(pkg, character.only = TRUE))
  }
}
detachAllPackages()
# Set working directory
setwd("#set your working directory here")
## Setup dependencies
# Load libraries
library('DT'); packageVersion("DT")
library(phyloseq); packageVersion("phyloseq")
#library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(iNEXT); packageVersion("iNEXT")
library(ape); packageVersion("ape")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(cowplot); packageVersion("cowplot")
#library(plyr); packageVersion("plyr")
#library(sjmisc); packageVersion("sjmisc")
library(data.table); packageVersion("data.table")
#library(metacoder);packageVersion("metacoder")
#library(hilldiv);packageVersion("hilldiv")
library(stringr); packageVersion("stringr")
library(ggh4x); packageVersion("ggh4x")
#library(car);packageVersion("car")
#library(lme4);packageVersion("lme4")
library(RColorBrewer);packageVersion("RColorBrewer")
library(ggpubr);packageVersion("ggpubr")
library(ggside);packageVersion("ggside")
library(ggtext);packageVersion("ggtext")
library(forcats);packageVersion("forcats")
library(microbiomeutilities);packageVersion("microbiomeutilities")

#Load data generated from 16S HiFi sequencing
OTU<-read.csv("asvtable.csv", header=TRUE, row.names=1)
TAX<-read.csv("taxonomy.csv", header=TRUE, row.names=1)
META<-read.csv("metadata.csv", header=TRUE, row.names=1)
class(OTU)
class(TAX)
OTU_data=otu_table(as.matrix(OTU), taxa_are_rows=TRUE)
TAX_data=tax_table(as.matrix(TAX))
META_data=sample_data(META)
physeq<-merge_phyloseq(phyloseq(OTU_data, TAX_data), META_data)
physeq<-t(physeq)
# identifying the ASVs with zero counts across the dataset
any(taxa_sums(physeq) == 0)
# total num of ASVs with zero counts across the dataset
sum(taxa_sums(physeq) == 0)
# number of asv features in phyloseq object
ntaxa(physeq)
# removing the ASVs with zero counts across the dataset (if required)
ps_ab = physeq
ps_ab = prune_taxa(taxa_sums(ps_ab) > 0, ps_ab)
any(sample_sums(ps_ab) == 0)
ps_ab

# Remove cyanobacteria, mitochondria and chloroplast
ps_ab <- ps_ab %>%
  subset_taxa(
    Phylum != "p__Cyanobacteria" &
      Family  != "f__mitochondria" &
      Class   != "c__Chloroplast"
  )

ps_ab
#Calculate alpha diversity matrics
hmp.div <- estimate_richness(ps_ab, split = TRUE)
datatable(hmp.div)
# get the metadata out as a separate object
hmp.meta <- meta(ps_ab)
# Add the rownames as a new column for easy integration later.
hmp.meta$sam_name <- rownames(hmp.meta)
# Add the rownames to the diversity table
hmp.div$sam_name <- rownames(hmp.div)
# merge these two data frames into one
div.df <- merge(hmp.div,hmp.meta, by = "sam_name")
# check the tables
colnames(div.df)
# convert phyloseq object into a long data format.
div.df2 <- div.df[, c("Term", "Timepoint", "Probiotic", "Observed", "Chao1", "Shannon", "Simpson")]
colnames(div.df2)
div_df_melt <- reshape2::melt(div.df2)
head(div_df_melt)
# Select each subset of variables
div_df_melt_ob <- div_df_melt %>%
  filter(variable == "Observed")
div_df_melt_chao1 <- div_df_melt %>%
  filter(variable == "Chao1")
div_df_melt_shannon <- div_df_melt %>%
  filter(variable == "Shannon")
div_df_melt_simpson <- div_df_melt %>%
  filter(variable == "Simpson")

#### Observed
#Re-order x-axis for term
div_df_melt_ob$Term <- factor(div_df_melt_ob$Term, levels=c("P ", "F "))
#Re-order x-axis for timepoint
div_df_melt_ob$Timepoint <- factor(div_df_melt_ob$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
#Re-order x-axis for probiotic
div_df_melt_ob$Probiotic <- factor(div_df_melt_ob$Probiotic, levels=c("N", "Y", "F"))
# add sample info as number of samples for observed species in "N", "Y" and "F" probiotic groups for each timepoint "1W", "4W", "4M" and "1Y"
div_df_melt_ob <- div_df_melt_ob %>% 
    group_by(Probiotic) %>% 
    mutate(n_pro= as.character(n())) %>% 
    mutate(label_pro= factor(paste0(Probiotic, '\n n = ', n_pro))) %>% 
    ungroup() %>% 
    select(-n_pro)
# Now use this data frame to plot 
p_ob <- ggboxplot(div_df_melt_ob, x = "Probiotic", y = "value",
                  color = "Probiotic",
                  fill = "Probiotic",
                  legend = "right",
                  scales = "free",
                  add = "dotplot",
                  alpha = 0.2, # Transparency of points
                  ggtheme = theme_bw(),
                  palette = c("#cd5b45", "#008a8a", "#E7B800", "#009acd")) + scale_x_discrete(labels= levels(div_df_melt_ob$label_pro)) # add sample size info in the x-axis

p_ob <- p_ob + 
    facet_wrap(~Timepoint) # Creates separate panels for each level of timepoint
p_ob
# Add p-values to the plot
p1_ob <- p_ob +
    stat_compare_means(aes(group = Probiotic),label = "p.format", method = "kruskal.test", 
                       label.y = 0.98*max(div_df_melt_ob$value)) # Add global p-value
p1_ob
# run multiple comparison wilcox test if Kruskal-Wallis test is significant and add p-values to the plot
# p1_ob_mc <- p1_ob + stat_compare_means(aes(group = Probiotic), label = "p.signif", method = "wilcox.test",
#                                  p.adjust.method = "bonferroni",  # other options: "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none"
#                                  label.y = c(0.95*max(div_df_melt_ob$value), 0.92*max(div_df_melt_ob$value), 0.89*max(div_df_melt_ob$value)),
#                                  comparisons = list(c("N", "Y"), c("N", "F"), c("Y", "F"))) # Add pairwise p-values
# p1_ob_mc
# run a multiple comparison test (Dunn's test) if the Kruskal-Wallis test is significant and add super script as "a"/"b" (not p-value) to the plot
# p1_ob_mct <- p1_ob + stat_compare_means(aes(group = Probiotic), label = after_stat(p.signif), method = "wilcox.test",
#                                  p.adjust.method = "bonferroni",  # other options: "holm", "hochberg", "hommel", "BH", "BY", "fdr", "none"
#                                  label.y = c(0.95*max(div_df_melt_ob$value), 0.92*max(div_df_melt_ob$value), 0.89*max(div_df_melt_ob$value)),
#                                  comparisons = list(c("N", "Y"), c("N", "F"), c("Y", "F"))) # Add pairwise p-values  
# p1_ob_mct
# Customize the theme
p2_ob <- p1_ob + theme_bw() +
    theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.75),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10)
    ) +
    labs(
        x = "Probiotic",
        y = "Observed Species"
    )
p2_ob
# annotate the plot with "a"/"b" for significant differences if required
p2_ob <- p2_ob + annotate("text", x=1, y=550, label= "b") + 
    annotate("text", x = 2, y=550, label = "b") +
    annotate("text", x = 3, y=550, label = "a")
p2_ob

#### Shannon
#Re-order x-axis for term
div_df_melt_shannon$Term <- factor(div_df_melt_shannon$Term, levels=c("P ", "F "))
#Re-order x-axis for timepoint
div_df_melt_shannon$Timepoint <- factor(div_df_melt_shannon$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
#Re-order x-axis for probiotic
div_df_melt_shannon$Probiotic <- factor(div_df_melt_shannon$Probiotic, levels=c("N", "Y", "F"))
# add sample info as number of samples in each group
div_df_melt_shannon <- div_df_melt_shannon %>% 
    group_by(Probiotic) %>% 
    mutate(n_pro= as.character(n())) %>% 
    mutate(label_pro= factor(paste0(Probiotic, '\n n = ', n_pro))) %>% 
    ungroup() %>% 
    select(-n_pro)
# Now use this data frame to plot
p_shannon <- ggboxplot(div_df_melt_shannon, x = "Probiotic", y = "value",
                       color = "Probiotic",
                       fill = "Probiotic",
                       legend = "right",
                       scales = "free",
                       add = "dotplot",
                       alpha = 0.2, # Transparency of points
                       ggtheme = theme_bw(),
                       palette = c("#cd5b45", "#008a8a", "#E7B800", "#009acd")) + scale_x_discrete(labels= levels(div_df_melt_shannon$label_pro)) # add sample size info in the x-axis
p_shannon <- p_shannon +
    facet_wrap(~Timepoint) # Creates separate panels for each level of timepoint
p_shannon
# Add p-values to the plot
p1_shannon <- p_shannon +
    stat_compare_means(aes(group = Probiotic),label = "p.format", method = "kruskal.test", 
                       label.y = 0.98*max(div_df_melt_shannon$value)) # Add global p-value
p1_shannon
# Customize the theme
p2_shannon <- p1_shannon + theme_bw() +
    theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.75),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10)
    ) +
    labs(
        x = "Probiotic",
        y = "Shannon diversity"
    )
p2_shannon
# annotate the plot with "a"/"b" for significant differences if required
p2_shannon <- p2_shannon + annotate("text", x=1, y=6.5, label= "b") + 
    annotate("text", x = 2, y=6.5, label = "b") +
    annotate("text", x = 3, y=6.5, label = "a")
p2_shannon

#### Simpson
#Re-order x-axis for term
div_df_melt_simpson$Term <- factor(div_df_melt_simpson$Term, levels=c("P ", "F "))
#Re-order x-axis for timepoint
div_df_melt_simpson$Timepoint <- factor(div_df_melt_simpson$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
#Re-order x-axis for probiotic
div_df_melt_simpson$Probiotic <- factor(div_df_melt_simpson$Probiotic, levels=c("N", "Y", "F"))
# add sample info as number of samples in each group
div_df_melt_simpson <- div_df_melt_simpson %>% 
    group_by(Probiotic) %>% 
    mutate(n_pro= as.character(n())) %>% 
    mutate(label_pro= factor(paste0(Probiotic, '\n n = ', n_pro))) %>% 
    ungroup() %>% 
    select(-n_pro)
# Now use this data frame to plot
p_simpson <- ggboxplot(div_df_melt_simpson, x = "Probiotic", y = "value",
                       color = "Probiotic",
                       fill = "Probiotic",
                       legend = "right",
                       scales = "free",
                       add = "dotplot",
                       alpha = 0.2, # Transparency of points
                       ggtheme = theme_bw(),
                       palette = c("#cd5b45", "#008a8a", "#E7B800", "#009acd")) + scale_x_discrete(labels= levels(div_df_melt_simpson$label_pro)) # add sample size info in the x-axis
p_simpson <- p_simpson +
    facet_wrap(~Timepoint) # Creates separate panels for each level of timepoint
p_simpson
# Add p-values to the plot
p1_simpson <- p_simpson +
    stat_compare_means(aes(group = Probiotic),label = "p.format", method = "kruskal.test", 
                       label.y = 0.71*max(div_df_melt_simpson$value)) # Add global p-value
p1_simpson
p1_simpson
# Customize the theme
p2_simpson <- p1_simpson + theme_bw() +
    theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle =  90, vjust = 0.75),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10)
    ) +
    labs(
        x = "Probiotic",
        y = "Simpson diversity"
    )
p2_simpson
# annotate the plot with "a"/"b" for significant differences if required
p2_simpson <- p2_simpson + annotate("text", x=1, y=0.65, label= "b") + 
    annotate("text", x = 2, y=0.65, label = "b") +
    annotate("text", x = 3, y=0.65, label = "a")  
p2_simpson

##### Chao1
#Re-order x-axis for term
div_df_melt_chao1$Term <- factor(div_df_melt_chao1$Term, levels=c("P ", "F "))
#Re-order x-axis for timepoint
div_df_melt_chao1$Timepoint <- factor(div_df_melt_chao1$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
#Re-order x-axis for probiotic
div_df_melt_chao1$Probiotic <- factor(div_df_melt_chao1$Probiotic, levels=c("N", "Y", "F"))
# add sample info as number of samples in each group
div_df_melt_chao1 <- div_df_melt_chao1 %>% 
    group_by(Probiotic) %>% 
    mutate(n_pro= as.character(n())) %>% 
    mutate(label_pro= factor(paste0(Probiotic, '\n n = ', n_pro))) %>% 
    ungroup() %>% 
    select(-n_pro)
# Now use this data frame to plot and add sample size info in the x-axis
p_chao <- ggboxplot(div_df_melt_chao1, x = "Probiotic", y = "value",
                    color = "Probiotic",
                    fill = "Probiotic",
                    legend = "right",
                    scales = "free",
                    add = "dotplot",
                    alpha = 0.2, # Transparency of points
                    ggtheme = theme_bw(),
                    palette = c("#cd5b45", "#008a8a", "#E7B800", "#009acd")) + scale_x_discrete(labels= levels(div_df_melt_chao1$label_pro)) # add sample size info in the x-axis

p_chao <- p_chao +
    facet_wrap(~Timepoint) # Creates separate panels for each level of timepoint
p_chao
# Add p-values to the plot
p1_chao <- p_chao +
    stat_compare_means(aes(group = Probiotic),label = "p.format", method = "kruskal.test", 
                       label.y = 0.98*max(div_df_melt_chao1$value)) # Add global p-value
p1_chao

# Customize the theme
p2_chao <- p1_chao + theme_bw() +
    theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.75),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10)
    ) +
    labs(
        x = "Probiotic",
        y = "Chao1 diversity"
    )
p2_chao
# annotate the plot with "a"/"b" for significant differences if required
p2_chao <- p2_chao + annotate("text", x=1, y=550, label= "b") + 
    annotate("text", x = 2, y=550, label = "b") +
    annotate("text", x = 3, y=550, label = "a")

# Remove legend and margin for observed, shannon and simpson plots
p2_ob <- p2_ob + theme(legend.position = "none",
                       unit(c(0.01,0,0.01,0.01), "cm"),
                       axis.text = element_blank(), axis.ticks.length = unit(0, "mm"))
p2_shannon <- p2_shannon + theme(legend.position = "none",
                                 unit(c(0.01,0,0.01,0.01), "cm"),
                                 axis.text = element_blank(), axis.ticks.length = unit(0, "mm"))
p2_simpson <- p2_simpson + theme(legend.position = "none",
                                 unit(c(0.01,0,0.01,0.01), "cm"),
                                 axis.text = element_blank(), axis.ticks.length = unit(0, "mm"))
p2_chao <- p2_chao + theme(legend.position = "none",
                           unit(c(0.01,0,0.01,0.01), "cm"),
                           axis.text = element_blank(), axis.ticks.length = unit(0, "mm"))                      
# Combine all the plots with cowplot
final_plot <- plot_grid(p2_ob, p2_shannon, p2_simpson, p2_chao, nrow = 1, align = 'v')
final_plot

# save figure in pdf format
ggsave("fig_alpha-div.pdf", plot = final_plot, device = "pdf", path = "path/to/save/figure",
       scale = 1, width = 12.0, height = 7.0, units = "in",
       dpi = 500)
# Save figure in tiff format
ggsave("fig_alpha-div.tiff", plot = final_plot, device = "tiff", path = "path/to/save/figure",
       scale = 1, width = 14.0, height = 7.0, units = "in",
       dpi = 500)
