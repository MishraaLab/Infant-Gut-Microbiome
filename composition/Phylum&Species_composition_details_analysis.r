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
# Set colour palette
col_blind18 <- c("azure2", "cadetblue", "darkgoldenrod1", "dodgerblue4", "mediumpurple4", "lightsteelblue", "indianred3", "mistyrose3", "slategray4", "olivedrab4", "navyblue", "navajowhite2", "lightsteelblue3", "mediumpurple4", "lightgoldenrod1", "deepskyblue4", "cornsilk3", "black")

col_blind5 <- c("cadetblue", "darkgoldenrod1",
                "dodgerblue4", "indianred3", "lightsteelblue",
                "palevioletred3", "slategray4", "olivedrab4", "navyblue", "navajowhite2",  "mediumpurple4", "lightgoldenrod1", "deepskyblue4", "cornsilk3", "black")

cbPalette <- c("burlywood3", "chocolate", "cadetblue3", "cadetblue4")
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
ntaxa(ps)
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

# threshold filtering or rarefaction
ps_filt1 <- prune_samples(sample_sums(ps_ab)>=3000, ps_ab) 
ps_filt1

# threshold filtering - removing very low abundance ASVs
ps_filt2 = filter_taxa(ps_filt1, function(x) sum(x) > 10 , TRUE) 
ps_filt2

ps_p <- ps_filt2 %>%
  subset_taxa(
    Phylum %in% c("p__Actinomycetota", "p__Bacillota", "p__Pseudomonadota", "p__Campylobacterota", "p__Fusobacteriota", "p__Verrucomicrobiota", "p__Desulfobacterota", "p__Unclassified"))
ps_p
# Phylum level relative abundance
pmelt_biol<-tax_glom(ps_p,taxrank = "Phylum") %>%
  transform_sample_counts(function(x) x*100/sum(x)) %>% 
  psmelt ()  %>% 
  group_by(Sample, Term, Probiotic, Timepoint, Phylum)%>%
  #summarize_at("Abundance", mean)
  dplyr::summarise(across("Abundance", ~mean(.x)))
#define abundance cutoff
cutoff<-1
pmelt_biol_filt<-pmelt_biol %>% 
  filter(mean(Abundance) >cutoff)
#cleaning the names
pmelt_biol_filt$Phylum <- str_replace_all(pmelt_biol_filt$Phylum, "p__", " ")
# Agrregate data
dat.agr_filt = aggregate(Abundance~Term+Timepoint+Probiotic+Phylum, data=pmelt_biol_filt, FUN=mean)
dat.agr_filt$Timepoint <- factor(dat.agr_filt$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
dat.agr_filt$Term <- factor(dat.agr_filt$Term, levels=c("P ", "F "))
dat.agr_filt$Probiotic <- factor(dat.agr_filt$Probiotic, levels=c("N", "Y", "F"))
# Plot
fig_p <- ggplot(dat.agr_filt,aes(x=Probiotic,y=Abundance,fill = fct_reorder(Phylum, Abundance))) +
  facet_nested(~Timepoint, scales = "free", switch = "y") +
  theme_bw() +
  geom_bar(stat = "identity", position="fill", colour="white") +
  scale_fill_manual(name = "Phylum", values = col_blind18)+
  theme_bw()+
  ylab("Relative Abundance (%)")+
  xlab("") +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text.x = element_text(face = "bold.italic"))+
  theme(plot.title = element_text(size=10, hjust = 0.5))+
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=10, face="bold"),
        legend.title = element_text(face = "bold"), 
        legend.text.align = 0)
fig_p

# Species level relative abundance
ps_ab_sp_trans <- tax_glom(ps_ab,taxrank = "Species") %>%
  transform_sample_counts(function(x) x*100/sum(x))
# Species level relative abundance- top 50
top50_species <- names(head(sort(colSums(otu_table(ps_ab_sp_trans)), decreasing = TRUE), 50))
top50_species
# Combine phyloseq with top50 species
comb_sp_top50 <- prune_taxa(top50_species, ps_ab_sp_trans)
comb_sp_top50
# Melting phyloseq object
comb_sp_top50_pmelt <- psmelt(comb_sp_top50) 
head(comb_sp_top50_pmelt, 5)
# Cleaning the names
comb_sp_top50_pmelt$Species<- str_replace_all(comb_sp_top50_pmelt$Species, "s__", " ")
# Aggregate data
dat.agr_filt = aggregate(Abundance~Term+Timepoint+Probiotic+Species, data=comb_sp_top50_pmelt, FUN=mean)
dat.agr_filt$Timepoint <- factor(dat.agr_filt$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
dat.agr_filt$Term <- factor(dat.agr_filt$Term, levels=c("P ", "F "))
dat.agr_filt$Probiotic <- factor(dat.agr_filt$Probiotic, levels=c("N", "Y", "F"))
#Set gen color
sp_col_t50 <- c ("#CB86B0", "#90509D", "#551A8B", "#FFCD62", "#FFB01A", "#EE8807", "#E4E6AF", "#B0C073", "#7C9938", "#A5D4E2", "#2D9FA2", "#098B81", "#7896CA", "#3C4BAA", "#00008B", "#FFD92F", "#FF7F00", "#C74B1D", "#77CCCC", "#44AAAA", "#117777", "#DDDD77", "#AAAA44", "#777711", "#DDAA77", "#AA7744", "#774411", "#f4e6d2", "#E5C494", "#d6a256", "#c5dcee", "#7aadd7", "#3274a9", "#c4cec0", "#8a9e82", "#3c5e30", "#ffa07a", "#cd8162", "#8b0a50", "#F2F2F2", "#cdc9c9", "#8b8b83", "#cd3333", "#88CCAA", "#75ebce", "#11775e", "#f39c9c", "#bb6a6b", "#8b3a3a", "grey20")
# Plot
fig_s_t50 <- ggplot(dat.agr_filt,aes(x=Probiotic,y=Abundance,fill = fct_reorder(Species, Abundance))) +
  facet_nested(~Timepoint, scales = "free", switch = "y") +
  theme_bw() +
  geom_bar(stat = "identity", position="fill", colour="white") +
  scale_fill_manual(name = "Species", values = sp_col_t50)+
  theme_bw()+
  ylab("Relative Abundance (%)")+
  xlab("") +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text.x = element_text(face = "bold.italic"))+
  theme(plot.title = element_text(size=10, hjust = 0.5))+
  theme(axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size=10, face="bold"),
        legend.margin = margin(0, 0, 0, 0), # No margin around legend box
        legend.spacing.x = unit(0.1, "cm"), # Small horizontal spacing
        legend.spacing.y = unit(0.1, "cm"), # Small vertical spacing
        legend.key.size = unit(0.5, "cm"), # Smaller legend keys
        legend.title = element_text(face = "bold"), 
        legend.text.align = 0)
fig_s_t50

#Combine phylum and species plots
fig_p_s_t50 <- plot_grid(fig_p + theme(legend.position="none"),
                         fig_s_t50 + theme(legend.position="none"),
                         nrow=1, align = 'v', axis = 'lr', rel_widths = c(1,1.05))
fig_p_s_t50
#saving the plot (medium definition)
dev.copy(tiff,'Fig_Phylum_Species_top50.tiff', width=8.5, height=5.5, units="in", res=300)
dev.off()
# Extract legends
legend_p <- get_legend(fig_p + theme(legend.position="bottom",
                                     legend.text = element_text(size=10)))
legend_s <- get_legend(fig_s_t50 + theme(legend.position="bottom",
                                         legend.text = element_text(size=9.5, face="italic")))
# Combine legends
legends <- plot_grid(legend_p, legend_s, nrow=2)
legends
#saving the legend label (medium definition)
dev.copy(tiff,'Fig_Phylum_Species_top50_levels.tiff', width=12.5, height=7, units="in", res=300)
dev.off()
# Final plot with legends
final_fig <- plot_grid(fig_p_s_t50, legends, nrow=2, rel_heights = c(1,2.7))
final_fig
# Save the final figure
ggsave("Fig2_Phylum_Species_top50.tiff", final_fig, width = 10, height = 5.5, dpi = 300)

###################################################### Phylum level relative abundance - preterm & full-term separate ###########################################################

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
ntaxa(ps)
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

# threshold filtering or rarefaction
ps_filt1 <- prune_samples(sample_sums(ps_ab)>=3000, ps_ab) 
ps_filt1

# threshold filtering - removing very low abundance ASVs
ps_filt2 = filter_taxa(ps_filt1, function(x) sum(x) > 10 , TRUE) 
ps_filt2

ps_p <- ps_filt2 %>%
  subset_taxa(
    Phylum %in% c("p__Actinomycetota", "p__Bacillota", "p__Pseudomonadota", "p__Campylobacterota", "p__Fusobacteriota", "p__Verrucomicrobiota", "p__Desulfobacterota", "p__Unclassified"))
ps_p

####### Preterm phylum abundance
# Subset to include only samples belonging to 'pre-term' without rarefaction
subset_pre_ps_rare_phy <- subset_samples(ps_p, Term %in% c("P "))
subset_pre_ps_rare_phy
# Subset to include only samples belonging to 'pre-term' with rarefaction
subset_pre_ps_norare_phy <- subset_samples(ps_ab, Term %in% c("P "))
subset_pre_ps_norare_phy
# Phylum level relative abundance -rare
pmelt_biol_rare<-tax_glom(subset_pre_ps_rare_phy,taxrank = "Phylum") %>%
  transform_sample_counts(function(x) x*100/sum(x)) %>% 
  psmelt ()  %>% 
  group_by(Sample, Term, Probiotic, Timepoint, Phylum)%>%
  #summarize_at("Abundance", mean)
  dplyr::summarise(across("Abundance", ~mean(.x)))

# Phylum level relative abundance - no rare
pmelt_biol_no_rare<-tax_glom(subset_pre_ps_norare_phy,taxrank = "Phylum") %>%
  transform_sample_counts(function(x) x*100/sum(x)) %>% 
  psmelt ()  %>% 
  group_by(Sample, Term, Probiotic, Timepoint, Phylum)%>%
  #summarize_at("Abundance", mean)
  dplyr::summarise(across("Abundance", ~mean(.x)))
#define abundance cutoff - rare
cutoff<-1
pmelt_biol_filt_rare<-pmelt_biol_rare %>% 
  filter(mean(Abundance) >cutoff)
#define abundance cutoff - no rare
cutoff<-1
pmelt_biol_filt_no_rare<-pmelt_biol_no_rare %>% 
  filter(mean(Abundance) >cutoff)
#cleaning the names - rare
pmelt_biol_filt_rare$Phylum <- str_replace_all(pmelt_biol_filt_rare$Phylum, "p__", "")
#cleaning the names - no rare
pmelt_biol_filt_no_rare$Phylum <- str_replace_all(pmelt_biol_filt_no_rare$Phylum, "p__", "")
# Set colour palette
col_blind18 <- c("azure2", "cadetblue", "mediumpurple4",  "lightsteelblue", "black", "indianred3", "dodgerblue4", "mistyrose3", "slategray4", "olivedrab4", "navyblue", "navajowhite2", "lightsteelblue3", "lightgoldenrod1", "deepskyblue4", "cornsilk3")
# Select phylum level relative abundance separately
#dat.agr_melt_phy_preterm <- pmelt_biol_filt %>%
#  filter(Term == "P ")
#dat.agr_melt_phy_fullterm <- pmelt_biol_filt %>%
#  filter(Term == "F ")
# Agrregate phylum data - rare
dat.agr_melt_phy_pre_rare = aggregate(Abundance~Term+Timepoint+Probiotic+Phylum, data=pmelt_biol_filt_rare, FUN=mean)
# Agrregate phylum data - no rare
dat.agr_melt_phy_pre_no_rare = aggregate(Abundance~Term+Timepoint+Probiotic+Phylum, data=pmelt_biol_filt_no_rare, FUN=mean)
# Re-order factors - preterm for rare
dat.agr_melt_phy_pre_rare$Timepoint <- factor(dat.agr_melt_phy_pre_rare$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
dat.agr_melt_phy_pre_rare$Term <- factor(dat.agr_melt_phy_pre_rare$Term, levels=c("P ", "F "))
dat.agr_melt_phy_pre_rare$Probiotic <- factor(dat.agr_melt_phy_pre_rare$Probiotic, levels=c("N", "Y", "F"))
# Re-order factors - preterm for no rare
dat.agr_melt_phy_pre_no_rare$Timepoint <- factor(dat.agr_melt_phy_pre_no_rare$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
dat.agr_melt_phy_pre_no_rare$Term <- factor(dat.agr_melt_phy_pre_no_rare$Term, levels=c("P ", "F "))
dat.agr_melt_phy_pre_no_rare$Probiotic <- factor(dat.agr_melt_phy_pre_no_rare$Probiotic, levels=c("N", "Y", "F"))

# Plot phylum - rare
fig_pre_phy_rare <- ggplot(dat.agr_melt_phy_pre_rare,aes(x=Probiotic,y=Abundance,fill = fct_reorder(Phylum, Abundance))) +
  facet_nested(~Timepoint, scales = "free", switch = "y") +
  theme_bw() +
  geom_bar(stat = "identity", position="fill", colour="white") +
  scale_fill_manual(name = "Phylum", values = col_blind18)+
  ylab("Relative Abundance (%)")+
  ggtitle("Preterm") +
  xlab("") +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text.x = element_text(face = "bold.italic"))+
  theme(plot.title = element_text(size=10, hjust = 0.5))+
  theme(axis.title = element_text(face = "bold"),
        legend.margin = margin(0, 0, 0, 0), # No margin around legend box
        legend.spacing.x = unit(0.1, "cm"), # Small horizontal spacing
        legend.spacing.y = unit(0.1, "cm"), # Small vertical spacing
        legend.key.size = unit(0.5, "cm"), # Smaller legend keys
        legend.title = element_text(face = "bold"), 
        legend.text.align = 0) +
  theme(plot.margin = unit(c(0.1,0,0.1,0.1), "cm"))
fig_pre_phy_rare
# Set colour palette for no rare
col_blind18 <- c("olivedrab4", "cadetblue", "mediumpurple4",  "lightsteelblue", "black", "dodgerblue4", "indianred3", "mistyrose3", "slategray4", "azure2", "navyblue", "navajowhite2", "lightsteelblue3", "lightgoldenrod1", "deepskyblue4", "cornsilk3")
# Plot phylum - no rare
fig_pre_phy_no_rare <- ggplot(dat.agr_melt_phy_pre_no_rare,aes(x=Probiotic,y=Abundance,fill = fct_reorder(Phylum, Abundance))) +
    facet_nested(~Timepoint, scales = "free", switch = "y") +
    theme_bw() +
    geom_bar(stat = "identity", position="fill", colour="white") +
    scale_fill_manual(name = "Phylum", values = col_blind18)+
    ylab("Relative Abundance (%)")+
    ggtitle("Preterm") +
    xlab("") +
    theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
              unit(0.1, "lines")) +
    theme(strip.text.x = element_text(face = "bold.italic"))+
    theme(plot.title = element_text(size=10, hjust = 0.5))+
    theme(axis.title = element_text(face = "bold"),
          legend.margin = margin(0, 0, 0, 0), # No margin around legend box
          legend.spacing.x = unit(0.1, "cm"), # Small horizontal spacing
          legend.spacing.y = unit(0.1, "cm"), # Small vertical spacing
          legend.key.size = unit(0.5, "cm"), # Smaller legend keys
          legend.title = element_text(face = "bold"), 
          legend.text.align = 0) +
    theme(plot.margin = unit(c(0.1,0,0.1,0.1), "cm"))
fig_pre_phy_no_rare
# (Without rarefaction genrated actual composition by retaining Bacteriodota)

######### Phylum abundance - full-term
# Subset to include only samples belonging to 'pre-term' with rarefaction
subset_full_ps_rare_phy <- subset_samples(ps_p, Term %in% c("F "))
subset_full_ps_rare_phy
# Subset to include only samples belonging to 'pre-term' without rarefaction
subset_full_ps_norare_phy <- subset_samples(ps_ab, Term %in% c("F "))
subset_full_ps_norare_phy
# Phylum level relative abundance -rare
pmelt_biol_rare<-tax_glom(subset_full_ps_rare_phy,taxrank = "Phylum") %>%
  transform_sample_counts(function(x) x*100/sum(x)) %>% 
  psmelt ()  %>% 
  group_by(Sample, Term, Probiotic, Timepoint, Phylum)%>%
  #summarize_at("Abundance", mean)
  dplyr::summarise(across("Abundance", ~mean(.x)))

# Phylum level relative abundance - no rare
pmelt_biol_no_rare<-tax_glom(subset_full_ps_norare_phy,taxrank = "Phylum") %>%
  transform_sample_counts(function(x) x*100/sum(x)) %>% 
  psmelt ()  %>% 
  group_by(Sample, Term, Probiotic, Timepoint, Phylum)%>%
  #summarize_at("Abundance", mean)
  dplyr::summarise(across("Abundance", ~mean(.x)))
#define abundance cutoff - rare
cutoff<-1
pmelt_biol_filt_rare<-pmelt_biol_rare %>% 
  filter(mean(Abundance) >cutoff)
#define abundance cutoff - no rare
cutoff<-1
pmelt_biol_filt_no_rare<-pmelt_biol_no_rare %>% 
  filter(mean(Abundance) >cutoff)
#cleaning the names - rare
pmelt_biol_filt_rare$Phylum <- str_replace_all(pmelt_biol_filt_rare$Phylum, "p__", "")
#cleaning the names - no rare
pmelt_biol_filt_no_rare$Phylum <- str_replace_all(pmelt_biol_filt_no_rare$Phylum, "p__", "")
# Select phylum level relative abundance separately
#dat.agr_melt_phy_preterm <- pmelt_biol_filt %>%
#  filter(Term == "P ")
#dat.agr_melt_phy_fullterm <- pmelt_biol_filt %>%
#  filter(Term == "F ")
# Agrregate phylum data - rare
dat.agr_melt_phy_full_rare = aggregate(Abundance~Term+Timepoint+Probiotic+Phylum, data=pmelt_biol_filt_rare, FUN=mean)
# Agrregate phylum data - no rare
dat.agr_melt_phy_full_no_rare = aggregate(Abundance~Term+Timepoint+Probiotic+Phylum, data=pmelt_biol_filt_no_rare, FUN=mean)
# Re-order factors - preterm for rare
dat.agr_melt_phy_full_rare$Timepoint <- factor(dat.agr_melt_phy_full_rare$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
dat.agr_melt_phy_full_rare$Term <- factor(dat.agr_melt_phy_full_rare$Term, levels=c("P ", "F "))
dat.agr_melt_phy_full_rare$Probiotic <- factor(dat.agr_melt_phy_full_rare$Probiotic, levels=c("N", "Y", "F"))
# Re-order factors - preterm for no rare
dat.agr_melt_phy_full_no_rare$Timepoint <- factor(dat.agr_melt_phy_full_no_rare$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
dat.agr_melt_phy_full_no_rare$Term <- factor(dat.agr_melt_phy_full_no_rare$Term, levels=c("P ", "F "))
dat.agr_melt_phy_full_no_rare$Probiotic <- factor(dat.agr_melt_phy_full_no_rare$Probiotic, levels=c("N", "Y", "F"))
# Set colour palette
col_blind18 <- c("azure2", "cadetblue", "mediumpurple4",  "lightsteelblue", "black", "indianred3", "dodgerblue4", "mistyrose3", "slategray4", "olivedrab4", "navyblue", "navajowhite2", "lightsteelblue3", "lightgoldenrod1", "deepskyblue4", "cornsilk3")
# Plot phylum - rare
fig_full_phy_rare <- ggplot(dat.agr_melt_phy_full_rare,aes(x=Probiotic,y=Abundance,fill = fct_reorder(Phylum, Abundance))) +
  facet_nested(~Timepoint, scales = "free", switch = "y") +
  theme_bw() +
  geom_bar(stat = "identity", position="fill", colour="white") +
  scale_fill_manual(name = "Phylum", values = col_blind18)+
  ylab("Relative Abundance (%)")+
  ggtitle("Full-term") +
  xlab("") +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text.x = element_text(face = "bold.italic"))+
  theme(plot.title = element_text(size=10, hjust = 0.5))+
  theme(axis.title = element_text(face = "bold"),
        legend.margin = margin(0, 0, 0, 0), # No margin around legend box
        legend.spacing.x = unit(0.1, "cm"), # Small horizontal spacing
        legend.spacing.y = unit(0.1, "cm"), # Small vertical spacing
        legend.key.size = unit(0.5, "cm"), # Smaller legend keys
        legend.title = element_text(face = "bold"), 
        legend.text.align = 0) +
  theme(plot.margin = unit(c(0.1,0,0.1,0.1), "cm"))
fig_full_phy_rare
# Set colour palette for no rare
col_blind18 <- c("olivedrab4", "cadetblue", "mediumpurple4", "black", "lightsteelblue", "dodgerblue4", "indianred3", "slategray4", "mistyrose3", "azure2", "navyblue", "navajowhite2", "lightsteelblue3", "lightgoldenrod1", "deepskyblue4", "cornsilk3")
# Plot phylum - no rare
fig_full_phy_no_rare <- ggplot(dat.agr_melt_phy_full_no_rare,aes(x=Probiotic,y=Abundance,fill = fct_reorder(Phylum, Abundance))) +
    facet_nested(~Timepoint, scales = "free", switch = "y") +
    theme_bw() +
    geom_bar(stat = "identity", position="fill", colour="white") +
    scale_fill_manual(name = "Phylum", values = col_blind18)+
    ylab("Relative Abundance (%)")+
    ggtitle("Full-term") +
    xlab("") +
    theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
              unit(0.1, "lines")) +
    theme(strip.text.x = element_text(face = "bold.italic"))+
    theme(plot.title = element_text(size=10, hjust = 0.5))+
    theme(axis.title = element_text(face = "bold"),
          legend.margin = margin(0, 0, 0, 0), # No margin around legend box
          legend.spacing.x = unit(0.1, "cm"), # Small horizontal spacing
          legend.spacing.y = unit(0.1, "cm"), # Small vertical spacing
          legend.key.size = unit(0.5, "cm"), # Smaller legend keys
          legend.title = element_text(face = "bold"), 
          legend.text.align = 0) +
    theme(plot.margin = unit(c(0.1,0,0.1,0.1), "cm"))
fig_full_phy_no_rare
# (Without rarefaction genrated actual composition by retaining Bacteriodota)
fig_pre_phy_no_rare | fig_full_phy_no_rare

######################################################## Species-level abundance - pre vs full-term ##########################################################
####### Preterm species abundance
# Subset to include only samples belonging to 'pre-term' without rarefaction
subset_pre_ps_rare_sp <- subset_samples(ps_p, Term %in% c("P "))
subset_pre_ps_rare_sp
# Subset to include only samples belonging to 'pre-term' with rarefaction
subset_pre_ps_norare_sp <- subset_samples(ps_ab, Term %in% c("P "))
subset_pre_ps_norare_sp
# Species level relative abundance -rare
trans_biol_sp_pre_rare<-tax_glom(subset_pre_ps_rare_sp,taxrank = "Species") %>%
  transform_sample_counts(function(x) x*100/sum(x)) 
# Species level relative abundance - no rare
trans_biol_sp_pre_no_rare<-tax_glom(subset_pre_ps_norare_sp,taxrank = "Species") %>%
  transform_sample_counts(function(x) x*100/sum(x)) 
# Species level relative abundance- top 50 rare
top50_species_pre_rare <- names(head(sort(colSums(otu_table(trans_biol_sp_pre_rare)), decreasing = TRUE), 50))
top50_species_pre_rare
# Species level relative abundance- top 50 no rare
top50_species_pre_no_rare <- names(head(sort(colSums(otu_table(trans_biol_sp_pre_no_rare)), decreasing = TRUE), 50))
top50_species_pre_no_rare
# Combine phyloseq with top50 species rare
comb_sp_top50_pre_rare <- prune_taxa(top50_species_pre_rare, trans_biol_sp_pre_rare)
comb_sp_top50_pre_rare
# Combine phyloseq with top50 species no rare
comb_sp_top50_pre_no_rare <- prune_taxa(top50_species_pre_no_rare, trans_biol_sp_pre_no_rare)
comb_sp_top50_pre_no_rare
# Melting phyloseq object - rare
comb_sp_top50_pmelt_pre_rare <- psmelt(comb_sp_top50_pre_rare)
head(comb_sp_top50_pmelt_pre_rare, 5)
# Melting phyloseq object - no rare
comb_sp_top50_pmelt_pre_no_rare <- psmelt(comb_sp_top50_pre_no_rare)
head(comb_sp_top50_pmelt_pre_no_rare, 5)
# Cleaning the names - rare 
comb_sp_top50_pmelt_pre_rare$Species<- str_replace_all(comb_sp_top50_pmelt_pre_rare$Species, "s__", "")
# Cleaning the names - no rare 
comb_sp_top50_pmelt_pre_no_rare$Species<- str_replace_all(comb_sp_top50_pmelt_pre_no_rare$Species, "s__", "")
# Aggregate data - rare
dat.agr_filt_sp_pre_rare = aggregate(Abundance~Term+Timepoint+Probiotic+Species, data=comb_sp_top50_pmelt_pre_rare, FUN=mean)
# Aggregate data - no rare
dat.agr_filt_sp_pre_no_rare = aggregate(Abundance~Term+Timepoint+Probiotic+Species, data=comb_sp_top50_pmelt_pre_no_rare, FUN=mean)
# Re-order names - rare
dat.agr_filt_sp_pre_rare$Timepoint <- factor(dat.agr_filt_sp_pre_rare$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
#dat.agr_filt_sp_pre_rare$Term <- factor(dat.agr_filt_sp_pre_rare$Term, levels=c("P ", "F "))
dat.agr_filt_sp_pre_rare$Probiotic <- factor(dat.agr_filt_sp_pre_rare$Probiotic, levels=c("N", "Y", "F"))
# Re-order names - no rare
dat.agr_filt_sp_pre_no_rare$Timepoint <- factor(dat.agr_filt_sp_pre_no_rare$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
#dat.agr_filt_sp_pre_no_rare$Term <- factor(dat.agr_filt_sp_pre_no_rare$Term, levels=c("P ", "F "))
dat.agr_filt_sp_pre_no_rare$Probiotic <- factor(dat.agr_filt_sp_pre_no_rare$Probiotic, levels=c("N", "Y", "F"))
# Set sp color - pre_rare
sp_col_t50 <- c ("#CB86B0", "#90509D", "#551A8B", "#FFCD62", "#FFB01A", "#EE8807", "#E4E6AF", "#B0C073", "#7C9938", "#A5D4E2", "#2D9FA2", "#098B81", "#7896CA", "#3C4BAA", "#00008B", "#FFD92F", "#FF7F00", "#C74B1D", "#77CCCC", "#44AAAA", "#117777", "#DDDD77", "#AAAA44", "#777711", "#DDAA77", "#AA7744", "#774411", "#f4e6d2", "#E5C494", "#d6a256", "#c5dcee", "#7aadd7", "#3274a9", "#c4cec0", "#8a9e82", "#3c5e30", "#ffa07a", "#cd8162", "#F2F2F2", "#8b0a50", "#cdc9c9",  "#cd3333", "#8b8b83", "#f39c9c", "#88CCAA", "#75ebce", "#bb6a6b", "#11775e", "#8b3a3a", "grey20")  
# Plot species - rare
fig_sp_pre_rare <- ggplot(dat.agr_filt_sp_pre_rare,aes(x=Probiotic,y=Abundance,fill = fct_reorder(Species, Abundance))) +
  facet_nested(~Timepoint, scales = "free", switch = "y") +
  theme_bw() +
  geom_bar(stat = "identity", position="fill", colour="white") +
  scale_fill_manual(name = "Species", values = sp_col_t50)+
  ylab("Relative Abundance (%)")+
  ggtitle("Preterm") +
  xlab("") +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text.x = element_text(face = "bold.italic"))+
  theme(plot.title = element_text(size=10, hjust = 0.5))+
  theme(axis.title = element_text(face = "bold"),
        legend.text = element_text(size=9.5, face="italic"), # Small legend italic font
        legend.margin = margin(0, 0, 0, 0), # No margin around legend box
        legend.spacing.x = unit(0.1, "cm"), # Small horizontal spacing
        legend.spacing.y = unit(0.1, "cm"), # Small vertical spacing
        legend.key.size = unit(0.5, "cm"), # Smaller legend keys
        legend.title = element_text(face = "bold"), 
        legend.text.align = 0) +
  theme(plot.margin = unit(c(0.1,0,0.1,0.1), "cm"))
fig_sp_pre_rare
# Set sp color - no rare
sp_col_t50 <- c ("#CB86B0", "#90509D", "#551A8B", "#FFCD62", "#FFB01A", "#EE8807", "#E4E6AF", "#B0C073", "#7C9938", "#A5D4E2", "#2D9FA2", "#098B81", "#7896CA", "#3C4BAA", "#00008B", "#FFD92F", "#FF7F00", "#C74B1D", "#77CCCC", "#44AAAA", "#117777", "#DDDD77", "#AAAA44", "#777711", "#DDAA77", "#AA7744", "#774411", "#f4e6d2", "#E5C494", "#d6a256", "#c5dcee", "#7aadd7", "#3274a9", "#c4cec0", "#8a9e82", "#3c5e30", "#ffa07a", "#cd8162", "#F2F2F2", "#8b0a50", "#cdc9c9",  "#cd3333", "#8b8b83", "#f39c9c", "#88CCAA", "#75ebce", "#bb6a6b", "#11775e", "#8b3a3a", "grey20")  
# Plot species - pre no_rare
fig_sp_pre_no_rare <- ggplot(dat.agr_filt_sp_pre_no_rare,aes(x=Probiotic,y=Abundance,fill = fct_reorder(Species, Abundance))) +
  facet_nested(~Timepoint, scales = "free", switch = "y") +
  theme_bw() +
  geom_bar(stat = "identity", position="fill", colour="white") +
  scale_fill_manual(name = "Species", values = sp_col_t50)+
  ylab("Relative Abundance (%)")+
  ggtitle("Preterm") +
  xlab("") +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text.x = element_text(face = "bold.italic"))+
  theme(plot.title = element_text(size=10, hjust = 0.5))+
  theme(axis.title = element_text(face = "bold"),
        legend.text = element_text(size=9.5, face="italic"), # Small legend italic font
        legend.margin = margin(0, 0, 0, 0), # No margin around legend box
        legend.spacing.x = unit(0.1, "cm"), # Small horizontal spacing
        legend.spacing.y = unit(0.1, "cm"), # Small vertical spacing
        legend.key.size = unit(0.5, "cm"), # Smaller legend keys
        legend.title = element_text(face = "bold"), 
        legend.text.align = 0) +
  theme(plot.margin = unit(c(0.1,0,0.1,0.1), "cm"))
fig_sp_pre_no_rare

####### Full-term species abundance
# Subset to include only samples belonging to 'full-term' without rarefaction
subset_full_ps_rare_sp <- subset_samples(ps_p, Term %in% c("F "))
subset_full_ps_rare_sp
# Subset to include only samples belonging to 'pre-term' with rarefaction
subset_full_ps_norare_sp <- subset_samples(ps_ab, Term %in% c("F "))
subset_full_ps_norare_sp
# Species level relative abundance -rare
trans_biol_sp_full_rare<-tax_glom(subset_full_ps_rare_sp,taxrank = "Species") %>%
  transform_sample_counts(function(x) x*100/sum(x)) 
# Species level relative abundance - no rare
trans_biol_sp_full_no_rare<-tax_glom(subset_full_ps_norare_sp,taxrank = "Species") %>%
  transform_sample_counts(function(x) x*100/sum(x)) 
# Species level relative abundance- top 50 rare
top50_species_full_rare <- names(head(sort(colSums(otu_table(trans_biol_sp_full_rare)), decreasing = TRUE), 50))
top50_species_full_rare
# Species level relative abundance- top 50 no rare
top50_species_full_no_rare <- names(head(sort(colSums(otu_table(trans_biol_sp_full_no_rare)), decreasing = TRUE), 50))
top50_species_full_no_rare
# Combine phyloseq with top50 species rare
comb_sp_top50_full_rare <- prune_taxa(top50_species_full_rare, trans_biol_sp_full_rare)
comb_sp_top50_full_rare
# Combine phyloseq with top50 species no rare
comb_sp_top50_full_no_rare <- prune_taxa(top50_species_full_no_rare, trans_biol_sp_full_no_rare)
comb_sp_top50_full_no_rare
# Melting phyloseq object - rare
comb_sp_top50_pmelt_full_rare <- psmelt(comb_sp_top50_full_rare)
head(comb_sp_top50_pmelt_full_rare, 5)
# Melting phyloseq object - no rare
comb_sp_top50_pmelt_full_no_rare <- psmelt(comb_sp_top50_full_no_rare)
head(comb_sp_top50_pmelt_full_no_rare, 5)
# Cleaning the names - rare 
comb_sp_top50_pmelt_full_rare$Species<- str_replace_all(comb_sp_top50_pmelt_full_rare$Species, "s__", "")
# Cleaning the names - no rare 
comb_sp_top50_pmelt_full_no_rare$Species<- str_replace_all(comb_sp_top50_pmelt_full_no_rare$Species, "s__", "")
# Aggregate data - rare
dat.agr_filt_sp_full_rare = aggregate(Abundance~Term+Timepoint+Probiotic+Species, data=comb_sp_top50_pmelt_full_rare, FUN=mean)
# Aggregate data - no rare
dat.agr_filt_sp_full_no_rare = aggregate(Abundance~Term+Timepoint+Probiotic+Species, data=comb_sp_top50_pmelt_full_no_rare, FUN=mean)
# Re-order names - rare
dat.agr_filt_sp_full_rare$Timepoint <- factor(dat.agr_filt_sp_full_rare$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
#dat.agr_filt_sp_full_rare$Term <- factor(dat.agr_filt_sp_full_rare$Term, levels=c("P ", "F "))
dat.agr_filt_sp_full_rare$Probiotic <- factor(dat.agr_filt_sp_full_rare$Probiotic, levels=c("N", "Y", "F"))
# Re-order names - no rare
dat.agr_filt_sp_full_no_rare$Timepoint <- factor(dat.agr_filt_sp_full_no_rare$Timepoint, levels=c("1W", "4W", "4M", "1Y"))
#dat.agr_filt_sp_full_no_rare$Term <- factor(dat.agr_filt_sp_full_no_rare$Term, levels=c("P ", "F "))
dat.agr_filt_sp_full_no_rare$Probiotic <- factor(dat.agr_filt_sp_full_no_rare$Probiotic, levels=c("N", "Y", "F"))
# Set sp color - pre_rare
sp_col_t50 <- c ("#CB86B0", "#90509D", "#551A8B", "#FFCD62", "#FFB01A", "#EE8807", "#E4E6AF", "#B0C073", "#7C9938", "#A5D4E2", "#2D9FA2", "#098B81", "#7896CA", "#3C4BAA", "#00008B", "#FFD92F", "#FF7F00", "#C74B1D", "#77CCCC", "#44AAAA", "#117777", "#DDDD77", "#AAAA44", "#777711", "#DDAA77", "#AA7744", "#774411", "#f4e6d2", "#E5C494", "#d6a256", "#c5dcee", "#7aadd7", "#3274a9", "#c4cec0", "#8a9e82", "#3c5e30", "#ffa07a", "#cd8162", "#F2F2F2", "#8b0a50", "#cdc9c9",  "#cd3333", "#8b8b83", "#f39c9c", "#88CCAA", "#75ebce", "#bb6a6b", "#11775e", "#8b3a3a", "grey20")  
# Plot species - rare
fig_sp_full_rare <- ggplot(dat.agr_filt_sp_full_rare,aes(x=Probiotic,y=Abundance,fill = fct_reorder(Species, Abundance))) +
  facet_nested(~Timepoint, scales = "free", switch = "y") +
  theme_bw() +
  geom_bar(stat = "identity", position="fill", colour="white") +
  scale_fill_manual(name = "Species", values = sp_col_t50)+
  ylab("Relative Abundance (%)")+
  ggtitle("Full-term") +
  xlab("") +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text.x = element_text(face = "bold.italic"))+
  theme(plot.title = element_text(size=10, hjust = 0.5))+
  theme(axis.title = element_text(face = "bold"),
        legend.text = element_text(size=9.5, face="italic"), # Small legend italic font
        legend.margin = margin(0, 0, 0, 0), # No margin around legend box
        legend.spacing.x = unit(0.1, "cm"), # Small horizontal spacing
        legend.spacing.y = unit(0.1, "cm"), # Small vertical spacing
        legend.key.size = unit(0.5, "cm"), # Smaller legend keys
        legend.title = element_text(face = "bold"), 
        legend.text.align = 0) +
  theme(plot.margin = unit(c(0.1,0,0.1,0.1), "cm"))
fig_sp_full_rare
# Set sp color - no rare
sp_col_t50 <- c ("#CB86B0", "#90509D", "#551A8B", "#FFCD62", "#FFB01A", "#EE8807", "#E4E6AF", "#B0C073", "#7C9938", "#A5D4E2", "#2D9FA2", "#098B81", "#7896CA", "#3C4BAA", "#00008B", "#FFD92F", "#FF7F00", "#C74B1D", "#77CCCC", "#44AAAA", "#117777", "#DDDD77", "#AAAA44", "#cdc9c9", "#DDAA77", "#AA7744", "#774411", "#f4e6d2", "#E5C494", "#d6a256", "#c5dcee", "#7aadd7", "#3274a9", "#8b0a50", "#c4cec0", "#8a9e82", "#11775e",  "#3c5e30", "#cd8162", "#cd3333", "#88CCAA", "#777711", "#75ebce", "#bb6a6b", "#8b8b83", "#8b3a3a", "grey20", "#F2F2F2", "#f39c9c", "#999999")
# Plot species - pre no_rare
fig_sp_full_no_rare <- ggplot(dat.agr_filt_sp_full_no_rare,aes(x=Probiotic,y=Abundance,fill = fct_reorder(Species, Abundance))) +
    facet_nested(~Timepoint, scales = "free", switch = "y") +
    theme_bw() +
    geom_bar(stat = "identity", position="fill", colour="white") +
    scale_fill_manual(name = "Species", values = sp_col_t50)+
    ylab("Relative Abundance (%)")+
    ggtitle("Full-term") +
    xlab("") +
    theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
              unit(0.1, "lines")) +
    theme(strip.text.x = element_text(face = "bold.italic"))+
    theme(plot.title = element_text(size=10, hjust = 0.5))+
    theme(axis.title = element_text(face = "bold"),
          legend.text = element_text(size=9.5, face="italic"), # Small legend italic font
          legend.margin = margin(0, 0, 0, 0), # No margin around legend box
          legend.spacing.x = unit(0.1, "cm"), # Small horizontal spacing
          legend.spacing.y = unit(0.1, "cm"), # Small vertical spacing
          legend.key.size = unit(0.5, "cm"), # Smaller legend keys
          legend.title = element_text(face = "bold"), 
          legend.text.align = 0) +
    theme(plot.margin = unit(c(0.1,0,0.1,0.1), "cm"))
fig_sp_full_no_rare

# (Without rarefaction genrated actual composition by retaining Bacteriodota)
# Combine plots - phylum 
fig_pre_phy_no_rare | fig_full_phy_no_rare
# Combine plots - species
fig_sp_pre_no_rare | fig_sp_full_no_rare

# Remove y-axis label, text, and ticks from full-term plot - phylum
fig_full_phy_no_rare_clean_y <- fig_full_phy_no_rare + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
fig_full_phy_no_rare_clean_y
#Combine pre&full-term plots - phylum
fig_p_f_abund_phy <- plot_grid(fig_pre_phy_no_rare + theme(legend.position="none"),
fig_full_phy_no_rare_clean_y + theme(legend.position="none"),
nrow=1, align = 'h', axis = 'lr', rel_widths = c(1.5,0.75))
fig_p_f_abund_phy
# Saving the plot (medium definition)
dev.copy(tiff,'Fig_Phylum_pre&full_separate_without_legend.tiff', width=4.0, height=5.5, units="in", res=300)
dev.off()
# Remove y-axis label, text, and ticks from full-term plot - species
fig_sp_full_no_rare_clean_y <- fig_sp_full_no_rare + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
fig_sp_full_no_rare_clean_y
#Combine pre&full plots - species
fig_p_f_abund_sp <- plot_grid(fig_sp_pre_no_rare + theme(legend.position="none"),
                               fig_sp_full_no_rare_clean_y + theme(legend.position="none"),
                               nrow=1, align = 'h', axis = 'lr', rel_widths = c(1.5,0.75))
fig_p_f_abund_sp

# Combine phylum and species plots
fig_p_f_abund_phy | fig_p_f_abund_sp
# Saving the plot legend (medium definition)
dev.copy(tiff,'Fig_Phylum_Species_top50_pre&full_separate_without_legend.tiff', width=11.2, height=6.2, units="in", res=300)
dev.off()
# Combine phylum and species plots with legends
fig_p_f_abund_sp_with_legend <- plot_grid(fig_sp_pre_no_rare + theme(legend.position="right"),
                              fig_sp_full_no_rare + theme(legend.position="right"),
                              nrow=1, align = 'h', axis = 'lr', rel_widths = c(1.0,0.95))
fig_p_f_abund_sp_with_legend
# Saving the plot legend (medium definition with legends)
dev.copy(tiff,'Fig_Phylum_Species_top50_pre&full_separate_with_legend.tiff', width=18.5, height=5.5, units="in", res=300)
dev.off()

# Extract legends from pre & full term plots - phylum
legend_fig_pre_phy_no_rare <- get_legend(fig_pre_phy_no_rare + theme(legend.position="right",
                                                                     legend.text = element_text(size=9.5)))
legend_fig_full_phy_no_rare <- get_legend(fig_full_phy_no_rare + theme(legend.position="right",
                                                                       legend.text = element_text(size=9.5)))
# Combine legends
legends_phy <- plot_grid(legend_fig_pre_phy_no_rare, legend_fig_full_phy_no_rare, nrow=2)
legends_phy
# Saving the legend label (medium definition)
dev.copy(tiff,'Fig_phy_levels_pre&full_separate.tiff', width=6.5, height=7, units="in", res=300)
dev.off()
# Extract legends from pre & full term plots - species
legend_fig_sp_pre_no_rare <- get_legend(fig_sp_pre_no_rare + theme(legend.position="right",
                                                             legend.text = element_text(size=9)))
legend_fig_sp_full_no_rare <- get_legend(fig_sp_full_no_rare + theme(legend.position="right",
                                                                 legend.text = element_text(size=9)))
# Combine legends - species
legends_sp <- plot_grid(legend_fig_sp_pre_no_rare, legend_fig_sp_full_no_rare, nrow=2)
legends_sp                                                                 
# Saving the legend label (medium definition)
dev.copy(tiff,'Fig_sp_levels_pre&full_separate.tiff', width=12.5, height=7, units="in", res=300)
dev.off()

############### Not required (saved for future work) #################
# Remove y-axis label, text, and ticks from preterm plot - phylum
fig_pre_phy_no_rare_clean_y <- fig_pre_phy_no_rare + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
fig_pre_phy_no_rare_clean_y
# Remove y-axis label, text, and ticks from full-term plot - phylum
fig_full_phy_no_rare_clean_y <- fig_full_phy_no_rare + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
fig_full_phy_no_rare_clean_y

# Remove y-axis label, text, and ticks from preterm plot - species
fig_sp_pre_no_rare_clean_y <- fig_sp_pre_no_rare + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
# Remove y-axis label, text, and ticks from full-term plot - species
fig_sp_full_no_rare_clean_y <- fig_sp_full_no_rare + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
# Remove y-axis label, text, and ticks from full-term plot
fig_s_t50_full_clean_y <- fig_s_t50_full + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
fig_s_t50_full_clean_y
