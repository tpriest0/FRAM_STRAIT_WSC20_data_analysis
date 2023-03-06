################### 

##### Script containing the steps used to process and create supplementary
##### figures for the paper "Spatial heterogeneity in carbohydrates and their 
##### utilisation by microbial communities in the high North Atlantic"

###################
# Set this variable to download datafiles directory
datafiles_dir = ('D:/ownCloud/PhD/Arctic Bacterioplankton/Fram Strait/MSM95/r_script/data_files')

# Define working directory
setwd(datafiles_dir)

# Create directory ready for figure output
output_figures_dir = "figures_output"
output_figures_dir_path <- file.path(datafiles_dir, output_figures_dir)
dir.create(output_figures_dir_path)

####################
##### Supplementary Figure S1
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggh4x)

# Import monosaccharide data
monosach_wide=read.table("FRAM_WSC20_monosaccharides_r.csv",header=TRUE, sep = ",",as.is=TRUE)

# For this figure, we will not plot the 100m and 200m data for stations S1 and S6
# So remove these and then convert the dataframe to long format and fix
# the depth in the correct order
monosach_long = monosach_wide %>%
  filter(Depth == "SRF" | 
           Depth == "BML") %>%
  reshape2::melt(., id.vars=c("Sample","Location","Station","Depth"),
                 variable.name="Compound", value.name="Concentration") %>%
  arrange(Location,Station,Compound) %>%
  mutate(Depth = fct_relevel(Depth, c("SRF","BML")))

# Fix factor order
monosach_long$Sample <- factor(monosach_long$Sample,levels=unique(monosach_long$Sample))
monosach_long$Depth <- factor(monosach_long$Depth,levels=unique(monosach_long$Depth))
monosach_long$Station <- factor(monosach_long$Station,levels=unique(monosach_long$Station))
monosach_long$Compound <- factor(monosach_long$Compound,levels=unique(monosach_long$Compound))

# Plot Supplementary_Figure_S1
Supplementary_Figure_S1 <- ggplot(monosach_long, aes(x=Concentration, y=Compound)) + 
  geom_bar(aes(fill = Compound), colour = "black", stat = "identity", position = "dodge") + 
  labs(y="Sample", x="Concentration µg L") + 
  scale_fill_brewer(palette = "Accent") + 
  facet_nested(Location+Station~Depth, scales = "fixed") + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white", color="black"), #
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour="gray90"),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        strip.placement = "outside", 
        panel.spacing = unit(0, "lines"), 
        strip.background = element_rect(colour="black", fill="white"), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size=14),
        legend.title = element_blank(), 
        axis.title.x = element_text(size=16), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.text.y = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + 
  guides(colour = guide_legend(override.aes = list(size=5), ncol=1))
Supplementary_Figure_S1

pdf(file="figures_output/Supplementary_Figure_S1.pdf", width=8, height=8)
Supplementary_Figure_S1
dev.off()

#######################

#####

### Supplementary Figure S2

# Import polysaccharide data
polysacch_wide=read.table("FRAM_WSC20_polysaccharides_r.csv",header=TRUE, 
                          sep = ",",as.is=TRUE,check.names=F)

# Import dataframe connecting antibody name to epitope
antibody_to_structure=read.table("antibodies_to_structures.txt",
                                 header=TRUE, sep = "\t", as.is=TRUE,
                                 encoding="UTF-8")

# Need to reformat the dataframe to long format 
# and then average the four values for each epitope in each sample for each solvent
polys_data_mean_per_solvent_long =  polysacch_wide %>%
  subset(., select=-c(Extraction_name)) %>%
  reshape2::melt(., id.vars=c("Station","Depth","Sample","Location","Solvent"),
                 variable.name="Antibody", value.name="Signal") %>%
  aggregate(Signal~Station+Depth+Sample+Location+Antibody+Solvent, data=., FUN=mean)

# Determine total epitope signal for each sample by summing those from the three
# solvents and then combine dataframe with information on epitope structure
# for each antibody
# ALSO, will only show SRF and BML depths, so must remove the 100m and 200m 
# samples for stations S1 and S6
polys_data_total_per_sample_long = polys_data_mean_per_solvent_long %>%
  aggregate(Signal~Station+Depth+Sample+Location+Antibody, data=., FUN=sum) %>%
  left_join(antibody_to_structure, by="Antibody") %>% 
  filter(Depth == "SRF" | 
           Depth == "BML") %>%
  arrange(Location)

# In the main manuscript section, we included the figure (Figure 3) with the 
# 12 most abundant epitopes. For Supplementary Figure S2, we will plot in a 
# similar fashion but show all epitopes. 

# Define factor orders
polys_data_total_per_sample_long$Station <- 
  factor(polys_data_total_per_sample_long$Station, 
         levels=unique(c("S1","S6","S8","S11","S25","S26","S36","S41","S44")))
polys_data_total_per_sample_long$Depth <- 
  factor(polys_data_total_per_sample_long$Depth, 
         levels=unique(c("SRF","BML")))

# Plot Supplementary_Figure_S2
Supplementary_Figure_S2 = 
  ggplot(polys_data_total_per_sample_long) +
  geom_bar(aes(x=ifelse(Depth=="BML",-Signal,NA), 
               y=Station,
               fill = Epitope), 
           stat = "identity", position = "identity",
           width = 1, colour = "black") +  
  geom_bar(aes(x=ifelse(Depth=="SRF",Signal,NA), 
               y=Station,
               fill = Epitope), 
           stat = "identity", position = "identity",
           width = 1, colour = "black") +
  scale_y_discrete(limits=rev) + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(17)) + 
  facet_wrap2(Epitope~., axes = "all", remove_labels = "all") + 
  labs(fill = "Epitope", x = "Epitope abundance (Antibody signal intensity)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(colour = "black", size = 10), 
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8), 
        legend.position = "right", 
        legend.text = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 12),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_blank()) + 
  guides(fill=guide_legend(aes(size=4), ncol = 1))
Supplementary_Figure_S2

# Export Supplementary_Figure_S2
# The figure was further annotated in Inkscape (location added to x-axis and SRF
# and BML labels added above plots as well as reordering of legend)
pdf(file="figures_output/Supplementary_Figure_S2.pdf", width=14, height=8)
Supplementary_Figure_S2
dev.off()



#######################

#####

### Supplementary Figure S3

# Import polysaccharide data
polysacch_wide=read.table("FRAM_WSC20_polysaccharides_r.csv",header=TRUE, 
                          sep = ",",as.is=TRUE,check.names=F)

# Import dataframe connecting antibody name to epitope
antibody_to_structure=read.table("antibodies_to_structures.txt",
                                 header=TRUE, sep = "\t", as.is=TRUE,
                                 encoding="UTF-8")

# Need to reformat the dataframe to long format 
# and then average the four values for each epitope in each sample for each solvent
polys_data_mean_per_solvent_long =  polysacch_wide %>%
  subset(., select=-c(Extraction_name)) %>%
  reshape2::melt(., id.vars=c("Station","Depth","Sample","Location","Solvent"),
                 variable.name="Antibody", value.name="Signal") %>%
  aggregate(Signal~Station+Depth+Sample+Location+Antibody+Solvent, data=., FUN=mean)

# Determine total epitope signal for each sample by summing those from the three
# solvents
# Then determine the mean epitope abundance for the two locations - 
# above slope and open ocean and then combine dataframe with information on
# epitope structure for each antibody
# ALSO, will only show SRF and BML depths, so must remove the 100m and 200m 
# samples for stations S1 and S6
polys_data_mean_per_epitope_per_location_long = polys_data_mean_per_solvent_long %>%
  aggregate(Signal~Station+Depth+Sample+Location+Antibody, data=., FUN=sum) %>%
  aggregate(Signal~Antibody+Depth+Location, data=., FUN=mean) %>%
  left_join(antibody_to_structure, by="Antibody") %>% 
  filter(Depth == "SRF" | 
           Depth == "BML") %>%
  arrange(Location)

# Define factor orders
polys_data_mean_per_epitope_per_location_long$Depth <- 
  factor(polys_data_mean_per_epitope_per_location_long$Depth, 
         levels=unique(c("SRF","BML")))

# Plot Supplementary_Figure_S3
Supplementary_Figure_S3 = 
  ggplot(polys_data_mean_per_epitope_per_location_long) +
  geom_bar(aes(x=Signal, y=Epitope, fill=Location),
           stat = "identity", position = "dodge",
           width = 0.8, colour = "black") +  
  scale_fill_manual(values=c("Above slope" = "#009292", "Open ocean" = "#b66dff"))  + 
  scale_y_discrete(limits=rev) + 
  facet_grid(.~Depth, scales="fixed") + 
  labs(fill = "Epitope", x = "Relative epitope abundance") + 
  theme_bw() + 
  theme(axis.text.x = element_text(colour = "black", size = 10), 
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8), 
        legend.position = "right", 
        legend.text = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 12),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "white"), 
        strip.text.x = element_text(colour = "black", size = 12))
Supplementary_Figure_S3

## Export Supplementary_Figure_S3
pdf(file="figures_output/Supplementary_Figure_S3.pdf", width=10, height=6)
Supplementary_Figure_S3
dev.off()



#######################

#####

### Supplementary Figure S4

## Import monosaccharide data
monosach_wide=read.table("FRAM_WSC20_monosaccharides_r.csv",header=TRUE, sep = ",",as.is=TRUE)

## Process monosaccharide data ready for plotting
# Convert the dataframe to long format, filter for stations s1 and s6 and fix
# the depth in the correct order
monosach_s1_s6_long = monosach_wide %>%
  filter(Station == "S1" | Station == "S6") %>%
  reshape2::melt(., id.vars=c("Sample","Location","Station","Depth"),
                 variable.name="Compound", value.name="Concentration") %>%
  filter(Concentration != "NA") %>%
  mutate(Depth = case_when(
    Depth == "m100" ~ "100m",
    Depth == "m200" ~ "200m",
    TRUE ~ Depth)) %>%
  mutate(Depth = as.factor(Depth)) %>%
  arrange(Compound = fct_relevel(Compound, c("Arabinose","Fucose",
                                             "Galactosamine","Galactose",
                                             "Glucosamine","Glucose",
                                             "Xylose")))

## Plot Supplementary_Figure_S4a

# Fix factor order
monosach_s1_s6_long$Depth <- factor(monosach_s1_s6_long$Depth,
                                    levels=c("SRF","BML","100m","200m"))
monosach_s1_s6_long$Compound <- factor(monosach_s1_s6_long$Compound,
                                       levels=unique(monosach_s1_s6_long$Compound))

# Plot Supplementary_Figure_S4a
Supplementary_Figure_S4a <- ggplot(monosach_s1_s6_long,
                                   aes(x=Concentration, y=Compound)) + 
  geom_bar(aes(fill = Compound), colour = "black", stat = "identity") + 
  labs(y="Sample", x="Concentration µg L") + 
  scale_fill_brewer(palette = "Accent") + 
  facet_nested(Station~as.factor(Depth), scales = "fixed") + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white", color="black"), #
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour="gray90"),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        strip.placement = "outside", 
        panel.spacing = unit(0, "lines"), 
        strip.background = element_rect(colour="black", fill="white"), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size=14),
        legend.title = element_blank(), 
        axis.title.x = element_text(size=16), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.text.y = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + 
  guides(fill = guide_legend(override.aes = list(size=5), ncol=1))
Supplementary_Figure_S4a


## 

## Import polysaccharide data
polysacch_wide=read.table("FRAM_WSC20_polysaccharides_r.csv",header=TRUE, 
                          sep = ",",as.is=TRUE,check.names=F)

# Import dataframe connecting antibody name to epitope
antibody_to_structure=read.table("antibodies_to_structures.txt",
                                 header=TRUE, sep = "\t", as.is=TRUE,
                                 encoding="UTF-8")

## Process polysaccharide data ready for plotting
# Reformat the dataframe to long format 
# then average the four values for each epitope in each sample for each solvent
# then sum up these mean values for each epitope in each sample
# and then filter for stations S1 and S6
polys_s1_s6_long =  polysacch_wide %>%
  subset(., select=-c(Extraction_name)) %>%
  reshape2::melt(., id.vars=c("Station","Depth","Sample","Location","Solvent"),
                 variable.name="Antibody", value.name="Signal") %>%
  aggregate(Signal~Station+Depth+Sample+Location+Antibody+Solvent, data=., FUN=mean) %>%
  aggregate(Signal~Station+Depth+Sample+Location+Antibody, data=., FUN=sum) %>%
  left_join(antibody_to_structure, by="Antibody") %>%
  filter(Station == "S1" | Station == "S6") %>%
  arrange(Station,Antibody) %>%
  filter(Signal != "NA") %>%
  mutate(Depth = case_when(
    Depth == "m100" ~ "100m",
    Depth == "m200" ~ "200m",
    TRUE ~ Depth)) %>%
  mutate(Depth = as.factor(Depth))

## Plot Supplementary_Figure_S4b

# Fix factor order
polys_s1_s6_long$Depth <- factor(polys_s1_s6_long$Depth,
                                    levels=c("SRF","BML","100m","200m"))

Supplementary_Figure_S4b <- ggplot(polys_s1_s6_long,
                                   aes(x=Signal, y=Epitope)) + 
  geom_bar(aes(fill = Epitope), colour = "black", stat = "identity") + 
  scale_y_discrete(limits=rev) + 
  labs(y="Sample", x="Relative epitope abundance") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(17)) + 
  facet_nested(Station~as.factor(Depth), scales = "fixed") + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white", color="black"), #
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour="gray90"),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        strip.placement = "outside", 
        panel.spacing = unit(0, "lines"), 
        strip.background = element_rect(colour="black", fill="white"), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size=14),
        legend.title = element_blank(), 
        axis.title.x = element_text(size=16), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.text.y = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + 
  guides(fill = guide_legend(override.aes = list(size=5), ncol=1))
Supplementary_Figure_S4b

## Combine and export
pdf(file="figures_output/Supplementary_Figure_S4.pdf", width=14, height=8)
Supplementary_Figure_S4a/Supplementary_Figure_S4b
dev.off()




#######################

#####

### Supplementary Figure S5

library(factoextra)
library(amap)
library(ggdendro)

## Import MASH distance data
metag_mash_dist_long=read.table("arctic_metagenomes_all_mash_distances.csv",
                                header=TRUE, sep = ",",as.is=TRUE)

# Import metagenome metadata
metag_metadata=read.table("arctic_metagenomes_metadata.csv",
                                header=TRUE, sep = ",",as.is=TRUE)

# Reformat to wide and set as distance matrix
metag_mash_dist_matrix <- metag_mash_dist_long %>%
  spread(ref, mash_dist, fill=0) %>%
  column_to_rownames(var="query") %>% 
  as.matrix(.)

## Perform hierarchical clustering analysis
# We will test what is the best clustering method

# Test which clustering method best represents the variation in the data
# First perform hierarchical clustering with different algorithms
metag_mash_clust_comp <- hclust(Dist(
  metag_mash_dist_matrix, method="correlation"),"complete")	
metag_mash_clust_avg <- hclust(Dist(
  metag_mash_dist_matrix, method="correlation"),"average")
metag_mash_clust_sing <- hclust(Dist(
  metag_mash_dist_matrix, method="correlation"),"single")
metag_mash_clust_ward <- hclust(Dist(
  metag_mash_dist_matrix, method="correlation"),"ward.D")	

# Calculate cophenetic correlation for each method
metag_mash_clust_comp_coph <- cophenetic(metag_mash_clust_comp)
metag_mash_clust_avg_coph <- cophenetic(metag_mash_clust_avg)
metag_mash_clust_single_coph <- cophenetic(metag_mash_clust_sing)
metag_mash_clust_ward_coph <- cophenetic(metag_mash_clust_ward)

# Print cophenetic values
paste("Complete cophenetic correlation: ", cor(Dist(
  metag_mash_dist_matrix), metag_mash_clust_comp_coph))
paste("Average cophenetic correlation: ", cor(Dist(
  metag_mash_dist_matrix), metag_mash_clust_avg_coph))
paste("Single cophenetic correlation: ", cor(Dist(
  metag_mash_dist_matrix), metag_mash_clust_single_coph))
paste("Ward cophenetic correlation: ", cor(Dist(
  metag_mash_dist_matrix), metag_mash_clust_ward_coph))

# Based on cophenetic correlation values, average (UPGMA) is the best 
# hierarchical clustering algorithmn to use
metag_mash_clust_avg <- hclust(Dist(
  metag_mash_dist_matrix, method="correlation"),"average")

# Extract dendrogram data
metag_mash_clust_avg_dend <- as.dendrogram(metag_mash_clust_avg)
metag_mash_clust_avg_dend_data <- dendro_data(metag_mash_clust_avg_dend, 
                                               type="rectangle")

# Extract segment information
metag_mash_clust_avg_dend_data_segments <- 
  metag_mash_clust_avg_dend_data$segments

# Extract labels information and combine with metadata table ready for annotation
metag_mash_clust_avg_dend_data_info = metag_mash_clust_avg_dend_data$labels %>%
  mutate(Sample = label) %>%
  subset(., select=-c(label)) %>%
  left_join(metag_metadata, by="Sample")

# Extract only terminal segment information in order to colour according to project
metag_mash_clust_avg_dend_data_segments_ends <- 
  metag_mash_clust_avg_dend_data_segments %>%
  filter(yend == 0) %>% # filter for terminal dendrogram ends
  left_join(metag_mash_clust_avg_dend_data$labels, by = "x")# labels contains the row names from dist_matrix (i.e., sample_name)

# reformat data
colnames(metag_mash_clust_avg_dend_data_segments_ends) <- 
  gsub("label","Sample", colnames(metag_mash_clust_avg_dend_data_segments_ends))

# Combine tables together
metag_mash_clust_avg_dend_data_complete <- 
  left_join(metag_mash_clust_avg_dend_data_segments_ends, 
            metag_metadata, by="Sample")

## Plot Supplementary_Figure_S5
Supplementary_Figure_S5 <- 
  ggplot(data=metag_mash_clust_avg_dend_data_complete) + 
  geom_segment(data=metag_mash_clust_avg_dend_data_segments, 
               aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(aes(x=x, y=y.x, xend=xend, yend=yend, colour=Month), linewidth = 2) + 
  geom_point(aes(x=xend, y=-0.05, shape=Project), size=3) +
  scale_colour_brewer(palette = "Paired") + 
  scale_shape_manual(values=c(23,22,25,21)) + 
  #geom_text(data = fram_mags_sig_functional_cluster_avg_dend_data_complete, 
  #          aes(x=xend, y=yend, label=Species_rep), hjust=0, nudge_y = 0.03) + ##this line can be used to add MAG names to dendrogram segments 
  scale_y_reverse(limits = c(1.5, -0.1)) +
  coord_flip() + 
  labs(y = "MASH Distance", colour = "Month of sampling", shape = "Sampling campaign") + # flipped x and y coordinates for aesthetic reasons
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  guides(shape = guide_legend(override.aes = list(size = 4))) + 
  theme_bw() + 
  theme(legend.position="right", 
        legend.text=element_text(size=12, colour = "black"), 
        legend.title=element_text(size=14, colour = "black"), 
        axis.title.x = element_text(colour = "black", size = 14), 
        axis.text.x = element_text(colour = "black", size = 12), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

## Export
pdf(file="figures_output/Supplementary_Figure_S5.pdf", width=10, height=8)
Supplementary_Figure_S5
dev.off()


###################

#####

### Supplementary Figure S6

# Import SC-RBP data on number of genomes and species

rbp_genes_num_genomes_and_species = read.table("FRAM_WSC20_SC_RBP_num_genomes_and_species.txt", 
                                               header = TRUE, sep = "\t",as.is=TRUE)

rbp_genes_num_genomes_and_species



###################

#####

### Supplementary Figure S7

# Eukaryotic composition






###################

#####

### Supplementary Figure S8

# Import SC-RBP abundance data
rbp_species_abund=read.table("FRAM_WSC20_SC_RBP_L6_species_abund.txt", 
                             header = TRUE, sep = "\t", as.is=TRUE)

# Import SC-RBP taxonomy information
rbp_species_taxonomy=read.table("FRAM_WSC20_SC_RBP_L6_species_taxonomy.txt", 
                                header = TRUE, sep = "\t", as.is=TRUE)

# The rbp species abund file we imported already includes relative
# abundance. But for sanity check, we will recalculate 
rbp_species_rel_abund_long = rbp_species_abund %>%
  reshape2::dcast(Species~Sample, value.var="Raw_abund", data=.,) %>%
  replace(is.na(.),0) %>%
  unique() %>%
  tibble::column_to_rownames(., var="Species") %>%
  vegan::decostand(., MARGIN=2, method="total") %>%
  tibble::rownames_to_column(., var="Species") %>%
  reshape2::melt(., variable.name="Sample", value.name="Rel_abund")

# Combine abundance and taxonomy data and then generate a genus and class
# relative abundance table
rbp_family_and_class_rel_abund_long = rbp_species_rel_abund_long %>%
  left_join(rbp_species_taxonomy, by="Species") %>% 
  subset(., select=c(GTDB_Family,GTDB_Class,Sample,Rel_abund)) %>%
  unique() %>%
  aggregate(Rel_abund~GTDB_Family+GTDB_Class+Sample, data=., FUN=sum)

# Identify those genera that reach >1% relative abundance in at least 
# one sample
rbp_family_and_class_of_interest = rbp_family_and_class_rel_abund_long %>% 
  filter(., Rel_abund >= 0.04) %>% 
  subset(., select=c(GTDB_Family,GTDB_Class)) %>%
  filter(GTDB_Family != "Not_assigned") %>%
  unique() %>%
  arrange(GTDB_Class,GTDB_Family)

# Create long format relative abundance table with above-selected genera
rbp_family_and_class_of_interest_rel_abund_long = 
  rbp_family_and_class_rel_abund_long %>%
  filter(GTDB_Family %in% rbp_family_and_class_of_interest$GTDB_Family)

# Rather than grouping everything else under 'other', we will divide the other
# fraction based on the relative proportions of each class, to provide more
# information on the remaining community fraction

# Firstly, sum up relative abundances for classes captured in our 
# genera of interest dataframe
rbp_classes_of_interest_rel_abund = rbp_family_and_class_of_interest_rel_abund_long %>%
  aggregate(Rel_abund~GTDB_Class+Sample, data=., FUN=sum) %>%
  rename(., class_captured_abund = Rel_abund)

# Now calculate the per sample total relative abundance for each class,
# then combine with previous information and determine 'other' fraction at
# class and sample level and then combine with the genus of interest relative
# abundance dataframe
rbp_reps_family_of_interest_rel_abund_long_with_class_other = 
  rbp_family_and_class_rel_abund_long %>%
  aggregate(Rel_abund~GTDB_Class+Sample, data=., FUN=sum) %>%
  filter(GTDB_Class %in% rbp_classes_of_interest_rel_abund$GTDB_Class) %>%
  left_join(rbp_classes_of_interest_rel_abund, by=c("GTDB_Class","Sample")) %>%
  rename(class_total_rel_abund = Rel_abund) %>%
  mutate(Rel_abund = class_total_rel_abund-class_captured_abund) %>%
  mutate(GTDB_Family = paste0("z",GTDB_Class,"_other")) %>%
  subset(., select=c(GTDB_Family,Sample,Rel_abund,GTDB_Class)) %>%
  rbind(., rbp_family_and_class_of_interest_rel_abund_long) 

# Now determine the remaining community fraction not captured and label
# as 'Other/Not_assigned' 
rbp_family_of_interest_rel_abund_final_df = 
  rbp_reps_family_of_interest_rel_abund_long_with_class_other %>%
  aggregate(Rel_abund~Sample, data=., FUN=sum) %>%
  mutate(Rel_abund = 1-Rel_abund) %>%
  mutate(GTDB_Family = "zOther/Not_assigned") %>%
  mutate(GTDB_Class = "zOther/Not_assigned") %>%
  rbind(., rbp_reps_family_of_interest_rel_abund_long_with_class_other) %>%
  arrange(GTDB_Class,GTDB_Family)

### Process data for figure 4b - relative proportion of transcription of RBP_L6
# Repeat the above steps using the TPM input table rather
# than relative abundance input table

# Import SC-RBP TPM data
rbp_species_TPM=read.table("FRAM_WSC20_SC_RBP_L6_species_TPM.txt", 
                           header = TRUE, sep = "\t", as.is=TRUE)

# The rbp species abund file we imported already includes relative
# abundance. But for sanity check, we will recalculate 
rbp_species_rel_TPM_long = rbp_species_TPM %>%
  reshape2::dcast(Species~Sample, value.var="TPM", data=.,) %>%
  replace(is.na(.),0) %>%
  unique() %>%
  tibble::column_to_rownames(., var="Species") %>%
  vegan::decostand(., MARGIN=2, method="total") %>%
  tibble::rownames_to_column(., var="Species") %>%
  reshape2::melt(., variable.name="Sample", value.name="Rel_TPM")

# Combine abundance and taxonomy data and then generate a genus and class
# relative abundance table
rbp_family_and_class_rel_tpm_long = rbp_species_rel_TPM_long %>%
  left_join(rbp_species_taxonomy, by="Species") %>% 
  subset(., select=c(GTDB_Family,GTDB_Class,Sample,Rel_TPM)) %>%
  unique() %>%
  aggregate(Rel_TPM~GTDB_Family+GTDB_Class+Sample, data=., FUN=sum)

# Identify those genera that reach >1% relative abundance in at least 
# one sample
rbp_family_and_class_of_interest_tpm = rbp_family_and_class_rel_tpm_long %>% 
  filter(., Rel_TPM >= 0.04) %>% 
  subset(., select=c(GTDB_Family,GTDB_Class)) %>%
  filter(GTDB_Family != "Not_assigned") %>%
  unique() %>%
  arrange(GTDB_Class,GTDB_Family)

# Create long format relative abundance table with above-selected genera
rbp_family_and_class_of_interest_rel_tpm_long = 
  rbp_family_and_class_rel_tpm_long %>%
  filter(GTDB_Family %in% rbp_family_and_class_of_interest_tpm$GTDB_Family)

# Rather than grouping everything else under 'other', we will divide the other
# fraction based on the relative proportions of each class, to provide more
# information on the remaining community fraction

# Firstly, sum up relative abundances for classes captured in our 
# genera of interest dataframe
rbp_classes_of_interest_rel_tpm = rbp_family_and_class_of_interest_rel_tpm_long %>%
  aggregate(Rel_TPM~GTDB_Class+Sample, data=., FUN=sum) %>%
  rename(., class_captured_tpm = Rel_TPM)

# Now calculate the per sample total relative abundance for each class,
# then combine with previous information and determine 'other' fraction at
# class and sample level and then combine with the genus of interest relative
# abundance dataframe
rbp_reps_family_of_interest_rel_tpm_long_with_class_other = 
  rbp_family_and_class_rel_tpm_long %>%
  aggregate(Rel_TPM~GTDB_Class+Sample, data=., FUN=sum) %>%
  filter(GTDB_Class %in% rbp_classes_of_interest_rel_tpm$GTDB_Class) %>%
  left_join(rbp_classes_of_interest_rel_tpm, by=c("GTDB_Class","Sample")) %>%
  rename(class_total_rel_tpm = Rel_TPM) %>%
  mutate(Rel_TPM = class_total_rel_tpm-class_captured_tpm) %>%
  mutate(GTDB_Family = paste0("z",GTDB_Class,"_other")) %>%
  subset(., select=c(GTDB_Family,Sample,Rel_TPM,GTDB_Class)) %>%
  rbind(., rbp_family_and_class_of_interest_rel_tpm_long) 

# Now determine the remaining community fraction not captured and label
# as 'Other/Not_assigned' 
rbp_family_of_interest_rel_tpm_final_df = 
  rbp_reps_family_of_interest_rel_tpm_long_with_class_other %>%
  aggregate(Rel_TPM~Sample, data=., FUN=sum) %>%
  mutate(Rel_TPM = 1-Rel_TPM) %>%
  mutate(GTDB_Family = "zOther/Not_assigned") %>%
  mutate(GTDB_Class = "zOther/Not_assigned") %>%
  rbind(., rbp_reps_family_of_interest_rel_tpm_long_with_class_other) %>%
  arrange(GTDB_Class,GTDB_Family)

## Determine colour pallete
# Bar plots will be coloured based on family
# Combine list of family from reduced lists (rel abund and rel TPM) and set
# colours

taxa_list_1 <- rbp_family_of_interest_rel_tpm_final_df %>%
  subset(., select=c("GTDB_Family","GTDB_Class"))

taxa_list_2 <- rbp_family_of_interest_rel_abund_final_df %>%
  subset(., select=c("GTDB_Family","GTDB_Class"))

rbp_family_list_for_plotting = rbind(taxa_list_1,taxa_list_2) %>%
  arrange(GTDB_Class,GTDB_Family) %>%
  subset(., select="GTDB_Family") %>%
  unique()

rbind(taxa_list_1,taxa_list_2) %>%
  arrange(GTDB_Class,GTDB_Family)i
# How many colours are needed?
dim(rbp_family_list_for_plotting)

#define custom color scale
myColors <- brewer.pal(3, "Spectral")
names(myColors) <- levels(iris$Species)
custom_colors <- scale_colour_manual(name = "Species Names", values = myColors)

colour_list <- colorRampPalette(brewer.pal(12,"Paired"))(19)
names(colour_list) <- factor(c(rbp_family_list_for_plotting$GTDB_Family))
family_colour_palette <- scale_fill_manual(name = "Family", values = colour_list)

## Plot Supplementary_Figure_S8a
rbp_family_of_interest_rel_abund_final_df$Sample <- 
  factor(rbp_family_of_interest_rel_abund_final_df$Sample,
         levels = c("S10_SRF", "S10_BML",
                    "S8_SRF", "S8_BML",
                    "S25_SRF", "S25_BML",
                    "S26_SRF", "S26_BML"))

rbp_family_of_interest_rel_abund_final_df$GTDB_Family <- 
  factor(rbp_family_of_interest_rel_abund_final_df$GTDB_Family,
         levels = unique(c(rbp_family_of_interest_rel_abund_final_df$GTDB_Family)))

Supplementary_Figure_S8a <- ggplot(rbp_family_of_interest_rel_abund_final_df, 
                    aes(x=Sample, y=Rel_abund*100)) + 
  geom_bar(aes(fill=GTDB_Family), stat="identity", position="stack") + 
  family_colour_palette + 
  labs(colour = "Class") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12, colour = "black", angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")) + 
  guides(fill=guide_legend(override.aes = list(size=5), ncol=1))

## Plot Supplementary_Figure_S8b
rbp_family_of_interest_rel_tpm_final_df$Sample <- 
  factor(rbp_family_of_interest_rel_tpm_final_df$Sample,
         levels = c("S10_SRF", "S25_SRF", "S25_BML","S26_BML"))

rbp_family_of_interest_rel_tpm_final_df$GTDB_Family <- 
  factor(rbp_family_of_interest_rel_tpm_final_df$GTDB_Family,
         levels = unique(c(rbp_family_of_interest_rel_tpm_final_df$GTDB_Family)))

Supplementary_Figure_S8b <- ggplot(rbp_family_of_interest_rel_tpm_final_df, 
                                   aes(x=Sample, y=Rel_TPM*100)) + 
  geom_bar(aes(fill=GTDB_Family), stat="identity", position="stack") + 
  family_colour_palette + 
  labs(colour = "Class") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12, colour = "black", angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"))  + 
  guides(fill=guide_legend(override.aes = list(size=5), ncol=1))


### Combine Supplementary_Figure_S8a and Supplementary_Figure_S8b and export to PDF
# The colours were further optimised in Inkscape
library(patchwork)
pdf(file="figures_output/Supplementary_Figure_S8.pdf", height=8, width=16)
Supplementary_Figure_S8a|Supplementary_Figure_S8b+plot_layout(guides = "collect")
dev.off()

###################

#####

### Supplementary Figure S9
# Import files

# Import cazyme and taxa abundance data
comm_cazys_count_with_taxa_long = read.table("FRAM_WSC20_cazyme_genes_all_with_taxonomy.txt", 
                                             header = TRUE, sep = "\t",as.is=TRUE)

# Import the cazyme transcription data
comm_cazys_tpm_with_taxa_long = read.table("FRAM_WSC20_cazyme_genes_all_with_taxonomy_and_TPM.txt", 
                                           header = TRUE, sep = "\t",as.is=TRUE)

# Import metadata about samples
metagenome_meta = read.table("metagenome_meta.csv", 
                             header = TRUE, sep = ",",as.is=TRUE)

## Convert raw cazyme gene counts into bray-curtis dissimilarity matrix
# CAZyme gene counts will be normalised by the number of microbial genomes
# estimated in each sample to provide 'counts per microbial genome'
# Then they will be converted into relative proportions of sample total
# and then a bray curtis dissimilarity matrix
comm_cazys_norm_bray = comm_cazys_count_with_taxa_long %>%
  subset(., select=c(Sample,Gene,Count)) %>%
  aggregate(Count~Gene+Sample, data=., FUN=sum) %>%
  left_join(metagenome_rbp_data, by="Sample") %>%
  mutate(Norm_count = Count/Num_genomes) %>%
  dcast(Gene~Sample, value.var = "Norm_count", data=.) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var="Gene") %>%
  decostand(., method="total", MARGIN=2) %>%
  t() %>%
  vegdist(., method="bray", na.rm=F)

## The same process is repeated for the cazyme gene TPM values
comm_cazys_tpm_norm_bray = comm_cazys_tpm_with_taxa_long %>%
  subset(., select=c(Sample,Gene,TPM)) %>%
  aggregate(TPM~Gene+Sample, data=., FUN=sum) %>%
  left_join(metagenome_rbp_data, by="Sample") %>%
  mutate(Norm_tpm = TPM/Num_genomes) %>%
  dcast(Gene~Sample, value.var = "Norm_tpm", data=.) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var="Gene") %>%
  decostand(., method="total", MARGIN=2) %>%
  t() %>%
  vegdist(., method="bray")


## Perform hierarchical clustering analysis on metagenome composition
# We will test what is the best clustering method

# Test which clustering method best represents the variation in the data
# First perform hierarchical clustering with different algorithms
comm_cazys_norm_bray_comp <- hclust(comm_cazys_norm_bray, "complete")	
comm_cazys_norm_bray_avg<- hclust(comm_cazys_norm_bray,"average")
comm_cazys_norm_bray_sing <- hclust(comm_cazys_norm_bray,"single")
comm_cazys_norm_bray_ward <- hclust(comm_cazys_norm_bray,"ward.D")	

# Calculate cophenetic correlation for each method
comm_cazys_norm_bray_comp_coph <- cophenetic(comm_cazys_norm_bray_comp)
comm_cazys_norm_bray_avg_coph <- cophenetic(comm_cazys_norm_bray_avg)
comm_cazys_norm_bray_sing_coph <- cophenetic(comm_cazys_norm_bray_sing)
comm_cazys_norm_bray_ward_coph <- cophenetic(comm_cazys_norm_bray_ward)

# Print cophenetic values
paste("Complete cophenetic correlation: ", cor(comm_cazys_norm_bray, 
                                               comm_cazys_norm_bray_comp_coph))
paste("Average cophenetic correlation: ", cor(comm_cazys_norm_bray, 
                                              comm_cazys_norm_bray_avg_coph))
paste("Single cophenetic correlation: ", cor(comm_cazys_norm_bray, 
                                             comm_cazys_norm_bray_sing_coph))
paste("Ward cophenetic correlation: ", cor(comm_cazys_norm_bray, 
                                           comm_cazys_norm_bray_ward_coph))

# Based on cophenetic correlation values, average (UPGMA) is the best 
# hierarchical clustering algorithmn to use
comm_cazys_norm_bray_clust_avg <- hclust(comm_cazys_norm_bray,"average")

# Extract dendrogram data
comm_cazys_norm_bray_clust_avg_dend <- as.dendrogram(comm_cazys_norm_bray_clust_avg)
comm_cazys_norm_bray_clust_avg_dend_data <- dendro_data(comm_cazys_norm_bray_clust_avg_dend, 
                                              type="rectangle")

# Extract segment information
comm_cazys_norm_bray_clust_avg_dend_data_segments <- 
  comm_cazys_norm_bray_clust_avg_dend_data$segments

# Extract labels information and combine with metadata table ready for annotation
comm_cazys_norm_bray_clust_avg_dend_data_info = comm_cazys_norm_bray_clust_avg_dend_data$labels %>%
  mutate(Sample = label) %>%
  subset(., select=-c(label)) %>%
  left_join(metagenome_meta, by="Sample")

# Extract only terminal segment information in order to colour according to project
comm_cazys_norm_bray_clust_avg_dend_data_segments_ends <- 
  comm_cazys_norm_bray_clust_avg_dend_data_segments %>%
  filter(yend == 0) %>% # filter for terminal dendrogram ends
  left_join(comm_cazys_norm_bray_clust_avg_dend_data$labels, by = "x")# labels contains the row names from dist_matrix (i.e., sample_name)

# reformat data
colnames(comm_cazys_norm_bray_clust_avg_dend_data_segments_ends) <- 
  gsub("label","Sample", colnames(comm_cazys_norm_bray_clust_avg_dend_data_segments_ends))

# Combine tables together
comm_cazys_norm_bray_clust_avg_dend_data_complete <- 
  left_join(comm_cazys_norm_bray_clust_avg_dend_data_segments_ends, 
            metagenome_meta, by="Sample")

# Create colour palette for stations
station_colours <- c("S8" = "#004949", "S10" = "#ffb6db",
                     "S25" = "#6db6ff", "S26" = "#924900")

comm_cazys_norm_bray_clust_avg_dend_data_complete

## Plot Supplementary_Figure_S9a
Supplementary_Figure_S9a <- 
  ggplot(data=comm_cazys_norm_bray_clust_avg_dend_data_complete) + 
  geom_segment(data=comm_cazys_norm_bray_clust_avg_dend_data_segments, 
               aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(aes(x=x, y=y.x, xend=xend, yend=yend, colour=Station), linewidth = 2) + 
  geom_point(aes(x=xend, y=-0.02, shape=Depth_cat), size=3) +
  scale_colour_manual(values = station_colours) + 
  labs(y = "Height", colour = "Station", shape = "Depth") + # flipped x and y coordinates for aesthetic reasons
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  guides(shape = guide_legend(override.aes = list(size = 4))) + 
  theme_bw() + 
  theme(legend.position="right", 
        legend.text=element_text(size=12, colour = "black"), 
        legend.title=element_text(size=14, colour = "black"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10))
Supplementary_Figure_S9a

## Supplementary Figure S9b

## Repeat the hierarchical clustering steps and plotting for TPM data

# Test which clustering method best represents the variation in the data
# First perform hierarchical clustering with different algorithms
comm_cazys_tpm_norm_bray_comp <- hclust(comm_cazys_tpm_norm_bray, "complete")	
comm_cazys_tpm_norm_bray_avg<- hclust(comm_cazys_tpm_norm_bray,"average")
comm_cazys_tpm_norm_bray_sing <- hclust(comm_cazys_tpm_norm_bray,"single")
comm_cazys_tpm_norm_bray_ward <- hclust(comm_cazys_tpm_norm_bray,"ward.D")	

# Calculate cophenetic correlation for each method
comm_cazys_tpm_norm_bray_comp_coph <- cophenetic(comm_cazys_tpm_norm_bray_comp)
comm_cazys_tpm_norm_bray_avg_coph <- cophenetic(comm_cazys_tpm_norm_bray_avg)
comm_cazys_tpm_norm_bray_sing_coph <- cophenetic(comm_cazys_tpm_norm_bray_sing)
comm_cazys_tpm_norm_bray_ward_coph <- cophenetic(comm_cazys_tpm_norm_bray_ward)

# Print cophenetic values
paste("Complete cophenetic correlation: ", cor(comm_cazys_tpm_norm_bray, 
                                               comm_cazys_tpm_norm_bray_comp_coph))
paste("Average cophenetic correlation: ", cor(comm_cazys_tpm_norm_bray, 
                                              comm_cazys_tpm_norm_bray_avg_coph))
paste("Single cophenetic correlation: ", cor(comm_cazys_tpm_norm_bray, 
                                             comm_cazys_tpm_norm_bray_sing_coph))
paste("Ward cophenetic correlation: ", cor(comm_cazys_tpm_norm_bray, 
                                           comm_cazys_tpm_norm_bray_ward_coph))

# Based on cophenetic correlation values, average (UPGMA) is the best 
# hierarchical clustering algorithmn to use
comm_cazys_tpm_norm_bray_clust_avg <- hclust(comm_cazys_tpm_norm_bray,"average")

# Extract dendrogram data
comm_cazys_tpm_norm_bray_clust_avg_dend <- as.dendrogram(comm_cazys_tpm_norm_bray_clust_avg)
comm_cazys_tpm_norm_bray_clust_avg_dend_data <- dendro_data(comm_cazys_tpm_norm_bray_clust_avg_dend, 
                                                        type="rectangle")

# Extract segment information
comm_cazys_tpm_norm_bray_clust_avg_dend_data_segments <- 
  comm_cazys_tpm_norm_bray_clust_avg_dend_data$segments

# Extract labels information and combine with metadata table ready for annotation
comm_cazys_tpm_norm_bray_clust_avg_dend_data_info = comm_cazys_tpm_norm_bray_clust_avg_dend_data$labels %>%
  mutate(Sample = label) %>%
  subset(., select=-c(label)) %>%
  left_join(metagenome_meta, by="Sample")

# Extract only terminal segment information in order to colour according to project
comm_cazys_tpm_norm_bray_clust_avg_dend_data_segments_ends <- 
  comm_cazys_tpm_norm_bray_clust_avg_dend_data_segments %>%
  filter(yend == 0) %>% # filter for terminal dendrogram ends
  left_join(comm_cazys_tpm_norm_bray_clust_avg_dend_data$labels, by = "x")# labels contains the row names from dist_matrix (i.e., sample_name)

# reformat data
colnames(comm_cazys_tpm_norm_bray_clust_avg_dend_data_segments_ends) <- 
  gsub("label","Sample", colnames(comm_cazys_tpm_norm_bray_clust_avg_dend_data_segments_ends))

# Combine tables together
comm_cazys_tpm_norm_bray_clust_avg_dend_data_complete <- 
  left_join(comm_cazys_tpm_norm_bray_clust_avg_dend_data_segments_ends, 
            metagenome_meta, by="Sample")

# Create colour palette for stations
station_colours <- c("S8" = "#004949", "S10" = "#ffb6db",
                     "S25" = "#6db6ff", "S26" = "#924900")

## Plot Supplementary_Figure_S9a
Supplementary_Figure_S9b <- 
  ggplot(data=comm_cazys_tpm_norm_bray_clust_avg_dend_data_complete) + 
  geom_segment(data=comm_cazys_tpm_norm_bray_clust_avg_dend_data_segments, 
               aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(aes(x=x, y=y.x, xend=xend, yend=yend, colour=Station), linewidth = 2) + 
  geom_point(aes(x=xend, y=-0.02, shape=Depth_cat), size=3) +
  scale_colour_manual(values = station_colours) + 
  labs(y = "Height", colour = "Station", shape = "Depth") + # flipped x and y coordinates for aesthetic reasons
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  guides(shape = guide_legend(override.aes = list(size = 4))) + 
  theme_bw() + 
  theme(legend.position="right", 
        legend.text=element_text(size=12, colour = "black"), 
        legend.title=element_text(size=14, colour = "black"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10))
Supplementary_Figure_S9b

## Combine Supplementary_Figure_S9a and Supplementary_Figure_S9b and export to PDF
pdf(file="figures_output/Supplementary_Figure_S9.pdf", height=6, width=10)
Supplementary_Figure_S9a|Supplementary_Figure_S9b+plot_layout(guides = "collect")
dev.off()



###################

#####

### Supplementary Figure S10
# Import files

# Import cazyme and taxa abundance data
comm_cazys_count_with_taxa_long = read.table("FRAM_WSC20_cazyme_genes_all_with_taxonomy.txt", 
                                             header = TRUE, sep = "\t",as.is=TRUE)

# Import the cazyme transcription data
comm_cazys_tpm_with_taxa_long = read.table("FRAM_WSC20_cazyme_genes_all_with_taxonomy_and_TPM.txt", 
                                           header = TRUE, sep = "\t",as.is=TRUE)

# Import metadata about samples
metagenome_meta = read.table("metagenome_meta.csv", 
                             header = TRUE, sep = ",",as.is=TRUE)

# Import estimation of number of genomes from single-copy RBP genes
metagenome_rbp_data = read.table("FRAM_WSC20_SC_RBP_num_genomes_and_species.txt", 
                                 header = TRUE, sep = "\t",as.is=TRUE)

## Convert raw cazyme gene counts into normalized counts
# CAZyme gene counts will be normalised by the number of microbial genomes
# estimated in each sample to provide 'counts per microbial genome'
# Then they will be converted into relative proportions of sample total
comm_cazys_norm_long = comm_cazys_count_with_taxa_long %>%
  subset(., select=c(Sample,Gene,Count)) %>%
  aggregate(Count~Gene+Sample, data=., FUN=sum) %>%
  left_join(metagenome_rbp_data, by="Sample") %>%
  mutate(Norm_count = Count/Num_genomes) %>%
  dcast(Gene~Sample, value.var = "Norm_count", data=.) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var="Gene") %>%
  decostand(., method="total", MARGIN=2) %>%
  rownames_to_column(., var="Gene") %>%
  reshape2::melt(., variable.name="Sample", value.name="Norm_count")

## The same process is repeated for the cazyme gene TPM values
comm_cazys_tpm_norm_long = comm_cazys_tpm_with_taxa_long %>%
  subset(., select=c(Sample,Gene,TPM)) %>%
  aggregate(TPM~Gene+Sample, data=., FUN=sum) %>%
  left_join(metagenome_rbp_data, by="Sample") %>%
  mutate(Norm_tpm = TPM/Num_genomes) %>%
  dcast(Gene~Sample, value.var = "Norm_tpm", data=.) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var="Gene") %>%
  decostand(., method="total", MARGIN=2) %>%
  rownames_to_column(., var="Gene") %>%
  melt(., variable.name="Sample", value.name="Norm_tpm")

# Now we have both cazyme gene count and TPM data normalised sufficiently

# We want to see which cazyme gene families vary the most across sample
# in terms of transcription
# To do this we will determine the min, max and then the fold change across samples
# For this, we chose to exclude S25_BML sample when calculating the range, 
# as it is so different from the others, it masks changes
# In addition, we will select those with the largest change but also only those
# that reach >1% rel proportion of transcription in a sample
comm_cazys_large_range_tpm_selected_wide = comm_cazys_tpm_norm_long %>%
  dcast(Gene~Sample, data=., value.var="Norm_tpm") %>%
  nest(-Gene, -S25_BML) %>% # mask the columns you don't want to include
  mutate(min_norm_tpm = map(data, min), # determine minimum row value
         max_norm_tpm = map(data, max)) %>% # determine maximum row value
  unnest(cols = c(data, min_norm_tpm, max_norm_tpm)) %>% 
  as.data.frame() %>%
  filter(max_norm_tpm > 0.01) %>%
  mutate(change_norm_tpm = case_when(
    min_norm_tpm > 0 ~ max_norm_tpm/min_norm_tpm,
    TRUE ~ max_norm_tpm-min_norm_tpm)) %>%
  arrange(desc(change_norm_tpm)) %>%
  head(., 15)

# Above, we created a dataframe containing the 15 cazyme gene families with the 
# largest range in normalised TPM values across samples
# Let's create a list of the gene families of interest in the order that we want
# to plot them in (descending in terms of tpm value range size)
gene_families_of_interest_for_range_in_order = comm_cazys_large_range_tpm_selected_wide %>%
  arrange(desc(max_norm_tpm)) %>%
  subset(., select="Gene")

# Convert TPM dataframe to long format, tidy and merge with metadata info
# ready for plotting 
comm_cazys_large_range_tpm_selected_long = 
  comm_cazys_large_range_tpm_selected_wide %>% 
  subset(., select=-c(min_norm_tpm, max_norm_tpm, change_norm_tpm)) %>%
  melt(., variable.name="Sample", value.name="Norm_tpm") %>%
  arrange(match(Gene, gene_families_of_interest_for_range_in_order$Gene)) %>%
  left_join(metagenome_meta, by="Sample")

# Now let's create a dataframe containing the normalised counts for the same 
# gene families across the metagenome samples and combine with metadata
# info ready for plotting. Also we will remove the samples which we do not have 
# metatranscriptomes for
comm_cazys_large_range_count_selected_long = comm_cazys_norm_long %>%
  filter(Gene %in% gene_families_of_interest_for_range_in_order$Gene) %>% 
  filter(Sample != "S10_BML",
         Sample != "S8_SRF",
         Sample != "S8_BML",
         Sample != "S26_SRF") %>%
  subset(., select=c(Gene,Sample,Norm_count)) %>%
  arrange(match(Gene, gene_families_of_interest_for_range_in_order$Gene)) %>%
  left_join(metagenome_meta, by="Sample")

# Fix gene order
comm_cazys_large_range_count_selected_long$Gene <- 
  factor(comm_cazys_large_range_count_selected_long$Gene, 
         levels=unique(comm_cazys_large_range_count_selected_long$Gene))

# Create colour palette for stations
station_colours <- c("S8" = "#004949", "S10" = "#ffb6db",
                     "S25" = "#6db6ff", "S26" = "#924900")

## Plotting Supplementary_Figure_S10_part1
Supplementary_Figure_S10_part1 <- 
  ggplot(comm_cazys_large_range_count_selected_long,
         aes(x=ifelse(Norm_count>0,Norm_count*100,NA), y=Gene)) + 
  geom_boxplot(stat = "boxplot") +
  geom_point(aes(colour = Station, shape = Depth_cat), 
             size=3, position = position_dodge(width = 0)) + 
  scale_shape_manual(values=c(19,17)) + 
  scale_colour_manual(values = station_colours) + 
  scale_x_continuous(trans = 'log10', limits=c(0.1,NA),
                     breaks=c(0.1,1,10,100)) + 
  scale_y_discrete(limits=rev) + 
  labs(x = "Proportion of sample CAZyme gene count (%)",
       shape = "Depth layer") + 
  theme_bw() + 
  theme(legend.position = "left",
        legend.title = element_text(size=12, colour="black"),
        legend.text=element_text(size=12, colour="black"), 
        axis.text.x= element_text(size=10, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        axis.title.x = element_text(size=12, colour="black"),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white", color="black"), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(0.5,0.01,0.5,0.5, "cm")) + 
  guides(colour=guide_legend(override.aes = list(size=4), ncol=1)) + 
  guides(shape=guide_legend(override.aes = list(size=4), ncol=1))

# Now repeat but for the transcription
comm_cazys_large_range_tpm_selected_long$Gene <- factor(
  comm_cazys_large_range_tpm_selected_long$Gene, levels=unique(
    comm_cazys_large_range_tpm_selected_long$Gene))

Supplementary_Figure_S10_part2 <- 
  ggplot(comm_cazys_large_range_tpm_selected_long,
         aes(x=ifelse(Norm_tpm>0,Norm_tpm*100,NA), y=Gene)) + 
  geom_boxplot(stat = "boxplot") +
  geom_point(aes(colour = Station, shape = Depth_cat), 
             size=3, position = position_dodge(width = 0)) + 
  scale_colour_manual(values = station_colours) + 
  scale_shape_manual(values=c(19,17)) +
  scale_x_continuous(trans = 'log10', limits=c(0.1,NA),
                     breaks=c(0.01,0.1,1,10)) +
  scale_y_discrete(limits=rev) + 
  labs(x = "Proportion of sample CAZyme TPM (%)",
       shape = "Depth layer") + 
  theme_bw() + 
  theme(legend.position = "none",
        legend.title = element_text(size=12, colour="black"),
        legend.text=element_text(size=12, colour="black"), 
        axis.text.x= element_text(size=10, colour="black"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=12, colour="black"),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white", color="black"), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(0.5,0.5,0.5,0.01, "cm"))

## Supplementary_Figure_S10_part3
# We want to plot the composition of taxa for the 
# genes that were identified as varying in transcription a lot (in Figure 5)

# Firstly, create a dataframe with the normalised, proportional transcription
# for each of the gene families of interest
comm_cazys_large_range_genes_of_interest_sum_tpm_per_sample = 
  comm_cazys_large_range_tpm_selected_long %>%
  subset(., select=c(Gene,Sample,Norm_tpm)) %>%
  mutate(gene_total_tpm = Norm_tpm) %>%
  subset(., select=-c(Norm_tpm))

# Let's go back to the original dataframe and determine the rel
# proportion of transcription for each gene family in each sample by each genera
comm_cazys_norm_tpm_genus_long = comm_cazys_tpm_with_taxa_long %>%
  subset(., select=c(Sample,Gene,GTDB_Genus,GTDB_Class,TPM)) %>%
  aggregate(TPM~Gene+GTDB_Genus+GTDB_Class+Sample, data=., FUN=sum) %>%
  left_join(metagenome_rbp_data, by="Sample") %>%
  mutate(Norm_tpm = TPM/Num_genomes) %>%
  mutate(Gene_Genus_Class = paste0(Gene,"%%%",GTDB_Class, "%%%", GTDB_Genus)) %>%
  aggregate(Norm_tpm~Gene_Genus_Class+Sample, data=., FUN=sum) %>%
  dcast(Gene_Genus_Class~Sample, value.var = "Norm_tpm", data=.) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var="Gene_Genus_Class") %>%
  decostand(., method="total", MARGIN=2) %>%
  rownames_to_column(., var="Gene_Genus_Class") %>%
  separate(Gene_Genus_Class, c("Gene","GTDB_Class","GTDB_Genus"), "%%%") %>%
  melt(., id.vars=c("Gene","GTDB_Class","GTDB_Genus"), variable.name="Sample", value.name="Norm_tpm") 

# Subset the dataframe based on the gene families of interest 
# then determine per gene family relative proportion of transcription for each 
# genus
comm_cazys_large_range_tpm_norm_genes_of_interest_top_genera_long = comm_cazys_norm_tpm_genus_long %>%
  filter(Gene %in% comm_cazys_large_range_genes_of_interest_sum_tpm_per_sample$Gene) %>%
  left_join(comm_cazys_large_range_genes_of_interest_sum_tpm_per_sample,
            by=c("Gene","Sample")) %>%
  mutate(Norm_rel_tpm_per_fam = (Norm_tpm/gene_total_tpm)*100) %>%
  subset(., select=c(Gene,Sample,GTDB_Class,GTDB_Genus,Norm_rel_tpm_per_fam)) %>%
  filter(GTDB_Genus != "Not_assigned",
         GTDB_Genus != "Cutibacterium") %>%
  arrange(Gene,desc(Norm_rel_tpm_per_fam)) %>%
  filter(Norm_rel_tpm_per_fam != "NaN" &
           Norm_rel_tpm_per_fam > 0) %>%
  group_by(Gene,Sample) %>%
  slice_head(n=1) %>%
  as.data.frame() %>%
  arrange(match(Gene, gene_families_of_interest_for_range_in_order$Gene), GTDB_Class, GTDB_Genus)

comm_cazys_large_range_tpm_norm_genes_of_interest_top_genera_long %>%
  subset(., select=c(GTDB_Class,GTDB_Genus)) %>%
  unique() %>%
  arrange(GTDB_Class,GTDB_Genus)

# Fix factor order
comm_cazys_large_range_tpm_norm_genes_of_interest_top_genera_long$Gene <- 
  factor(comm_cazys_large_range_tpm_norm_genes_of_interest_top_genera_long$Gene, 
         levels = unique(comm_cazys_large_range_tpm_norm_genes_of_interest_top_genera_long$Gene))
comm_cazys_large_range_tpm_norm_genes_of_interest_top_genera_long$GTDB_Genus <- 
  factor(comm_cazys_large_range_tpm_norm_genes_of_interest_top_genera_long$GTDB_Genus, 
         levels = unique(comm_cazys_large_range_tpm_norm_genes_of_interest_top_genera_long$GTDB_Genus))
comm_cazys_large_range_tpm_norm_genes_of_interest_top_genera_long$Sample <- 
  factor(comm_cazys_large_range_tpm_norm_genes_of_interest_top_genera_long$Sample, 
         levels = unique(c("S10_SRF","S25_SRF","S25_BML","S26_BML")))

# Determine number of colours neeeded
cazy_genus_tpm_colourcount=
  length(unique(comm_cazys_large_range_tpm_norm_genes_of_interest_top_genera_long$GTDB_Genus))

## Plot Supplementary_Figure_S10_part3
Supplementary_Figure_S10_part3 <-  
  ggplot(comm_cazys_large_range_tpm_norm_genes_of_interest_top_genera_long, aes(x=Norm_rel_tpm_per_fam, y=Gene)) +
  geom_bar(aes(fill=GTDB_Genus), position="stack", stat="identity") + 
  scale_y_discrete(limits=rev) + 
  scale_x_continuous(limits=c(0,100)) + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(cazy_genus_tpm_colourcount)) + 
  facet_grid(.~Sample, scales = "free_x") +
  labs(x = "Relative proportion of per gene family TPM (%)", fill = "Genera") + 
  theme(legend.key=element_blank(),
        axis.text.x = element_text(colour = "black", size = 10), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
        legend.position = "right", 
        panel.grid.major.y = element_line(colour = "grey95"),
        panel.grid.major.x = element_line(colour = "grey95"),
        strip.background.x = element_rect(fill = "white", colour="black",
                                          linewidth = 1),
        strip.text.x = element_text(colour = "black", size = 14)) +
  guides(fill=guide_legend(override.aes = list(size=5), ncol=1))

### Combine two parts together and export to PDF
pdf(file="figures_output/Supplementary_Figure_S10.pdf", height=8, width=18)
Supplementary_Figure_S10_part1+Supplementary_Figure_S10_part2+
  Supplementary_Figure_S10_part3+plot_layout(ncol=4, widths=c(1,1,2))
dev.off()


###################

#####

### Supplementary Figure S11

# Import cazyme and taxa abundance data
comm_cazys_count_with_taxa_long = read.table("FRAM_WSC20_cazyme_genes_all_with_taxonomy.txt", 
                                             header = TRUE, sep = "\t",as.is=TRUE)

# Import the cazyme transcription data
comm_cazys_tpm_with_taxa_long = read.table("FRAM_WSC20_cazyme_genes_all_with_taxonomy_and_TPM.txt", 
                                           header = TRUE, sep = "\t",as.is=TRUE)

# Import metadata about samples
metagenome_meta = read.table("metagenome_meta.csv", 
                             header = TRUE, sep = ",",as.is=TRUE)

# Import estimation of number of genomes from single-copy RBP genes
metagenome_rbp_data = read.table("FRAM_WSC20_SC_RBP_num_genomes_and_species.txt", 
                                 header = TRUE, sep = "\t",as.is=TRUE)

## Convert raw cazyme gene counts into normalized counts
# CAZyme gene counts per genus will be normalised by the number of microbial genomes
# estimated in each sample to provide 'counts per microbial genome'
# Then they will be converted into relative proportions of sample total
comm_cazys_rel_prop_genera_long = comm_cazys_count_with_taxa_long %>%
  subset(., select=c(Sample,GTDB_Class,GTDB_Genus,Count)) %>%
  aggregate(Count~GTDB_Class+GTDB_Genus+Sample, data=., FUN=sum) %>%
  filter(GTDB_Genus != "Cutibacterium") %>%
  left_join(metagenome_rbp_data, by="Sample") %>%
  mutate(Norm_count = Count/Num_genomes) %>%
  mutate(Class_Genus = paste0(GTDB_Class,"%%%",GTDB_Genus)) %>%
  dcast(Class_Genus~Sample, value.var = "Norm_count", data=.) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var="Class_Genus") %>%
  decostand(., method="total", MARGIN=2) %>%
  rownames_to_column(., var="Class_Genus") %>%
  separate(Class_Genus, c("GTDB_Class","GTDB_Genus"), "%%%") %>%
  reshape2::melt(., id.vars=c("GTDB_Class","GTDB_Genus"), 
                 variable.name="Sample", value.name="Rel_prop")
  
# Identify genera who contribute >1% in any sample
# and then retrieve their relative proportions from above created dataframe
comm_cazys_genera_of_interest = comm_cazys_rel_prop_genera_long %>%
  filter(Rel_prop > 0.02) %>%
  filter(GTDB_Genus != "Not_assigned") %>%
  subset(., select=c(GTDB_Genus))
  
# Retrieve the genera of interest relative proportions for all samples
comm_cazys_rel_prop_top_genera_long = comm_cazys_rel_prop_genera_long %>%
  filter(GTDB_Genus %in% comm_cazys_genera_of_interest$GTDB_Genus)

# Sum up the relative proportions from the refined dataframe above at the 
# Class level
comm_cazys_rel_prop_top_class_long = 
  comm_cazys_rel_prop_top_genera_long %>%
  aggregate(Rel_prop~GTDB_Class+Sample, data=., FUN=sum) %>%
  rename(., class_rel_prop = Rel_prop)

# To provide more information, we will calculate the 'Other' fraction
# of the classes that were represented by the top contributing genera
# For this, we first calculate total rel prop for class in each sample
# then substract the proportion represented in the above created genera of 
# interest dataframe
comm_cazys_class_of_interest_other_rel_prop_long = 
  comm_cazys_rel_prop_genera_long %>%
  aggregate(Rel_prop~GTDB_Class+Sample, data=., FUN=sum) %>%
  filter(GTDB_Class %in% comm_cazys_rel_prop_top_genera_long$GTDB_Class) %>%
  rename(., class_total_rel_rel_prop = Rel_prop) %>%
  left_join(comm_cazys_rel_prop_top_class_long, by=c("GTDB_Class","Sample")) %>%
  replace(is.na(.),0) %>%
  mutate(Rel_prop = class_total_rel_rel_prop-class_rel_prop) %>%
  mutate(GTDB_Genus = paste0("z",GTDB_Class,"_other")) %>%
  subset(., select=c(GTDB_Class,GTDB_Genus,Sample,Rel_prop))

# Combine the dataframes together and merge with metadata ready for plotting
comm_cazys_genera_of_interest_for_plotting_long = 
  rbind(comm_cazys_class_of_interest_other_rel_prop_long,
      comm_cazys_rel_prop_top_genera_long) %>%
  arrange(GTDB_Class,GTDB_Genus) %>%
  left_join(metagenome_meta, by="Sample")

# Identify classes
comm_cazys_genera_of_interest_for_plotting_long %>%
  subset(., select=c(GTDB_Class)) %>%
  unique()

# Define class-level custom colour palette
# Here we will retain the same colours used for classes in other figures in the
# manuscript (e.g. Figure 4) and will add additional ones where necessary
class_colours <- c("Acidimicrobiia" = "#ffff6d", "Actinomycetia" = "#ff6db6", 
                   "Alphaproteobacteria" = "#6db6ff", "Bacteroidia" = "#db6d00", 
                   "Cyanobacteriia" = "#009292", "Gammaproteobacteria" = "#920000",
                   "Marinisomatia" = "#924900", "Phycisphaerae" = "#004949", 
                   "Poseidoniia" = "#490092", "Verrucomicrobiae" = "#b66dff",
                   "zOther/Not_assigned" = "#808080")

### Supplementary_Figure_S11a
# Set correct factor levels
comm_cazys_genera_of_interest_for_plotting_long$Sample <- 
  factor(comm_cazys_genera_of_interest_for_plotting_long$Sample,
         levels = c("S10_SRF", "S10_BML",
                    "S8_SRF", "S8_BML",
                    "S25_SRF", "S25_BML",
                    "S26_SRF", "S26_BML"))

comm_cazys_genera_of_interest_for_plotting_long$GTDB_Genus <- 
  factor(comm_cazys_genera_of_interest_for_plotting_long$GTDB_Genus,
         levels = unique(c(comm_cazys_genera_of_interest_for_plotting_long$GTDB_Genus)))

Supplementary_Figure_S11a <- ggplot(comm_cazys_genera_of_interest_for_plotting_long, 
                    aes(x=Sample, y=GTDB_Genus)) + 
  geom_point(aes(colour=GTDB_Class, size=Rel_prop*100), stat="identity", position="identity") + 
  scale_colour_manual(values=class_colours) +  
  scale_size(range = c(1, 10), breaks=c(0.5,1,2,4,8,16), name="Relative Abundance (%)") + 
  labs(colour = "Class") + 
  scale_y_discrete(limits=rev) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 11, colour = "black", angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")) + 
  guides(colour=guide_legend(override.aes = list(size=5), ncol=1))
Supplementary_Figure_S11a


### Supplementary Figure S11b
## The same process is repeated for the cazyme gene TPM values

# Normalize TPM values by number of microbial genomes and then convert
# to relative proportions 
comm_cazys_rel_tpm_genera_long = comm_cazys_tpm_with_taxa_long %>%
  subset(., select=c(Sample,GTDB_Class,GTDB_Genus,TPM)) %>%
  aggregate(TPM~GTDB_Class+GTDB_Genus+Sample, data=., FUN=sum) %>%
  filter(GTDB_Genus != "Cutibacterium") %>%
  left_join(metagenome_rbp_data, by="Sample") %>%
  mutate(Norm_TPM = TPM/Num_genomes) %>%
  mutate(Class_Genus = paste0(GTDB_Class,"%%%",GTDB_Genus)) %>%
  dcast(Class_Genus~Sample, value.var = "Norm_TPM", data=.) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var="Class_Genus") %>%
  decostand(., method="total", MARGIN=2) %>%
  rownames_to_column(., var="Class_Genus") %>%
  separate(Class_Genus, c("GTDB_Class","GTDB_Genus"), "%%%") %>%
  reshape2::melt(., id.vars=c("GTDB_Class","GTDB_Genus"), 
                 variable.name="Sample", value.name="Rel_TPM")

# Identify genera who contribute >1% in any sample
# and then retrieve their relative proportions from above created dataframe
comm_cazys_tpm_genera_of_interest = comm_cazys_rel_tpm_genera_long %>%
  filter(Rel_TPM > 0.02) %>%
  filter(GTDB_Genus != "Not_assigned") %>%
  subset(., select=c(GTDB_Genus))

# Retrieve the genera of interest relative proportions for all samples
comm_cazys_rel_tpm_top_genera_long = comm_cazys_rel_tpm_genera_long %>%
  filter(GTDB_Genus %in% comm_cazys_tpm_genera_of_interest$GTDB_Genus)

# Sum up the relative proportions from the refined dataframe above at the 
# Class level
comm_cazys_rel_tpm_top_class_long = 
  comm_cazys_rel_tpm_top_genera_long %>%
  aggregate(Rel_TPM~GTDB_Class+Sample, data=., FUN=sum) %>%
  rename(., class_rel_tpm = Rel_TPM)

# To provide more information, we will calculate the 'Other' fraction
# of the classes that were represented by the top contributing genera
# For this, we first calculate total rel prop for class in each sample
# then substract the proportion represented in the above created genera of 
# interest dataframe
comm_cazys_class_of_interest_other_rel_tpm_long = 
  comm_cazys_rel_tpm_genera_long %>%
  aggregate(Rel_TPM~GTDB_Class+Sample, data=., FUN=sum) %>%
  filter(GTDB_Class %in% comm_cazys_rel_tpm_top_genera_long$GTDB_Class) %>%
  rename(., class_total_rel_tpm = Rel_TPM) %>%
  left_join(comm_cazys_rel_tpm_top_class_long, by=c("GTDB_Class","Sample")) %>%
  replace(is.na(.),0) %>%
  mutate(Rel_TPM = class_total_rel_tpm-class_rel_tpm) %>%
  mutate(GTDB_Genus = paste0("z",GTDB_Class,"_other")) %>%
  subset(., select=c(GTDB_Class,GTDB_Genus,Sample,Rel_TPM))

# Combine the dataframes together and merge with metadata ready for plotting
comm_cazys_genera_of_interest_tpm_for_plotting_long = 
  rbind(comm_cazys_class_of_interest_other_rel_tpm_long,
        comm_cazys_rel_tpm_top_genera_long) %>%
  arrange(GTDB_Class,GTDB_Genus) %>%
  left_join(metagenome_meta, by="Sample")

# Identify classes to make sure no new levels need to be added to the 
# previously defined class colour palette
comm_cazys_genera_of_interest_tpm_for_plotting_long %>%
  subset(., select=c(GTDB_Class)) %>%
  unique()

### Supplementary_Figure_S11a
# Set correct factor levels
comm_cazys_genera_of_interest_tpm_for_plotting_long$Sample <- 
  factor(comm_cazys_genera_of_interest_tpm_for_plotting_long$Sample,
         levels = c("S10_SRF","S25_SRF", "S25_BML","S26_BML"))

comm_cazys_genera_of_interest_tpm_for_plotting_long$GTDB_Genus <- 
  factor(comm_cazys_genera_of_interest_tpm_for_plotting_long$GTDB_Genus,
         levels = unique(c(comm_cazys_genera_of_interest_tpm_for_plotting_long$GTDB_Genus)))

Supplementary_Figure_S11b <- ggplot(comm_cazys_genera_of_interest_tpm_for_plotting_long, 
                                    aes(x=Sample, y=GTDB_Genus)) + 
  geom_point(aes(colour=GTDB_Class, size=Rel_TPM*100), stat="identity", position="identity") + 
  scale_colour_manual(values=class_colours) +  
  scale_size(range = c(1, 10), breaks=c(0.5,1,2,4,8,16), name="Relative Abundance (%)") + 
  labs(colour = "Class") + 
  scale_y_discrete(limits=rev) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 11, colour = "black", angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")) + 
  guides(colour=guide_legend(override.aes = list(size=5), ncol=1))
Supplementary_Figure_S11b

### Combine two parts together and export to PDF
pdf(file="figures_output/Supplementary_Figure_S11.pdf", height=8, width=16)
Supplementary_Figure_S11a|Supplementary_Figure_S11b
dev.off()


###################

#####

### Supplementary Figure S12
library(ggpmisc)
library(ggrepel)

# Import cazyme and taxa abundance data
comm_cazys_count_with_taxa_long = read.table("FRAM_WSC20_cazyme_genes_all_with_taxonomy.txt", 
                                             header = TRUE, sep = "\t",as.is=TRUE)

# Import the cazyme transcription data
comm_cazys_tpm_with_taxa_long = read.table("FRAM_WSC20_cazyme_genes_all_with_taxonomy_and_TPM.txt", 
                                           header = TRUE, sep = "\t",as.is=TRUE)


## Convert raw cazyme gene counts into normalized counts
# CAZyme gene counts will be normalised by the number of microbial genomes
# estimated in each sample to provide 'counts per microbial genome'
# Then they will be converted into relative proportions of sample total
comm_cazys_norm_long = comm_cazys_count_with_taxa_long %>%
  subset(., select=c(Sample,Gene,Count)) %>%
  aggregate(Count~Gene+Sample, data=., FUN=sum) %>%
  left_join(metagenome_rbp_data, by="Sample") %>%
  mutate(Norm_count = Count/Num_genomes) %>%
  dcast(Gene~Sample, value.var = "Norm_count", data=.) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var="Gene") %>%
  decostand(., method="total", MARGIN=2) %>%
  rownames_to_column(., var="Gene") %>%
  reshape2::melt(., variable.name="Sample", value.name="Norm_count")

## The same process is repeated for the cazyme gene TPM values
comm_cazys_tpm_norm_long = comm_cazys_tpm_with_taxa_long %>%
  subset(., select=c(Sample,Gene,TPM)) %>%
  aggregate(TPM~Gene+Sample, data=., FUN=sum) %>%
  left_join(metagenome_rbp_data, by="Sample") %>%
  mutate(Norm_tpm = TPM/Num_genomes) %>%
  dcast(Gene~Sample, value.var = "Norm_tpm", data=.) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(., var="Gene") %>%
  decostand(., method="total", MARGIN=2) %>%
  rownames_to_column(., var="Gene") %>%
  melt(., variable.name="Sample", value.name="Norm_tpm")
comm_cazys_tpm_norm_long

comm_cazys_norm_and_tpm_norm_long = comm_cazys_tpm_norm_long %>%
  left_join(comm_cazys_norm_long, by=c("Gene","Sample"))

## Plot Supplementary Figure S12
switch=T
Supplementary_Figure_S12_base = ggplot(data=comm_cazys_norm_and_tpm_norm_long, aes(x=Norm_count*100,
                                                   y=Norm_tpm*100)) +
  facet_wrap(Sample~., scales="free") + 
  stat_poly_line() + 
  stat_poly_eq() + 
  geom_point(size = 2) +
  #{if(Switch)geom_text_repel(aes(label=ifelse(Norm_tpm>(Norm_count*2),Gene,NA)), 
  #                           stat="identity", nudge_y=0.75),}
  #geom_text_repel(aes(label=ifelse(Norm_tpm>(Norm_count*2),Gene,NA)), 
  #                stat="identity", nudge_y=0.75) + 
  labs(x = "Relative proportion of CAZyme gene count (%)",
       y = "Relative proportion of CAZyme gene transcription (%)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black"),
        strip.background = element_rect(fill = "white", colour = "black"), 
        strip.text.x = element_text(size = 12, colour = "black"))

# We are interested to see which of those are disproportionately transcribed
# compared to their abundance
# So filter out those whose relative transcription is more than double
# the relative proportion of gene count AND whose transcription is >1% of sample
# total
text_labels_to_add <- comm_cazys_norm_and_tpm_norm_long %>%
  filter(Norm_tpm > (Norm_count*2) & 
           Norm_tpm > 0.01)

# Add text labels
Supplementary_Figure_S12_complete = 
  Supplementary_Figure_S12 + geom_text_repel(data=text_labels_to_add,
                                           aes(x=Norm_count*100, y=Norm_tpm*100,
                                               label=Gene), nudge_y=2)

### Export to PDF
pdf(file="figures_output/Supplementary_Figure_S12.pdf", height=8, width=8)
Supplementary_Figure_S12_complete
dev.off()
