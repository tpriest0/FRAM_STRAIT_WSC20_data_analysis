################### 

##### Script encapsulating the processing, analysis and visualisation 
##### of data presented in the main body text of the paper  
##### paper "Spatial heterogeneity in carbohydrates and their utilisation by
##### microbial communities in the high North Atlantic"

###################

##### Using the script
### The script has been sectioned based on Figure number. The sections for 
### each figure is independent of all others. Although this does involve a lot
### replication of some processing steps, it means that you don't have to run
### the whole script from beginning to end if you are just interested in a 
### specific plot.
### You can use ctrl+F to jump to the figure section of interest

####################

#### Setting up

# Set this variable to download datafiles directory
datafiles_dir = ('data_files')

# Define working directory
setwd(datafiles_dir)

# Create directory ready for figure output
output_figures_dir = "figures_output"
output_figures_dir_path <- file.path(datafiles_dir, output_figures_dir)
dir.create(output_figures_dir_path)

####################

#####
# Figure 1
#####

### Plotting CTD data

# CTD profiles were plotted in R (below section of script) and then overlaid
# with the bathymetric map in Inkscape

# The following section of script that incorporates the plotting of 
# Temperature, salinity and fluorescence needs to be repeated for 
# each sampling station. 
# The result will be one output PDF file for each station
# Therefore, after each station, you need to modify the input filename, 
# output filename and change the y-intercept value for the 
# BML sampling depth (values provided below)

# This is the most cumbersome plotting but unfortunately comes first.
# If you are not interested in recreating the CTD profiles then skip onto the
# section for Figure 2.

# Load libraries
library(oce)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)

# Import the data using the read.ctd function from oce package
ctd=read.ctd("ctd_profiles/MSM95_10_1processE.cnv", columns=list(Average_sound_velocity=list(
  name="avgsvCM",unit=list(unit=expression(),scale="m/s"))))

# process the ctd file by only retaining the downcast and normalizing decimal places
ctd = ctd %>%
  ctdTrim(method = "downcast")%>%
  ctdDecimate(p = 0.2)

# Here we are first extracting the data information (measured variables) from the 
# ctd profile and then converting it to a dataframe.
# Next, we are extracting date, time, lat and long information from the metadata
# section of the ctd profile
# Then we are keeping only the measured variables of interest
ctd.tb = ctd@data %>%
  as.data.frame() %>%
  mutate(datetime =  ctd@metadata$startTime,
         longitude = ctd@metadata$longitude,
         latitude = ctd@metadata$latitude) %>%
  separate(datetime, c("date","time"), sep=" ") %>%
  select(date, time, longitude, latitude, pressure, temperature, salinity, fluorescence)

# For generating the CTD plots, we plot each of the three variables
# separately. The horizontal line that represents the BML sampling depth
# is included in the temperature profile plot. The y-intercept value needs
# to be changed for each station accordingly. Here are the necessary values
# S1 = 30
# S6 = 30
# S8 = 35
# S10 = 35
# S11 = 38
# S25 = 23
# S26 = 38
# S36 = 38
# S41 = 38
# S44 = 38

# Plot temperature profile
station_temperature = ggplot(data = ctd.tb%>%na.omit(), 
       aes(x = temperature, y = pressure)) +
  geom_path(col = "black", linewidth=1.5) +
  scale_y_reverse() +
  scale_x_continuous(position = "top") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 1),
        axis.title = element_text(size = 14, colour = 1)) +
  labs(x = expression(~Temperature~(degree~C)), y = "Depth (m)") + 
  geom_hline(yintercept = 35, colour = "purple", linewidth=1.5) ### This needs to be changed for each station - 
### it represents the BML sampling depth

# Plot salinity profile
station_salinity = ggplot(data = ctd.tb%>%na.omit(), 
       aes(x = salinity, y = pressure)) +
  geom_path(col = "blue", linewidth=1.5) +
  scale_y_reverse() +
  scale_x_continuous(position = "bottom") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 1),
        axis.title = element_text(size = 14, colour = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = expression(~Salinity~(psu)), y = "Depth (m)")

# Plot fluorescence profile
station_fluorescence = ggplot(data = ctd.tb%>%na.omit(), 
                          aes(x = fluorescence, y = pressure)) +
  geom_path(col = "darkolivegreen2", linewidth=1.5) +
  scale_y_reverse() +
  scale_x_continuous(position = "bottom") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 1),
        axis.title = element_text(size = 14, colour = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = expression(~Fluorescence~(mg~m3)), y = "Depth (m)")

# Combine the three plots into one grid and then export
# The three plots were then overlaid in inkscape to generate a single plot - 
# this is why we removed the background panel grid for salinity and fluorescence
# to make for easy dragging over the temperature plot.

pdf(file="figures_output/FRAM_WSC20_S10_ctd_profiles.pdf", height=4, width=12)
plot(station_temperature+station_salinity+station_fluorescence)
dev.off()

# Once a single plot for each station had been generated, these were laid
# out in a grid format and placed over the bathymetric map in Inkscape to 
# generate the final Figure 1

###################

#####
# Figure 2
#####

# Load libraries
library(dplyr)
library(reshape2)
library(ggplot2)
library(patchwork)
library(ggh4x)
library(RColorBrewer)
library(ggrepel)

### Figure 2a

# Import data
mono_data_raw = read.table("FRAM_WSC20_monosaccharides_r.txt",
                           header=TRUE, sep = "\t", as.is=TRUE)

# Import sample metadata dataframe
sample_metadata=read.table("sample_metadata.txt",
                           header=TRUE, sep = "\t", as.is=TRUE,
                           encoding="UTF-8")

# Reformat data to long format
mono_data_long = mono_data_raw %>%
  left_join(sample_metadata, by="Sample") %>%
  subset(., select=c("Sample","Station","Location","Depth_cat","Fucose",
                     "Galactosamine","Arabinose","Glucosamine","Galactose",
                     "Glucose","Xylose")) %>%
  reshape2::melt(., id.vars=c("Sample","Location","Station","Depth_cat"),
                 variable.name="Compound", value.name="Concentration")

# Calculate total monosaccharide conc. per sample
mono_data_sample_total_long = mono_data_long %>%
  aggregate(Concentration~Station+Depth_cat+Location, data=., FUN=sum) %>%
  arrange(Location, Station) %>%
  filter(Depth_cat != "m100", 
         Depth_cat != "m200")

# Define factor levels
mono_data_sample_total_long$Station <- 
  factor(mono_data_sample_total_long$Station, 
         levels=unique(c("S1","S6","S8","S11","S25","S26","S36","S41","S44")))

mono_data_sample_total_long$Depth_cat <- 
  factor(mono_data_sample_total_long$Depth_cat, 
         levels=unique(c("BML","SRF")))

# Plot Figure 2a
Figure_2a <- ggplot(mono_data_sample_total_long, aes(x=Concentration, y=Depth_cat)) +  
  geom_bar(fill = "#AAAAAA", stat="identity", position="identity") + 
  facet_nested(Location+Station ~ ., scales="fixed") + 
  labs(x = "Concentration (µg/L)", y = "Depth layer") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        strip.background.y = element_rect(fill = "white", colour = "black"),
        strip.text.y = element_text(size = 12, colour = "black"))
Figure_2a

### Figure 2b

# Reformat dataframe from wide to long and calculate total monosaccharide
# concentration per depth and location (in relation to shelf)
mono_data_depth_location_total_long = mono_data_long %>%
  filter(Depth_cat != "m100", 
         Depth_cat != "m200") %>%
  aggregate(Concentration~Depth_cat+Location, data=., FUN=sum) %>%
  rename(SUM_conc = Concentration)

# Reformat data to long, find the sum of each monosaccharide per station
# and depth and then combine with previous and calculate relative proportions
mono_data_rel_prop_per_depth_and_location = mono_data_long %>%
  filter(Depth_cat != "m100", 
         Depth_cat != "m200") %>%
  aggregate(Concentration~Compound+Depth_cat+Location, data=., FUN=sum) %>%
  left_join(mono_data_depth_location_total_long) %>%
  mutate(Rel_concentration = round((Concentration/SUM_conc)*100, 2)) %>% 
  mutate(label_ypos = cumsum(Rel_concentration) - 0.5*Rel_concentration)

# Define factor order
mono_data_rel_prop_per_depth_and_location$Depth_cat <- 
  factor(mono_data_rel_prop_per_depth_and_location$Depth_cat, 
         levels=unique(c("SRF","BML")))

# Define monosaccharide colour palette
mono_colour_palette <- c("Arabinose" = "#db6d00",
                         "Galactosamine" = "#490092",
                         "Galactose" = "#6db6ff",
                         "Glucosamine" = "#ffff6d",
                         "Glucose" = "#009292",
                         "Fucose" = "#ffb6db",
                         "Xylose" = "#920000")

# Plot Figure 2b
Figure_2b <- ggplot(mono_data_rel_prop_per_depth_and_location, 
       aes(x = " ", y = Rel_concentration, fill = Compound)) +
  geom_col(color = "black") + 
  geom_text_repel(aes(label = paste0(Rel_concentration, "%")), 
                position = position_stack(vjust = 0.5)) + 
  coord_polar(theta = "y", start = 0) + 
  scale_fill_manual(values=mono_colour_palette) +
  theme_bw() + 
  facet_nested_wrap(.~Location+Depth_cat) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, colour = "black"),
        strip.background.x = element_rect(fill = "white", colour = "black"),
        strip.text.x = element_text(size = 12, colour = "black"))
Figure_2b

# Export Figure 2
# NB : the value labels of the piecharts were further modified in Inkscape
pdf(file = "figures_output/Figure_2.pdf", heigh=6, width=12)
Figure_2a|Figure_2b
dev.off()

###################

#####
# Figure 3
#####

# Load libraries
# Import data from each extraction solvent, reformat and combine

# Import polysaccharide antibody signal data
polys_data_wide=read.table("FRAM_WSC20_polysaccharide_r_not_averaged.txt",
                           header=TRUE, sep = "\t", as.is=TRUE,
                           encoding="UTF-8", check.names = F)

# Import dataframe connecting antibody name to epitope
antibody_to_structure=read.table("antibodies_to_structures.txt",
                           header=TRUE, sep = "\t", as.is=TRUE,
                           encoding="UTF-8")

# Import sample metadata dataframe
sample_metadata=read.table("sample_metadata.txt",
                                 header=TRUE, sep = "\t", as.is=TRUE,
                                 encoding="UTF-8")

### Process data files
# Reformat polysaccharide data to long format and then determine mean 
# values for each antibody in each sample based off the four replicates
polys_data_mean_per_solvent_long =  polys_data_wide %>%
  reshape2::melt(., id.vars=c("Array_sample","Sample","Extraction_solvent"),
                        variable.name="Antibody", value.name="Signal") %>%
  aggregate(Signal~Sample+Extraction_solvent+Antibody, data=., FUN=mean) %>%
  filter(., !grepl("control", Sample))

# Export table with averaged signals
write.table(polys_data_mean_per_solvent_long, file="FRAM_WSC20_polysaccharide_r_averaged.txt",
            sep="\t")

# Determine total epitope signal for each sample by summing those from the three
# solvents and then combine dataframe with information on epitope structure
# for each antibody and combine with metadata
polys_data_total_per_sample_long = polys_data_mean_per_solvent_long %>%
  aggregate(Signal~Sample+Antibody, data=., FUN=sum) %>%
  left_join(antibody_to_structure, by="Antibody") %>%
  left_join(sample_metadata, by="Sample") %>%
  arrange(Location)

# As there are 17 epitopes in total, it is too many for discerning on a plot
# Therefore, we will filter the 12 with the highest abundance values 
# (antibody binding signals) to plot in the main manuscript figure
# (A complete plot will be made for supplementary figure)
antibodies_for_plotting = polys_data_total_per_sample_long %>%
  aggregate(Signal~Antibody+Epitope_antibody, data=., FUN=sum) %>%
  arrange(desc(Signal)) %>%
  head(., 12) %>%
  subset(., select="Antibody")

# Extract those antibodies from the above created dataframe ready for plotting
polys_data_total_per_sample_long_refined = polys_data_total_per_sample_long %>%
  filter(Antibody %in% antibodies_for_plotting$Antibody)

# Define factor orders
polys_data_total_per_sample_long_refined$Station <- 
  factor(polys_data_total_per_sample_long_refined$Station, 
         levels=unique(c("S1","S6","S8","S11","S25","S26","S36","S41","S44")))
polys_data_total_per_sample_long_refined$Depth_cat <- 
  factor(polys_data_total_per_sample_long_refined$Depth_cat, 
         levels=unique(c("SRF","BML")))

# Determine the frequency of detection of each epitope
polys_data_total_per_sample_long_refined %>%
  mutate(., presabs = case_when(
    Signal > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  aggregate(presabs~Epitope_antibody+Antibody, data=., FUN=sum) %>%
  mutate(det_freq = (presabs/22)*100)

# Plot barplot faceted by depth levels
Figure_3 = polys_data_total_per_sample_long_refined %>%
  filter(., Depth_cat == "SRF" | 
           Depth_cat == "BML") %>%
  ggplot(.,) +
  geom_bar(aes(x=ifelse(Depth_cat=="SRF",-Signal,NA), 
               y=Station,
               fill = Epitope_antibody), 
           stat = "identity", position = "identity",
           width = 1, colour = "black") +
  geom_bar(aes(x=ifelse(Depth_cat=="BML",Signal,NA), 
               y=Station,
               fill = Epitope_antibody), 
           stat = "identity", position = "identity",
           width = 1, colour = "black") +  
  scale_y_discrete(limits=rev) + 
  scale_fill_brewer(palette = "Paired") + 
  facet_wrap2(Epitope_antibody~., axes = "all", remove_labels = "all") + 
  labs(fill = "Epitope", x = "Epitope abundance (Antibody signal intensity)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(colour = "black", size = 10), 
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8), 
        legend.position = "bottom", 
        legend.text = element_text(colour = "black", size = 12),
        legend.title = element_text(colour = "black", size = 12),
        panel.grid.major.x = element_line(colour = "grey95"),
        panel.grid.major.y = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_blank()) + 
  guides(fill=guide_legend(ncol = 4))
Figure_3

# Export Figure 3
# Figure was further annotated, e.g. with "times of sampling" and inclusion of greek 
# letters for alpha and beta, in inkscape 
pdf(file="figures_output/Figure_3.pdf", height=8, width=10)
Figure_3
dev.off()


###################

library(ggplot2)
library(RColorBrewer)
library(dplyr)

#####
# Figure 4
#####

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
rbp_genus_and_class_rel_abund_long = rbp_species_rel_abund_long %>%
  left_join(rbp_species_taxonomy, by="Species") %>% 
  subset(., select=c(GTDB_Genus,GTDB_Class,Sample,Rel_abund)) %>%
  unique() %>%
  aggregate(Rel_abund~GTDB_Genus+GTDB_Class+Sample, data=., FUN=sum)
  
# Identify those genera that reach >1% relative abundance in at least 
# one sample
rbp_genus_and_class_of_interest = rbp_genus_and_class_rel_abund_long %>% 
  filter(., Rel_abund >= 0.025) %>% 
  subset(., select=c(GTDB_Genus,GTDB_Class)) %>%
  filter(GTDB_Genus != "Not_assigned") %>%
  unique() %>%
  arrange(GTDB_Class,GTDB_Genus)

# Create long format relative abundance table with above-selected genera
rbp_genus_and_class_of_interest_rel_abund_long = 
  rbp_genus_and_class_rel_abund_long %>%
  filter(GTDB_Genus %in% rbp_genus_and_class_of_interest$GTDB_Genus)

# Rather than grouping everything else under 'other', we will divide the other
# fraction based on the relative proportions of each class, to provide more
# information on the remaining community fraction

# Firstly, sum up relative abundances for classes captured in our 
# genera of interest dataframe
rbp_classes_of_interest_rel_abund = rbp_genus_and_class_of_interest_rel_abund_long %>%
  aggregate(Rel_abund~GTDB_Class+Sample, data=., FUN=sum) %>%
  rename(., class_captured_abund = Rel_abund)

# Now calculate the per sample total relative abundance for each class,
# then combine with previous information and determine 'other' fraction at
# class and sample level and then combine with the genus of interest relative
# abundance dataframe
rbp_reps_genus_of_interest_rel_abund_long_with_class_other = 
  rbp_genus_and_class_rel_abund_long %>%
  aggregate(Rel_abund~GTDB_Class+Sample, data=., FUN=sum) %>%
  filter(GTDB_Class %in% rbp_classes_of_interest_rel_abund$GTDB_Class) %>%
  left_join(rbp_classes_of_interest_rel_abund, by=c("GTDB_Class","Sample")) %>%
  rename(class_total_rel_abund = Rel_abund) %>%
  mutate(Rel_abund = class_total_rel_abund-class_captured_abund) %>%
  mutate(GTDB_Genus = paste0("z",GTDB_Class,"_other")) %>%
  subset(., select=c(GTDB_Genus,Sample,Rel_abund,GTDB_Class)) %>%
  rbind(., rbp_genus_and_class_of_interest_rel_abund_long) 

# Now determine the remaining community fraction not captured and label
# as 'Other/Not_assigned' 
rbp_genera_of_interest_rel_abund_final_df = 
  rbp_reps_genus_of_interest_rel_abund_long_with_class_other %>%
  aggregate(Rel_abund~Sample, data=., FUN=sum) %>%
  mutate(Rel_abund = 1-Rel_abund) %>%
  mutate(GTDB_Genus = "zOther/Not_assigned") %>%
  mutate(GTDB_Class = "zOther/Not_assigned") %>%
  rbind(., rbp_reps_genus_of_interest_rel_abund_long_with_class_other) %>%
  arrange(GTDB_Class,GTDB_Genus)

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
rbp_genus_and_class_rel_tpm_long = rbp_species_rel_TPM_long %>%
  left_join(rbp_species_taxonomy, by="Species") %>% 
  subset(., select=c(GTDB_Genus,GTDB_Class,Sample,Rel_TPM)) %>%
  unique() %>%
  aggregate(Rel_TPM~GTDB_Genus+GTDB_Class+Sample, data=., FUN=sum)

# Identify those genera that reach >1% relative abundance in at least 
# one sample
rbp_genus_and_class_of_interest_tpm = rbp_genus_and_class_rel_tpm_long %>% 
  filter(., Rel_TPM >= 0.025) %>% 
  subset(., select=c(GTDB_Genus,GTDB_Class)) %>%
  filter(GTDB_Genus != "Not_assigned") %>%
  unique() %>%
  arrange(GTDB_Class,GTDB_Genus)

# Create long format relative abundance table with above-selected genera
rbp_genus_and_class_of_interest_rel_tpm_long = 
  rbp_genus_and_class_rel_tpm_long %>%
  filter(GTDB_Genus %in% rbp_genus_and_class_of_interest_tpm$GTDB_Genus)

# Rather than grouping everything else under 'other', we will divide the other
# fraction based on the relative proportions of each class, to provide more
# information on the remaining community fraction

# Firstly, sum up relative abundances for classes captured in our 
# genera of interest dataframe
rbp_classes_of_interest_rel_tpm = rbp_genus_and_class_of_interest_rel_tpm_long %>%
  aggregate(Rel_TPM~GTDB_Class+Sample, data=., FUN=sum) %>%
  rename(., class_captured_tpm = Rel_TPM)

# Now calculate the per sample total relative abundance for each class,
# then combine with previous information and determine 'other' fraction at
# class and sample level and then combine with the genus of interest relative
# abundance dataframe
rbp_reps_genus_of_interest_rel_tpm_long_with_class_other = 
  rbp_genus_and_class_rel_tpm_long %>%
  aggregate(Rel_TPM~GTDB_Class+Sample, data=., FUN=sum) %>%
  filter(GTDB_Class %in% rbp_classes_of_interest_rel_tpm$GTDB_Class) %>%
  left_join(rbp_classes_of_interest_rel_tpm, by=c("GTDB_Class","Sample")) %>%
  rename(class_total_rel_tpm = Rel_TPM) %>%
  mutate(Rel_TPM = class_total_rel_tpm-class_captured_tpm) %>%
  mutate(GTDB_Genus = paste0("z",GTDB_Class,"_other")) %>%
  subset(., select=c(GTDB_Genus,Sample,Rel_TPM,GTDB_Class)) %>%
  rbind(., rbp_genus_and_class_of_interest_rel_tpm_long) 

# Now determine the remaining community fraction not captured and label
# as 'Other/Not_assigned' 
rbp_genera_of_interest_rel_tpm_final_df = 
  rbp_reps_genus_of_interest_rel_tpm_long_with_class_other %>%
  aggregate(Rel_TPM~Sample, data=., FUN=sum) %>%
  mutate(Rel_TPM = 1-Rel_TPM) %>%
  mutate(GTDB_Genus = "zOther/Not_assigned") %>%
  mutate(GTDB_Class = "zOther/Not_assigned") %>%
  rbind(., rbp_reps_genus_of_interest_rel_tpm_long_with_class_other) %>%
  arrange(GTDB_Class,GTDB_Genus)

### Determine colour pallete
# Bar plots will be coloured based on genera
# Combine list of genera from reduced lists (rel abund and rel TPM) and set
# colours#

taxa_list_1 <- rbp_genera_of_interest_rel_tpm_final_df %>%
  subset(., select=c("GTDB_Genus","GTDB_Class"))

taxa_list_2 <- rbp_genera_of_interest_rel_abund_final_df %>%
  subset(., select=c("GTDB_Genus","GTDB_Class"))

rbp_class_list_for_plotting = rbind(taxa_list_1,taxa_list_2) %>%
  arrange(GTDB_Class,GTDB_Genus) %>%
  subset(., select="GTDB_Class") %>%
  unique()

# View class list and create custom colour palette
head(rbp_class_list_for_plotting)

# Define custom palette
rbp_class_colours <- c("Actinomycetia" = "#ff6db6", "Alphaproteobacteria" = "#6db6ff",
                       "Bacteroidia" = "#db6d00", "Cyanobacteriia" = "#009292",
                       "Gammaproteobacteria" = "#920000", "Poseidoniia" = "#490092",
                       "Verrucomicrobiae" = "#b66dff",
                       "zOther/Not_assigned" = "#808080")

### Plot Figure 4a
# Set correct factor levels
rbp_genera_of_interest_rel_abund_final_df$Sample <- 
  factor(rbp_genera_of_interest_rel_abund_final_df$Sample,
         levels = c("S10_SRF", "S10_BML",
                    "S8_SRF", "S8_BML",
                    "S25_SRF", "S25_BML",
                    "S26_SRF", "S26_BML"))

rbp_genera_of_interest_rel_abund_final_df$GTDB_Genus <- 
  factor(rbp_genera_of_interest_rel_abund_final_df$GTDB_Genus,
         levels = unique(c(rbp_genera_of_interest_rel_abund_final_df$GTDB_Genus)))

# Plot figure 4a
Figure_4a <- ggplot(rbp_genera_of_interest_rel_abund_final_df, 
                   aes(x=Sample, y=GTDB_Genus)) + 
  geom_point(aes(colour=GTDB_Class, size=Rel_abund*100), stat="identity", position="identity") + 
  scale_colour_manual(values=rbp_class_colours) +  
  scale_size(range = c(1, 10), breaks=c(0.5,1,2,4,8,16), name="Relative Abundance (%)") + 
  labs(colour = "Class") + 
  scale_y_discrete(limits=rev) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, colour = "black", angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")) + 
  guides(colour=guide_legend(override.aes = list(size=5), ncol=1))

### Plot Figure 4b
rbp_genera_of_interest_rel_tpm_final_df$Sample <- 
  factor(rbp_genera_of_interest_rel_tpm_final_df$Sample,
         levels = c("S10_SRF", "S25_SRF", "S25_BML","S26_BML"))

rbp_genera_of_interest_rel_tpm_final_df$GTDB_Genus <- 
  factor(rbp_genera_of_interest_rel_tpm_final_df$GTDB_Genus,
         levels = unique(c(rbp_genera_of_interest_rel_tpm_final_df$GTDB_Genus)))

# Plot figure 4b
Figure_4b <- ggplot(rbp_genera_of_interest_rel_tpm_final_df, 
       aes(x=Sample, y=GTDB_Genus)) + 
  geom_point(aes(colour=GTDB_Class, size=Rel_TPM*100), stat="identity", position="identity") + 
  scale_colour_manual(values=rbp_class_colours) +  
  scale_size(range = c(1, 10), breaks=c(0.5,1,2,4,8,16), name="Relative proportion of TPM (%)") + 
  labs(colour = "Class") + 
  scale_y_discrete(limits=rev) + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 10, colour = "black", angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")) + 
  guides(colour=guide_legend(override.aes = list(size=5), ncol=1))


### Combine Figure 4a and Figure4b and export to PDF
# The colours were further optimised in Inkscape
library(patchwork)
pdf(file="figures_output/Figure_4.pdf", height=8, width=16)
Figure_4a|Figure_4b+plot_layout(guides = "collect")
dev.off()

###################



###################

#####
# Figure 5
#####
library(reshape2)
library(tibble)
library(vegan)
library(tidyverse)
library(tidyr)
library(RColorBrewer)

# Import files

# Import cazyme and taxa abundance data
comm_cazys_count_with_taxa_long = read.table("FRAM_WSC20_cazyme_genes_all_with_taxonomy.txt", 
                                     header = TRUE, sep = "\t",as.is=TRUE)

# Import the cazyme transcription data
comm_cazys_tpm_with_taxa_long = read.table("FRAM_WSC20_cazyme_genes_all_with_taxonomy_and_TPM.txt", 
                                     header = TRUE, sep = "\t",as.is=TRUE)

# Import metadata about samples
sample_metadata = read.table("sample_metadata.txt", 
                             header = TRUE, sep = "\t",as.is=TRUE)

# Import estimation of number of genomes from single-copy RBP genes
metagenome_rbp_data = read.table("FRAM_WSC20_SC_RBP_num_genomes_and_species.txt", 
                                 header = TRUE, sep = "\t",as.is=TRUE)

## Create wide format cazyme count unnormalised dataframe
# This table is provided as supplementary table S6 in the manuscript
comm_cazys_count_wide = comm_cazys_count_with_taxa_long %>%
  subset(., select=c(Sample,Gene,Count)) %>%
  aggregate(Count~Gene+Sample, data=., FUN=sum) %>%
  dcast(Gene~Sample, value.var = "Count", data=.)

write.table(comm_cazys_count_wide, file="FRAM_WSC20_cazyme_genes_count_wide.txt",
            sep="\t")

## Create wide format cazyme TPM unnormalised dataframe
# This table is provided as supplementary table S7 in the manuscript
comm_cazys_tpm_wide = comm_cazys_tpm_with_taxa_long %>%
  subset(., select=c(Sample,Gene,TPM)) %>%
  aggregate(TPM~Gene+Sample, data=., FUN=sum) %>%
  dcast(Gene~Sample, value.var = "TPM", data=.)

write.table(comm_cazys_tpm_wide, file="FRAM_WSC20_cazyme_genes_TPM_wide.txt",
            sep="\t")

## Convert raw cazyme gene counts into normalized counts
# CAZyme gene counts will be normalised by the number of microbial genomes
# estimated in each sample to provide 'counts per microbial genome'
# Then they will be converted into relative proportions of sample total
comm_cazys_norm_wide = comm_cazys_count_with_taxa_long %>%
  subset(., select=c(Sample,Gene,Count)) %>%
  aggregate(Count~Gene+Sample, data=., FUN=sum) %>%
  left_join(metagenome_rbp_data, by="Sample") %>%
  mutate(Norm_count = Count/Num_genomes) %>%
  dcast(Gene~Sample, value.var = "Norm_count", data=.) %>%
  replace(is.na(.), 0) %>%
  tibble::column_to_rownames(., var="Gene") %>%
  decostand(., method="total", MARGIN=2) %>%
  rownames_to_column(., var="Gene")

# Convert the table to long format
comm_cazys_norm_long = comm_cazys_norm_wide %>%
  reshape2::melt(., variable.name="Sample", value.name="Norm_count")

## The same process is repeated for the cazyme gene TPM values
comm_cazys_tpm_norm_wide = comm_cazys_tpm_with_taxa_long %>%
  subset(., select=c(Sample,Gene,TPM)) %>%
  aggregate(TPM~Gene+Sample, data=., FUN=sum) %>%
  left_join(metagenome_rbp_data, by="Sample") %>%
  mutate(Norm_tpm = TPM/Num_genomes) %>%
  dcast(Gene~Sample, value.var = "Norm_tpm", data=.) %>%
  replace(is.na(.), 0) %>%
  tibble::column_to_rownames(., var="Gene") %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="Gene")

# Convert the table to long format
comm_cazys_tpm_norm_long = comm_cazys_tpm_norm_wide %>%
  melt(., variable.name="Sample", value.name="Norm_tpm")

# Now we have both cazyme gene count and TPM data normalised sufficiently
# We would like to visualise those that are transcribed the most across
# samples and additionaly see their abundance in the metagenomes and the 
# taxonomic group contributing the most to transcription in each sample
# As sample S25_BML is so different from the others, we will mask it 
# when determining the maximum transcribed gene families
comm_cazys_most_transcribed_selected_wide = comm_cazys_tpm_norm_long %>%
  dcast(Gene~Sample, data=., value.var="Norm_tpm") %>%
  nest(-Gene, -S25_BML) %>% # mask the columns you don't want to include
  mutate(max_norm_tpm = map(data, max)) %>% # determine maximum row value
  unnest(cols = c(data, max_norm_tpm)) %>% 
  as.data.frame() %>%
  arrange(desc(max_norm_tpm)) %>%
  head(., 15)

# Create a list of the gene families of interest in the order that we want
# to plot them in (descending in terms of normalised, relative transcription)
gene_families_of_interest_in_order = comm_cazys_most_transcribed_selected_wide %>%
  arrange(desc(max_norm_tpm)) %>%
  subset(., select="Gene")

# Convert TPM dataframe to long format, tidy and merge with metadata info
# ready for plotting 
comm_cazys_most_transcribed_tpm_selected_long = 
  comm_cazys_most_transcribed_selected_wide %>% 
  subset(., select=-c(max_norm_tpm)) %>%
  melt(., variable.name="Sample", value.name="Norm_tpm") %>%
  arrange(match(Gene, gene_families_of_interest_in_order$Gene)) %>%
  left_join(sample_metadata, by="Sample")

# Now let's create a dataframe containing the normalised counts for the same 
# gene families across the metagenome samples and combine with metadata
# info ready for plotting. Also we will remove the samples which we do not have 
# metatranscriptomes for
comm_cazys_selected_norm_count_long = comm_cazys_norm_long %>%
  filter(Gene %in% gene_families_of_interest_in_order$Gene) %>% 
  filter(Sample != "S10_BML",
         Sample != "S8_SRF",
         Sample != "S8_BML",
         Sample != "S26_SRF") %>%
  subset(., select=c(Gene,Sample,Norm_count)) %>%
  arrange(match(Gene, gene_families_of_interest_in_order$Gene)) %>%
  left_join(sample_metadata, by="Sample")

# Fix gene order
comm_cazys_selected_norm_count_long$Gene <- 
  factor(comm_cazys_selected_norm_count_long$Gene, 
         levels=unique(comm_cazys_selected_norm_count_long$Gene))

# Create colour palette for stations
station_colours <- c("S8" = "#004949", "S10" = "#ffb6db",
                     "S25" = "#6db6ff", "S26" = "#924900")

### Plotting Figure 5a
# Plot will be performed in two halves, first normalised counts
Figure_5a_part1 <- 
  ggplot(comm_cazys_selected_norm_count_long,
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
comm_cazys_most_transcribed_tpm_selected_long$Gene <- factor(
  comm_cazys_most_transcribed_tpm_selected_long$Gene, levels=unique(
    comm_cazys_most_transcribed_tpm_selected_long$Gene))

Figure_5a_part2 <- 
  ggplot(comm_cazys_most_transcribed_tpm_selected_long,
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

### Figure_5a_part3
# We want to plot the composition of taxa for the 
# genes that were identified as varying in transcription a lot (in Figure 5)

# Firstly, create a dataframe with the normalised, proportional transcription
# for each of the gene families of interest
comm_cazys_genes_of_interest_sum_tpm_per_sample = 
  comm_cazys_most_transcribed_tpm_selected_long %>%
  subset(., select=c(Gene,Sample,Norm_tpm)) %>%
  mutate(gene_total_tpm = Norm_tpm) %>%
  subset(., select=-c(Norm_tpm))

# Let's go back to the original dataframe and determine the relative
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
  tibble::column_to_rownames(., var="Gene_Genus_Class") %>%
  decostand(., method="total", MARGIN=2) %>%
  tibble::rownames_to_column(., var="Gene_Genus_Class") %>%
  separate(Gene_Genus_Class, c("Gene","GTDB_Class","GTDB_Genus"), "%%%") %>%
  melt(., id.vars=c("Gene","GTDB_Class","GTDB_Genus"), variable.name="Sample", value.name="Norm_tpm") 

# Subset the dataframe based on the gene families of interest 
# then determine per gene family relative proportion of transcription for each 
# genus
comm_cazys_tpm_norm_genes_of_interest_top_genera_long = comm_cazys_norm_tpm_genus_long %>%
  filter(Gene %in% comm_cazys_genes_of_interest_sum_tpm_per_sample$Gene) %>%
  left_join(comm_cazys_genes_of_interest_sum_tpm_per_sample,
            by=c("Gene","Sample")) %>%
  mutate(Norm_rel_tpm_per_fam = (Norm_tpm/gene_total_tpm)*100) %>%
  subset(., select=c(Gene,Sample,GTDB_Class,GTDB_Genus,Norm_rel_tpm_per_fam)) %>%
  filter(GTDB_Genus != "Not_assigned",
         GTDB_Genus != "Cutibacterium") %>%
  arrange(Gene,desc(Norm_rel_tpm_per_fam)) %>%
  filter(Norm_rel_tpm_per_fam != "NaN" &
           Norm_rel_tpm_per_fam > 0) %>%
  group_by(Gene,Sample) %>%
  slice_head(n=2) %>%
  as.data.frame() %>%
  arrange(match(Gene, gene_families_of_interest_in_order$Gene), GTDB_Class, GTDB_Genus)

# Fix factor order
comm_cazys_tpm_norm_genes_of_interest_top_genera_long$Gene <- 
  factor(comm_cazys_tpm_norm_genes_of_interest_top_genera_long$Gene, 
         levels = unique(comm_cazys_tpm_norm_genes_of_interest_top_genera_long$Gene))
comm_cazys_tpm_norm_genes_of_interest_top_genera_long$GTDB_Genus <- 
  factor(comm_cazys_tpm_norm_genes_of_interest_top_genera_long$GTDB_Genus, 
         levels = unique(comm_cazys_tpm_norm_genes_of_interest_top_genera_long$GTDB_Genus))
comm_cazys_tpm_norm_genes_of_interest_top_genera_long$Sample <- 
  factor(comm_cazys_tpm_norm_genes_of_interest_top_genera_long$Sample, 
         levels = unique(c("S10_SRF","S25_SRF","S25_BML","S26_BML")))

# Determine number of colours neeeded
cazy_genus_tpm_colourcount=
  length(unique(comm_cazys_tpm_norm_genes_of_interest_top_genera_long$GTDB_Genus))

# Plot Figure 5b
Figure_5a_part3 <-  
  ggplot(comm_cazys_tpm_norm_genes_of_interest_top_genera_long, aes(x=Norm_rel_tpm_per_fam, y=Gene)) +
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
  guides(fill=guide_legend(override.aes = list(size=5), ncol=2))

Figure_5a_part3

### Combine two parts together and export to PDF
library(patchwork)
pdf(file="figures_output/Figure_5.pdf", height=6, width=18)
Figure_5a_part1+Figure_5a_part2+Figure_5a_part3+plot_layout(ncol=4, widths=c(1,1,2))
dev.off()

###################

###################

#####
# Figure 6
#####


## Figure 6 consists of a phylogenetic tree along with relative abundance
## and relative proportion of transcription information for representative MAGs 
## The Tree was produced using FastTree and beautified in iToL. 
## Below is the processing of the data to generate relative abundance
## and relative proportion of transcription


###################

###################

#####
# Figure 7
#####

# Import new libraries
library(stringr)

mag_summary_stats = read.table("FRAM_WSC20_MAG_reps_summary_statistics_and_taxonomy.csv", 
                               header = TRUE, sep = ",", as.is=TRUE)

mag_reps_cazyme_genes_with_tpm = read.table("FRAM_WSC20_MAG_reps_cazyme_genes_with_TPM.txt", 
                                                header = TRUE, sep = "\t", as.is=TRUE)

mag_reps_peptidase_genes_with_tpm = read.table("FRAM_WSC20_MAG_reps_peptidase_genes_with_TPM.txt", 
                                            header = TRUE, sep = "\t", as.is=TRUE)

mag_reps_sulfatase_genes_with_tpm = read.table("FRAM_WSC20_MAG_reps_sulfatase_genes_with_TPM.txt", 
                                            header = TRUE, sep = "\t", as.is=TRUE)

mag_reps_tbdt_genes_with_tpm = read.table("FRAM_WSC20_MAG_reps_TBDT_genes_with_TPM.txt", 
                                            header = TRUE, sep = "\t", as.is=TRUE)

# Firstly, we will determine the number of genes in each group for each MAG
mag_reps_cazyme_gene_count = mag_reps_cazyme_genes_with_tpm %>%
  subset(., select=c(MAG_name)) %>%
  group_by(MAG_name) %>%
  count(., name="Count") %>%
  as.data.frame(.) %>%
  mutate(Gene_group="CAZymes")

mag_reps_peptidase_gene_count = mag_reps_peptidase_genes_with_tpm %>%
  subset(., select=c(MAG_name)) %>%
  group_by(MAG_name) %>%
  count(., name="Count") %>%
  as.data.frame(.) %>%
  mutate(Gene_group="Peptidases")

mag_reps_sulfatase_gene_count = mag_reps_sulfatase_genes_with_tpm %>%
  subset(., select=c(MAG_name)) %>%
  group_by(MAG_name) %>%
  count(., name="Count") %>%
  as.data.frame(.) %>%
  mutate(Gene_group="Sulfatases")

mag_reps_tbdt_gene_count = mag_reps_tbdt_genes_with_tpm %>%
  subset(., select=c(MAG_name)) %>%
  group_by(MAG_name) %>%
  count(., name="Count") %>%
  as.data.frame(.) %>%
  mutate(Gene_group="TBDTs")

# Combine the counts of each gene group and then divide by the genome size in Mbp
# To get counts per Mbp
mag_reps_gene_counts_long <- rbind(mag_reps_cazyme_gene_count,mag_reps_peptidase_gene_count,
      mag_reps_sulfatase_gene_count,mag_reps_tbdt_gene_count) %>%
  left_join(mag_summary_stats, by="MAG_name") %>%
  mutate(Value=Count/(Genome_size/1000000)) %>%
  subset(., select=c(MAG_name,Gene_group,Value)) %>%
  mutate(Sample="Count_per_Mbp")

# Now we will determine the total TPM for each gene group for each MAG
# in each sample
# First CAZymes
mag_reps_cazyme_gene_tpm_per_sample = mag_reps_cazyme_genes_with_tpm %>%
  subset(., select=-c(Gene_name)) %>%
  reshape2::melt(., id.vars=c("MAG_name","Gene"), variable.name="Sample", 
                 value.name="Value") %>%
  subset(., select=-c(Gene)) %>%
  group_by(MAG_name,Sample) %>%
  summarise(., across(everything(), sum)) %>%
  as.data.frame(.) %>%
  mutate(Gene_group="CAZymes")

# Repeat for peptidases
mag_reps_peptidase_gene_tpm_per_sample = mag_reps_peptidase_genes_with_tpm %>%
  subset(., select=-c(Gene_name,Gene_1)) %>%
  reshape2::melt(., id.vars=c("MAG_name","Gene_2"), variable.name="Sample", 
                 value.name="Value") %>%
  subset(., select=-c(Gene_2)) %>%
  group_by(MAG_name,Sample) %>%
  summarise(., across(everything(), sum)) %>%
  as.data.frame(.) %>%
  mutate(Gene_group="Peptidases")

# Repeat for sulfatases
mag_reps_sulfatase_gene_tpm_per_sample = mag_reps_sulfatase_genes_with_tpm %>%
  subset(., select=-c(Gene_name,Gene_1)) %>%
  reshape2::melt(., id.vars=c("MAG_name","Gene_2"), variable.name="Sample", 
                 value.name="Value") %>%
  subset(., select=-c(Gene_2)) %>%
  group_by(MAG_name,Sample) %>%
  summarise(., across(everything(), sum)) %>%
  as.data.frame(.) %>%
  mutate(Gene_group="Sulfatases")

# Repeat for TBDTs
mag_reps_tbdt_gene_tpm_per_sample = mag_reps_tbdt_genes_with_tpm %>%
  subset(., select=-c(Gene_name)) %>%
  reshape2::melt(., id.vars=c("MAG_name","Gene"), variable.name="Sample", 
                 value.name="Value") %>%
  subset(., select=-c(Gene)) %>%
  group_by(MAG_name,Sample) %>%
  summarise(., across(everything(), sum)) %>%
  as.data.frame(.) %>%
  mutate(Gene_group="TBDTs")

# Combine TPM dataframes and the gene count dataframe generated
# before. Then add taxonomic information to the table from the mag summary statistics
# file
mag_reps_gene_counts_and_tpm_long = 
  rbind(mag_reps_cazyme_gene_tpm_per_sample,mag_reps_peptidase_gene_tpm_per_sample,
      mag_reps_sulfatase_gene_tpm_per_sample,mag_reps_tbdt_gene_tpm_per_sample,
      mag_reps_gene_counts_long) %>%
  left_join(mag_summary_stats, by="MAG_name") %>%
  arrange(GTDB_Class,GTDB_Genus,MAG_name) %>%
  mutate(Sample = str_replace_all(Sample, "_TPM", ""))

# Set factor orders
mag_reps_gene_counts_and_tpm_long$GTDB_Class <- 
  factor(mag_reps_gene_counts_and_tpm_long$GTDB_Class,
         levels=unique(mag_reps_gene_counts_and_tpm_long$GTDB_Class))

mag_reps_gene_counts_and_tpm_long$MAG_name <- 
  factor(mag_reps_gene_counts_and_tpm_long$MAG_name,
         levels=unique(mag_reps_gene_counts_and_tpm_long$MAG_name))

mag_reps_gene_counts_and_tpm_long$Sample <- 
  factor(mag_reps_gene_counts_and_tpm_long$Sample,
         levels = c("Count_per_Mbp", "S10_SRF", "S25_SRF", "S25_BML", "S26_BML"))

func_gene_col <- c("Sulfatases" = "#ffe680ff", "Peptidases" = "#5a2ca0ff", 
                   "CAZymes" = "#8dd35fff", "TBDTs" = "#784421ff")

# Plot stacked bar chart
Figure_7 = ggplot(mag_reps_gene_counts_and_tpm_long, 
                  aes(x=Value, y=MAG_name)) + 
  geom_bar(aes(fill=Gene_group), colour = "black", stat="identity", 
           position="stack", show.legend = TRUE) + 
  scale_fill_manual(values = func_gene_col) + 
  theme_bw() + 
  facet_grid(GTDB_Class~Sample, scales = "free", space = "free_y") + 
  theme(panel.background = element_rect(fill = "white", color="black"), #
        panel.grid.major.y = element_line(colour = "gray", linetype = "dashed", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        strip.placement = "outside", 
        panel.spacing = unit(0, "lines"), 
        strip.background = element_rect(colour="black", fill="white"), 
        strip.text.x = element_text(size = 12, colour="black"),
        strip.text.y.right = element_text(size=10, colour="black", angle=0),
        legend.text = element_text(size=12),
        legend.title = element_blank(), 
        legend.position = "bottom", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size = 7, colour = "black")) + 
  guides(fill=guide_legend(size = 5, ncol=4))
Figure_7

### Export Figure 7
pdf(file="figures_output/Figure_7.pdf", height=8, width=12)
Figure_7
dev.off()

###################


#####
# Figure 8
#####

# Import new libraries
library(lemon)

# Import MAG cazyme gene TPM values
MAG_reps_cazy_genes_TPM_wide = 
  read.table("FRAM_WSC20_MAG_reps_cazyme_genes_with_TPM.txt", header = TRUE, 
             sep = "\t", as.is=TRUE)

# Import MAG RBP gene TPM values
MAG_reps_rbp_genes_TPM = 
  read.table("FRAM_WSC20_MAG_reps_RBP_genes_with_TPM.txt", header = TRUE, 
             sep = "\t", as.is=TRUE)

# Import MAG cazyme gene relative proportion of transcript values
MAG_reps_cazy_family_trans_counts_long = 
  read.table("FRAM_WSC20_MAG_reps_cazyme_per_family_transcript_counts.txt", header = TRUE, 
             sep = "\t", as.is=TRUE)

# Import taxonomic information for the MAGs
MAG_reps_summary_info = 
  read.table("FRAM_WSC20_MAG_reps_summary_statistics_and_taxonomy.csv", header = TRUE, 
             sep = ",", as.is=TRUE, check.names=F)

# Import CAZyme gene name to predicted substrate target information
cazyme_to_substrate = 
  read.table("cazyme_gene_to_target_substrate.csv", header = TRUE, 
             sep = ",", as.is=TRUE, check.names=F)

# Generate a refined MAG to taxonomy table
MAG_reps_taxa = MAG_reps_summary_info %>%
  subset(., select=c(MAG_name,GTDB_Class,GTDB_Family,GTDB_Genus,
                     SILVA_Genus))

# Determine the average RBP gene TPM for each MAG in each sample
MAG_reps_rbp_gene_avg_TPM = MAG_reps_rbp_genes_TPM %>%
  subset(., select=-c(Gene_name)) %>% 
  reshape2::melt(id.vars=c("MAG_name","Gene"), variable.name = "Sample", 
                 value.name="RBP_TPM") %>%
  aggregate(RBP_TPM~MAG_name+Sample, data=., FUN=mean)

### Filter dataset to plot only selected MAGs:

# Step 1: select only those MAGs with the highest total CAZyme TPM in each 
# sample
MAG_selected_list = mag_reps_cazyme_genes_with_tpm %>%
  subset(., select=-c(Gene_name)) %>%
  reshape2::melt(., id.vars=c("MAG_name","Gene"), variable.name="Sample", 
                 value.name="Value") %>%
  subset(., select=-c(Gene)) %>%
  group_by(MAG_name,Sample) %>%
  summarise(., across(everything(), sum)) %>%
  group_by(Sample) %>%
  arrange(desc(Value)) %>%
  slice_head(n=6) %>%
  as.data.frame(.) %>%
  subset(., select=c(MAG_name)) %>%
  unique()
dim(MAG_selected_list)


# Step 2: Select only those CAZyme gene families that are transcribed
# at equal to or higher levels than ribosomal protein genes, e.g. they are 
# more transcribed than the housekeeping/MAG average
MAG_reps_selected_cazy_families = MAG_reps_cazy_genes_TPM_wide %>%
  subset(., select=-c(Gene_name)) %>%
  filter(MAG_name %in% MAG_selected_list$MAG_name) %>%
  reshape2::melt(id.vars=c("MAG_name","Gene"), variable.name = "Sample", 
                 value.name="CAZyme_TPM") %>%
  left_join(MAG_reps_rbp_gene_avg_TPM, by=c("MAG_name","Sample")) %>%
  filter(CAZyme_TPM >= RBP_TPM & CAZyme_TPM != 0)

# Using the table generated above, now extract the MAGs and Genes from the 
# MAG_reps_cazy_family_trans_counts_long table so we can plot relative
# proportion of transcription for each gene family
MAG_reps_and_selected_cazys_long = MAG_reps_cazy_family_trans_counts_long %>%
  filter(MAG_name %in% MAG_reps_selected_cazy_families$MAG_name & 
           Gene %in% MAG_reps_selected_cazy_families$Gene) %>%
  left_join(MAG_reps_taxa, by="MAG_name") %>%
  arrange(GTDB_Class,GTDB_Family,GTDB_Genus) %>%
  mutate(Sample = fct_relevel(Sample, c("S10_SRF","S25_SRF","S25_BML","S26_BML"))) %>%
  mutate(MAG_name_new = case_when(
    SILVA_Genus == GTDB_Genus ~ paste0(MAG_name," (", GTDB_Genus,")"),
    GTDB_Genus == "Unassigned" & 
      SILVA_Genus != "NA" ~ paste0(MAG_name," (", SILVA_Genus,")"),
    GTDB_Genus == "Unassigned" ~ paste0(MAG_name," (",GTDB_Family,")"),
    GTDB_Genus != "Unassigned" & 
      SILVA_Genus == "NA" ~ paste0(MAG_name, " (", GTDB_Genus,")"),
    TRUE ~ paste0(MAG_name," (", GTDB_Genus," & ",SILVA_Genus,")"))) %>%
  mutate(MAG_name_new = str_replace(MAG_name_new, "FRAM_WSC20_", "")) %>%
  mutate(MAG_name_new = str_replace(MAG_name_new, " & NA",""))

# Add in predicted substrate information for the gene families
# The table imported above - cazyme_to_substrate - contains this information
# It was prepared by searching CAZyDB for each CAZyme
# gene family and conservatively designating predicted substrate targets
MAG_reps_and_selected_cazys_for_plotting = MAG_reps_and_selected_cazys_long %>%
  left_join(cazyme_to_substrate, by="Gene") %>%
  arrange(Target_broad,GTDB_Class,GTDB_Family,GTDB_Genus)

# Fix factor variables
MAG_reps_and_selected_cazys_for_plotting$MAG_name_new <- factor(
  MAG_reps_and_selected_cazys_for_plotting$MAG_name_new, levels=unique(
    MAG_reps_and_selected_cazys_for_plotting$MAG_name_new))

MAG_reps_and_selected_cazys_for_plotting$Gene <- factor(
  MAG_reps_and_selected_cazys_for_plotting$Gene, levels=unique(
    MAG_reps_and_selected_cazys_for_plotting$Gene))

# Determine number of colours needed
MAG_reps_and_selected_cazys_for_plotting$MAG_name_new %>%
  unique

# Plot Figure 8
library(RColorBrewer)
Figure_8 = ggplot(MAG_reps_and_selected_cazys_for_plotting, 
       aes(x=Rel_transcript_per_fam, y=Gene)) + 
  geom_bar(aes(fill=MAG_name_new), position="stack", stat="identity") +
  facet_grid(.~Sample, scales="fixed") + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(23)) + 
  labs(x = "Relative proportion of whole community gene family transcription (%)",
       fill = "MAG name") + 
  theme_bw() + 
  theme(strip.background.x = element_rect(fill = "white", colour = "black"),
        strip.text.x = element_text(colour = "black", size = 14),
        legend.text = element_text(colour = "black", size = 10),
        legend.title = element_text(colour = "black", size = 12),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.title.x = element_text(colour = "black", size = 12)) + 
  guides(fill=guide_legend(size = 5, ncol=1))

### Export Figure 8
pdf(file="figures_output/Figure_8.pdf", height=8, width=12)
Figure_8
dev.off()

