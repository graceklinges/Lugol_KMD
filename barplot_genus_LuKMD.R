library("phyloseq")
library("data.table")
library("plyr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("reshape2")
library("indicspecies")
library("ggnetwork")
library("ape")
library("microbiome")
library("ggthemes")
library("cowplot")
library("ggsignif")

rm(list=ls())
setwd("~/Desktop/Lugol_KMD")

load(file = "ps_rel.RData") #renamed, clr transformed data
nsamples(ps_rel)

mapfile = "LuKMD_metadata.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps_rel) <- map

#ps_rel <- subset_samples(ps_rel, project == "KMD")
ps_rel <- subset_samples(ps_rel, Species_coral == "APAL")
ps_rel <- subset_samples(ps_rel, Timepoint != "NA")


pseq <- aggregate_rare(ps_rel, level = "Genus", detection = 2/100, prevalence = 10/100)
pseq@sam_data$project <- factor(pseq@sam_data$project, levels = c("Lugol", "KMD")) 
#pseq@sam_data$Exposure_basic <- factor(pseq@sam_data$Exposure_basic, levels = c("Baseline","Nutrients 3 weeks","No Disease Exposure", "Disease Exposed", "Disease Exposed - Unaffected", "No Disease Exposure - Mortality")) 

#pseq@sam_data$Replicate <- as.factor(pseq@sam_data$Replicate) 

ps_rel_genus_melt <- pseq %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()
head(ps_rel_genus_melt)

# Get families with mean relative abundance >0.01 across all samples 
gen_sum <- ps_rel_genus_melt %>% group_by(Genus) %>% dplyr::summarise(Aver = mean(Abundance))
gen_sub <- gen_sum[which(gen_sum$Aver > 0.01),]
names <- gen_sub$Genus
names

# Replace minor taxa with <0.01 abundance with NA
ps_rel_genus_melt$genus <- ps_rel_genus_melt$Genus

#look at your top taxa (names command above) and copy them below, replacing the taxa I have written

#All APAL
ps_rel_genus_melt$genus[ps_rel_genus_melt$genus != "[Caedibacter] taeniospiralis group" & 
                          ps_rel_genus_melt$genus != "Acinetobacter" &
                          ps_rel_genus_melt$genus != "Ascidiimonas" &
                          ps_rel_genus_melt$genus != "Candidatus Amoebophilus" &
                          ps_rel_genus_melt$genus != "Eilatimonas" &
                          ps_rel_genus_melt$genus != "Family_BD2-7" &
                          ps_rel_genus_melt$genus != "Family_Caulobacteraceae" &
                          ps_rel_genus_melt$genus != "Family_Rhodobacteraceae" &
                          ps_rel_genus_melt$genus != "MBIC10086" &
                          ps_rel_genus_melt$genus != "P3OB-42" &
                          ps_rel_genus_melt$genus != "Pleurocapsa PCC-7319" &
                          ps_rel_genus_melt$genus != "Thalassotalea"] <- NA

#APAL KMD
ps_rel_genus_melt$genus[ps_rel_genus_melt$genus != "[Caedibacter] taeniospiralis group" & 
                          ps_rel_genus_melt$genus != "Candidatus Amoebophilus" &
                          ps_rel_genus_melt$genus != "Eilatimonas" &
                          ps_rel_genus_melt$genus != "Family_BD2-7" &
                          ps_rel_genus_melt$genus != "Family_Rhodobacteraceae" &
                          ps_rel_genus_melt$genus != "MBIC10086" &
                          ps_rel_genus_melt$genus != "P3OB-42" &
                          ps_rel_genus_melt$genus != "Pelagibius" &
                          ps_rel_genus_melt$genus != "Pleurocapsa PCC-7319" &
                          ps_rel_genus_melt$genus != "Pseudomonas" &
                          ps_rel_genus_melt$genus != "Synechococcus PCC-7336" &
                          ps_rel_genus_melt$genus != "Thalassotalea" &
                          ps_rel_genus_melt$genus != "Vibrio"] <- NA

#OFAV KMD
ps_rel_genus_melt$genus[ps_rel_genus_melt$genus != "Candidatus Amoebophilus" & 
                          ps_rel_genus_melt$genus != "Eilatimonas" &
                          ps_rel_genus_melt$genus != "Ekhidna" &
                          ps_rel_genus_melt$genus != "Family_Terasakiellaceae" &
                          ps_rel_genus_melt$genus != "Fodinicurvata" &
                          ps_rel_genus_melt$genus != "Marine Methylotrophic Group 3" &
                          ps_rel_genus_melt$genus != "Pelagibius"] <- NA

#OFAV lugol
ps_rel_genus_melt$genus[ps_rel_genus_melt$genus != "Candidatus Amoebophilus" & 
                          ps_rel_genus_melt$genus != "Aestuariicella" &
                          ps_rel_genus_melt$genus != "Family_Rhodobacteraceae" &
                          ps_rel_genus_melt$genus != "Family_Terasakiellaceae" &
                          ps_rel_genus_melt$genus != "Family_Unknown Family" &
                          ps_rel_genus_melt$genus != "Fodinicurvata" &
                          ps_rel_genus_melt$genus != "Marine Methylotrophic Group 3" &
                          ps_rel_genus_melt$genus != "Woeseia" &
                          ps_rel_genus_melt$genus != "Pelagibius"] <- NA

#OFAV KMD and lugol
ps_rel_genus_melt$genus[ps_rel_genus_melt$genus != "Candidatus Amoebophilus" & 
                          ps_rel_genus_melt$genus != "Family_Terasakiellaceae" &
                          ps_rel_genus_melt$genus != "Family_Unknown Family" &
                          ps_rel_genus_melt$genus != "Marine Methylotrophic Group 3"] <- NA

#now go into ps_rel_genus_melt and check that only those taxa have names in the "genus" column and everything else says NA
names

# plot
#rename breaks and labels using the top taxa names above
#double check you have added a break for each taxa you kept from "names" or things will get weird
#you will need to pick a color for each taxon and add its hex code to "values"
#breaks and labels should be the same unless you want to rename any of your taxa in the label, i.e. if there are underscores. "labels" is what will show in the figure legend
#change x category to one of your variable names from your metadata file (I suggest narrow categories such as genotype, replicate, tank), change facet_grid to a different variable (I suggest broad categories that you only have a couple of e.g. treatment, species)

#apal KMD
bar_genus = ggplot(ps_rel_genus_melt, aes(x=Treatment, y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = reorder(genus, Abundance))) +
  scale_fill_manual(values=c("#660066", "#277B00", "#FF8889", "#CC3300","#00ffff", "#6495ed", "#8CD09F", "#AC84DD","#BF5B9A", "#AAE247", "#003E94", "#FF9200", "#65BFFF", "#DADADA"), 
                    breaks = c("[Caedibacter] taeniospiralis group", "Candidatus Amoebophilus", "Eilatimonas", "Family_BD2-7", "Family_Rhodobacteraceae", "MBIC10086" ,"P3OB-42", "Pelagibius", "Pleurocapsa PCC-7319", "Pseudomonas", "Synechococcus PCC-7336", "Thalassotalea", "Vibrio"),
                    labels = c("[Caedibacter] taeniospiralis group", "Candidatus Amoebophilus", "Eilatimonas", "Family_BD2-7", "Family_Rhodobacteraceae", "MBIC10086" ,"P3OB-42", "Pelagibius", "Pleurocapsa PCC-7319", "Pseudomonas", "Synechococcus PCC-7336", "Thalassotalea", "Vibrio"),
                    na.value = "#DADADA")  +
  facet_grid(Genotype~Timepoint, scales = "free_x", space = "free_x") +
  theme_bw() +
  #theme(legend.text = element_text(face = "italic")) +
  ylab("Relative Abundance") +
  xlab("Treatment By Genotype") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  labs(fill = "Genus", size = "Genus")
bar_genus
ggsave(filename="apal_trt_by_timepoint_KMD.png", plot=bar_genus, device="png",  width = 7, height = 4, dpi=500)

#apal KMD and Lugol
bar_genus = ggplot(ps_rel_genus_melt, aes(x=trt, y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = reorder(genus, Abundance))) +
  scale_fill_manual(values=c("#006BFF",  "#ffd700", "#660066", "#277B00", "#FF8889", "#CC3300", "#97DE00","#00ffff", "#6495ed", "#8CD09F","#BF5B9A", "#FF9200",  "#DADADA"), 
                    breaks = c("[Caedibacter] taeniospiralis group", "Acinetobacter", "Ascidiimonas", "Candidatus Amoebophilus", "Eilatimonas", "Family_BD2-7", "Family_Caulobacteraceae", "Family_Rhodobacteraceae", "MBIC10086" ,"P3OB-42","Pleurocapsa PCC-7319", "Thalassotalea"),
                    labels = c("Caedibacter taeniospiralis group", "Acinetobacter", "Ascidiimonas", "Candidatus Amoebophilus", "Eilatimonas", "Family BD2-7", "Family Caulobacteraceae", "Family Rhodobacteraceae", "MBIC10086" ,"P3OB-42","Pleurocapsa PCC-7319", "Thalassotalea"),
                    na.value = "#DADADA")  +
  facet_grid(project~Timepoint, scales = "free_x", space = "free_x") +
  theme_bw() +
  #theme(legend.text = element_text(face = "italic")) +
  ylab("Relative Abundance") +
  xlab("Treatment By Genotype") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  labs(fill = "Genus", size = "Genus")
bar_genus
ggsave(filename="apal_trt_by_timepoint_both.png", plot=bar_genus, device="png",  width = 7, height = 4, dpi=500)

#ofav kmd
bar_genus = ggplot(ps_rel_genus_melt, aes(x=Treatment, y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = reorder(genus, Abundance))) +
  scale_fill_manual(values=c("#6BC80E", "#00DABC","#FFDA75", "#D60000","#8400C2", "#FF8310", "#EE00F1", "#ccff33"), 
                    breaks = c("Candidatus Amoebophilus","Eilatimonas","Ekhidna", "Family_Terasakiellaceae", "Fodinicurvata", "Marine Methylotrophic Group 3", "Pelagibius"),
                    labels = c("Candidatus Amoebophilus","Eilatimonas","Ekhidna", "Family_Terasakiellaceae", "Fodinicurvata", "Marine Methylotrophic Group 3", "Pelagibius"),
                    na.value = "#DADADA")  +
  facet_grid(Genotype~Timepoint, scales = "free_x", space = "free_x") +
  theme_bw() +
  #theme(legend.text = element_text(face = "italic")) +
  ylab("Relative Abundance") +
  xlab("Treatment") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  labs(fill = "Genus", size = "Genus")
bar_genus
ggsave(filename="ofav_trt_by_timepoint_KMD.png", plot=bar_genus, device="png",  width = 7, height = 4, dpi=500)

#ofav lugol
bar_genus = ggplot(ps_rel_genus_melt, aes(x=Treatment, y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = reorder(genus, Abundance))) +
  scale_fill_manual(values=c("#667CFF", "#6BC80E", "#00ffff", "#D60000", "#0EC8C8", "#8400C2", "#FF8310", "#EE00F1", "#ccff33"), 
                    breaks = c("Aestuariicella","Candidatus Amoebophilus", "Family_Rhodobacteraceae","Family_Terasakiellaceae", "Family_Unknown Family", "Fodinicurvata", "Marine Methylotrophic Group 3", "Pelagibius", "Woeseia"),
                    labels = c("Aestuariicella","Candidatus Amoebophilus", "Family_Rhodobacteraceae","Family_Terasakiellaceae",  "Order Gammaproteobacteria Incertae Sedis", "Fodinicurvata", "Marine Methylotrophic Group 3", "Pelagibius", "Woeseia"),
                    na.value = "#DADADA")  +
  facet_grid(Genotype~Timepoint, scales = "free_x", space = "free_x") +
  theme_bw() +
  #theme(legend.text = element_text(face = "italic")) +
  ylab("Relative Abundance") +
  xlab("Treatment") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  labs(fill = "Genus", size = "Genus")
bar_genus
ggsave(filename="ofav_trt_by_timepoint_lugol.png", plot=bar_genus, device="png",  width = 7, height = 4, dpi=500)


#ofav all
bar_genus = ggplot(ps_rel_genus_melt, aes(x=trt, y=Abundance)) + 
  geom_bar(stat="identity", position="fill", aes(fill = reorder(genus, Abundance))) +
  scale_fill_manual(values=c("#6BC80E", "#D60000", "#0EC8C8", "#FF8310", "#EE00F1", "#ccff33"), 
                    breaks = c("Candidatus Amoebophilus","Family_Terasakiellaceae", "Family_Unknown Family", "Marine Methylotrophic Group 3"),
                    labels = c("Candidatus Amoebophilus","Family Terasakiellaceae", "Order Gammaproteobacteria Incertae Sedis", "Marine Methylotrophic Group 3"),
                    na.value = "#DADADA")  +
  facet_grid(project~Timepoint, scales = "free_x", space = "free_x") +
  theme_bw() +
  #theme(legend.text = element_text(face = "italic")) +
  ylab("Relative Abundance") +
  xlab("Treatment") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  labs(fill = "Genus", size = "Genus")
bar_genus
ggsave(filename="ofav_trt_by_timepoint_both.png", plot=bar_genus, device="png",  width = 8, height = 4, dpi=500)
