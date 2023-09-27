# libraries
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
setwd("~/Lugol_KMD")

load(file = "ps_rel.RData") #renamed, clr transformed data
nsamples(ps_rel)

mapfile = "LuKMD_metadata.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps_rel) <- map

ps_rel <- subset_samples(ps_rel, project == "Lugol")
ps_rel <- subset_samples(ps_rel, Species_coral == "APAL")
ps_rel <- subset_samples(ps_rel, Timepoint != "NA")

ps_rel = subset_samples(ps_rel, species_time != "APAL_Control_T4" | project != "Lugol")
ps_rel = subset_samples(ps_rel, species_time != "APAL_Lugol_T4")

pseq <- aggregate_rare(ps_rel, level = "Genus", detection = 2/100, prevalence = 10/100)
pseq@sam_data$project <- factor(pseq@sam_data$project, levels = c("Lugol", "KMD")) 

ps_rel_genus_melt <- pseq %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()
#head(ps_rel_genus_melt)

# summarize by Genus, sort, get top taxa
# Get families with mean relative abundance >0.01 across all samples 
genus <- ps_rel_genus_melt %>% group_by(Genus) %>% dplyr::summarise(Aver = mean(Abundance))
top <- genus[which(genus$Aver > 0.01),]
names <- top$Genus
names
ps_melt_top <- ps_rel_genus_melt[ps_rel_genus_melt$Genus %in% names,]


# plot
#apal KMD
pallette = c("#660066", "#277B00", "#FF8889", "#CC3300","#00ffff", "#6495ed", "#8CD09F", "#AC84DD","#BF5B9A", "#AAE247", "#003E94", "#FF9200", "#65BFFF", "#DADADA")
ps_melt_top$Genus <- factor(ps_melt_top$Genus, levels = c("[Caedibacter] taeniospiralis group", "Candidatus Amoebophilus", "Eilatimonas", "Family_BD2-7", "Family_Rhodobacteraceae", "MBIC10086" ,"P3OB-42", "Pelagibius", "Pleurocapsa PCC-7319", "Pseudomonas", "Synechococcus PCC-7336", "Thalassotalea", "Vibrio"))
nozero<- subset(ps_melt_top, Abundance > 0)
p <- ggplot(nozero, aes(x = trt, y = Genus)) +
  geom_point(aes(size=Abundance, fill = Genus), position = position_jitter(width = 0.25, height = 0), shape = 21) + #shape = 21 tells it to do colored circles
  facet_grid(Genotype~Timepoint, scales = "free_x", space = "free_x") + #splits into resistant v susceptible categories
  theme_facet() +
  ylab("Relative Abundance") +
  xlab("Treatment By Genotype") +
  labs(fill = "Genus", size = "Genus") +
  scale_y_discrete(limits=rev) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_fill_manual(values = pallette, 
                    breaks = c("[Caedibacter] taeniospiralis group", "Candidatus Amoebophilus", "Eilatimonas", "Family_BD2-7", "Family_Rhodobacteraceae", "MBIC10086" ,"P3OB-42", "Pelagibius", "Pleurocapsa PCC-7319", "Pseudomonas", "Synechococcus PCC-7336", "Thalassotalea", "Vibrio"),
                    labels = c("[Caedibacter] taeniospiralis group", "Candidatus Amoebophilus", "Eilatimonas", "Family_BD2-7", "Family_Rhodobacteraceae", "MBIC10086" ,"P3OB-42", "Pelagibius", "Pleurocapsa PCC-7319", "Pseudomonas", "Synechococcus PCC-7336", "Thalassotalea", "Vibrio")) +
                    guides(fill = guide_legend(override.aes = list(size=4))) + #makes legend circles bigger
  labs(fill = "Genus", size = "Relative \nAbundance")
p

ggsave(filename="rel/apal_bubble_KMD.png", plot=p, device="png",  width = 9, height = 6, dpi=500)

#apal Lugol
pallette = c("#006BFF", "#ffd700", "#277B00", "#000080", "#ff0000","#00ffff", "#6495ed" , '#843BFF', "#BF5B9A", "#FF9200")
ps_melt_top$Genus <- factor(ps_melt_top$Genus, levels = c("Acinetobacter", "Ascidiimonas", "Candidatus Amoebophilus", "Corynebacterium", "Family_Caulobacteraceae", "Family_Rhodobacteraceae", "MBIC10086" ,"Nannocystis",  "Pleurocapsa PCC-7319", "Thalassotalea"))
nozero<- subset(ps_melt_top, Abundance > 0)
p <- ggplot(nozero, aes(x = trt, y = Genus)) +
  geom_point(aes(size=Abundance, fill = Genus), position = position_jitter(width = 0.25, height = 0), shape = 21) + #shape = 21 tells it to do colored circles
  facet_grid(Genotype~Timepoint, scales = "free_x", space = "free_x") + #splits into resistant v susceptible categories
  theme_facet() +
  ylab("Relative Abundance") +
  xlab("Treatment By Genotype") +
  labs(fill = "Genus", size = "Genus") +
  scale_y_discrete(limits=rev) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_fill_manual(values = pallette, 
                    breaks = c("Acinetobacter", "Ascidiimonas", "Candidatus Amoebophilus", "Corynebacterium", "Family_Caulobacteraceae", "Family_Rhodobacteraceae", "MBIC10086" ,"Nannocystis",  "Pleurocapsa PCC-7319", "Thalassotalea"),
                    labels = c("Acinetobacter", "Ascidiimonas", "Candidatus Amoebophilus", "Corynebacterium", "Family_Caulobacteraceae", "Family_Rhodobacteraceae", "MBIC10086" ,"Nannocystis",  "Pleurocapsa PCC-7319", "Thalassotalea")) +
  guides(fill = guide_legend(override.aes = list(size=4))) + #makes legend circles bigger
  labs(fill = "Genus", size = "Relative \nAbundance")
p
ggsave(filename="rel/apal_bubble_Lugol.png", plot=p, device="png",  width = 7, height = 6, dpi=500)

#ofav Lugol
pallette = c("#667CFF", "#6BC80E", "#00ffff", "#D60000", "#0EC8C8", "#8400C2", "#FF8310", "#EE00F1", "#ccff33")
ps_melt_top$Genus <- factor(ps_melt_top$Genus, levels = c("Aestuariicella","Candidatus Amoebophilus", "Family_Rhodobacteraceae","Family_Terasakiellaceae", "Family_Unknown Family", "Fodinicurvata", "Marine Methylotrophic Group 3", "Pelagibius", "Woeseia"))
nozero<- subset(ps_melt_top, Abundance > 0)
p <- ggplot(nozero, aes(x = trt, y = Genus)) +
  geom_point(aes(size=Abundance, fill = Genus), position = position_jitter(width = 0.25, height = 0), shape = 21) + #shape = 21 tells it to do colored circles
  facet_grid(Genotype~Timepoint, scales = "free_x", space = "free_x") + #splits into resistant v susceptible categories
  theme_facet() +
  ylab("Relative Abundance") +
  xlab("Treatment By Genotype") +
  labs(fill = "Genus", size = "Genus") +
  scale_y_discrete(limits=rev) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_fill_manual(values = pallette, 
                    breaks = c("Aestuariicella","Candidatus Amoebophilus", "Family_Rhodobacteraceae","Family_Terasakiellaceae", "Family_Unknown Family", "Fodinicurvata", "Marine Methylotrophic Group 3", "Pelagibius", "Woeseia"),
                    labels = c("Aestuariicella","Candidatus Amoebophilus", "Family_Rhodobacteraceae","Family_Terasakiellaceae",  "Order Gammaproteobacteria Incertae Sedis", "Fodinicurvata", "Marine Methylotrophic Group 3", "Pelagibius", "Woeseia")) +
  guides(fill = guide_legend(override.aes = list(size=4))) + #makes legend circles bigger
  labs(fill = "Genus", size = "Relative \nAbundance")
p
ggsave(filename="rel/ofav_bubble_Lugol.png", plot=p, device="png",  width = 9, height = 6, dpi=500)

#ofav KMD
pallette = c("#6BC80E", "#00DABC","#FFDA75", "#D60000","#8400C2", "#FF8310", "#EE00F1", "#ccff33")
ps_melt_top$Genus <- factor(ps_melt_top$Genus, levels = c("Candidatus Amoebophilus","Eilatimonas","Ekhidna", "Family_Terasakiellaceae", "Fodinicurvata", "Marine Methylotrophic Group 3", "Pelagibius"))
nozero<- subset(ps_melt_top, Abundance > 0)
p <- ggplot(nozero, aes(x = trt, y = Genus)) +
  geom_point(aes(size=Abundance, fill = Genus), position = position_jitter(width = 0.25, height = 0), shape = 21) + #shape = 21 tells it to do colored circles
  facet_grid(Genotype~Timepoint, scales = "free_x", space = "free_x") + #splits into resistant v susceptible categories
  theme_facet() +
  ylab("Relative Abundance") +
  xlab("Treatment By Genotype") +
  labs(fill = "Genus", size = "Genus") +
  scale_y_discrete(limits=rev) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_fill_manual(values = pallette, 
                    breaks = c("Candidatus Amoebophilus","Eilatimonas","Ekhidna", "Family_Terasakiellaceae", "Fodinicurvata", "Marine Methylotrophic Group 3", "Pelagibius"),
                    labels = c("Candidatus Amoebophilus","Eilatimonas","Ekhidna", "Family_Terasakiellaceae", "Fodinicurvata", "Marine Methylotrophic Group 3", "Pelagibius")) +
  guides(fill = guide_legend(override.aes = list(size=4))) + #makes legend circles bigger
  labs(fill = "Genus", size = "Relative \nAbundance")
p
ggsave(filename="rel/ofav_bubble_KMD.png", plot=p, device="png",  width = 9, height = 6, dpi=500)

