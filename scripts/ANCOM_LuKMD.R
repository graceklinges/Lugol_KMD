library(exactRankTests)
library(dplyr)
library(ggplot2)
library(compositions)
library(phyloseq)
library(robustbase)
library(ggpubr)
library("gridExtra")
library("cowplot")
rm(list=ls())

source("~/ancom_v2.1.R") #sourced from https://github.com/FrederickHuangLin/ANCOM/tree/master/scripts
setwd("~/Lugol_KMD")

# load unrarefied, renamed data for differential abundance analysis
load(file = "ps_pruned.RData")
sample_sums(ps)

mapfile = "LuKMD_metadata.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

#subset dataset (by whatever variables you are interested, could be treatment, site, coral species, etc.)
ps <- subset_samples(ps, Species_coral == "OFAV")
ps <- subset_samples(ps, project == "Lugol")
#ps <- subset_samples(ps, Timepoint == "T4" | Timepoint == "T0")

ps <- subset_samples(ps, Timepoint == "T4")
ps@sam_data$project <- factor(ps@sam_data$project, levels = c("Lugol", "KMD")) 

# filter to keep only most abundant taxa
ps = filter_taxa(ps, function(x) sum(x > 10) > (0.1*length(x)), TRUE)
# agglomerate to genus
ps <- tax_glom(ps, "Genus")
ntaxa(ps) #included in figure legend/plots

OTUdf <- as.data.frame(t(otu_table(ps)))
metadf <- read.delim(file = "LuKMD_metadata.txt")
colnames(metadf)[1] <- "Sample.ID"

#levels
metadf$T0_merge <- factor(metadf$T0_merge, levels = c("T0","Control_T4", "Lugol_T4")) 
metadf$T0_merge <- factor(metadf$T0_merge, levels = c("T0","Control_T4", "KoralMD_T4")) 

metadf$Treatment <- factor(metadf$Treatment, levels = c("Control", "Lugol")) 


# Step 1: Data preprocessing
feature_table = OTUdf; meta_data = metadf; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)

#zero_cut: Taxa with proportion of zeroes greater than zero_cut are not included in the analysis.
#out_cut: observations with proportion of mixture distribution less than out_cut (or greater than 1- out_cut) 
#will be detected as outliers

feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM
main_var = "Treatment"; p_adj_method = "fdr"; alpha = 0.05
adj_formula = NULL ; rand_formula = NULL #random formula can be any additional variable that might affect your data that you don't want to specifically look at, but rather ignore. e.g. replicate, tank effect

res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)

resdf <- as.data.frame(res$out)

# add taxonomy
tax<-as(tax_table(ps),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
tax<-as.data.frame(tax)
colnames(tax)
tax$taxa_id <- rownames(tax)
rownames(tax) <- NULL
dim(tax)
tax <- tax[,c(8,1:7)] #taxID
keep <- tax[tax$taxa_id %in% rownames(feature_table), ]
tax <- keep
rm(keep)

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.6"], label = "W[0.6]")

#test figure
fig = res$fig +  
   geom_hline(yintercept = cut_off["detected_0.6"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
              size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig  
 
# # specialized plot
figdf <- as.data.frame(fig$data)
figdf <- cbind(figdf,tax, by = "taxa.id")
figdf <- figdf[,-6]
figdf$group <- as.factor(figdf$group)
#can rename levels here or in inkscape
levels(figdf$group) <- c("Control T3 vs T0", "Lugol T3 vs T0")
levels(figdf$group) <- c("Control T4 vs T0", "KoralMD T4 vs T0")

# order genus
x = tapply(figdf$y, figdf$Genus, function(x) max(x))
x = sort(x, TRUE)
figdf$Genus = factor(as.character(figdf$Genus), levels=names(x))
figdf$col_genus <- figdf$Genus

#below here is relevant for plotting, I'm changing what taxa I am keeping in the final dataset based on the contrasts that were significant- look at resdf and figdf
#these are just examples, showing which taxa I kept for each contrast
#OFAV T0_merge KMD
figdf$col_genus[figdf$col_genus != "Pseudomonas" &
                  figdf$col_genus != "Vibrio" &
                  figdf$col_genus != "Thalassotalea"] <- NA
#APAL T0_merge KMD
figdf$col_genus[figdf$col_genus != "Family_Caulobacteraceae" &
                  figdf$col_genus != "Litoribrevibacter" &
                  figdf$col_genus != "Pseudomonas" &
                  figdf$col_genus != "Pseudoalteromonas" &
                  figdf$col_genus != "Thalassotalea" &
                  figdf$col_genus != "Enterobacter" &
                  figdf$col_genus != "Vibrio" &
                  figdf$col_genus != "Family_Enterobacteriaceae" &
                  figdf$col_genus != "Shimia" &
                  figdf$col_genus != "Acanthopleuribacter" &
                  figdf$col_genus != "Bermanella"] <- NA

#APAL t0 vs t3 Lugol
figdf$col_genus[figdf$col_genus != "Litoribrevibacter"] <- NA

#APAL Lugol T4
figdf$col_genus[figdf$col_genus != "Aliiglaciecola" &
                  figdf$col_genus != "Lentisphaera"] <- NA

#OFAV Lugol T4
figdf$col_genus[figdf$col_genus != "Allofrancisella" &
                figdf$col_genus != "MBIC10086" &
                  figdf$col_genus != "Family_Terasakiellaceae" &
                  figdf$col_genus != "Maritimimonas" &
                  figdf$col_genus != "Methylobacterium-Methylorubrum" &
                  figdf$col_genus != "Family_Enterobacteriaceae"] <- NA

#OFAV Lugol T0 vs T4
figdf$col_genus[figdf$col_genus != "Allofrancisella" &
                  figdf$col_genus != "Family_Enterobacteriaceae" &
                  figdf$col_genus != "Enterobacter" &
                  figdf$col_genus != "Hafnia-Obesumbacterium" &
                  figdf$col_genus != "Serratia" &
                  figdf$col_genus != "Family_Yersiniaceae" &
                  figdf$col_genus != "Family_Caulobacteraceae" &
                  figdf$col_genus != "Candidatus Amoebophilus" &
                  figdf$col_genus != "Family_Terasakiellaceae"] <- NA



#add new factor
figdf$col_genus <- factor(figdf$col_genus, levels = c(levels(figdf$col_genus), "Other"))
# convert NAs to other
figdf$col_genus[is.na(figdf$col_genus)] = "Other"

colors = c('#BE0032', '#F3C300', '#831DCB', '#37A337','#0067A5', '#848482')
apal_colors = c("#FF8700", "#848482","#BE0032","#831DCB", "#AA2384", "#F3C300", "#58B8EC","#E27AB7", "#AAE247", "#FF0087", '#848482')
ofav_colors = c("#D703E2", "#19A4EB", "#D60000", "#19EB91","#F56299", "#AA2384", "#2A6AFB", "#848482")
ofav_colors2 = c("#D703E2", "#19A4EB", "#D60000", "#19EB91","#F56299", "#AA2384", "#2A6AFB","#FB782A", "#FFCE00", "#848482")
# Taxa above the dashed line are significant with the null‐hypothesis rejected 80% of the time.
ggfig <- ggplot(figdf, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 2) +
  #facet_grid(~group) +
  theme_bw() +
  ylab("W statistic") +
  xlim(-3, 3) +
  theme(axis.title.x=element_blank()) +
  #theme(legend.position = "none") +
  scale_color_manual(name = "Genus", values = ofav_colors2) +
  geom_hline(yintercept = cut_off["detected_0.6"], linetype = "dashed") +
  annotate("text", label = "W = 0.6", size = 3.5, x = -1.3, y = 150, color = "red")
ggfig
ggsave(filename="ofav_lugol.svg", plot=ggfig, device="svg",  width = 5.75, height = 4, dpi=500)

ggfig <- ggplot(figdf, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 2) +
  facet_grid(~group) +
  theme_bw() +
  ylab("W statistic") +
  xlim(-3, 3) +
  theme(axis.title.x=element_blank()) +
  #theme(legend.position = "none") +
  scale_color_manual(name = "Genus", values = apal_colors) +
  geom_hline(yintercept = cut_off["detected_0.6"], linetype = "dashed") +
  annotate("text", label = "W = 0.6", size = 3.5, x = -1.4, y = 47, color = "red")
ggfig
ggsave(filename="ancom/apal_lugolT3.svg", plot=ggfig, device="svg",  width = 5, height = 4, dpi=500)

ggfig <- ggplot(figdf, aes(x = x, y = y, color = col_genus)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(size = 2) +
  facet_grid(~group) +
  theme_bw() +
  ylab("W statistic") +
  xlim(-3, 3) +
  theme(axis.title.x=element_blank()) +
  #theme(legend.position = "none") +
  scale_color_manual(name = "Genus", values = apal_colors) +
  geom_hline(yintercept = cut_off["detected_0.6"], linetype = "dashed") +
  annotate("text", label = "W = 0.6", size = 3.5, x = -1.4, y = 47, color = "red")
ggfig

