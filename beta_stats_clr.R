library('phyloseq'); packageVersion('phyloseq')
library('vegan'); packageVersion('vegan')
library("usedist")
library("ggplot2")
library("microbiome")
library(multcompView)

rm(list=ls())
setwd("~/Lugol_KMD")

pairwise.adonis.dm <- function(x,factors,stratum=NULL,p.adjust.m="fdr",perm=999){
  
  library(vegan)
  if(class(x) != "dist") stop("x must be a dissimilarity matrix (dist object)")
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  
  for(elem in 1:ncol(co)){
    sub_inds <- factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))
    resp <- as.matrix(x)[sub_inds,sub_inds]
    ad = adonis(as.dist(resp) ~
                  
                  factors[sub_inds], strata=stratum[sub_inds], permutations=perm);
    
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
    
  }
  
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


load(file = "ps_clr.RData") #renamed, clr transformed data
nsamples(clr)

mapfile = "LuKMD_metadata.txt"
map = import_qiime_sample_data(mapfile)
sample_data(clr) <- map

#clr <- subset_samples(clr, Species_coral == "OFAV")
clr <- subset_samples(clr, project == "KMD")
#clr <- subset_samples(clr, Timepoint != "T4")
#clr <- subset_samples(clr, Timepoint == "T4" | Timepoint == "T0")

clr <- subset_samples(clr, sample.id != "OF3-41-T3-L")

nsamples(clr)

clr@sam_data$T0_merge<- factor(clr@sam_data$T0_merge, levels = c("T0", "Control_T1", "Control_T2", "Control_T3", "Control_T4", "Lugol_T1", "Lugol_T2", "Lugol_T3", "Lugol_T4")) 
clr@sam_data$Description<- factor(clr@sam_data$Description, levels = c("Pre-Treatment","1 Day Post-Treatment","After 6 treatments ", "1 Month Washout", "2 Months Washout"))


### Calculate distance matrices
# only calculating euclidean distance for clr transformed data because we have negative values, and with the clr transform we've moved our data into "real space"
clr_euc <- phyloseq::distance(clr, method = "euclidean")

### PERMANOVAs

# PERMANOVA's with Adonis -  Permutational Multivariate Analysis of Variance
# Make a data frame from the sample_data
sampledf <- data.frame(sample_data(clr))
# Adonis 
adonis(dist ~ variable, data = metadata) #significantly dif by project
adonis(clr_euc ~ Timepoint*Treatment, data = sampledf) #interaction not sig, only timepoint sig

species_time_KMD <- pairwise.adonis.dm(clr_euc, sample_data(clr)$species_time, p.adjust.m = "fdr")
write.table(species_time_KMD, file = "beta_stats_species_time_KMD.txt", sep = "\t")
#major takeaways: KMD: both apal and ofav not sig different by treatment at T4. Ofav did not change in community composition over time, Apal did
#Lugol: communities at T3 were not significantly different by treatment at T3. Ofav did not change in community composition over time, Apal did
pairwise.adonis.dm(clr_euc, sample_data(clr)$Treatment, p.adjust.m = "fdr")
#sig for ofav, NS for apal


### PERMDISPR
anova(betadisper(clr_euc, sampledf$Treatment, bias.adjust = TRUE)) #NS
anova(betadisper(clr_euc, sampledf$Timepoint, bias.adjust = TRUE)) #NS for KMD

p.adjust(permutest(betadisper(clr_euc, sampledf$T0_merge, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr') #all NS KMD
#explore further


#PCA ordination
ord_clr <- phyloseq::ordinate(clr, "RDA") #RDA without constraints = PCA
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n") #screen plot plateaus quickly
head(ord_clr$CA$eig)
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))    

#clr PCA figure
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
PCA <- phyloseq::plot_ordination(clr, ord_clr, type="samples", color="Description", shape="Treatment") + 
  geom_point(size = 2) +
  coord_fixed(((clr2 / clr1))*2) +
  stat_ellipse(aes(group = Timepoint, linetype = project)) +
  labs(color = "Timepoint") +
  labs(shape = "Treatment") +
  scale_color_manual(values = c("#ffa500", "#00ff7f", "#00bfff", "#0000ff", "#ff1493", "#000000", "#FF0000"))
PCA
ggsave(filename="KMD_PCA-OFAV.png", plot=PCA, width = 6, height = 6, device="png", dpi=700)

clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
PCA <- phyloseq::plot_ordination(clr, ord_clr, type="samples", color="T0_merge") + 
  geom_point(size = 2) +
  coord_fixed(((clr2 / clr1))*2) +
  stat_ellipse(aes(group = Species_coral, linetype = Species_coral)) +
  labs(color = "Timepoint") +
  labs(shape = "Treatment") +
  scale_color_manual(values = c("#3786E2", "#EF0000", "#8F8C8C", "#0000ff", "#ff1493", "#000000", "#FF0000"))
PCA
ggsave(filename="Lugol_both_species.png", plot=PCA, width = 6, height = 6, device="png", dpi=700)


#dispersion
adonis(clr_euc ~ T0_merge, data = sampledf) 
dispr <- betadisper(clr_euc, sampledf$T0_merge, bias.adjust = TRUE)
dispr
#Average distance to median:
#  Ammonium  Combined      Ctrl   Nitrate Phosphate        T0 
#21.74     29.20     29.45     33.83     29.97     31.58 

#Average distance to median:
#  A3    A6    C3    C6 CTRL3 CTRL6    N3    N6    P3    P6    T0 
#24.64 18.15 37.58 21.38 27.85 30.24 36.82 32.01 35.49 24.87 31.58 

plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "") #quick dispersion plot
permutest(dispr)

# extract distance to centroid
clr <- subset_samples(clr, Species_coral == "APAL")
clr_euc <- phyloseq::distance(clr, method = "euclidean")
sampledf <- data.frame(sample_data(clr))

disp <- betadisper(clr_euc, sampledf$species_time_wT0, bias.adjust = TRUE)
dispd <- as.data.frame(disp$distances)
dispd <- cbind(dispd, sample_data(clr))
colnames(dispd)[1] <- "distance"

dispersion_stats <- p.adjust(permutest(betadisper(clr_euc, sampledf$species_time_wT0, 
                                                  bias.adjust = TRUE), 
                                       pairwise=TRUE)$pairwise$permuted, method = 'fdr')

write.table(dispersion_stats, file = "beta/dispersion_stats_KMD_APAL.txt", sep = "\t")
write.table(dispersion_stats, file = "beta/dispersion_stats_KMD_OFAV.txt", sep = "\t")
write.table(dispersion_stats, file = "beta/dispersion_stats_Lugol_APAL.txt", sep = "\t")
write.table(dispersion_stats, file = "beta/dispersion_stats_Lugol_OFAV.txt", sep = "\t")


#manually turned this into a triangular matrix
dispersion_table <- read.table("beta/dispersion_lugol_OFAV_matrix.txt", header = TRUE, sep = "\t", row.names = 1)
myletters<-multcompLetters(dispersion_table,compare="<=",threshold=0.05,Letters=letters)

myCols <- c("#EF476F", "#06D6A0", "#FFC107", "#843BFF", "#3BFF51")

dispersion_plot <- ggplot(dispd, aes(x=Timepoint, y=distance)) +
  facet_grid(Species_coral~Treatment)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = Description), 
             position = position_jitter(width = .25, height = 0)) +
  ylab("Distance-to-centroid") +
  theme_bw() +
  #theme(axis.title.x = element_blank(),
    #    legend.position = "none") +
  labs(color = "Description") +
  scale_colour_manual(values = myCols) +
  stat_summary(geom = 'text', label = c("a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a"), fun.y = max, vjust = -1, size = 3.5) +
  #ggtitle("Dispersion Over Time: OFAV") +
  #scale_y_continuous(expand = expand_scale(mult = c(.1))) +
  ylim(0,300)
dispersion_plot
ggsave(filename="beta/disp_KMD.svg", plot=dispersion_plot, height = 4.5, width = 8, device="svg", dpi=500)

