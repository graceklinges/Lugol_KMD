library(phyloseq)
library(multcomp)
library(picante)
library(ggplot2 )
library(tidyverse)
library(plyr)
library(dplyr)
library(ggplot2)

setwd("~/Lugol_KMD")
rm(list=ls())

load(file = "ps_unpruned.RData") #renamed, unpruned ps object

mapfile = "LuKMD_metadata.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

ps = subset_samples(ps, project == "KMD")
ps = subset_samples(ps, species_time != "APAL_Control_T4")
ps = subset_samples(ps, species_time != "APAL_Lugol_T4")#for lugol only, removing APAL T4
#ps = subset_samples(ps, Species_coral == "APAL")
nsamples(ps)
#functions
asinTransform <- function(p) { asin(sqrt(p)) }

#mymat
library(multcompView)

tri.to.squ<-function(x)
{
  rn<-row.names(x)
  cn<-colnames(x)
  an<-unique(c(cn,rn))
  myval<-x[!is.na(x)]
  mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
  for(ext in 1:length(cn))
  {
    for(int in 1:length(rn))
    {
      if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
      mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
    }
    
  }
  return(mymat)
}


# Estimate faith's phylogenetic diversity 
estimate_pd <- function(phylo){
  
  # Error if input is not of class phylo
  if(class(phylo) != "phyloseq"){
    stop("Input file is not of class 'phyloseq'.")
  }
  
  # Error if no class phy_tree
  if(!(.hasSlot(phylo, "phy_tree"))){
    stop("Could not find tree slot in phylo object.")
  }
  
  # Transpose if needed
  # Adapted from phyloseq/vegan import
  OTU <- phyloseq::otu_table(phylo)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  
  # Get matrix version of OTU table
  otutable <- as(OTU, "matrix")
  
  # Get phylogenetic tree from pyloseq object
  tree <- phyloseq::phy_tree(phylo)
  
  # Print status message
  message("Calculating Faiths PD-index...")
  
  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }
  
  # Calculate Faith's PD-index
  pdtable <- picante::pd(otutable, tree, include.root = F)
  
  # Return data frame of results
  return(pdtable)
}

# Calculate standard error
sderr <- function(x) {sd(x)/sqrt(length(x))}
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sderr(x[[col]]), na.rm=TRUE)
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

normality.plots <- function(x) {
  par(mfrow=c(2,2))
  hist(residuals(x), main = "Histogram", xlab = "Values")
  boxplot(residuals(x), main = "Boxplot", ylab = "Values")
  qqnorm(residuals(x))
  qqline(residuals(x))
  plot(density(residuals(x)), main = "Kernel Density Estimate")
}

#Initialize matrices to store alpha diversity estimates
nsamp = nsamples(ps)
observed <- matrix(nrow = nsamp)
row.names(observed) <- sample_names(ps)
simpson <- matrix(nrow =nsamp)
row.names(simpson) <- sample_names(ps)
shannon <- matrix(nrow =nsamp)
row.names(shannon) <- sample_names(ps)

#Calculate statistics
# Options for measures = ("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
# Calculate observed
obs <- as.numeric(as.matrix(estimate_richness(ps, measures = "Observed")))
observed[ ,] <- obs
colnames(observed) [1] <- "observed"
# Calculate simpson
simp <- as.numeric(as.matrix(estimate_richness(ps, measures = "Simpson")))
simpson[ ,] <- simp
colnames(simpson) [1] <- "simpson"
# Calculate shannon
shan <- as.numeric(as.matrix(estimate_richness(ps, measures = "Shannon")))
shannon[ ,] <- shan
colnames(shannon) [1] <- "shannon"

#Combine our estimates for observed, simpson, and Shannon indices into one dataframe
alpha <- cbind(observed,simpson,shannon)
head(alpha)
# Add the sample metadata into this dataframe
s <- data.frame(sample_data(ps))

alphadiv <- cbind(alpha, s)
head(alphadiv)
# alphadiv <- alphadiv[,-5]
# head(alphadiv)
levels(alphadiv$Description) <- c("Pre-Treatment","1 Day Post-Treatment", "After 6 treatments", "1 Month Washout", "2 Months Washout")

odiv <- subset(alphadiv, Species_coral == "OFAV")
adiv <- subset(alphadiv, Species_coral == "APAL")

#plots --> decided to just do stats within each species
shannon_trt_time <- pairwise.wilcox.test(odiv$shannon, odiv$species_time_wT0, p.adjust.method = 'fdr')
mymat<-tri.to.squ(shannon_trt_time[["p.value"]])
write.table(mymat, file = "alpha/ofav_shan_stats_KMD.txt", sep = "\t")
ofav_letters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)
ofav_letters$Letters

shannon_trt_time <- pairwise.wilcox.test(adiv$shannon, adiv$species_time_wT0, p.adjust.method = 'fdr')
mymat<-tri.to.squ(shannon_trt_time[["p.value"]])
write.table(mymat, file = "alpha/apal_shan_stats_KMD.txt", sep = "\t")
apal_letters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)
apal_letters$Letters

myCols <- c("#EF476F", "#06D6A0", "#FFC107", "#843BFF", "#3BFF51")
#goes with stats from KW_simpsons_No_level_weeks_noT0, use alpha_div_sig_letters to get letters
A <- ggplot(alphadiv, aes(x=Timepoint, y=shannon)) +
  facet_grid(Species_coral~Treatment)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = Description), 
             position = position_jitter(width = .25, height = 0)) +
  ylab("Shannon Diversity Index") +
  theme_bw() +
  #theme(aspect.ratio = 1.8) +
  theme(
    #legend.position = c(2,0.3),
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 12)) +
  scale_colour_manual(values = myCols) +
  stat_summary(geom = 'text', label = c("a","a", "a", "a", "a", "a", "a", "a", "a", "a", "a" , "a", "a", "a", "a", "a", "a", "a", "a", "a"), fun = max, vjust = -1, size = 3.5) + #stats with rick
    ylim(3,8) +
  #ggtitle("Shannon Diversity by Treatment Weeks") +
theme(plot.title = element_text(hjust = 0.5))
A
#merged with lugol plot in Illustrator
ggsave(filename="alpha/alphadiv_trtweeks_KMD.svg", plot=A, height = 4, width = 8, device="svg", dpi=500)

#load 
A <- ggplot(alphadiv, aes(x=Timepoint, y=shannon)) +
  facet_grid(Species_coral~Treatment)+
  geom_boxplot(outlier.shape = NA, color = "gray35") +
  geom_point(aes(color = Description), 
             position = position_jitter(width = .25, height = 0)) +
  ylab("Shannon Diversity Index") +
  theme_bw() +
  #theme(aspect.ratio = 1.8) +
  theme(
    #legend.position = c(2,0.3),
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 12)) +
  scale_colour_manual(values = myCols) +
  stat_summary(geom = 'text', label = c("a","a", "a", "a", "a", "a", "a", "a", "ab" , "a", "bc", "abc", "ab", "abc", "abc", "bc", "abc", "c"), fun = max, vjust = -1, size = 3.5) + #stats with rick
  ylim(3,8) +
  #ggtitle("Shannon Diversity by Treatment Weeks") +
  theme(plot.title = element_text(hjust = 0.5))
A
ggsave(filename="alpha/alphadiv_trtweeks_Lugol.svg", plot=A, height = 4, width = 8, device="svg", dpi=500)


