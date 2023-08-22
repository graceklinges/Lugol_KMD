library(dada2); packageVersion("dada2")
library(phyloseq)

setwd("~/Desktop/KoralMD/")

# Location of raw reads
#If they are all in the same folder set the same path for both
path <- "koralmd"

# Double check you've set the right paths - should list your read 1 files
list.files(path)

# Assumes forward and reverse fastq filenames have format: ***_R1.fastq.gz and ***_R2.fastq.gz, change to match your files
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

# Plot read quality- you expect read quality to drop off towards the end
plotQualityProfile(fnFs[1:4])
#plotQualityProfile(fnRs[1:4])
#
#reads look super good
#In gray-scale is a heat map of the frequency of each quality score at each 
#base position. The mean quality score at each position is shown by the green 
#line, and the quartiles of the quality score distribution by the orange lines.
#The red line shows the scaled proportion of reads that extend to at least that
#position (this is more useful for other sequencing technologies, as Illumina 
#reads are typically all the same length, hence the flat red line).
#
length(fnFs) 
length(fnRs) #both 200, correct

#this is just a good idea to record
write.table(sample.names, file = "names.txt", sep="\t")

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim, change parameters how you see fit. truncLen usually based off of plotQualityProfile
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
write.table(out, file = "filter_out.txt", sep="\t")

#keep <- out[,"reads.out"] > 100 # Or other cutoff. 
#minimum read depth was 170,759 lol
#filtFs <- filtFs[keep]
#filtRs <- filtRs[keep]
#filtFs <- filtFs[!is.na(filtFs)] #this removes any F reads when the corresponding R read had less than 100 reads
#filtRs <- filtRs[!is.na(filtRs)]

sample.names <- names(filtFs)

set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# If you like, visualize these error rates. You want the estimated error rates (black lines) to be a good 
# fit to the observed error rates (points) and for the error rates drop with increased quality
# plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

# Note: if this loop doesn't work, you've done something wrong. You've probably made an error earlier in the code, 
# probably with the location of your forward and reverse reads and it's generated duplicate sample names
# Normal to get warning message that there are duplicate sequences
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE) # core sample interference algorithm
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE) # core sample interference algorithm
  merger <- mergePairs(ddF, derepF, ddR, derepR) # Merge paired end reads
  mergers[[sam]] <- merger
}


# Make sequence table from merged reads
st.all <- makeSequenceTable(mergers) # Normal to get warning message saying the sequences being tabled vary in length

ncol(st.all) #576923
nrow(st.all) #200

names_st.all <- row.names(st.all)
write.table(names_st.all, file = "names_merged.txt", sep="\t")

table(nchar(getSequences(st.all)))

# Remove any ASVs that are considerably off target length, depends on your primers
seqtab_trimmed <- st.all[,nchar(colnames(st.all)) %in% seq(290,294)]

# Inspect distribution of read lengths after removal of off-target reads
table(nchar(getSequences(seqtab_trimmed)))

# Remove chimeric sequences
seqtab <- removeBimeraDenovo(seqtab_trimmed, method="consensus", multithread=TRUE, verbose = T) #Identified 194157 bimeras out of 223860 input sequences.
sum(st.all)-sum(seqtab) # How many chimeras were removed? 32840844
sum(seqtab)/sum(st.all) #percent chimeras 0.6114076
ncol(seqtab) #number of ASVs left 90155

names_seqtab <- row.names(seqtab)
identical(names_st.all,names_seqtab) #were there samples lost with removing off target reads or chimeras
write.table(names_seqtab, file = "names_final.txt", sep="\t")

getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=names_st.all, dada2_input=out[,1],
                          filtered=out[,2],
                          nonchim=rowSums(seqtab),
                          final_perc_reads_retained=round(rowSums(seqtab)/out[,1]*100, 1))

summary_tab #sometimes this doesn't work, don't worry about it too much, just tells you how many reads lost with filtering
write.table(summary_tab, file = "reads_lost.txt", sep="\t")

# Save chimera-free ASV table as downstream tasks may cause R to crash
saveRDS(seqtab, "seqtab-KMD.rds")

readRDS("seqtab.rds")
# Assign taxonomy based on silva reference database at species (100%) level, you must have the appropriate Silva database downloaded
tax_silva <- assignTaxonomy(seqtab, "/Volumes/Grace\ External\ 2/Databases/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)

library(phyloseq)
# Export sequence table with genus and species assignments as phyloseq objects
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(tax_silva))

# Save as RDS objects
saveRDS(ps, file = "ps_object.rds")
save(ps, file = "ps_full.RData")
