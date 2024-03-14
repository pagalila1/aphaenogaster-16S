# ============================================================
'R code for Microbiome analysis using DADA2

Christina Pavloudi
christina.pavloudi@embrc.eu
https://cpavloud.github.io/mysite/

	Copyright (C) 2024 Christina Pavloudi
  
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
  
    This script is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.'

# =============================================================


#The ShortRead and ggplot2 packages are available from Bioconductor
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")

setwd("/Volumes/Drive/Data/")

#Define the following path variable so that it points to the extracted directory on your machine:
path <- getwd() #CHANGE ME to the directory containing the fastq files after unzipping. (path = getwd)
path

#List the files in your current/working directory
fns <- list.files(path)
fns

# Forward and reverse fastq filenames
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names.Rs <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)

#Visualize the quality profile of the forward reads:
plot1 <- plotQualityProfile(fnFs[[1]])
plot2 <- plotQualityProfile(fnFs[[2]])
#Visualize the quality profile of the reverse reads:
plot3 <- plotQualityProfile(fnRs[[1]])
plot4 <- plotQualityProfile(fnRs[[2]])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(5,8), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
head(out)
#The standard filtering parameters are starting points, not set in stone. 
#If you want to speed up downstream computation, consider tightening maxEE. 
#If too few reads are passing the filter, consider relaxing maxEE, 
#perhaps especially on the reverse reads (eg. maxEE=c(2,5)), and reducing 
#the truncLen to remove low quality tails. #The maxEE parameter sets the maximum number 
#of expected errors allowed in a read, 
#Remember though, when choosing #truncLen for paired-end reads you must maintain overlap 
#after truncation in order to merge them later.

# Save this output as RDS file:
saveRDS(out, "filter_and_trim_out.rds")

############# paused here on Friday ################

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
plot_errF <- plotErrors(errF, nominalQ=TRUE)
plot_errR <- plotErrors(errR, nominalQ=TRUE)


# save error calculation as RDS files:
saveRDS(errF, "errF.rds")
saveRDS(errR, "errR.rds")

#readRDS()

#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
#derepRs <- derepFastq(filtRs, verbose=TRUE)

derepR1s <- derepFastq("/Volumes/Drive/Data/filtered/16S-PCR-neg_R_filt.fastq", verbose=TRUE)
derepR2s <- derepFastq("/Volumes/Drive/Data/filtered/A48-F1-G_R_filt.fastq", verbose=TRUE)
derepR3s <- derepFastq("/Volumes/Drive/Data/filtered/A48-F10-G_R_filt.fastq", verbose=TRUE)
derepR4s <- derepFastq("/Volumes/Drive/Data/filtered/A48-F2-G_R_filt.fastq", verbose=TRUE)
derepR5s <- derepFastq("/Volumes/Drive/Data/filtered/A48-F3-G-a_R_filt.fastq", verbose=TRUE)
derepR6s <- derepFastq("/Volumes/Drive/Data/filtered/A48-F3-G-b_R_filt.fastq", verbose=TRUE)
derepR7s <- derepFastq("/Volumes/Drive/Data/filtered/A48-F4-G-a_R_filt.fastq", verbose=TRUE)
derepR8s <- derepFastq("/Volumes/Drive/Data/filtered/A48-F4-G-b_R_filt.fastq", verbose=TRUE) 
derepR9s <- derepFastq("/Volumes/Drive/Data/filtered/A48-F5-G_R_filt.fastq", verbose=TRUE)
derepR10s <- derepFastq("/Volumes/Drive/Data/filtered/A48-F6-G_R_filt.fastq", verbose=TRUE)   
derepR11s <- derepFastq("/Volumes/Drive/Data/filtered/A48-F7-G_R_filt.fastq", verbose=TRUE)
derepR12s <- derepFastq("/Volumes/Drive/Data/filtered/A48-F8-G_R_filt.fastq", verbose=TRUE)
derepR13s <- derepFastq("/Volumes/Drive/Data/filtered/A48-F9-G_R_filt.fastq", verbose=TRUE)    
derepR14s <- derepFastq("/Volumes/Drive/Data/filtered/A48-N10-G_R_filt.fastq", verbose=TRUE)  
derepR15s <- derepFastq("/Volumes/Drive/Data/filtered/A48-N3-G_R_filt.fastq", verbose=TRUE)
derepR16s <- derepFastq("/Volumes/Drive/Data/filtered/A48-N4-G_R_filt.fastq", verbose=TRUE)
derepR17s <- derepFastq("/Volumes/Drive/Data/filtered/A48-N5-G_R_filt.fastq", verbose=TRUE)
derepR18s <- derepFastq("/Volumes/Drive/Data/filtered/A48-N6-G_R_filt.fastq", verbose=TRUE)   
derepR19s <- derepFastq("/Volumes/Drive/Data/filtered/A48-N7-G_R_filt.fastq", verbose=TRUE)
derepR20s <- derepFastq("/Volumes/Drive/Data/filtered/A48-N8-G_R_filt.fastq", verbose=TRUE)
derepR21s <- derepFastq("/Volumes/Drive/Data/filtered/A48-N9-G_R_filt.fastq", verbose=TRUE)
derepR22s <- derepFastq("/Volumes/Drive/Data/filtered/A51-F10-G_R_filt.fastq", verbose=TRUE) 
derepR23s <- derepFastq("/Volumes/Drive/Data/filtered/A51-F2-G_R_filt.fastq", verbose=TRUE)
derepR24s <- derepFastq("/Volumes/Drive/Data/filtered/A51-F3-G_R_filt.fastq", verbose=TRUE)   
derepR25s <- derepFastq("/Volumes/Drive/Data/filtered/A51-F8-G_R_filt.fastq", verbose=TRUE)
derepR26s <- derepFastq("/Volumes/Drive/Data/filtered/A51-F9-G_R_filt.fastq", verbose=TRUE)
derepR27s <- derepFastq("/Volumes/Drive/Data/filtered/A51-N1-G_R_filt.fastq", verbose=TRUE)
derepR28s <- derepFastq("/Volumes/Drive/Data/filtered/A51-N2-G_R_filt.fastq", verbose=TRUE)
derepR29s <- derepFastq("/Volumes/Drive/Data/filtered/A51-N3-G_R_filt.fastq", verbose=TRUE)
derepR30s <- derepFastq("/Volumes/Drive/Data/filtered/A51-N5-G_R_filt.fastq", verbose=TRUE)  
derepR31s <- derepFastq("/Volumes/Drive/Data/filtered/A51-N6-G_R_filt.fastq", verbose=TRUE)
derepR32s <- derepFastq("/Volumes/Drive/Data/filtered/A53-F1-G_R_filt.fastq", verbose=TRUE)
derepR33s <- derepFastq("/Volumes/Drive/Data/filtered/A53-F10-G_R_filt.fastq", verbose=TRUE)
derepR34s <- derepFastq("/Volumes/Drive/Data/filtered/A53-F2-G_R_filt.fastq", verbose=TRUE)
derepR35s <- derepFastq("/Volumes/Drive/Data/filtered/A53-F3-G_R_filt.fastq", verbose=TRUE)
derepR36s <- derepFastq("/Volumes/Drive/Data/filtered/A53-F4-G_R_filt.fastq", verbose=TRUE)
derepR37s <- derepFastq("/Volumes/Drive/Data/filtered/A53-F5-G_R_filt.fastq", verbose=TRUE)
derepR38s <- derepFastq("/Volumes/Drive/Data/filtered/A53-F6-G_R_filt.fastq", verbose=TRUE)
derepR39s <- derepFastq("/Volumes/Drive/Data/filtered/A53-F7-G_R_filt.fastq", verbose=TRUE)
derepR40s <- derepFastq("/Volumes/Drive/Data/filtered/A53-F8-G_R_filt.fastq", verbose=TRUE)
derepR41s <- derepFastq("/Volumes/Drive/Data/filtered/A53-F9-G_R_filt.fastq", verbose=TRUE)
derepR42s <- derepFastq("/Volumes/Drive/Data/filtered/A53-N1-G_R_filt.fastq", verbose=TRUE)
derepR43s <- derepFastq("/Volumes/Drive/Data/filtered/A53-N10-G_R_filt.fastq", verbose=TRUE)
derepR44s <- derepFastq("/Volumes/Drive/Data/filtered/A53-N2-G_R_filt.fastq", verbose=TRUE)
derepR45s <- derepFastq("/Volumes/Drive/Data/filtered/A53-N3-G_R_filt.fastq", verbose=TRUE)
derepR46s <- derepFastq("/Volumes/Drive/Data/filtered/A53-N4-G_R_filt.fastq", verbose=TRUE)
derepR47s <- derepFastq("/Volumes/Drive/Data/filtered/A53-N5-G_R_filt.fastq", verbose=TRUE)
derepR48s <- derepFastq("/Volumes/Drive/Data/filtered/A53-N6-G_R_filt.fastq", verbose=TRUE)
derepR49s <- derepFastq("/Volumes/Drive/Data/filtered/A53-N7-G_R_filt.fastq", verbose=TRUE)
derepR50s <- derepFastq("/Volumes/Drive/Data/filtered/A53-N8-G_R_filt.fastq", verbose=TRUE)
derepR51s <- derepFastq("/Volumes/Drive/Data/filtered/A53-N9-G_R_filt.fastq", verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# save dereplication as RDS files:
saveRDS(derepFs, "derepF.rds")
saveRDS(derepRs, "derepR.rds")

#Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=FALSE, pool=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE, pool=TRUE)
#Inspecting the returned dada-class object:
dadaFs[[1]]
dadaRs[[1]]

# Save sequence-variant inference as RDS files which may be uploaded in case R crashes (or save session, see above): 
saveRDS(dadaFs, "dadaFs.rds")
saveRDS(dadaRs, "dadaRs.rds")

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#save mergers
saveRDS(mergers,"mergers.rds")


#Constructing the sequence table
#We can now construct a sequence table of the samples that is analogous to the OTU table
#produced by classical methods.
#Construct sequence table:
seqtab <- makeSequenceTable(mergers)
#Check the dimensions of the sequence table
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#The sequence table is a matrix with rows corresponding to (and named by) the samples, 
#and columns corresponding to (and named by) the sequence variants. 
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, 
#and may be worth removing, e.g.
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)])
#This is analogous to cutting a band in-silico to get amplicons of the targeted length.

# Save sequence table
saveRDS(seqtab, "seqtab.rds")

#Remove chimeras
#The core dada method removes substitution and indel errors, 
#but chimeras remain. Fortunately, the accuracy of the sequences after denoising makes 
#identifying chimeras easier than it is when dealing with fuzzy OTUs: 
#all sequences which can be exactly reconstructed as a bimera (two-parent chimera) from more abundant sequences.
#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. Here chimeras make up about XX% of the inferred sequence variants, 
#but those variants account for only about XX% of the total sequence reads.
#Most of your reads should remain after chimera removal 
#(it is not uncommon for a majority of sequence variants to be removed though). 
#If most of your reads were removed as chimeric, upstream processing may need to be revisited. 
#In almost all cases this is caused by primer sequences with ambiguous nucleotides that 
#were not removed prior to beginning the DADA2 pipeline.

# Save table with the non-chimeric sequences as rds-file:
saveRDS(seqtab.nochim, "seqtab_nochim.rds")


#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
#This is a great place to do a last sanity check. 
#Outside of filtering (depending on how stringent you want to be) 
#there should no step in which a majority of reads are lost. 
#If a majority of reads failed to merge, you may need to revisit the  
#truncLen parameter used in the filtering step and make sure that the 
#truncated reads span your amplicon. If a majority of reads were removed 
#as chimeric, you may need to revisit the removal of primers, as the ambiguous 
#nucleotides in unremoved primers interfere with chimera identification.
write.table(track,"track.txt",sep="\t",col.names = NA)

#Assign taxonomy
#It is common at this point, especially in 16S/18S/ITS amplicon sequencing, 
#to classify sequence variants taxonomically. The DADA2 package provides a native implementation 
#of the RDPs naive Bayesian classifier for this purpose. 
#The assignTaxonomy function takes a set of sequences and a training set of taxonomically classified sequences, 
#and outputs the taxonomic assignments with at least minBoot bootstrap confidence.
#Assign taxonomy:

taxa_silva <- assignTaxonomy(seqtab.nochim, "/Volumes/Drive/silva_nr99_v138.1_train_set.fa.gz")
taxa_silva_revcomp <- assignTaxonomy(seqtab.nochim, "/Volumes/Drive/silva_nr99_v138.1_train_set.fa.gz", tryRC=TRUE)

taxa_GTDB <- assignTaxonomy(seqtab.nochim, "/Volumes/Drive/GTDB_bac-arc_ssu_r86.fa.gz")
taxa_GTDB_revcomp <- assignTaxonomy(seqtab.nochim, "/Volumes/Drive//GTDB_bac-arc_ssu_r86.fa.gz", tryRC=TRUE)

#taxa_silva
#colnames(taxa_silva)
#taxa_silva_revcomp
#colnames(taxa_silva_revcomp)

taxa_silva <- addSpecies(taxa_silva, "/Volumes/Drive/silva_species_assignment_v138.1.fa.gz")
taxa_silva_revcomp <- addSpecies(taxa_silva_revcomp, "/Volumes/Drive/silva_species_assignment_v138.1.fa.gz")

taxa_GTDB <- addSpecies(taxa_GTDB, "/Volumes/Drive/GTDB_dada2_assignment_species.fa.gz")
taxa_GTDB_revcomp <- addSpecies(taxa_GTDB_revcomp, "/Volumes/Drive/GTDB_dada2_assignment_species.fa.gz")

colnames(taxa_silva) <- c("domain", "phylum", "class", "order", "family", "genus", "species")
unname(head(taxa_silva))
#head(taxa_silva)
colnames(taxa_silva_revcomp) <- c("domain", "phylum", "class", "order", "family", "genus", "species")
unname(head(taxa_silva_revcomp))
#head(taxa_silva_revcomp)
colnames(taxa_GTDB) <- c("domain", "phylum", "class", "order", "family", "genus", "species", "subspecies")
unname(head(taxa_GTDB))
colnames(taxa_GTDB_revcomp) <- c("domain", "phylum", "class", "order", "family", "genus", "species", "subspecies")
unname(head(taxa_GTDB_revcomp))


#Save sequence file:
sink("seqs.fa");cat(paste(">","ASV_",seq(1:dim(seqtab.nochim)[2]),"\n",paste(colnames(seqtab.nochim),"\n",sep=""),sep=""),sep="");sink()
seqtab.final<-seqtab.nochim
colnames(seqtab.final)<-paste("ASV_",seq(1:dim(seqtab.nochim)[2]),sep="")
seqtab.final
#Save the sequence table:
write.csv(t(seqtab.final),"seq_table.csv",quote=F)

#Save the sequence taxonomy:
write.table(data.frame(row.names=paste("ASV_",seq(1:dim(seqtab.nochim)[2]),sep=""),unname(taxa_silva)),"seq_Taxonomy_silva.csv",sep=",",col.names=F,quote=F,na="")
write.table(data.frame(row.names=paste("ASV_",seq(1:dim(seqtab.nochim)[2]),sep=""),unname(taxa_silva_revcomp)),"seq_Taxonomy_silva_revcomp.csv",sep=",",col.names=F,quote=F,na="")
write.table(data.frame(row.names=paste("ASV_",seq(1:dim(seqtab.nochim)[2]),sep=""),unname(taxa_GTDB)),"seq_Taxonomy_GTDB.csv",sep=",",col.names=F,quote=F,na="")
write.table(data.frame(row.names=paste("ASV_",seq(1:dim(seqtab.nochim)[2]),sep=""),unname(taxa_GTDB_revcomp)),"seq_Taxonomy_GTDB_revcomp.csv",sep=",",col.names=F,quote=F,na="")

#Come out of R and check the files. DADA2 is done. 
