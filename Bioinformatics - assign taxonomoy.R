# Corinne Klohmann 
# cak268@uw.edu

# eDNA samples 
# 16s illumina MiSeq 

# Bioinformatics 
# Dada2 Pipeline 
# step 2

# load packages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)
install_phyloseq(branch = "devel")
install.packages("devtools")
install_phyloseq(branch = "github")
install_github("phyloseq/joey711")


# Assign taxonomy
# It is common at this point, especially in 16S/18S/ITS amplicon sequencing, 
# to assign taxonomy to the sequence variants. 
# The DADA2 package provides a native implementation of the naive 
# Bayesian classifier method for this purpose. The assignTaxonomy function 
# takes as input a set of sequences to be classified and a 
# training set of reference sequences with known taxonomy, 
# and outputs taxonomic assignments with at least minBoot bootstrap confidence.
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v132_train_set.fa.gz", 
                       multithread=TRUE)
#Extensions: The dada2 package also implements a method to make species level 
# assignments based on exact matching between ASVs and sequenced reference 
# strains. Recent analysis suggests that exact matching (or 100% identity) 
# is the only appropriate way to assign species to 16S gene fragments. 
# Currently, species-assignment training fastas are available for the 
# Silva and RDP 16S databases. To follow the optional species addition step, 
# download the silva_species_assignment_v132.fa.gz file, and place it in 
# the directory with the fastq files.
taxa <- addSpecies(taxa, "~/tax/silva_species_assignment_v132.fa.gz")

# Let’s inspect the taxonomic assignments:

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
##      Kingdom    Phylum          Class         Order           Family          
## [1,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
## [2,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
## [3,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
## [4,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
## [5,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
## [6,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
##      Genus         Species
## [1,] NA            NA     
## [2,] NA            NA     
## [3,] NA            NA     
## [4,] NA            NA     
## [5,] "Bacteroides" NA     
## [6,] NA            NA

# Evaluate accuracy
#Evaluating DADA2’s accuracy on the mock community:

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), 
    "sample sequences present in the Mock community.\n")
## DADA2 inferred 20 sample sequences present in the Mock community.
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), 
    "were exact matches to the expected reference sequences.\n")
## Of those, 20 were exact matches to the expected reference sequences.

#### end regular pipeline 

# Bonus: Handoff to phyloseq
# The phyloseq R package is a powerful framework for further analysis 
# of microbiome data. We now demonstrate how to straightforwardly import 
# the tables produced by the DADA2 pipeline into phyloseq. We’ll also add 
# the small amount of metadata we have – the samples are named by the 
# gender (G), mouse subject number (X) and the day post-weaning (Y) 
# it was sampled (eg. GXDY).

#Import into phyloseq:

library(phyloseq); packageVersion("phyloseq")
## [1] '1.30.0'
library(Biostrings); packageVersion("Biostrings")
## [1] '2.54.0'
library(ggplot2); packageVersion("ggplot2")
## [1] '3.3.0'
theme_set(theme_bw())

# We can construct a simple sample data.frame from the information 
# encoded in the filenames. 
# Usually this step would instead involve reading the sample data in from a file.
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
# We now construct a phyloseq object directly from the dada2 outputs.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

# we’ll store the DNA sequences of our ASVs in the refseq slot of the phyloseq 
# object, and then rename our taxa to a short string. That way, the short 
# new taxa names will appear in tables and plots, and we can still recover 
# the DNA sequences corresponding to each ASV as needed with refseq(ps).
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 232 taxa and 19 samples ]
## sample_data() Sample Data:       [ 19 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 232 taxa by 7 taxonomic ranks ]
## refseq()      DNAStringSet:      [ 232 reference sequences ]

# Visualize alpha-diversity:

plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")

# Ordinate:

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

# Bar plot:

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")