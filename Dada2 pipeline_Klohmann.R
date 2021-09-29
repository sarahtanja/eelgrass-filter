# Corinne Klohmann 
# cak268@uw.edu

# eDNA samples 
# 16s illumina MiSeq 

# Bioinformatics 
# Dada2 Pipeline 

#install package 
install.packages("devtools")
library(devtools)
devtools::install_github("benjjneb/dada2", ref="v1.16")
# load dada2 package 
library(dada2); packageVersion("dada2")

#Define the following path variable so that it points 
#to the extracted directory on your machine: 
path <- "/Users/corinneklohmann/Documents/GitHub/Bioinformatics-/fastq_files" 
# CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq 
# and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# inspect read quality profiles
# We start by visualizing the quality profiles of the forward reads:
plotQualityProfile(fnFs[1:2])
# n gray-scale is a heat map of the frequency of 
# each quality score at each base position.

# We generally advise trimming the last few nucleotides to avoid 
# less well-controlled errors that can arise there. 
# These quality profiles do not suggest that any additional trimming 
# is needed. We will truncate the forward reads at position 240 
# (trimming the last 10 nucleotides).

# Now we visualize the quality profile of the reverse reads:
plotQualityProfile(fnRs[1:2])
# The reverse reads are of significantly worse quality, 
# especially at the end, which is common in Illumina sequencing.

# This isn’t too worrisome, as DADA2 incorporates quality information 
# into its error model which makes the algorithm robust to lower quality 
# sequence, but trimming as the average qualities crash will improve 
# the algorithm’s sensitivity to rare sequence variants. 
# Based on these profiles, we will truncate the reverse reads at 
# position 160 where the quality distribution crashes.

# Filter and trim
# Assign the filenames for the filtered fastq.gz files.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


# We’ll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), 
# truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the 
# maximum number of “expected errors” allowed in a read, 
# which is a better filter than simply averaging quality scores.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

##                               reads.in reads.out
## F3D0_S188_L001_R1_001.fastq       7793      7113
## F3D1_S189_L001_R1_001.fastq       5869      5299
## F3D141_S207_L001_R1_001.fastq     5958      5463
## F3D142_S208_L001_R1_001.fastq     3183      2914
## F3D143_S209_L001_R1_001.fastq     3178      2941
## F3D144_S210_L001_R1_001.fastq     4827      4312

# Alternatives: Zymo Research has recently developed a tool called 
# Figaro that can help you choose DADA2 truncation length parameters: 
# https://github.com/Zymo-Research/figaro#figaro


# Learn the Error Rates

# The DADA2 algorithm makes use of a parametric error model (err) 
# and every amplicon dataset has a different set of error rates. 
# The learnErrors method learns this error model from the data, 
# by alternating estimation of the error rates and inference of 
# sample composition until they converge on a jointly consistent solution. 

errF <- learnErrors(filtFs, multithread=TRUE)
## 33514080 total bases in 139642 reads from 20 samples will 
#be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
## 22342720 total bases in 139642 reads from 20 samples will be used for 
# learning the error rates.

# It is always worthwhile, as a sanity check if nothing else, to visualize 
# the estimated error rates:

plotErrors(errF, nominalQ=TRUE)
# The error rates for each possible transition (A→C, A→G, …) are shown. 
# Points are the observed error rates for each consensus quality score.

# Sample Inference
# apply core samp inference algorithm to filtered and trimmed sequence data.

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

## Sample 1 - 7113 reads in 1979 unique sequences.
## Sample 2 - 5299 reads in 1639 unique sequences.
## Sample 3 - 5463 reads in 1477 unique sequences.
## Sample 4 - 2914 reads in 904 unique sequences.
## Sample 5 - 2941 reads in 939 unique sequences.
## Sample 6 - 4312 reads in 1267 unique sequences.
## Sample 7 - 6741 reads in 1756 unique sequences.
## Sample 8 - 4560 reads in 1438 unique sequences.
## Sample 9 - 15637 reads in 3590 unique sequences.

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

## Sample 1 - 7113 reads in 1660 unique sequences.
## Sample 2 - 5299 reads in 1349 unique sequences.
## Sample 3 - 5463 reads in 1335 unique sequences.
## Sample 4 - 2914 reads in 853 unique sequences.
## Sample 5 - 2941 reads in 880 unique sequences.
## Sample 6 - 4312 reads in 1286 unique sequences.
## Sample 7 - 6741 reads in 1803 unique sequences.

#Inspecting the returned dada-class object:
  
dadaFs[[1]]
## dada-class: object describing DADA2 denoising results
## 128 sequence variants were inferred from 1979 input unique sequences.
## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# Merge paired reads
# Now merge forward and reverse reads together to get full denoised sequences
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
##                                                                                                                                                                                                                                                       sequence
## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
##   abundance forward reverse nmatch nmismatch nindel prefer accept
## 1       579       1       1    148         0      0      1   TRUE
## 2       470       2       2    148         0      0      2   TRUE
## 3       449       3       4    148         0      0      1   TRUE
## 4       430       4       3    148         0      0      2   TRUE
## 5       345       5       6    148         0      0      1   TRUE
## 6       282       6       5    148         0      0      2   TRUE


# Construct sequence table
# construct an amplicon sequence variant table (ASV) table, 
# a higher-resolution version of the OTU table produced by traditional methods
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
## [1]  20 293
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
## 
## 251 252 253 254 255 
##   1  88 196   6   2

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
## [1]  20 232
sum(seqtab.nochim)/sum(seqtab)
## [1] 0.964263

# Track reads through the pipeline
# final check of our progress
# look at number of reads that made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
##        input filtered denoisedF denoisedR merged nonchim
## F3D0    7793     7113      6996      6978   6551    6539
## F3D1    5869     5299      5227      5239   5025    5014
## F3D141  5958     5463      5339      5351   4973    4850
## F3D142  3183     2914      2799      2833   2595    2521
## F3D143  3178     2941      2822      2868   2553    2519
## F3D144  4827     4312      4146      4224   3622    3483

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















