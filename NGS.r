####
#### Next generation sequencing (NGS) represents new technology for obtaining massively-
#### parallel reads of genomic sequences. NGS has replaced Sanger sequencing as the 
#### current state-of-the-art in sequencing technology. NGS allows for much greater 
#### throughput at much lower cost, compared to Sanger sequencing. While NGS is 
#### fundamentally about obtain genomic sequences, it can be applied in a variety of 
#### ways, including transcription analysis, epigenetics, variant discovery, and meta-
#### genomics. A good overview of NGS technology is in the 'metzger_10.pdf' article on 
#### eCampus.
####

####
#### Sequence Read Archive (SRA) database access. SRA is a database of short sequence 
#### reads from various NGS technologies. 
####

## We will use the 'SRAdb' package.
source("http://bioconductor.org/biocLite.R")
biocLite("SRAdb")
library(SRAdb)

## We need the SQLite file that manages database queries to the SRA database. We download 
## and uncompress the file using the 'getSRAdbFile' function. This takes a while, because 
## the file is quite large.
sqlFile <- getSRAdbFile()

## We then create a connection between the SQLite file and the SRA database.
sraCon <- dbConnect(SQLite(), sqlFile)

## We obtain a list of all the available tables in the database with the 'dbListTables' 
## function. The available fields for a given table are obtained with 'dbListFields'.
dbListTables(sraCon)
dbListFields(sraCon, "study")

## The dbGetQuery function uses SQL query language to search metadata (such as titles,  
## abstracts, and study accession IDs). Here, we search for studies with keywords equal 
## to 'embryo'. The 'paste' function is used to break the search string over two lines.
##
## Here's some information on the one study that is returned:
##
##   http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17621
embryo <- dbGetQuery(sraCon, paste("select study_accession,study_title from study where", 
  "study_description like'%embryo'", sep = " "))
embryo

## The getSRA function enables "free text" searches of the SRA database. Free text search 
## means that the actual content of the files are searched, as opposed to just the files'
## metadata; the dbGetQuery function only searches metadata. Here, we search for database 
## entries with the word "brain" in them. We ask for output variables that are relevant 
## to the "runs" involved as well as the studies. We get about 23k results.
brain <- getSRA(search_terms = "brain", out_types = c('run','study'), sraCon)
head(brain)

## Can combine terms using logical operators. This time we search for database entries 
## that have either "alzheimers" or "epilepsy" in them, and we ask for output variables 
## that are relevant to the samples. We get about 12k results.
alz_epi <- getSRA(search_terms ='Alzheimers OR Epilepsy', out_types = 'sample', sraCon)

####
#### Downloading SRA data. Above, we conducted queries of the SRA database. Once we have 
#### identified a data set of interest, we can download its FASTQ file.
####
#### FASTQ files. FASTQ is a file format that is widely used in sequencing. A FASTQ file 
#### for a sequence consists of four lines: (1) the sequence name; (2) the sequence; (3) 
#### optional information about the sequence; (4) confidence scores for each called base. 
#### Confidence is quantified in terms of something called the Phred score. In sequencing, 
#### the probability that any given called base is correct can be estimated. The Phred 
#### score transforms this estimated probability to make it be on a more manageable scale 
#### by log transformation. Specifically, with P the estimated probability of a base call 
#### being correct, the Phred score is Q = -10 log_10 P. As an example, if P = 1e-5 
#### (0.00001, corresponding to an accuracy of 99.999%), then Q = 50. The Phred scores in 
#### a FASTQ file are represented by ASCII characters, to make the data representation 
#### more compact. Examples of ASCII character representations:
####
####   http://www.somewhereville.com/?p=1508
####

## The first two sample IDs returned in our 'Alzheimers OR Epilepsy' search above were 
## ERS354366 and SRS266589. Given these study IDs, sraConvert can be used to obtain other 
## descriptors: submission ID, study ID, experiment ID, and run ID.
sra_conversion <- sraConvert(c('ERS354366','SRS266589'), sra_con = sraCon)

## Given any of a submission ID, study ID, experiment ID, or run ID, the getSRAfile 
## function is for downloading the corresponding FASTQ files. Unfortunately, this 
## function does not run for me. 
getSRAfile(c("SRR351672", "SRR351673"), sraCon, fileType ='fastq')

## Alternatively, we can use getSRAinfo to obtain an FTP address for downloading a .sra 
## file. A .sra file can be converted into a FASTQ file using the SRA Toolkit:
##
##   http://www.ncbi.nlm.nih.gov/books/NBK158900/
##
## I used the SRA Toolkit to convert the data from one "run" (sample) of the experiment 
## with SRA code SRX100465 into a FASTQ file. Searching for this experiment code on the 
## SRA site, we find that these are ChIP-Seq experiments. We are considering only one of 
## 797 runs in this experiment. To keep the FASTQ file size down (otherwise, it would be 
## about 5GB in size), I only extracted the first 10k reads. The resulting file is on 
## eCampus, as "SRR351672.fastq".
sra_info <- getSRAinfo("SRX100465", sraCon, sraType = "sra")

####
#### Reading FASTQ files. FASTQ is the standard file format for storing short sequence 
#### reads from NGS applications. A FASTQ file contains the sequence and meta data for 
#### all the reads in a single run (run = sample). For each read, there are four lines of 
#### text. Here's an example entry for one read:
####
#### @SRR351672.751 ILLUMINA-EAS45_1:1:1:1089:18662:0:1:1 length=36
#### ACTGCGGCAAGGAATATAGAACATGAGTCTGGAACA
#### +SRR351672.751 ILLUMINA-EAS45_1:1:1:1089:18662:0:1:1 length=36
#### CCCCCBCCCDCBCCBC@CBCCCCCCCCCCCCCCCC?
####
#### The first line contains an identifier for the read as well as platform-specific data 
#### on where on the device the read came from. This read came from an Illumina 
#### instrument. The second line is the sequence. The third line is optional and may 
#### contain additional description of the read. The fourth line contains quality codes 
#### for each of the bases in the read.
####
####   http://en.wikipedia.org/wiki/FASTQ_format
####

## We'll use the ShortRead package for working with FASTQ files.
source("http://bioconductor.org/biocLite.R")
biocLite("ShortRead")
library(ShortRead)  

## Read the example FASTQ file from above. This creates an object of class ShortReadQ. We 
## can access the elements of each read entry using the sread (for reading the sequences 
## themselves), id (for reading the identifying information), and quality (for reading 
## the quality scores) functions. The first 750 or so reads appear to be of poor quality; 
## the first few hundred, in particular, have very few actual calls, with the majority of 
## the attempted calls shown as just "N" (for unknown).
fastq <- readFastq("../data/NGS_Data/SRR351672.fastq")
sread(fastq)[-(1:750)]
id(fastq)[-(1:750)]
quality(fastq)[-(1:750)]

## You can also use the generic readLines function to read lines of a text file. Here are 
## the first four lines of our FASTQ file, corresponding to the entries for the first 
## read. It returns a character vector of length four.
readLines("../data/NGS_Data/SRR351672.fastq", 4)

####
#### Alignment data. Unless we are doing whole genome sequencing (assembling a genome 
#### from scratch), we will want to align our reads with a reference genome. For example, 
#### in RNA-Seq, we sequence RNA fragments and align those with the genome, so that we 
#### know which genes the fragments belong to. Since reference genomes are typically very 
#### large, alignment of NGS reads is a computationally-intensive task. There are many 
#### standalone programs for alignment; BWA and Bowtie are widely used.
####
####   http://en.wikibooks.org/wiki/Next_Generation_Sequencing_%28NGS%29/Alignment
####
#### The standard file representation of aligned NGS data is the sequence alignment map 
#### (SAM) file format; a BAM file is the binary version of a SAM file. SAM files include 
#### a variety of variables for each aligned read, including the sequence itself, the 
#### location in the genome where it was deemed to align, and a quality score for the 
#### alignment (a "CIGAR" code).
####
####   http://genome.sph.umich.edu/wiki/SAM
####

## We can use the Rsamtools package to import and manipulate SAM / BAM files. The scanBam 
## function loads a BAM file. We'll use an example BAM file from the Rsamtools package.
library(Rsamtools)

## This locates the example BAM file included in the extdata folder of the Rsamtools 
## package installation. There are 3,307 records in this BAM file.
fl <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
countBam(fl)

## The scanBam function returns a list of lists, one entry per BAM file. Since we only 
## read one BAM file, we just take the first element of the list. 
bam <- scanBam(fl)[[1]]
names(bam)

## We can read in just a subset of the attributes if we like.
param <- ScanBamParam(what = c("rname", "strand", "pos", "qwidth", "seq"))
bam_1 <- scanBam(fl, param = param)[[1]]
lapply(bam_1, head)

## Can convert the BAM file object into a data frame. This can be easily searched for 
## aligned reads of the desired criteria.
bam_df <- do.call("DataFrame", bam)
ii <- (1:nrow(bam_df))[bam_df$pos <= 100 & bam_df$qwidth == 35]
ii <- ii[!is.na(ii)]
bam_df[ii, ]

####
#### Analyzing RNA-Seq data with edgeR. RNA-Seq data (and data from related technologies) 
#### consists of read counts for thousands of "tags", where ideally each tag corresponds 
#### to a gene. A typical application is to identify tags that differ in abundance 
#### between comparison groups. Consider a two-class scenario, where there are n_i 
#### samples in the ith comparison group. The (ij)th sample "library" consists of m_{ij} 
#### reads, y_{ijk} of which correspond to tag k. It is natural to model these counts via 
#### a Binomial model or a Poisson approximation to the Binomial. In the Binomial case, 
#### we would say the counts are Binomial with probabilities lambda_{ik}. In the Poisson 
#### case, we would say the counts are Poisson with means m_{ij} lambda_{ik}. Our 
#### interest is in testing the null hypothesis that lambda_{1k} = lambda_{2k} = ... = 
#### lambda_{Gk}, where G is the number of comparison groups. 
####
#### One challenge with using the Binomial or Poisson models for NGS count data is that 
#### the data are frequently "overdispersed" with respect to these models. That is, their 
#### variation exceeds that expected by the Binomial or Poisson distributions. A commonly-
#### used solution to this problem is to use a Poisson model that allows for over-
#### dispersion. One such model is the Negative Binomial, which is equivalent to a 
#### hierarchical Poisson model that allows the Poisson rate parameters to vary according 
#### to a Gamma distribution. 
####
#### A more general challenge in differential expression analysis is high variability in 
#### feature-specific dispersion parameter estimates. For example, in microarray-based 
#### differential expression analysis, suppose we are comparing just two groups. We might 
#### use a two-sample t-test with equal variances on the intensities for feature k. To 
#### compute the t statistic, we have to estimate a pooled standard deviation, s_k. The 
#### s_k are known to be highly variable, and this can effect the precision of the 
#### overall analysis. A solution that is commonly used is to use a "moderated" version 
#### of the test statistic, where the dispersion parameters are all computed together and 
#### shrunken toward a common value. This results in dispersion parameter estimates that 
#### are much more stable, and this has been shown to lead to more precise conclusions. 
#### The moderated test statistic approach was originally developed in the context of 
#### microarray-based differential expression analysis. When we talk about microarrays 
#### later this semester, we will use the limma package for differential expression 
#### analysis with moderated test statistics.
####
####   http://rileylab.bio.umb.edu/sites/g/files/g1314676/f/201403/ebayes.pdf
####
#### The edgeR package implements a Negative Binomial model with moderated test statistics
#### for the analysis of count-based tag data (most often applied to RNA-Seq data and 
#### differential expression analysis). 
####
####   http://bioinformatics.oxfordjournals.org.lib-ezproxy.tamu.edu:2048/content/26/1/139.long
####

## We need the edgeR library. We will also install the goseq library so we can use one of 
## its example data sets.
biocLite("edgeR")
biocLite("goseq")
library(edgeR)
library(goseq)

## Read example data from goseq. Here is the description of the data from the goseq 
## vignette:
##
## "This experiment examined the effects of androgen stimulation on a human prostate 
## cancer cell line, LNCaP. The data set includes more than 17 million short cDNA reads 
## obtained for both the treated and untreated cell line and sequenced on Illumina’s 1G 
## genome analyzer. For each sample we were provided with the raw 35 bp RNA-seq reads 
## from the authors. For the untreated prostate cancer cells (LNCaP cell line) there were 
## 4 lanes totaling 10 million 35 bp reads. For the treated cells there were 3 lanes 
## totaling 7 million 35 bp reads. All replicates were technical replicates. Reads were 
## mapped to NCBI version 36.3 of the human genome using bowtie. Any read that  mapped to 
## multiple locations was discarded. Using the ENSEMBL 54 annotation from biomart, each 
## mapped read was associated with an ENSEMBL gene. This was done by associating any read 
## that overlapped with any part of the gene (not just the exons) with that gene. Reads 
## that did not correspond to genes were discarded."
lncap <- read.table(system.file("extdata", "Li_sum.txt", package='goseq'), sep = '\t', 
  header = TRUE, stringsAsFactors = FALSE, row.names = 1)
head(lncap)

## There are many ENSEMBL genes that had no counts in any of the samples. Remove them 
## from the analysis.
n_count <- drop(as.matrix(lncap) %*% rep(1, 7))
lncap <- lncap[n_count > 0, ]
m <- nrow(lncap)
n <- ncol(lncap)

## Now create a DGElist object for holding the counts along with the library sizes (the 
## total number of counts) and the group membership information for each sample. The 
## first four samples (run on the first four "lanes" of an Illumina flow cell device) are 
## control samples, and the last three samples are treated samples. The DGEList object is 
## a list with two components, one for the counts and one containing the sample 
## information.
lncap_dge <- DGEList(lncap, lib.size = colSums(lncap), group = 
  factor(rep(c("Control", "Treatment"), times = c(4, 3))))
lncap_dge
lncap_norm<- calcNormFactors(lncap_dge)
## Estimate the "common dispersion parameter." This is the number toward which all gene-
## specific dispersion estimates will be shrunken.
lncap_disp <- estimateCommonDisp(lncap_dge)

## edgeR does an exact test via the exactTest function, where a p-value is computed for 
## each feature as the probability under the NB model of observing counts as or more 
## extreme as the ones observed for that feature. Thus, in this case, because there are 
## 22,743 genes, we get 22,743 p-values. The topTags function then extracts the top DE 
## genes by p-value rank. Also reported in the output is the estimated false discovery 
## rate (FDR) for each gene. We will discuss FDR more later in the semester. For now, 
## just note that FDR provides a confidence measure appropriate for conducting thousands 
## of hypothesis tests simultaneously. It is defined as the expected proportion of tests 
## with p-values equal to or less than that for a particular gene are expected to be null 
## (not differentially expressed). If it is small, it means we can be confident that the 
## gene is differentially expressed.
lncap_test <- exactTest(lncap_disp)
lncap_out <- topTags(lncap_test, n = 100, sort.by = "PValue")

head(lncap_out)

## Look at the read counts for the selected genes.
sig_genes <- rownames(lncap_out$table)
lncap[sig_genes, ]

## Just out of curiosity, what would happen if we simulated read counts for which there 
## is *no* differential expression? We'll allow the different genes to have their own 
## means, but the means do not differ between comparison groups. We'll also randomly 
## generate the dispersion parameters. The edgeR p-value distribution does not look 
## uniform. I also try doing a simple two-sample t-test. Its p-values also do not look 
## uniform.
Y <- matrix(NA, nrow = m, ncol = n)
for(i in 1:m) {
  mu_i <- round(runif(1, 0, 500), 0)
  phi_i <- runif(1, 0, 0.5)
  
  Y[i, ] <- rnbinom(n, mu = mu_i, size = phi_i)
}

Y_dge <- DGEList(Y, lib.size = colSums(Y), group = 
  factor(rep(c("Control", "Treatment"), times = c(4, 3))))
Y_disp <- estimateCommonDisp(Y_dge)
Y_test <- exactTest(Y_disp)

hist(Y_test$table$PValue, main = "edgeR P-Value Distribution (Should Be Uniform)")

## Removing the genes with very low read counts helps with the p-values near to or equal 
## to one, but there is still a peak near zero.
ii_keep <- (1:m)[rowSums(Y) >= 10]
hist(Y_test$table$PValue[ii_keep])

## T-test p-values also look non-uniform. But, at least they don't have a spike near zero.
T_pvals <- rep(NA, m)
for(i in 1:m) {
  T_pvals[i] <- t.test(Y[i, 1:4], Y[i, 5:7], var.equal = TRUE)$p.value
}

hist(T_pvals, main = "T-Test P-Value Distribution")

####
#### DESeq2 is another widely-used algorithm and suite of software for the differential 
#### analysis of NGS count data. DESeq2 is similar to edgeR, in that it also uses 
#### "moderated" Negative Binomial models for its analysis. While they are similar, 
#### DESeq2 and edgeR will give slightly different results, although there will typically 
#### be a lot of overlap between their results (e.g., between the lists of genes called 
#### differentially expressed). It is claimed that DESeq2 will tend to call more genes 
#### differentially expressed than edgeR. The differences between the two methods can 
#### most likely be attributed to differences in normalization techniques used, as well 
#### as minor differences in the models used. Both methods can handle arbitrary numbers 
#### of comparison groups and timecourse studies. The use of DESeq2 to conduct 
#### differential expression analysis of RNA-Seq data is summarized in a Bioconductor 
#### workflow:
####
####   http://www.bioconductor.org/help/workflows/rnaseqGene/
####
#### The workflow shows how to construct counts for genes, DNA binding sites, etc., 
#### starting from SAM / BAM alignment files. We will proceed with the already-prepared 
#### count data we used with edgeR above.
####

biocLite("DESeq2")
library(DESeq2)

## The main object class for DESeq2 is DESeqDataSet. This is analogous to edgeR's DGEList.
col_data <- data.frame("Group" = factor(rep(c("Control", "Treatment"), times = c(4, 3))), 
  row.names = colnames(lncap))
lncap_deseq <- DESeqDataSetFromMatrix(countData = lncap, colData = col_data, 
  design = ~ Group)
  
## For exploratory analysis of NGS counts, DESeq2 recommends applying its 'rlog' 
## transformation. This is a modified version of a log2 transformation, intended to make 
## the features have approximately the same variance. Untransformed read counts have 
## variances that increase with the mean. Thus, for example, PCA analysis of the raw 
## counts would likely depend almost entirely on the features with the highest counts. 
## The rlog function returns an object of class SummarizedExperiment. This is a general 
## class for representing -omics data, where the rows correspond to features and the 
## columns correspond to samples. We can extract the numeric values for each feature with 
## the assay function.
lncap_rld <- rlog(lncap_deseq)
head(assay(lncap_rld))

## One use of the rlog-transformed data is the clustering of samples. We would expect, 
## for example, that samples from the same comparison groups are more similar than samples 
## from different comparison groups. We can also get insight into any possibly outlying 
## samples (perhaps of poor quality). We'll use the heatmap.2 function from the gplots 
## package to create a clustering heatmap. We'll use the RColorBrewer package to create 
## some fancy colors to use in the heatmap. The dist function computes Euclidean 
## distances between the rlog-transformed data. The hclust function does hierarchical 
## clustering, given the sample distances.
##
## As expected, the first four control samples (lanes) cluster together, and the last 
## three treatment samples cluster together.
biocLite("gplots")
library(gplots)
library(RColorBrewer)

sampleDists <- dist(t(assay(lncap_rld)))
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
hc <- hclust(sampleDists)
heatmap.2(as.matrix(sampleDists), Rowv = as.dendrogram(hc), symm = TRUE, trace = "none", 
  col = colors, margins = c(2,10), labCol = FALSE)
  
## Another way to visualize sample differences is with a principal component analysis 
## (PCA). As we saw earlier in the semester, PCA spreads out the sample data (with 
## dimensionality equal to the number of features) on a lower dimension, typically 2D. 
## DESeq2 has a built-in function for this: plotPCA. Here, it is clear that, by far, the 
## biggest source of difference between the samples is the control vs. treatment 
## distinction.
plotPCA(lncap_rld, intgroup = "Group")

## Differential expression (DE) analysis is done with the DESeq function. Note that some 
## p-values are reported as NA. This can happen for two reasons: (1) there is an outlying 
## value; (2) the overall mean count is too low. The padj column contains FDR estimates.
lncap_de <- DESeq(lncap_deseq)
lncap_de_out <- results(lncap_de)
summary(lncap_de_out)

## The lncap_de_out object is of class DataFrame. The function mcols extracts the "meta 
## data" for the columns of a DataFrame object.
mcols(lncap_de_out, use.names = TRUE)

## We can extract the top DE genes as ordered by log fold change, among genes with small 
## FDR. We first show the top DE genes as ordered by negative log fold change (expression 
## in the control group is greater than that in the treated group). Then, the top DE 
## genes as ordered by positive log fold change (expression in the treated group is 
## greater than that in the control group).
lncap_de_out_sig <- subset(lncap_de_out, padj < 0.1)
head(lncap_de_out_sig[order(lncap_de_out_sig$log2FoldChange), ])
head(lncap_de_out_sig[order(-lncap_de_out_sig$log2FoldChange), ])

## We can look at the raw counts for a particular gene (say, one with a small FDR). The 
## plotCounts function creates a plot to visualize the normalized version of the counts.
## DESeq normalizes in a way that factors in the library size (total number of reads) for 
## each sample. Counts in a sample with a library size much larger than those of the 
## other samples are shrunken to more comparable with the counts for the other samples.
top_gene <- rownames(lncap_de_out)[which.min(lncap_de_out$padj)]
assay(lncap_deseq)[match(top_gene, rownames(assay(lncap_deseq))), ]
plotCounts(lncap_de, gene = top_gene, intgroup = "Group")

## We can also make an MA plot using DESeq2's plotMA function. Without the ylim argument, 
## some of the log fold changes are truncated and shown as a flat line at the extreme y-
## axis values (for some reason). Here, we circle the point for the gene with the lowest 
## FDR and label it with its Ensembl gene name.
plotMA(lncap_de_out, ylim = c(-5, 5))
with(lncap_de_out[top_gene, ], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, top_gene, pos = 2, col = "dodgerblue")
})

## As we'll see when we discuss FDRs in more detail, the distribution of the p-values is 
## of interest. We are interested in the p-values themselves (not the "adjusted" p-
## values from DESeq, which are actually FDR estimates). The peak in the middle is un-
## expected. It turns out that this is due to genes with very small counts.
hist(lncap_de_out$pvalue)

quantile(lncap_de_out$baseMean)
hist(lncap_de_out$pvalue[lncap_de_out$baseMean > 0.65])

####
#### Yet another method for analyzing NGS data is limma. The limma method and software 
#### was originally developed for gene expression microarrays, where the response is an
#### essentially continuous intensity value. With NGS data, the response is a count. 
#### However, using a particular transformation (voom), NGS count data can be turned into 
#### continuous data, and all the mature tools of limma can be employed. In particular, 
#### limma allows for random effects modeling and gene set testing, while the negative 
#### binomial-based methods do not. We will use limma more extensively when we discuss
#### microarrays.
#### 

## We will use an example data set from the pasilla package. These are RNA-Seq data in 
## which two groups were compared in Drosophila melanogaster (the common fruit fly). The 
## control group had no treatment administered to it. In the treatment group RNAi (RNA 
## interference) was used to inhibit the expression of the pasilla (PS) gene. The 
## interest was in finding genes that are differentially expressed under the two 
## conditions. 
##
## We will use the DESeq package (a precursor to DESeq2) to extract the count information 
## from the pasilla data. The limma package was already installed as a dependency of a 
## previous Bioconductor package that we used.
biocLite(c("DESeq", "pasilla"))
library(limma)
library(DESeq)
library(pasilla)

## The pasillaGenes object contains gene counts for the 7 samples. It is an object of 
## class CountDataSet, which is defined by the DESeq package. DESeq's count function can 
## be used to extract the gene-specific read counts. Note that the first three samples 
## are from the treatment group, and the last four samples are from the control group. 
## We'll rename the columns to be shorter.
data(pasillaGenes)
pasillaGenes

ps_counts <- counts(pasillaGenes)
head(ps_counts)
colnames(ps_counts) <- c(paste("T", 1:3, sep = "_"), paste("C", 1:4, sep = "_"))

## Limma and voom require a design matrix for representation of which samples correspond 
## to which comparison groups. A design matrix is n x p, where n is the number of samples 
## and p is the number of parameters to be estimated. In the context of simple linear 
## regression, the design matrix would have n rows and two columns. The first column 
## would be all ones, for the intercept, and the second column would contain the 
## covariate values (the values of x in y = beta_0 + beta_1 x + epsilon). In our example, 
## we just have two comparison groups, so x is binary.
des_mat <- cbind(1, c(1, 1, 1, 0, 0, 0, 0))

## Now do the voom transformation and pass the result to limma's lmFit function for 
## fitting a linear model to each gene. The linear model is the basis for our differential 
## expression test (we're essentially testing whether the beta_1 from above = 0, for each 
## gene). The eBayes function does the actual test, implementing an Empirical Bayes model 
## that we will discuss in more detail when we get to microarrays.
ps_voom <- voom(ps_counts, design = des_mat, plot = FALSE)
ps_fit <- lmFit(ps_voom, design = des_mat)
ps_eb <- eBayes(ps_fit)

## The topTable function can then be used to extract the genes with small estimated FDRs.
## There are 445 genes with FDR < 0.05.
ps_FDRs <- topTable(ps_eb, n = nrow(ps_counts), coef = 2)
ps_DE <- rownames(ps_FDRs[which(ps_FDRs$adj.P.Val < 0.05), ])
length(ps_DE)
head(ps_FDRs[ps_DE, ])

## We can further filter the results by requiring both that FDR < 0.05 and at least a 
## two-fold difference in mean expression (so, abs(log(fold change)) > log(2)). There are 
## 297 such genes.
rownames(ps_FDRs[which(ps_FDRs$adj.P.Val < 0.05 & abs(ps_FDRs$logFC) > log(2)), ])

## Here is a "volcano" plot for the results. A volcano plot puts log fold change on the 
## x-axis and -log10(p-value) on the y-axis. Instead of p-values, we use estimated FDRs.
## Here, we color-code the points based on whether they correspond to genes with FDR less 
## than 0.05. The red points would be called differentially expressed.
clr <- rep("black", nrow(ps_FDRs))
clr[which(ps_FDRs$adj.P.Val < 0.05)] <- "red"
plot(ps_FDRs$logFC, -log10(ps_FDRs$adj.P.Val), col = clr, xlab = "log fold change", 
  ylab = "-log(FDR)")
abline(h = -log10(0.05), col = "blue", lwd = 2)

####
#### GO enrichment of RNA-Seq data. As we have seen earlier in the semester, GO enrichment
#### analysis involves inspecting a set of "interesting" features for gene ontology (GO) 
#### terms that present in greater frequency than would be expected by chance. 
####
#### There is a unique challenge when performing GO enrichment on RNA-Seq (or any count-
#### based) data. Longer reads tend to have higher read counts, even if those reads are 
#### expressed at the same levels as shorter reads. Because of the higher read count, 
#### there is more statistical power to detect DE in genes with longer reads. Similarly, 
#### there is more statistical power to detect DE in genes with high abundance. When 
#### carrying out GO enrichment analysis, GO terms that are involved in genes with long 
#### reads and / or genes that are DE will tend to be preferentially selected. It is 
#### therefore desirable to somehow adjust for these biases when doing the GO enrichment 
#### analysis.
####
#### The GOseq method creates a "probability weighting function" (PWF) that estimates the 
#### probability that a gene with reads of a particular length will be selected as DE. 
#### For assessing the statistical significance of a GO term's enrichment, genes are 
#### randomly selected as DE, using the PWF as a weighting function in the selection 
#### process. The GO enrichment algorithm is run on each random selection of DE genes, 
#### and a p-value is calculated as the proportion of randomization-based enrichment 
#### statistics that exceed our observed statistic.
####
####   http://genomebiology.com/2010/11/2/r14
####

## We will use the goseq and edgeR packages.
library(goseq)
library(edgeR)

## We will again use the LNCaP data set from the goseq package. We will also, again, 
## remove genes with zero counts.
lncap <- read.table(system.file("extdata", "Li_sum.txt", package='goseq'), sep = '\t', 
  header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  
n_count <- drop(as.matrix(lncap) %*% rep(1, 7))
lncap <- lncap[n_count > 0, ]

## Again create a DGEList object for use by edgeR. Recall that the first four samples are 
## control samples, and the last three samples are treatment samples.
lncap_dge <- DGEList(lncap, lib.size = colSums(lncap), group = 
  factor(rep(c("Control", "Treatment"), times = c(4, 3))))

## Now carry out and edgeR differential expression analysis to extract "interesting" 
## genes from the list of all genes. Here, we extract the top 100 genes, as ordered by 
## their p-values.
lncap_disp <- estimateCommonDisp(lncap_dge)
lncap_test <- exactTest(lncap_disp)
lncap_out <- topTags(lncap_test, n = 100, sort.by = "PValue")

## To create the PWF, we need a vector of ones and zeros for each of the genes, with ones 
## indicating which genes were called DE.
lncap_DE <- rep(0, nrow(lncap_test$table))
names(lncap_DE) <- rownames(lncap_test$table)
lncap_DE[match(rownames(lncap_out), rownames(lncap_test$table))] <- 1
table(lncap_DE)

## Now use the nullp function to create the PWF. The genome argument is a code for which 
## genome to use. For a list of supported genomes, use supportedGenomes. The id argument 
## is a code for which gene identifier to use in the genome. FOr a list of supported gene 
## IDs, use supportedGeneIDs. We use the hg19 genome, which is one human genome annotation 
## list, and the ensGene, which is the set of Ensembl gene identifiers. See the GOseq 
## publication (reference above) for interpreting the plot that is automatically made by 
## nullp.
supportedGenomes()
supportedGeneIDs()
lncap_PWF <- nullp(lncap_DE, genome = "hg19", id = "ensGene")
head(lncap_PWF)

## Use the PWF data as input to goseq to perform the actual GO enrichment analysis. The 
## top GO term is GO:0006576. This is a biological process described the amigo site we've 
## used earlier as "cellular biogenic amine metabolic process". 
lncap_GO_enrich <- goseq(lncap_PWF, genome = "hg19", id = "ensGene")
head(lncap_GO_enrich)

####
#### KEGG enrichment of RNA-Seq data. The goseq package can again be used.
####

## We do the same several initial steps from the GO enrichment code above, through the 
## computation of the PWF. Then, to perform KEGG enrichment, we need the org.Hs.eg.db 
## package to (1) translate the Ensembl gene identifiers from the LNCaP data set into 
## Entrez IDs, then (2) obtain all the KEGG IDs associated with these Entrez IDs, and 
## finally to (3) map the KEGG IDs for each Entrez ID back to the corresponding Ensembl 
## IDs.
library(org.Hs.eg.db)

en_eg <- as.list(org.Hs.egENSEMBL2EG)
eg_kegg <- as.list(org.Hs.egPATH)

## This code goes to each Ensembl ID, pulls out the Entrez ID(s) it corresponds to, then 
## returns the KEGG IDs that correspond to those Entrez IDs.
en_kegg_f <- function(eg, eg_kegg) {
  unlist(eg_kegg[eg], use.names = FALSE)
}
en_kegg <- lapply(en_eg, en_kegg_f, eg_kegg)

## We can now use the goseq function. The top KEGG pathway is 00565, which is described 
## as "ether lipid metabolism".
lncap_kegg_enrich <- goseq(lncap_PWF, gene2cat = en_kegg)
head(lncap_kegg_enrich)

####
#### DNA methylation analysis. DNA methylation happens when a methyl group is added to 
#### either a C or an A in a DNA sequence. Methylation of a C is known to reduce the 
#### expression of the corresponding gene. Methylation typically occurs at "CpG" sites; a 
#### CpG site is location in a DNA sequence where a C is next to a G (this is not the 
#### same thing as a C being complementary to a G). Methylation is important to many 
#### biological functions, including cancer progression. Methylation is an example of 
#### "epigenetics," in that it regulates expression in a way that does not require any 
#### change to the actual DNA sequence. 
####
#### The methyAnalysis package that we use below is actually designed for use with 
#### specialized microarrays for methylation sites. For each of many genomic regions in 
#### which methylation is known to occur, there are one or more probes specific to CpG 
#### sites in that region. At each probe, we get different light intensities for genome 
#### fragments at that site that are methylated and not methylated. The output intensity 
#### for a given probe is a relative intensity, comparing the methylated intensity to the 
#### non-methylated intensity. Higher relative intensities correspond to probes where 
#### there is substantial methylation. A region can then be deemed interesting if it has 
#### many methylated probes inside it. A common analysis task is then to identify 
#### differentially methylated regions (DMRs), comparing two or more treatment / control 
#### groups.
####
#### While methyAnalysis is specifically for microarrays , it can be applied to NGS data, 
#### where sequencing is done of DNA regions where methylation is known to occur. However, 
#### this is not ideal, since NGS data are count-based, not continuous like microarray 
#### intensities. The edgeR, DESeq2, and limma packages could be used instead, with some 
#### extra attention paid to how the data are represented.
####
#### One limitation of methyAnalysis is that it can only compare two groups.
####

## Install both the methyAnalysis package and the TxDb.Hsapiens.UCSC.hg19.knownGene 
## package (for annotating the results of our DMR analysis).
biocLite("methyAnalysis", "TxDb.Hsapiens.UCSC.hg19.knownGene")

library(methyAnalysis)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

## We will use an example data set from the methyAnalysis package. There are two 
## comparison groups, named "Type1" and "Type2", with four samples each. The samples are 
## randomly-selected cell lines from two tissues. There are 4,243 features, where each 
## feature is a genomic region in which methylation is known to occur.
data(exampleMethyGenoSet)
slotNames(exampleMethyGenoSet)
head(locData(exampleMethyGenoSet))
dim(exprs(exampleMethyGenoSet))
pData(exampleMethyGenoSet)
str(exampleMethyGenoSet)

## We use the detectDMR.slideWin function to compute p-values (and FDRs) for individual 
## probes. This function smooths the intensities at nearby CpG sites, in an effort to 
## reduce the noisy nature of the data. Then, a two-sample t-test is carried out at each 
## site. The identifySigDMR function selects DNA regions as interesting if they contain 
## probes with small enough FDRs.
DMR_out <- detectDMR.slideWin(exampleMethyGenoSet, sampleType = c(rep("Type1", 4), 
  rep("Type2", 4)))
head(DMR_out)
DMR_sig_out <- identifySigDMR(DMR_out)

## We can obtain annotation information for the regions with differential methylation 
## from UCSC using the TxDb.Hsapiens.UCSC.hg19.knownGene pacakge. 
DMR_annot <- annotateDMRInfo(DMR_sig_out, 'TxDb.Hsapiens.UCSC.hg19.knownGene')
DMR_annot

## One example gene that showed up as differentially methylated in the annotation is SIM2.
## The heatmapByChromosome function plots the relative methylation intensities for the 
## probes located in the SIM2 region.
heatmapByChromosome(exampleMethyGenoSet, gene = '6493', genomicFeature = 
  'TxDb.Hsapiens.UCSC.hg19.knownGene', includeGeneBody = TRUE)

####
#### ChIP-Seq analysis. Methylation is one "epigenetic" mechanism by which gene expression
#### is regulated. Another epigenetic mechanism involves the binding of certain proteins 
#### to the DNA. For example, a transcription factor is a protein that binds to DNA (at a 
#### transcription factor binding site) and regulates the expression level of nearby 
#### genes. Histones are another class of proteins that bind to DNA and regulate its 
#### expression. ChIP-Seq is a technology for locating where in the genome specific 
#### proteins bind. In ChIP-Seq, a protein of interest is mixed with DNA and allowed to 
#### bind to its designated spots along the DNA sequence. The DNA is then fragmented, and 
#### the fragments to which the protein has bound are extracted by immunoprecipitation. 
#### The proteins are detached from the fragments, and NGS is used to sequence the 
#### fragments. Counting the number of reads for each fragment and mapping them back to 
#### the genome allows us to quantify the proclivity of the protein to bind to different 
#### genomic locations. The two main applications of ChIP-Seq are (1) detecting locations 
#### of protein binding (finding read count peaks that are statistically significant) and 
#### (2) finding binding sites for which protein binding frequency differs between 
#### comparison groups.
####

## We will use the chipseq package for our analysis and the TxDb.Mmusculus.UCSC.mm9.knownGene
## package for annotating our results.
biocLite(c("chipseq", "TxDb.Mmusculus.UCSC.mm9.knownGene"))
library(chipseq)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)

## We will use the example dataset cstest from the chipseq package. We have two samples 
## (the two samples were actually run on two lanes of the same device). One sample is for 
## a mouse to which the transcription factor CTCF has been bound. The other sample is for 
## a mouse to which the protein for the reporter gene GFP has been bound. The data we 
## have are a subset of the above, restricted to chromosomes 10, 11, and 12.
##
## The cstest object is of class GRangesList, which is basically of list of GRanges 
## objects. The GRanges class is for storing genomic ranges, including chromosome numbers 
## and strand directions. In our example, we have one range for each aligned read. The 
## actual read sequences have been removed for efficiency purposes. Once we know which 
## genome we are looking at and a particular range, we know the sequence.
data(cstest)
cstest
names(cstest)
elementLengths(cstest)

## While the DNA fragments that are sequenced after immunoprecipitation may be 100's of 
## bp long, NGS reads are typically shorter. After sequencing and alignment to a 
## reference genome, we have read counts associated with specific genomic ranges. There 
## is interest in extending these ranges, so that we are more likely to capture the 
## actual binding sites. The estimate.mean.fraglen function uses one of a handful of 
## available algorithms for estimating the average DNA fragment length associated with 
## our reads. The average fragment length for the CTCF mouse is about 190 nucleotides 
## (nt), while that for the GFP mouse is about 250 nt. We will extend all CTCF reads to 
## be 200 nt long and all the GFP reads to be 250 nt long. The resize function is 
## designed to do this.
estimate.mean.fraglen(cstest$ctcf)
estimate.mean.fraglen(cstest$gfp)

ctcf_ext <- resize(cstest$ctcf, width = 200)
gfp_ext <- resize(cstest$gfp, width = 250)

## The coverage function computes the number of times each base along the genome is 
## covered by a read interval. For example, in cov_ctcf, we are told that on chromosome 
## 12, no reads aligned to any of the first 3,018,214 bases, then one read (of extended 
## length 200 nt) aligned, followed by 24,958 bases with no reads, etc. Wherever reads 
## overlapped, we have counts > 1.
cov_ctcf <- coverage(ctcf_ext)
cov_gfp <- coverage(gfp_ext)

cov_ctcf[[12]]

## The slice function can be used to identify contiguous segments of non-zero coverage. 
## Any binding sites would be expected to lie within such segments. Here, we extract 
## segments with >= 1 overlapping read. For example, there are 87,956 "islands" on 
## chromosome 10 of the CTCF mouse. 
ctcf_islands <- slice(cov_ctcf, lower = 1)
gfp_islands <- slice(cov_gfp, lower = 1)

ctcf_islands[[10]]
length(ctcf_islands[[10]])

## The viewSums function returns the number of reads contained in each island, and the 
## viewMaxs function returns the maximum coverage depth within each island. For example, 
## there are 2,400 / 200 = 12 reads in the first island on chromosome 10 of the CTCF 
## mouse, and the maximum coverage depth in that island is 11.
viewSums(ctcf_islands)
viewMaxs(ctcf_islands)
viewSums(gfp_islands)
viewMaxs(gfp_islands)

## On chromosome 10, there are 95 islands with 12 reads contained in them. There are 99 
## islands with max read depth of 11.
ctcf_nread_tab <- table(viewSums(ctcf_islands) / 200)
ctcf_depth_tab <- table(viewMaxs(ctcf_islands))
head(ctcf_nread_tab, 10)
head(ctcf_depth_tab, 10)

## Create a plot that shows the number of intervals by read depth. The x-axis is depth, 
## and the y-axis is log of interval counts. There are peaks on the left side of the 
## plots, with the counts leveling off on the right, corresponding to greater read depths.
islandDepthPlot(cov_ctcf)
islandDepthPlot(cov_gfp)

## The chipseq package allows for a statistical test of whether a peak is "significant." 
## This is done by assuming a Poisson distribution on the background (not significant) 
## island depths. The peakCutoff function returns an estimated depth cutoff required for 
## a desired FDR. We'll use cutoffs of 6 and 8 for CTCF and GFP, respectively.
peakCutoff(cov_ctcf, fdr = 0.01)
peakCutoff(cov_gfp, fdr = 0.01)

## Extract the islands with depths at or above our cutoffs from above. Shown, for example, 
## is the result for chromosome 10.
peaks_ctcf <- slice(cov_ctcf, lower = 6)
peaks_gfp <- slice(cov_gfp, lower = 8)

peaks_ctcf[[10]]

## The coverageplot function shows the peaks for each island. Below, we show the most 
## significant peak on chromosome 10 for the CTCF mouse. You can also display a figure 
## like this that differentiates between the strands to which each read aligned; see the 
## chipseq vignette for an example.
depths_ctcf <- viewMaxs(peaks_ctcf)
views_ctcf <- Views(cov_ctcf, ranges(peaks_ctcf))
order_ctcf <- order(depths_ctcf$chr10, decreasing = TRUE)

coverageplot(views_ctcf$chr10[order_ctcf[1]])

## We can also test for differential peaks between the two mice using diffPeakSummary. 
## A peak is called differential by chipseq by a simple (non-statistical) rule; see the 
## chipseq vignette for details. Because a formal statistical method is not used, no p-
## values or FDRs are returned.

## The diffPeakSummary computes summary statistics for comparing peaks at each location 
## for the two mice. We get read counts and maximum read depth for each location, for 
## each mouse. We can then employ a simple cutoff rule for extracting "differential" 
## peaks, peaks that are substantially different. Our cutoff rule does not employ a 
## formal statistical test, so we do not get p-values or FDRs for the differential peaks.
diff_peaks <- diffPeakSummary(peaks_ctcf, peaks_gfp)
head(data.frame(diff_peaks))

diff_peaks <- within(diff_peaks, 
  {
    diffs <- sums2 - sums1
    resids <- (diffs - median(diffs)) / mad(diffs)
    up <- resids > 2
    down <- resids < -2
    change <- ifelse(up, "up", ifelse(down, "down", "flat"))
  }
)
diff_peaks <- diff_peaks[diff_peaks$change != "flat", ]

## Now that we have some peaks of interest, we can check for whether they fall within 
## "promoter" regions of known genes in the mouse genome. A promoter region is a region 
## of the DNA sequence that is "upstream" of a gene. Transcription factors typically bind 
## upstream, so that they can regulate the expression of their target gene. "Upstream" 
## would mean that the peak is to the left of the gene on the positive strand and to the 
## right of the gene on the negative strand. Here, we obtain the location information for 
## all the known transcripts in the mouse genome and create candidate promoter regions by 
## extending 1000 bases in either direction of the starting location for each transcript.
gregions <- transcripts(TxDb.Mmusculus.UCSC.mm9.knownGene)
promoters <- flank(gregions, 1000, both = TRUE)

## Now count the number of peaks that fall within the promoter regions. There are 48 
## differential peaks (present in CTCF at greater abundance than GTF) that overlap with 
## promoter regions.
diff_peaks$inPromoter <- diff_peaks %over% promoters
xtabs(~ inPromoter + change, diff_peaks)

## We can view our peak summaries in a genome browser like that provided by UCSC. For 
## this, we need the rtracklayer package, which was installed by one of our already-used 
## packages. The browserView function then opens a UCSC link that displays the results, 
## together with annotation for each DNA region.
library(rtracklayer)

session <- browserSession()
genome(session) <- "mm9"
session$gfpCov <- cov_gfp
session$gfpPeaks <- peaks_gfp
session$ctcfCov <- cov_ctcf
session$ctcfPeaks <- peaks_ctcf

peaks_ctcf_sort <- as(peaks_ctcf, "GRanges")[order_ctcf]
view <- browserView(session, peaks_ctcf_sort[1], full = c("gfpCov", "ctcfCov"))

####
#### Aside: Here is a picture of what overdispersion looks like for the Poisson dist'n.
####

## The rnbinom function can work with either of the two parameterizations we've seen for 
## the negative binomial model. By default, it uses the r, p parameterization, but we can 
## access the specification used by edgeR (in terms of the mean mu and a dispersion 
## parameter phi) via the 'mu' and 'size' arguments.
n <- 1000

mu <- 10
phi <- 0.5

Y_P <- rpois(n, mu)
Y_NB <- rnbinom(n, mu = mu, size = phi)

mean(Y_P); var(Y_P)
mean(Y_NB); var(Y_NB)

hist(Y_NB, col = "pink", xlim = c(0, 50), ylim = c(0, 0.15), prob = TRUE, main = "")
hist(Y_P, col = "lightblue", prob = TRUE, add = TRUE)
legend(30, 0.12, legend = c("Poisson", "Negative Binomial"), lwd = rep(5, 2), 
  col = c("lightblue", "pink"), bty = "n")
