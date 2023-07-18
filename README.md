# Rbiotools
R scripts for bioinformatics

```
RNAseq
	--rdata	RData file path
        --wd	Working directory path
        --name	Experiment name
        --input	Aligner output files path (RSEM, STAR or kallisto)
        --mapper	One of RSEM, Kallisto, Salmon, HTseq-count, Counts
        --design	DESeq2 design formula (e.g. ~Genotype+treatment)
        --exp	Experimental design file (Control first)
        --reduced	DESeq2 reduced design formula (Default = ~1)
        --process	Process:
          1 = Run settings
          2 = DESeq2
          3 = LRT
          4 = Transformation (not parallelized)
          5 = Heatmaps
          6 = Dispersion estimates, PCAs and PCoAs
          7 = Optimize K-means
          8 = K-means clustering
          9 = Hierarchial clustering
          10 = MA-plot
          11 = DEG heatmap
          12 = goseq GO analysis
          13 = Venn diagram
          14 = topGO analysis
          15 = Word cloud
          16 = Variable heatmap and report
          (e.g. 2,5,6,7,8,9,10,11,12)
        --remove_isoforms	Remove isoform suffix from gene IDs
        --isoforms	Use isoforms instead of genes
        --vst	Use VST instead of rlog
        --lengths	File with extention .gtf, .rds, .RData, or a tab delimited file
        --tx2gene	A transcript to gene mapping file (Kallisto/Salmon)
        --alpha	Cut off for p values
        --FDR	Cut off for FDR
        --k	Number of clusters for K-means clustering
        --seed	Seed value for random number generator
        --heatmap_no_clust	Cluster rows in heatmaps
        --GO_file	Path to GO annotation file
        --noiseq	Perform NOISeq correction (ARSyNseq: counts | proportion | FDR)
        --venn_go	0 = None, 1 = Up-regulated, 2 = Down-regulated, 3 = Both (Default = 0; Multiple input possible)
        --3dpca	1 = Single factors, 2 = Multiple factors, 3 = DEGs (Default = None; All = 123)
        --ensembl	Select Ensembl species ID (enter 0 for options)
        --t	Number of compute threads
        --arg	Additional R arguments (multiple arguments in separate flags)
```

```
ChIPSeq

	--rdata	RData file path
	--wd	Working directory path
	--name	Experiment name (all data will output to a directory by that name in the working directory)
	--annotation	Annotation file path (GFF or GTF file)
	--xls	MACS2 Excel output path (if applicable)
	--paired	Designates that the data is paired-end (if applicable); default = FALSE (single-end)
	--t	Number of compute threads
	--arg	Additional R arguments (multiple arguments in separate flags)
```

- First set-up in interactive R session, then run in shell (optional).

```
MethylSeq

	--rdata	RData file path
	--wd	Working directory path
	--name	Experiment name (all data will output to a directory by that name in the working directory)
	--bismark	Bismark output BAM files path
	--import	Import BAM files (Y), do not import BAM files (N) or only import BAM files (O); default = N
	--context	Methylation context: CpG, CHG or CHH
	--savedb	Data will be output into flat files (for reduced RAM usage); default = FALSE
	--t	Number of compute threads
	--arg	Additional R arguments (multiple arguments in separate flags)
```

- First set-up in interactive R session, then run in shell (optional).

```
3D_PCA_run

	-d | --dir      DESeq2 output directory
        -f | --factors  Factors for 3D PCA (e.g. factor1,factor2,factor3; Default = All factors)
        -l | --list     List all available factors for a specific analysis
        -p | --pca      Types of PCA to return (1 = Single factor, 2 = Multi-factor, 3 = DEGs; e.g. 13; Default = 123)
        -r | --rscript  Path to PCA_3D.R (Defaults bash script directory)
        -h | --help     Display help
```

- The DESeq2 directory should include PCA_data.RData, created by RNAseq.R
