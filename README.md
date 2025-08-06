# Rbiotools
R scripts for bioinformatics


## RNAseq
### Analyze RNA sequencing data
```
	--rdata	RData file path
        --wd	Working directory path
        --name	Experiment name
        --input	Aligner output files path (RSEM, STAR or kallisto)
        --mapper	One of RSEM, Kallisto, Salmon, HTseq-count, FeatureCounts, Counts
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


## ChIPSeq
### Analyze MACS3 output of ChIP sequencing
```
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


## MethylSeq
### Analyze Bismark output of Bisulfite sequencing
```
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

## WGCNA
### Perform WGCNA clustering on RNAseq data
```
                           --rdata	RData file path
                           --wd	Working directory path
                           --name	Experiment name
                           --factors	Experimantal factors to test(Seperated by comma; Default = All factors)
                           --power	Scale Free Topology Model Fit
                           --TOM	Specify sft value and create Topology overlap matrix (TOM)
                           --blockwise	Perform blockwise calculation for TOM
                           --minmod	Minimum module size
                           --maxblock	Maximum block size (only used in blockwise)
                           --min_exp	Minimum expression for all samples of a gene
                           --threads	Number of CPU threads
                           --arg	Additional R arguments (multiple arguments in separate flags)
```
                           

## dups
### Create duplcate plot from FASTQ alignment
```
	--bam	Path to duplicate marked / unmarked BAM file
	--gtf	GTF file used in the alignment
	--stranded	Strandedness of the FASTQ file (0 = unstranded [default]; 1 = stranded; 2 = reverse)
	--verbose	Verbose
```

## Ks
### Ks ridge plot from MCScanX output
```
                         --gff	Path to folder with gff files
                         --kaks	Path to folder with ka/ks collinearity files
                         --order	Order of species by kaks file names (without extention, comma separated)
                         --color	Species group number (for color designation, e.g. '1,1,2,2,2,3,3')
                         --ks_cutoff	Cutoff for Ks plot (Default = 2)
                         --regex	Regex to rename files to samples
```
