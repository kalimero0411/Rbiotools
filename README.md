# Rbiotools
R scripts for bioinformatics

```
RNAseq

	--rdata	RData file path
	--wd	Working directory path
	--name	Experiment name (all data will output to a directory by that name in the working directory)
	--input	Aligner output files path (RSEM,STAR or kallisto)
	--design	DESeq2 design formula (e.g. ~Genotype+treatment)
	--reduced	DESeq2 reduced design formula (Default = ~1)
	--process	Process: 
		1 = Run settings (interactive only)
		2 = DESeq2
		3 = LRT
		4 = LRT-DESeq
		5 = Transformation (not parallelized)
		6 = Heatmaps
		7 = Dispersion estimates	 PCAs and PCoAs
		8 = Optimize K-means
		9 = K-means clustering
		10 = Hierarchial clustering
		11 = MA-plot
		12 = DEG heatmap
		13 = goseq GO analysis
		14 = topGO analysis
		15 = Word cloud
		16 = Venn diagram
		17 = Variable heatmap and report
		(e.g. 2,5,6,7,8,9,10,11,12)
	--remove_isoforms	Remove isoform suffix from gene IDs
	--k	Number of clusters for K-means clustering
	--seed	Seed value for random number generator
	--heatmap_no_clust	Cluster rows in heatmaps
	--GO_file	Path to GO annotation file
	--t	Number of compute threads
	--arg	Additional R arguments (multiple arguments in separate flags)
```
	
- The script creates a directory with the experiment name in the provided working directory.
- Create 3D PCA videos with 3D\_PCA_run (shell).


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
