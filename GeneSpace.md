Written by Li Lei/2022/10/28

This tutorial is for Miguel Campos

Here are the resources about GENESPACE:

- Publication: 
[GENESPACE tracks regions of interest and gene copy number variation across multiple genomes](https://elifesciences.org/articles/78526)
- Github: https://github.com/jtlovell/GENESPACE 

1) Install the GENESPACE
Install R from the recent release in [Can](https://www.r-project.org/) >4.0

2) Install orthofinder

```
conda create -n orthofinder
conda activate orthofinder
conda install -c bioconda orthofinder 
```
 Since if intall the orthofinder, the conda will directly install diamond0.9 version.
So you need to install diamond=2.0.15 via below command line

``` 
conda config --add channels conda-forge
conda install -c bioconda diamond=2.0.15
```
3) Install MCScanX should be installed from [github](https://github.com/wyp1125/MCScanX)

```
git clone https://github.com/wyp1125/MCScanX.git

cd MCScanx
make
```
4) Open an interactive R session
```
open -na rstudio
```
5) Install GeneSpace in Rstudio
```
devtools::install_github("jtlovell/GENESPACE@v0.9.3")
library(GENESPACE)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))

```
6) Run GENESPACE

```
library(GENESPACE)

wd <- "~/Desktop/genespace_runs/liTroubleshoot3June2022"

rawAnnotationDir <- "~/Desktop/genespace_runs/rawGenomeRepo/"

mcscanDir <- "/Users/lovell/Desktop/programs/MCScanX/"

nThreads <- 4

path2of <- "orthofinder" # since orthofinder is in the path via conda

gi <- list.files(file.path(rawAnnotationDir, "LiBrachy"))
gpar <- init_genespace(
  wd = wd,
  nCores = nThreads, 
  orthofinderInBlk = F, 
  speciesIDs = rep("LiBrachy", length(gi)),
  genomeIDs = gi,
  versionIDs = gi,
  path2orthofinder = path2of, 
  path2mcscanx = mcscanDir, 
  ploidy = rep(1, length(gi)),
  rawGenomeDir = rawAnnotationDir)

parse_phytozome(gpar) #This is very important and with this command line, I've fixed the bus error!!!!!

gpar <- run_orthofinder(gpar)

gpar <- synteny(gpar)

gpar <- run_orthofinder(gpar)

```
Then it should be have the below things coming out

```
Synteny Parameters have not been set! Setting to defaults
	Running 'defualt' genespace orthofinder method 
	############################################################
	Cleaning out orthofinder directory and prepping run
	Calculating blast results and running OrthoFinder 
	################################################## 
	##################################################
	
	
	OrthoFinder version 2.5.4 Copyright (C) 2014 David Emms
	
	2022-06-03 10:01:46 : Starting OrthoFinder 2.5.4
	4 thread(s) for highly parallel tasks (BLAST searches etc.)
	1 thread(s) for OrthoFinder algorithm
	
	Checking required programs are installed
	----------------------------------------
	Test can run "mcl -h" - ok
	Test can run "fastme -i /Users/lovell/Desktop/genespace_runs/liTroubleshoot3June2022/orthofinder/Results_Jun03/WorkingDirectory/SimpleTest.phy -o /Users/lovell/Desktop/genespace_runs/liTroubleshoot3June2022/orthofinder/Results_Jun03/WorkingDirectory/SimpleTest.tre" - ok
	
	Dividing up work for BLAST for parallel processing
	--------------------------------------------------
	2022-06-03 10:01:48 : Creating diamond database 1 of 11
	2022-06-03 10:01:48 : Creating diamond database 2 of 11
	2022-06-03 10:01:48 : Creating diamond database 3 of 11
	2022-06-03 10:01:48 : Creating diamond database 4 of 11
	2022-06-03 10:01:48 : Creating diamond database 5 of 11
	2022-06-03 10:01:49 : Creating diamond database 6 of 11
	2022-06-03 10:01:49 : Creating diamond database 7 of 11
	2022-06-03 10:01:49 : Creating diamond database 8 of 11
	2022-06-03 10:01:49 : Creating diamond database 9 of 11
	2022-06-03 10:01:49 : Creating diamond database 10 of 11
	2022-06-03 10:01:49 : Creating diamond database 11 of 11
	
	Running diamond all-versus-all
	------------------------------
	Using 4 thread(s)
	2022-06-03 10:01:49 : This may take some time....
	2022-06-03 10:01:50 : Done 0 of 121
	2022-06-03 10:17:36 : Done 10 of 121
	2022-06-03 10:29:55 : Done 20 of 121
	2022-06-03 10:43:57 : Done 30 of 121
	2022-06-03 10:57:10 : Done 40 of 121
	2022-06-03 11:09:41 : Done 50 of 121
	2022-06-03 11:24:56 : Done 60 of 121
	2022-06-03 11:37:12 : Done 70 of 121
	2022-06-03 11:52:14 : Done 80 of 121
	2022-06-03 12:04:00 : Done 90 of 121
	2022-06-03 12:18:15 : Done 100 of 121
	2022-06-03 12:29:31 : Done 110 of 121
	2022-06-03 12:44:51 : Done all-versus-all sequence search
	
	Running OrthoFinder algorithm
	-----------------------------
	2022-06-03 12:44:52 : Initial processing of each species
	2022-06-03 12:45:31 : Initial processing of species 0 complete
	2022-06-03 12:46:06 : Initial processing of species 1 complete
	2022-06-03 12:46:45 : Initial processing of species 2 complete
	2022-06-03 12:47:20 : Initial processing of species 3 complete
	2022-06-03 12:47:51 : Initial processing of species 4 complete
	2022-06-03 12:48:20 : Initial processing of species 5 complete
	2022-06-03 12:48:55 : Initial processing of species 6 complete
	2022-06-03 12:49:31 : Initial processing of species 7 complete
	2022-06-03 12:50:05 : Initial processing of species 8 complete
	2022-06-03 12:50:41 : Initial processing of species 9 complete
	2022-06-03 12:51:19 : Initial processing of species 10 complete
	2022-06-03 12:51:49 : Connected putative homologues
	2022-06-03 12:51:55 : Written final scores for species 0 to graph file
	2022-06-03 12:52:00 : Written final scores for species 1 to graph file
	2022-06-03 12:52:06 : Written final scores for species 2 to graph file
	2022-06-03 12:52:11 : Written final scores for species 3 to graph file
	2022-06-03 12:52:15 : Written final scores for species 4 to graph file
	2022-06-03 12:52:20 : Written final scores for species 5 to graph file
	2022-06-03 12:52:25 : Written final scores for species 6 to graph file
	2022-06-03 12:52:31 : Written final scores for species 7 to graph file
	2022-06-03 12:52:36 : Written final scores for species 8 to graph file
	2022-06-03 12:52:41 : Written final scores for species 9 to graph file
	2022-06-03 12:52:46 : Written final scores for species 10 to graph file
	2022-06-03 12:53:49 : Ran MCL
	
	Writing orthogroups to file
	---------------------------
	OrthoFinder assigned 332270 genes (93.1% of total) to 34634 orthogroups. Fifty percent of all genes were in orthogroups with 11 or more genes (G50 was 11) and were contained in the largest 11061 orthogroups (O50 was 11061). There were 14687 orthogroups with all species present and 8890 of these consisted entirely of single-copy genes.
	
	2022-06-03 12:54:03 : Done orthogroups
	
	Analysing Orthogroups
	=====================
	
	Calculating gene distances
	--------------------------
	2022-06-03 12:59:37 : Done
	2022-06-03 12:59:39 : Done 0 of 25720
	2022-06-03 12:59:45 : Done 1000 of 25720
	2022-06-03 12:59:50 : Done 2000 of 25720
	2022-06-03 12:59:56 : Done 3000 of 25720
	2022-06-03 13:00:01 : Done 4000 of 25720
	2022-06-03 13:00:06 : Done 5000 of 25720
	2022-06-03 13:00:12 : Done 6000 of 25720
	2022-06-03 13:00:17 : Done 7000 of 25720
	2022-06-03 13:00:23 : Done 8000 of 25720
	2022-06-03 13:00:28 : Done 9000 of 25720
	2022-06-03 13:00:34 : Done 10000 of 25720
	2022-06-03 13:00:39 : Done 11000 of 25720
	2022-06-03 13:00:44 : Done 12000 of 25720
	2022-06-03 13:00:50 : Done 13000 of 25720
	2022-06-03 13:00:55 : Done 14000 of 25720
	2022-06-03 13:01:01 : Done 15000 of 25720
	2022-06-03 13:01:06 : Done 16000 of 25720
	2022-06-03 13:01:12 : Done 17000 of 25720
	2022-06-03 13:01:17 : Done 18000 of 25720
	2022-06-03 13:01:22 : Done 19000 of 25720
	2022-06-03 13:01:28 : Done 20000 of 25720
	2022-06-03 13:01:33 : Done 21000 of 25720
	2022-06-03 13:01:38 : Done 22000 of 25720
	2022-06-03 13:01:44 : Done 23000 of 25720
	2022-06-03 13:01:49 : Done 24000 of 25720
	2022-06-03 13:01:55 : Done 25000 of 25720
	
	Inferring gene and species trees
	--------------------------------
	
	14687 trees had all species present and will be used by STAG to infer the species tree
	
	Best outgroup(s) for species tree
	---------------------------------
	2022-06-03 13:07:16 : Starting STRIDE
	2022-06-03 13:07:21 : Done STRIDE
	Observed 289 well-supported, non-terminal duplications. 286 support the best root and 3 contradict it.
	Best outgroup for species tree:
	  BmexP, BmexU
	
	Reconciling gene trees and species tree
	---------------------------------------
	Outgroup: BmexP, BmexU
	2022-06-03 13:07:21 : Starting Recon and orthologues
	2022-06-03 13:07:21 : Starting OF Orthologues
	2022-06-03 13:07:22 : Done 0 of 25720
	2022-06-03 13:07:32 : Done 1000 of 25720
	2022-06-03 13:07:38 : Done 2000 of 25720
	2022-06-03 13:07:43 : Done 3000 of 25720
	2022-06-03 13:07:47 : Done 4000 of 25720
	2022-06-03 13:07:50 : Done 5000 of 25720
	2022-06-03 13:07:54 : Done 6000 of 25720
	2022-06-03 13:07:57 : Done 7000 of 25720
	2022-06-03 13:08:00 : Done 8000 of 25720
	2022-06-03 13:08:04 : Done 9000 of 25720
	2022-06-03 13:08:07 : Done 10000 of 25720
	2022-06-03 13:08:11 : Done 11000 of 25720
	2022-06-03 13:08:15 : Done 12000 of 25720
	2022-06-03 13:08:19 : Done 13000 of 25720
	2022-06-03 13:08:24 : Done 14000 of 25720
	2022-06-03 13:08:29 : Done 15000 of 25720
	2022-06-03 13:08:33 : Done 16000 of 25720
	2022-06-03 13:08:36 : Done 17000 of 25720
	2022-06-03 13:08:40 : Done 18000 of 25720
	2022-06-03 13:08:43 : Done 19000 of 25720
	2022-06-03 13:08:46 : Done 20000 of 25720
	2022-06-03 13:08:49 : Done 21000 of 25720
	2022-06-03 13:08:51 : Done 22000 of 25720
	2022-06-03 13:08:53 : Done 23000 of 25720
	2022-06-03 13:08:54 : Done 24000 of 25720
	2022-06-03 13:08:55 : Done 25000 of 25720
	2022-06-03 13:08:56 : Done OF Orthologues
	
	Writing results files
	=====================
	2022-06-03 13:09:04 : Done orthologues
	
	Results:
	    /Users/lovell/Desktop/genespace_runs/liTroubleshoot3June2022/orthofinder/Results_Jun03/
	
	CITATION:
	 When publishing work that uses OrthoFinder please cite:
	 Emms D.M. & Kelly S. (2019), Genome Biology 20:238
	
	 If you use the species tree in your work then please also cite:
	 Emms D.M. & Kelly S. (2017), MBE 34(12): 3267-3278
	 Emms D.M. & Kelly S. (2018), bioRxiv https://doi.org/10.1101/267914

```
Then run below commandline in Rstudio

```
gpar <- synteny(gsParam = gpar)

ripdat <- plot_riparianHits(gpar)
pg <- pangenome(gpar)
```
It should give your the plot
