pkgname <- "SoupX"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SoupX')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("SoupChannel")
### * SoupChannel

flush(stderr()); flush(stdout())

### Name: SoupChannel
### Title: Construct a SoupChannel object
### Aliases: SoupChannel

### ** Examples

#Load droplet and count tables
if(requireNamespace("Seurat", quietly = TRUE)) {
  tod = Seurat::Read10X(system.file('extdata','toyData','raw_gene_bc_matrices','GRCh38',
                                    package='SoupX'))
  toc = Seurat::Read10X(system.file('extdata','toyData','filtered_gene_bc_matrices','GRCh38',
                                    package='SoupX'))
  #Default calculates soup profile
  sc = SoupChannel(tod,toc)
  names(sc)
  #This can be suppressed
  sc = SoupChannel(tod,toc,calcSoupProfile=FALSE)
  names(sc)
}



cleanEx()
nameEx("adjustCounts")
### * adjustCounts

flush(stderr()); flush(stdout())

### Name: adjustCounts
### Title: Remove background contamination from count matrix
### Aliases: adjustCounts

### ** Examples

# Ensure contamination fraction is calculated first
sc <- scToy
# Set contamination fraction manually for the toy dataset
sc <- setContaminationFraction(sc, 0.1)
out = adjustCounts(sc)
#Return integer counts only
out = adjustCounts(sc,roundToInt=TRUE)



cleanEx()
nameEx("autoEstCont")
### * autoEstCont

flush(stderr()); flush(stdout())

### Name: autoEstCont
### Title: Automatically calculate the contamination fraction
### Aliases: autoEstCont

### ** Examples

#Use less specific markers
scToy = autoEstCont(scToy,tfidfMin=0.8)
#Allow large contamination fractions to be allocated
scToy = autoEstCont(scToy,forceAccept=TRUE)
#Be quiet
scToy = autoEstCont(scToy,verbose=FALSE,doPlot=FALSE)



cleanEx()
nameEx("benchmark_soupx")
### * benchmark_soupx

flush(stderr()); flush(stdout())

### Name: benchmark_soupx
### Title: Benchmark SoupX Performance
### Aliases: benchmark_soupx

### ** Examples

# Benchmark core SoupX operations
if(requireNamespace("Matrix", quietly = TRUE)) {
  # Create a small test dataset
  toc <- Matrix::Matrix(rpois(1000, 5), nrow = 100, ncol = 10, sparse = TRUE)
  tod <- Matrix::Matrix(rpois(2000, 3), nrow = 100, ncol = 20, sparse = TRUE)
  rownames(toc) <- rownames(tod) <- paste0("Gene", 1:100)
  colnames(toc) <- paste0("Cell", 1:10)
  colnames(tod) <- paste0("Droplet", 1:20)
  
  sc <- SoupChannel(tod, toc)
  results <- benchmark_soupx(sc, iterations = 2, verbose = FALSE)
  print(results)
}



cleanEx()
nameEx("calculateContaminationFraction")
### * calculateContaminationFraction

flush(stderr()); flush(stdout())

### Name: calculateContaminationFraction
### Title: Calculate the contamination fraction
### Aliases: calculateContaminationFraction

### ** Examples

#Common gene list in real world data
geneList = list(HB=c('HBB','HBA2'))
#Gene list appropriate to toy data
geneList = list(CD7 = 'CD7')
ute = estimateNonExpressingCells(scToy,geneList)
sc = calculateContaminationFraction(scToy,geneList,ute)



cleanEx()
nameEx("estimateNonExpressingCells")
### * estimateNonExpressingCells

flush(stderr()); flush(stdout())

### Name: estimateNonExpressingCells
### Title: Calculate which cells genuinely do not express a particular gene
###   or set of genes
### Aliases: estimateNonExpressingCells

### ** Examples

#Common gene list in real world data
geneList = list(HB=c('HBB','HBA2'))
#Gene list appropriate to toy data
geneList = list(CD7 = 'CD7')
ute = estimateNonExpressingCells(scToy,geneList)



cleanEx()
nameEx("estimateSoup")
### * estimateSoup

flush(stderr()); flush(stdout())

### Name: estimateSoup
### Title: Get expression profile of soup
### Aliases: estimateSoup

### ** Examples

#Load droplet and count tables
tod = Seurat::Read10X(system.file('extdata','toyData','raw_gene_bc_matrices','GRCh38',
                                  package='SoupX'))
toc = Seurat::Read10X(system.file('extdata','toyData','filtered_gene_bc_matrices','GRCh38',
                                  package='SoupX'))
#Suppress calculation of soup profile automatically on load
sc = SoupChannel(tod,toc,calcSoupProfile=FALSE)
#Retain table of droplets
sc = estimateSoup(sc,keepDroplets=TRUE)
#Or use non-default values
sc = estimateSoup(sc,soupRange=c(60,100))



cleanEx()
nameEx("generateSoupXReport")
### * generateSoupXReport

flush(stderr()); flush(stdout())

### Name: generateSoupXReport
### Title: Generate Comprehensive SoupX Report
### Aliases: generateSoupXReport

### ** Examples

# Generate a comprehensive report
if(requireNamespace("ggplot2", quietly = TRUE)) {
  # Use the toy dataset
  data(scToy)
  
  # Generate report
  report_path <- generateSoupXReport(scToy, format = "html")
  cat("Report generated at:", report_path, "\n")
}



cleanEx()
nameEx("load10X")
### * load10X

flush(stderr()); flush(stdout())

### Name: load10X
### Title: Load a collection of 10X data-sets
### Aliases: load10X

### ** Examples

sc = load10X(system.file('extdata','toyData',package='SoupX'))



cleanEx()
nameEx("plotChangeMap")
### * plotChangeMap

flush(stderr()); flush(stdout())

### Name: plotChangeMap
### Title: Plot maps comparing corrected/raw expression
### Aliases: plotChangeMap

### ** Examples

out = adjustCounts(scToy)
gg = plotChangeMap(scToy,out,'S100A9')



cleanEx()
nameEx("plotContaminationByCluster")
### * plotContaminationByCluster

flush(stderr()); flush(stdout())

### Name: plotContaminationByCluster
### Title: Plot Contamination Fraction by Cluster
### Aliases: plotContaminationByCluster

### ** Examples

# Plot contamination by cluster
if(requireNamespace("ggplot2", quietly = TRUE)) {
  # Use the toy dataset
  data(scToy)
  
  # Set some dummy clusters
  scToy$metaData$clusters <- rep(c("A", "B"), each = nrow(scToy$metaData)/2)
  
  # Plot contamination by cluster
  gg <- plotContaminationByCluster(scToy)
  print(gg)
}



cleanEx()
nameEx("plotGeneContamination")
### * plotGeneContamination

flush(stderr()); flush(stdout())

### Name: plotGeneContamination
### Title: Plot Gene-Specific Contamination Analysis
### Aliases: plotGeneContamination

### ** Examples

# Plot gene contamination analysis
if(requireNamespace("ggplot2", quietly = TRUE)) {
  # Use the toy dataset
  data(scToy)
  
  # Analyze specific genes
  gg <- plotGeneContamination(scToy, c("CD7", "LTB", "S100A9"))
  print(gg)
}



cleanEx()
nameEx("plotMarkerDistribution")
### * plotMarkerDistribution

flush(stderr()); flush(stdout())

### Name: plotMarkerDistribution
### Title: Plots the distribution of the observed to expected expression
###   for marker genes
### Aliases: plotMarkerDistribution

### ** Examples

gg = plotMarkerDistribution(scToy,list(CD7='CD7',LTB='LTB'))



cleanEx()
nameEx("plotMarkerMap")
### * plotMarkerMap

flush(stderr()); flush(stdout())

### Name: plotMarkerMap
### Title: Plot ratio of observed to expected counts on reduced dimension
###   map
### Aliases: plotMarkerMap

### ** Examples

gg = plotMarkerMap(scToy,'CD7')



cleanEx()
nameEx("plotQualityControl")
### * plotQualityControl

flush(stderr()); flush(stdout())

### Name: plotQualityControl
### Title: Generate Quality Control Plots for SoupX Analysis
### Aliases: plotQualityControl

### ** Examples

# Generate QC plots
if(requireNamespace("ggplot2", quietly = TRUE)) {
  # Use the toy dataset
  data(scToy)
  
  # Generate QC plots
  qc_plots <- plotQualityControl(scToy)
  
  # Access individual plots
  print(qc_plots$soup_profile)
  print(qc_plots$contamination_distribution)
}



cleanEx()
nameEx("quickMarkers")
### * quickMarkers

flush(stderr()); flush(stdout())

### Name: quickMarkers
### Title: Gets top N markers for each cluster
### Aliases: quickMarkers

### ** Examples

#Calculate markers of clusters in toy data
mrks = quickMarkers(scToy$toc,scToy$metaData$clusters)
## Not run: 
##D #Calculate markers from Seurat (v3) object
##D mrks = quickMarkers(srat@assays$RNA@count,srat@active.ident)
## End(Not run)



cleanEx()
nameEx("setClusters")
### * setClusters

flush(stderr()); flush(stdout())

### Name: setClusters
### Title: Sets clustering for SoupChannel
### Aliases: setClusters

### ** Examples

sc = load10X(system.file('extdata','toyData',package='SoupX'))
mDat = read.table(system.file('extdata','toyData','metaData.tsv',package='SoupX'),sep='\t')
sc = setClusters(sc,mDat$res.1)



cleanEx()
nameEx("setContaminationFraction")
### * setContaminationFraction

flush(stderr()); flush(stdout())

### Name: setContaminationFraction
### Title: Manually set contamination fraction
### Aliases: setContaminationFraction

### ** Examples

sc = load10X(system.file('extdata','toyData',package='SoupX'))
sc = setContaminationFraction(sc,0.1)



cleanEx()
nameEx("setDR")
### * setDR

flush(stderr()); flush(stdout())

### Name: setDR
### Title: Manually set dimension reduction for a channel
### Aliases: setDR

### ** Examples

sc = load10X(system.file('extdata','toyData',package='SoupX'))
mDat = read.table(system.file('extdata','toyData','metaData.tsv',package='SoupX'),sep='\t')
sc = setDR(sc,mDat[,c('tSNE_1','tSNE_2')])



cleanEx()
nameEx("setSoupProfile")
### * setSoupProfile

flush(stderr()); flush(stdout())

### Name: setSoupProfile
### Title: Set soup profile
### Aliases: setSoupProfile

### ** Examples

#Suppose only table of counts is available
toc = Seurat::Read10X(system.file('extdata','toyData','filtered_gene_bc_matrices','GRCh38',
                                  package='SoupX'))
#Suppress calculating soup profile automatically
sc = SoupChannel(toc,toc,calcSoupProfile=FALSE)
#And add manually
rowSums = Matrix::rowSums
soupProf = data.frame(row.names = rownames(toc),est=rowSums(toc)/sum(toc),counts=rowSums(toc))
sc = setSoupProfile(sc,soupProf)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
