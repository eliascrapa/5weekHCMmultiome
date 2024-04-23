#devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())
#ArchR::installExtraPackages()
suppressPackageStartupMessages(library(ArchR))
library("xlsx")
library(tidyverse)
library(hdf5r)
addArchRGenome("mm10")
set.seed(1)
library(here)
library(ggplot2)
addArchRGenome("mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
#HDF5_USE_FILE_LOCKING=FALSE
#RHDF5_USE_FILE_LOCKING=FALSE
set.seed(1)
library(here)
library(ggplot2)
set.seed(1)
library(hexbin)
library(pheatmap)

# setwd("~/ArchR/All/")
setwd(here::here())


#Get Input Fragment Files
inputFiles <- getInputFiles(c("~/ArchR/MHC/TNTdata","~/ArchR/MHC/data","~/ArchR/CTR-TNT/data","~/ArchR/CTR-MHC/data"))
names(inputFiles) <- c("TNT","MHC","CTR-TNT","CTR-MHC")

#Create ArrowFile
ArrowFiles <- createArrowFiles(
 inputFiles = inputFiles,
 sampleNames = names(inputFiles),
 addTileMat = TRUE,
 addGeneScoreMat = TRUE,
 subThreading = F
)

#ArchRProject #Add scRNA
#Import scRNA
listcs17<-list.files(pattern=".arrow")
proj<-ArchRProject(ArrowFiles=listcs17, outputDirectory="cs17-out", copyArrows=FALSE)
names<-gsub(".arrow", "", listcs17)
files<-paste("~/ArchR/", names, "/filtered_feature_bc_matrix.h5", sep="")
seRNA <- import10xFeatureMatrix(input = files, names = names)
seRNAcombined<-cbind(assay(seRNA[[1]]), assay(seRNA[[2]]), assay(seRNA[[3]]), assay(seRNA[[4]]))
seRNA2<-SummarizedExperiment(assays=list(counts=seRNAcombined), rowRanges= rowRanges(seRNA[[1]]))

#Add scRNA to ArchRProject
proj<-addGeneExpressionMatrix(input=proj, seRNA=seRNA2)

#Filter Cells
proj <-subsetArchRProject(
  ArchRProj = proj,
  cells = getCellNames(proj[which(proj$TSSEnrichment > 5 & proj$nFrags > 1000 & !is.na(proj$Gex_nUMI))]),
  outputDirectory = "New",
  dropCells = FALSE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

#Doublet Filtration. Currently disabled just for tutorial. If you want to filter doublets uncomment below.
proj <- addDoubletScores(proj)
proj <- filterDoublets(proj)

#LSI-ATAC
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.3, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  iterations = 4,
  name = "LSI_ATAC"
)

#LSI-RNA
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)

#Combined Dims
proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

#Harmonize Samples
 proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "LSI_Combined",
  name = "Harmony2",
  groupBy = "Sample"
)




#UMAPs
proj <- addUMAP(proj, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "Harmony", name = "UMAP_Harmony", minDist = 0.8, force = TRUE)

#Add Clusters:
proj <- addClusters(proj, reducedDims = "Harmony", name = "Clusters", resolution = 0.8, force = TRUE)


#Calculation Cells per Cluster and Identity
{
#Show number of cells per cluster
cellspercluster <- as.data.frame(table(proj$Clusters))
colnames(cellspercluster) <- c("Cluster","Cellnumber")
cellspercluster

#Confusion Martix + Plot to see which samples are in which cluster
cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))
cM

#plot this confusionmatrix
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
pClusterHM <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
pClusterHM

plotPDF(pClusterHM, name = "Confusion-Matrix-Cluster-vs-Samples_Heatmap", addDOC = FALSE)
cMdf <- as.data.frame(cM)
cMdf <-  cbind(rownames(cMdf),cMdf)
colnames(cMdf)[1] <- "Cluster"
cMall <- left_join(cMdf,cellspercluster, by = "Cluster")
write.xlsx(cMall, file = paste0("~/ArchR/All/Confusion-Matrix-ClusterperSample.xlsx"))

rm(cMdf)
rm(cellspercluster)
rm(cMall)
rm(pClusterHM)
rm(cM)
}
------
#Save Normalized, PCAed, Harmonyiesed, and UMAP proected and clusted data
saveArchRProject(ArchRProj = proj, outputDirectory = "allout", load = FALSE)


#Plot Embedding
p1 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_ATAC", size = 0.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_RNA", size = 0.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Combined", size = 0.1, labelAsFactors=F, labelMeans=F)
p4 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Harmony", size = 0.05, labelAsFactors=T, labelMeans=T)

#Save Plot
plotPDF(p1, p2, p3, p4, name = "UMAP-scATAC-scRNA-Combined", addDOC = FALSE)


#Adding Columns to identify subsets of data
{
proj <- addCellColData(ArchRProj = proj, data = paste0(proj@cellColData$Clusters,"_x_",proj@cellColData$Sample), name = "ClusterBySample", cells = getCellNames(proj), force = TRUE)

cellgroups <- plyr::mapvalues(
  x = proj@cellColData$Sample, 
  from = c("CTR-MHC", "MHC", "CTR-TNT", "TNT"), 
  to = c("CTR", "MHC", "CTR", "TNT")
)
proj <- addCellColData(ArchRProj = proj, data = paste0(proj@cellColData$Clusters,"_x_",cellgroups), name = "ClusterOneControl", cells = getCellNames(proj), force = TRUE)

rm(cellgroups)

celltypes <- plyr::mapvalues(
  x = proj@cellColData$Clusters, 
  from = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25"), 
  to = c("Leukocyte", "Leukocyte", "Leukocyte", "Leukocyte", "Leukocyte", 
         "Cardiomyocyte","Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", "EpicardialCell",
         "EpicardialCell", "EpicardialCell", "EndothelialCell", "EndothelialCell", "EndothelialCell",
         "EndothelialCell", "EndothelialCell", "VSMC", "VSMC", "VSMC", 
         "Fibroblast", "Fibroblast", "Fibroblast", "undeterm", "undeterm")
)
proj <- addCellColData(ArchRProj = proj, data = paste0(celltypes,"_x_",cellgroups), name = "CelltypesPerGroups", cells = getCellNames(proj), force = TRUE)

rm(cellgroups)

proj <- addCellColData(ArchRProj = proj, data = paste0(proj@cellColData$Sample), name = "Samples", cells = getCellNames(proj), force = TRUE)

}

#Pseudobulking
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "ClusterOneControl",force = FALSE)

getAvailableMatrices(proj)

proj <- addImputeWeights(proj, reducedDims = "Harmony")

#Exporting Embeddings for Seurat
df<- proj@embeddings@listData[["UMAP_Harmony"]]@listData[["df"]]
df$ArchR_Clusters <- proj@cellColData@listData[["Clusters"]]
df$ClusterBySample <- proj@cellColData@listData[["ClusterBySample"]]
rownames(df) <- gsub("#","_",rownames(df))

row.names(df) <- lapply(row.names(df), function(x) paste(x, sep="_", "1"))
df <- cbind(rownames(df), data.frame(df, row.names=NULL))
df$`rownames(df)` <- gsub(".*#","",df$`rownames(df)`)

markerGenes  <- c(
  "Col5a1", "Col3a1", "Col1a1", #Fibroblasts
 "MHRT", "Tnnt2", #Cardimyocyte
 "Nppa", "Myl7", "Sln", "Myh6", #atrial Myocyte
 "PAX5", "MS4A1", "EBF1", "MME", "Cd79a", #B-Cell Trajectory
 "Jchain", #"Igk-C", "Igh-6", #Plasmacyte
 "CD14", "CEBPB", "MPO", #Monocytes
 "C1qa", "Vsig4",#Macrophages
 "S100a9","S100a8", #Neurophils
 "IRF8", "Postn", "TCF21","Mfap5", "Col3a1", "Col5a1","Pdgfra", #Fibroblasts
 "H19", "Igf2", "Vwf", #Endocardial Cell
 "CD3D", "CD8A", "TBX21", "IL7R", "CD8A", "Ccl5", #Adult CD8+ T-Cell
 "Cd209a", "H2-Ab1", #Dendritic  Cells
 "Pecam1", "Itgb1", "CD34", "Lyve1", "Plvap", "Vtn", #Endothelial Cells Lyve1, Plvap and Vtn high clusters described
 "Gpihbp1", "Cd36", #Capillary endothelial cell_Gpihbp1 high-C2
 "Acta2", "Tagln", "Mylk", "Rgs5", #Smooth muscle_cells - RGs5 and ACTA2 high types possible
 "Hba-a1", "Alas2", #Erythroid Cells
 "Retnla", "Lyz1", #Eosinophil
 "Car3", "Cfd", "Adipoq" # Adipose Cell
)

#project markergenes with Genexpression and Genescore on Clusters
{
p5 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneExpressionMatrix", 
    name = markerGenes, 
  continuousSet = "horizonExtra",
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(proj)
)

p5

p6 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  continuousSet = "horizonExtra",
  name = markerGenes, 
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(proj)
)

p6

plotPDF(p5, p6, name = "PLot-Embedding by Geneexpression and Genescore", addDOC = FALSE)



# Geneexpressionmarkers for Each Cluster
{
markersGEC <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneExpressionMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

  markerList <- getMarkers(markersGEC, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
  
  for (n in 1:25) {
    ClustermarkerList <- paste0("markerList$C",n)
    file <-  paste0("~/ArchR/All/cs17-out/Plots/Pseudobulk/","ClusterMarker_C",n,"_.csv") 
    
    #ClustermarkerList <- as.data.frame(Clustermarker)
    #sink("output.txt")
    #eval(parse(text=ClustermarkerList))
    #sink()
    
    write.csv(eval(parse(text=ClustermarkerList)), file = file)
    
 
  }

heatmapGEC <- markerHeatmap(
  seMarker = markersGEC, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE,
  returnMat=TRUE
)

heatmapGEC

ComplexHeatmap::draw(heatmapGEC, heatmap_legend_side = "bot", annotation_legend_side = "bot")
}


#PeakCAlling with MACS2 and adding MarkerPeaks
{
#pathToMacs2 <- "/System/Volumes/Data/Users/XXX/opt/anaconda3/bin/macs2"
pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "ClusterOneControl",
  pathToMacs2 = pathToMacs2
)

#Adding PeakMatrix - Only one Peakmatrix in one ArchR project possible
proj<-addPeakMatrix(
  ArchRProj = proj,
  ceiling = 4,
  binarize = FALSE,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addPeakMatrix")
)

#get Infos of Object
getPeakSet(proj)
table(proj$Clusters)
getAvailableMatrices(proj)

}

#Markerpeak Analysis
#Identify marker peaks
{
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "ClusteroneControl",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks

#show only relevant marker peaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")
markerList

}

#draw marker peaks as heatmap
{
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 3",
  transpose = TRUE
)

options(bitmapType ='cairo')

heatmapPeaks 

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", ArchRProj = proj, width = 8, height = 6, addDOC = FALSE, use_raster=TRUE, res=72)


#draw marker peaks as MA plot
pma1 <- plotMarkers(seMarker = markersPeaks, name = "TNT", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma1
#draw marker peaks as Vulcanoplot
pv1 <- plotMarkers(seMarker = markersPeaks, name = "TNT", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
pv1
#Plot
plotPDF(pma1, pv1, name = "Markers-MA-Volcano-TNT", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)


}


# Pairwise Testing Between Groups 11.3 + 12.1
{
  for (n in 1:25){
    #Geneexpression TNTvsCT
    usegroup <- paste0("C",n,"_x_TNT") 
    file <-  paste0("//Users/eliascrapa/ArchR/All/additional/motifs2/","Peakdiffereneces_C",n,"_x_TNT","vs_CTR-Markers-MA-Volcano")  
    motifsenriched <- paste0("PlotsMarkerMotifsOneControl/","Markers-Motifs-Enriched_C",n,"_x_TNT","vs_allCTR")
    
    motifsUpfile <- paste0("//Users/eliascrapa/ArchR/All/additional/motifs2/","Table_Markermotifs_up_Enriched_C",n,"_x_TNT","vs_allCTR.csv")
    motifsDofile <- paste0("//Users/eliascrapa/ArchR/All/additional/motifs2/","Table_Markermotifs_down_Enriched_C",n,"_x_TNT","vs_allCTR.csv")
    
    markerTest <- getMarkerFeatures(
      ArchRProj = proj, 
      useMatrix = "PeakMatrix",
      groupBy = "ClusterOneControl",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = usegroup,
      maxCells = 5000,
     
    )
    
    
    # pma2 <- plotMarkers(seMarker = markerTest, name = usegroup, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
    #pma2
    
    # pv2 <- plotMarkers(seMarker = markerTest, name = usegroup, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
    #pv2
    
    #plotPDF(pma2, pv2, name = file , width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
    
  }   

# Motif Plotting Colorful - https://github.com/GreenleafLab/ArchR/issues/550
{
  PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
      }

library(seqLogo)
##Chose which Motif to plot
plottingmotif <- paste0("Egr3_183")

#MHC Cluster9 Top10
tf <- c("Smarcc1_843","Bach1_108", "Fos_104", "Jund_135", "Pgr_849", "Bach2_119", "Nr3c1_673", "Mef2a_640","Mef2d_842","Mef2b_882")
#TNT Cluster9 Top10
tf <- c("Smarcc1_843","Bach1_108", "Fos_104", "Jund_135", "Bach2_119", "Fosb_98", "Junb_127", "Fosl1_107","Nfe2_132","Nfe2l2_101")
for (plottingmotif in tf) {
plottingmotiffile <- paste0("/Users/eliascrapa/ArchR/All/additional/Motif_",plottingmotif,".png")

png(plottingmotiffile, 490, 350)
seqLogo(PWMatrixToProbMatrix(getPeakAnnotation(proj, "Motif")$motifs[[plottingmotif]]))
dev.off()
}
}


#Motif Enrichment in Differential Peaks 
  
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force=TRUE)

{ 

  enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 1"
  )
  
  markersPeaks
  enrichMotifs
  
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 20, transpose = TRUE)
  
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  
  plotPDF(heatmapEM, name = "Reanalysis/Motifs-Enriched-Marker-Heatmap_Clusters", width = 20, height = 6, ArchRProj = proj, addDOC = FALSE)
}

#ChromVAR Deviatons Enrichment with ArchR
#Motif Deviations
{
  if("Motif" %ni% names(proj@peakAnnotation)){
    proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
  }
  
  proj<- addBgdPeaks(proj)
  
  proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = "Motif",
    force = TRUE
  )
  

  
  plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
  
  plotVarDev
  
  plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
  
  #extract a subset of motifs for downstream analysis
  motifs <- c("GATA4", "TCF21", "Gm4881", "Sfpi1", "Mef2a", "Mef2d", "Mef2b", "Elf1", "Bcl11a", "Spic")
  markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
  markerMotifs
  
  #remove motifs that we don't want
  #markerMotifs <- grep("z:", markerMotifs, value = TRUE)
  #markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
  #markerMotifs
  
  #Plot ChromVar deviation scores for each cluster
  
  p <- plotGroups(ArchRProj = proj, 
                  groupBy = "Sample", 
                  colorBy = "MotifMatrix", 
                  name = markerMotifs,
                  imputeWeights = getImputeWeights(proj)
                  
  )
  
  #Cowplot
  
  p2 <- lapply(seq_along(p), function(x){
    if(x != 1){
      p[[x]] + guides(color = FALSE, fill = FALSE) + 
        theme_ArchR(baseSize = 6) +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
        theme(
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank()
        ) + ylab("")
    }else{
      p[[x]] + guides(color = FALSE, fill = FALSE) + 
        theme_ArchR(baseSize = 6) +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
        theme(
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank()
        ) + ylab("")
    }
  })
  do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
  
  plotPDF(p, name = "Plot-Groups-Deviations-w-ImputationbyClusterbySample", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

  #overlay the z-scores on our UMAP embedding
  
  p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(proj)
  )
  
  #another Cowplot 
  p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6.5) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
      )
  })
  do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
  
  plotPDF(p2, name = "Plot-z-scores of MotifMatrix on our UMAP embedding2", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
  
  
  # TF deviation z-scores compare to the inferred gene expression via gene scores
  
  markerRNA <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
  markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
  markerRNA
  
  p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = sort(markerRNA), 
    embedding = "UMAP_Harmony",
    imputeWeights = getImputeWeights(proj),
    
    )
    
    p2 <- lapply(p, function(x){
      x + guides(color = FALSE, fill = FALSE) + 
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()
        )
    })
    do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
    
    plotPDF(p2, name = "Plot-z-scores compare to the inferred gene expression via gene scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
    
    
  #Similarly, because we previously linked our scATAC-seq data with corresponding scRNA-seq data,
  #we can plot the linked gene expression for each of these TFs on the UMAP embedding.
    
    markerRNA <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
    markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
    markerRNA
  
    p <- plotEmbedding(
      ArchRProj = proj, 
      colorBy = "GeneExpressionMatrix", 
      name = sort(markerRNA), 
      embedding = "UMAP_Harmony",
      continuousSet = "blueYellow",
      imputeWeights = getImputeWeights(proj)
    )
  
    
    p2 <- lapply(p, function(x){
      x + guides(color = FALSE, fill = FALSE) + 
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()
        )
    })
    do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
    
    
    plotPDF(p2, name = "Plot scatac linked gene expression for each of these TFs on the UMAP embedding", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
    
    
  }

#13.2 ArchR and Custom Deviations - No Encodedataset


###DIFFERENTIAL EXPRESSION WITH MOTIFMATRIX vis ChromVAR Deviatons Enrichment

# Pairwise Testing Between Groups 11.3 + 12.1
{
  for (n in 9:9){
    #Geneexpression TNTvsCT
    usegroup <- paste0("C",n,"_x_TNT") 
    bdggroups <- paste0("C",n,"_x_CTR")
    file <-  paste0("/Reanalysis/ChromvarDeviation/","ChromvarMotifdiffereneces_C",n,"_x_TNT","vs_CTR-Markers-MA-Volcano")  
    motifsenriched <- paste0("/Reanalysis/ChromvarDeviation/","ChromvarMotifs-Enriched_C",n,"_x_TNT","vs_CTR")
    motifsUpfile <- paste0("//Users/eliascrapa/ArchR/All/additional/CHROMVARdevDeviation/","ChromvarTable_Markermotifs_up_Enriched_C",n,"_x_TNT","vs_CTR.csv")
    motifsDofile <- paste0("//Users/eliascrapa/ArchR/All/additional/CHROMVARdevDeviation/","ChromvarTable_Markermotifs_down_Enriched_C",n,"_x_TNT","vs_CTR.csv")
    
    #usegroup <- paste0("C",n,"_x_MHC") 
    #bdggroups <- paste0("C",n,"_x_CTR-MHC")
    #file <-  paste0("/Pairwise/","Peakdiffereneces_C",n,"_x_MHC","vs_CTR-Markers-MA-Volcano")  
    #motifsenriched <- paste0("/Pairwise/New/","Markers-Motifs-Enriched_C",n,"_x_MHC","vs_CTR")
    #motifsUpfile <- paste0("//Users/eliascrapa/ArchR/All/cs17-out/Plots//Pairwise/New/","Table_Markermotifs_up_Enriched_C",n,"_x_MHC","vs_CTR.csv")
    #motifsDofile <- paste0("//Users/eliascrapa/ArchR/All/cs17-out/Plots/Pairwise/New/","Table_Markermotifs_down_Enriched_C",n,"_x_MHC","vs_CTR.csv")
    
    
    
    diffMotif <- getMarkerFeatures(
      ArchRProj = proj, 
      testMethod = "wilcoxon",
      useGroups = usegroup,
      bgdGroups = bdggroups,
      binarize = FALSE,
      useMatrix = "MotifMatrix",
      groupBy = "ClusterOneControl",
      useSeqnames="z"
    )
    
    diffMotifList <- getMarkers(diffMotif, cutOff = "FDR <= 0.01 & MeanDiff >= 0.1")
    diffMotifList
    diffMotifList$C9_x_TNT
    
    diffMotif
    motifs <- c("Smarcc1","Bach1", "Fos", "Jund", "Bach2", "Fosb", "Junb", "Fosl1","Nfe2","Nfe2l2")
    diffMotifs <- getFeatures(proj, select = NULL, useMatrix = "MotifMatrix")
    diffMotifs
    
    diffMotifs <- grep("z:", diffMotif, value = TRUE)
    markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
    diffMotifs
    
    #pma2 <- plotMarkers(seMarker = diffMotif, name = usegroup, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
    #pma2
    
  }
    # pv2 <- plotMarkers(seMarker = markerTest, name = usegroup, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
    #pv2
    
    #plotPDF(pma2, pv2, name = file , width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
    
    
    
    
    # MOTIFMATRIX Differential Analysis Part 2
    
    #proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force=TRUE)
    
    #save.image(file='After_Motifannotations.RData')
    
    
    motifsUp <- peakAnnoEnrichment(
      seMarker = diffMotif,
      ArchRProj = proj,
      peakAnnotation = "Motif",
      cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
    )
    
    motifsUp
    
    df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
    df <- df[order(df$mlog10Padj, decreasing = TRUE),]
    df$rank <- seq_len(nrow(df))
    
    head(df)
    
    #write Dataframe with Enriched Motifs Up
    write.csv(df, file = motifsUpfile)
    write.csv(df, file = motifsDofile)
    
    
    ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
      geom_point(size = 1) +
      ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
      ) + theme_ArchR() + 
      ylab("-log10(P-adj) Motif Enrichment") + 
      xlab("Rank Sorted TFs Enriched") +
      scale_color_gradientn(colors = paletteContinuous(set = "comet"))
    
    ggUp
    
    # analyses for the peaks that are more accessible in second Cluster
    motifsDo <- peakAnnoEnrichment(
      seMarker = markerTest,
      ArchRProj = proj,
      peakAnnotation = "Motif",
      cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
    )
    
    motifsDo
    
    df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
    df <- df[order(df$mlog10Padj, decreasing = TRUE),]
    df$rank <- seq_len(nrow(df))
    
    head(df)
    
    #write Dataframe with Enriched Motifs Up
    write.csv(df, file = motifsDofile)
    
    ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
      geom_point(size = 1) +
      ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
      ) + theme_ArchR() + 
      ylab("-log10(FDR) Motif Enrichment") +
      xlab("Rank Sorted TFs Enriched") +
      scale_color_gradientn(colors = paletteContinuous(set = "comet"))
    
    ggDo
    
    
    #plotPDF(ggUp, ggDo, name = motifsenriched, width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
    
  }
  
  
}



#14.1 Motif Footprinting & 14.2 Normalization of Footprints for Tn5 Bias

{
  
  motifPositions <- getPositions(proj)
  
  motifPositions
  
  motifs <- c("Smarcc1","Bach1", "Fos", "Jund", "Bach2", "Fosb", "Junb", "Fosl1","Nfe2","Nfe2l2")
  markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
  markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
  markerMotifs
  
  # Pseudo bulking - already done before - proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")
  
  seFoot <- getFootprints(
    ArchRProj = proj, 
    positions = motifPositions[markerMotifs], 
    groupBy = "Clusters"
  )
  
  # Subtracting the Tn5 Bias
  
  plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj, 
    normMethod = "Subtract",
    plotName = "Reanalysis/Footprints-Subtract-Bias",
    addDOC = FALSE,
    smoothWindow = 5
  )
  
  plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj, 
    normMethod = "Divide",
    plotName = "Footprints-Divide-Bias",
    addDOC = FALSE,
    smoothWindow = 5
  )
  
  plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj, 
    normMethod = "None",
    plotName = "Footprints-No-Normalization",
    addDOC = FALSE,
    smoothWindow = 5
  )
  
  
}


#14.3 Feature Footprinting
{
  
  seTSS <- getFootprints(
    ArchRProj = proj, 
    positions = GRangesList(TSS = getTSS(proj)), 
    groupBy = "Clusters",
    flank = 2000
  )
  
  plotFootprints(
    seFoot = seTSS,
    ArchRProj = proj, 
    normMethod = "None",
    plotName = "TSS-No-Normalization",
    addDOC = FALSE,
    flank = 2000,
    flankNorm = 100
  )
  
  
}



#15.2 Co-accessibility with ArchR

{
  
  proj <- addCoAccessibility(
    ArchRProj = proj,
    reducedDims = "Harmony"
  )
  
  cA <- getCoAccessibility(
    ArchRProj = proj,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
  )
  
  cA
  
  metadata(cA)[[1]]
  
  #either resolution 1
  cA <- getCoAccessibility(
    ArchRProj = proj,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = TRUE
  )
  
  
  cA[[1]]
  
  # or decrease the resolution of our loops to 1000. Can help with over-plotting of co-accessibility interactions
  cA <- getCoAccessibility(
    ArchRProj = proj,
    corCutOff = 0.5,
    resolution = 1000,
    returnLoops = TRUE
  )
  
  cA[[1]]
  
  #even fewere if resolution set to 10000
  
  cA <- getCoAccessibility(
    ArchRProj = proj,
    corCutOff = 0.5,
    resolution = 10000,
    returnLoops = TRUE
  )
  
  cA[[1]]
  
  
  
  idxSample <- BiocGenerics::which(proj$Samples %in% "CTR-MHC" %in% "CTR-TNT")
  cellsSample <- proj$cellNames[idxSample]
  projCTR <- proj[cellsSample, ]
  
  #15.2.1 Plotting browser tracks of Co-accessibility
  
  markerGenes  <- c(
    "Tgfb2", "Tgfb1", "Tgfb3", #Fibroblasts
    "Vegfa", "Vegfc", #Cardimyocyte
    "Fgf1", "Pdgfd", "Pdgfb", "Postn", #B-Cell Trajectory
    "Ptprm", "Angpt1", "Lama2" #Monocytes
  )
  
  
  markerGenes  <- c(
    "Col5a1", "Col3a1", #Fibroblasts
    "MHRT", "Tnnt2", #Cardimyocyte
    "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
    "CD14", "CEBPB", "MPO", #Monocytes
    "IRF8", "Postn", "TCF21", #Fibroblasts
    "CD3D", "CD8A", "TBX21", "IL7R", #TCells
    "Pecam1", "Itgb1", "CD34" #Endothelial Cells
    "Hbb"
  )
  
  markerGenes  <- c(
    "Col5a1", "Col3a1", #Fibroblasts
    "MHRT", "Tnnt2", #Cardimyocyte
    "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
    "CD14", "CEBPB", "MPO", #Monocytes
    "IRF8", "Postn", "TCF21", #Fibroblasts
    "CD3D", "CD8A", "TBX21", "IL7R", #TCells
    "Pecam1", "Itgb1", "CD34", #Endothelial Cells
    "Nrg1", "Mecom", "Pde3a", "Prkg1","Acer2", #Cluster16,
    "Ptprg", "Flt1", "Zbtb20", 
    "Slc39a11", "Hbb-bs","Tnnt2","Prkg1","Dmd","Sorbs2", #TNTC9
    "Sorbs2", "Pdlim5","Celf2","Ctnna3","Pcdh7", #MHCC9
    "Ptprg", "Tshz2", "Plcb1", "Magi1", "Pitpnc1", "Tcf4", # C14TNT
    "Ptprg","Flt1","Pkp4","Zbtb16","Plcb4", #MHC C14
    "Zbtb16", "Pcdh9", "Cfh", "Celf2", "Fmo2", #MHC_C21
    "Nrg1","Cfh","Col8a1","Slit3","Pde3a","Heg1","Nfia","Lama2", # C21TNT
    "Col8a1", "Celf2","Cfh", "Lama2", "Zbtb20","Gab2","Mast4", "Rora", #C22 TNT
    "Zbtb16","Celf2","Errfi1","Rora","Pcdh9","Zbtb20","Lhfp","Cfh","Gpc6" #C22 MHC
  )
  
  p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "ClusterOneControl",
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(proj, corCutOff = 0.5,
                               resolution = 10000,
                               returnLoops = TRUE)
  )
  
  grid::grid.newpage()
  grid::grid.draw(p$Postn)
  
  plotPDF(plotList = p, 
          name = "Reanalysis/ClusteroneControl-BrianBlack-Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
          ArchRProj = proj, 
          addDOC = FALSE, width = 5, height = 5)
  
  
}


#15.3 Peak2GeneLinkage with ArchR

{


proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "Harmony",
  useMatrix = "GeneExpressionMatrix",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  k = 20,
  knnIteration = 500,
  overlapCutoff = 0.8,
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.4,
  seed = 1,
  threads = max(floor(getArchRThreads()/2), 1),
  verbose = TRUE,
  logFile = createLogFile("addPeak2GeneLinks")
)


  
  p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE
  )
  
  p2g
  
  metadata(p2g)[[1]]
  
  p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = TRUE
  )
  
  p2g[[1]]
  
  #decrease the resolution of these links by setting resolution = 1000
  
  p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = 0.45,
    resolution = 1000,
    returnLoops = TRUE
  )
  
  p2g[[1]]
  
  # or again reduce it to 10000
  p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = 0.45,
    resolution = 10000,
    returnLoops = TRUE
  )
  
  p2g[[1]]
  
  markerGenes  <- c(
    "Col5a1", "Col3a1", #Fibroblasts
    "MHRT", "Tnnt2", #Cardimyocyte
    "Nppa", "Myl7", "Sln", "Myh6", #atrial Myocyte
    "PAX5", "MS4A1", "EBF1", "MME", "Cd79a", #B-Cell Trajectory
    "Jchain", #"Igk-C", "Igh-6", #Plasmacyte
    "CD14", "CEBPB", "MPO", #Monocytes
    "C1qa", "Vsig4",#Macrophages
    "S100a9","S100a8", #Neurophils
    "IRF8", "Postn", "TCF21","Mfap5", "Col3a1", "Col5a1","Pdgfra", #Fibroblasts
    "H19", "Igf2", "Vwf", #Endocardial Cell
    "CD3D", "CD8A", "TBX21", "IL7R", "CD8A", "Ccl5", #Adult CD8+ T-Cell
    "Cd209a", "H2-Ab1", #Dendritic  Cells
    "Pecam1", "Itgb1", "CD34", "Lyve1", "Plvap", "Vtn", #Endothelial Cells Lyve1, Plvap and Vtn high clusters described
    "Gpihbp1", "Cd36", #Capillary endothelial cell_Gpihbp1 high-C2
    "Acta2", "Tagln", "Mylk", "Rgs5", #Smooth muscle_cells - RGs5 and ACTA2 high types possible
    "Hba-a1", "Alas2", #Erythroid Cells
    "Retnla", "Lyz1", #Eosinophil
    "Car3", "Cfd", "Adipoq" # Adipose Cell
  )
  
  markerGenes  <- c("Prkg1","Acer2", #Cluster16,
    "Ptprg", "Flt1", "Zbtb20", 
    "Slc39a11", "Hbb-bs","Tnnt2","Prkg1","Dmd","Sorbs2", #TNTC9
    "Sorbs2", "Pdlim5","Celf2","Ctnna3","Pcdh7", #MHCC9
    "Ptprg", "Tshz2", "Plcb1", "Magi1", "Pitpnc1", "Tcf4", # C14TNT
    "Ptprg","Flt1","Pkp4","Zbtb16","Plcb4", #MHC C14
    "Zbtb16", "Pcdh9", "Cfh", "Celf2", "Fmo2", #MHC_C21
    "Nrg1","Cfh","Col8a1","Slit3","Pde3a","Heg1","Nfia","Lama2", # C21TNT
    "Col8a1", "Celf2","Cfh", "Lama2", "Zbtb20","Gab2","Mast4", "Rora", "Postn", #C22 TNT
    "Zbtb16","Celf2","Errfi1","Rora","Pcdh9","Zbtb20","Lhfp","Cfh","Gpc6" #C22 MHC 
  )
  p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(proj, corCutOff = 0.45,
                              resolution = 1000,
                              returnLoops = TRUE)
  )
  
  grid::grid.newpage()
  grid::grid.draw(p$CDPostn)
  
  plotPDF(plotList = p, 
          name = "Reanalysis/BrianBlack-Plot-Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
          ArchRProj = proj, 
          addDOC = FALSE, width = 5, height = 5)
  
  #Plot Comparison Heatmaps of scATAC and scRNA
  
  p <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "Clusters")
  
  p
  
  plotPDF(p, name = "Plot Heatmaps of scRNA and scATAQ beside each other", width = 8, height = 8, ArchRProj = proj, addDOC = FALSE)
  
}


#15.4 Identification of Positive TF-Regulators
{
  projTNT <- proj[proj@cellColData@listData[["Samples"]] == "TNT" ]
  projMHC <- proj[proj@cellColData@listData[["Samples"]] == "MHC" ]
  
group  <- list()
group  <- c("TNT","MHC")
  
for (u in group) {
    
      for ( n in 9:9){
    
  idxSampleMUT <- BiocGenerics::which(proj$ClusterOneControl %in% paste0("C",n,"_x_",u))
  idxSampleCTR <- BiocGenerics::which(proj$ClusterOneControl %in% paste0("C",n,"_x_CTR"))
  
  cellsSampleMUT <- proj$cellNames[idxSampleMUT]
  cellsSampleCTR <- proj$cellNames[idxSampleCTR]
 
  
  
 # idxSample <- BiocGenerics::which(proj$ClusterOneControl %in% paste0("C",n,"_x_",u))
  #cellsSample <- proj$cellNames[idxSample]
  #projsub <- proj[cellsSample, ]
  
  
  
 # seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "ClusterOneControl")
  
#  seGroupMotif
  
#  seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

#  rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
#    rowMaxs(assay(seZ) - assay(seZ)[,x])
#  }) %>% Reduce("cbind", .) %>% rowMaxs
  

 # lapply(seq_len(ncol(seZ)), function(x){     #https://github.com/GreenleafLab/ArchR/issues/354
 #   rowMaxs(assay(seZ) - assay(seZ)[,x])
 # }) %>% Reduce("cbind", .) %>% rowMaxs  
  
  corGSM_MM <- correlateMatrices(
    ArchRProj = proj[c(cellsSampleMUT,cellsSampleCTR),],
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "Harmony"
  )
  
  corGSM_MM
  
  corGIM_MM <- correlateMatrices(
    ArchRProj = proj[,c(idxSampleMUT,idxSampleCTR)],
    useMatrix1 = "GeneExpressionMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "Harmony"
  )
  
  corGIM_MM
  
  corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
  corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
  
  # Identify Positive TF Regulators
  corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
  corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
  corGSM_MM$TFRegulator <- "NO"
  corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
  sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
  
  
  pGS <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
    geom_point() + 
    theme_ArchR() +
    geom_vline(xintercept = 0, lty = "dashed") + 
    scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
    xlab("Correlation To Gene Score") +
    ylab("Max TF Motif Delta") +
    scale_y_continuous(
      expand = c(0,0), 
      limits = c(0, max(corGSM_MM$maxDelta)*1.05)
    )
  
  pGS
  
  TFregulators <- as.data.frame(sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1]))
  write_csv(TFregulators, file = paste0("C",n,"_x_",u,"-TFregulators_GeneScore.csv"))
  
  # same analysis for the correlations derived from our GeneexpressionMatrix
  
  corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
  corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
  corGIM_MM$TFRegulator <- "NO"
  corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01  & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
  sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
  
  TFregulators <- as.data.frame(sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1]))
  write_csv(TFregulators, file = paste0("/Users/eliascrapa/ArchR/All/cs17-out/Plots/Reanalysis/C",n,"_x_",u,"-TFregulators_GeneExpression.csv"))
  
  pGE <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
    geom_point() + 
    theme_ArchR() +
    geom_vline(xintercept = 0, lty = "dashed") + 
    scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
    xlab("Correlation To Gene Expression") +
    ylab("Max TF Motif Delta") +
    scale_y_continuous(
      expand = c(0,0), 
      limits = c(0, max(corGIM_MM$maxDelta)*1.05)
    )
  
  pGE
  
  
  plotPDF(pGS, pGE, name = paste0("Reanalysis/C",n,"_x_",u,"-Cluster_Dot Plots GeneScore and GeneExpression based pos TF-Regulators"), width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
  write.csv( corGSM_MM, file = paste0("/Users/eliascrapa/ArchR/All/cs17-out/Plots/Reanalysis/C",n,"_x_",u,"-Cluster_Identification of Positive TF-Regulators - Genescore.csv"))
  write.csv( corGIM_MM, file = paste0("/Users/eliascrapa/ArchR/All/cs17-out/Plots/Reanalysis/C",n,"_x_",u,"-Sample_Identification of Positive TF-Regulators - Geneexpression.csv"))
  
  
                  }
  
  
              }
  
  
  
  }

# Chapter 16 Trajectory Analysis with ArchR
{
# Plotting old and new names
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP_Harmony")

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP_Harmony")

ggAlignPlots(p1, p2 type = "h")

trajectory <- c("C22", "C21", "C20")
trajectory

proj <- addTrajectory(
  ArchRProj = proj, 
  name = "Fibroblasts", 
  groupBy = "Clusters",
  trajectory = trajectory, 
  embedding = "UMAP_Harmony", 
  force = TRUE
)

head(proj$Fibroblasts[!is.na(proj$Fibroblasts)])

p <- plotTrajectory(proj, trajectory = "Fibroblasts", colorBy = "cellColData", name = "Fibroblasts", embedding = "UMAP_Harmony")

p[[1]]

plotPDF(p, name = "Plot-Fibroblas-Traj-UMAP.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

}


cM <- confusionMatrix(proj$Clusters, proj$predictedGroup)
labelOld <- (proj$Clusters)
labelOld
#  [1] "Cluster11" "Cluster2"  "Cluster12" "Cluster1"  "Cluster8"  "Cluster4" 
#  [7] "Cluster9"  "Cluster5"  "Cluster7"  "Cluster14" "Cluster3"  "Cluster10"
# [13] "Cluster6"  "Cluster13"

labelNew <- proj$Clusters[apply(proj$Clusters, 1, which.max)]
labelNew

proj$Clusters

#SdafsAVE Projectevniroment
save.image(file='All Steps finished.RData')


#Print Session Info
sessionInfo()
