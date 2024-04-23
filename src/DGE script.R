suppressPackageStartupMessages({
  library(Seurat)
  library(venn)
  library(dplyr)
  library(cowplot)
  library(ggplot2)
  library(pheatmap)
  library(enrichR)
  library(rafalib)
  library(naturalsort)
  library(openxlsx)
})


setwd("C:\\Users\\tthottakara\\OneDrive - UCSF\\DATA\\SingleCellSeq\\Analysis\\All")

#Idents(object = all) <- all@meta.data$ArchR_Clusters

Idents(object = all) <- "ArchR_Clusters"
my_levels <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25")
Idents(all) <- factor(Idents(all), levels= my_levels)

Idents(object = all) <- "ClusterbySample"
my_levels <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25")
Idents(all) <- factor(Idents(all), levels= my_levels)
write.xlsx(table(all@active.ident),file="cellspercluster_ClusterbySample.xlsx")

Idents(object = all) <- all@meta.data[["ClusterOneControl"]]
my_levels <-naturalsort(levels(factor(all@meta.data[["ClusterOneControl"]])))
Idents(all) <- factor(Idents(all), levels= my_levels)


#Idents(object = all) <- "clusternames"
Idents(object = all) <- "ArchR_Clusters"
my_levels <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25")
Idents(all) <- factor(Idents(all), levels= my_levels)

# my_levels <- c("C1", "Macrophage_C2", "Monocytes_C3", "CD8TCells_C4", "BCells_C5", "atrialCM_C6", 
#                "CM2_C7", "FGF13+CM_C8", "ventricularCM_C9", "C10", "C11", 
#                "C12", "ET1_C13", "mainET2_C14", "ET3_C15", "venousET4_C16", 
#                "lymphaticET5_C17", "VSMC1_C18", "VSMC2_C19", "mainVSMC3_C20", "FB1_C21",
#                "FB2_C22", "FB3_C23", "CM5?_C24", "C25")
# Idents(all) <- factor(Idents(all), levels= my_levels)

#Idents(object = all) <- "ArchR_Clusters"
new.cluster.ids <- c("Leukocyte", "Leukocyte", "Leukocyte", "Leukocyte", "Leukocyte", "Cardiomyocyte", 
               "Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", "EpicardialCell", "EpicardialCell", 
               "EpicardialCell", "EndothelialCell", "EndothelialCell", "EndothelialCell", "EndothelialCell", 
               "EndothelialCell", "VSMC", "VSMC", "VSMC", "Fibroblast",
               "Fibroblast", "Fibroblast", "undeterm", "undeterm")
names(new.cluster.ids) <- levels(all)
all <- RenameIdents(all, new.cluster.ids)

all@meta.data$Celltype <- all@active.ident
all@meta.data$CelltypesbyGroup <- paste0(all@meta.data$Celltype, "_x_", all@meta.data$groups)


Idents(object = all) <- all@meta.data$CelltypesbyGroup
Idents(object = all) <- all@meta.data$Celltype 

data.dir <- './NewDGE'
dir.create(data.dir)
setwd(data.dir)

# plot this clustering
png("ClusterPlots.png", width = 1920, height = 960)
plot_grid(ncol = 3, DimPlot(all, label = T, reduction= 'harmony_archr' ) + NoAxes(), DimPlot(all, group.by = "orig.ident",  reduction= 'harmony_archr') + 
            NoAxes(), DimPlot(all, group.by = "ArchR_Clusters",  reduction= 'harmony_archr') + NoAxes())
dev.off()


DefaultAssay(all) <- "RNA"

# Compute differentiall expression
markers_genes <- FindAllMarkers(all, logfc.threshold = 0.2, test.use = "DESeq2", 
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 1000, 
                                assay = "RNA", slot = "counts")

marker_genes_bp <- markers_genes
markers_genes <- markers_genes[!grepl("^mt-",markers_genes$gene),]

#We can now select the top 25 up regulated genes for plotting.

top25 <- markers_genes %>% group_by(cluster) %>% top_n(-25, p_val_adj)
factor(top25, levels= my_levels)
top25

write.xlsx(markers_genes,"ClusterMarkerGenes_DESseq.xlsx")
write.csv(top25,"ClusterMarkerGenes_DEsSeq_Top25.csv")
#We can now select the top 25 up regulated genes for plotting.

mypar(2, 5, mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
  BarPlot <- paste0("BarPLot_ClusterGenes_C_",i,".png")
  
  png(BarPlot)
  barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F), horiz = T, 
          las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
  abline(v = c(0, 0.25), lty = c(1, 2))
  dev.off()
}
  
  

# We can visualize them as a heatmap. Here we are selecting the top 5.

top5 <- markers_genes %>% group_by(cluster) %>% slice_min(p_val_adj, n = 5, with_ties = FALSE) 

levels(all) <- my_levels
# create a scale.data slot for the selected genes
all <- ScaleData(all, features = as.character(unique(top5$gene)), assay = "RNA")
png("HearMapClusterGenes_Top3.png", width = 2400, height = 1440)

gg <- DoHeatmap(all, features = as.character(unique(top5$gene)), 
          assay = "RNA", group.by = "ArchR_Clusters") + scale_fill_gradientn(colors = c("blue", "black", "yellow")) 
ggdev.off()

ggsave(plot=gg, filename ="HeatMapClusterGenes_Top5_300.pdf", device="pdf", dpi="print", width = 30, height = 25, unit="cm")

#Another way is by representing the overall group expression and detection rates in a dot-plot.

png("DotPlotClusters.png",width = 960, height = 1440)
DotPlot(all, features = rev(as.character(unique(top5$gene))), "ArchR_Clusters", 
        assay = "RNA") + coord_flip() 
dev.off()



#We can also plot a violin plot for each gene.

# take top 3 genes per cluster/
top3 <- top5 %>% group_by(cluster) %>%  slice_min(p_val_adj, n = 3, with_ties = FALSE)

png("DotPlotClusters_Top3.png",width = 960, height = 1440)
DotPlot(all, features = rev(as.character(unique(top3$gene))), "ArchR_Clusters", 
        assay = "RNA",  scale = FALSE, col.max = 1) + coord_flip() 
dev.off()

# take top 3 genes per cluster/
top1 <- top5 %>% group_by(cluster) %>%  slice_min(p_val_adj, n = 1, with_ties = FALSE)

png("DotPlotClusters_Top1.png",width = 960, height = 1440)
DotPlot(all, features = rev(as.character(unique(top1$gene))), "ArchR_Clusters", 
        assay = "RNA",  scale = FALSE, col.max = 1) + coord_flip() 
dev.off()


# set pt.size to zero if you do not want all the points to hide the violin
# shapes, or to a small value like 0.1

png("VlnPlots.png", width = 2880, height = 2160)
VlnPlot(all, features = as.character(unique(top3$gene)), ncol = 5, 
        assay = "RNA", pt.size = 0)#, group.by = "ArchR_Clusters"
dev.off()



##RE-Do Differential Expresion with grouped Celltypes
Clusternumber <- levels(Idents(all))

#CompareClusterGroups <- function(n,m,bggroup,testgroup){
for (n in Clusternumber){
  
  #Geneexpression TNTvsCT
  #  value <- paste0(groups[[n]])
  n="Leukocyte"
  #n=2L
  m=n
  #n=4L
  #m=3L
  #bggroup="NS-MyHC"
  testgroup = "TnT"
  bggroup = "CTR"
  
  ident1<- paste0(n,"_x_",testgroup)  
  ident2<- paste0(m,"_x_",bggroup)
  filem <-  here("RNA", "DGE_NewClusters",(paste0("1_RNA_ClusterMarker_C",n,"_",testgroup,"_vs_C",m,"_",bggroup,".xlsx")))
  Volcano <-  here("RNA", "DGE_NewClusters",(paste0("Volcanoplot_C",n,"_",testgroup,"_vs_C",m,"_",bggroup,".png")))
  
  
  # find all markers of cluster 9
  
  cluster.markers <- FindMarkers(all, ident.1= ident1, ident.2 = ident2 , min.cells.feature = 5, test.use = "DESeq2", slot = "counts", assay = "RNA")
  names <- rownames(cluster.markers)
  
  
  
  #png device
  png(Volcano)
  
  print(EnhancedVolcano(cluster.markers,
                        lab = rownames(cluster.markers),
                        x = 'avg_log2FC',
                        y = 'p_val',
                        pCutoff = 0.05,
                        FCcutoff = 0.5,
                        xlim = -1.5,
                        ylim = 0.1)
        
  )
  # Close device
  dev.off()
  
  
  
  rownames(cluster.markers) <- NULL
  cluster.markers <- cbind(names,cluster.markers)
  
  #cluster.markers$rownames <- rownames(cluster.markers)
  write.xlsx(cluster.markers, file = filem)
  #head(cluster22.markers, n = 25, wt = avg_log2FC)   
  
} 

sapply(Clusternumber,CompareClusterGroups(n,m=n,bggroup="NS-TnT",testgroup = "S-TnT"))
sapply(Clusternumber,CompareClusterGroups(n,m=n,bggroup="NS-MyHC",testgroup = "S-MyHC"))









##RE-Do Differential Expresion with grouped Celltypes
Idents(object = all) <- all@meta.data$Celltype 
Clusternumber <- levels(Idents(all))
Idents(object = all) <- all@meta.data$CelltypesbyGroup


for (n in Clusternumber){
  #Geneexpression TNTvsCT
  #  value <- paste0(groups[[n]])
  ident1<- paste0(n,"_x_TNT")  
  ident2<- paste0(n,"_x_CTR")
  filem <-  paste0("C://Users//tthottakara//OneDrive - UCSF//DATA//SingleCellSeq//Analysis//All//NewDGE","Diffexpr_",ident1,"_vs_",ident2,".xlsx")  
  #Volcano <- paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/","Volcanoplot_TNT_vs_CONTROLS_C_",n,".png")
  
  #motifsenriched <- paste0("/Pairwise/New/","Markers-Motifs-Enriched_C",n,"_x_TNT","vs_CTR")
  #motifsUpfile <- paste0("//Users/eliascrapa/ArchR/All/cs17-out/Plots//Pairwise/New/","Table_Markermotifs_up_Enriched_C",n,"_x_TNT","vs_CTR.csv")
  #motifsDofile <- paste0("//Users/eliascrapa/ArchR/All/cs17-out/Plots/Pairwise/New/","Table_Markermotifs_down_Enriched_C",n,"_x_TNT","vs_CTR.csv")
  #print(value)
  
  # find all markers of cluster 9
  
  cluster.markers <- FindMarkers(all, ident.1= ident1, ident.2 = ident2 , min.cells.feature = 5, test.use = "DESeq2", slot = "counts", assay = "RNA")
  names <- rownames(cluster.markers)
  
  
  
  # png device
  #  png(Volcano)
  
  #  print(EnhancedVolcano(cluster.markers,
  #                        lab = rownames(cluster.markers),
  #                        x = 'avg_log2FC',
  #                        y = 'p_val',
  #                        #pCutoff = 0.05,
  #                       FCcutoff = 0.5,
  #                        xlim = -1.5,
  #                        ylim = 0.1)
  
  #  )
  # Close device
  #dev.off()
  
  
  
  rownames(cluster.markers) <- NULL
  cluster.markers <- cbind(names,cluster.markers)
  
  #cluster.markers$rownames <- rownames(cluster.markers)
  write.xlsx(cluster.markers, file = filem)
  #head(cluster22.markers, n = 25, wt = avg_log2FC)   
  
  
}



# Assuming 'object' is your Seurat object and 'CelltypesbyGroup' is the column name containing the cell type labels
barcodes <- rownames(all@meta.data)
labels <- all@meta.data$CelltypesbyGroup

# Combine the barcodes and labels into a data frame
df <- data.frame(Barcode = barcodes, CelltypesbyGroup = labels)
df$Barcode <- gsub("_","#",df$Barcode)

# Write the data frame to a CSV file
write.csv(df, file = "barcodes_labels.csv", row.names = FALSE)




genos <- c("TNT","MHC")
for (k in genos){ 
  for (i in 2:25){
    
    
    sample <- paste0("C",i,"_x_",k)
    file <- paste0("Marker_C",i,"_",k,"vsCTR")
    csv <- paste0("Marker",k,"vs CTR.csv")
    controlcluster <- paste0("C",i,"_x_CTR")
    
    markersPeaks <- getMarkerFeatures(
      ArchRProj = proj, 
      useMatrix = "PeakMatrix",
      groupBy = "ClusterOneControl",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = sample,
      bgdGroups = controlcluster
    )
    
    
    markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")
    markerList
    
    #markerpeaks for Sample
    write.csv(markerList$sample,file = csv)
    
    
    #draw marker peaks as MA plot
    pma1 <- plotMarkers(seMarker = markersPeaks, name = sample, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
    pma1
    
    
    #draw marker peaks as Vulcanoplot
    pv1 <- plotMarkers(seMarker = markersPeaks, name = sample, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
    pv1
    
    plotPDF(pma1, pv1, name = file, width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
    
    
  }
}

