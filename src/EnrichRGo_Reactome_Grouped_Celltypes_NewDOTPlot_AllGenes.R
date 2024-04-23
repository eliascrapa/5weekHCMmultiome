#https://ycl6.github.io/GO-Enrichment-Analysis-Demo/4_enrichR.html
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")

#library(enrichplot)
#library("clusterProfiler")
library("enrichplot")
library("org.Mm.eg.db")
library("ggplot2")
library(ggpubr)
library(openxlsx)
library(tidyverse)
library("enrichR")
library(cowplot)
library(here)
library(DOSE)
library(patchwork)
library(scales)

#setwd("C:/Users/tthottakara/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/NewDGE/TNTvsMHC/comparingsame/")
#setwd("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/MP/TnT")
setwd("~/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/EC/")
setwd("~/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/DGE/DotPlots")
setwd("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/EC/Geneset")
dotplot_function <- function(data) {
  df <- data
 #df <- dnEnriched_go[["GO_Biological_Process_2023"]]
  
  df$Term <- gsub(" \\(.*\\)", "", df$Term)
  
  df <- df %>% dplyr::filter(Adjusted.P.value <0.05) 
  
  df <- df %>% slice_min(n = 10, order_by = Adjusted.P.value, with_ties = FALSE)
  print(df$Term)
  
  Significant <- as.numeric(sub("/\\d+$", "", as.character(df$Overlap)))
  
  df$Generatio <- sapply(df$Overlap, function(x) eval(parse(text = x)))
  df$Term <- gsub(" WP\\d+", "", df$Term)
  df$Term <- factor(df$Term, levels = rev(df$Term[order(df$Adjusted.P.value)]))

  title <- paste0(group," Cluster ",n)
  
  try(ggplot(df, aes_string(x="Adjusted.P.value", y="Term", size="Significant", color="Generatio")) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = "Gene ratio",
                           guide=guide_colorbar(reverse=TRUE)) +
    ylab(NULL) + ggtitle(title) + theme_dose(12) +
    scale_size(range=c(3, 8), name="Gene Count") +
    guides(size  = guide_legend(order = 1),
           color = guide_colorbar(order = 2)) +
    theme(axis.text.y = element_text(size = 16, face = "bold")) +
    scale_x_reverse(labels = label_number(accuracy = 0.001)) + 
    xlab("adj. P value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.text.y=element_text(lineheight=0.8)) +
    scale_y_discrete(labels=function(x) str_wrap(x, width=30)))
}


test_ggplot_error <- function(ggplot_object) {
  tryCatch({
    # Try to print the ggplot object
    print(ggplot_object)
    # If no error occurs, return FALSE
    return(TRUE)
  }, error = function(e) {
    # If an error occurs, return TRUE
    return(FALSE)
  })
}

#########CHOOOSE group
# group = "Common Genes"
# upgroup <-  paste0("Enriched in Upregulated Common Genes")
# downgroup <-  paste0("Enriched in Downregulated Common Genes")
# group = "Up in MyHC Genes"
# upgroup <-  paste0("Enriched in Upregulated Genes in MyHC")
# downgroup <-  paste0("Enriched in Downregulated Genes in MyHC")
group = "TnT"
upgroup <-  paste0("Enriched in Upregulated Genes in TnT")
downgroup <-  paste0("Enriched in Downregulated Genes in TnT")
#1 is TnT, 2 is MyHC, 3 is Common Genes
q=1
termnumber = 10

#data.dir <- paste0(here("RNA","DGE_negbionomial",paste0("EnrichRGo_Reactome_",termnumber,'\\',group,'\\')))
#dir.create(data.dir)
#setwd(data.dir)

#celltypelevels <- c("Cardiomyocyte","EndothelialCell","EpicardialCell","Fibroblast","Leukocyte","undeterm","VSMC")

#filem <-  paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsCTR_C",n,".rnk")  
#title <- paste0("Venn Diagram Cluster ",n)

#TNT Export RNK File for Genesetenrichment

# clusternumbers <- c(1,2,3,4,5.6,8,9,12,13,14,15,16,17,19,20,21,22,23)
# clusternumbers <- c(5,6,8,9,12,13,14,15,16,17,19,20,21,22,23)
# clusternumbers <- c(22,23)
# clusternumbers <- c(9,12,14,16,17,19,20,21,22)
# clusternumbers <- c(19,20,21,22)
# clusternumbers <- c(22)

clusternumbers <- c(2,9,14,16,20,22)
#clusternumbers <- c(9)
#clusternumbers <- c(21,22)
clusternumbers <- c(22)
for (n in clusternumbers){
  
  #n="Cardiomyocyte"
  
  #n=2L
  m=n
  #q=1L
  #skip_to_next <- FALSE

#source <- paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsMHC/RNA_ClusterMarker_TNT_vs_MHC_C",n,".xlsx")
source <- paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_TNT_vs_CONTROLS_C",n,".xlsx")
#source <- paste0("~/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_MHC_vs_CONTROLS_C",n,".xlsx")
#source <- paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/DGE/Venn_Extraction/Venn Diagram Cluster ",n,"-Genelists.xlsx")
#source <- paste0("/Users/eliascrapa/Library/CloudStorage/Box-Box/SinglecellSeq/Analysis_Swimming/RNA/DGE_negbionomial/RNA_ClusterMarker_C",n,"_S-MyHC_vs_C",m,"_NS-MyHCnegbinomial.xlsx")
#source <- paste0("/Users/eliascrapa/Library/CloudStorage/Box-Box/SinglecellSeq/Analysis_Swimming/RNA/DGE_negbionomial/RNA_ClusterMarker_C",n,"_S-WT_vs_C",m,"_NS-WTnegbinomial.xlsx")
#upgroup <-  paste0("Enriched in Upregulated Genes ",group," ")
#downgroup <-  paste0("Enriched in Downregulated Genes ",group," ")

TNTvsMHC <- openxlsx::read.xlsx(source, q)

sapply(TNTvsMHC, class) 
#TNTvsMHC[,3:6] <- sapply(TNTvsMHC[,3:6],as.numeric)



#remove_rownames(TNT_logFC_sel)
#TNT_logFC_sel <- dplyr::arrange(TNTvsMHC, desc(logFC_A))


TNT_logFC_sel <-TNTvsMHC 
up.idx <- TNT_logFC_sel$names[which(TNT_logFC_sel$avg_log2FC > 0 & TNT_logFC_sel$p_val_adj < 0.01)]
# FDR < 0.05 and logFC > 0
dn.idx <- TNT_logFC_sel$names[which(TNT_logFC_sel$avg_log2FC < 0 & TNT_logFC_sel$p_val_adj < 0.01)]
# FDR < 0.05 and logFC > 0

#up.idx <- TNT_logFC_sel$Genes[TNT_logFC_sel$Trend == "UP"]
#dn.idx <- TNT_logFC_sel$Genes[TNT_logFC_sel$Trend == "DOWN"]

length(up.idx)
length(dn.idx)

if (length(up.idx) >0) {
  
  if (length(dn.idx) >0) {
  
# all.genes <- TNT_logFC_sel$names
# up.genes <- TNT_logFC_sel[up.idx,]$names
# dn.genes <- TNT_logFC_sel[dn.idx,]$names

up.genes <- up.idx
dn.genes <- dn.idx


# Use fromType = "ENSEMBL" if your input identifier is Ensembl gene ID
up.genes.df = as.data.frame(up.genes)
dn.genes.df = as.data.frame(dn.genes)

colnames(dn.genes.df) <- "names"
colnames(up.genes.df) <- "names"

head(up.genes.df, 10)
head(dn.genes.df, 10)

#load databases
dbs <- listEnrichrDbs()
dbs <- dbs[order(dbs$libraryName),]

dbs
class(dbs)
dim(dbs)
head(dbs)

#Check which databases are available
dbs$libraryName
dbs[grep("2021",dbs$libraryName),]$libraryName 

dbs_go <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
dbs_pw <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse")#, "BioPlanet_2019")
dbs_dd <- c("Reactome_2022")
dbs_pa <- c("Panther_2016")

#Enrichment of GO Terms

upEnriched_go <- enrichr(genes = up.genes.df$names, databases = dbs_go)
Sys.sleep(2)
dnEnriched_go <- enrichr(genes = dn.genes.df$names, databases = dbs_go)

#class(upEnriched_go)
#names(upEnriched_go)

# View top 5 terms in the first element of the list
#head(upEnriched_go[[1]], 15)


#Pathways
Sys.sleep(2)
upEnriched_pw <- enrichr(genes = up.genes.df$names, databases = dbs_pw)
Sys.sleep(2)
dnEnriched_pw <- enrichr(genes = dn.genes.df$names, databases = dbs_pw)

#head(upEnriched_pw[[3]], 10)

#Diseases/Drugs analysis
Sys.sleep(2)
upEnriched_dd <- enrichr(genes = up.genes.df$names, databases = dbs_dd)
Sys.sleep(2)
dnEnriched_dd <- enrichr(genes = dn.genes.df$names, databases = dbs_dd)

#head(dnEnriched_dd[[2]], 10)

Sys.sleep(2)
upEnriched_pa <- enrichr(genes = up.genes.df$names, databases = dbs_pa)
Sys.sleep(2)
dnEnriched_pa <- enrichr(genes = dn.genes.df$names, databases = dbs_pa)


#upEnriched_go[[1]][order(upEnriched_go[["GO_Molecular_Function_2023"]][["Adjusted.P.value"]]),] %>% dplyr::filter(upEnriched_go[["GO_Molecular_Function_2023"]][["Adjusted.P.value"]] < "0.25")


plotlist <- list()


#GOdatabase <- paste0(gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_go)[l]))))


#GO-Terms
GOdatabase <- paste0("Cluster ",n," - ",gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_go)[1]))))
goup  <- try(plotEnrich(upEnriched_go[[1]], showTerms = termnumber, numChar = 50, y = "Count", orderBy = "P.Value", title = upgroup)+ theme(text=element_text(size=14)))
godown <- try(plotEnrich(dnEnriched_go[[1]], showTerms = termnumber, numChar = 50, y = "Count", orderBy = "P.Value", title = downgroup)+ theme(text=element_text(size=14))) 
try(write.xlsx(upEnriched_go[[1]], file = paste0(group,"_",GOdatabase," - ",upgroup,"Cluster",n,".xlsx")))
try(write.xlsx(dnEnriched_go[[1]], file = paste0(group,"_",GOdatabase," - ",downgroup,"Cluster",n,".xlsx")))
plot <- ggarrange(goup,godown, labels =c("A","B"))
p1 <-annotate_figure(plot, top = text_grob(GOdatabase, color = "red", face = "bold", size = 17)) 
p1
GOdatabase <- paste0("Cluster ",n," - ",gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_go)[2]))))
goup  <- try(plotEnrich(upEnriched_go[[2]], showTerms = termnumber, numChar = 50, y = "Count", orderBy = "P.Value", title = upgroup)+ theme(text=element_text(size=14)))
godown <- try(plotEnrich(dnEnriched_go[[2]], showTerms = termnumber, numChar = 50, y = "Count", orderBy = "P.Value", title = downgroup)+ theme(text=element_text(size=14)))  
try(write.xlsx(upEnriched_go[[2]], file = paste0(group,"_",GOdatabase," - ",upgroup,"Cluster",n,".xlsx")))
try(write.xlsx(dnEnriched_go[[2]], file = paste0(group,"_",GOdatabase," - ",downgroup,"Cluster",n,".xlsx")))


plot <- ggarrange(goup,godown, labels =c("A","B"))

plot <- ggarrange(goup,godown, labels =c("A","B"))
p2<- annotate_figure(plot, top = text_grob(GOdatabase, color = "red", face = "bold", size = 17)) 

GOdatabase <- paste0("Cluster ",n," - ",gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_go)[3]))))
goup  <- try(plotEnrich(upEnriched_go[[3]], showTerms = termnumber, numChar = 50, y = "Count", orderBy = "P.Value", title = upgroup)+ theme(text=element_text(size=14)))
godown <- try(plotEnrich(dnEnriched_go[[3]], showTerms = termnumber, numChar = 50, y = "Count", orderBy = "P.Value", title = downgroup)+ theme(text=element_text(size=14)))  
try(write.xlsx(upEnriched_go[[3]], file = paste0(group,"_",GOdatabase," - ",upgroup,"Cluster",n,".xlsx")))
try(write.xlsx(dnEnriched_go[[3]], file = paste0(group,"_",GOdatabase," - ",downgroup,"Cluster",n,".xlsx")))

plot <- ggarrange(goup,godown, labels =c("A","B"))
p3<- annotate_figure(plot, top = text_grob(GOdatabase, color = "red", face = "bold", size = 17)) 




#Pathwayanalysis

Pathwaydatabase <- paste0("Cluster ",n," - ",gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_pw)[1]))))
pwup  <- try(plotEnrich(upEnriched_pw[[1]], showTerms = termnumber, numChar = 50, y = "Count", orderBy = "P.Value", title = upgroup)+ theme(text=element_text(size=14)))
pwdown <- try(plotEnrich(dnEnriched_pw[[1]], showTerms = termnumber, numChar = 50, y = "Count", orderBy = "P.Value", title = downgroup)+ theme(text=element_text(size=14))) 
try(write.xlsx(upEnriched_pw[[1]], file = paste0(group,"_",Pathwaydatabase," - ",upgroup,"Cluster",n,".xlsx")))
try(write.xlsx(dnEnriched_pw[[1]], file = paste0(group,"_",Pathwaydatabase," - ",downgroup,"Cluster",n,".xlsx")))

plot <- ggarrange(pwup,pwdown, labels =c("A","B"))
p4 <-annotate_figure(plot, top = text_grob(Pathwaydatabase, color = "red", face = "bold", size = 17)) 

Pathwaydatabase <- paste0("Cluster ",n," - ",gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_pw)[2]))))
pwup  <- try(plotEnrich(upEnriched_pw[[2]], showTerms = termnumber, numChar = 50, y = "Count", orderBy = "P.Value", title = upgroup)+ theme(text=element_text(size=14)))
pwdown <- try(plotEnrich(dnEnriched_pw[[2]], showTerms = termnumber, numChar = 50, y = "Count", orderBy = "P.Value", title = downgroup)+ theme(text=element_text(size=14)))  
try(write.xlsx(upEnriched_pw[[2]], file = paste0(group,"_",Pathwaydatabase," - ",upgroup,"Cluster",n,".xlsx")))
try(write.xlsx(dnEnriched_pw[[2]], file = paste0(group,"_",Pathwaydatabase," - ",downgroup,"Cluster",n,".xlsx")))

plot <- ggarrange(pwup,pwdown, labels =c("A","B"))
p5 <-annotate_figure(plot, top = text_grob(Pathwaydatabase, color = "red", face = "bold", size = 17)) 

Pathwaydatabase <- paste0("Cluster ",n," - ",gsub("_"," ",(gsub('_[[:digit:]]+', '', names(upEnriched_dd)[1]))))
ddup  <- try(plotEnrich(upEnriched_dd[[1]], showTerms = termnumber, numChar = 50, y = "Count", orderBy = "P.Value", title = upgroup)+ theme(text=element_text(size=14))) 
dddown <- try(plotEnrich(dnEnriched_dd[[1]], showTerms = termnumber, numChar = 50, y = "Count", orderBy = "P.Value", title = downgroup)+ theme(text=element_text(size=14)))  
try(write.xlsx(upEnriched_dd[[1]], file = paste0(group,"_",Pathwaydatabase," - ",upgroup,"Cluster",n,".xlsx")))
try(write.xlsx(dnEnriched_dd[[1]], file = paste0(group,"_",Pathwaydatabase," - ",downgroup,"Cluster",n,".xlsx")))
 
plot <- ggarrange(ddup,dddown, labels =c("A","B"))
p6 <-annotate_figure(plot, top = text_grob(Pathwaydatabase, color = "red", face = "bold", size = 17)) 



plotlist <- list(p1,p2,p3,p4,p5,p6)

ggexport(plotlist = plotlist, filename = paste0(group,"Cluster ",n," - GO and Pathway Analysis - Enriched in ",upgroup,"Cluster",n,".pdf"), width=20, height=11, res=300)

rm(up,dn,upk,dnk,upw,dnw)
# Example usage
up <- dotplot_function(upEnriched_go[["GO_Biological_Process_2023"]]%>% dplyr::filter(!grepl("cancer", Term, ignore.case = TRUE)))+ ggtitle(paste0(" Enriched in Upregulated Genes"))
dn <- dotplot_function(dnEnriched_go[["GO_Biological_Process_2023"]]%>% dplyr::filter(!grepl("cancer", Term, ignore.case = TRUE)))+ ggtitle(paste0(" Enriched in Downregulated Genes"))

#together <- plot_grid(up, dn, ncol= 2)
if (test_ggplot_error(up) && test_ggplot_error(dn)) {
  # Both up and dn have values, so plot them with ncol = 2
  together <- plot_grid(up, dn, ncol = 2)
  
} else if (test_ggplot_error(up)) {
  # Only up has a value, so plot it alone
  together <- plot_grid(up,NULL, ncol = 2)
} else if (test_ggplot_error(dn)) {
  # Only dn has a value, so plot it alone
  together <- plot_grid(NULL,dn, ncol = 2)
} else {
  # Neither up nor dn have values, so do not create a plot
  together <- NULL
}

if (is.null(together)) {
  } else {
    together <- try(together + plot_annotation(title = paste0("GO-BP ",group)) & theme(plot.title = element_text(hjust = 0.5, face = "bold")))
    ggsave(together, file=paste0(group," DotPlot_Go_Biological_Process Cluster", n,".pdf"), units="in", width=13, height=8, dpi=600,device="pdf")
  
}



upk <- dotplot_function(upEnriched_pw[["KEGG_2019_Mouse"]]%>% dplyr::filter(!grepl("cancer", Term, ignore.case = TRUE)))+ ggtitle(paste0(" Enriched in Upregulated Genes"))
dnk <- dotplot_function(dnEnriched_pw[["KEGG_2019_Mouse"]]%>% dplyr::filter(!grepl("cancer", Term, ignore.case = TRUE)))+ ggtitle(paste0(" Enriched in Downregulated Genes"))

if (test_ggplot_error(upk) && test_ggplot_error(dnk)) {
  # Both up and dn have values, so plot them with ncol = 2
  togetherk <- plot_grid(upk, dnk, ncol = 2)
} else if (test_ggplot_error(upk)) {
  # Only up has a value, so plot it alone
  togetherk <- plot_grid(upk,NULL, ncol = 2)
} else if (test_ggplot_error(dnk)) {
  # Only dn has a value, so plot it alone
  togetherk <- plot_grid(NULL,dnk, ncol = 2)
} else {
  # Neither up nor dn have values, so do not create a plot
  togetherk <- NULL
}

if (is.null(togetherk)) {
} else {
  togetherk <- try(togetherk + plot_annotation(title = paste0("KEGG ", group)) & theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  ggsave(togetherk, file=paste0(group," DotPlot_KEGG Cluster", n,".pdf"), units="in", width=13, height=8, dpi=600,device="pdf")
   
}



upw <- dotplot_function(upEnriched_pw[["WikiPathways_2019_Mouse"]] %>% dplyr::filter(!grepl("cancer", Term, ignore.case = TRUE)))+ ggtitle(paste0(" Enriched in Upregulated Genes"))
dnw <- dotplot_function(dnEnriched_pw[["WikiPathways_2019_Mouse"]] %>% dplyr::filter(!grepl("cancer", Term, ignore.case = TRUE)))+ ggtitle(paste0(" Enriched in Downregulated Genes"))

if (test_ggplot_error(upw) && test_ggplot_error(dnw)) {
  # Both up and dn have values, so plot them with ncol = 2
  togetherkw <- plot_grid(upw, dnw, ncol = 2)
} else if (test_ggplot_error(upw)) {
  # Only up has a value, so plot it alone
  togetherkw <- plot_grid(upw,NULL, ncol = 2)
} else if (test_ggplot_error(dnw)) {
  # Only dn has a value, so plot it alone
  togetherkw <- plot_grid(NULL,dnw, ncol = 2)
} else {
  # Neither up nor dn have values, so do not create a plot
  togetherkw <- NULL
}

if (is.null(togetherkw)) {
} else {
  togetherkw <- try(togetherkw + plot_annotation(title = paste0("Wikipathway ", group)) & theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  ggsave(togetherkw, file=paste0(group," DotPlot_WikiP Cluster", n,".pdf"), units="in", width=13, height=8, dpi=600,device="pdf")
  
}



upcc <- dotplot_function(upEnriched_go[["GO_Cellular_Component_2023"]] %>% dplyr::filter(!grepl("cancer", Term, ignore.case = TRUE)))+ ggtitle(paste0(" Enriched in Upregulated Genes"))
dncc <- dotplot_function(dnEnriched_go[["GO_Cellular_Component_2023"]] %>% dplyr::filter(!grepl("cancer", Term, ignore.case = TRUE)))+ ggtitle(paste0(" Enriched in Downregulated Genes"))

if (test_ggplot_error(upcc) && test_ggplot_error(dncc)) {
  # Both up and dn have values, so plot them with ncol = 2
  togethercc <- plot_grid(upcc, dncc, ncol = 2)
} else if (test_ggplot_error(upcc)) {
  # Only up has a value, so plot it alone
  togethercc <- plot_grid(upcc,NULL, ncol = 2)
} else if (test_ggplot_error(dncc)) {
  # Only dn has a value, so plot it alone
  togethercc <- plot_grid(NULL,dncc, ncol = 2)
} else {
  # Neither up nor dn have values, so do not create a plot
  togethercc <- NULL
}

if (is.null(togethercc)) {
} else {
  togethercc <- try(togethercc + plot_annotation(title = paste0("GO Cellular Components ", group)) & theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  ggsave(togethercc, file=paste0(group," GO Cellular Components - Cluster", n,".pdf"), units="in", width=13, height=8, dpi=600,device="pdf")
  
}

upk <- dotplot_function(upEnriched_pw[["KEGG_2019_Mouse"]]%>% dplyr::filter(!grepl("cancer", Term, ignore.case = TRUE)))+ ggtitle(paste0(" Enriched in KEGG"))
upw <- dotplot_function(upEnriched_pw[["WikiPathways_2019_Mouse"]] %>% dplyr::filter(!grepl("cancer", Term, ignore.case = TRUE)))+ ggtitle(paste0(" Enriched in WikiPathway"))

if (test_ggplot_error(upw) && test_ggplot_error(dnw)) {
  # Both up and dn have values, so plot them with ncol = 2
  togetherkw <- plot_grid(upk, upw, ncol = 2)
} else if (test_ggplot_error(upk)) {
  # Only up has a value, so plot it alone
  togetherkw <- plot_grid(upk,NULL, ncol = 2)
} else if (test_ggplot_error(upw)) {
  # Only dn has a value, so plot it alone
  togetherkw <- plot_grid(NULL,upw, ncol = 2)
} else {
  # Neither up nor dn have values, so do not create a plot
  togetherkw <- NULL
}

if (is.null(togetherkw)) {
} else {
  togetherkw <- try(togetherkw + plot_annotation(title = paste0(group)) & theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  ggsave(togetherkw, file=paste0(group," UP_Combined_DotPlot_WikiandKEGG Cluster", n,".pdf"), units="in", width=13, height=8, dpi=600,device="pdf")
}




uppa <- dotplot_function(upEnriched_pa[["Panther_2016"]]%>% dplyr::filter(!grepl("cancer", Term, ignore.case = TRUE)))+ ggtitle(paste0(" Enriched in Upregulated Genes"))
dnpa <- dotplot_function(dnEnriched_pa[["Panther_2016"]] %>% dplyr::filter(!grepl("cancer", Term, ignore.case = TRUE)))+ ggtitle(paste0(" Enriched in Downregulated Genes"))

if (test_ggplot_error(dnpa) && test_ggplot_error(uppa)) {
  # Both up and dn have values, so plot them with ncol = 2
  togetherpa <- plot_grid(uppa, dnpa, ncol = 2)
} else if (test_ggplot_error(uppa)) {
  # Only up has a value, so plot it alone
  togetherpa <- plot_grid(uppa,NULL, ncol = 2)
} else if (test_ggplot_error(dnpa)) {
  # Only dn has a value, so plot it alone
  togetherpa <- plot_grid(NULL,dnpa, ncol = 2)
} else {
  # Neither up nor dn have values, so do not create a plot
  togetherpa <- NULL
}

if (is.null(togetherpa)) {
} else {
  togetherpa <- try(togetherpa + plot_annotation(title = paste0(group)) & theme(plot.title = element_text(hjust = 0.5, face = "bold")))
  ggsave(togetherpa, file=paste0(group," Panther2016 Cluster", n,".pdf"), units="in", width=13, height=8, dpi=600,device="pdf")
}

#together <- try(together + plot_annotation(title = paste0("KEGG ",group)) & theme(plot.title = element_text(hjust = 0.5, face = "bold")))


#d1 <- plotEnrich(upEnriched_dd[[1]], showTerms = 10, numChar = 30, y = "Count", orderBy = "Combined.Score", title = upgroup)
#d2 <- plotEnrich(dnEnriched_dd[[1]], showTerms = 10, numChar = 30, y = "Count", orderBy = "Combined.Score", title = downgroup)
#d1+d2


} else {
  
print("No Down Genes")
  
}
  
  
} else { print("No Up Genes")
  
}

}



writeLines(capture.output(sessionInfo()), paste0(here("RNA","DGE",paste0("EnrichGo_Reactome_09292023Comparison.txt"))))





##########

# 
# 
# dotplot_function_old <- function(data) {
#   df <- data
#   
#   
#   df$Term <- gsub(" \\(.*\\)", "", df$Term)
#   
#   df$Term <- substr(df$Term, 1, 35)
#   df$Term[duplicated(df$Term)] <- paste0(df$Term[duplicated(df$Term)], " ")
#   df$Term[duplicated(df$Term)] <- paste0(df$Term[duplicated(df$Term)], "  ")
#   
#   
#   df <- df %>% dplyr::dplyr::filter(Adjusted.P.value <0.05) 
#   
#   df <- df %>% slice_min(n = 10, order_by = Adjusted.P.value)
#   print(df$Term)
#   
#   Significant <- as.numeric(sub("/\\d+$", "", as.character(df$Overlap)))
#   
#   df$Generatio <- sapply(df$Overlap, function(x) eval(parse(text = x)))
#   
#   df$Term <- factor(df$Term, levels = rev(df$Term[order(df$Adjusted.P.value)]))
#   title <- paste0(group," Cluster ",n)
#   
#   
#   ggplot(df, aes_string(x="Adjusted.P.value", y="Term", size="Significant", color="Generatio")) +
#     geom_point() +
#     scale_color_continuous(low="red", high="blue", name = "Generatio",
#                            guide=guide_colorbar(reverse=TRUE)) +
#     ylab(NULL) + ggtitle(title) + theme_dose(12) +
#     scale_size(range=c(3, 8), name="Gene Count") +
#     guides(size  = guide_legend(order = 1),
#            color = guide_colorbar(order = 2)) +
#     theme(axis.text.y = element_text(size = 12, face = "bold")) +
#     scale_x_reverse() + 
#     xlab("adj. P value") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
#   
# }
