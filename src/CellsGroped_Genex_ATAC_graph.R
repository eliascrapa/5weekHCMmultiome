#Make Venn Diagramm
library(openxlsx)
library(GOplot)
library(ggplot2)
library(tidyverse)

#saveRDS(celltypelevels, "celltypelevels.Rds")

#setwd("C:/Users/tthottakara/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/NewDGE/")
setwd("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/NewDGE/")


celltypelevels <- readRDS("celltypelevels.Rds")
data.dir <- './TNTvsMHC'
dir.create(data.dir)
setwd(data.dir)
#VEnnDiagramm MHC and TNT
my_list <-list()
for (n in celltypelevels){
  
  skip_to_next <- FALSE
  
  
  
    
    tnt_input <- paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/NewDGE/NewDGEDiffexpr_",n,"_x_TNT_vs_",n,"_x_CTR.xlsx")
    mhc_input <- paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/NewDGE/NewDGEDiffexpr_",n,"_x_MHC_vs_",n,"_x_CTR.xlsx")
    #filem <-  paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsCTR_C",n,".rnk")  
    title <- paste0("Venn Diagram Cluster ",n)
    vennoutput <- paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/NewDGE/Venn/Venndiagram_TNTvsMHC_Cluster_",n,".pdf")
    
  # tnt_input <- paste0("C:/Users/tthottakara/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/NewDGE/NewDGEDiffexpr_",n,"_x_TNT_vs_",n,"_x_CTR.xlsx")
  # mhc_input <- paste0("C:/Users/tthottakara/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/NewDGE/NewDGEDiffexpr_",n,"_x_MHC_vs_",n,"_x_CTR.xlsx")
  # #filem <-  paste0("~/OneDrive - University of California, San Francisco/DATA/SingleCellSeq/Analysis/All/DGE/TNTvsCTR_C",n,".rnk")  
  # title <- paste0("Venn Diagram Cluster ",n)
  # vennoutput <- paste0("C:/Users/tthottakara/OneDrive - UCSF/DATA/SingleCellSeq/Analysis/All/NewDGE/Venn/Venndiagram_TNTvsMHC_Cluster_",n,".pdf")
  # 
  TNT <- read.xlsx(tnt_input, 1)
  MHC <- read.xlsx(mhc_input, 1)
  
  sapply(TNT, class) 
  sapply(MHC, class) 
  #TNT[,3:6] <- sapply(TNT[,3:6],as.numeric)
  #MHC[,3:6] <- sapply(MHC[,3:6],as.numeric)
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  #TNTvsCTR_C8 <- transform(TNTvsCTR_C8, avg_log2FC = as.numeric(avg_log2FC))
  
  TNT <- dplyr::filter(TNT, p_val_adj < 0.05)
  MHC <- dplyr::filter(MHC, p_val_adj < 0.05)
  TNT_sel <-dplyr::select(TNT, names, avg_log2FC)
  MHC_sel <-dplyr::select(MHC, names, avg_log2FC)
 
  
  df <- cbind(length(TNT_sel$avg_log2FC[TNT_sel$avg_log2FC>=0]),length(TNT_sel$avg_log2FC[TNT_sel$avg_log2FC<0]))
  df <- rbind(df,(cbind(length(MHC_sel$avg_log2FC[MHC_sel$avg_log2FC>=0]),length(MHC_sel$avg_log2FC[MHC_sel$avg_log2FC<0]))))
  df <- cbind(c("TNT","MHC"),n,df)
  colnames(df) <- c("Genotype", "Cluster", "Upregulated", "Downregulated")
  
  
  colnames(df) <- c("Genotype", "Cluster", "Upregulated", "Downregulated")
  
  my_list[[length(my_list) + 1]] <- df
  
}

#saveRDS(my_list, file="GeneexpressionPerGEnotypeacrossCelltypes.Rds")

df2 <- as.data.frame(do.call(rbind, (my_list)))
df2$Upregulated <- as.numeric(df2$Upregulated)
df2$Downregulated <- as.numeric(df2$Downregulated)
df2$Cluster <- as.factor(df2$Cluster)


df2$Genotype <- factor(df2$Genotype, levels = c("TNT","MHC"))
df2$Cluster <- factor(df2$Cluster, levels =  rev(c("Cardiomyocyte","Fibroblast","EndothelialCell","Leukocyte","VSMC","EpicardialCell","undeterm")))
df2$DownregulatedNeg <- -df2$Downregulated

df2_long <- df2 %>% gather(key = "Regulation", value = "Value", Upregulated, DownregulatedNeg)

gex <- ggplot(df2_long, aes(fill = Genotype, y = Value, x = Cluster)) +
  geom_bar(position = position_dodge(width = 0.5), stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("#CC9900", "#006666", "#006699"), labels=c("TnT","MyHC")) +
  xlab("Grouped Celltypes") +
  ylab("Number of Up- and Downregulated Genes") +
  ggtitle(paste0("Differential Geneexpression per Celltype")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill = "Groups") +
  theme(legend.title.align = 0.5, legend.text = element_text(size = 8)) +
  theme(legend.position = c(0.86, 0.15))+
  geom_hline(yintercept = 0, color = "black")+
  theme(text=element_text(size=20),axis.text.y=element_text(size=30))

ggsave(filename = paste0("DiffGenexpressionPerCelltype.pdf"), gex, device = "pdf", width = 6, height = 5)

file <- paste0("DiffGenexpressionPerCelltype.pdf")
system2('open', args = c('-a Preview.app', file), wait = FALSE)





####################################################

ATAC <- read.xlsx("C:\\Users\\tthottakara\\OneDrive - UCSF\\DATA\\SingleCellSeq\\Analysis\\All\\NewDGE\\ATAC\\CelltypesPerSample_DiffFeatures.xlsx",1)
colnames(ATAC) <- c("Number","Group","Upregulated","Downregulated","Celltype","Genotype")
#ATAC <- ATAC %>% mutate_at('Number','Upregulated','Downregulated', as.numeric) %>% str()

ATAC  <-  ATAC %>% filter(Genotype =="TNT"|Genotype =="MHC")

ATAC$Genotype <- factor(ATAC$Genotype, levels = c("TNT","MHC"))
ATAC$Celltype <- factor(ATAC$Celltype, levels =  rev(c("Cardiomyocyte","Fibroblast","EndothelialCell","Leukocyte","VSMC","EpicardialCell","undeterm")))
ATAC$DownregulatedNeg <- -ATAC$Downregulated

ATAC_long <- ATAC %>% gather(key = "Regulation", value = "Value", Upregulated, DownregulatedNeg)

max_val <- max(abs(ATAC_long$Value))

ATACplot <-ggplot(ATAC_long, aes(fill = Genotype, y = Value, x = Celltype)) +
  geom_bar(position = position_dodge(width = 0.5), stat = "identity") +
  coord_flip(ylim = c(-max_val, max_val)) +
  scale_fill_manual(values = c("#CC9900", "#006666", "#006699"), labels=c("TnT","MyHC")) +
  xlab("Grouped Celltypes") +
  ylab("Number of Up- and Downregulated Features") +
  ggtitle(paste0("Differential Feature Accessibility per Celltype")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill = "Groups") +
  theme(legend.title.align = 0.5, legend.text = element_text(size = 8)) +
  theme(legend.position = c(0.86, 0.15))+
  geom_hline(yintercept = 0, color = "black")+
  scale_y_continuous(breaks = c(-60000,-40000,-20000,0,20000,40000,60000))+
  theme(text=element_text(size=20),axis.text.y=element_text(size=30))



ggsave(filename = paste0("ATACPerCelltype.pdf"), ATACplot, device = "pdf", width = 6, height = 5)
file <- paste0("ATACPerCelltype.pdf")
system2('open', args = c('-a Preview.app', file), wait = FALSE)







ggplot(ATAC, aes(Cluster), ylim(-95000, 95000)) + 
  geom_bar(data = ATAC, 
           aes(y = Upregulated, fill = Genotype), stat = "identity", position = "dodge") +
  geom_bar(data = ATAC, 
           aes(y = -Downregulated, fill = Genotype), stat = "identity", position = "dodge") + 
  geom_hline(yintercept = 0,colour = "grey90")+
  theme_minimal()+
  scale_fill_brewer(palette="Paired")


gg2 <- last_plot() + 
  geom_text(data = ATAC, 
            aes(Cluster, Upregulated, group=Genotype, label=Upregulated),
            position = position_dodge(width=0.9), vjust = -1.3, size=2.3) +
  geom_text(data = ATAC, 
            aes(Cluster, -Downregulated, group=Genotype, label=Downregulated),
            position = position_dodge(width=0.9), vjust = +2.2, size=2.3) +
  coord_cartesian(ylim = c(-95000, 50000))


gg2
ggsave(file="ATAC-PLot.pdf",gg2, width=12, height=9, units="in")


Cluster






  
 # df2 %>%
  ggplot(df2, aes(fill=Genotype, y=Upregulated, x= Cluster)) + 
  geom_bar(position=position_dodge(width = 0.5), stat="identity") +
  coord_flip() +
  #scale_fill_brewer(palette="Set1") +
  #scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99"))+ 
  scale_fill_manual(values=c("#CC9900","#006666", "#006699"))+ 
  xlab("Upregulated Genes") + ylab("Relative Regulon Activity") +
  ggtitle(paste0("Differential Genexpression")) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(fill = "Groups") + theme(legend.title.align=0.5, legend.text = element_text(size=8)) + #, legend.key.width = unit(1, 'in'))  
  theme(legend.position = c(0.9, 0.15))



p <- ggplot(df2, aes(fill=Genotype, y=Upregulated, x= Cluster)) + geom_bar(position=position_dodge(width = 0.5), stat="identity")
print(p)




CXregplot <- CX_choRegulators %>%
  mutate(CellType = fct_relevel(CellType, paste0("CTR C",n," - ",name), paste0("TnT C",n," - ",name), paste0("MyHC C",n," - ",name))) %>%
  ggplot( aes(fill=CellType, y=RelativeActivity, x= reorder(forcats::fct_rev(Regulon), RelativeActivity))) + 
  geom_bar(position=position_dodge(width = 0.5), stat="identity") +
  coord_flip() +
  #scale_fill_brewer(palette="Set1") +
  #scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99"))+ 
  scale_fill_manual(values=c("#CC9900", "#006666", "#006699"))+ 
  xlab("Active Regulons") + ylab("Relative Regulon Activity") +
  ggtitle(paste0("Regulon Activity Comparison")) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(fill = "Groups") + theme(legend.title.align=0.5, legend.text = element_text(size=8)) + #, legend.key.width = unit(1, 'in'))  
  theme(legend.position = c(0.86, 0.15))
ggsave(filename = paste0(selectedClustering,"_C",n,"_TF-",paste0(chosenTFs, collapse ="-"),"_overlapping.pdf"), CXregplot, device = "pdf", width = 6, height = 5)
file <- paste0(selectedClustering,"_C",n,"_TF-",paste0(chosenTFs, collapse ="-"),"_overlapping.pdf")
system2('open', args = c('-a Preview.app', file), wait = FALSE)
