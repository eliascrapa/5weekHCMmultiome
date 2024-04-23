library(ggplot2)
library(gridExtra)
library(here)
library(openxlsx)
library("readxl")
library(tidyverse)
library(stringr)
library(cowplot)
library(ggh4x)
library(GOplot)
library(pheatmap)
library("org.Mm.eg.db") 
library(viridis)
library(grid)


setwd("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/MP/")
setwd("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/Fib/")
setwd("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/EC/")
setwd("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/EC-C16/")
setwd("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/")
setwd("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/Fib-C21/")
setwd("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/CM-C8/")
setwd("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/VSMC/")

ipanumber <- c("2","6","9","14","16","20","22")
wrap_text <- function(x, width = 20) {str_wrap(x, width)}

save_pheatmap_pdf <- function(x, filename, width=4, height=12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


for (n in ipanumber){
q <- 1
cluster <- "20"
#pathwaysource <- read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/DGE/GenesOfPathways/Pathways.xlsx",q)
ipadataTNT <- read_excel(paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/IPA/CanonicalPathway/TNT_C",cluster,"_CanonicalPathway.xls"),1, skip=1)
ipadataMHC <- read_excel(paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/IPA/CanonicalPathway/MHC_C",cluster,"_CanonicalPathway.xls"),1, skip=1)

diffgeneTNT <-  openxlsx::read.xlsx(paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_TNT_vs_CONTROLS_C",cluster,".xlsx"),q)
diffgeneMHC <-  openxlsx::read.xlsx(paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/DGE/RNA_ClusterMarker_MHC_vs_CONTROLS_C",cluster,".xlsx"),q)

regulons <- read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/Regulons.xlsx",q)
allregulons <- read.table("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/Step2_regulonTargetsInfo.tsv")

kegg <- openxlsx::read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/KEGG_2019_Mouse.xlsx",q, colNames = FALSE)
reactome <- openxlsx::read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/Reactome_2022.xlsx",q, colNames = FALSE)
wikipathways <- openxlsx::read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/WikiPathways_2019_Mouse.xlsx",q, colNames = FALSE)
PANTHERpathways <- openxlsx::read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/EC/DB/PANTHER_Pathways_2014Dec17_gene_set_library_crisp.xlsx",q, colNames = FALSE)
panther2016 <- openxlsx::read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/EC/DB/Panther_2016_.xlsx",q, colNames = FALSE)

amigo <- read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/AMIGO/CardiacHypertrophyAmigo.xlsx", colNames = FALSE)
amigo <- amigo %>% dplyr::select(3,5,10,12,13,14,15)
colnames(amigo) <- c("gene","spGO-term","fullname","datasource","organism","GO-category","GO-Code")


amigo_TGF <- openxlsx::read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/Fib/DB/AmiGO-TGFbetaPathways.xlsx",q, colNames = FALSE)
amigo_TGF <- amigo_TGF %>% dplyr::select(3,5,10,12,13,14,15)
colnames(amigo_TGF) <- c("gene","spGO-term","fullname","datasource","organism","GO-category","GO-Code")


amigo_EC <- openxlsx::read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/EC/DB/EC_GO_list.xlsx",q, colNames = FALSE)
amigo_EC <- amigo_EC %>% dplyr::select(3,5,10,12,13,14,15)
colnames(amigo_EC) <- c("gene","spGO-term","fullname","datasource","organism","GO-category","GO-Code")

amigo_PDGF <- openxlsx::read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/VSMC/PDGF_EC__VOrlageGO_list.xlsx",q, colNames = FALSE)
amigo_PDGF <- amigo_PDGF %>% dplyr::select(3,5,10,12,13,14,15)
colnames(amigo_PDGF) <- c("gene","spGO-term","fullname","datasource","organism","GO-category","GO-Code")


ILTNF <- openxlsx::read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/MP/DB/ILTNFgenes.xlsx",q, colNames = TRUE)

amigo_phag <- openxlsx::read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/MP/DB/phagosome_maturation-GO-0090382.xlsx",q, colNames = FALSE)
amigo_phag <- amigo_phag %>% dplyr::select(3,9,11,12,13,14)
colnames(amigo_phag) <- c("gene","fullname","datasource","organism","GO-category","GO-Code")

  
#AMIGO Hypertrophy Genelist
hypertrophy <- amigo %>% filter(organism =="NCBITaxon:10116"|organism =="NCBITaxon:10090"|organism =="NCBITaxon:9606") # ratus norvegicus, mus musculus, homo sapiens
# Assuming your dataframe is called hypertrophy and the column is gene
#hypertrophy2 <- distinct(hypertrophy, tolower(gene), .keep_all = TRUE) %>% filter(datasource != "miRNA" )
hypertrophy <-  hypertrophy %>% dplyr::filter(datasource != "miRNA" ) %>% dplyr::select(6,1) %>% group_by(`GO-category`)
hypertrophy <-hypertrophy %>% group_by(`GO-category`) %>% mutate(row = row_number()) %>% spread(key = row, value = gene) %>% ungroup()
hypertrophy <- as.data.frame(t(hypertrophy))
colnames(hypertrophy) <- hypertrophy[1,]
hypertrophy <- hypertrophy[-1,]


#AMIGO Hypertrophy Genelist
PDGF  <- amigo_PDGF  %>% filter(organism =="NCBITaxon:10116"|organism =="NCBITaxon:10090"|organism =="NCBITaxon:9606") # ratus norvegicus, mus musculus, homo sapiens
# Assuming your dataframe is called hypertrophy and the column is gene
#hypertrophy2 <- distinct(hypertrophy, tolower(gene), .keep_all = TRUE) %>% filter(datasource != "miRNA" )
PDGF  <-  PDGF  %>% dplyr::filter(datasource != "miRNA" ) %>% dplyr::select(6,1) %>% group_by(`GO-category`)
PDGF  <-PDGF  %>% group_by(`GO-category`) %>% mutate(row = row_number()) %>% spread(key = row, value = gene) %>% ungroup()
PDGF  <- as.data.frame(t(PDGF ))
colnames(PDGF ) <- PDGF [1,]
PDGF  <- PDGF [-1,]

PDGF[,2:3]

write.xlsx(PDGF[,2:3],"PDGF-Genes.xlsx")
PDGFGenes <- read.xlsx("PDGF-Genes.xlsx")[,1]

#AMIGO Hypertrophy Genelist
EC <- amigo_EC %>% filter(organism =="NCBITaxon:10116"|organism =="NCBITaxon:10090"|organism =="NCBITaxon:9606") # ratus norvegicus, mus musculus, homo sapiens
# Assuming your dataframe is called hypertrophy and the column is gene
#hypertrophy2 <- distinct(hypertrophy, tolower(gene), .keep_all = TRUE) %>% filter(datasource != "miRNA" )
EC <-  EC %>% dplyr::filter(datasource != "miRNA" ) %>% dplyr::select(6,1) %>% group_by(`GO-category`)
EC <-EC %>% group_by(`GO-category`) %>% mutate(row = row_number()) %>% spread(key = row, value = gene) %>% ungroup()
EC <- as.data.frame(t(EC))
colnames(EC) <- EC[1,]
EC <- EC[-1,]



#AMIGO Hypertrophy Genelist
TGFb <- amigo_TGF %>% dplyr::filter(organism =="NCBITaxon:10116"|organism =="NCBITaxon:10090"|organism =="NCBITaxon:9606") # ratus norvegicus, mus musculus, homo sapiens
# Assuming your dataframe is called hypertrophy and the column is gene
TGFb <-  TGFb %>% dplyr::filter(datasource != "miRNA" ) %>% dplyr::select(6,1) %>% group_by(`GO-category`)
TGFb <-TGFb %>% group_by(`GO-category`) %>% mutate(row = row_number()) %>% spread(key = row, value = gene) %>% ungroup()
TGFb <- as.data.frame(t(TGFb))
colnames(TGFb) <- TGFb[1,]
TGFb <- TGFb[-1,]

phag <- amigo_phag %>% dplyr::filter(organism =="NCBITaxon:10116"|organism =="NCBITaxon:10090"|organism =="NCBITaxon:9606") # ratus norvegicus, mus musculus, homo sapiens
# Assuming your dataframe is called hypertrophy and the column is gene
phag <-  phag %>% dplyr::filter(datasource != "miRNA" ) %>% dplyr::select(5,1) %>% group_by(`GO-category`)
phag <-phag %>% group_by(`GO-category`) %>% mutate(row = row_number()) %>% spread(key = row, value = gene) %>% ungroup()
phag <- as.data.frame(t(phag))
colnames(phag) <- phag[1,1]
phag <- phag[-1,]


#######IPA-DATA
#Extract Molecules and put them in new columns
ipadataTNT$Molecules <- strsplit(ipadataTNT$Molecules, ",")
ipadataTNT <- ipadataTNT %>% unnest_wider(Molecules, names_sep = ",")
ipadataMHC$Molecules <- strsplit(ipadataMHC$Molecules, ",")
ipadataMHC <- ipadataMHC %>% unnest_wider(Molecules, names_sep = ",")
#Summarize Genes that define a pathway
GenesforPathway <- left_join(ipadataTNT %>% dplyr::select (-c(2,3,4)),ipadataMHC %>% dplyr::select (-c(2,3,4)),by="Ingenuity Canonical Pathways")
GenesforPathway <- t(GenesforPathway)
colnames(GenesforPathway) <- GenesforPathway[1,]
GenesforPathway <- as.data.frame(GenesforPathway[2:nrow(GenesforPathway),])
ipa <- GenesforPathway
rm(GenesforPathway)


######KEGG & REACTOME
kegg <- t(kegg)
colnames(kegg) <- kegg[1,]
kegg <- kegg[-1,]
kegg <- as.data.frame(kegg)
reactome <- t(reactome)
colnames(reactome) <- reactome[1,]
reactome <- reactome[-1,]
reactome <- as.data.frame(reactome)
wikipathways <- t(wikipathways)
colnames(wikipathways) <- wikipathways[1,]
wikipathways <- wikipathways[-1,]
wikipathways <- as.data.frame(wikipathways)
PANTHERpathways <- t(PANTHERpathways)
colnames(PANTHERpathways) <- PANTHERpathways[1,]
PANTHERpathways <- PANTHERpathways[-1,]
PANTHERpathways <- as.data.frame(PANTHERpathways)

panther2016 <- t(panther2016)
colnames(panther2016) <- panther2016[1,]
panther2016 <- panther2016[-1,]
panther2016 <- as.data.frame(panther2016)



panther2016
#######REGULONS
colnames(allregulons) <- allregulons[1,]
allregulons <- allregulons[-1,]

# Group the allregulons dataframe by 'TF' column
grouped <- split(allregulons, allregulons$TF)
#grouped <- split(regulons, regulons$TF)
######MAke Regulon TF and Genelist
# Split the dataframe based on 'TF'
split_regulon_genes <- split(allregulons, allregulons$TF)
# Extract 'gene' column from each dataframe
gene_list <- lapply(split_regulon_genes, function(regulon_genes) regulon_genes$gene)
# Convert the list to a data frame
regulon_genes <- do.call(rbind, lapply(gene_list, `length<-`, max(lengths(gene_list))))
# Set the row names to be the TF names
rownames(regulon_genes) <- names(gene_list)
regulon_genes <- as.data.frame(regulon_genes)
regulon_genes$TF <- rownames(regulon_genes)
regulon_genes <- regulon_genes %>% dplyr::select(TF,everything())
regulon_genes <- as.data.frame(t(regulon_genes))
write.xlsx(regulon_genes,file=paste0("Cluster ",cluster,"RegulonswithGenes.xlsx"))



#####CALCULATE the number of genes for Each TF and Each Pathway
TFgenenumber <- data.frame('TF' = names(grouped), 'Regulongenes' = sapply(grouped, function(x) length(unique(x$gene))))
ipagenenumber <- data.frame('Pathway' = colnames(ipa), 'Pathwaygenes' = sapply(ipa, function(x) sum(!is.na(x))))
reactomegenenumber <- data.frame('Pathway' = colnames(reactome), 'Pathwaygenes' = sapply(reactome, function(x) sum(!is.na(x))))
kegggenenumber <- data.frame('Pathway' = colnames(kegg), 'Pathwaygenes' = sapply(kegg, function(x) sum(!is.na(x))))
hypertrophygenenumber <- data.frame('Pathway' = colnames(hypertrophy), 'Pathwaygenes' = sapply(hypertrophy, function(x) sum(!is.na(x))))
TGFbgenenumber <- data.frame('Pathway' = colnames(TGFb), 'Pathwaygenes' = sapply(TGFb, function(x) sum(!is.na(x))))
wikipathwaysgenenumber <- data.frame('Pathway' = colnames(wikipathways), 'Pathwaygenes' = sapply(wikipathways, function(x) sum(!is.na(x))))
PANTHERpathwaysgenenumber <- data.frame('Pathway' = colnames(PANTHERpathways), 'Pathwaygenes' = sapply(PANTHERpathways, function(x) sum(!is.na(x))))
panther2016genenumber <- data.frame('Pathway' = colnames(panther2016), 'Pathwaygenes' = sapply(panther2016, function(x) sum(!is.na(x))))

ECgenenumber <- data.frame('Pathway' = colnames(EC), 'Pathwaygenes' = sapply(EC, function(x) sum(!is.na(x))))

#phaggenenumber <- data.frame('Pathway' = colnames(phag), 'Pathwaygenes' = sapply(phag, function(x) sum(!is.na(x))))


###########CHOOSE GENESET TO USE
whatpathway <- "panther2016"
#Chose which Pathwaylist to USE
GenesforPathway <- get(whatpathway)#kegg
pathwaygenenumber <- get(paste0(whatpathway,"genenumber"))
###########CHOOSE GENESET TO USE

#Get#Intersection between Regulons and Genelist
# Initialize an empty list to store the results
# results <- list()
# # Iterate over each group in the grouped dataframe
# for (name in names(grouped)) {
#   group <- grouped[[name]]
# 
#   # Iterate over each column in the GenesofPathway dataframe
#   for (colname in colnames(GenesforPathway)) {
#     # Find the intersection of genes between the group and the GenesofPathway column
# 
#     group$gene <- tolower(group$gene)  # Convert 'name1' to upper case
#     GenesforPathway[[colname]] <- tolower(GenesforPathway[[colname]])  # Convert 'name2' to upper case
# 
# 
#     intersect_genes <- intersect(group$gene, GenesforPathway[[colname]])
#     #intersect_genes <- rbind(name,colname,intersect_genes)
# 
#     # Store the result in the list
#     results[[paste(name, colname, sep = "_")]] <- intersect_genes
#   }
# }




# Create a new list to store the results
results <- list()

#Retrive all Protein Coding Genes for Mus musculus for Fisher's Exact test
library(biomaRt)
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
#mm_genes <- getBM(attributes = "external_gene_name", mart = mart)
mm_genes <- getBM(attributes = c("external_gene_name"), filters = 'biotype', values = 'protein_coding', mart = mart)
mm_genes <- as.vector(mm_genes$external_gene_name)



# Iterate over each group
for (name in names(grouped)) {
  group <- grouped[[name]]
  
  # Iterate over each column in the GenesofPathway dataframe
  for (colname in colnames(GenesforPathway)) {
    # Find the intersection of genes between the group and the GenesofPathway column
    group$gene <- tolower(group$gene)
    GenesforPathway[[colname]] <- tolower(GenesforPathway[[colname]])
    intersect_genes <- intersect(group$gene, GenesforPathway[[colname]])
    # Create a contingency table
    # Number of genes in both the TF and the pathway
    both <- length(intersect_genes)
    # Number of genes in the TF but not in the pathway
    tf_only <- length(setdiff(group$gene, intersect_genes))
    # Number of genes in the pathway but not in the TF
    pathway_only <- length(setdiff(GenesforPathway[[colname]], intersect_genes))
    # Number of genes in neither the TF nor the pathway
    neither <- length(setdiff(mm_genes, union(group$gene, GenesforPathway[[colname]])))
    
    # Create the contingency table
    contingency_table <- matrix(c(both, tf_only, pathway_only, neither), nrow = 2, byrow = TRUE)
    dimnames(contingency_table) <- list(c("In TF", "Not in TF"), c("In Pathway", "Not in Pathway"))
    # Perform the Fisher's exact test
    fisher_test_results <- fisher.test(contingency_table)
    # Store the result in the list
    results[[paste(name, colname, sep = "_")]] <- list(Intersect_Genes = intersect_genes, Fisher_p_value = fisher_test_results$p.value)
  }
}


saveRDS(results, file=paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA",whatpathway," Cluster",cluster,"Resultsfiles.Rds"))




# Convert the results list to a data frame
results_df <- data.frame(TF_Column = names(results), 
                         Fisher_p_value = sapply(results, function(x) x$Fisher_p_value),
                         Intersect_Genes = I(lapply(results, function(x) x$Intersect_Genes)))

# Find the maximum number of Intersect_Genes
max_genes <- max(sapply(results_df$Intersect_Genes, length))

# Convert Intersect_Genes list into multiple columns
genes_df <- do.call(rbind, lapply(results_df$Intersect_Genes, function(x) {
  c(x, rep(NA, max_genes - length(x)))
}))

# Combine the dataframes
results_df <- cbind(results_df[, c("TF_Column", "Fisher_p_value")], genes_df)
#results_df <- data.frame(TF_Column = names(results), Intersect_Genes = I(results), Fisher_p_value = I(results[[Fisher_p_value]]))
#results_df$Intersectgenes <- results_df %>% unnest_wider(Intersect_Genes, names_sep = ",")
results_df$Geneintersect <- apply(results_df, 1, function(x) sum(!is.na(x))-1)
results_df <- results_df %>%filter(Geneintersect >2)
results_df <- results_df %>% separate(TF_Column,into=c("TF","Pathway"), sep="_")
results_df <- left_join(results_df,TFgenenumber, by="TF")
results_df <- left_join(results_df,pathwaygenenumber, by="Pathway")
results_df <-  results_df %>% mutate(Regulongeneratio=Geneintersect/Regulongenes)
results_df <-  results_df %>% mutate(Pathwaygeneratio=Geneintersect/Pathwaygenes)
results_df <-  results_df %>% mutate(Combinedratio=Regulongeneratio*Pathwaygeneratio) 
results_df <- results_df %>% mutate(CombinedMetric = Combinedratio * -log10(Fisher_p_value)) %>% arrange(desc(CombinedMetric))
results_df <-  results_df %>% dplyr::select(TF,Pathway,Geneintersect,Regulongenes,Regulongeneratio,Pathwaygenes,Pathwaygeneratio,Combinedratio,Fisher_p_value,CombinedMetric,everything()) %>% arrange( desc(CombinedMetric))
openxlsx::write.xlsx(results_df, file=paste0(whatpathway," - Cluster ",cluster,"_
                                   Allregulonsmateched.xlsx"))


######################SAVE respective results
# ipa_results_df <- results_df
# saveRDS(ipa_results_df, file=paste0("IPA_Results_df.Rds"))
ipa_results_df <- readRDS(file=paste0("../IPA_Results_df.Rds"))

# kegg_results_df <- results_df
# saveRDS(kegg_results_df, file=paste0("KEGG_Results_df.Rds"))
kegg_results_df <- readRDS(file=paste0("../KEGG_Results_df.Rds"))
#results_df <- kegg_results_df

# TGFb_results_df <- results_df
# saveRDS(TGFb_results_df, file=paste0("TGFb_Results_df.Rds"))
TGFb_results_df <- readRDS(file=paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/Fib/TGFb_Results_df.Rds"))
#results_df <- kegg_results_df

# wikipathways_results_df <- results_df
# saveRDS(wikipathways_results_df, file=paste0("wikipathways_Results_df.Rds"))
wikipathways_results_df <- readRDS(file=paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/MP/wikipathways_Results_df.Rds"))
#results_df <- kegg_results_df

# PANTHERpathways_results_df <- results_df
# saveRDS(PANTHERpathways_results_df, file=paste0("PANTHERpathways_Results_df.Rds"))
PANTHERpathways_results_df <- readRDS(file=paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/EC/PANTHERpathways_Results_df.Rds"))
#results_df <- kegg_results_df


# EC_results_df <- results_df
# saveRDS(EC_results_df, file=paste0("EC_Results_df.Rds"))
EC_results_df <- readRDS(file=paste0("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/EC/EC_Results_df.Rds"))
#results_df <- kegg_results_df



#TF for Fibroblast Cluster C22
cho <- results_df %>% dplyr::filter(TF %in% c("Foxp2","Srebf1","Atf5","Creb3l2","Egr2","Foxo1","Srebf2")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
write.xlsx(cho, file=paste0(whatpathway," - Cluster ",cluster,"_Chosenregulonsmateched.xlsx"))



#TF for CardiomyocyteCluster C9
cho <- results_df %>% dplyr::filter(TF %in% c("Gabpa","Atf6","Bach1","Smarca4","Jun","Fos","Smarcc2")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
write.xlsx(cho, file=paste0(whatpathway," - Cluster ",cluster,"_Chosenregulonsmateched.xlsx"))


#TF for CardiomyocyteCluster C2
cho <- results_df %>% dplyr::filter(TF %in% c("Irf5","Irf2","Irf4","Rel","Pbx3","Irf9","Stat2","Mitf","Nr3c1")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
write.xlsx(cho, file=paste0(whatpathway," - Cluster ",cluster,"_Chosenregulonsmateched.xlsx"))


cho_filtered <- PANTHERpathways_results_df %>% dplyr::filter(TF %in% c("Gabpa","Atf6","Bach1","Smarca4","Jun","Fos","Smarcc2")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)

#########################################################################################################################################
###GO-Chordplot Regulons per Pathways
######################################################################################################################################
#############CARDIOMYOCYTES######################
#SORT OUT ALL NON CARDIAC Terms in KEGG
whatpathway <- "kegg"
cho_filtered <- kegg_results_df %>% dplyr::filter(TF %in% c("Gabpa","Atf6","Bach1","Smarca4","Jun","Fos","Smarcc2")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique(cho$Pathway)#[c(-4,-6,-7,-10,-11,-13,-14,-15,-16,-20,-21,-23,-26,-27,-28,-30,-34,-35,-36,-37,-38,-40,-41,-42,-43,-44,-47,-48,-50,-56,-58,-59,-60,-61,-63,-64,-68,-72,-73,-78,-81,-86,-89,-94,-95,-97,-99,-100,-102,-104,-106,-110,-113,-116,-120,-121,-122,-125,-126,-130,-132,-134,-135,-138,-139,-141,-149,-152,-154,-158,-160,-161)]
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ] %>% filter(Fisher_p_value < 0.01)
#interpathways <- read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/interestingpathways.xlsx")
#cho_filtered <- cho_filtered %>% dplyr::filter(Pathway %in% interpathways[,1])
#write.xlsx(cho_filtered,file=paste0(whatpathway,"Chosen_Regulons_p0.01.xlsx"))
#SORT OUT ALL NON CARDIAC Terms in IPA

whatpathway <- "ipa"
cho <- ipa_results_df %>% dplyr::filter(TF %in% c("Gabpa","Atf6","Bach1","Smarca4","Jun","Fos","Smarcc2")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique(cho$Pathway)[c(-6,-7,-8,-9,-10,-12,-14,-31,-33,-40,-47,-51,-86,-96,-98,-107,-108,-110,-141,-180,-347,-348,-384)]
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
interpathways <- read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/ipa_interestingpathways.xlsx")
cho_filtered <- cho_filtered %>% filter(Pathway %in% interpathways[,1])
#write.xlsx(as.data.frame(unique(cho_filtered$Pathway)),file=paste0(whatpathway,"Chosen_Regulons_p0.01.xlsx"))


########################FIBROBLASTS
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. in Fibroblasts
whatpathway <- "TGFb"
cho <- TGFb_results_df %>% dplyr::filter(TF %in% c("Foxp2","Srebf1","Atf5","Creb3l2","Egr2","Foxo1","Srebf2")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique(cho$Pathway)#[c(-6,-7,-8,-9,-10,-12,-14,-31,-33,-40,-47,-51,-86,-96,-98,-107,-108,-110,-141,-180,-347,-348,-384)]
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]#%>% filter(Fisher_p_value < 0.01)
#interpathways <- read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/ipa_interestingpathways.xlsx")
cho_filtered <- cho_filtered %>% filter(Pathway %in% interpathways[,1])
#write.xlsx(as.data.frame(unique(cho_filtered$Pathway)),file=paste0(whatpathway,"Chosen_Regulons_p0.01.xlsx"))

########################FIBROBLASTS
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. in Fibroblasts
whatpathway <- "kegg"
cho <- kegg_results_df %>% dplyr::filter(TF %in% c("Foxp2","Srebf1","Atf5","Creb3l2","Egr2","Foxo1","Srebf2")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique(cho$Pathway)[c(-6,-7,-8,-9,-10,-12,-14,-31,-33,-40,-47,-51,-86,-96,-98,-107,-108,-110,-141,-180,-347,-348,-384)]
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho_filtered[!grepl('diabetic|cancer|Cancer|Prion|leukemia|carcinogen|Bacterial|Hepatitis|Longevity|Oocyte|Parathyroid|Amyotrophic|Amphetamine|Estrogen|carbon|Circadian|Thyroid', cho_filtered$Pathway), ]
#interpathways <- read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/ipa_interestingpathways.xlsx")
#cho_filtered <- cho_filtered %>% filter(Pathway %in% interpathways[,1])
#write.xlsx(as.data.frame(unique(cho_filtered$Pathway)),file=paste0(whatpathway,"Chosen_Regulons_p0.01.xlsx"))


#WNTIPA
whatpathway <- "ipa"
cho <- ipa_results_df %>% dplyr::filter(TF %in% c("Foxp2","Srebf1","Atf5","Creb3l2","Egr2","Foxo1","Srebf2","Tcf7l1","Nr3c1")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique(cho$Pathway)#[c(-6,-7,-8,-9,-10,-12,-14,-31,-33,-40,-47,-51,-86,-96,-98,-107,-108,-110,-141,-180,-347,-348,-384)]
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
#interpathways <- read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/ipa_interestingpathways.xlsx")
#cho_filtered <- cho_filtered %>% filter(Pathway %in% interpathways[,1])
cho_filtered <- cho_filtered %>% dplyr::filter(str_detect(Pathway, fixed("WNT")))





########################FIBROBLASTS
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. in Fibroblasts
whatpathway <- "ipa"
cho <- ipa_results_df %>% dplyr::filter(TF %in% c("Foxp2","Srebf1","Atf5","Creb3l2","Egr2","Foxo1","Srebf2")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique(cho$Pathway)#[c(-6,-7,-8,-9,-10,-12,-14,-31,-33,-40,-47,-51,-86,-96,-98,-107,-108,-110,-141,-180,-347,-348,-384)]
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho_filtered[!grepl('Melanoma|Glioblastoma|Leukemia|diabetic|cancer|Cancer|Prion|leukemia|carcinogen|Bacterial|Hepatitis|Longevity|Oocyte|Parathyroid|Amyotrophic|Amphetamine|Estrogen|carbon|Circadian|Thyroid', cho_filtered$Pathway), ]


########################Macrophages
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. in Macrophages
whatpathway <- "kegg"
cho <- kegg_results_df %>% dplyr::filter(TF %in% c("Irf5","Irf2","Irf4","Rel","Pbx3","Irf9","Stat2","Mitf","Nr3c1")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique((cho$Pathway)[c(-6,-7,-8,-9,-10,-12,-14,-31,-33,-40,-47,-51,-86,-96,-98,-107,-108,-110,-141,-180,-347,-348,-384)])
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho_filtered[!grepl('Asthma|Hematopoietic|Platelet|cancer|Hepatitis|Herpes|herpesvirus|papillomavirus|Leishmaniasis|Pertussis|Toxoplasmosis|Tuberculosis|carcinogenesis|Measeles|virus|Thyroid', cho_filtered$Pathway), ]
#cho_filtered <- cho_filtered[grepl('Phago', cho_filtered$Pathway), ]
#interpathways <- read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/ipa_interestingpathways.xlsx")
#cho_filtered <- cho_filtered %>% filter(Pathway %in% interpathways[,1])
#write.xlsx(as.data.frame(unique(cho_filtered$Pathway)),file=paste0(whatpathway,"Chosen_Regulons_p0.01.xlsx"))

########################Macrophages
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. in Macrophages
whatpathway <- "wikipathways"
cho <- wikipathways_results_df %>% dplyr::filter(TF %in% c("Irf5","Irf2","Irf4","Rel","Pbx3","Irf9","Stat2","Mitf","Nr3c1"))  %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique((cho$Pathway)[c(-1,-2,-3,-5,-6,-11,-13,-17,-20,-39,-50,-64,-68,-74,-75,-76,-78,-84,-99,-101,-103,-105,-109,-116,-117,-120,-129,-134, -137, -140, -142, -146, -150, -152, -154, -160, -164, -165, -167, -174, -184)])
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho_filtered[!grepl('Cancer', cho_filtered$Pathway), ]
cho_filtered <- cho_filtered[grepl('Phago', cho_filtered$Pathway), ]

########################Macrophages
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. in Macrophages
whatpathway <- "ipa"
cho <- ipa_results_df %>% dplyr::filter(TF %in% c("Irf5","Irf2","Irf4","Rel","Pbx3","Irf9","Stat2","Mitf","Nr3c1"))  %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique((cho$Pathway))#[c(-1,-2,-3,-5,-6,-11,-13,-17,-20,-39,-50,-64,-68,-74,-75,-76,-78,-84,-99,-101,-103,-105,-109,-116,-117,-120,-129,-134, -137, -140, -142, -146, -150, -152, -154, -160, -164, -165, -167, -174, -184)])
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho_filtered[!grepl('Influenza|diabetic|cancer|Cancer|Prion|leukemia|carcinogen|Bacterial|Hepatitis|Longevity|Oocyte|Parathyroid|Amyotrophic|Amphetamine|Estrogen|carbon|Circadian|Thyroid', cho_filtered$Pathway), ]


########################Macrophages
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. in Endothelials
whatpathway <- "ipa"
cho <- ipa_results_df %>% dplyr::filter(TF %in% c("Klf4","Bcl6b","Klf2","Irf1","Tcf4","Foxo1","Pbx3","Msx1"))  %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique((cho$Pathway))#[c(-1,-2,-3,-5,-6,-11,-13,-17,-20,-39,-50,-64,-68,-74,-75,-76,-78,-84,-99,-101,-103,-105,-109,-116,-117,-120,-129,-134, -137, -140, -142, -146, -150, -152, -154, -160, -164, -165, -167, -174, -184)])
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho_filtered[!grepl('Axonal|Senescence|Influenza|diabetic|cancer|Cancer|Prion|leukemia|carcinogen|Bacterial|Hepatitis|Longevity|Oocyte|Parathyroid|Amyotrophic|Amphetamine|Estrogen|carbon|Circadian|Thyroid', cho_filtered$Pathway), ]

########################EndothelialCells
#SORT OUT ALL NON CARDIAC Terms in TEndothelialCells
whatpathway <- "kegg"
cho <- kegg_results_df %>% dplyr::filter(TF %in% c("Klf4","Bcl6b","Klf2","Irf1","Tcf4","Foxo1","Pbx3","Msx1")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique((cho$Pathway)[c(-6,-7,-8,-9,-10,-12,-14,-31,-33,-40,-47,-51,-86,-96,-98,-107,-108,-110,-141,-180,-347,-348,-384)])
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho_filtered[!grepl('Inlfuenza|Asthma|Hematopoietic|Platelet|cancer|Hepatitis|Herpes|herpesvirus|papillomavirus|Leishmaniasis|Pertussis|Toxoplasmosis|Tuberculosis|carcinogenesis|Measeles|virus|Thyroid', cho_filtered$Pathway), ]
#cho_filtered <- cho_filtered[grepl('Phago', cho_filtered$Pathway), ]
#interpathways <- read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/ipa_interestingpathways.xlsx")
#cho_filtered <- cho_filtered %>% filter(Pathway %in% interpathways[,1])
#write.xlsx(as.data.frame(unique(cho_filtered$Pathway)),file=paste0(whatpathway,"Chosen_Regulons_p0.01.xlsx"))

########################EndothelialCells
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. in Macrophages
whatpathway <- "wikipathways"
cho <- wikipathways_results_df %>% dplyr::filter(TF %in% c("Klf4","Bcl6b","Klf2","Irf1","Tcf4","Foxo1","Pbx3","Msx1"))   %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique((cho$Pathway)[c(-1,-2,-3,-5,-6,-11,-13,-17,-20,-39,-50,-64,-68,-74,-75,-76,-78,-84,-99,-101,-103,-105,-109,-116,-117,-120,-129,-134, -137, -140, -142, -146, -150, -152, -154, -160, -164, -165, -167, -174, -184)])
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho_filtered[!grepl('Cancer', cho_filtered$Pathway), ]



########################EndothelialCells
#SORT OUT ALL NON CARDIAC Terms in TEndothelialCells
whatpathway <- "kegg"
cho <- kegg_results_df %>% dplyr::filter(TF %in% c("Klf4","Bcl6b","Klf2","Irf1","Tcf4","Foxo1","Pbx3","Msx1")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique((cho$Pathway)[c(-6,-7,-8,-9,-10,-12,-14,-31,-33,-40,-47,-51,-86,-96,-98,-107,-108,-110,-141,-180,-347,-348,-384)])
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho_filtered[!grepl('Inlfuenza|Asthma|Hematopoietic|Platelet|cancer|Hepatitis|Herpes|herpesvirus|papillomavirus|Leishmaniasis|Pertussis|Toxoplasmosis|Tuberculosis|carcinogenesis|Measeles|virus|Thyroid', cho_filtered$Pathway), ]
#cho_filtered <- cho_filtered[grepl('Phago', cho_filtered$Pathway), ]
#interpathways <- read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/ipa_interestingpathways.xlsx")
#cho_filtered <- cho_filtered %>% filter(Pathway %in% interpathways[,1])
write.xlsx(as.data.frame(unique(cho_filtered$Pathway)),file=paste0("C21-Fib ",whatpathway,"Chosen_Regulons_p0.01.xlsx"))

########################EndothelialCells
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. in Macrophages
whatpathway <- "wikipathways"
cho <- wikipathways_results_df %>% dplyr::filter(TF %in% c("Klf4","Bcl6b","Klf2","Irf1","Tcf4","Foxo1","Pbx3","Msx1"))   %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique((cho$Pathway)[c(-1,-2,-3,-5,-6,-11,-13,-17,-20,-39,-50,-64,-68,-74,-75,-76,-78,-84,-99,-101,-103,-105,-109,-116,-117,-120,-129,-134, -137, -140, -142, -146, -150, -152, -154, -160, -164, -165, -167, -174, -184)])
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho_filtered[!grepl('Cancer', cho_filtered$Pathway), ]




######################## REG FIBROBLASTS C21
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. in Fibroblasts
whatpathway <- "ipa"
cho <- ipa_results_df %>% dplyr::filter(TF %in% c("Foxp2","Gli3","Tef","Nr3c1","Pbx1")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique(cho$Pathway)#[c(-6,-7,-8,-9,-10,-12,-14,-31,-33,-40,-47,-51,-86,-96,-98,-107,-108,-110,-141,-180,-347,-348,-384)]
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho_filtered[!grepl('Melanoma|Glioblastoma|Leukemia|diabetic|cancer|Cancer|Prion|leukemia|carcinogen|Bacterial|Hepatitis|Longevity|Oocyte|Parathyroid|Amyotrophic|Amphetamine|Estrogen|carbon|Circadian|Thyroid', cho_filtered$Pathway), ]


######################## REG FIBROBLASTS C21
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. in Macrophages
whatpathway <- "kegg"
cho <- kegg_results_df %>% dplyr::filter(TF %in% c("Foxp2","Gli3","Tef","Nr3c1","Pbx1")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique(cho$Pathway)#[c(-6,-7,-8,-9,-10,-12,-14,-31,-33,-40,-47,-51,-86,-96,-98,-107,-108,-110,-141,-180,-347,-348,-384)])
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
#cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)
cho_filtered <- cho_filtered[!grepl('Asthma|Hematopoietic|Platelet|cancer|Hepatitis|Herpes|herpesvirus|papillomavirus|Leishmaniasis|Pertussis|Toxoplasmosis|Tuberculosis|carcinogenesis|Measeles|virus|Thyroid', cho_filtered$Pathway), ]
#cho_filtered <- cho_filtered[grepl('Phago', cho_filtered$Pathway), ]
#interpathways <- read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/ipa_interestingpathways.xlsx")
#cho_filtered <- cho_filtered %>% filter(Pathway %in% interpathways[,1])
#write.xlsx(as.data.frame(unique(cho_filtered$Pathway)),file=paste0(whatpathway,"Chosen_Regulons_p0.01.xlsx"))



########################Endothelial Cell GO LSIT
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. inEC
whatpathway <- "EC"
cho <- EC_results_df %>% dplyr::filter(TF %in% c("Klf4","Bcl6b","Klf2","Irf1","Tcf4","Foxo1","Pbx3","Msx1")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique(cho$Pathway)#[c(-6,-7,-8,-9,-10,-12,-14,-31,-33,-40,-47,-51,-86,-96,-98,-107,-108,-110,-141,-180,-347,-348,-384)]
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)



########################Endothelial Cell PANTHERpathways
#SORT OUT ALL NON CARDIAC Terms in TGFbeta. inEC
whatpathway <- "PANTHERpathways"
cho <- PANTHERpathways_results_df %>% dplyr::filter(TF %in% c("Klf4","Bcl6b","Klf2","Irf1","Tcf4","Foxo1","Pbx3","Msx1")) %>% arrange(desc(CombinedMetric)) %>% filter(Fisher_p_value < 0.05)
values_to_keep <- unique(cho$Pathway)#[c(-6,-7,-8,-9,-10,-12,-14,-31,-33,-40,-47,-51,-86,-96,-98,-107,-108,-110,-141,-180,-347,-348,-384)]
cho_filtered <- cho[cho$Pathway %in% values_to_keep, ]%>% filter(Fisher_p_value < 0.01)




############################
# Create a matrix for the links and make GOPlots
cho_filtered <- cho_filtered[!grepl('Eosinophils|Osteoarthritis|Gastrin|Neutrophils|Human|Neurons|Brain|Adipogenesis|Glioma|Arthritis|Glioblastoma|Melanoma|Leukemia|Opioid|Adenocarcinoma|Lupus|Axonal|Influenza|diabetic|cancer|Cancer|Prion|leukemia|carcinogen|Bacterial|Hepatitis|Longevity|Oocyte|Parathyroid|Amyotrophic|Amphetamine|Estrogen|carbon|Circadian|Thyroid', cho_filtered$Pathway), ]#cho_filtered <- cho_filtered[!grepl('Gastrin|Neutrophils|Neurons|Brain|Adipogenesis|Glioma|Arthritis|Glioblastoma|Melanoma|Leukemia|Opioid|Adenocarcinoma|Lupus|Axonal|Senescence|Influenza|diabetic|cancer|Cancer|Prion|leukemia|carcinogen|Bacterial|Hepatitis|Longevity|Oocyte|Parathyroid|Amyotrophic|Amphetamine|Estrogen|carbon|Circadian|Thyroid', cho_filtered$Pathway), ]

#first_three_values <- c("Foxp2","Srebf2","Creb3l2")
first_three_values <- unique(cho_filtered$TF)[1:6]
#first_three_values <- c("Klf2", "Klf4","Bcl6b","Irf1","Tcf4")
#first_three_values <- c("Klf2", "Klf4","Tcf4")
#first_three_values <- c("Mitf", "Irf2","Irf5","Irf9","Stat2","Rel")
#first_three_values <-c("Foxp2","Nr3c1","Creb3l2","Tcf7l1","Srebf2")


#first_three_values <- cho_filtered
cho_short <- cho_filtered %>% filter(TF %in% first_three_values) %>%group_by(TF) %>% slice_head(n = 10)
links <- cho_short %>% slice_max(n=40, order_by = CombinedMetric) %>% dplyr::select(TF,Pathway,Geneintersect)#,Combinedratio,Fisher_p_value) 
#links$Pathway <- wrap_text(links$Pathway, threshold = 10)
links$Pathway <-str_wrap(links$Pathway, width = 25, indent = 0, exdent = 5, whitespace_only = TRUE)
matrix <- as.data.frame(spread(links, key = TF, value = Geneintersect, fill = 0))
#matrix <- matrix[2:3,]
rownames(matrix) <- matrix[,1] 
matrix <- as.matrix(matrix %>% dplyr::select(-1))
matrix

#GOChord(matrix, space = 0.05 , gene.space = 0.3, gene.size = 3,ribbon.col = c("#E69F00","#009E73","#0072B2","grey"), border.size = 0.1)

# tiff(paste0("circos_",whatpathway," Cluster ",cluster,".tiff"), width = 300, height = 300, units = 'mm', res = 600, compression = "lzw")
# GOChord(matrix, space = 0.05 , gene.space = 0.3, gene.size = 5,ribbon.col = c("#E69F00","#009E73","#0072B2","grey"), border.size = 0.1)#,lfc.max=50)
# dev.off()

pdf(paste0("circos_",whatpathway," Cluster ",cluster,".pdf"), width = 12, height = 12)
GOChord(matrix, title=paste0(whatpathway," Cluster ",cluster), space = 0.05 , gene.space = 0.4, gene.size = 4,ribbon.col = c("#E69F00","#009E73","#0072B2","grey"), border.size = 0.05)#,lfc.max=50)
dev.off()

pdf(paste0("circos_",whatpathway,"Cluster ",cluster,".pdf"), width = 12, height = 12)
GOChord(matrix, title=paste0(whatpathway," Cluster ",cluster), space = 0.05 , gene.space = 0.4, gene.size = 4,ribbon.col = c("#E69F00","#009E73"), border.size = 0.05)#,lfc.max=50)
dev.off()

pdf(paste0("circos_",whatpathway,"WNT-Signaling Cluster ",cluster,".pdf"), width = 12, height = 12)
GOChord(matrix, title=paste0(whatpathway," Cluster ",cluster), space = 0.05 , gene.space = 0.4, gene.size = 4,ribbon.col = c("#E69F00","#009E73","#0072B2"), border.size = 0.05)#,lfc.max=50)
dev.off()


pdf(paste0("NEW_circos_",whatpathway,"Cluster ",cluster,".pdf"), width = 12, height = 12)
GOChord(matrix, title=paste0(whatpathway," Cluster ",cluster), space = 0.05 , gene.space = 0.4, gene.size = 4,ribbon.col = c("#440154FF","#31688EFF","#FDE725FF"), border.size = 0.05)#,lfc.max=50)
dev.off()


pdf(paste0("circos_",whatpathway,"Cluster ",cluster,".pdf"), width = 12, height = 12)
GOChord(matrix, title=paste0(whatpathway," Cluster ",cluster),space = 0.05 , gene.space = 0.3, gene.size = 4,ribbon.col = c("#440154FF","#31688EFF","#35B779FF","#FDE725FF"), border.size = 0.05)#,lfc.max=50)
dev.off()

pdf(paste0("circos_",whatpathway," Phagocytosiscentered Cluster ",cluster,".pdf"), width = 12, height = 12)
GOChord(matrix, title=paste0(whatpathway," Cluster ",cluster),space = 0.05 , gene.space = 0.3, gene.size = 5,ribbon.col = c("#440154FF","#31688EFF","#35B779FF","#FDE725FF"), border.size = 0.05)#,lfc.max=50)
dev.off()

# pdf(paste0("circos_",whatpathway,"ChosenTFs Cluster ",cluster,".pdf"), width = 12, height = 12)
# GOChord(matrix, title=paste0(whatpathway," Cluster ",cluster),space = 0.05 , gene.space = 0.3, gene.size = 4,ribbon.col = inferno(n = ncol(matrix)), border.size = 0.01)#,lfc.max=50)
# dev.off()

pdf(paste0("circos_",whatpathway,"ChosenTFs ",cluster,".pdf"), width = 12, height = 12)
GOChord(matrix, title=paste0(whatpathway," Cluster ",cluster),space = 0.05 , gene.space = 0.3, gene.size = 4,ribbon.col = plasma(n = ncol(matrix)), border.size = 0.01)#,lfc.max=50)
dev.off()

pdf(paste0("circos_",whatpathway,"ChosenTFs Cluster ",cluster,".pdf"), width = 12, height = 12)
GOChord(matrix, title=paste0(whatpathway," Cluster ",cluster),space = 0.05 , gene.space = 0.3, gene.size = 4,ribbon.col = magma(n = ncol(matrix)), border.size = 0.05)#,lfc.max=50)
dev.off()

pdf(paste0("circos_",whatpathway,"Cluster ",cluster,".pdf"), width = 12, height = 12)
GOChord(matrix, title=paste0(whatpathway," Cluster ",cluster),space = 0.05 , gene.space = 0.3, gene.size = 4,ribbon.col = viridis(n = ncol(matrix)), border.size = 0.05)#,lfc.max=50)
dev.off()

pdf(paste0("circos_",whatpathway," Cluster ",cluster,".pdf"), width = 12, height = 12)
GOChord(matrix, title=paste0(whatpathway," Cluster ",cluster),space = 0.05 , gene.space = 0.3, gene.size = 3,ribbon.col = c("#1B9E77","#7570B3","#D95F02","grey","#440154FF","#31688EFF"), border.size = 0.02)#,lfc.max=50)
dev.off()
##########################################################################################################################################################
##Produce Dataframe for DotPlot per Process#####
################################################################################################################################################


###########CHOOSE GENESET TO USE
whatpathway <- ""
#Chose which Pathwaylist to USE
GenesforPathway <- get(whatpathway)#kegg
pathwaygenenumber <- get(paste0(whatpathway,"genenumber"))
###########CHOOSE GENESET TO USE

#chosenpathway <-c("Gabpa","Atf6","Bach1","Smarca4","Jun","Fos","Smarcc2")
chosenpathway <-c("Irf5","Irf2","Rel","Pbx3","Irf9","Stat2","Mitf","Nr3c1")

chosenpathway <-c("Foxp2","Srebf1","Atf5","Creb3l2","Egr2","Foxo1","Srebf2","")
chosenpathway <- c("Pulmonary Fibrosis Idiopathic Signaling Pathway","Hepatic Fibrosis Signaling Pathway")
chosenpathway <- c("EGF receptor signaling pathway","Angiogenesis","PDGF signaling pathway","Integrin signalling pathway","VEGF signaling pathway","Wnt signaling pathway")
# 
# chosenpathway <- c("Cardiac Hypertrophy Signaling (Enhanced)","Autophagy","Protein Ubquitination Pathway","FGF Signaling")
# chosenpathway <- c("positive regulation of cardiac muscle hypertrophy","negative regulation of cardiac muscle hypetrophy")
# chosenpathway <- factor(c("Ubiquitin mediated proteolysis","Autophagy","Insulin signaling pathway","AMPK signaling pathway"), levels=c("Ubiquitin mediated proteolysis","Autophagy","Insulin signaling pathway","AMPK signaling pathway"))
# 
# chosenpathway <- c("nfat", "mef", "gata4", "srf", "foxo", "mitf", "YY1")

#Match PathwayGenes with Differential Geneexpression
# chosenpathway <- c("Cardiac Hypertrophy Signaling","Autophagy")
chosenpathway <- colnames(GenesforPathway)
#chosenpathway <- unique(c(TNTnames$`Ingenuity Canonical Pathways`,MHCnames$`Ingenuity Canonical Pathways`))

chosenpathway <- c("Foxp2","Srebf1","Atf5","Creb3l2","Egr2","Foxo1","Srebf2","Pbx1","Sox9","Nr3c1","Tcf7l1","Sp3","Cebpd")
chosenpathway <-c("Klf4","Bcl6b","Klf2","Irf1","Tcf4","Foxo1","Pbx3","Msx1")


chosenpathway <-c("Foxp2","Gli3","Tef","Nr3c1","Pbx1")
chosenpathway <- c("Pulmonary Fibrosis Idiopathic Signaling Pathway")


rm(datas)
allselected <- data.frame()
for (pathway in chosenpathway) {
  matchMHC <- dplyr::intersect(tolower(diffgeneMHC$names),tolower(GenesforPathway[[pathway]]))
  #matchMHC <- dplyr::intersect(tolower(diffgeneMHC$names),tolower([[pathway]]))
  selected_MHC <- diffgeneMHC[tolower(diffgeneMHC$names) %in% matchMHC, ]
  selected_MHC <- selected_MHC %>% dplyr::select(2,4,7) %>% mutate(`Group` = "MyHC")
  
  matchTNT <- dplyr::intersect(tolower(diffgeneTNT$names),tolower(GenesforPathway[[pathway]]))
  selected_TNT <- diffgeneTNT[tolower(diffgeneTNT$names) %in% matchTNT, ]
  selected_TNT <- selected_TNT %>% dplyr::select(2,4,7) %>% mutate(`Group` = "TnT")
  
  allselected <- rbind(allselected, rbind(selected_MHC,selected_TNT) %>% mutate(`Process` = pathway))
}

#allselected <- allselected %>% mutate(`Process` = "Hepatic & Pulmonary Fibrosis")

#Automatically Extract TopIPA terms from enrichment
###deleted 09/22/2023 Commit did
# names <- rbind(TNTnames,MHCnames)
# names <- names[!duplicated(names$Process), ]
# allselected_joined <-left_join(allselected,names, by="Process")
# datas <-allselected_joined

#allselected$Process <- factor(allselected$Process,levels=c("Ubiquitin mediated proteolysis","Autophagy","Insulin signaling pathway","AMPK signaling pathway"))
datas <- allselected

datas$Group <- factor(datas$Group,levels=c("TnT","MyHC"))

# Find the maximum absolute value of avg_log2FC
max_val <- max(abs(min(datas$avg_log2FC, na.rm = TRUE)), abs(max(datas$avg_log2FC, na.rm = TRUE)))
# Set the default theme with base_size = 14
theme_set(theme_gray(base_size = 18))
# Define a wrapping function


levels(wrap_text(datas$Process))
# Apply the wrapping function to your variable
#datas$namearrow_wrapped <- factor(wrap_text(datas$Process),levels=c("transforming growth\nfactor beta receptor\nsignaling pathway","positive regulation\nof transforming\ngrowth factor beta\nreceptor signaling\npathway","negative regulation\nof transforming\ngrowth factor beta\nreceptor signaling\npathway"))
datas$namearrow_wrapped <- factor(wrap_text(datas$Process))

datas$namearrow_wrapped <- gsub(" WP\\d+", "", datas$namearrow_wrapped)

#datas$namearrow_wrapped <-datas$Process
# Your ggplot code
p <- ggplot(datas, aes(x=avg_log2FC, y=rnorm(nrow(datas)))) +
  geom_point(aes(colour=p_val_adj < 0.01)) +
  scale_color_manual(values=c("black", "red")) +
  facet_grid(namearrow_wrapped ~ Group) +
  force_panelsizes(rows = 1, cols = unit(c(5, 5, 12), "cm"), TRUE) +
  #facet_grid2(namearrow ~ Group, scales = "free", independent = "all")+
  geom_vline(xintercept=0, color="black")+
  xlab("log2FC") +
  ylab("") + 
  theme(panel.background = element_rect(fill="white"),
        axis.title.y=element_blank(), # removes y-axis title
        axis.text.y=element_blank(), # removes y-axis text
        axis.ticks.y=element_blank(), # removes y-axis ticks
        legend.position="none", # removes legend
        text = element_text(face = "bold"), # makes all text bold
        strip.placement = "outside", # places facet labels outside the plot area
        strip.text.y = element_text(angle = 0),
        strip.clip = "off") + # rotates facet labels
  scale_x_continuous(limits = c(-max_val, max_val)) # sets x-axis limits


p <- p+ ggtitle(paste0("Differentially Regulated Genes \nby ",whatpathway," Process of Cluster ",cluster))
#ggsave(p, file=paste0("GenesofPathways ",whatpathway," Cluster ",cluster,"_",paste0(chosenpathway,collapse="-"),".pdf"), units="in", width=10, height=8, dpi=600,device="pdf")
p
ggsave(p, file=paste0("GenesofPathways ",whatpathway," Cluster ",cluster,"_","ScatterPlot",".pdf"), units="in", width=9, height=7, dpi=600,device="pdf")

}

paste(my_vector, collapse = " ")
##########################################################################################################################################################
#Heatmap for chosen Pathways
##########################################################################################################################################################


# save_pheatmap_pdf <- function(x, filename, width=3.5, height=12, margin = 0.2) {
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
# 
# 
# 
#   # Create a viewport with the desired margins.
#   viewport <- grid::grid.viewport(width = width, height = height, margin = margin)
# 
#   # Draw the heatmap within the viewport.
#   grid::grid.newpage()
#   grid::grid.draw(grid::grid.viewport.draw(viewport, x$gtable))
#   dev.off()
# }
###########CHOOSE GENESET TO USE
whatpathway <- "regulon_genes"
#Chose which Pathwaylist to USE
GenesforPathway <- get(whatpathway)#kegg
pathwaygenenumber <- get(paste0(whatpathway,"genenumber"))
###########CHOOSE GENESET TO USE


##Extract Gene-Groups like TFs to get gene names
# Create a vector of prefixes
prefixes <- c("Gata", "Nfat", "Jun", "Fos", "Smarc", "Bach", "MAF", "Gabpa","Atf")
prefixes <- c("Col", "Tgf", "Tcf", "Nfk", "Creb", "Fox", "Egr", "Sreb","Atf","Sp","Pbx1","Sox9")
prefixes <- c("Col", "Tgf")
prefixes <- "MMP|TIMP|LOXL|HSP"
chosenpathway <- "ColandTGFGenes"
chosenpathway <- "MMPs&TIMPs&LOX"

# Create a pattern that matches 'Phago' or 'phago' anywhere in the string
pattern <- paste0("(?i)(", paste(prefixes, collapse = "|"), ")")
GenesforPathway <- diffgeneTNT$names
indices <- grep(pattern, GenesforPathway, perl=TRUE)
chosengenes <- GenesforPathway[indices]
prefixes <- c("Creb", "Fox", "Egr", "Sreb","Atf","Pbx1","Sox9")
prefixes <- c("Tnf")
prefixes <- c("Klf2","Klf4")

chosenpathway <- "Collagen and TGF Genes"
# Create a pattern that matches any of the prefixes
prefixes <- c("Ar","Pge","Fiz","Spp","TGF","VEG,Ym","Nos")



# Your Term
prefixes <- "EGFR|Integrin-mediated|Wnt Signaling Pathway and Pluripotency WP723|MAPK signaling pathway WP493|Adhesion"
chosenpathway <- gsub("\\|", "&", prefixes)
#chosenpathway <- "MMPs&TIMPs&LOX"
# Create a pattern that matches 'Phago' or 'phago' anywhere in the string
pattern <- paste0("(?i)(", paste(prefixes, collapse = "|"), ")")
# Use grep() with perl=TRUE to enable case-insensitive matching
indices <- grep(pattern, colnames(GenesforPathway), perl=TRUE)
# Extract the matching names
chosenpathway <- colnames(GenesforPathway)[indices]
chosenpathway 
# Order the chosenpathway vector
chosenpathway <- chosenpathway[chosenpathway %>% naturalsort::naturalorder()]


chosengenes <- as.vector(read.xlsx("/Users/eliascrapa/Library/CloudStorage/OneDrive-UCSF/DATA/SingleCellSeq/Analysis/All/CrossATACIPA/DB/TFGenelist.xlsx")[,1])
chosenpathway <- c("negative regulation of cardiac muscle hypetrophy","negative regulation of cardiac muscle hypertrophy in response to stress","positive regulation of cardiac muscle hypertrophy","positive regulation of cardiac muscle hypertrophy in response to stress")
chosenpathway <- c("positive regulation of cardiac muscle hypertrophy","positive regulation of cardiac muscle hypertrophy in response to stress")
chosenpathway <- c("negative regulation of cardiac muscle hypetrophy","negative regulation of cardiac muscle hypertrophy in response to stress")
chosenpathway <- c("Cardiac Hypertrophy Signaling (Enhanced)")
chosenpathway <- c("Pulmonary Fibrosis Idiopathic Signaling Pathway")
chosenpathway <- c("Bach1")
chosenpathway <- c("Bach2")
chosenpathway <- c("Atf6")
chosenpathway <- c("Jun")
chosenpathway <- c("Gabpa")
chosenpathway <- c("Smarca4")
chosenpathway <- c("Smarcc2")
chosenpathway <-"negative regulation of transforming growth factor beta receptor signaling pathway"
chosenpathway <- "positive regulation of transforming growth factor beta receptor signaling pathway"
chosenpathway <- "transforming growth factor beta receptor signaling pathway"  

chosenpathway <- c("EGFR1 Signaling Pathway WP572","Integrin-mediated Cell Adhesion WP6","Focal Adhesion WP85","Focal Adhesion-PI3K-Akt-mTOR-signaling pathway WP2841")


chosenpathway <- colnames(TGFb)
whatgenes <- "GOTERM"
whatgenes <- chosenpathway

chosenpathway <- "TFGenes Fib Regulons_2"
chosenpathway <- "TFGenes EC Regulons"
chosenpathway <- "TFGenes EC Regulons"
chosenpathway <- "ColandTGFGenes"
chosenpathway <- "Hepatic & Pulmonary Fibrosis Genes"
chosenpathway <- "RoselleChosenGenes"
chosenpathway <- "TNFGenes"
chosenpathway <- "MerTK"
chosenpathway <- "GO Phagosome Maturation"


chosenpathway <- c("Foxp2","Srebf1","Creb3l2","Egr2","Foxo1","Srebf2","Pbx1","Sox9","Nr3c1","Tcf7l1","Sp3","Crebbp","Atf5")
chosenpathway <-c("Irf5","Irf2","Irf4","Rel","Pbx3","Irf9","Stat2","Mitf","Nr3c1")
#chosengenes <-c("Irf5","Irf2","Pbx3","Irf9","Stat2","Mitf","Nr3c1","Rel")
chosengenes <-ILTNF$Approved.symbol
chosenpathway <- "Stiffness related DEGs in HUVEC"
chosengenes <- intersect(tolower(c("ADAMTSL1","LAMB3","SMOC1","STC1","ZNF703","ZSCAN31","ID1","KLF10","CXCL12","GADD45B","KRT7","INSR-A","EPHA5","TGFb2")),tolower(diffgeneTNT$names))
chosenpathway <- c("Cebpd","Cebpb","Cebpa")
chosenpathway <-c("Klf4","Bcl6b","Klf2","Irf1","Tcf4","Foxo1","Msx1","Pbx3")
chosenpathway <-c("Klf4","Klf2")
chosenpathway <- colnames(GenesforPathway)

WNTsig <- vector()
for (pathway in chosenpathway) {
  WNTsig <- c(WNTsig,GenesforPathway[[pathway]])
}
chosengenes <- WNTsig
chosenpathway <- ("Endocytosis")

#Make Gene-List for Fibroblalsstx
chosenpathway <- as.data.frame(c(tolower(GenesforPathway[["Pulmonary Fibrosis Idiopathic Signaling Pathway"]]),tolower(GenesforPathway[["Hepatic Fibrosis Signaling Pathway"]])))
colnames(chosenpathway) <- "Hepatic & Pulmonary Fibrosis Genes"
chosenpathway$`TGF-β Signaling` <- GenesforPathway$`TGF-β Signaling`
chosenpathway$`WNT Signaling`<- WNTsig

GenesforPathway <- chosenpathway 
chosenpathway<- colnames(GenesforPathway)

chosenpathway<-c("Pbx1","Atf6","Bach1","Smarca4","Jun","Smarcc2","Foxp1","Nr3c1")

chosenpathway <- "Myofibroblastmarker"
chosengenes <- c("Postn","Col1a1","Col3a1","ACTA2","PDGFRA","PDGFRB","PCOLCE","INHBA","FN1","Timp1", "Thbs1", "Thbs4", "Ccn2","Meox1","Col8a1")


chosenpathway <- "EC14genesR1"
chosengenes <- c("Agtr1a", "AT2R", "Nos", "NOX")

chosenpathway <- c("Clathrin-mediated Endocytosis Signaling","Caveolar-mediated Endocytosis Signaling")

PDGFGenes <- read.xlsx("PDGF-Genes.xlsx")[,3]
chosengenes <- PDGFGenes
chosenpathway <- ("PDGF Signaling Genes")

chosenpathway <- colnames(PDGF)[2]
chosengenes <- PDGF[2]
chosenpathway <- colnames(PDGF)[3]
chosengenes <- PDGF[3]

PDGFGenes <- read.xlsx("PDGF-Genes.xlsx")[,1]
chosengenes <- PDGFGenes
chosenpathway <- ("PDGF Receptor Activity")

prefixes <- "Agt|Nos|Nox"
prefixes <- paste0(c("ERK","MEK","Raf","Ras","Rap1","Sos","Grb2","Shc","Fyn","Cav-1","C3g","Crk","CAS","PAX","FAK","c-Src","Integrin"),collapse="|")


chosenpathway <- gsub("\\|", "&", prefixes)


#chosenpathway <- "MMPs&TIMPs&LOX"
# Create a pattern that matches 'Phago' or 'phago' anywhere in the string
#pattern <- paste0( paste(prefixes, collapse = "|"), ")","(?i)(")
pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")","(?i)")
GenesforPathway <- diffgeneTNT$names
indices <- grep(pattern, GenesforPathway, perl=TRUE)
chosengenes <- GenesforPathway[indices]


chosengenes <- c(EC$`cellular response to fluid shear stress`, "Nos3","Nostrin")
chosenpathway <- "Sheer Stress induces genes"


# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

#Extract mousdgenes from Input list
chosenpathway <- "Sheer Stress induces genes_graphlist"
chosengenes <- paste0(c("ERK","MEK","Raf","Ras","Rap1","Sos","Grb2","Shc","Fyn","Cav-1","C3g","Crk","CAS","PAX","FAK","c-Src","Integrin","Nesprin2", "AP-2", "GEM", "TFID", "FOS", "LDLR", "LaminA", "Stat1","Stat","Stat3","Stat5","Stat6","BCL2L1","IRF1","IFNG","IL4","CCND2","E2F1","KLF4","KLF5","Sp1","Emerin","Nfe2","Sirt1","HDAC5","Klf2","Erk5","Mef2","AMPK","MEF2A","MEF2C","BMPR2","IL1R1","NOS3","CAV1","PDGFB","PTK2","CDH5","AKT3","KDR","PECAM1","MAP2K5","MAP3K5"))

# Use the ensembl dataset
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart=ensembl)
# Get the mouse genes
mouse_genes <- getBM(attributes=c('mgi_symbol'), filters = 'mgi_symbol', values = chosengenes, mart = ensembl)
mouse_genes

chosenpathway <- ("Klf2/Klf4")

for (pathway in chosenpathway) {
  ####CHOOSE if Pathwaylist or Genelist
  #filteringgenes <- tolower(GenesforPathway[[pathway]])
  filteringgenes <- tolower(chosengenes)
  
  matchMHC <- dplyr::intersect(tolower(diffgeneMHC$names),filteringgenes)
  selected_MHC_nopcutoff <- diffgeneMHC[tolower(diffgeneMHC$names) %in% matchMHC, ] %>% mutate(Group = "MyHC")
  selected_MHC <- selected_MHC_nopcutoff %>% dplyr::select(2,4,7) %>% filter(p_val_adj < 0.01) %>% dplyr::select(-"p_val_adj")
  
  matchTNT <- dplyr::intersect(tolower(diffgeneTNT$names),filteringgenes)
  selected_TNT_nopcutoff <- diffgeneTNT[tolower(diffgeneTNT$names) %in% matchTNT, ] %>% mutate(Group = "TnT")
  selected_TNT <- selected_TNT_nopcutoff %>% dplyr::select(2,4,7)  %>% filter(p_val_adj < 0.01) %>% dplyr::select(-"p_val_adj")
 
  
  merged_df_nopcutoff <- rbind(selected_MHC_nopcutoff, selected_TNT_nopcutoff) %>% arrange(desc(names)) %>% dplyr::select(names,everything())
  openxlsx::write.xlsx(merged_df_nopcutoff,file=paste0(str_replace_all(pathway, "/", "_"),"_","_Cluster ",cluster,"_","HeatmapGenexpressionTable_Log2Fold.xlsx"))
  
  
  
  # Assuming selected_MHC and selected_TNT are your dataframes
  merged_df <- try(merge(selected_TNT,selected_MHC, by = "names", all = TRUE) %>%
    setNames(c("Gene","TnT","MyHC")) %>%
    replace_na(list(MyHC = 0, TnT = 0)) %>%
    arrange(desc(MyHC)) %>%
    dplyr::slice(1:30) %>%
    tibble::column_to_rownames(var = "Gene"))
  breaksList = seq(0, max(merged_df$MyHC,merged_df$TnT)+0.1, by = 0.05)
  merged_df <-as.matrix(merged_df)
 
  heatmap <- pheatmap(merged_df, main = paste0(wrap_text(pathway)), display_numbers = F, cluster_rows = FALSE, cluster_cols = FALSE,scale = "none", cellwidth = 18, cellheight = 13, fontsize_row = 16, breaks = unique(breaksList), color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)))
  
  
  if (nrow(merged_df) >0 ){
    save_pheatmap_pdf(heatmap, paste0(str_replace_all(pathway, "/", "_"),"_MyHC_Heatmedr
                                      ap_Cluster ",cluster,".pdf"))
  } else {"Dataframe without any value"}
   
  
   merged_df <- try(merge(selected_TNT,selected_MHC, by = "names", all = TRUE) %>%
     setNames(c("Gene","TnT","MyHC")) %>%
     replace_na(list(MyHC = 0, TnT = 0)) %>%
     arrange(desc(TnT)) %>%
     dplyr::slice(1:30) %>%
     tibble::column_to_rownames(var = "Gene"))
   breaksList = seq(0, max(merged_df$MyHC,merged_df$TnT)+0.1, by = 0.05)
   merged_df <-as.matrix(merged_df)
   
   heatmap <- pheatmap(merged_df, main = paste0(wrap_text(pathway)), display_numbers = F, cluster_rows = FALSE, cluster_cols = FALSE,scale = "none", cellwidth = 18, cellheight = 13, fontsize_row = 16, breaks = unique(breaksList), color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)))
   

   if (nrow(merged_df) >0 ){
     save_pheatmap_pdf(heatmap, paste0(str_replace_all(pathway, "/", "_"),"_TnT_Heatmap_Cluster ",cluster,".pdf"))
   } else {"Dataframe without any value"}
  
  
   
   
   
   
}


library(grid)

# Your existing code to generate the heatmap
heatmap <- pheatmap(merged_df, main = paste0(wrap_text(pathway)), display_numbers = F, cluster_rows = FALSE, cluster_cols = FALSE,scale = "none", cellwidth = 18, cellheight = 13, fontsize_row = 16, breaks = unique(breaksList), color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)))

# Add vertical text
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1, 3), "null"))))
print(heatmap$gtable, vp = viewport(layout.pos.col = 2))

# Increase font size and move text closer to the plot
grid.text(paste0(wrap_text(pathway)), x = 0.8, rot = 90, gp = gpar(fontsize = 20), vp = viewport(layout.pos.col = 1))

###################################################################################################################
#JOIN all GO-Terms to have an nice List of Cardiac Hypertrophy Genes
# Check each vector to make sure it is valid
hypertrophy_vectors <- list(
  unique(c(hypertrophy$`cardiac muscle hypertrophy`, hypertrophy$`cardiac muscle hypertrophy in response to stress`, hypertrophy$`regulation of cardiac muscle hypertrophy`, hypertrophy$`regulation of cardiac muscle hypertrophy in response to stress`)),
  unique(c(hypertrophy$`negative regulation of cardiac muscle hypertrophy in response to stress`, hypertrophy$`negative regulation of cardiac muscle hypetrophy`)),
  unique(c(hypertrophy$`positive regulation of cardiac muscle hypertrophy in response to stress`, hypertrophy$`positive regulation of cardiac muscle hypertrophy`))
)

# Get the maximum length of the elements in the list
max_length <- max(lengths(hypertrophy_vectors))

# Fill each element with NA up to the maximum length
hypertrophy_vectors_padded <- lapply(hypertrophy_vectors, function(x) {
  # If the element is shorter than the maximum length, add NA
  if (length(x) < max_length) {
    x <- c(x, rep(NA, max_length - length(x)))
  }
  
  # If the element is longer than the maximum length, remove the excess elements
  if (length(x) > max_length) {
    x <- x[1:max_length]
  }
  
  return(x)
})

# Add the names to the list elements


# Print the padded list
hypertrophy_joined <- as.data.frame(hypertrophy_vectors_padded)

colnames(hypertrophy_joined) <- c("Cardiac Muscle Hypertrophy",
                                       "negative regulation of cardiac muscle hypertrophy",
                                       "positive regulation of cardiac muscle hypertrophy")

GenesforPathway <- hypertrophy_joined

summary(tolower(hypertrophy_joined$`Cardiac Muscle Hypertrophy`) %in% tolower(ipa$`Cardiac Hypertrophy Signaling (Enhanced)`))
