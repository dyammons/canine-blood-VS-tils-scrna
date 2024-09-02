#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")

### Analysis note: 

############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin B cell preprocessing   ######## <<<<<<<<<<<<<<
############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set output params
outName <- "bcell"

#subset on the B cells & complete subset analysis
table(seu.obj$conSense)[grepl("B",names(table(seu.obj$conSense)))]
Idents(seu.obj) <- "conSense"
seu.obj.sub <- subset(seu.obj, 
                          ident = c("Activated B cell", "B cell", "Class switched B cell", "Immature B cell", 
                                    "Naive B cell","Plasma cell")
                         )

table(seu.obj.sub$conSense)
min(table(seu.obj.sub$name)) > 100
seu.sub.list <- SplitObject(seu.obj.sub, split.by = "name")
seu.sub.list <- seu.sub.list[-c(3,5,12,14,16)]#too few cells - innacurate clsutering

seu.obj <- indReClus(seu.obj = NULL, outDir = "./output/s2/", subName = "20230507_bcell_bloodANDtils", preSub = T, seu.list = seu.sub.list,
                      vars.to.regress = "percent.mt",nfeatures = 2000
                       )

clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "20230507_bcell_bloodANDtils", test_dims = c(40,35,30), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "20230507_bcell_bloodANDtils", final.dims = 30, final.res = 0.4, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.5, n.neighbors = 60, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                      )

########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin B cell analysis   ######## <<<<<<<<<<<<<<
########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in processed data
seu.obj <- readRDS("../output/s3/20230507_bcell_bloodANDtils_res0.4_dims30_dist0.5_neigh60_S3.rds")
seu.obj$cellSource <- ifelse(grepl("tils",seu.obj$name),"TILs","Blood")
outName <- "bcell"

#check consenseous
seu.obj$ct <- ifelse(seu.obj$cellSource == "TILs", paste0("TIL;", seu.obj$celltype.l3), paste0("Blood;",seu.obj$celltype.l3_pbmc))
table(seu.obj$ct, seu.obj$clusterID_sub) %>% melt() %>% separate(Var.1, sep = ";", c("source", "ct")) %>% filter(source == "Blood") %>% group_by(Var.2) %>% mutate(pct = value/sum(value)) %>% filter(value == max(value))
table(seu.obj$ct, seu.obj$clusterID_sub) %>% melt() %>% separate(Var.1, sep = ";", c("source", "ct")) %>% filter(source == "TIL") %>% group_by(Var.2) %>% mutate(pct = value/sum(value)) %>% filter(value == max(value))

Idents(seu.obj) <- "cellSource"
seu.obj.ds <- subset(seu.obj, downsample = min(table(seu.obj$cellSource)))
table(seu.obj.ds$ct, seu.obj.ds$clusterID_sub) %>% sweep(., 2, colSums(.), `/`) %>% melt() %>% separate(Var.1, sep = ";", c("source", "ct")) %>% group_by(source, Var.2) %>% filter(value == max(value))

#set metadata
colz.df <- read.csv("./metaData/majorGroups.csv")
colz.df <- colz.df[colz.df$majorID2 == "b" | colz.df$majorID2 == "plasma", ]

table(seu.obj$clusterID_sub,seu.obj$conSense)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "#F8766D", "1" = "#A3A500", 
                                   "2" = "#00BF7D", "3" = "#00B0F6",
                                   "4" = "#E76BF3")
                       )
seu.obj$sub_colz <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Naive (c0)", "1" = "Class_switched (c1)", 
                                   "2" = "Plasma_cell (c2)", "3" = "Immature (c3)",
                                   "4" = "Plasma_cell (c2)")
                       )
seu.obj$majorID_sub <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "0", "1" = "1", 
                                   "2" = "2", "3" = "3",
                                   "4" = "2")
                       )
seu.obj$clusID_new <- Idents(seu.obj)


### Supp data - Generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, outName = outName, returnViln = F,
                     outDir = paste0("../output/viln/", outName, "/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
                    )

### Fig 4a - UMAP by clusterID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusID_new",
              cols = colz.df$colour[c(1,4,5,2)],
              pt.size = 0.5,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, 
                alpha = 0.8, labCol = colz.df$labCol[c(1,4,5,2)]) + NoLegend() + theme(axis.title = element_blank(),
                                                                                                               panel.border = element_blank())
ggsave(paste("../output/", outName, "/", "rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig 4b - skew plot for abundance analysis
p <- skewPlot(seu.obj, groupBy = "majorID_sub", yAxisLabel = "Percent of B cells",
              dout = paste0("../output/", outName), outName = outName)
ggsave(paste0("../output/", outName, "/", "skewPlot.png"), width = 6, height = 4)


### Fig 4c - DGE analysis

#exclude plasma cells from DGE analysis
seu.obj <- subset(seu.obj, invert = T,
                  subset = majorID_sub == "Plasma_cell (c2)"
                  )

#load in the tumor and pal signatures to exlude from DE analysis
pal_feats = c('TIMP1', 'NAA10', 'ENSCAFG00000037735', 'GP6', 'SEC11C', 'FTL', 'NRGN', 'ACOT7', 'VCL', 'RSU1', 'ITGB1', 'H3-3A', 'RABGAP1L', 'SELP', 'SH3GLB1', 'ACTB', 'ENSCAFG00000008221', 'TLN1', 'GSN', 'AMD1', 'TREM2', 'SH3BGRL2', 'MYH9', 'PLEK', 'ENSCAFG00000042554', 'RAP1B', 'ENSCAFG00000004260', 'NAP1L1', 'PPBP', 'RASA3', 'ITGA2B', 'EIF1', 'ACTG1', 'C9H17orf64', 'JMJD6', 'CCL14', 'GNG11', 'IGF2BP3', 'TBXAS1', 'VDAC3', 'MARCHF2', 'TPM4', 'TKT', 'FTH1.1', 'FERMT3', 'RTN3', 'PRKAR2B', 'SVIP', 'ENSCAFG00000030286', 'ADA', 'MYL9', 'TUBB1', 'TUBA1B', 'METTL7A', 'THBS1', 'SERF2', 'PIF1', 'B2M', 'GAS2L1', 'YWHAH', 'HPSE', 'ATG3', 'ENSCAFG00000015217', 'ITGA6','RGS18', 'SUB1', 'LGALS1', 'CFL1', 'BIN2', 'CAT', 'RGS10', 'MGST3', 'TMBIM6', 'PFN1', 'CD63', 'RALBP1', 'GNAS', 'SEPTIN7', 'TPT1', 'UBB', 'ATF4', 'BBLN', 'MTDH', 'ENSCAFG00000017655','FYB1', 'ENO1', 'GABARAP', 'SSR4', 'MSN', 'ENSCAFG00000011134', 'ENSCAFG00000046637', 'COX8A', 'DLA-64', 'CD47', 'VASP', 'DYNLRB1', 'DLA88', 'SMDT1', 'ATP5PF','ELOB', 'ENSCAFG00000029155', 'ARPC3', 'VPS28', 'LRRFIP1', 'SRP14', 'ABRACL', 'ENSCAFG00000043577', 'ENSCAFG00000042598')
tumor.sig <- read.csv("./metaData/tumorSig.csv", header = T)$x

createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name",
         outDir = paste0("../output/", outName, "/pseudoBulk/"), 
         grepTerm = "tils", grepLabel = c("TILs", "Blood"), featsTOexclude = c(pal_feats,tumor.sig), lowFilter = T, dwnSam = F)

p_volc <- pseudoDEG(inDir = paste0("../output/", outName, "/pseudoBulk/"), metaPWD = paste0("../output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"),
                    outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = outName,
                    padj_cutoff = 0.01, lfcCut = 0.58,  idents.1_NAME = "TILs", idents.2_NAME = "Blood", title = "TILS vs Blood", 
                    fromFile = T, returnVolc = T, filterTerm = "^ENSCAF", mkDir = T)

p <- prettyVolc(plot = p_volc[[1]], rightLab = NULL, leftLab = NULL, rightCol = "red", leftCol = "blue", arrowz = F
                    ) + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
                              axis.line = element_blank())
ggsave(paste0("../output/", outName, "/", "volcPlot.png"), width = 7, height = 7)


### Fig 2d - GO GSEA of DEGs
p <- plotGSEA(pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv"),
              geneList = NULL, category = "C5", species = "dog", termsTOplot = 10, upOnly = T, trunkTerm = T,
              pvalueCutoff = 0.05, subcategory = NULL, 
              saveRes = paste0("../output/", outName, "/c5_", outName, "_res.csv")) + theme(axis.title=element_text(size = 16))
p <- p + scale_x_continuous(limits = c(-80,ceiling(max(p$data$x_axis)*1.05)), 
                            breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
                            name = "log10(p.adj)") + ggtitle("Gene ontology") + theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste0("../output/", outName, "/", "gseaPlot_1.png"), width = 7, height = 7)


### Fig 2e - Reactome GSEA of DEGs
p <- plotGSEA(pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv"),
         geneList = NULL, category = "C2", species = "dog", termsTOplot = 10, upOnly = T, trunkTerm = T,
                     pvalueCutoff = 0.05, subcategory = "CP:REACTOME", saveRes = paste0("../output/", outName, "/c2_", outName, "_res.csv")
                    ) + theme(axis.title=element_text(size = 16))
p <- p + scale_x_continuous(limits = c(-30,ceiling(max(p$data$x_axis)*1.05)), breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),name = "log10(p.adj)") + ggtitle("Reactome") + theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste0("../output/", outName, "/", "gseaPlot_2.png"), width = 7, height = 7)


### Fig 2f - Split UMAP of selected DEGs
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$cellSource)))

features <- c("LYZ", "BTLA", "ADD3", "IGHM", "FOS", "FOSB", "FAS", "IGFLR1", "EGFL7")

p <- FeaturePlot(seu.obj.sub,features = features, pt.size = 0.1, split.by = "cellSource", order = T, cols = c("lightgrey","darkblue"), by.col = F) + labs(x = "UMAP1", y = "UMAP2") & theme(axis.text = element_blank(),
                                                                                                                                                axis.title.y.right = element_text(size = 11),
                                                                                                                                            
                                                                                                                           axis.ticks = element_blank(),
                                                                                                                           axis.title = element_blank(),
                                                                                                                           axis.line = element_blank(),
                                                                                                                           plot.title = element_text(size=11),
                                                                                                                           title = element_blank(),
                                                                                                                           plot.margin = unit(c(0, 0, 0, 0), "cm")
                                                                                                                          ) 
ggsave(paste("../output/", outName, "/", "splitFeats.png", sep = ""), width = 12, height = 4)


res.df <- read.csv(paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv"))
res.df <- res.df[!grepl("^ENS", res.df$gene), ]
geneList_UP <- res.df %>% filter(padj < 0.1) %>% filter(log2FoldChange > 1) %>% pull(gene)
geneList_DWN <- res.df %>% filter(padj < 0.1) %>% filter(log2FoldChange < -1) %>% pull(gene)

seu.obj$cellSource <- as.factor(seu.obj$cellSource)
p <- splitDot(
    seu.obj = seu.obj, groupBy = "majorID_sub", splitBy = "cellSource", buffer = 150,
    namedColz = setNames(c("#F8766D", "#00BFC4"),  c("Blood", "TILs")), 
    geneList_UP = geneList_UP[1:20], geneList_DWN = geneList_DWN[1:20], geneColz = c("red", "blue")
)
ggsave(plot = p, paste0("../output/", outName, "/", outName, "_splitDot.png"), width = 11, height = 7)


## too few to do within cluster de
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end B cell analysis   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


