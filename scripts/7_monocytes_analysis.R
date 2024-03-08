#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")

### Analysis note: 

################################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin monocyte preprocessing   ######## <<<<<<<<<<<<<<
################################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set output params
outName <- "mono"

#subset on the monocytes & complete subset analysis
seu.obj.sub <- subset(seu.obj, 
                          idents = c("IFN-TAM","M-MDSC" ,"Monocyte, CD4+","Monocyte, CD4-","Monocyte, IFN signature",
                                     "TIM")
                         )
table(seu.obj.sub$conSense)
table(seu.obj.sub$orig.ident)

seu.sub.list <- SplitObject(seu.obj.sub, split.by = "orig.ident")

seu.obj <- indReClus(seu.obj = NULL, outDir = "./output/s2/", subName = "20230507_mono_bloodANDtils", preSub = T, seu.list = seu.sub.list,
                      vars.to.regress = "percent.mt"
                       )

clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "20230507_mono_bloodANDtils", test_dims = c(40,35,30), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "20230507_mono_bloodANDtils", final.dims = 40, final.res = 0.5, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.6, n.neighbors = 75, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                      )

############################################ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin monocyte analysis   ######## <<<<<<<<<<<<<<
############################################ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in processed data
seu.obj <- readRDS("../output/s3/20230507_mono_bloodANDtils_res0.5_dims40_dist0.6_neigh75_S3.rds")
seu.obj$cellSource <- ifelse(grepl("tils",seu.obj$name),"TILs","Blood")
outName <- "mono"

#check consenseous
seu.obj$ct <- ifelse(seu.obj$cellSource == "TILs", paste0("TIL;", seu.obj$celltype.l3), paste0("Blood;",seu.obj$celltype.l3_pbmc))
table(seu.obj$ct, seu.obj$clusterID_sub) %>% melt() %>% separate(Var.1, sep = ";", c("source", "ct")) %>% filter(source == "Blood") %>% group_by(Var.2) %>% mutate(pct = value/sum(value)) %>% filter(value == max(value))
table(seu.obj$ct, seu.obj$clusterID_sub) %>% melt() %>% separate(Var.1, sep = ";", c("source", "ct")) %>% filter(source == "TIL") %>% group_by(Var.2) %>% mutate(pct = value/sum(value)) %>% filter(value == max(value))

Idents(seu.obj) <- "cellSource"
seu.obj.ds <- subset(seu.obj, downsample = min(table(seu.obj$cellSource)))
table(seu.obj.ds$ct, seu.obj.ds$clusterID_sub) %>% sweep(., 2, colSums(.), `/`) %>% melt() %>% separate(Var.1, sep = ";", c("source", "ct")) %>% group_by(source, Var.2) %>% filter(value == max(value))

#set metadata
colz.df <- read.csv("./metaData/majorGroups.csv")
colz.df <- colz.df[colz.df$majorID2 == "mono", ]

table(seu.obj$clusterID_sub,seu.obj$conSense)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "CD4neg_monocyte (c0)", "1" = "CD4Pos_monocyte (c1)", 
                                   "2" = "CD4neg_monocyte (c0)", "3" = "CD4neg_monocyte (c0)",
                                   "4" = "IFN_monocyte (c2)", "5" = "M-MDSC (c3)")
                       )
seu.obj$majorID_sub <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "0", "1" = "1", 
                                   "2" = "0", "3" = "0",
                                   "4" = "2", "5" = "3")
                       )
seu.obj$clusID_new <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "#F8766D", "1" = "#B79F00", 
                                   "2" = "#00BA38", "3" = "#00BFC4",
                                   "4" = "#619CFF", "5" = "#F564E3")
                       )
seu.obj$sub_colz <- Idents(seu.obj)


### Supp data - Generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "clusID_new", numOfFeats = 24, outName = outName,
                     outDir = paste0("../output/viln/", outName, "/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )

### Fig 6a - UMAP by clusterID_sub
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusID_new",
              cols = colz.df$colour[c(2,5,4,3)],
              pt.size = 0.5,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, 
                alpha = 0.8, labCol = colz.df$labCol[c(2,5,4,3)]) + NoLegend() + theme(axis.title = element_blank(),
                                                                                                               panel.border = element_blank())
ggsave(paste("../output/", outName, "/", "rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig 6b - skew plot for abundance analysis
p <- skewPlot(seu.obj, groupBy = "majorID_sub", outDir = paste0("../output/", outName), outName = outName)
ggsave(paste0("../output/", outName, "/", "skewPlot.png"), width = 6, height = 4)


### Fig 6c - DGE analysis

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


### Fig 6d - GO GSEA of DEGs
p <- plotGSEA(pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv"),
              geneList = NULL, category = "C5", species = "dog", termsTOplot = 10, upOnly = T, trunkTerm = T,
              pvalueCutoff = 0.05, subcategory = NULL, 
              saveRes = paste0("../output/", outName, "/c5_", outName, "_res.csv")) + theme(axis.title=element_text(size = 16))
p <- p + scale_x_continuous(limits = c(-20,ceiling(max(p$data$x_axis)*1.05)), 
                            breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
                            name = "log10(p.adj)") + ggtitle("Gene ontology") + theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste0("../output/", outName, "/", "gseaPlot_1.png"), width = 7, height = 7)


### Fig 6e - Reactome GSEA of DEGs
p <- plotGSEA(pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv"),
         geneList = NULL, category = "C2", species = "dog", termsTOplot = 10, upOnly = T, trunkTerm = T,
                     pvalueCutoff = 0.05, subcategory = "CP:REACTOME", saveRes = paste0("../output/", outName, "/c2_", outName, "_res.csv")
                    ) + theme(axis.title=element_text(size = 16))
p <- p + scale_x_continuous(limits = c(-20,ceiling(max(p$data$x_axis)*1.05)), breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),name = "log10(p.adj)") + ggtitle("Reactome") + theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste0("../output/", outName, "/", "gseaPlot_2.png"), width = 7, height = 7)


### Fig 6f - Split UMAP of selected DEGs
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$cellSource)))

features <- c("LTF", "IL16","CCR2", "CCL7", "IL1A",  "PTGES","OSM", "CD80","C1QC")

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

### Generate supplemental figure for mono vs macrophage
#load in annotated tumor data
seu.obj.tumor <- readRDS(file = "./output/s3/canine_naive_n6_annotated.rds")

### Select cells to integrate (all pbmcs from OS dogs and only immune cells from tumor)
Idents(seu.obj.tumor) <- "celltype.l3"
table(seu.obj.tumor@active.ident)
seu.obj.tams <- subset(seu.obj.tumor,
                          idents = c("ANGIO_TAM","LA-TAM_C1QC_hi","LA-TAM_SPP2_hi","TAM_ACT","TAM_INT", "CD4+_TIM", "CD4-_TIM", "IFN-TAM")
                         )

seu.obj.bl <- subset(seu.obj, subset = cellSource == "Blood")

seu.sub.list <- c(SplitObject(seu.obj.bl, split.by = "name"),SplitObject(seu.obj.tams, split.by = "name"))


seu.obj <- indReClus(seu.obj = NULL, outDir = "./output/s2/", subName = "20231014_macMono_bloodANDtils", preSub = T, seu.list = seu.sub.list,
                      vars.to.regress = "percent.mt",nfeatures = 2000
                       )


clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "20231014_macMono_bloodANDtils", test_dims = c(40,35,30), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "20231014_macMono_bloodANDtils", final.dims = 35, final.res = 0.3, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.5, n.neighbors = 60, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                      )

seu.obj$macMono <- ifelse(!is.na(seu.obj$cellSource), "Monocyte", ifelse(seu.obj$celltype.l3 == "CD4+_TIM" | seu.obj$celltype.l3 == "CD4-_TIM", "TIM", "TAM"))

#complete hc
seu.obj$type <- paste0(seu.obj$macMono,"_",seu.obj$name)

metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@data)) #use log noralized count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL                                                                       

M <- (1- cor(t(clusAvg_expression),method="pearson"))/2
hc <- hclust(as.dist(M),method="complete")

ggtree(as.phylo(hc)) + geom_tiplab(offset = 0.003) + xlim(NA,0.1) #+ geom_tippoint(shape = 21,size = 8,alpha = 1, colour="black", fill = "white") + geom_tiplab(aes(label = seq(1:43)-1,colour=c("black"),offset = -0.001))
ggsave(paste("./output/", outName, "/", subName, "/",subName, "_hc.png", sep = ""), width = 8, height = 8)

# Transpose count matrix and calculate distances using dist()
sampleDists <- clusAvg_expression %>% dist()

# Convert distance dataframe to matrix
sampleDistMatrix <- as.matrix(sampleDists)

# Choose continuous palette from RColorBrewer
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Plot sample distance heatmap with ComplexHeatmap
png(file = paste0("./output/", outName, "/", subName, "/",subName, "_samDist.png"), width=4000, height=4000, res=400)
par(mfcol=c(1,1))         

ht <- Heatmap(sampleDistMatrix, #name = "mat", #col = col_fun,
              name = "Sample distance",
              cluster_rows = hc,
              row_title = "",
              row_title_gp = gpar(fontsize = 24),
              col=colors,#viridis(option = "magma",100),
              cluster_columns = hc,
              column_title = "",
              column_title_gp = gpar(fontsize = 24),
              column_title_side = "bottom",
              column_names_rot = 45,
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
              
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                   grid.text(sprintf("%.0f", small_mat[i, j]), x, y, gp = gpar(fontsize = 10))
#               }
             )

draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),heatmap_legend_side = "top")

dev.off()

########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end monocyte analysis   ######## <<<<<<<<<<<<<<
########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

