#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")

### Analysis note: 

################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin CD4 T cell preprocessing   ######## <<<<<<<<<<<<<<
################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set output params
outName <- "cd4"

#subset on the CD4 T cells & complete subset analysis
table(seu.obj$conSense)[grepl("CD4",names(table(seu.obj$conSense)))]
Idents(seu.obj) <- "conSense"
seu.obj.sub <- subset(seu.obj, 
                          ident = c("CD4+ Naive", "CD4+ T reg", "CD4+ TCM", "CD4+ TEM", 
                                    "CD4+ TEM, Th1-like","CD4+ TEM, Th17-like", "CD4+ TEM, Th2-like",
                                   "CD4_act","CD4_ex","CD4_naive","CD4_reg")
                         )

table(seu.obj.sub$conSense)
min(table(seu.obj.sub$orig.ident)) > 100
seu.sub.list <- SplitObject(seu.obj.sub, split.by = "orig.ident")

seu.obj <- indReClus(seu.obj = NULL, outDir = "../output/s2/", subName = "20230507_cd4_bloodANDtils", preSub = T, seu.list = seu.sub.list,
                      vars.to.regress = "percent.mt",nfeatures = 2000
                       )

clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = "20230507_cd4_bloodANDtils", test_dims = 40, algorithm = 3, prefix = "integrated_snn_res.")
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "20240305_cd4_bloodANDtils", final.dims = 40, final.res = 0.4, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.6, n.neighbors = 75, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                      )


############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin CD4 T cell analysis   ######## <<<<<<<<<<<<<<
############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in processed data
seu.obj <- readRDS("../output/s3/20230507_cd4_bloodANDtils_res0.7_dims40_dist0.6_neigh75_S3.rds")
seu.obj$cellSource <- ifelse(grepl("tils",seu.obj$name),"TILs","Blood")
outName <- "cd4"

#fix nomencalture error
seu.obj$celltype.l3 <- as.factor(seu.obj$celltype.l3)
levels(seu.obj$celltype.l3)[2] <- "CD4_fh"

tmpColz <- gg_color_hue(8)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "1", "1" = "2", 
                                   "2" = "0", "3" = "0",
                                   "4" = "3", "5" = "4",
                                   "6" = "5", "7" = "0")
                       )
seu.obj$clusID_new <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "TCM (c1)", "1" = "Treg/Tfh (c2)", 
                                   "2" = "Naive (c0)", "3" = "Naive (c0)",
                                   "4" = "TEM_Th2-like (c3)", "5" = "TEM (c4)",
                                   "6" = "TEM_Th1-like (c5)", "7" = "Naive (c0)")
                       )
seu.obj$majorID_sub <- Idents(seu.obj)
seu.obj$majorID_sub <- factor(seu.obj$majorID_sub, levels = c("TCM (c1)", "Treg/Tfh (c2)", 
                                                              "Naive (c0)",
                                                              "TEM_Th2-like (c3)", "TEM (c4)",
                                                              "TEM_Th1-like (c5)"
                                                             )[c(3,1,5,6,4,2)])

Idents(seu.obj) <- "clusID_new"
seu.obj <- RenameIdents(seu.obj, c("0" = "#009DA5", "1" = "#148B7D", 
                                   "2" = "#0E3F5C", "3" = "#D4F3A3",
                                   "4" = "#9CD6BA", "5" = "#2E4F79")
                       )
seu.obj$sub_colz <- Idents(seu.obj)


#generate violin plots for each cluster
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = subName,
                     outDir = paste0("../output/", outName, "/viln/"), outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS")
                    )


### Fig 2a - UMAP by clusID_new
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusID_new",
              cols = levels(seu.obj$sub_colz),
              pt.size = 0.5,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, 
                alpha = 0.8, labCol = c("black","black","white","black","black","white")) + NoLegend() + theme(axis.title = element_blank(),
                                                                                                               panel.border = element_blank())
ggsave(paste("../output/", outName, "/", "rawUMAP.png", sep = ""), width = 7, height = 7)


### Supp Fig 2a - split UMAP by cell source with original labels
seu.obj$ct <- ifelse(seu.obj$cellSource == "TILs", paste0("TIL_", seu.obj$celltype.l3), paste0("Blood_",seu.obj$celltype.l3_pbmc))
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "ct",
              split.by = "cellSource",
              pt.size = 0.5,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi)
ggsave(paste("../output/", outName, "/", "supp_label.png", sep = ""), width = 16, height = 7)

#check consenseous
table(seu.obj$ct, seu.obj$clusterID_sub) %>% sweep(., 2, colSums(.), `/`) %>% melt() %>% group_by(Var.2) %>% filter(value == max(value))


### Supp Fig 2b - CXCL13 expression by CD4 T cells
seu.obj$ct <- paste0(ifelse(seu.obj$cellSource == "TILs", "TIL_","Blood_"), seu.obj$majorID_sub)

cxcl13.df <- FetchData(seu.obj, vars = c("cellSource", "CXCL13"))
cxcl13.df <- cxcl13.df[cxcl13.df$CXCL13 > 0, ]
cxcl13.df <- cxcl13.df %>% group_by(cellSource) %>% summarize(N = n())
cxcl13.df$pct <- c(cxcl13.df$N[1]/table(seu.obj$cellSource)[1]*100, cxcl13.df$N[2]/table(seu.obj$cellSource)[2]*100)
cxcl13.df$label <- c(paste0(cxcl13.df$N[1],"/",table(seu.obj$cellSource)[1]), paste0(cxcl13.df$N[2],"/",table(seu.obj$cellSource)[2]))

p1 <- ggplot(cxcl13.df, aes(x = cellSource, y = pct)) +
stat_summary(fun = mean, geom = "bar", fill = "grey", width = 0.7) +
geom_text(aes(label=label), vjust=-0.25) + 
scale_y_continuous(limits = c(0, 10), expand = expansion(mult = c(0, 0))) + 
theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.line = element_line(colour="black"),
    legend.key = element_rect(colour = NA, fill = NA),
    axis.ticks = element_line(colour="black"),
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(angle = 90, vjust = 2),
    axis.title.x = element_blank(),
    axis.text = element_text(face = "bold"),
      legend.background = element_rect(colour = "transparent", fill = NA),
    legend.box.background = element_rect(colour = "transparent", fill = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.background = element_blank()
  ) + labs(y = "% CXCL13+ cells") 
ggsave(paste("../output/", outName, "/", "barchart.png", sep = ""), width = 3, height = 2)


# split cluster map
set.seed(12)
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$cellSource)))

pi <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              split.by = "cellSource",
              pt.size = 0.5,
              label = F,
              label.box = F,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste("../output/", outName, "/", "ID_UMAP.png", sep = ""), width = 8, height = 5)


### Fig 2b - skew plot for abundance analysis
p <- skewPlot(seu.obj)
ggsave(paste0("../output/", outName, "/", "skewPlot.png"), width = 6, height = 4)


### Fig 2c - DGE analysis

#load in the tumor and pal signatures to exlude from DE analysis
pal_feats = c('TIMP1', 'NAA10', 'ENSCAFG00000037735', 'GP6', 'SEC11C', 'FTL', 'NRGN', 'ACOT7', 'VCL', 'RSU1', 'ITGB1', 'H3-3A', 'RABGAP1L', 'SELP', 'SH3GLB1', 'ACTB', 'ENSCAFG00000008221', 'TLN1', 'GSN', 'AMD1', 'TREM2', 'SH3BGRL2', 'MYH9', 'PLEK', 'ENSCAFG00000042554', 'RAP1B', 'ENSCAFG00000004260', 'NAP1L1', 'PPBP', 'RASA3', 'ITGA2B', 'EIF1', 'ACTG1', 'C9H17orf64', 'JMJD6', 'CCL14', 'GNG11', 'IGF2BP3', 'TBXAS1', 'VDAC3', 'MARCHF2', 'TPM4', 'TKT', 'FTH1.1', 'FERMT3', 'RTN3', 'PRKAR2B', 'SVIP', 'ENSCAFG00000030286', 'ADA', 'MYL9', 'TUBB1', 'TUBA1B', 'METTL7A', 'THBS1', 'SERF2', 'PIF1', 'B2M', 'GAS2L1', 'YWHAH', 'HPSE', 'ATG3', 'ENSCAFG00000015217', 'ITGA6','RGS18', 'SUB1', 'LGALS1', 'CFL1', 'BIN2', 'CAT', 'RGS10', 'MGST3', 'TMBIM6', 'PFN1', 'CD63', 'RALBP1', 'GNAS', 'SEPTIN7', 'TPT1', 'UBB', 'ATF4', 'BBLN', 'MTDH', 'ENSCAFG00000017655','FYB1', 'ENO1', 'GABARAP', 'SSR4', 'MSN', 'ENSCAFG00000011134', 'ENSCAFG00000046637', 'COX8A', 'DLA-64', 'CD47', 'VASP', 'DYNLRB1', 'DLA88', 'SMDT1', 'ATP5PF','ELOB', 'ENSCAFG00000029155', 'ARPC3', 'VPS28', 'LRRFIP1', 'SRP14', 'ABRACL', 'ENSCAFG00000043577', 'ENSCAFG00000042598')
tumor.sig <- 

createPB(seu.obj = seu.obj, groupBy = "allCells", comp = "cellSource", biologicalRep = "name",
         outDir = paste0("../output/", outName, "/pseudoBulk/"), 
         grepTerm = "tils", grepLabel = c("TILs", "Blood"), featsTOexclude = c(pal_feats,tumor.sig), lowFilter = T, dwnSam = F
        )
p_volc <- pseudoDEG(metaPWD = paste0("../output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"),
          padj_cutoff = 0.01, lfcCut = 0.58, outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = subName, idents.1_NAME = "TILs", idents.2_NAME = "Blood",
          inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "TILS vs Blood", fromFile = T, meta = NULL, pbj = NULL, returnVolc = T, paired = F, pairBy = "", 
          minimalOuts = F, saveSigRes = T, filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T
                     )
p <- prettyVolc(plot = p_volc[[1]], rightLab = NULL, leftLab = NULL, rightCol = "red", leftCol = "blue", arrowz = F
                    ) + theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
                              axis.line = element_blank())
ggsave(paste0("../output/", outName, "/", "volcPlot.png"), width =7, height = 7)


### Fig 2d - GO GSEA of DEGs
p <- plotGSEA(pwdTOgeneList = "../output/cd4/pseudoBulk/allCells/cd4_cluster_allCells_all_genes.csv",
              geneList = NULL, category = "C5", species = "dog", termsTOplot = 10, upOnly = T, 
              pvalueCutoff = 0.05, subcategory = NULL, 
              saveRes = paste0("../output/", outName, "/gsea_res/c5_cd4_res.csv")) + theme(axis.title=element_text(size = 16))
p <- p + scale_x_continuous(limits = c(-12,ceiling(max(p$data$x_axis)*1.05)), 
                            breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
                            name = "log10(p.adj)") + ggtitle("Gene ontology") + theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste0("../output/", outName, "/", "gseaPlot_1.png"), width =7, height = 7)


### Fig 2e - Reactome GSEA of DEGs
p <- plotGSEA(pwdTOgeneList = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/tils_naive6_w_all_OS_PBMC/cd4/pseudoBulk/allCells/cd4_cluster_allCells_all_genes.csv",
         geneList = NULL, category = "C2", species = "dog", termsTOplot = 10, upOnly = T,
                     pvalueCutoff = 0.05, subcategory = "CP:REACTOME", saveRes = paste0("../output/", outName, "/gsea_res/c2_cd4_res.csv")
                    ) + theme(axis.title=element_text(size = 16))
p <- p + scale_x_continuous(limits = c(-12,ceiling(max(p$data$x_axis)*1.05)), breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),name = "log10(p.adj)") + ggtitle("Reactome") + theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste0("../output/", outName, "/", "gseaPlot_2.png"), width =7, height = 7)

### Fig 2f - Split UMAP of selected DEGs
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$cellSource)))

features <- c("SELL", "LEF1", "TMEM154", "TNFRSF18", "TNFRSF4", "HAVCR1", "CXCL13", "LAG3","IL4I1")

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


########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end CD4 T cell analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


