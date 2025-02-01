#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")

### Analysis note: 

################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin CD8 T cell preprocessing   ######## <<<<<<<<<<<<<<
################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set output params
outName <- "cd8"

#subset on the CD8 T cells & complete subset analysis
table(seu.obj$conSense)[grepl("CD8",names(table(seu.obj$conSense)))]
Idents(seu.obj) <- "conSense"
seu.obj.sub <- subset(
    seu.obj, 
    ident = c(
        "CD8_SPP1_hi", "CD8+ Memory","CD8+ Effector", "CD8+ gd T cell",
        "CD8+ Naive", "NK cell", "NK T cell", "NK cell","CD8_eff",
        "CD8_ex", "NK"
    )
)

table(seu.obj.sub$conSense)
min(table(seu.obj.sub$orig.ident)) > 100
seu.sub.list <- SplitObject(seu.obj.sub, split.by = "orig.ident")

seu.obj <- indReClus(seu.obj = NULL, outDir = "../output/s2/", subName = "20230507_cd8Friends_bloodANDtils", preSub = T, seu.list = seu.sub.list,
                      vars.to.regress = "percent.mt",nfeatures = 2000
                       )

clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = "20230507_cd8Friends_bloodANDtils", test_dims = c(40,35,30), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "20230507_cd8Friends_bloodANDtils", final.dims = 40, final.res = 0.4, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                      )

#remove suspect cell clusters -- they were expressing CSF1
Idents(seu.obj) <- "clusterID_sub"
seu.obj.sub <- subset(seu.obj, invert = T,
                          ident = c(3,4,5)
                         )
table(seu.obj.sub$clusterID_sub)
DefaultAssay(seu.obj.sub) <- "integrated"
clusTree(seu.obj = seu.obj.sub, dout = "../output/clustree/", outName = "20230507_cd8Friends_bloodANDtils", test_dims = c(40), algorithm = 3, prefix = "integrated_snn_res.")


seu.obj <- dataVisUMAP(seu.obj = seu.obj.sub, outDir = "../output/s3/", outName = "20230507_cd8Friends_bloodANDtils", final.dims = 40, final.res = 0.6, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                      )

############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin CD8 T cell analysis   ######## <<<<<<<<<<<<<<
############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in processed data
seu.obj <- readRDS("../output/s3/20230507_cd8Friends_bloodANDtils_res0.6_dims40_dist0.3_neigh30_S3.rds")
seu.obj$cellSource <- ifelse(grepl("tils",seu.obj$name),"Tumor","Blood")
outName <- "cd8"
#Check consenseous
seu.obj$ct <- ifelse(
    seu.obj$cellSource == "Tumor", 
    paste0("TIL;", seu.obj$celltype.l3), 
    paste0("Blood;",seu.obj$celltype.l3_pbmc)
)
#Plot % in blood
table(seu.obj$ct, seu.obj$clusterID_sub) %>% 
    melt() %>% separate(Var.1, sep = ";", c("source", "ct")) %>% 
    filter(source == "Blood") %>%
    group_by(Var.2) %>% 
    mutate(pct = value/sum(value)) %>% 
    filter(value == max(value))
#Plot % in TIL
table(seu.obj$ct, seu.obj$clusterID_sub) %>% 
    melt() %>% 
    separate(Var.1, sep = ";", c("source", "ct")) %>% 
    filter(source == "TIL") %>% 
    group_by(Var.2) %>% 
    mutate(pct = value/sum(value)) %>% 
    filter(value == max(value))
#Plot % in both
Idents(seu.obj) <- "cellSource"
seu.obj.ds <- subset(seu.obj, downsample = min(table(seu.obj$cellSource)))
table(seu.obj.ds$ct, seu.obj.ds$clusterID_sub) %>% 
    sweep(., 2, colSums(.), `/`) %>% 
    melt() %>% separate(Var.1, sep = ";", c("source", "ct")) %>% 
    group_by(source, Var.2) %>% 
    filter(value == max(value))

#set metadata
colz.df <- read.csv("./metaData/majorGroups.csv")
colz.df <- colz.df[colz.df$majorID2 == "cyto", ]
colz.df$colour

table(seu.obj$clusterID_sub,seu.obj$conSense)

tmpColz <- gg_color_hue(5)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "#0066A5", "1" = "#99BFEF", 
                                   "2" = "#D0E4FF", "3" = "#B3B7CF",
                                   "4" = "#3267AD") #5D6CAE
                       )
seu.obj$dcColz <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "CD8_eff_1 (c0)", "1" = "CD8_mem_1 (c1)", 
                                   "2" = "CD8_eff_2 (c2)", "3" = "CD8_naive (c3)",
                                   "4" = "CD8_mem_2 (c4)")
                       )
seu.obj$majorID_sub <- Idents(seu.obj)
seu.obj$majorID_sub <- factor(
    seu.obj$majorID_sub, 
    levels = c("CD8_eff_1 (c0)", "CD8_mem_1 (c1)", "CD8_eff_2 (c2)",
               "CD8_naive (c3)", "CD8_mem_2 (c4)")[c(4,1,3,2,5)]
)
#Export annotations
write.csv(seu.obj@meta.data["majorID_sub"], file = paste0("../output/annotations/", outName, ".csv"))


### Supp data - Generate violin plots for each cluster
vilnPlots(
    seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, 
    outName = outName, outDir = paste0("../output/viln/", outName, "/"), 
    outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), 
    returnViln = F
)

### Use gene expression patterns to further support classifications
### Fig sup - Use enrichment scoring help ID cells
#load in gene lists as a named list
modulez <- list(
    "Naïve" = c("CCR7", "LEF1", "SELL", "TCF7"),
    "Effector" = c("IL2RA", "TNFRSF8", "CD69", "TNFRSF4", "ICOS", "KLRG1", 
                 "HAVCR2", "TBX21", "IFNG", "IL2", "PRF1", "GZMB", "GZMA", 
                 "CCL3", "CCL4", "CCL5"), 
    "CD8 REG" = c("KLRG1", "PDCD1", "LAG3", "FOXP3", "EGR2"), 
    "CD8 REG-MEM" = c("CD69", "CD101", "PDCD1", "CXCR6", "CCR8", "ITGA1", 
                      "ITGAE", "SELPLG", "RUNX3"), 
    "CD8 TCM" = c("CD27", "CD28", "CXCR3", "TBX21", "EOMES", "IFNG", "IL2"), 
    "CD8 TEM" = c("CD44", "KLRG1", "EOMES", "TBX21", "GZMK", "IFNG", "PRF1", 
                  "IL2"),
    "Exhausted" = c("BTLA", "CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT"),
    "Costim" = c("ICOS", "CD226", "SLAMF1", "TNFRSF14", "TNFRSF25", "TNFRSF9")
)

#run module score
seu.obj <- AddModuleScore(
    seu.obj, features = modulez, name = "_score"
)
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)
#plot the results of enrichment scores
features <- names(modulez)
ecScores <- majorDot(
    seu.obj = seu.obj, groupBy = "clusterID_sub", features = rev(features)
) + 
    coord_flip() + 
    theme(
        plot.margin = margin(3, 0, 3, 0, "pt"),
        axis.text.y=element_text(size=10),
        axis.title = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        axis.text.x = element_text(angle=0, hjust = 0.5)
    ) + 
    scale_y_discrete(position = "right") + 
    scale_colour_continuous(name="Enrichment score", type = "viridis")
ggsave(paste("../output/", outName, "/", outName, "_dotPlot_ecScores.png", sep = ""), width = 6, height = 4)

#plot indivdual members of each term
modulez <- c(list("Module score" = names(modulez)), modulez)
labelz <- as.data.frame(names(modulez))
colnames(labelz) <- "labz"
labelz$modLen <- unname(unlist(lapply(modulez, length)))
cntr <- 0
plots <- lapply(modulez, function(x){
    cntr <<- cntr+1
    labz.df <- labelz[cntr,]

    majorDot(
        seu.obj = seu.obj, groupBy = "clusterID_sub", features = rev(unname(unlist(x)))
    ) + 
        theme(
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "right",
            legend.direction = "vertical",
            axis.title = element_blank(),
            plot.margin = margin(3, 0, 3, 0, "pt")
        ) + 
        scale_colour_distiller(palette = "RdYlBu", name='Average\nexpression', limits = c(-2.5,2.5)) + 
        geom_text(
            data = labz.df, aes(label = labz, y = 5.85, x = (modLen+1)/2),
            angle = 270, vjust = 0.5, hjust=0.5, size = 12*0.36
        ) + 
        coord_flip(ylim = c(1, 5.75), clip = "off") + 
        annotate(
            "segment", x = -Inf, y = 5.5, xend = Inf, yend = 5.5, 
            lineend = "round", linejoin = "bevel", linetype ="solid", 
            colour = "grey70", alpha = 0.7, size = 0.5
        )
})

#get all the plots together
    patch <- area()
    nrow <- length(modulez)
    ncol <- 1
    counter=0
    for (i in 1:nrow) {
        for (x in 1:ncol) {
                patch <- append(patch, area(t = i, l = x, b = i, r = x))
        }
    }

#change to the color of the module scores for visual distinction & plot final
plots$`Module score` <- plots$`Module score` + 
    theme(axis.text.x = element_text(angle=0, hjust = 0.5)) + 
    scale_y_discrete(position = "right") + 
    scale_colour_viridis() + 
    guides(color = guide_colorbar(title = 'Module\nscore'), limits = c(-2.5,2.5))
p <- Reduce( `+`, plots ) +
    plot_layout(
        guides = "collect", design = patch, 
        height = unname(unlist(lapply(modulez, length)))/sum(unname(unlist(lapply(modulez, length))))
    )
ggsave(paste("../output/", outName, "/", outName, "_dotPlot_ecScores_2.png", sep = ""), width = 2.1, height = 7.4, scale = 2)

### Fig supp - Reviewer requested support for exhaustion claim
p <- VlnPlot(
    seu.obj, group.by = "majorID_sub", 
    features = "Exhausted", 
    split.by = "cellSource"
)
p <- p +
    theme(
        axis.title = element_blank()
    ) +
    ggtitle("Exhaustion enrichment score")
ggsave(paste0("../output/", outName, "/", outName, "viln.png"), width = 6, height = 3, scale = 1.25)
#run stats
stat_df <- seu.obj@meta.data %>%
    select(name, majorID_sub, cellSource, Exhausted) %>%
    group_by(name, majorID_sub) %>%
    mutate(
        sum_score = mean(Exhausted)
    ) %>%
    ungroup() %>%
    distinct(sum_score, .keep_all = TRUE)

statz <- compare_means(sum_score ~ cellSource, group.by = "majorID_sub", stat_df) %>%
    select(-p.format, -.y.) %>%
    mutate(
        p.signif = symnum(p.adj, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                          symbols = c("****", "***", "**", "*", "ns"))
    )
stat_df %>%
    group_by(cellSource, majorID_sub) %>%
    summarize(MEAN = mean(sum_score))

### Fig 3a - UMAP by clusterID_sub
pi <- DimPlot(
    seu.obj, 
    reduction = "umap", 
    group.by = "clusterID_sub",
    cols = levels(seu.obj$dcColz),
    pt.size = 0.5,
    label = TRUE,
    label.box = TRUE,
    shuffle = TRUE
)
pi <- cusLabels(
    plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8, 
    labCol = c("white","black","black","black","white")
) + 
    NoLegend() + 
    theme(
            axis.title = element_blank(),
            panel.border = element_blank()
    )
ggsave(paste("../output/", outName, "/", "rawUMAP.png", sep = ""), width = 7, height = 7)

### Supp Fig 3a - UMAP highlighting NK cells
seu.obj$nk <- ifelse(seu.obj$celltype.l3_pbmc == "NK cell" | 
                     seu.obj$celltype.l3 == "NK", "NK", "Other")
Idents(seu.obj) <- "nk"
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              cells.highlight = WhichCells(seu.obj, idents = "NK"),
              cols.highlight = "#f22e9d",
              pt.size = 0.5,
              label = F,
              label.box = F,
              shuffle = F
)
pi <- formatUMAP(plot = pi)
ggsave(paste("../output/", outName, "/", "highlight_NK.png", sep = ""), width = 7, height = 7)


### Fig 3b - skew plot for abundance analysis
p <- skewPlot(
    seu.obj, groupBy = "majorID_sub", yAxisLabel = "Percent of CD8 T cells",
    dout = paste0("../output/", outName), outName = outName, 
    sampleRep = "name", grepTerm = "tils", grepRes = c("Tumor","Blood")
)
ggsave(paste0("../output/", outName, "/", "skewPlot.png"), width = 6, height = 4)

### Fig 3c - DGE analysis
#load in the tumor and pal signatures to exlude from DE analysis
pal_feats = c(
    "TIMP1", "NAA10", "ENSCAFG00000037735", "GP6", "SEC11C", "FTL", 
    "NRGN", "ACOT7", "VCL", "RSU1", "ITGB1", "H3-3A", "RABGAP1L", 
    "SELP", "SH3GLB1", "ACTB", "ENSCAFG00000008221", "TLN1", "GSN", 
    "AMD1", "TREM2", "SH3BGRL2", "MYH9", "PLEK", "ENSCAFG00000042554", 
    "RAP1B", "ENSCAFG00000004260", "NAP1L1", "PPBP", "RASA3", "ITGA2B", 
    "EIF1", "ACTG1", "C9H17orf64", "JMJD6", "CCL14", "GNG11", "IGF2BP3", 
    "TBXAS1", "VDAC3", "MARCHF2", "TPM4", "TKT", "FTH1.1", "FERMT3", 
    "RTN3", "PRKAR2B", "SVIP", "ENSCAFG00000030286", "ADA", "MYL9", 
    "TUBB1", "TUBA1B", "METTL7A", "THBS1", "SERF2", "PIF1", "B2M", 
    "GAS2L1", "YWHAH", "HPSE", "ATG3", "ENSCAFG00000015217", "ITGA6", 
    "RGS18", "SUB1", "LGALS1", "CFL1", "BIN2", "CAT", "RGS10", "MGST3", 
    "TMBIM6", "PFN1", "CD63", "RALBP1", "GNAS", "SEPTIN7", "TPT1", 
    "UBB", "ATF4", "BBLN", "MTDH", "ENSCAFG00000017655", "FYB1", 
    "ENO1", "GABARAP", "SSR4", "MSN", "ENSCAFG00000011134", "ENSCAFG00000046637", 
    "COX8A", "DLA-64", "CD47", "VASP", "DYNLRB1", "DLA88", "SMDT1", 
    "ATP5PF", "ELOB", "ENSCAFG00000029155", "ARPC3", "VPS28", "LRRFIP1", 
    "SRP14", "ABRACL", "ENSCAFG00000043577", "ENSCAFG00000042598"
)
tumor.sig <- read.csv("./metaData/tumorSig.csv", header = T)$x

createPB(
    seu.obj = seu.obj, groupBy = "allCells", 
    comp = "cellSource", biologicalRep = "name",
    outDir = paste0("../output/", outName, "/pseudoBulk/"), 
    grepTerm = "tils", grepLabel = c("Tumor", "Blood"), 
    featsTOexclude = c(pal_feats,tumor.sig), 
    lowFilter = T, dwnSam = F
)
p_volc <- pseudoDEG(
    inDir = paste0("../output/", outName, "/pseudoBulk/"), 
    metaPWD = paste0("../output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"),
    outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = outName,
    padj_cutoff = 0.01, lfcCut = 0.58, strict_lfc = T,
    idents.1_NAME = "Tumor", idents.2_NAME = "Blood", title = "TILS vs Blood", 
    fromFile = T, returnVolc = T, filterTerm = "^ENSCAF", mkDir = T
)
p <- prettyVolc(
    plot = p_volc[[1]], rightLab = "Tumor", leftLab = "Blood", 
    rightCol = "red", leftCol = "blue", arrowz = T
) + 
    theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_blank(),
        legend.position = c(0.10, 0.85),
    )
ggsave(paste0("../output/", outName, "/", "volcPlot.png"), width = 7, height = 7)

### Fig 3d - GO GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", 
                           outName, "_cluster_allCells_all_genes.csv"),
    geneList = NULL, category = "C5", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = NULL,
    lolli = T,
    saveRes = paste0("../output/", outName, "/c5_", outName, "_res.csv")
) + 
    theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)
    )
p <- p + scale_x_continuous(
    limits = c(-22,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Gene ontology") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.07),
        legend.background = element_rect(fill = 'transparent', colour = NA),
        legend.key = element_rect(fill = 'transparent', colour = NA)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_1.png"), width = 7, height = 7)

### Fig 3e - Reactome GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", 
                           outName, "_cluster_allCells_all_genes.csv"),
    geneList = NULL, category = "C2", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "CP:REACTOME", 
    lolli = T,
    saveRes = paste0("../output/", outName, "/c2_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-8,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Reactome") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_2.png"), width = 7, height = 7)

### Fig 2e - ImmuneSigDB GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", 
                           outName, "_cluster_allCells_all_genes.csv"),
    geneList = NULL, category = "C7", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
    lolli = T, filterTerm = "_CD8_",
    saveRes = paste0("../output/", outName, "/c7_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-30,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("ImmuneSigDB") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_3.png"), width = 7, height = 7)


### Fig 2f - Split UMAP of selected DEGs
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$cellSource)))

features <- c("TRIM22", "LEF1", "PTGDR", "CD53", "CX3CR1", "HAVCR2", "TNFSF9", "LAG3","TNFRSF18")
p <- FeaturePlot(
    seu.obj.sub, features = features, pt.size = 0.1, split.by = "cellSource", 
    order = T, cols = c("lightgrey","darkblue"), by.col = F
) + 
    labs(x = "UMAP1", y = "UMAP2") & 
    theme(
        axis.text = element_blank(),
        axis.title.y.right = element_text(size = 11),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(size=11),
        title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) 
ggsave(paste("../output/", outName, "/", "splitFeats.png", sep = ""), width = 12, height = 4)

### Fig sup - Split dotplot of selected DEGs
res.df <- read.csv(paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv"))
res.df <- res.df[!grepl("^ENS", res.df$gene), ]
geneList_UP <- res.df %>% filter(padj < 0.1) %>% filter(log2FoldChange > 1) %>% pull(gene)
geneList_DWN <- res.df %>% filter(padj < 0.1) %>% filter(log2FoldChange < -1) %>% pull(gene)

seu.obj$cellSource <- as.factor(seu.obj$cellSource)
p <- splitDot(
    seu.obj = seu.obj, groupBy = "majorID_sub", splitBy = "cellSource", buffer = 125,
    namedColz = setNames(c("#F8766D", "#00BFC4"),  c("Blood", "Tumor")), 
    geneList_UP = c(geneList_UP[1:20]), geneList_DWN = c(geneList_DWN[1:20], "CD53"), 
    geneColz = c("red", "blue")
)
p <- p +
    theme(
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.justification = 'center'
    )
ggsave(plot = p, paste0("../output/", outName, "/", outName, "_splitDot.png"), width = 12, height = 5)


# Complete DE analysis in each subcluser
createPB(
    seu.obj = seu.obj, groupBy = "majorID_sub", comp = "cellSource", 
    biologicalRep = "name", outDir = paste0("../output/", outName, "/pseudoBulk/"), 
    min.cell = 10, grepTerm = "tils", grepLabel = c("Tumor", "Blood"), 
    featsTOexclude = c(pal_feats,tumor.sig), lowFilter = T, dwnSam = F
)
pseudoDEG(
    inDir = paste0("../output/", outName, "/pseudoBulk/"), 
    metaPWD = paste0("../output/", outName, "/pseudoBulk/",
                     "majorID_sub_deg_metaData.csv"),
    outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = outName,
    padj_cutoff = 0.01, lfcCut = 0.58,  idents.1_NAME = "Tumor", 
    idents.2_NAME = "Blood", title = "Tumor vs Blood", fromFile = T, 
    returnVolc = F, filterTerm = "^ENSCAF", mkDir = T
)


### Supp fig xx -- heatmap of sig DEGs
# Load in the cluster DE results
files <- paste0("../output/", outName, "/pseudoBulk/", 
                levels(seu.obj$majorID_sub), "/", outName, "_cluster_", 
                levels(seu.obj$majorID_sub), "_all_genes.csv")

df.list <- lapply(files, read.csv, header = T)
res.df <- do.call(rbind, df.list)
# Convert to matrix
cnts_mat <- res.df  %>% 
    mutate(
        direction = ifelse(log2FoldChange > 0, "Up", "Down")
    ) %>% 
    group_by(gs_base, direction) %>% 
    summarize(nRow = n()) %>% 
    pivot_wider(names_from = gs_base, values_from = nRow) %>% 
    as.matrix() %>% 
    t()
colnames(cnts_mat) <- cnts_mat[1,]
cnts_mat <- cnts_mat[-c(1),]
class(cnts_mat) <- "numeric"

# Order by number of total # of DEGs
cnts_mat[is.na(cnts_mat)] <- 0
orderList <- rev(rownames(cnts_mat)[order(rowSums(cnts_mat))])
cnts_mat <- cnts_mat[match(orderList, rownames(cnts_mat)),]        
anno_mat <- cnts_mat
cnts_mat[,1] <- -cnts_mat[,1]

# Generate plot
ht <- Heatmap(
    cnts_mat,#name = "mat", #col = col_fun,
    name = "# of DEGs",
    cluster_rows = F,
    row_title = "Cell type",
    col = circlize::colorRamp2(c(min(cnts_mat), 0,max(cnts_mat)), colors = c("blue", "white", "red")),
    cluster_columns = F,
    show_column_names = TRUE,
    column_title_side = "top",
    column_names_rot = 0,
    column_names_gp = gpar(fontsize = 14, col = "black"),
    column_names_centered = TRUE,
    heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                              labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if(cnts_mat[i, j] < -40) {
          grid.text(sprintf("%.0f", as.matrix(anno_mat)[i, j]), x, y, gp = gpar(fontsize = 14, col = "white"))
      } else if(cnts_mat[i, j] > -40) {
          grid.text(sprintf("%.0f", as.matrix(anno_mat)[i, j]), x, y, gp = gpar(fontsize = 14, col = "black"))
      }
    })
png(file = paste0("../output/", outName, "/", outName, "_deg_heatmap.png"), width=1500, height=2000, res=400)
par(mfcol=c(1,1))     
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),show_heatmap_legend = FALSE)
dev.off()

clus_colz <- levels(seu.obj$dcColz)
names(clus_colz) <- levels(seu.obj$majorID_sub)
cond_colz <- gg_color_hue(2)
names(cond_colz) <- c("Blood","Tumor")

genez <- res.df %>% 
    filter(!grepl("^ENS", gene)) %>%
    mutate(
        grp = paste0(gs_base, ifelse(log2FoldChange > 0, "_UP", "_DNN"))
    ) %>%
    group_by(grp) %>%
    top_n(-15, padj) %>%
    pull(gene)
res.df <- res.df[res.df$gene %in% c(genez, features), ]

res.df$gs_base <- toupper(res.df$gs_base)
ht <- sigDEG_heatmap(
    seu.obj = seu.obj, groupBy = "majorID_sub", splitBy = "cellSource", forceCleanPlot = T, 
    dge_res = res.df, lfc_thres = 1, cond_colz = cond_colz, clus_colz = clus_colz,
    font_colz = c("white", rep("black", 3), "white"),
    saveName = paste0("../output/", outName, "/", "splitHeat.png"),
    ht_height = 4750, ht_width = 3000
)

### Fig supp - GSEA for DEGs in each cluster
# Run GSEA analysis for Reactome and GO
gsea_helper <- list(c("C5", NULL), c("C2", "CP:REACTOME"))
lapply(gsea_helper, function(db){
    lapply(levels(seu.obj$majorID_sub), function(group){
        group <- gsub("/", "-", group)
        if(is.na(db[2])){
            subcategory <- NULL
        } else {
            subcategory <- db[2]
        }
        p <- plotGSEA(
            pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/", group, "/",
                                   outName, "_cluster_", group, "_all_genes.csv"),
            geneList = NULL, category = db[1], species = "dog", termsTOplot = 10, 
            upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, 
            subcategory = subcategory, 
            lolli = T,
            saveRes = paste0("../output/", outName, "/", db[1], "_", outName, 
                             "_", group, "_res.csv")
        ) + theme(
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 12)
        )
        if(! is.null(p$data)){
            p <- p + scale_x_continuous(
                limits = c(-12,ceiling(max(p$data$x_axis)*1.05)), 
                breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
                name = "-log10(p.adj)"
            ) + 
                ggtitle("Reactome") + 
                theme(
                    plot.title = element_text(size = 20, hjust = 0.5),
                    legend.position = c(0.05, 1.05)
                )
            ggsave(paste0("../output/", outName, "/", group, "_gseaPlot_", db[1], ".png"), width = 7, height = 7)
        }
    })
})

#organize gsea results
df.list <- lapply(c("C2", "C5"), function(database){
    res.df.list <- lapply(levels(seu.obj$majorID_sub), function(group){
        group <- gsub("/", "-", group)
        inFile <- paste0("../output/", outName, "/", database, "_",
                                      outName, "_", group, "_res.csv")
        if(file.exists(inFile)){
            df <- read.csv(inFile, header = T)
            df$CellType <- group
            df$Database <- database
            message(paste("Saving", group, "and", database))
            return(df)
        }
    })
    res.df <- do.call(rbind, res.df.list)
    return(res.df)
})
res.df <- do.call(rbind, df.list)

#data cleaning and conversion to matrix
df <- res.df %>% 
    mutate(
        row = row_number(),
        Description = paste0(
            "[", ifelse(Database == "C2", "R", "G"), "] ", Description
        )
    ) %>% 
    select(CellType, Description, x_axis, row) %>%
    pivot_wider(names_from = CellType, values_from = x_axis) %>% 
    as.data.frame()

mat <- aggregate(
    df[ , 3:ncol(df)], by = df[1], 
    FUN = function(x){sum(as.numeric(x), na.rm = TRUE)}
) %>%
    column_to_rownames("Description") %>%
    as.matrix()

maxCol <- colnames(mat)[max.col(abs(mat), ties.method = "first")]
rowOrderPos <- lapply(colnames(mat), function(x){
    posi <- which(maxCol == x)[rev(order(abs(mat[which(maxCol == x), x])))]
    #create df with position and color
    if(length(posi > 0)){
        posi.df <- as.data.frame(matrix(
            c(
                posi,
                ifelse(grepl("\\[R\\]", rownames(mat)[posi]), "grey50", "black")
            ), 
            nrow = length(posi), ncol = 2, byrow = FALSE,
            dimnames = list(1:length(posi), c("posi", "color"))
        ))
    } else {

            posi.df <- data.frame(matrix(ncol = 2, nrow = 0))
            colnames(posi.df) <- c("posi", "color")
    }
    posi.df$posi <- as.numeric(posi.df$posi)
    return(posi.df)
})
position <- unlist(lapply(1:ncol(mat), function(i){na.omit(rowOrderPos[[i]]$posi[1:7])}))
labelz <- rownames(mat)[c(position)]
mat <- mat[do.call(rbind, rowOrderPos)$posi, ]
position <- match(labelz, rownames(mat))

ra_right <- rowAnnotation(key_feats = anno_mark(
    at = position, labels = rownames(mat)[position], 
    labels_gp = gpar(fontsize = 9, col = do.call(rbind, rowOrderPos)$color[position])))

row_split <- unlist(lapply(1:length(rowOrderPos), function(i){rep(i, nrow(rowOrderPos[[i]]))}))
#create the heatmap
ht <- Heatmap(
    mat,
    name = "Signed log10(FDR)",
    col = circlize::colorRamp2(c(0, 4), c("#EEEEEE", "red")),
    cluster_rows = F,
    row_title_gp = gpar(fontsize = 24),
    row_title = NULL,
    row_split = row_split,
    border = TRUE,
    show_row_names = F,
    cluster_columns = F,
    show_column_names = T,
    right_annotation = ra_right,
    column_split = unlist(lapply(colnames(mat), function(x){strsplit(x, "-_-")[[1]][2]})),
    heatmap_legend_param = list(
        direction = "horizontal",
        legend_width = unit(6, "cm")
    ),
    column_title = NULL
)
png(file = paste0("../output/", outName, "/",  outName, "_heatTest.png"), width = 2750, height = 4500, res = 400)
par(mfcol = c(1, 1))   
draw(
    ht, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "top"
)
dev.off()

#rename metadata
seu.obj$celltype <- seu.obj$majorID_sub
seu.obj$clusterID <- seu.obj$clusterID_sub

# Export data for UCSC cell browser
ExportToCB_cus(
    seu.obj = seu.obj, dataset.name = "CD8", outDir = "../output/cb_input/", 
    markers = "../output/supplementalData/supplemental_data_4.csv", 
    metaSlots = c(
        "cluster","gene","avg_log2FC","p_val_adj"
    ),
    reduction = "umap",  
    colsTOkeep = c(
        "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase",
        "name", "clusterID", "celltype"
    ),
    skipEXPR = F, test = F,
    feats = c(
        "CD44", "IL2", "TBX21", "SELPLG", "TNFRSF8", "FOXP3", "LAG3", 
        "CCR8", "TNFRSF9", "CD226", "EOMES", "HAVCR2", "CCL5", "CXCR6", 
        "KLRG1", "RUNX3", "IFNG", "IL2", "ICOS", "EOMES"
    )
)

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end CD8 T cell analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


