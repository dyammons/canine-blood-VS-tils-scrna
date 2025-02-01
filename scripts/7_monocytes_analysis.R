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
Idents(seu.obj) <- "conSense"
seu.obj.sub <- subset(
    seu.obj, idents = c("IFN-TAM","M-MDSC" ,"Monocyte, CD4+",
                        "Monocyte, CD4-","Monocyte, IFN signature", "TIM")
)
table(seu.obj.sub$conSense)
table(seu.obj.sub$orig.ident)

seu.sub.list <- SplitObject(seu.obj.sub, split.by = "orig.ident")

seu.obj <- indReClus(
    seu.obj = NULL, outDir = "./output/s2/", subName = "20230507_mono_bloodANDtils", 
    preSub = T, seu.list = seu.sub.list, vars.to.regress = "percent.mt"
)

clusTree(
    seu.obj = seu.obj, dout = "./output/clustree/", 
    outName = "20230507_mono_bloodANDtils", test_dims = c(40,35,30), 
    algorithm = 3, prefix = "integrated_snn_res."
)

seu.obj <- dataVisUMAP(
    seu.obj = seu.obj, outDir = "./output/s3/", 
    outName = "20230507_mono_bloodANDtils", final.dims = 40, final.res = 0.5, 
    stashID = "clusterID_sub", 
    algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.6, 
    n.neighbors = 75, assay = "integrated", saveRDS = T
)

############################################ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin monocyte analysis   ######## <<<<<<<<<<<<<<
############################################ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in processed data
seu.obj <- readRDS("../output/s3/20230507_mono_bloodANDtils_res0.5_dims40_dist0.6_neigh75_S3.rds")
seu.obj$cellSource <- ifelse(grepl("tils",seu.obj$name), "Tumor", "Blood")
outName <- "mono"

#check consenseous
seu.obj$ct <- ifelse(seu.obj$cellSource == "Tumor", 
                     paste0("TIL;", seu.obj$celltype.l3), 
                     paste0("Blood;",seu.obj$celltype.l3_pbmc))
table(seu.obj$ct, seu.obj$clusterID_sub) %>%
    melt() %>%
    separate(Var.1, sep = ";", c("source", "ct")) %>%
    filter(source == "Blood") %>%
    group_by(Var.2) %>%
    mutate(pct = value/sum(value)) %>%
    filter(value == max(value))
table(seu.obj$ct, seu.obj$clusterID_sub) %>%
    melt() %>%
    separate(Var.1, sep = ";", c("source", "ct")) %>%
    filter(source == "TIL") %>%
    group_by(Var.2) %>% 
    mutate(pct = value/sum(value)) %>%
    filter(value == max(value))
Idents(seu.obj) <- "cellSource"
seu.obj.ds <- subset(seu.obj, downsample = min(table(seu.obj$cellSource)))
    table(seu.obj.ds$ct, seu.obj.ds$clusterID_sub) %>%
    sweep(., 2, colSums(.), `/`) %>%
    melt() %>%
    separate(Var.1, sep = ";", c("source", "ct")) %>%
    group_by(source, Var.2) %>%
    filter(value == max(value))

#set metadata
colz.df <- read.csv("./metaData/majorGroups.csv")
colz.df <- colz.df[colz.df$majorID2 == "mono", ]

table(seu.obj$clusterID_sub,seu.obj$conSense)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "CD4-_monocyte (c0)", 
                                   "1" = "CD4+_monocyte (c1)", 
                                   "2" = "CD4-_monocyte (c0)", 
                                   "3" = "CD4-_monocyte (c0)",
                                   "4" = "IFN_monocyte (c2)", 
                                   "5" = "M-MDSC (c3)")
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

#Export annotations
write.csv(seu.obj@meta.data["majorID_sub"], file = paste0("../output/annotations/", outName, ".csv"))

### Supp data - Generate violin plots for each cluster
vilnPlots(
    seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, 
    outName = outName, returnViln = F, 
    outDir = paste0("../output/viln/", outName, "/"), 
    outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
)

### Fig 6a - UMAP by clusterID_sub
pi <- DimPlot(
    seu.obj, 
    reduction = "umap", 
    group.by = "clusID_new",
    cols = colz.df$colour[c(2,5,4,3)],
    pt.size = 0.5,
    label = TRUE,
    label.box = TRUE,
    shuffle = TRUE
)
pi <- cusLabels(
    plot = pi, shape = 21, size = 10, textSize = 6,
    alpha = 0.8, labCol = colz.df$labCol[c(2,5,4,3)]
) + 
    NoLegend() + 
    theme(
        axis.title = element_blank(),
        panel.border = element_blank()
    )
ggsave(paste("../output/", outName, "/", "rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig 6b - skew plot for abundance analysis
p <- skewPlot(
    seu.obj, groupBy = "majorID_sub", yAxisLabel = "Percent of monocytes",
    dout = paste0("../output/", outName), outName = outName, 
    sampleRep = "name", grepTerm = "tils", grepRes = c("Tumor","Blood")
)
ggsave(paste0("../output/", outName, "/", "skewPlot.png"), width = 6, height = 4)


### Fig 6c - DGE analysis
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

### Fig 6d - GO GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", 
                           outName, "_cluster_allCells_all_genes.csv"),
    category = "C5", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = F, pvalueCutoff = 0.05, subcategory = NULL,
    lolli = T,
    saveRes = paste0("../output/", outName, "/c5_", outName, "_res.csv")
) + 
    theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)
    )
p <- p + scale_x_continuous(
    limits = c(-20,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Gene ontology") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_1.png"), width = 7, height = 7)

### Fig 6e - Reactome GSEA of DEGs
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
    limits = c(-20,ceiling(max(p$data$x_axis)*1.05)), 
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
    category = "C7", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
    lolli = T, filterTerm = "MONOCYTE|MACROPHAGE",
    saveRes = paste0("../output/", outName, "/c7_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-50, ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0, ceiling(max(p$data$x_axis)*1.05)/2, ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    labs(title = "TIM vs Blood monocyte",
         subtitle = "(ImmuneSigDB)") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_3.png"), width = 7, height = 7)


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

### Investigate similarities between mono and dc gsea results
dc1 <- read.csv("../output/dc/c2_dc_res.csv")
mono1 <- read.csv("../output/mono/c2_mono_res.csv")
interSECT <- length(intersect(dc1$ID, mono1$ID))
dc1 <- length(dc1$ID) - interSECT
mono1 <- length(mono1$ID) - interSECT
cat(paste0(
    "unique dc: ", dc1, "\n",
    "unique mono: ", mono1, "\n",
    "overlap: ", interSECT, "\n"
))
dc2 <- read.csv("../output/dc/c5_dc_res.csv")
mono2 <- read.csv("../output/mono/c5_mono_res.csv")
interSECT <- length(intersect(dc2$ID, mono2$ID))
dc2 <- length(dc2$ID) - interSECT
mono2 <- length(mono2$ID) - interSECT
cat(paste0(
    "unique dc: ", dc2, "\n",
    "unique mono: ", mono2, "\n",
    "overlap: ", interSECT, "\n"
))

### Fig 6f - Split UMAP of selected DEGs
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$cellSource)))

features <- c("LTF", "IL16","CCR2", "CCL3", "IL1A", "OSM", "CD274", "PTGES", "C1QC")

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



res.df <- read.csv(paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv"))
res.df <- res.df[!grepl("^ENS", res.df$gene), ]
geneList_UP <- res.df %>% filter(padj < 0.1) %>% filter(log2FoldChange > 1) %>% pull(gene)
geneList_DWN <- res.df %>% filter(padj < 0.1) %>% filter(log2FoldChange < -1) %>% pull(gene)

seu.obj$cellSource <- as.factor(seu.obj$cellSource)
p <- splitDot(
    seu.obj = seu.obj, groupBy = "majorID_sub", splitBy = "cellSource", buffer = 150,
    namedColz = setNames(c("#F8766D", "#00BFC4"),  c("Blood", "Tumor")), 
    geneList_UP = unique(c(geneList_UP[1:20], "CCL7", "IL1A",  "PTGES", "OSM", 
                           "CD80", "CD86","C1QC")),
    geneList_DWN = unique(c(geneList_DWN[1:20], "LTF", "IL16","CCR2")), 
    geneColz = c("red", "blue")
)
p <- p +
    theme(
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.justification = 'center'
    )
ggsave(plot = p, paste0("../output/", outName, "/", outName, "_splitDot.png"), width = 13, height = 5)


# p <- splitDot(
#     seu.obj = seu.obj, groupBy = "majorID_sub", splitBy = "cellSource", buffer = 150,
#     namedColz = setNames(c("#F8766D", "#00BFC4"),  c("Blood", "Tumor")), 
#     geneList_UP = c("CXCL10", "CXCL16", "CCL5", "CCL8", "CCL19", "CCL7"),
#     geneList_DWN = NULL, 
#     geneColz = c("red", "blue")
# )
# p <- p +
#     theme(
#         legend.box = "horizontal",
#         legend.direction = "horizontal",
#         legend.position = "bottom",
#         legend.justification = 'center'
#     )
# ggsave(plot = p, paste0("../output/", outName, "/", outName, "_cc_splitDot.png"), width = 7, height = 5)
VlnPlot(
    seu.obj, features = c("CXCL10", "CXCL16", "CCL5", "CCL8", "CCL19", "CCL7"), 
    split.by = "cellSource", group.by = "majorID_sub", stack = TRUE, flip = TRUE
) + theme(axis.title.x = element_blank())
ggsave(paste0("../output/", outName, "/", outName, "_cc_vilnPlot.png"), width = 7, height = 5)


createPB(seu.obj = seu.obj, groupBy = "majorID_sub", comp = "cellSource", biologicalRep = "name",
         outDir = paste0("../output/", outName, "/pseudoBulk/"), min.cell = 10,
         grepTerm = "tils", grepLabel = c("TILs", "Blood"), featsTOexclude = c(pal_feats,tumor.sig), lowFilter = T, dwnSam = F)

pseudoDEG(inDir = paste0("../output/", outName, "/pseudoBulk/"), metaPWD = paste0("../output/", outName, "/pseudoBulk/majorID_sub_deg_metaData.csv"),
                    outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = outName,
                    padj_cutoff = 0.01, lfcCut = 0.58,  idents.1_NAME = "TILs", idents.2_NAME = "Blood", title = "TILS vs Blood", 
                    fromFile = T, returnVolc = F, filterTerm = "^ENSCAF", mkDir = T)


### Supp fig xx -- heatmap of sig DEGs
files <- paste0("../output/", outName, "/pseudoBulk/", levels(seu.obj$majorID_sub), "/", outName, "_cluster_", levels(seu.obj$majorID_sub), "_all_genes.csv")

df.list <- lapply(files, read.csv, header = T)
res.df <- do.call(rbind, df.list)

# write.csv(res.df, file = paste0("../output/supplementalData/", timepoint, "_degs.csv"))

cnts_mat <- res.df  %>% 
    mutate(
        direction = ifelse(log2FoldChange > 0, "Up", "Down")
    ) %>% 
    group_by(gs_base, direction) %>% 
    summarize(nRow = n()) %>% 
    pivot_wider(names_from = gs_base, values_from = nRow) %>% 
    as.matrix() %>% t()

colnames(cnts_mat) <- cnts_mat[1,]
cnts_mat <- cnts_mat[-c(1),]
class(cnts_mat) <- "numeric"

#order by number of total # of DEGs
cnts_mat[is.na(cnts_mat)] <- 0
orderList <- rev(rownames(cnts_mat)[order(rowSums(cnts_mat))])
cnts_mat <- cnts_mat[match(orderList, rownames(cnts_mat)),]        

anno_mat <- cnts_mat
cnts_mat[,1] <- -cnts_mat[,1]

png(file = paste0("../output/", outName, "/", outName, "_deg_heatmap.png"), width=1750, height=2000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(cnts_mat,#name = "mat", #col = col_fun,
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
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),show_heatmap_legend = FALSE)
dev.off()

clus_colz <- colz.df$colour[c(2,5,4,3)]
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
    font_colz = c("black", "white", rep("black", 2)),
    saveName = paste0("../output/", outName, "/", "splitHeat.png"),
    ht_height = 4750, ht_width = 3000
)


### Generate supplemental figure for mono vs macrophage
#load in annotated tumor data
seu.obj.tumor <- readRDS(file = "../output/s3/canine_naive_n6_annotated.rds")

### Select cells to integrate (all pbmcs from OS dogs and only immune cells from tumor)
Idents(seu.obj.tumor) <- "celltype.l3"
table(seu.obj.tumor@active.ident)
seu.obj.tams <- subset(
    seu.obj.tumor,
    idents = c(
        "ANGIO_TAM","LA-TAM_C1QC_hi","LA-TAM_SPP2_hi","TAM_ACT","TAM_INT", 
        "CD4+_TIM", "CD4-_TIM", "IFN-TAM"
    )
)
seu.obj.bl <- subset(seu.obj, subset = cellSource == "Blood")
seu.sub.list <- c(
    SplitObject(seu.obj.bl, split.by = "name"),
    SplitObject(seu.obj.tams, split.by = "name")
)
seu.obj <- indReClus(
    seu.obj = NULL, outDir = "./output/s2/", 
    subName = "20231014_macMono_bloodANDtils", preSub = T, 
    seu.list = seu.sub.list, vars.to.regress = "percent.mt",nfeatures = 2000
)
clusTree(
    seu.obj = seu.obj, dout = "./output/clustree/",
    outName = "20231014_macMono_bloodANDtils", test_dims = c(40,35,30), 
    algorithm = 3, prefix = "integrated_snn_res."
)

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "20231014_macMono_bloodANDtils", final.dims = 35, final.res = 0.3, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.5, n.neighbors = 60, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                      )



### Load in the data
seu.obj <- readRDS("../output/s3/20231014_macMono_bloodANDtils_res0.3_dims35_dist0.5_neigh60_S3.rds")
seu.obj$macMono <- ifelse(!is.na(seu.obj$cellSource), "Monocyte", 
                          ifelse(
                              seu.obj$celltype.l3 == "CD4+_TIM" | 
                              seu.obj$celltype.l3 == "CD4-_TIM" |
                              seu.obj$celltype.l3 == "IFN-TAM",
                              "TIM", "TAM"
                          )
                         )
seu.obj$name2 <- paste0(seu.obj$macMono, "_", seu.obj$name)
#Plot integrated data UMAP
### Fig 6a - UMAP by clusterID_sub
pi <- DimPlot(
    seu.obj, 
    reduction = "umap", 
    group.by = "macMono",
    split.by = "macMono",
    pt.size = 0.5#,
#     label = TRUE,
#     label.box = TRUE,
#     shuffle = TRUE
)
ggsave(paste("../output/", outName, "/", "rawUMAP.png", sep = ""), width = 10, height = 4)

seu.obj.sub <- subset(seu.obj, invert = T, subset = macMono == "TIM")
seu.obj.sub$macMono <- as.factor(seu.obj.sub$macMono)
seu.obj.sub$allCells <- "allCells"
seu.obj.sub$allCells <- as.factor(seu.obj.sub$allCells)
createPB(
    seu.obj = seu.obj.sub, groupBy = "allCells", 
    comp = "cellSource", biologicalRep = "name",
    outDir = paste0("../output/", outName, "/pseudoBulk/"), 
    grepTerm = "N", grepLabel = c("TAM", "Monocyte"), 
    featsTOexclude = c(pal_feats,tumor.sig), 
    lowFilter = T, dwnSam = F
)
p_volc <- pseudoDEG(
    inDir = paste0("../output/", outName, "/pseudoBulk/"), 
    metaPWD = paste0("../output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"),
    outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = outName,
    padj_cutoff = 0.01, lfcCut = 0.58, strict_lfc = T,
    idents.1_NAME = "TAM", idents.2_NAME = "Monocyte", title = "TILS vs Blood", 
    fromFile = T, returnVolc = T, filterTerm = "^ENSCAF", mkDir = T
)
p <- prettyVolc(
    plot = p_volc[[1]], rightLab = "TAM", leftLab = "Blood monocyte", 
    rightCol = "red", leftCol = "blue", arrowz = T
) + 
    theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_blank(),
        legend.position = c(0.10, 0.85)
    )
ggsave(paste0("../output/", outName, "/", "volcPlot.png"), width = 7, height = 7)

### Fig 6e - Reactome GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/",
                           "TAM_vs_Monocyte_all_genes.csv"),
    geneList = NULL, category = "C2", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "CP:REACTOME", 
    lolli = T,
    saveRes = paste0("../output/", outName, "/c2_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-40,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Reactome") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_TAM_2.png"), width = 7, height = 7)


### Fig 2e - ImmuneSigDB GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/",
                           "TAM_vs_Monocyte_all_genes.csv"),
    geneList = NULL, category = "C7", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
    lolli = T, filterTerm = "MONOCYTE|MACROPHAGE",
    saveRes = paste0("../output/", outName, "/c7_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-50,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    labs(title = "TAM vs Blood monocyte",
         subtitle = "(ImmuneSigDB)") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_TAM_3.png"), width = 7, height = 7)



#complete hc
seu.obj$type <- paste0(seu.obj$macMono, "_", seu.obj$name)

metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@data)) #use log noralized count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL                                                                       

M <- (1- cor(t(clusAvg_expression),method="pearson"))/2
hc <- hclust(as.dist(M),method="complete")

ggtree(as.phylo(hc)) + geom_tiplab(offset = 0.003) + xlim(NA,0.1) #+ geom_tippoint(shape = 21,size = 8,alpha = 1, colour="black", fill = "white") + geom_tiplab(aes(label = seq(1:43)-1,colour=c("black"),offset = -0.001))
ggsave(paste("../output/", outName, "/_hc.png", sep = ""), width = 8, height = 8)

# Transpose count matrix and calculate distances using dist()
sampleDists <- clusAvg_expression %>% dist()

# Convert distance dataframe to matrix
sampleDistMatrix <- as.matrix(sampleDists)

# Choose continuous palette from RColorBrewer
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Plot sample distance heatmap with ComplexHeatmap
png(file = paste0("../output/", outName, "/samDist.png"), width=4000, height=4000, res=400)
par(mfcol=c(1,1))         

ht <- Heatmap(sampleDistMatrix, #name = "mat", #col = col_fun,
              name = "Sample distance",
              cluster_rows = T,#hc,
              row_title = "",
              row_title_gp = gpar(fontsize = 24),
              col=colors,#viridis(option = "magma",100),
              cluster_columns = T,#hc,
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


### Core cross species TIMs signature
#Compare to human BrCA TEMo vs Mono
hu_degs <- read.csv("./metaData/brCA_TEMo_VS_Mo.csv") %>%
    column_to_rownames(var = "gene.name")

hu_degs <- read.csv("./metaData/lung_TIM_vs_Mono_Zilionis.csv") %>%
    column_to_rownames(var = "gene")

genez <- orthogene::convert_orthologs(
    gene_df = rownames(hu_degs), gene_output = "columns", 
    input_species = "human", output_species = "dog", 
    non121_strategy = "drop_both_species"
)
genez <- genez %>% 
    left_join(rownames_to_column(hu_degs), by = c("input_gene" = "rowname"))

degs_dog <- read.csv(paste0(
    "../output/", outName, "/pseudoBulk/allCells/", outName, 
    "_cluster_allCells_all_genes.csv"
))
degs_dog <- degs_dog[degs_dog$log2FoldChange > 0, ]
genez_up <- genez[genez$FDR < 0.01 & genez$fold_change > 1.5, ]$ortholog_gene
overlap_tims <- intersect(genez_up, degs_dog$gene)

hu_only <- genez_up[! genez_up %in% overlap_tims]
can_only <- degs_dog$gene
can_only <- genez_up[! can_only %in% overlap_tims]


gene_res <- list(
    "Human" = genez_up,
    "Canine" = degs_dog$gene
)

#Make the plot
VennDiagram::venn.diagram(
    x = gene_res,
    filename = paste0("../output/", outName, "/venn.png"),
    output = TRUE,
    imagetype = "png",
    fontfamily = "sans",
    cat.fontfamily = "sans",
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.dist = c(0.05, 0.05),
    margin = 0.05,
    disable.logging = TRUE,
    height = 2250, 
    width = 2250
)

### Fig 2e - ImmuneSigDB GSEA of DEGs
p <- plotGSEA(
    geneList = hu_only, category = "C7", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
    lolli = T, filterTerm = "MONOCYTE|MACROPHAGE",
    saveRes = paste0("../output/", outName, "/c7_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-160, ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0, ceiling(max(p$data$x_axis)*1.05)/2, ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    labs(title = "Human TIMs vs Blood monocyte",
         subtitle = "(ImmuneSigDB)") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_huCan_1.png"), width = 7, height = 7)

### Fig 2e - ImmuneSigDB GSEA of DEGs
p <- plotGSEA(
    geneList = can_only, category = "C7", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
    lolli = T, filterTerm = "MONOCYTE|MACROPHAGE",
    saveRes = paste0("../output/", outName, "/c7_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-160, ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0, ceiling(max(p$data$x_axis)*1.05)/2, ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    labs(title = "Canine TIMs vs Blood monocyte",
         subtitle = "(ImmuneSigDB)") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_huCan_1.png"), width = 7, height = 7)



### Fig 2e - ImmuneSigDB GSEA of DEGs
p <- plotGSEA(
    geneList = genez_up, category = "C7", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
    lolli = T, filterTerm = "MONOCYTE|MACROPHAGE",
    saveRes = paste0("../output/", outName, "/c7_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-160, ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0, ceiling(max(p$data$x_axis)*1.05)/2, ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    labs(title = "Human NSCLC TIMs vs Blood monocyte",
         subtitle = "(ImmuneSigDB)") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_huCan_1.png"), width = 7, height = 7)



### Fig 2e - ImmuneSigDB GSEA of DEGs
p <- plotGSEA(
    geneList = overlap_tims, category = "C7", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
    lolli = T, filterTerm = "MONOCYTE|MACROPHAGE",
    saveRes = paste0("../output/", outName, "/c7_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-60, ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0, ceiling(max(p$data$x_axis)*1.05)/2, ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    labs(title = "Core TIM signature",
         subtitle = "(ImmuneSigDB)") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_huCan_2.png"), width = 7, height = 7)

#rename metadata
seu.obj$celltype <- seu.obj$majorID_sub
seu.obj$clusterID <- seu.obj$clusID_new

# Export data for UCSC cell browser
ExportToCB_cus(
    seu.obj = seu.obj, dataset.name = "Mono", outDir = "../output/cb_input/", 
    markers = "../output/supplementalData/supplemental_data_7.csv", 
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
        "LTF", "IL16","CCR2", "CCL3", "IL1A", "OSM", "CD274", "PTGES", "C1QC"
    )
)

########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end monocyte analysis   ######## <<<<<<<<<<<<<<
########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

library(CellChat)

##############################
### CELLCHAT PREPROCESSING ###
##############################


seu.obj$grouping <- seu.obj$macMono

#get 1:1 orthologues
cnts <- seu.obj@assays$RNA@data
cnts <- orthogene::convert_orthologs(
    gene_df = cnts,
    gene_input = "rownames", 
    gene_output = "rownames", 
    input_species = "ecaballus",
    output_species = "human",
    non121_strategy = "drop_both_species"
)
rownames(cnts) <- unname(rownames(cnts))
#prep metadata and create cellchat object
meta <- seu.obj@meta.data
cell.use <- rownames(meta)
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "grouping")
cellchat@idents <- factor(
    cellchat@idents, 
    levels = as.character(str_sort(levels(cellchat@idents), numeric = TRUE))
)
cellchat@DB <- CellChatDB.human
#run standard cellchat v1 workflow
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, "../output/cellchat/monocyte_cellChatobj.rds")


##################################
### END CELLCHAT PREPROCESSING ###
##################################


###############################
### BEGIN CELLCHAT ANALYSIS ###
###############################

#set output specifications
outName <- "cellchat"

#load in processed cellchat data
cellchat <- readRDS("../output/cellchat/monocyte_cellChatobj.rds")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
pathways <- cellchat@netP$pathways
pathwayz <- pathways
subName <- "cellCom"
lapply(pathwayz, function(pathway){
    
    #extract required plotting data
    lrData <- as.data.frame(cellchat@LR)
    net <- cellchat@net
    prob <- net$prob
    prob <- prob[,,rownames(lrData[lrData$LRsig.pathway_name == pathway,])]
    prob.sum <- apply(prob, c(1,2), sum)

    #identify which cell types are active in pathway
    idx1 <- which(Matrix::rowSums(prob.sum) == 0)
    idx2 <- which(Matrix::colSums(prob.sum) == 0)
    idx <- intersect(idx1, idx2)
    net <- prob.sum[-idx, ]
    net <- net[, -idx]
    cellTypeOFinterest <- rownames(net)

#     grey out cell types not involved
#     colz2 <- colz
#     colz2[names(colz2)[!names(colz2) %in% cellTypeOFinterest]] <- "grey"

    #save the plot
    outfile <- paste("../output/", outName, "/", subName, "/", pathway ,"_cell_cell_3.png", sep = "")
    png(file = outfile, width=2500, height=2500, res=400)
    par(mfrow=c(1,1), mar = c(0,0,0,0))
    gg7 <- netVisual_aggregate(
        cellchat, layout = "chord", signaling = pathway, #group = groupzNames, 
        remove.isolate = F, big.gap = 5
    ) # color.use = colz2, 
    dev.off()

    #get the active features in the pathway
    genez <- lapply(pathways, function(x){extractEnrichedLR(cellchat, signaling = x, geneLR.return = TRUE, enriched.only = T)[["geneLR"]]})
    names(genez) <- pathways

    Idents(seu.obj) <- "grouping"
    seu.obj.sub <- subset(seu.obj, idents = cellTypeOFinterest)
    
    #plot expression using Seurat function
    pi <- VlnPlot(
        object = seu.obj,
        pt.size = 0,
        same.y.lims = T,
        flip = T,
        group.by = "grouping",
        fill.by = "ident",
#         cols = colz2,
        stack = TRUE,
        combine = FALSE,
        features = unlist(genez[pathway])
    ) + NoLegend() + theme(axis.title.x = element_blank(),
                           axis.title.y.right = element_blank())
    
    ggsave(paste("../output/", outName, "/", subName, "/", pathway ,"_viln.png", sep = ""), height = 7, width = 7)
    
})
