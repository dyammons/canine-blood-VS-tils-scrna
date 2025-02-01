#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")

### Analysis note: 

################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin neutrophil preprocessing   ######## <<<<<<<<<<<<<<
################################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set output params
outName <- "neut"

#subset on the monocytes & complete subset analysis
Idents(seu.obj) <- "conSense"
seu.obj.sub <- subset(seu.obj, 
                          ident = c("Neutrophil", "PMN-MDSC")
                         )

table(seu.obj.sub$conSense)
min(table(seu.obj.sub$name)) > 100
seu.sub.list <- SplitObject(seu.obj.sub, split.by = "name")
seu.sub.list <- seu.sub.list[-c(3,11)] #remove samples with too few cells - innacurate clsutering

seu.obj <- indReClus(seu.obj = NULL, outDir = "./output/s2/", subName = "20230507_neuts_bloodANDtils", preSub = T, seu.list = seu.sub.list,
                      vars.to.regress = "percent.mt",nfeatures = 2000
                       )

clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "20230507_neuts_bloodANDtils", test_dims = c(40,35,30), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "20230507_neuts_bloodANDtils", final.dims = 40, final.res = 0.2, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.5, n.neighbors = 60, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                      )

############################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin neutrophil analysis   ######## <<<<<<<<<<<<<<
############################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in processed data
seu.obj <- readRDS("../output/s3/20230507_neuts_bloodANDtils_res0.2_dims40_dist0.5_neigh60_S3.rds")
seu.obj$cellSource <- ifelse(grepl("tils",seu.obj$name),"Tumor","Blood")
outName <- "neut"

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
colz.df <- colz.df[colz.df$majorID2 == "neut", ]

table(seu.obj$clusterID_sub,seu.obj$conSense)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "#F8766D", "1" = "#00BA38", 
                                   "2" = "#619CFF")
                       )
seu.obj$sub_colz <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Neutrophil (c0)", "1" = "PMN_MDSC (c1)", 
                                   "2" = "Neutrophil (c0)")
                       )
seu.obj$majorID_sub <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "0", "1" = "1", 
                                   "2" = "0")
                       )
seu.obj$clusID_new <- Idents(seu.obj)

#Export annotations
write.csv(seu.obj@meta.data["majorID_sub"], file = paste0("../output/annotations/", outName, ".csv"))

### Supp data - Generate violin plots for each cluster
vilnPlots(
    seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, 
    outName = outName, returnViln = F, outDir = paste0("../output/viln/", outName, "/"), 
    outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
)


pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "cellSource",
              split.by = "name",
              pt.size = 0.5,
              label = F,
              label.box = F,
              shuffle = F,
              ncol = 4
)
pi <- formatUMAP(plot = pi)
ggsave(paste("../output/", outName, "/", "supp_label.png", sep = ""), width = 16, height = 16)


### Fig 7a - UMAP by clusID_new
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusID_new",
              cols = c(colz.df$colour[1],"hotpink"),
              pt.size = 0.5,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, 
                alpha = 0.8, labCol = "black") + NoLegend() + theme(axis.title = element_blank(),
                                                                                                               panel.border = element_blank())
ggsave(paste("../output/", outName, "/", "rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig 7b - skew plot for abundance analysis
p <- skewPlot(
    seu.obj, groupBy = "majorID_sub", yAxisLabel = "Percent of neutrophils",
    dout = paste0("../output/", outName), outName = outName, 
    sampleRep = "name", grepTerm = "tils", grepRes = c("Tumor","Blood")
)
ggsave(paste0("../output/", outName, "/", "skewPlot.png"), width = 6, height = 4)

### Fig 7c - DGE analysis
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

### Fig 7d - GO GSEA of DEGs
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

### Fig 7e - Reactome GSEA of DEGs
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
    limits = c(-10,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Reactome") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_2.png"), width = 7, height = 7)

### Fig 7f - Split UMAP of selected DEGs
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$cellSource)))

features <- c("PLAUR", "PLAU", "IL18BP", "CXCL8", "OSM", "IL1R2", "IL22RA2", "IL1A", "CD274")

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

clus_colz <- c(colz.df$colour[1],"hotpink")
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
    font_colz = rep("black", 2),
    saveName = paste0("../output/", outName, "/", "splitHeat.png"),
    ht_height = 4000, ht_width = 2500
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

#obsolete
# res.df <- read.csv(paste0("../output/", outName, "/pseudoBulk/allCells/", outName, "_cluster_allCells_all_genes.csv"))
# res.df <- res.df[!grepl("^ENS", res.df$gene), ]
# geneList_UP <- res.df %>% filter(padj < 0.1) %>% filter(log2FoldChange > 1) %>% pull(gene)
# geneList_DWN <- res.df %>% filter(padj < 0.1) %>% filter(log2FoldChange < -1) %>% pull(gene)

# seu.obj$cellSource <- as.factor(seu.obj$cellSource)
# p <- splitDot(
#     seu.obj = seu.obj, groupBy = "majorID_sub", splitBy = "cellSource", buffer = 125,
#     namedColz = setNames(c("#F8766D", "#00BFC4"),  c("Blood", "Tumor")), 
#     geneList_UP = unique(c(geneList_UP[1:20], "PLAUR", "PLAU", "IL18BP", 
#                     "CXCL8", "OSM", "IL1R2", "IL22RA2", "IL1A", "CD274")),
#     geneList_DWN = geneList_DWN[1:4], geneColz = c("red", "blue")
# )
# p <- p +
#     theme(
#         legend.box = "horizontal",
#         legend.direction = "horizontal",
#         legend.position = "bottom",
#         legend.justification = 'center'
#     )
# ggsave(plot = p, paste0("../output/", outName, "/", outName, "_splitDot.png"), width = 10, height = 3.5)

# VlnPlot(seu.obj, features = c("CXCL8", "CCL5", "CCL7"), 
#         split.by = "cellSource", group.by = "majorID_sub") & theme(axis.title.x = element_blank())
# ggsave(paste0("../output/", outName, "/", outName, "_cc_vilnPlot.png"), width = 9, height = 4)

# #load in the tumor and pal signatures to exlude from DE analysis
# pal_feats = c('TIMP1', 'NAA10', 'ENSCAFG00000037735', 'GP6', 'SEC11C', 'FTL', 'NRGN', 'ACOT7', 'VCL', 'RSU1', 'ITGB1', 'H3-3A', 'RABGAP1L', 'SELP', 'SH3GLB1', 'ACTB', 'ENSCAFG00000008221', 'TLN1', 'GSN', 'AMD1', 'TREM2', 'SH3BGRL2', 'MYH9', 'PLEK', 'ENSCAFG00000042554', 'RAP1B', 'ENSCAFG00000004260', 'NAP1L1', 'PPBP', 'RASA3', 'ITGA2B', 'EIF1', 'ACTG1', 'C9H17orf64', 'JMJD6', 'CCL14', 'GNG11', 'IGF2BP3', 'TBXAS1', 'VDAC3', 'MARCHF2', 'TPM4', 'TKT', 'FTH1.1', 'FERMT3', 'RTN3', 'PRKAR2B', 'SVIP', 'ENSCAFG00000030286', 'ADA', 'MYL9', 'TUBB1', 'TUBA1B', 'METTL7A', 'THBS1', 'SERF2', 'PIF1', 'B2M', 'GAS2L1', 'YWHAH', 'HPSE', 'ATG3', 'ENSCAFG00000015217', 'ITGA6','RGS18', 'SUB1', 'LGALS1', 'CFL1', 'BIN2', 'CAT', 'RGS10', 'MGST3', 'TMBIM6', 'PFN1', 'CD63', 'RALBP1', 'GNAS', 'SEPTIN7', 'TPT1', 'UBB', 'ATF4', 'BBLN', 'MTDH', 'ENSCAFG00000017655','FYB1', 'ENO1', 'GABARAP', 'SSR4', 'MSN', 'ENSCAFG00000011134', 'ENSCAFG00000046637', 'COX8A', 'DLA-64', 'CD47', 'VASP', 'DYNLRB1', 'DLA88', 'SMDT1', 'ATP5PF','ELOB', 'ENSCAFG00000029155', 'ARPC3', 'VPS28', 'LRRFIP1', 'SRP14', 'ABRACL', 'ENSCAFG00000043577', 'ENSCAFG00000042598')
# tumor.sig <- read.csv("./metaData/tumorSig.csv", header = T)$x

# createPB(seu.obj = seu.obj, groupBy = "majorID_sub", comp = "cellSource", biologicalRep = "name",
#          outDir = paste0("../output/", outName, "/pseudoBulk/"), 
#          grepTerm = "tils", grepLabel = c("TILs", "Blood"), featsTOexclude = c(pal_feats,tumor.sig), lowFilter = T, dwnSam = F)

# pseudoDEG(inDir = paste0("../output/", outName, "/pseudoBulk/"), metaPWD = paste0("../output/", outName, "/pseudoBulk/majorID_sub_deg_metaData.csv"),
#                     outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = outName,
#                     padj_cutoff = 0.01, lfcCut = 0.58,  idents.1_NAME = "TILs", idents.2_NAME = "Blood", title = "TILS vs Blood", 
#                     fromFile = T, returnVolc = F, filterTerm = "^ENSCAF", mkDir = T)


# ### Supp fig xx -- heatmap of sig DEGs
# files <- paste0("../output/", outName, "/pseudoBulk/", levels(seu.obj$majorID_sub), "/", outName, "_cluster_", levels(seu.obj$majorID_sub), "_all_genes.csv")

# df.list <- lapply(files, read.csv, header = T)
# res.df <- do.call(rbind, df.list)

# # write.csv(res.df, file = paste0("../output/supplementalData/", timepoint, "_degs.csv"))

# cnts_mat <- res.df  %>% 
#     mutate(
#         direction = ifelse(log2FoldChange > 0, "Up", "Down")
#     ) %>% 
#     group_by(gs_base, direction) %>% 
#     summarize(nRow = n()) %>% 
#     pivot_wider(names_from = gs_base, values_from = nRow) %>% 
#     as.matrix() %>% t()

# colnames(cnts_mat) <- cnts_mat[1,]
# cnts_mat <- cnts_mat[-c(1),]
# class(cnts_mat) <- "numeric"

# #order by number of total # of DEGs
# cnts_mat[is.na(cnts_mat)] <- 0
# orderList <- rev(rownames(cnts_mat)[order(rowSums(cnts_mat))])
# cnts_mat <- cnts_mat[match(orderList, rownames(cnts_mat)),]        

# anno_mat <- cnts_mat
# cnts_mat[,1] <- -cnts_mat[,1]

# png(file = paste0("../output/", outName, "/", outName, "_deg_heatmap.png"), width=1500, height=1500, res=400)
# par(mfcol=c(1,1))         
# ht <- Heatmap(cnts_mat,#name = "mat", #col = col_fun,
#               name = "# of DEGs",
#               cluster_rows = F,
#               row_title = "Cell type",
#               col = circlize::colorRamp2(c(min(cnts_mat), 0,max(cnts_mat)), colors = c("blue", "white", "red")),
#               cluster_columns = F,
#               show_column_names = TRUE,
#               column_title_side = "top",
#               column_names_rot = 0,
#               column_names_gp = gpar(fontsize = 14, col = "black"),
#               column_names_centered = TRUE,
#               heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
#                                           labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
#               cell_fun = function(j, i, x, y, width, height, fill) {
#                   if(cnts_mat[i, j] < -10) {
#                       grid.text(sprintf("%.0f", as.matrix(anno_mat)[i, j]), x, y, gp = gpar(fontsize = 14, col = "white"))
#                   } else if(cnts_mat[i, j] > -10) {
#                       grid.text(sprintf("%.0f", as.matrix(anno_mat)[i, j]), x, y, gp = gpar(fontsize = 14, col = "black"))
#                   }
#               })
# draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),show_heatmap_legend = FALSE)
# dev.off()

# clus_colz <- c(colz.df$colour[1],"hotpink")
# names(clus_colz) <- levels(seu.obj$majorID_sub)
# cond_colz <- gg_color_hue(2)
# names(cond_colz) <- c("Blood","Tumor")

# genez <- res.df %>% 
#     filter(!grepl("^ENS", gene)) %>%
#     mutate(
#         grp = paste0(gs_base, ifelse(log2FoldChange > 0, "_UP", "_DNN"))
#     ) %>%
#     group_by(grp) %>%
#     top_n(-15, padj) %>%
#     pull(gene)
# res.df <- res.df[res.df$gene %in% c(genez, features), ]

# res.df$gs_base <- toupper(res.df$gs_base)
# ht <- sigDEG_heatmap(
#     seu.obj = seu.obj, groupBy = "majorID_sub", splitBy = "cellSource", forceCleanPlot = T, 
#     dge_res = res.df, lfc_thres = 1, cond_colz = cond_colz, clus_colz = clus_colz,
#     font_colz = rep("black", 2),
#     saveName = paste0("../output/", outName, "/", "splitHeat.png"),
#     ht_height = 4000, ht_width = 2500
# )

### Core cross species TIMs signature
#Compare to human BrCA TEMo vs Mono
hu_degs <- read.csv("./metaData/lung_TAN_vs_neut_Zilionis.csv") %>%
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
genez_up <- genez[genez$FDR < 0.01 & genez$fold_change > 1, ]$ortholog_gene
overlap <- intersect(genez_up, degs_dog$gene)

hu_only <- genez_up[! genez_up %in% overlap]
can_only <- degs_dog$gene
can_only <- genez_up[! can_only %in% overlap]


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
    lolli = T, filterTerm = "NEUTROPHIL",
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
    labs(title = "Human TAN vs Blood monocyte",
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
    lolli = T, filterTerm = "NEUTROPHIL",
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
    labs(title = "Canine TAN vs Blood monocyte",
         subtitle = "(ImmuneSigDB)") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_huCan_2.png"), width = 7, height = 7)



### Fig 2e - ImmuneSigDB GSEA of DEGs
p <- plotGSEA(
    geneList = genez_up, category = "C7", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
    lolli = T, filterTerm = "NEUTROPHIL",
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
    labs(title = "Human NSCLC TANs vs Blood monocyte",
         subtitle = "(ImmuneSigDB)") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_huCan_3.png"), width = 7, height = 7)



### Fig 2e - ImmuneSigDB GSEA of DEGs
p <- plotGSEA(
    geneList = overlap, category = "C7", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
    lolli = T, filterTerm = "NEUTROPHIL",
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
    labs(title = "Core TAN signature",
         subtitle = "(ImmuneSigDB)") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_huCan_4.png"), width = 7, height = 7)

#rename metadata
seu.obj$celltype <- seu.obj$majorID_sub
seu.obj$clusterID <- seu.obj$clusID_new

# Export data for UCSC cell browser
ExportToCB_cus(
    seu.obj = seu.obj, dataset.name = "Neut", outDir = "../output/cb_input/", 
    markers = "../output/supplementalData/supplemental_data_8.csv", 
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
        "PLAUR", "PLAU", "IL18BP", "CXCL8", "OSM", 
        "IL1R2", "IL22RA2", "IL1A", "CD274"
    )
)

########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end neut analysis   ######## <<<<<<<<<<<<<<
########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#Compile DE results in major pops and subclusters (supplemental_data_1.csv)
subsets <- c("cd4", "cd8", "bcell", "plasma","dc", "mono", "neut")
df.list <- lapply(seq_along(subsets), function(x){
    files <- list.dirs(
        path = paste0("../output/", subsets[x], "/pseudoBulk"), recursive = FALSE
    )
    dirNames <- unlist(lapply(files, function(x) {unlist(strsplit(x, split = "/"))[5]}))
    dirNames <- dirNames[! grepl("^\\.", dirNames)]
    df.list <- lapply(seq_along(dirNames), function(y){
        basePath <- paste0("../output/", subsets[x], "/pseudoBulk/", dirNames[y], "/")
        df <- read.csv(file = list.files(path = basePath, 
                                         pattern = "_all_genes.csv",
                                         full.names = T))
        df$gs_base <- "Tumor_VS_Blood"
        if(dirNames[y] == "allCells") {
            cOntrast <- paste0("All cells (", toupper(subsets[x]), ")")
        } else {
            cOntrast <- dirNames[y]
        }
        df$subset <- cOntrast
        return(df)
    })
    df <- do.call(rbind, df.list)
    return(df)
})
res <- do.call(rbind, df.list)
#append with tam vs monocyte
dat2 <- read.csv("../output/mono/TAM_vs_Monocyte_all_genes.csv")
dat2$subset <- "TAM_and_monocyte"
res <- rbind(res, dat2)
write.csv(res, file = "../output/supplementalData/supplemental_data_2.csv", row.names = F)

#Compile GSEA results (supplemental_data_2.csv)
read_gsea <- function(
    inFile = NULL,
    msigdb_reference = NULL,
    major_pop = NULL
){
    if(grepl(".csv", inFile) & file.exists(inFile)){
        df <- read.csv(file = inFile)
        if(nrow(df) > 0) {
            df$msigdb_reference <- msigdb_reference
            df$subset <- toupper(major_pop)
            print(inFile)
            df$celltype <- toupper(substr(
                inFile,
                nchar(paste0("../output/", major_pop, major_pop)) + 6,
                nchar(inFile) - 8
            ))
            if(nchar(df$celltype[1]) < 1){
                df$celltype <- paste0(toupper(major_pop), " (All cells)")
            }
            return(df) 
        } else {
            return(NULL)
        }
    } else {
        message("No file found for ", inFile)
    }

}
    
subsets <- c("cd4", "cd8", "bcell", "plasma", "dc", "mono", "neut")
df.list <- lapply(subsets, function(x){
    inFiles <- list.files(paste0("../output/", x, "/"), pattern = "_res.csv")
    c5Files <- c(
        paste0("../output/", x, "/", inFiles[grepl("C5.*_res.csv", inFiles)]),
        paste0("../output/", x, "/", "c5_", x, "_res.csv")
    )
    df.list_c5 <- lapply(c5Files, read_gsea, msigdb_reference = "C5", major_pop = x)
    df_c5 <- do.call(rbind, df.list_c5)

    c2Files <- c(
        paste0("../output/", x, "/", inFiles[grepl("C2.*_res.csv", inFiles)]),
        paste0("../output/", x, "/", "c2_", x, "_res.csv")
    )
    df.list_c2 <- lapply(c2Files, read_gsea, msigdb_reference = "C2", major_pop = x)
    df_c2 <- do.call(rbind, df.list_c2)
    
    c7Files <- c(
#         paste0("../output/", x, "/", inFiles[grepl("C7.*_res.csv", inFiles)]),
        paste0("../output/", x, "/", "c7_", x, "_res.csv")
    )
    df.list_c7 <- lapply(c7Files, read_gsea, msigdb_reference = "C7", major_pop = x)
    df_c7 <- do.call(rbind, df.list_c7)
    
    return(rbind(df_c2, df_c5, df_c7))
})
df_save <- do.call(rbind, df.list) %>%
    select(-X)
write.csv(df_save, file = "../output/supplementalData/supplemental_data_3.csv", row.names = F)


# Prep supp data
subsets <- c("cd4", "cd8", "bcell", "dc", "mono", "neut")
lapply(seq_along(subsets), function(x){
    df <- read.csv(file = paste0("../output/viln/", subsets[x], "/", subsets[x], "_gene_list.csv"), row.names = 1)
    df <- df[ , c(7, 6, 1, 5, 2:4)]
    write.csv(df, paste0("../output/supplementalData/supplemental_data_", x + 2, ".csv"), row.names = F)
})


#Compile supp tables, need to manually add tables 1-3
library(openxlsx)
subsets <- c("allCells", "cd4", "cd8", "bcell", "dc", "mono", "neut")
df.list <- lapply(seq_along(subsets), function(x){
    df <- read.csv(file = paste0("../output/", subsets[x], "/", subsets[x], "_skewPlot_stats.csv"))
})

names(df.list) <- paste0("supp_table_", seq_along(subsets) + 3)
write.xlsx(df.list, file = "../output/supplementalData/supplemental_tables.xlsx")