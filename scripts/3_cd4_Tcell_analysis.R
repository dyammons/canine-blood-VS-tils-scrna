#!/usr/bin/Rscript

### Author: Dylan Ammons
### Analysis code for [title here]

### Purpose: 
# Load in the clean, integrated object that contains canine pbmcs and tils to
# then subset on CD4 T cells and complete subcluster analysis. With the
# clustered data, complete differential gene expression and abundance analysis

### Analysis approach (2 steps):

## (1) Subcluster analysis - Load in cleaned object generated from 
# `2_allCells_analysis.R`, subset on CD4 T cells, then re-integrate data and
# repeat unsupervised clustering.

## (2) Plot figures (Figure 1) - Complete visualization and run 
# differential abundance analysis

############################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin CD4 T cell analysis   ######## <<<<<<<<<<<<<<<
############################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


### (1) Subcluster analysis

## Load custom functions, packages, and set output path
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")
outName <- "cd4"

## Load in data and subset on cell types of interest
seu.obj <- readRDS("../output/s3/bloodANDtils_filtered_allCells_S3.rds")
Idents(seu.obj) <- "conSense"
seu.obj.sub <- subset(
    seu.obj, 
    ident = c(
        "CD4+ Naive", "CD4+ T reg", "CD4+ TCM", "CD4+ TEM", 
        "CD4+ TEM, Th1-like","CD4+ TEM, Th17-like", "CD4+ TEM, Th2-like",
        "CD4_act","CD4_ex","CD4_naive","CD4_reg"
    )
)
#Check that all samples have at least 100 cells
min(table(seu.obj.sub$orig.ident)) > 100
#Split data by sample then reintegrate and cluster
seu.sub.list <- SplitObject(seu.obj.sub, split.by = "orig.ident")
seu.obj <- indReClus(
    seu.obj = NULL, outDir = "../output/s2/", 
    subName = "20230507_cd4_bloodANDtils", preSub = T, seu.list = seu.sub.list,
    vars.to.regress = "percent.mt", nfeatures = 2000
)
clusTree(
    seu.obj = seu.obj, dout = "../output/clustree/", 
    outName = "20230507_cd4_bloodANDtils", 
    test_dims = 40, algorithm = 3, prefix = "integrated_snn_res."
)
seu.obj <- dataVisUMAP(
    seu.obj = seu.obj, 
    outDir = "../output/s3/", outName = "20240305_cd4_bloodANDtils", 
    final.dims = 40, min.dist = 0.6, n.neighbors = 75,
    final.res = 0.7, algorithm = 3,
    prefix = "integrated_snn_res.",  stashID = "clusterID_sub", 
    assay = "integrated", saveRDS = T,
    features = c(
        "PTPRC", "CD3E", "CD8A", "GZMA", 
        "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
        "CD4", "MS4A1", "PPBP","HBM"
    )
)


### (2) Data visualization
#Load in processed data (in case not already loaded)
seu.obj <- readRDS("../output/s3/20230507_cd4_bloodANDtils_res0.7_dims40_dist0.6_neigh75_S3.rds")
outName <- "cd4"

#Correct metadata errors
seu.obj$cellSource <- ifelse(grepl("tils", seu.obj$name), "Tumor", "Blood")
seu.obj$celltype.l3 <- as.factor(seu.obj$celltype.l3)
levels(seu.obj$celltype.l3)[2] <- "CD4_fh"

#Manually join clusters based on clustree output to better fit annotations
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(
    seu.obj, 
    c(
        "0" = "1", "1" = "2", 
        "2" = "0", "3" = "0",
        "4" = "3", "5" = "4",
        "6" = "5", "7" = "1"
    )
)
seu.obj$clusID_new <- Idents(seu.obj)
seu.obj$clusID_new <- factor(seu.obj$clusID_new, levels = c(0, 1, 2, 3, 4, 5))
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
    melt() %>% 
    separate(Var.1, sep = ";", c("source", "ct")) %>% 
    group_by(source, Var.2) %>% 
    filter(value == max(value))
#Update the metadata with cell type annotations
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(
    seu.obj, 
    c(
        "0" = "TCM (c1)", "1" = "Treg/Tfh (c2)", 
        "2" = "Naive (c0)", "3" = "Naive (c0)",
        "4" = "TEM_Th2-like (c3)", "5" = "TEM (c4)",
        "6" = "TEM_Th1-like (c5)", "7" = "TCM (c1)"
    )
)
seu.obj$majorID_sub <- Idents(seu.obj)
seu.obj$majorID_sub <- factor(
    seu.obj$majorID_sub, 
    levels = c(
        "TCM (c1)", "Treg/Tfh (c2)", 
        "Naive (c0)", "TEM_Th2-like (c3)", 
        "TEM (c4)", "TEM_Th1-like (c5)"
    )[c(3,1,5,6,4,2)]
)
Idents(seu.obj) <- "clusID_new"
seu.obj <- RenameIdents(
    seu.obj, 
    c(
        "0" = "#0E3F5C", "1" = "#009DA5", 
        "2" = "#148B7D", "3" = "#D4F3A3",
        "4" = "#9CD6BA", "5" = "#2E4F79"
    )
)
seu.obj$sub_colz <- Idents(seu.obj)

#Export annotations
write.csv(seu.obj@meta.data["majorID_sub"], file = paste0("../output/annotations/", outName, ".csv"))

## Supplemental data - Generate gene lists for each cluster
vilnPlots(
    seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, 
    outName = outName, outDir = paste0("../output/viln/", outName, "/"), 
    outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), 
    returnViln = F
)

## Fig 2a - UMAP by clusID_new
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusID_new",
              cols = levels(seu.obj$sub_colz),
              pt.size = 0.5,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- cusLabels(
    plot = pi, shape = 21, size = 10, textSize = 6, 
    alpha = 0.8, labCol = c("white","black","black","black","black","white")
) + 
    NoLegend() + 
    theme(
        axis.title = element_blank(),
        panel.border = element_blank()
    )
ggsave(paste("../output/", outName, "/", "rawUMAP.png", sep = ""), width = 7, height = 7)

### Supp Fig 2a - split UMAP by cell source with original labels
seu.obj$ct <- ifelse(seu.obj$cellSource == "Tumor", 
                     paste0("TIL_", seu.obj$celltype.l3), 
                     paste0("Blood_",seu.obj$celltype.l3_pbmc))
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

# #Find centroid of each sample and cell source - not used
# as.data.frame(seu.obj@reductions$umap@cell.embeddings) %>%
#     rownames_to_column() %>%
#     left_join(rownames_to_column(seu.obj@meta.data), by = "rowname") %>%
#     group_by(name) %>%
#     summarize(
#         mean_x = mean(UMAP_1),
#         mean_y = mean(UMAP_2)
#     )
# as.data.frame(seu.obj@reductions$umap@cell.embeddings) %>%
#     rownames_to_column() %>%
#     left_join(rownames_to_column(seu.obj@meta.data), by = "rowname") %>%
#     group_by(cellSource) %>%
#     summarize(
#         mean_x = mean(UMAP_1),
#         mean_y = mean(UMAP_2)
#     )
# pi <- DimPlot(seu.obj, 
#               reduction = "umap", 
#               group.by = "name",
#               split.by = "cellSource",
#               pt.size = 0.5,
#               label = TRUE,
#               label.box = F,
#               shuffle = TRUE
# )
# pi <- formatUMAP(plot = pi)
# ggsave(paste("../output/", outName, "/", "supp_label.png", sep = ""), width = 16, height = 7)

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

# ### Add supp stacked bargraph for reviewer
# p <- stackedBar(
#     seu.obj = seu.obj, 
#     groupBy = "name", 
#     clusters = "majorID_sub",
#     flipOrder = F,
#     downSampleBy = NULL
# )
# ggsave(paste("../output/", outName, "/", "stacked.png", sep = ""), width = 10, height = 5)

### Use enrichment scoring help ID cells
#load in gene lists as a named list
modulez <- list(
    "Naïve" = c("CCR7", "LEF1", "SELL", "TCF7"),
    "TFH" = c("CXCR3", "CXCR5", "ICOS", "PDCD1", "BCL6", "MAF", "IRF4",
              "STAT3", "IL21"),
    "TH1" = c("KLRD1", "IFNGR1", "CXCR3", "CXCR6", "CCR1", "CCR5", "STAT1", 
              "STAT4", "TBX21", "TNF", "LTA", "IFNG", "IL2", "IL12RB1",
              "IL18R1", "TNFSF11", "HAVCR2"),
    "TH2" = c("CXCR4", "CCR4", "CCR8", "PTGDR2", "HAVCR1", "IL17RB", "IL33",
              "GATA3", "IRF4", "STAT6", "IL4", "IL5", "IL13", "AREG"), 
    "TREG" = c("CTLA4", "IL2RA", "FOXP3"),
    "Exhausted" = c("BTLA", "CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT"),
    "Costimulatory" = c("ICOS", "CD226", "SLAMF1", "TNFRSF14", "TNFRSF25", "TNFRSF9")
)
#run module score
seu.obj <- AddModuleScore(
    seu.obj, features = modulez, name = "_score"
)
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)
#plot the results of enrichment scores
features <- names(modulez)
ecScores <- majorDot(
    seu.obj = seu.obj, groupBy = "clusID_new", features = rev(features)
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
        seu.obj = seu.obj, groupBy = "clusID_new", features = rev(unname(unlist(x)))
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
            data = labz.df, aes(label = labz, y = 6.85, x = (modLen+1)/2),
            angle = 270, vjust = 0.5, hjust=0.5, size = 12*0.36
        ) + 
        coord_flip(ylim = c(1, 6.75), clip = "off") + 
        annotate(
            "segment", x = -Inf, y = 6.5, xend = Inf, yend = 6.5, 
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
ggsave(paste("../output/", outName, "/", outName, "_dotPlot_ecScores_2.png", sep = ""), width = 2.25, height = 7.25, scale = 2)

#Reviewer requested support for exhaustion claim

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

#Run stats on the enrichment scores
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
    
### dotPlot by keyFeats
p <- majorDot(
    seu.obj = seu.obj, groupBy = "majorID_sub",
    features = c(
        "ZNF536","LEF1","SELL",
        "CCR7", "TSHZ2",
        "CD38","IL7R","CD44",
        "TBX21", "EOMES","GZMK",
        "LGALS3","CCDC3","GATA3",
        "RORC","RORA","SMAD2","TGFB1","TNFRSF8",
        "FOXP3", "CTLA4",
        "CCL5","IFNG","BCL6","IRF1"
    )
) + 
    theme(
        legend.position = "bottom",
        axis.title.y = element_blank(),
        plot.margin = margin(7, 7, 0, 24, "pt")
    ) + 
    scale_y_discrete(position = "right")
ggsave(paste0("../output/", outName, "/", outName, "_helper_dots.png"), width = 8, height = 4)

p <- autoDot(
    seu.obj, inFile = "../output/viln/cd4/cd4_gene_list.csv", groupBy = "majorID_sub",
    MIN_LOGFOLD_CHANGE = 0.5, MIN_PCT_CELLS_EXPR_GENE = 0.1,
    filterTerm = "ENSCAFG", n_feat = 10
)
ggsave(paste0("../output/", outName, "/", outName, "_helper_dots.png"), width = 4, height = 8)

### Supp Fig 2b - CXCL13 expression by CD4 T cells
seu.obj$ct <- paste0(ifelse(seu.obj$cellSource == "TILs", "TIL_","Blood_"), 
                     seu.obj$majorID_sub)
cxcl13.df <- FetchData(seu.obj, vars = c("cellSource", "CXCL13"))
cxcl13.df <- cxcl13.df[cxcl13.df$CXCL13 > 0, ]
cxcl13.df <- cxcl13.df %>% group_by(cellSource) %>% summarize(N = n())
cxcl13.df$pct <- c(cxcl13.df$N[1]/table(seu.obj$cellSource)[1]*100, 
                   cxcl13.df$N[2]/table(seu.obj$cellSource)[2]*100)
cxcl13.df$label <- c(paste0(cxcl13.df$N[1],"/",table(seu.obj$cellSource)[1]), 
                     paste0(cxcl13.df$N[2],"/",table(seu.obj$cellSource)[2]))

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


### Fig 2b - skew plot for abundance analysis
p <- skewPlot(
    seu.obj, groupBy = "majorID_sub", yAxisLabel = "Percent of CD4 T cells",
    dout = paste0("../output/", outName), outName = outName, 
    sampleRep = "name", grepTerm = "tils", grepRes = c("Tumor","Blood")
)
ggsave(paste0("../output/", outName, "/", "skewPlot.png"), width = 6, height = 4)


seu.obj$majorID_sub <- as.factor(seu.obj$majorID_sub) # groupBy must be a factor
groupBy <- "majorID_sub"
refVal <- "name"
cellCnts <- table(seu.obj$majorID_sub, seu.obj$name)
sampleData <- as.data.frame(list(
    name = c(
        "bl_1", "bl_2", "bl_3", "bl_4", "bl_5", "bl_6", "bl_7", "bl_8", 
        "bl_9", "bl_10", "tils_5", "tils_1", "tils_2", "tils_3", "tils_4", 
        "tils_6"
    ),
    cellSource = c(rep("blood", 10), rep("tils", 6))
))

cellCnts <- cellCnts[rowSums(cellCnts == 0) == 0, ] # Only use clusters present in all horses
#set up edgeR
dge <- DGEList(counts = cellCnts, samples = sampleData)
cellSource <- factor(sampleData$cellSource, levels = c("blood","tils"))
design <- model.matrix(~1 + cellSource)
#run stats and extract res
dge <- estimateDisp(dge, design, trend = "none")
fit <- glmQLFit(dge, design, robust = TRUE, abundance.trend = FALSE)
res <- glmQLFTest(fit, coef = ncol(design))
summary(decideTests(res))

### Fig 2c - DGE analysis

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

### Fig 7d - GO GSEA of DEG
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
    limits = c(-12,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Gene ontology") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_1.png"), width = 7, height = 7)

### Fig 2e - Reactome GSEA of DEGs
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
    limits = c(-12,ceiling(max(p$data$x_axis)*1.05)), 
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
    lolli = T, filterTerm = "_CD4_",
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

features <- c("SELL", "LEF1", "TMEM154", "TNFRSF18", "TNFRSF4", "HAVCR1", "CXCL13", "LAG3","IL4I1")

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


res.df <- read.csv("../output/cd4/pseudoBulk/allCells/cd4_cluster_allCells_all_genes.csv")
res.df <- res.df[!grepl("^ENS", res.df$gene), ]
geneList_UP <- res.df %>% filter(padj < 0.1) %>% filter(log2FoldChange > 1) %>% pull(gene)
geneList_DWN <- res.df %>% filter(padj < 0.1) %>% filter(log2FoldChange < -1) %>% pull(gene)

seu.obj$cellSource <- as.factor(seu.obj$cellSource)
p <- splitDot(
    seu.obj = seu.obj, groupBy = "majorID_sub", splitBy = "cellSource",
    namedColz = setNames(c("#F8766D", "#00BFC4"),  c("Blood", "Tumor")), 
    geneList_UP = c(geneList_UP[1:20], "HAVCR1", "CXCL13", "LAG3"), geneList_DWN = c(geneList_DWN[1:20], "SELL", "LEF1", "TMEM154"), geneColz = c("red", "blue")
)
p <- p +
    theme(
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.justification = 'center'
    )
ggsave(plot = p, paste0("../output/", outName, "/", outName, "_splitDot.png"), width = 13, height = 6)


Idents(seu.obj) <- "majorID_sub"
seu.obj <- RenameIdents(seu.obj, "Treg/Tfh (c2)" = "Treg-Tfh (c2)")
seu.obj$majorID_sub <- Idents(seu.obj)
seu.obj$majorID_sub <- factor(seu.obj$majorID_sub, levels = c("Naive (c0)", "TCM (c1)", "Treg-Tfh (c2)", 
                                                              "TEM_Th2-like (c3)", "TEM (c4)", "TEM_Th1-like (c5)"))


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

clus_colz <- levels(seu.obj$sub_colz)
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
    font_colz = c("white", rep("black", 4), "white"),
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
seu.obj$clusterID <- seu.obj$clusID_new

# Export data for UCSC cell browser
ExportToCB_cus(
    seu.obj = seu.obj, dataset.name = "CD4", outDir = "../output/cb_input/", 
    markers = "../output/supplementalData/supplemental_data_1.csv", 
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
        "ZNF536","LEF1","SELL",
        "CCR7", "TSHZ2",
        "CD38","IL7R","CD44",
        "TBX21", "EOMES","GZMK",
        "LGALS3","CCDC3","GATA3",
        "RORC","RORA","SMAD2","TGFB1","TNFRSF8",
        "FOXP3", "CTLA4",
        "CCL5","IFNG","BCL6","IRF1"
    )
)


########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end CD4 T cell analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


