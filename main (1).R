library(Seurat)

dir <- "data/Synthetic Analyses of Single-Cell Transcriptomes from Multiple Brain Organoids and Fetal Brain/"
tmp <- Seurat::Read10X(data.dir = dir, gene.column = 1)


minibrains <- CreateSeuratObject(tmp,project = "minibrains",
                                 min.cells = 3, min.features = 200)


minibrains_norm <- NormalizeData(minibrains,
                                 normalization.method = "LogNormalize",
                                 scale.factor = 10000)

minibrains_norm_feat <- FindVariableFeatures(minibrains_norm, selection.method = "vst",
                             nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(minibrains_norm_feat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(minibrains_norm_feat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
