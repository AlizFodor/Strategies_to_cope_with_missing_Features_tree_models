### PREPROCESSING OF THE PANCREATIC DATA

library(scRNAseq)
library(SingleCellExperiment)
library(Matrix)
library(DESeq2)
library(VennDiagram)
library(grid)
library(Seurat)
library(ggplot2)
library(tidyr)
library(scuttle)


# Human pancreas cross-platform datasets
sce_baron   <- BaronPancreasData()        # inDrop
sce_seger   <- SegerstolpePancreasData()  # Smart-seq2
sce_muraro  <- MuraroPancreasData()       # CEL-Seq2


head(sce_baron)
colnames(colData(sce_baron))

table(colData(sce_baron)$label)
table(colData(sce_seger)$label)
table(colData(sce_muraro))


# ---- Binary labeling (alpha vs beta)

# Baron
keep_b <- colData(sce_baron)[["label"]] %in% c("alpha", "beta")
sce_baron_bin <- sce_baron[, keep_b]
y_baron <- ifelse(colData(sce_baron_bin)[["label"]]=="beta", 1, 0)
table(y_baron)

# Segerstolpe
keep_s <- colData(sce_seger)[["cell type"]] %in% c("alpha cell", "beta cell")
sce_seger_bin <- sce_seger[, keep_s]
y_seger <- ifelse(colData(sce_seger_bin)[["cell type"]]=="beta cell", 1, 0)
table(y_seger)

# Muraro
keep_m <- colData(sce_muraro)[["label"]] %in% c("alpha", "beta")
sce_muraro_bin <- sce_muraro[, keep_m]
y_muraro <- ifelse(colData(sce_muraro_bin)[["label"]]=="beta", 1, 0)
table(y_muraro)


# ---- Gene overlap 
g_baron <- rownames(sce_baron_bin)
g_seger <- rownames(sce_seger_bin)
g_muraro <- rownames(sce_muraro_bin)

# clean muraro gene names
g_muraro_raw <- g_muraro
g_muraro <- sub("__chr.*$", "", g_muraro_raw)

length(g_baron)
length(g_seger)
length(g_muraro)

length(intersect(g_baron, g_seger))
length(intersect(g_seger, g_muraro))
length(intersect(g_muraro, g_baron))

# --- compare groups for pseudo bulk 

table(colData(sce_baron_bin)$donor,
      colData(sce_baron_bin)$label)

libsize_cells <- colSums(assay(sce_baron_bin, "counts"))

summary(libsize_cells)
quantile(libsize_cells, probs = c(0,0.25,0.5,0.75,0.9,0.99,1))

hist(libsize_cells,
     breaks=50,
     main="Single-cell library sizes (Baron)",
     xlab="Total counts per cell")


# --- build pseudo bulk groups ( donor x celltype)

group_b <- paste(colData(sce_baron_bin)$donor, colData(sce_baron_bin)$label, sep="__")
counts_b <- assay(sce_baron_bin, "counts")
pb_baron <- rowsum(t(counts_b), group=group_b)
pb_baron <- t(pb_baron)

# metadata: 
donor <- sub("__.*","", colnames(pb_baron))
celltype <- sub(".*__","", colnames(pb_baron))

cd_baron <- data.frame(donor=donor, celltype=factor(celltype, levels=c("alpha", "beta")), row.names=colnames(pb_baron))

cd_baron$donor <- factor(cd_baron$donor)
cd_baron$celltype <- factor(cd_baron$celltype,levels = c("alpha","beta"))

# DESeq2 Baron 
dds_baron <- DESeqDataSetFromMatrix(countData = round(pb_baron), colData=cd_baron, design = ~ donor + celltype)

dds_baron <- dds_baron[rowSums(counts(dds_baron)) > 10,]
dds_baron <- DESeq(dds_baron)

res_baron <- results(dds_baron,contrast = c("celltype","beta","alpha"))

res_baron <- as.data.frame(res_baron)
res_baron$gene <- rownames(res_baron)
head(res_baron)

sig_baron <- subset(res_baron, padj < 0.05 & abs(log2FoldChange) > 0.5)

nrow(sig_baron)

sig_baron <- sig_baron[order(-abs(sig_baron$log2FoldChange)),]
baron_markers <- sig_baron[,c("gene","log2FoldChange","padj")]


# ----- Segerstolpe (Smart-seq2)

group_s  <- paste(colData(sce_seger_bin)$individual, colData(sce_seger_bin)[["cell type"]],sep="__")

counts_s <- assay(sce_seger_bin, "counts")
pb_seger <- rowsum(t(counts_s), group = group_s)
pb_seger <- t(pb_seger)

# metadata:
donor    <- sub("__.*", "", colnames(pb_seger))
celltype <- sub(".*__", "", colnames(pb_seger))

cd_seger <- data.frame(donor = donor, celltype = factor(celltype, levels = c("alpha cell", "beta cell")),
  row.names = colnames(pb_seger))

cd_seger$donor <- factor(cd_seger$donor)
cd_seger$celltype <- factor(cd_seger$celltype, levels = c("alpha cell", "beta cell"))

# DESeq2 Seger
dds_seger <- DESeqDataSetFromMatrix(countData = round(pb_seger), colData = cd_seger,
  design    = ~ donor + celltype)

dds_seger <- dds_seger[rowSums(counts(dds_seger)) > 10, ]
dds_seger <- DESeq(dds_seger)

res_seger <- results(dds_seger, contrast = c("celltype", "beta cell", "alpha cell"))

res_seger <- as.data.frame(res_seger)
res_seger$gene <- rownames(res_seger)
head(res_seger)

sig_seger <- subset(res_seger, padj < 0.05 & abs(log2FoldChange) > 0.5)


sig_seger <- sig_seger[order(-abs(sig_seger$log2FoldChange)), ]
seger_markers <- sig_seger[, c("gene","log2FoldChange","padj")]

nrow(sig_seger)



# ---- Muraro (CEL-Seq2)

label_col_m <- "label"
donor_col_m <- "donor"
table(colData(sce_muraro)[[label_col_m]])
table(colData(sce_muraro)[[donor_col_m]])
rownames(sce_muraro) <- sub("__chr.*$", "", rownames(sce_muraro))


# subset to alpha / beta
keep <- colData(sce_muraro)[[label_col_m]] %in% c("alpha", "beta")
sce_muraro_bin <- sce_muraro[, keep]

# build pseudo bulk groups
group_m  <- paste(colData(sce_muraro_bin)[[donor_col_m]],
                  colData(sce_muraro_bin)[[label_col_m]],
                  sep="__")

counts_m <- assay(sce_muraro_bin, "counts")
pb_muraro <- rowsum(t(counts_m), group = group_m)
pb_muraro <- t(pb_muraro)

# metadata:
donor    <- sub("__.*", "", colnames(pb_muraro))
celltype <- sub(".*__", "", colnames(pb_muraro))

cd_muraro <- data.frame(donor = donor, 
                        celltype = factor(celltype, levels = c("alpha", "beta")),
                        row.names = colnames(pb_muraro))

cd_muraro$donor    <- factor(cd_muraro$donor)
cd_muraro$celltype <- factor(cd_muraro$celltype, levels = c("alpha", "beta"))

# DESeq2 Muraro
dds_muraro <- DESeqDataSetFromMatrix(countData = round(pb_muraro),
                                     colData = cd_muraro,
                                     design = ~ donor + celltype)

dds_muraro <- dds_muraro[rowSums(counts(dds_muraro)) > 10, ]
dds_muraro <- DESeq(dds_muraro)

res_muraro <- results(dds_muraro, contrast = c("celltype", "beta", "alpha"))

res_muraro <- as.data.frame(res_muraro)
res_muraro$gene <- rownames(res_muraro)
head(res_muraro)

sig_muraro <- subset(res_muraro, padj < 0.05 & abs(log2FoldChange) > 0.5)

nrow(sig_muraro)

sig_muraro <- sig_muraro[order(-abs(sig_muraro$log2FoldChange)), ]
muraro_markers <- sig_muraro[, c("gene","log2FoldChange","padj")]



# ---- Common Marker genes

baron_set  <- unique(baron_markers$gene)
seger_set  <- unique(seger_markers$gene)
muraro_set <- unique(muraro_markers$gene)


shared_genes <- Reduce(intersect, list(
  baron_markers$gene,
  seger_markers$gene,
  muraro_markers$gene))

length(shared_genes)

sce_baron_pca  <- sce_baron_bin[shared_genes, ]
sce_seger_pca  <- sce_seger_bin[shared_genes, ]
sce_muraro_pca <- sce_muraro_bin[shared_genes, ]


sce_baron_pca  <- logNormCounts(sce_baron_pca)
sce_seger_pca  <- logNormCounts(sce_seger_pca)
sce_muraro_pca <- logNormCounts(sce_muraro_pca)

X_b <- t(assay(sce_baron_pca,"logcounts"))
X_s <- t(assay(sce_seger_pca,"logcounts"))
X_m <- t(assay(sce_muraro_pca,"logcounts"))


X_all <- rbind(X_b, X_s, X_m)

dataset <- c(
  rep("Baron",  nrow(X_b)),
  rep("Seger",  nrow(X_s)),
  rep("Muraro", nrow(X_m)))

label <- c(
  colData(sce_baron_pca)$label == "beta",
  colData(sce_seger_pca)[["cell type"]] == "beta cell",
  colData(sce_muraro_pca)$label == "beta")

pc <- prcomp(X_all, scale.=TRUE)



# Baron seurat object
seu_b <- CreateSeuratObject(counts = counts(sce_baron_pca))
seu_b <- SCTransform(seu_b, vst.flavor = "v2", verbose = FALSE)

# SCTransform output to use for PCA
# - scale.data: Pearson residuals (variance-stabilized) genes x cells
X_b_sct <- t(GetAssayData(seu_b, assay = "SCT", slot = "scale.data"))

seu_s <- CreateSeuratObject(counts = counts(sce_seger_pca))
seu_s <- SCTransform(seu_s, vst.flavor="v2", verbose=FALSE)
X_s_sct <- t(GetAssayData(seu_s, assay="SCT", slot="scale.data"))

seu_m <- CreateSeuratObject(counts = counts(sce_muraro_pca))
seu_m <- SCTransform(seu_m, vst.flavor="v2", verbose=FALSE)
X_m_sct <- t(GetAssayData(seu_m, assay="SCT", slot="scale.data"))

common_genes <- Reduce(intersect, list(colnames(X_b_sct), colnames(X_s_sct), colnames(X_m_sct)))
X_all_sct <- rbind(X_b_sct[,common_genes], X_s_sct[,common_genes], X_m_sct[,common_genes])

pc <- prcomp(X_all_sct, center=TRUE, scale.=FALSE)


var_expl <- pc$sdev^2 / sum(pc$sdev^2)

pc1_lab <- paste0("PC1 (", round(100 * var_expl[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(100 * var_expl[2], 1), "%)")


dataset <- c(
  rep("Baron",  nrow(X_b_sct)),
  rep("Seger",  nrow(X_s_sct)),
  rep("Muraro", nrow(X_m_sct)))

label <- c(
  colData(sce_baron_pca)$label == "beta",
  colData(sce_seger_pca)[["cell type"]] == "beta cell",
  colData(sce_muraro_pca)$label == "beta")

dataset <- factor(dataset, levels = c("Baron","Seger","Muraro"))

p <- plot(pc$x[,1], pc$x[,2],
     col = ifelse(label, "#D7191C", "#2C7BB6"),
     pch = ifelse(dataset=="Baron",4,
                  ifelse(dataset=="Seger",1,17)),
     cex = 0.8,
     xlab = pc1_lab,
     ylab = pc2_lab,
     main = "PCA of shared alpha and beta marker genes")

p + legend("topright",
       legend = c("alpha","beta","Baron (inDrop)","Segerstolpe (Smart-Seq2)","Muraro (CEL-Seq2)"),
       col    = c("#2C7BB6","#D7191C","black","black","black"),
       pch    = c(16,16,4,1,17),
       pt.cex = 1,
       bty = "n")


###### ------ SPARSITY FILTERING 

sparsity_b <- 1 - Matrix::colMeans(X_b > 0)
sparsity_s <- 1 - Matrix::colMeans(X_s > 0)
sparsity_m <- 1 - Matrix::colMeans(X_m > 0)

# convert to percentage
sparsity_b <- sparsity_b * 100
sparsity_s <- sparsity_s * 100
sparsity_m <- sparsity_m * 100

df_sparse <- data.frame(
  sparsity = c(sparsity_b, sparsity_s, sparsity_m),
  dataset  = factor(c(
    rep("Baron",  length(sparsity_b)),
    rep("Segerstolpe", length(sparsity_s)),
    rep("Muraro", length(sparsity_m)))))


df_sparse_scatter <- rbind(data.frame(gene = seq_along(sparsity_b),
                                      sparsity = sort(sparsity_b), 
                                      dataset = "Baron"),
                           data.frame(gene = seq_along(sparsity_s), 
                                      sparsity = sort(sparsity_s), 
                                      dataset = "Segerstolpe"),
                           data.frame(gene = seq_along(sparsity_m),
                                      sparsity = sort(sparsity_m),
                                      dataset = "Muraro"))

df_sparse_scatter$dataset <- factor(df_sparse_scatter$dataset, 
                                    levels = c("Baron","Segerstolpe","Muraro"))


genes <- colnames(X_b)

# sparsity vectors
sp_b <- sparsity_b[genes]
sp_s <- sparsity_s[genes]
sp_m <- sparsity_m[genes]

df_sparse <- data.frame(
  gene = genes, Baron = sp_b, Segerstolpe = sp_s, Muraro = sp_m)

# average sparsity to define order
df_sparse$mean_sparsity <- rowMeans(df_sparse[,2:4])

# order genes globally
df_sparse <- df_sparse[order(df_sparse$mean_sparsity), ]


df_long <- pivot_longer(
  df_sparse, cols = c("Baron","Segerstolpe","Muraro"), names_to = "dataset",
  values_to = "sparsity")

df_long$gene_index <- rep(seq_len(nrow(df_sparse)), times=3)



# ---- Filter genes according to sparsity

thr <- 50
  
n_b <- sum(sp_b < thr)
n_s <- sum(sp_s < thr)
n_m <- sum(sp_m < thr)

data.frame(Dataset = c("Baron","Segerstolpe","Muraro"), Genes_under_25pct_sparsity = c(n_b, n_s, n_m))



ggplot(df_long, aes(x = sparsity, color = dataset)) +
  stat_ecdf(linewidth = 1.2) +
  geom_vline(xintercept = 50, linetype = "dashed") +
  scale_color_manual(values = c(
    "Baron" = "#1b9e77",
    "Segerstolpe" = "#d95f02",
    "Muraro" = "#7570b3"
  )) +
  labs(
    x = "Gene sparsity (% zeros)",
    y = "Cumulative fraction of genes",
    title = ""
  ) +
  theme_bw(base_size = 13) +
  theme(legend.title = element_blank())




## Filter all datasets with genes <50% sparsity

keep_50 <- genes[sp_b < thr & sp_s < thr & sp_m < thr]

length(keep_50)


X_b_50 <- X_b[, keep_50, drop=FALSE]
X_s_50 <- X_s[, keep_50, drop=FALSE]
X_m_50 <- X_m[, keep_50, drop=FALSE]

X_b_50_df <- as.data.frame(as.matrix(X_b_50))
X_s_50_df <- as.data.frame(as.matrix(X_s_50))
X_m_50_df <- as.data.frame(as.matrix(X_m_50))


saveRDS(X_b_50_df, "Datasets/df_baron.rds")
saveRDS(X_s_50_df, "Datasets/df_seger.rds")
saveRDS(X_m_50_df, "Datasets/df_muraro.rds")

saveRDS(y_baron, "Datasets/y_baron.rds")
saveRDS(y_seger, "Datasets/y_seger.rds")
saveRDS(y_muraro, "Datasets/y_muraro.rds")






## --- FIX MURARO DUPLICATE GENE NAMES FIRST ----------------------------------

rownames(sce_muraro) <- sub("__chr.*$", "", rownames(sce_muraro))

counts_m_raw <- counts(sce_muraro)
counts_m_collapsed <- rowsum(as.matrix(counts_m_raw), group = rownames(sce_muraro))

sce_muraro <- SingleCellExperiment(
  assays  = list(counts = counts_m_collapsed),
  colData = colData(sce_muraro)
)

## rebuild Muraro binary object
keep_m <- colData(sce_muraro)[["label"]] %in% c("alpha", "beta")
sce_muraro_bin <- sce_muraro[, keep_m]
y_muraro <- ifelse(colData(sce_muraro_bin)[["label"]] == "beta", 1, 0)

## rebuild shared-gene SCEs
shared_genes <- Reduce(intersect, list(
  rownames(sce_baron_bin),
  rownames(sce_seger_bin),
  rownames(sce_muraro_bin)
))

sce_baron_pca  <- sce_baron_bin[shared_genes, ]
sce_seger_pca  <- sce_seger_bin[shared_genes, ]
sce_muraro_pca <- sce_muraro_bin[shared_genes, ]

## log-normalize for sparsity filtering
sce_baron_pca  <- logNormCounts(sce_baron_pca)
sce_seger_pca  <- logNormCounts(sce_seger_pca)
sce_muraro_pca <- logNormCounts(sce_muraro_pca)

X_b <- t(assay(sce_baron_pca,  "logcounts"))
X_s <- t(assay(sce_seger_pca,  "logcounts"))
X_m <- t(assay(sce_muraro_pca, "logcounts"))

## sparsity filtering
sparsity_b <- 100 * (1 - Matrix::colMeans(X_b > 0))
sparsity_s <- 100 * (1 - Matrix::colMeans(X_s > 0))
sparsity_m <- 100 * (1 - Matrix::colMeans(X_m > 0))

genes <- colnames(X_b)

sp_b <- sparsity_b[genes]
sp_s <- sparsity_s[genes]
sp_m <- sparsity_m[genes]

thr <- 50
keep_50 <- genes[sp_b < thr & sp_s < thr & sp_m < thr]

## --- ARTIFICIAL DOWNSAMPLING ------------------------------------------------

counts_b_50 <- round(as.matrix(counts(sce_baron_pca)[keep_50, , drop = FALSE]))
counts_s_50 <- round(as.matrix(counts(sce_seger_pca)[keep_50, , drop = FALSE]))
counts_m_50 <- round(as.matrix(counts(sce_muraro_pca)[keep_50, , drop = FALSE]))

downsample_binom <- function(mat, p = 0.5) {
  mat <- round(as.matrix(mat))
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  
  out <- matrix(
    rbinom(length(mat), size = as.vector(mat), prob = p),
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  
  out[is.na(out)] <- 0
  storage.mode(out) <- "integer"
  
  ## make absolutely sure Seurat gets unique feature names
  rownames(out) <- make.unique(rownames(out))
  
  out
}

counts_down_b <- downsample_binom(counts_b_50, p = 0.5)
counts_down_s <- downsample_binom(counts_s_50, p = 0.5)
counts_down_m <- downsample_binom(counts_m_50, p = 0.5)

## remove all-zero cells
keep_cells_b <- colSums(counts_down_b) > 0
keep_cells_s <- colSums(counts_down_s) > 0
keep_cells_m <- colSums(counts_down_m) > 0

counts_down_b <- counts_down_b[, keep_cells_b, drop = FALSE]
counts_down_s <- counts_down_s[, keep_cells_s, drop = FALSE]
counts_down_m <- counts_down_m[, keep_cells_m, drop = FALSE]

y_baron_down  <- y_baron[keep_cells_b]
y_seger_down  <- y_seger[keep_cells_s]
y_muraro_down <- y_muraro[keep_cells_m]

## Seurat + SCTransform
seu_b_down <- CreateSeuratObject(counts = counts_down_b)
seu_b_down <- SCTransform(seu_b_down, vst.flavor = "v2", verbose = FALSE)
X_b_down <- t(GetAssayData(seu_b_down, assay = "SCT", slot = "scale.data"))

seu_s_down <- CreateSeuratObject(counts = counts_down_s)
seu_s_down <- SCTransform(seu_s_down, vst.flavor = "v2", verbose = FALSE)
X_s_down <- t(GetAssayData(seu_s_down, assay = "SCT", slot = "scale.data"))

seu_m_down <- CreateSeuratObject(counts = counts_down_m)
seu_m_down <- SCTransform(seu_m_down, vst.flavor = "v2", verbose = FALSE)
X_m_down <- t(GetAssayData(seu_m_down, assay = "SCT", slot = "scale.data"))


saveRDS(X_b_down, "R_objects/X_down.rds")
saveRDS(X_s_down, "R_objects/X_down_s.rds")
saveRDS(X_m_down, "R_objects/X_down_m.rds")

common_genes_down <- Reduce(intersect, list(
  colnames(X_b_down),
  colnames(X_s_down),
  colnames(X_m_down)
))

X_b_use <- X_b_down[, common_genes_down, drop = FALSE]
X_s_use <- X_s_down[, common_genes_down, drop = FALSE]
X_m_use <- X_m_down[, common_genes_down, drop = FALSE]

X_all_down <- rbind(X_b_use, X_s_use, X_m_use)
saveRDS(X_all_down, "R_objects/X_all_down.rds")

dataset_down <- c(
  rep("Baron",  nrow(X_b_use)),
  rep("Seger",  nrow(X_s_use)),
  rep("Muraro", nrow(X_m_use))
)

label_down <- c(y_baron_down, y_seger_down, y_muraro_down)
celltype_down <- ifelse(label_down == 1, "Beta", "Alpha")

pc_down <- prcomp(X_all_down, center = TRUE, scale. = FALSE)

## --- PCA PLOT ---------------------------------------------------------------

var_expl_down <- pc_down$sdev^2 / sum(pc_down$sdev^2)

pc1_lab_down <- paste0("PC1 (", round(100 * var_expl_down[1], 1), "%)")
pc2_lab_down <- paste0("PC2 (", round(100 * var_expl_down[2], 1), "%)")

df_pca_down <- data.frame(
  PC1 = pc_down$x[, 1],
  PC2 = pc_down$x[, 2],
  Dataset = factor(dataset_down, levels = c("Baron", "Seger", "Muraro")),
  Celltype = celltype_down
)

plot_pca<- ggplot(df_pca_down, aes(x = PC1, y = PC2, color = Celltype, shape = Dataset)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Alpha" = "#2C7BB6", "Beta" = "#D7191C")) +
  scale_shape_manual(values = c("Baron" = 4, "Seger" = 16, "Muraro" = 17)) +
  labs(
    x = pc1_lab_down,
    y = pc2_lab_down,
    title = "PCA after downsampling (50% reads)"
  ) +
  theme_bw(base_size = 22) +
  theme(legend.title = element_blank())


ggsave(
  filename = "pca_pancre.png",
  plot = plot_pca,
  path = "Plots/",
  width = 3500,
  height = 3500,
  units = "px",
  dpi = 300
)



saveRDS(as.data.frame(X_b_use), "Datasets/df_baron_down_50.rds")
saveRDS(as.data.frame(X_s_use), "Datasets/df_seger_down_50.rds")
saveRDS(as.data.frame(X_m_use), "Datasets/df_muraro_down_50.rds")

saveRDS(y_baron_down,  "Datasets/y_baron_down_50.rds")
saveRDS(y_seger_down,  "Datasets/y_seger_down_50.rds")
saveRDS(y_muraro_down, "Datasets/y_muraro_down_50.rds")