#### SCRIPT FOR DATASET PREPROCESSING 
# datasets for two biological contexts: Pancreatic tissue profiling and NMD knockdown experiments
# The following script performs feature selection and normalization of the data

library(Seurat)
library(ggplot2)
library(ggpubr)
library(zeroSum) #install.packages("https://github.com/rehbergT/zeroSum/raw/main/zeroSum_2.0.7.zip", repos = NULL)
library(Matrix)
library(tidyverse)
library(patchwork)


#### NMD knockdown
# Knockdown of NMD pathway, two different measurement platforms: Smart-Seq2 and Perturb-Seq

load("Datasets/features_selection_64_32_32.rdata")
# This is the list of preselected features for the Smart-Seq2 dataset consisting of: 
# - 64 positive targets (downregulated genes in treatment groups)
# - 32 negative targets (upregulated genes in treatment groups)
# - 32 controls (constant gene expression across treatment groups)
# features are stored in features_anno$gene_name

features_64_32_32 <- unique(features_anno$gene_name)
length(features_64_32_32) # 128

# Smart-Seq2: 

SmartSeqSeurat <- readRDS("Datasets/SmartSeq_Seurat.rds")

# Experimental conditions
unique(SmartSeqSeurat$Treatment)
# WT_Luc: WT cells with UPF1 overexpressin 
# SQ_Luc: integrated mutant of UPF1 (impaired NMD)
# SQ_UPF1: UPF1 knockdown and integrated mutant of UPF1 (extreme impairment)

# we will only use WT_Luc and SQ_UPF1

smart_mat <- as.matrix(SmartSeqSeurat[["SCT"]]@data)
smart_sub <- t(smart_mat[rownames(smart_mat) %in% features_64_32_32, , drop=FALSE])
# cells = rows, coloumns = genes


# Perturb-Seq:

# contains Perturbseq data as seurat object: 
load("Datasets/Perturb_seq_UPF1KD_Ctrl.rdata")

perturb_mat <- as.matrix(Perturb_seq_Ctrl_UPF1)

perturb_sub <- perturb_mat[, colnames(perturb_mat) %in% features_anno$gene_ID, drop=FALSE]

row_map <- features_anno[, c("gene_ID", "gene_name")]
row_map <- row_map[!duplicated(row_map$gene_ID), ]

common_ids <- intersect(colnames(perturb_sub), row_map$gene_ID)
perturb_sub2 <- perturb_sub[,common_ids, drop = FALSE]
colnames(perturb_sub2) <- row_map$gene_name[match(common_ids, row_map$gene_ID)]

length(colnames(perturb_sub2))

dim(perturb_sub2)

### --- RESPONSE VECTORS --- ###

map_cond_P <- read.csv("Datasets/Conditions_perturb.csv", row.names = 1)
map_cond_S <- read.csv("Datasets/Conditions_Smart_seq3.csv", row.names = 1)

cat("Perturb conditions:\n"); print(table(map_cond_P$condition))
cat("SmartSeq conditions:\n"); print(table(map_cond_S$condition))


### --- Filter smartseq specifically
conditions <-c("SQ_UPF1", "WT_Luc")
map_cond_S <- map_cond_S[map_cond_S$condition %in% conditions, , drop=FALSE] # only keep these two conditions
table(map_cond_S)

cells_S <- intersect(rownames(smart_sub), rownames(map_cond_S))
smart_sub <- smart_sub[cells_S, , drop=FALSE]
map_cond_S <- map_cond_S[cells_S, , drop=FALSE]


### ---- Align feature space between Smartseq and Perturbseq dataset
common_genes <- intersect(colnames(smart_sub), colnames(perturb_sub2))
length(common_genes)

SMARTSEQ_FINAL <- as.data.frame(smart_sub[, common_genes])
PERTURB_FINAL  <- as.data.frame(perturb_sub2[, common_genes])

cat("Final SmartSeq dims: ", dim(SMARTSEQ_FINAL), "\n")
cat("Final PerturbSeq dims:", dim(PERTURB_FINAL), "\n")

### --- Create binary vector 1=NMD inactive 0=NMD active
map_cond_S$NMD_status <- ifelse(map_cond_S$condition=="WT_Luc", 1, 0)
map_cond_P$NMD_status <- ifelse(map_cond_P$condition=="Control", 1, 0)


### --- Align with expression matrices 
y_smart <- map_cond_S[rownames(SMARTSEQ_FINAL), "NMD_status", drop=TRUE]
y_perturb <- map_cond_P[rownames(PERTURB_FINAL), "NMD_status", drop=TRUE]

### --- Check
length(y_smart); length(y_perturb)


# ---- Downsample the Perturbseq Dataset 

set.seed(093)
table(y_perturb)

i_0 <- which(y_perturb==0)
i_1 <- sample(which(y_perturb==1), length(i_0)*2)

keep <- c(i_0, i_1)

PERTURB_FINAL <- PERTURB_FINAL[keep,]
y_perturb     <- y_perturb[keep]

table(y_perturb)
table(y_smart)

names(y_smart) <- rownames(SMARTSEQ_FINAL)
names(y_perturb) <- rownames(PERTURB_FINAL)

##### ---- Z-scaling 

# this is coloumn wise standardization: 
normalize <- function(x, mean, sd) { 
  return((x - mean) / sd) 
}

# Scaling of SmartSeq:
smart_mean <- colMeans(SMARTSEQ_FINAL)
smart_sd   <- apply(SMARTSEQ_FINAL, 2, sd)
smart_scaled   <- normalize(SMARTSEQ_FINAL, smart_mean, smart_sd)


# Scaling of PerturbSeq:
pert_mean <- colMeans(PERTURB_FINAL)
pert_sd   <- apply(PERTURB_FINAL, 2, sd)
pert_scaled <- normalize(PERTURB_FINAL, pert_mean, pert_sd)



## ========== PLOTTING ========== #

col_map <- c("Active" = "#2C7BB6", "Inactive" = "#D7191C")

# ---- SmartSeq PCA
pca_s <- prcomp(smart_scaled, center = FALSE, scale. = FALSE)
ve_s  <- (pca_s$sdev^2) / sum(pca_s$sdev^2)

pc1_lab_s <- paste0("PC1 (", round(ve_s[1]*100, 1), "%)")
pc2_lab_s <- paste0("PC2 (", round(ve_s[2]*100, 1), "%)")

df_s <- data.frame(
  PC1 = pca_s$x[,1],
  PC2 = pca_s$x[,2],
  NMD_status = factor(
    y_smart[rownames(smart_scaled)],
    levels = c(0,1),
    labels = c("Active", "Inactive")))

p_s <- ggplot(df_s, aes(PC1, PC2, color = NMD_status)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = col_map) +
  theme_bw(base_size = 13) +
  labs(
    title = "NMD Smart-Seq2 (source domain)",
    x = pc1_lab_s,
    y = pc2_lab_s
  ) +
  theme(
    legend.title = element_blank())

# ---- PerturbSeq PCA
pca_p <- prcomp(pert_scaled, center = FALSE, scale. = FALSE)
ve_p  <- (pca_p$sdev^2) / sum(pca_p$sdev^2)

pc1_lab_p <- paste0("PC1 (", round(ve_p[1]*100, 1), "%)")
pc2_lab_p <- paste0("PC2 (", round(ve_p[2]*100, 1), "%)")

df_p <- data.frame(PC1 = pca_p$x[,1], PC2 = pca_p$x[,2],
  NMD_status = factor(
    y_perturb[rownames(pert_scaled)],
    levels = c(0,1),
    labels = c("Active", "Inactive")))

p_p <- ggplot(df_p, aes(PC1, PC2, color = NMD_status)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = col_map) +
  theme_bw(base_size = 13) +
  labs(
    title = "NMD Perturb-Seq (external validation)",
    x = pc1_lab_p,
    y = pc2_lab_p
  ) +
  theme(
    legend.title = element_blank())

(p_s | p_p) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal")



# Combined PCA (cross-platform)

common_genes <- intersect(colnames(smart_scaled), colnames(pert_scaled))
smart_c <- smart_scaled[, common_genes, drop = FALSE]
pert_c  <- pert_scaled[,  common_genes, drop = FALSE]

X_all <- rbind(smart_c, pert_c)

meta <- data.frame(
  platform = factor(c(rep("SmartSeq", nrow(smart_c)),
      rep("PerturbSeq", nrow(pert_c))), levels = c("SmartSeq", "PerturbSeq")),
  NMD_status = c(y_smart[rownames(smart_c)], y_perturb[rownames(pert_c)]))

meta$NMD_status <- factor(meta$NMD_status,
                          levels = c(0,1),
                          labels = c("Active", "Inactive"))

pca_all <- prcomp(X_all, center = FALSE, scale. = FALSE)
var_expl <- (pca_all$sdev^2) / sum(pca_all$sdev^2)

pc1_lab_all <- paste0("PC1 (", round(var_expl[1] * 100, 1), "%)")
pc2_lab_all <- paste0("PC2 (", round(var_expl[2] * 100, 1), "%)")

df_all <- data.frame(
  PC1 = pca_all$x[,1],
  PC2 = pca_all$x[,2],
  platform = meta$platform,
  NMD_status = meta$NMD_status)

p_all <- ggplot(df_all, aes(PC1, PC2, color = NMD_status, shape = platform)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = col_map) +
  scale_shape_manual(values = c("SmartSeq" = 4, "PerturbSeq" = 16)) +
  theme_bw(base_size = 13) +
  labs(
    title = "Combined PCA (cross-platform comparison)",
    x = pc1_lab_all,
    y = pc2_lab_all
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

p_all


