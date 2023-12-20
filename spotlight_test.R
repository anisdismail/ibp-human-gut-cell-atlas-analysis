library(SPOTlight)
library(rhdf5)
Sys.setenv(RETICULATE_PYTHON = "/media/gambino/students_workdir/ibp/gautam/miniconda3/envs/r_kernel/bin/python")
library(reticulate)
use_condaenv(condaenv = "r_kernel", conda  = "/media/gambino/students_workdir/ibp/gautam/miniconda3/bin/conda")
library(anndata)
library(SpatialExperiment)
library(SingleCellExperiment)
library(scater)
library(scran)
library(zellkonverter)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggcorrplot)
library(Polychrome)

spe <- read10xVisium(
  samples = "/media/gambino/students_workdir/ibp/visium_data/Slide2-3/1345_start-EP/outs",
  sample_id = "1345_start-EP",
  type = "HDF5",
  data = "filtered",
  images = "hires",
  load = FALSE
)

rownames(spe) <- rowData(spe)$symbol

adata <- read_h5ad("/media/gambino/students_workdir/ibp/new_exploded_categories.h5ad")

adult_intestine_idx <- which((adata$obs)$Region == "SmallInt", (adata$obs)$Diagnosis == "Healthy adult")
count_matrix <- adata$X
new_obs <- adata$obs
counts_int <- count_matrix[adult_intestine_idx, ]
obs_int <- new_obs[adult_intestine_idx, ]
adult_intestine_markers <- scoreMarkers(t(counts_int), obs_int$category)

mgs_fil <- lapply(names(adult_intestine_markers), function(i) {
    x <- adult_intestine_markers[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.8, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

dec <- modelGeneVar(t(counts_int))
hvg <- getTopHVGs(dec, n = 3000)

markers_df <- read.csv("/media/gambino/students_workdir/ibp/new_exploded_categories.csv")
filtered_markers <- as_tibble(markers_df) %>% 
    filter(gene_name %in% rownames(spe)) %>% 
    filter(z_score > 1.645) %>% 
    group_by(cell_type)  %>%
    slice_max(z_score, n = 100)
#markers_new <- as_tibble(markers_df)  %>% 
#    group_by(cell_type)  %>%
#    slice_max(pval, n = 10)

idx <- split(seq(nrow(counts_int)), obs_int$category)
n_cells <- 500
cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
        n_cells <- n
    sample(i, n_cells)
})

new_counts <- counts_int[unname(unlist(cs_keep)), ]
new_obs <- obs_int[unname(unlist(cs_keep)), ]
groups_to_use <- new_obs$category

counts_final <- as(t(new_counts), "CsparseMatrix")

res <- SPOTlight(
    x = counts_final,
    y = spe,
    groups = groups_to_use,
    mgs = mgs_df,
    n_top = 50,
    hvg = hvg,
    weight_id = "z_score",
    group_id = "cell_type",
    gene_id = "gene_name",
    scale = TRUE,
    verbose = TRUE)

mat <- res$mat
mod <- res$NMF

write.table(mat, "/media/gambino/students_workdir/ibp/gautam/spotlight_matrix_res_scaled_50_markers.txt")
write.table(res[[2]][colnames(spe)], "/media/gambino/students_workdir/ibp/gautam/res_ss_scaled_50_markers.txt")

cor_plot <- plotCorrelationMatrix(mat)
ggsave(plot = cor_plot, file = "/media/gambino/students_workdir/ibp/gautam/cor_plot_100.png")

ct <- colnames(mat)
pal = createPalette(length(ct),  c("#ff0000", "#00ff00", "#0000ff"))
names(pal) <- ct
#mat[mat < 0.1] <- 0

# Define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
#paletteMartin <- c(
#    "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
#    "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
#    "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")

#pal <- colorRampPalette(paletteMartin)(length(ct))
#names(pal) <- ct

scatterpie <- plotSpatialScatterpie(
    x = spe,
    y = mat,
    cell_types = colnames(mat),
    img = TRUE,
    scatterpie_alpha = 1,
    pie_scale = 0.4,
    axis = "h",
    degrees = 270
    ) +
    scale_fill_manual(
        values = pal,
        breaks = names(pal))

ggsave(plot = scatterpie, file = "/media/gambino/students_workdir/ibp/gautam/scatterpie_50_scaled.png")

