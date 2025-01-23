# Step 3: Perform enrichment analysis on SL candidates  --------------------------------
p_load(tidyverse, GRETTA, rcompanion, doMC, broom, diptest)
registerDoMC(20)

# Paths & data ---------------------------------------------------------------
dir.create("~/Reproduce_data/") # clone github repo here.
dir.create("~/Reproduce_data/data/") 
dir.create("~/Reproduce_data/output/") 
setwd("~/Reproduce_data/") 

# Download DepMap 22Q2 source data from: https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2/, and deposit in ~/Reproduce_data/data/
gretta_data_dir <- paste0("/path/to/22Q2/data/")
gretta_output_dir <- paste0("/path/to/output/")
res_corrected <- read_csv(paste0(gretta_output_dir,"/lowCASP37_v_highCASP37_screen_result_perms_10000.csv"))

# Enrichment ---------------------------------------------------------------
SL_genes <- res_corrected %>% filter(Candidate, log2FC_by_median > 0) %>%
  mutate(ENTREZID = str_split_fixed(GeneNameID, "_",2)[,2])
All_genes <- res_corrected %>%
  mutate(ENTREZID = str_split_fixed(GeneNameID, "_", 2)[, 2])

# SL enrichment - gather all GO term mappings
SL_ego_bp <- enrichGO(
  gene = unique(SL_genes$ENTREZID),
  universe = unique(All_genes$ENTREZID),
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH", # nothing was enriched by qvalue (all q < 0.3)
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  # minGSSize = 3,
  # maxGSSize = 500,
  # qvalueCutoff  = 1,
  readable = TRUE
)

saveRDS(SL_ego_bp, paste0(gretta_output_dir, "/lowCASP37_v_highCASP37_SL_enrichment_GOBP.rds"))

SL_ego_bp@result %>%
  write_csv(paste0(gretta_output_dir, "/lowCASP37_v_highCASP37_SL_enrichment_GOBP.csv"))

# Calculate jaccard index ----------------------------------------------------------
# Functions
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return(intersection / union)
}

# SL candidates
SL_ego_bp_tbl <- SL_ego_bp@result %>%
  as_tibble() %>%
  filter(pvalue < 0.05)
jaccard_df <- expand_grid(GO1 = SL_ego_bp_tbl$ID, GO2 = SL_ego_bp_tbl$ID) %>%
  mutate(jaccard_index = NA)

All_res <- NULL
All_res <- foreach(i = 1:nrow(jaccard_df), .combine = bind_rows) %dopar% {
  if (i == 1) {
    cat(paste0("Processing ", i, " of ", nrow(jaccard_df)), "\n")
  } else if (i == nrow(jaccard_df)) {
    cat(paste0("Processing ", i, " of ", nrow(jaccard_df)), "\n")
  } else if (i %% 1000 == 0) {
    cat(paste0("Processing ", i, " of ", nrow(jaccard_df)), "\n")
  }

  GO1 <- jaccard_df$GO1[i]
  GO2 <- jaccard_df$GO2[i]

  GO1_genes <- SL_ego_bp_tbl %>%
    filter(ID == GO1) %>%
    pull(geneID) %>%
    str_split_fixed(., "/", Inf) %>%
    as.vector()
  GO2_genes <- SL_ego_bp_tbl %>%
    filter(ID == GO2) %>%
    pull(geneID) %>%
    str_split_fixed(., "/", Inf) %>%
    as.vector()

  # jaccard_df$jaccard_index[i] <- jaccard(GO1_genes, GO2_genes)
  tibble(GO1 = GO1, GO2 = GO2, jaccard_index = jaccard(GO1_genes, GO2_genes))
}
write_csv(All_res, paste0(gretta_output_dir, "/lowCASP37_v_highCASP37_SL_enrichment_GOBP_jaccard.csv"))

# Cluster by jaccard index -----------------------------------------------------------
p_load(cluster, factoextra) # clustering algorithms & visualization

# SL -----------------
jaccard_df <- read_csv(paste0(gretta_output_dir, "/lowCASP37_v_highCASP37_SL_enrichment_GOBP_jaccard.csv"))

jaccard_mat <- jaccard_df %>%
  pivot_wider(names_from = "GO2", values_from = "jaccard_index") %>%
  dplyr::select(-GO1) %>%
  as.matrix()
rownames(jaccard_mat) <- jaccard_df$GO1 %>% unique()

# Calculate optimal clustering
# Elbow method
df <- jaccard_mat
distance <- get_dist(df)
k2 <- kmeans(df, centers = 2, nstart = 25)
set.seed(123)
# function to compute total within-cluster sum of square
wss <- function(k) {
  kmeans(df, k, nstart = 10)$tot.withinss
}
# Compute and plot wss for k = 1 to k = 20
k.values <- 1:20
# extract wss for 2-20 clusters
wss_values <- map_dbl(k.values, wss)
pdf(paste0(gretta_output_dir, "/lowCASP37_v_highCASP37_SL_enrichment_GOBP_jaccard_reordered_elbow.pdf"), width = 6, height = 6)
plot(k.values, wss_values,
  type = "b", pch = 19, frame = FALSE,
  xlab = "Number of clusters K",
  ylab = "Total within-clusters sum of squares"
)
dev.off()

# Gap stat method
gap_stat <- clusGap(df,
  FUN = kmeans, nstart = 25,
  K.max = 20, B = 10000
)
pdf(paste0(gretta_output_dir, "/lowCASP37_v_highCASP37_SL_enrichment_GOBP_jaccard_reordered_gapstat.pdf"), width = 6, height = 6)
fviz_gap_stat(gap_stat)
dev.off()
n_clust <- with(gap_stat, maxSE(Tab[, "gap"], Tab[, "SE.sim"]))

# plot heatmap
col_fun <- colorRamp2(c(0, 1), c("white", "red"))
hmap <- Heatmap(
  jaccard_mat,
  col = col_fun,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  row_split = n_clust,
  column_split = n_clust,
  top_annotation = HeatmapAnnotation(foo = anno_block(
    gp = gpar(fill = 2),
    labels = c(1:n_clust),
    labels_gp = gpar(col = "white", fontsize = 10)
  )),
  left_annotation = rowAnnotation(foo = anno_block(
    gp = gpar(fill = 2),
    labels = c(1:n_clust),
    labels_gp = gpar(col = "white", fontsize = 10)
  )),
)

# Supplemental Figure S7B
pdf(paste0(gretta_output_dir, "/lowCASP37_v_highCASP37_SL_enrichment_GOBP_jaccard_reordered.pdf"), width = 7, height = 6)
draw(hmap)
dev.off()

hmap2 <- draw(hmap)
r_order <- row_order(hmap2) %>% unlist()
c_order <- column_order(hmap2) %>% unlist()
# rownames(jaccard_mat)[row_order(hmap2)]
cluster_annot <- NULL
for (i in 1:length(row_order(hmap2))) {
  # i <- 1
  len <- row_order(hmap2)[[i]] %>% length()
  cluster_annot <- c(cluster_annot, rep(i, len))
}

# output heatmap annotations
jaccard_mat[r_order, c_order] %>%
  as_tibble(rownames = "GO1") %>%
  left_join(., SL_ego_bp_tbl %>% dplyr::select(GO1 = "ID", Description, geneID)) %>%
  dplyr::select(Description, geneID, everything()) %>%
  mutate(cluster_id = cluster_annot) %>%
  dplyr::select(cluster_id, everything()) %>%
  write.csv(paste0(gretta_output_dir, "/lowCASP37_v_highCASP37_SL_enrichment_GOBP_jaccard_reordered.csv"))