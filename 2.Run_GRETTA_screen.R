# Step 2:Run GRETTA screen -------------------------------------------------------
p_load(tidyverse, GRETTA, rcompanion, doMC, broom, diptest)
registerDoMC(20)

# Paths & data ---------------------------------------------------------------
dir.create("~/Reproduce_data/") # clone github repo here.
dir.create("~/Reproduce_data/data/") 
dir.create("~/Reproduce_data/output/") 
setwd("~/Reproduce_data/") 
source("./function_screen_perms.R")

# Download DepMap 22Q2 source data from: https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2/, and deposit in ~/Reproduce_data/data/
gretta_data_dir <- paste0("/path/to/22Q2/data/")
gretta_output_dir <- paste0("/path/to/output/")
load(paste0(gretta_data_dir, "sample_annot.rda"))
prot_df <- read_csv(paste0(gretta_output_dir, "/cell_line_protein_percentile.csv"))

# Run screen ---------------------------------------------------------------
# 2. CASP3+7 low protein expressors vs high
mutant_lines <- prot_df %>% filter(group == "bottom 10%") %>% pull(DepMap_ID)
control_lines <- prot_df %>% filter(group == "top 10%") %>% pull(DepMap_ID)

res <- GI_screen(
  mutant_id = mutant_lines,
  control_id = control_lines,
  core_num = 10,
  data_dir = gretta_data_dir,
  output_dir = gretta_output_dir,
  filename = "lowCASP37_v_highCASP37_screen_result"
)

# Run perms for correction ----------------------------------------------------------
# There are intermediate files 
res_perms <- screen_perms(control_samples = control_lines, mutant_samples = mutant_lines, perms = 10000) %>% 
  write_csv(paste0(gretta_output_dir,"/lowCASP37_v_highCASP37_screen_perms_10000.csv"))
res <- read_csv(paste0(gretta_output_dir,"/lowCASP37_v_highCASP37_screen_result.csv"))
res_perms <- read_csv(paste0(gretta_output_dir,"/lowCASP37_v_highCASP37_screen_perms_10000.csv"))

# Correct for multiple testing using permutations (Table used in Figure 7I.)
res_corrected <- res %>%
  mutate(
    Pval_perm_corrected = map_dbl(Pval, .f = perm_pvalue, res_perms$Pval),
    log2FC_by_median = ifelse(log2FC_by_median %in% c(Inf, -Inf), 1, log2FC_by_median),
    Interaction_score = -log10(Pval_perm_corrected) * (log2FC_by_median),
    Interaction_score = ifelse(is.na(Interaction_score),0, Interaction_score),
    Candidate = case_when(
      (Pval_perm_corrected < 0.05) & (abs(Control_median) >= 0.5 | abs(Mutant_median) >= 0.5)  ~ TRUE,
      TRUE ~ FALSE
    )) %>% 
  filter(Interaction_score > 0) %>% # select synthetic lethal interactions
  select(GeneNameID, GeneNames, Control_median, Mutant_median, log2FC_by_median, Interaction_score, Pval_perm_corrected, Candidate) %>%
  arrange(-Candidate, -Interaction_score) %>%
  write_csv(paste0(gretta_output_dir,"/lowCASP37_v_highCASP37_screen_result_perms_10000.csv"))

# Plot Figure 7I-----------------------------------------------------------
top5_candidates <- res_corrected %>%
  filter(Candidate) %>%
  arrange(-Interaction_score) %>%
  pull(GeneNames) %>%
  .[1:5]
  
plot_df <- res_corrected %>% 
  arrange(-Interaction_score) %>%  
  mutate(
    Rank = 1:length(Interaction_score),
    lethal_score = case_when(
      Interaction_score < 0 ~ 0,
      Interaction_score == 0 ~ 0,
      Interaction_score > 0 ~ Interaction_score,
      TRUE ~ FALSE),
    Lethal_candidate = ifelse(Candidate & (log2FC_by_median > 0), TRUE, FALSE),
    label = case_when(
      GeneNames %in% top5_candidates ~ TRUE,
      TRUE ~ FALSE)
  )

pdf(paste0(gretta_output_dir,"lowCASP37_v_highCASP37_screen_lethal_interaction_score.pdf"), width = 7, height = 5)
ggplot(plot_df, aes(x = Rank, y = lethal_score)) + 
  geom_point(aes(color = Lethal_candidate, size = ifelse(Lethal_candidate, 3, 1))) + 
  scale_color_manual(values = c("grey", "red")) +
  scale_size_identity() + scale_x_reverse() +
  ggrepel::geom_label_repel(aes(label = ifelse(label, GeneNames, "")), max.overlaps = Inf, force = 5,
    box.padding = 3, direction = "both", 
    min.segment.length = unit(0, "lines"), segment.color = "grey50", color = "black") +
  ylab("Lethal genetic interaction score")+
  theme_light() + 
  theme(text = element_text(size = 12))
dev.off()
