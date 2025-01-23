# Step 1: Identify CASP3 and CASP7 low and high expressing cell lines
library(tidyverse, GRETTA)

# Paths & data ---------------------------------------------------------------------------------
dir.create("~/Reproduce_data/") # clone github repo here.
dir.create("~/Reproduce_data/data/") 
dir.create("~/Reproduce_data/output/") 
setwd("~/Reproduce_data/") 

# Download DepMap 22Q2 source data from: https://www.bcgsc.ca/downloads/ytakemon/GRETTA/22Q2/, and deposit in ~/Reproduce_data/data/
gretta_data_dir <- paste0("/path/to/22Q2/data/")
gretta_output_dir <- paste0("/path/to/output/")
load(paste0(gretta_data_dir, "sample_annot.rda"))

# Identify CASP3 and CASP7 low and high expressing cell lines -----------------------------
prot_df <-
  extract_prot(
    input_samples = sample_annot$DepMap_ID,
    input_genes = "CASP3",
    data_dir = gretta_data_dir
  ) %>%
  rename(CASP3_prot = protein_expr) %>%
  select(DepMap_ID, CASP3_prot) %>% distinct %>%
    left_join(., extract_prot(
      input_samples = sample_annot$DepMap_ID,
      input_genes = "CASP7",
      data_dir = gretta_data_dir
    )) %>%
  rename(CASP7_prot = protein_expr) %>%
  select(DepMap_ID, CASP3_prot, CASP7_prot) %>% distinct

# Percentile 
CASP3_10pct <- quantile(prot_df$CASP3_prot, c(.1)) 
CASP7_10pct <- quantile(prot_df$CASP7_prot, c(.1)) 
CASP3_90pct <- quantile(prot_df$CASP3_prot, c(.9)) 
CASP7_90pct <- quantile(prot_df$CASP7_prot, c(.9)) 

prot_df <- prot_df %>% mutate(
  bottom_pct10 = case_when(
    (CASP3_prot <= CASP3_10pct) & (CASP7_prot <= CASP7_10pct) ~ TRUE,
    TRUE ~ FALSE
  ),
  top_pct10 = case_when(
    (CASP3_prot >= CASP3_90pct) & (CASP7_prot >= CASP7_90pct) ~ TRUE,
    TRUE ~ FALSE
  ),
  group = case_when(
    bottom_pct10 ~ "bottom 10%",
    top_pct10 ~ "top 10%",
    TRUE ~ "NA"
  )
) 

# Figure 7H
pdf(paste0(gretta_output_dir, "/cell_line_protein_percentile.pdf"), width = 6, height = 5)
ggplot(prot_df, aes(x = CASP3_prot, y = CASP7_prot)) +
  geom_point(size = 2, aes(colour = group)) +
  scale_colour_manual(values = c("red", "grey","blue"))+
  theme_bw()+
  theme(text = element_text(size = 12))
dev.off()

# Figure S7A
prot_df %>% left_join(sample_annot %>% select(DepMap_ID, stripped_cell_line_name, disease)) %>%
  select(DepMap_ID, stripped_cell_line_name, disease, everything()) %>%
  write_csv(paste0(gretta_output_dir, "/cell_line_protein_percentile.csv"))