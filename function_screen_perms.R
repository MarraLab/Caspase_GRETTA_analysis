# Additional functions. These functions are used in the main script
screen_perms <- function(control_samples = NULL, mutant_samples = NULL, perms = NULL){
  # subset 
  select_dep <- dep %>%
  pivot_longer(cols = matches("\\d"), names_to = "GeneNameID", values_to = "DepProb") %>%
  filter(
    !is.na(DepProb)) %>%
  mutate(CellType = case_when(
    DepMap_ID %in% mutant_samples ~ "Mutant",
    DepMap_ID %in% control_samples ~ "Control",
    TRUE ~ "Others")) %>%
  filter(CellType != "Others") %>%
  mutate(CellType = fct_relevel(CellType, "Control", "Mutant"))

  # Begin single gene analysis -----------------------------------------------
  # Begin nested loop
  # 1:length(unique(dep$GeneNameID))
  All_res <- NULL
  All_res <- foreach(each = 1:perms, .combine = bind_rows) %dopar% {
  # All_res <- foreach(each = 1:10, .combine = bind_rows) %dopar% {
    if(each == 1){
      cat(paste0("Processing ", each, " of ", perms),"\n")
    } else if(each == perms){
      cat(paste0("Processing ", each, " of ", perms),"\n")
    } else if(each%%100 == 0){
      cat(paste0("Processing ", each, " of ", perms),"\n")
    }

    # Create randomly sampled DepProb dataframe
    dummy_geneID <- "A1BG_1"
    df <- select_dep %>% 
      mutate(DepProb_randomize = sample(DepProb, size = length(DepProb), replace = FALSE)) %>%
      filter(GeneNameID == dummy_geneID) %>%
      select(-DepProb) %>%
      rename(DepProb = DepProb_randomize)

    # Begin analysis - analysis was stripped down to only measure the essentials
    if(all(df$DepProb == 0)){
      populate <- rep(0,5)
    # } else if(all(df$DepProb < 0.01)) {
    #   populate <- rep(0,5)

    } else if(all(df$DepProb == 1)){
      populate <- rep(1,5)

    } else {
      # # MWU doesn't handle na or zero's well so
      # # FOR NOW remove zeros.
      # df <- df %>% filter(!is.na(DepProb)) %>%
      #   filter(DepProb != 0)

      stats <- df %>%
        group_by(CellType) %>%
        summarize(Median = median(DepProb, na.rm = TRUE),
                  Mean = mean(DepProb, na.rm = TRUE),
                  .groups = "drop")

      if((any(is.na(stats)) != TRUE) & (nrow(stats) == 2)){

        fit_pval <- wilcox.test(DepProb ~ CellType, df,
                            paired = F,
                            alternative = "two.sided",
                            conf.int = T,
                            na.action = "na.omit")$p.value

      } else if((any(is.na(stats)) == TRUE) & (nrow(stats) == 2)){
        populate <- rep(0,5)
      }

      if((any(is.na(stats)) != TRUE) & (nrow(stats) == 2)){
        populate <- as.numeric(c(unlist(stats)[-c(1,2)], fit_pval))
      } else {
        populate <- rep(0,5)
      }
    }
    geneID <- paste0("gene_",each)
    tibble(Result = c("Control_median", "Mutant_median", "Control_mean", "Mutant_mean","Pval")) %>%
        mutate(!!sym(geneID) := populate) %>%
        pivot_longer(-Result) %>%
        pivot_wider(names_from = Result, values_from = value) %>%
        rename(GeneNameID = name) %>%
        mutate(
          GeneNameID = paste0("GenePerm_",each),
          Mutant_group = "Mutant",
          Control_group = "Control")
  } # End of for loop
}

perm_pvalue <- function(p_value, perm_pvalues) {

  if(is.na(p_value)|p_value == Inf){
    return(NA_integer_)
  }

  x <- table(perm_pvalues <= p_value)

  if(length(x) < 2){
    if(names(x) == TRUE){
      res <- 1
    } else {
      res <- 1/(x[[1]]+1)
    }
  } else {
    res <- x[["TRUE"]]/(x[["FALSE"]]+ x[["TRUE"]]+1)
  }
  return(res)
}


