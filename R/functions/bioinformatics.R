library(tidyverse)
library(ggplot2)
library(vegan)



PatientDonorHeatmap <- function(metadata_with_tax, patient, project_filter, stages = c("Pre", "treatment_5", "treatment_10", "treatment_15", "treatment_21", "Post")) {
  
  check_if_placebo_or_NA <- filter(metadata_with_tax, id == patient, project == project_filter)$group[1] == "placebo"
  
  if (check_if_placebo_or_NA == TRUE | is.na(check_if_placebo_or_NA)) {
    stop("Sample either not in project or in the placebo group")
    
  }
  
  check_if_Pre_and_Post_exists <- length(unique(filter(metadata_with_tax, id == patient, x_axis %in% c("Pre", "Post"))$x_axis)) == 2
  
  if (check_if_Pre_and_Post_exists == FALSE) {
    stop("Patient lack pre or post information")
    
  }
  
  
  
  donors_used_for_patient <- metadata_with_tax %>%
    filter(id == patient) %>% 
    pivot_longer(batch_1:batch_3, names_to = "batch", values_to = "batch_number") %>% 
    filter(!is.na(batch_number)) %>% 
    mutate(donor_batch = paste0(donor, "_", batch_number)) %>% 
    distinct(donor_batch, .keep_all = T) %>% 
    separate(stage, into = c("variable", "stage_value"), sep = "\\_") %>% 
    mutate(stage_value = as.numeric(stage_value)) %>% 
    arrange(stage_value, batch) %>% 
    select(id, LibID, variable, stage_value, donor, donor_batch)
  
  unique_donors <- donors_used_for_patient %>% 
    distinct(donor, .keep_all = T) %>% 
    mutate(donor_treatment = if_else(stage_value == 5, 1,
                                     if_else(stage_value == 10, 2,
                                             if_else(stage_value == 15, 3, 4)))) %>% 
    mutate(donor_treatment = paste0("Donor_",donor_treatment)) %>% 
    select(donor, donor_treatment) %>% 
    rename(id = donor)
  
  
  donors <- metadata_with_tax %>% 
    filter(donor_batch %in% donors_used_for_patient$donor_batch) %>% 
    full_join(unique_donors, by = "id") %>% 
    mutate(
      donor_batch = as.factor(donor_batch),
      donor_batch = fct_relevel(.f = donor_batch, donors_used_for_patient$donor_batch),
      x_axis = donor_batch,
      facet = donor_treatment
    )
  
  
  patient_metadata <- metadata_with_tax %>% 
    filter(id == patient, x_axis %in% stages) %>% 
    mutate(
      facet = patient
    )
  
  
  final_df <- bind_rows(donors, patient_metadata)
  
  TopTaxa <- final_df %>% 
    group_by(clade_name) %>% 
    summarise(relative_abundance = mean(value)) %>% 
    arrange(desc(relative_abundance)) %>% 
    head(20)
  
  gg <- final_df %>% 
    filter(clade_name %in% TopTaxa$clade_name) %>% 
    mutate(
      x_axis = factor(x_axis, levels = c(stages, donors_used_for_patient$donor_batch)),
      facet = factor(facet, levels = c(patient, unique_donors$donor_treatment)),
      clade_name = as.factor(clade_name),
      clade_name = fct_relevel(clade_name, TopTaxa$clade_name),
      value = round(value, 2)
    ) %>%
    
    ggplot(aes(x = x_axis, y = fct_rev(clade_name), label = value, fill = value))+
    geom_tile(color = "white") +
    geom_text() +
    scale_fill_gradient(low = "white", high = "red") + 
    facet_grid(~ facet, scales ="free_x", space = "free_x") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = -0.005)
    )
  
  return(gg)  
}


abundance_filter <- function(df, abundance_threshold = 0.1, verbose = F) {
  n_before <- ncol(df)
  
  df_filter <- df %>% 
    rownames_to_column("id") %>% 
    pivot_longer(-id, names_to = "clade_name", values_to = "relative_abundance") %>% 
    mutate(relative_abundance = relative_abundance / 100) %>% 
    group_by(clade_name) %>% 
    summarise(relative_abundance = max(relative_abundance)) %>% 
    filter(relative_abundance >= abundance_threshold / 100)
  
  df_filtered <- df %>% 
    select(df_filter$clade_name)
  
  n_after <- ncol(df_filtered)
  if (verbose) {
    cat(paste0("Species before: ", n_before, "\n", 
               "Species after: ", n_after, "\n",
               "Removed species: ", n_before - n_after, "\n"))
  }
  
  return(df_filtered)
  
}

hellinger_transform <- function(df, transform = transform, method = "hellinger") {
  if (transform) {
    df <- vegan::decostand(df, method = method)
  } else {
    df <- df
  }
  return(df)
}



# beta_diversity_bray <- function(metadata_with_remission, t_taxdata, transform = F, abundance_threshold = 0.1){
# 
# set.seed(1)
# metadata_patients <- metadata_with_remission %>% 
#   filter(!is.na(remission))
# 
# 
# metadata_donor <- metadata_with_remission %>% 
#   filter(x_axis == "Donor")
# 
# 
# t_metaphlan_tib <- t_taxdata %>% 
#   rownames_to_column(var = "LibID")
# 
# beta_diversity <- tibble()
# 
# for (pt_ID in unique(metadata_patients$id)) {
#   
#   patient <- metadata_patients %>% 
#     filter(id == pt_ID)
#   
#   print(patient$id)
#   
#   patient_PRE_POST <- patient %>% 
#     filter(x_axis %in% c("Pre", "Post"))
#   
#   isPatientFMT = patient_PRE_POST$group[1] == "FMT"
#   
#   if (isPatientFMT) {
#     patient_donors <- patient %>% 
#       pivot_longer(batch_1:batch_3, names_to = "batch", values_to = "batch_number") %>% 
#       filter(!is.na(batch_number)) %>% 
#       mutate(donor_batch = paste0(donor, "_", batch_number))
#     
#     donors <- metadata_donor %>% 
#       filter(donor_batch %in% patient_donors$donor_batch)
#     
#     excluded_donors <- metadata_donor %>% 
#       filter(!id %in% donors$id) %>% 
#       mutate(id_number = parse_number(id))
#     
#   } else {
#     excluded_donors <- metadata_donor %>% 
#       mutate(id_number = parse_number(id))
#   }
#   
#   
#   random_donors <- excluded_donors %>% 
#     filter(id_number %in% sample(unique(excluded_donors$id_number), 4)) %>% 
#     arrange(id) %>% 
#     group_by(id_number) %>% 
#     mutate(n_samples = row_number()) %>% 
#     filter(n_samples == sample(1:n(),1)) %>% 
#     ungroup()
#   
#   if (isPatientFMT) {
#     
#     
#     metaphlan_donor <- t_metaphlan_tib %>%
#       filter(LibID %in% donors$LibID) %>% 
#       right_join(select(donors, id, LibID), Joining, by = "LibID") %>% 
#       group_by(id) %>% 
#       summarise(across(where(is.numeric), ~mean(.))) %>% 
#       column_to_rownames("id")
#     
#   }
#   
#   
#   
#   metaphlan_random_donor <- t_metaphlan_tib %>%
#     filter(LibID %in% random_donors$LibID) %>%
#     right_join(select(random_donors, id, LibID), Joining, by = "LibID") %>% 
#     select(-LibID) %>% 
#     column_to_rownames("id")
#   
#   metaphlan_patient <- t_metaphlan_tib %>%
#     filter(LibID %in% patient_PRE_POST$LibID) %>% 
#     column_to_rownames(var = "LibID")
#   
#   
#   metaphlan_random_donor_patient <- bind_rows(metaphlan_patient, metaphlan_random_donor) %>% 
#     abundance_filter(abundance_threshold = abundance_threshold) %>% 
#     hellinger_transform(transform = transform) %>%
#     as.matrix()
#   
#   
#   
#   beta_dist_random <- vegdist(metaphlan_random_donor_patient) %>% 
#     as.matrix()
#   
#   beta_tibble_random <- beta_dist_random %>%
#     as.data.frame() %>% 
#     rownames_to_column("donor_comparison") %>%
#     as_tibble() %>% 
#     filter(row_number() > 2) %>% 
#     select(1:3) %>% 
#     pivot_longer(-donor_comparison,names_to = "LibID", values_to = "bray_curtis") %>% 
#     mutate(comparison = "random")
#   
#   
#   
#   if (isPatientFMT) {
#     metaphlan_donor_patient <- bind_rows(metaphlan_patient, metaphlan_donor) %>% 
#       abundance_filter(abundance_threshold = abundance_threshold) %>% 
#       hellinger_transform(transform = transform) %>% 
#       as.matrix()
#     
#     
#     beta_dist_donor <- vegdist(metaphlan_donor_patient) %>% 
#       as.matrix()
#     
#     beta_tibble_donor <- beta_dist_donor %>%
#       as.data.frame() %>% 
#       rownames_to_column("donor_comparison") %>%
#       as_tibble() %>% 
#       filter(row_number() > 2) %>% 
#       select(1:3) %>% 
#       pivot_longer(-donor_comparison,names_to = "LibID", values_to = "bray_curtis") %>% 
#       mutate(comparison = "real_donor")
#     
#     beta_diversity <- bind_rows(beta_diversity,beta_tibble_donor)
#     
#   }
#   
#   
#   
#   beta_diversity <- bind_rows(beta_diversity,beta_tibble_random)
#   
#   
# }
# 
# metadata_beta_diversity <- left_join(metadata_patients, beta_diversity)
# 
# return(metadata_beta_diversity)
# 
# }



# beta_diversity_sorensen <- function(metadata_with_remission, t_taxdata){
#   set.seed(1)
#   metadata_patients <- metadata_with_remission %>% 
#     filter(!is.na(remission))
#   
#   
#   metadata_donor <- metadata_with_remission %>% 
#     filter(x_axis == "Donor")
#   
#   
#   t_metaphlan_tib <- t_taxdata %>% 
#     rownames_to_column(var = "LibID")
#   
#   beta_diversity <- tibble()
#   
#   for (pt_ID in unique(metadata_patients$id)) {
#     
#     patient <- metadata_patients %>% 
#       filter(id == pt_ID)
#     
#     patient_PRE_POST <- patient %>% 
#       filter(x_axis %in% c("Pre", "Post"))
#     
#     isPatientFMT = patient_PRE_POST$group[1] == "FMT"
#     
#     if (isPatientFMT) {
#       patient_donors <- patient %>% 
#         pivot_longer(batch_1:batch_3, names_to = "batch", values_to = "batch_number") %>% 
#         filter(!is.na(batch_number)) %>% 
#         mutate(donor_batch = paste0(donor, "_", batch_number))
#       
#       donors <- metadata_donor %>% 
#         filter(donor_batch %in% patient_donors$donor_batch)
#       
#       excluded_donors <- metadata_donor %>% 
#         filter(!id %in% donors$id) %>% 
#         mutate(id_number = parse_number(id))
#       
#     } else {
#       excluded_donors <- metadata_donor %>% 
#         mutate(id_number = parse_number(id))
#     }
#     
#     
#     random_donors <- excluded_donors %>% 
#       filter(id_number %in% sample(unique(excluded_donors$id_number), 4)) %>% 
#       arrange(id) %>% 
#       group_by(id_number) %>% 
#       mutate(n_samples = row_number()) %>% 
#       filter(n_samples == sample(1:n(),1)) %>% 
#       ungroup()
#     
#     if (isPatientFMT) {
#       
#       
#       metaphlan_donor <- t_metaphlan_tib %>%
#         filter(LibID %in% donors$LibID) %>% 
#         right_join(select(donors, id, LibID)) %>% 
#         group_by(id) %>% 
#         summarise(across(where(is.numeric), ~mean(.))) %>% 
#         column_to_rownames("id")
#       
#     }
#     
#     
#     
#     metaphlan_random_donor <- t_metaphlan_tib %>%
#       filter(LibID %in% random_donors$LibID) %>%
#       right_join(select(random_donors, id, LibID)) %>% 
#       select(-LibID) %>% 
#       column_to_rownames("id")
#     
#     metaphlan_patient <- t_metaphlan_tib %>%
#       filter(LibID %in% patient_PRE_POST$LibID) %>% 
#       column_to_rownames(var = "LibID")
#     
#     
#     metaphlan_random_donor_patient <- bind_rows(metaphlan_patient, metaphlan_random_donor) %>% 
#       as.matrix()
#     
#     
#     
#     beta_dist_random <- vegdist(metaphlan_random_donor_patient, binary = T) %>% 
#       as.matrix()
#     
#     beta_tibble_random <- beta_dist_random %>%
#       as.data.frame() %>% 
#       rownames_to_column("donor_comparison") %>%
#       as_tibble() %>% 
#       filter(row_number() > 2) %>% 
#       select(1:3) %>% 
#       pivot_longer(-donor_comparison,names_to = "LibID", values_to = "sorensen_coefficient") %>% 
#       mutate(comparison = "random")
#     
#     
#     
#     if (isPatientFMT) {
#       metaphlan_donor_patient <- bind_rows(metaphlan_patient, metaphlan_donor) %>% 
#         as.matrix()
#       
#       
#       beta_dist_donor <- vegdist(metaphlan_donor_patient, binary = T) %>% 
#         as.matrix()
#       
#       beta_tibble_donor <- beta_dist_donor %>%
#         as.data.frame() %>% 
#         rownames_to_column("donor_comparison") %>%
#         as_tibble() %>% 
#         filter(row_number() > 2) %>% 
#         select(1:3) %>% 
#         pivot_longer(-donor_comparison,names_to = "LibID", values_to = "sorensen_coefficient") %>% 
#         mutate(comparison = "real_donor")
#       
#       beta_diversity <- bind_rows(beta_diversity,beta_tibble_donor)
#       
#     }
#     
#     
#     
#     beta_diversity <- bind_rows(beta_diversity,beta_tibble_random)
#     
#     
#   }
#   
#   metadata_beta_diversity <- left_join(metadata_patients, beta_diversity)
#   return(metadata_beta_diversity)
#   
# }














