library(tidyverse)

stages_pre_post <- c("Pre", "Post", "Donor")
stages_all <- c("Pre", "treatment_5", "treatment_10", "treatment_15", "treatment_21", "Post", "followup_1m", "followup_3m", "followup_6m", "followup_12m", "Donor")
stages_labels <- c("Pre", paste0("Day ", c(5,10,16,28)), "Post", paste0("Follow-up ", c(1,3,6, 12), "M"), "Donor")


mutate_x_axis_factor <- function(df) {
  df %>% 
    mutate(x_axis = factor(x_axis, levels = c("Pre", "Post", "Donor", "followup_1m",  "followup_3m", "followup_6m",  "followup_12m")))
}


## Load data
load_metadata <- function() {
    #metadata <- read_delim("../data/metadata/csv/full_metadata.csv", delim = ";")
    metadata <- read_delim("data/metadata/csv/selected_metadata.csv", delim = ";")  
    return(metadata)
  
}

load_metaphlan <- function(taxonomic_level = "genus") {
  metaphlan <- read_delim(paste0("data/metaphlan3/MetaPhlAn3_Subsample_5000000_", taxonomic_level), delim = "\t")
  return(metaphlan)
  
}

### Metadata functions
selectMetadata <- function(metadata) {
  
  metadata_all <- metadata %>% 
  select(id, sample_barcode, LibID, stage, project, group, donor, batch_1:batch_3, fecal_donation_number, fecal_batch_date, pdai_score) %>% 
  arrange(LibID)
  
  return(metadata_all) 
}


IsolateProjectAndDonor <- function(metadata, project_filter = "MP"){
  
  metadata_donor_filter <- metadata %>% 
    pivot_longer(batch_1:batch_3, names_to = "batch", values_to = "batch_number") %>% 
    filter(!is.na(batch_number), project == project_filter) %>% 
    mutate(donor_batch = paste0(donor, "_", batch_number))
  
  metadata_in_MP_pre_post <- metadata %>% 
    mutate(donor_batch = paste0(id, "_", fecal_donation_number)) %>% 
    filter(project == project_filter | donor_batch %in% metadata_donor_filter$donor_batch) %>% 
    mutate(x_axis = if_else(stage == "inclusion", "Pre",
                            if_else(stage == "followup_30d", "Post", stage)))
  
  metadata_donors_used_in_trial <- metadata_in_MP_pre_post %>% 
    filter(project == "donor_batch") %>% 
    mutate(x_axis = "Donor",
           group = "FMT")
  
  metadata_patients_completing_trial <- metadata_in_MP_pre_post %>% 
    filter(
      project == project_filter
      #stage %in% c("inclusion", "treatment_5", "treatment_10", "treatment_15", "treatment_21", "followup_30d")
    ) 
  
  
  metadata_relevant <- bind_rows(metadata_donors_used_in_trial, metadata_patients_completing_trial)
  return(metadata_relevant)
}


isolateDonorBatchesUsed <- function(metadata, project_filter){
  patientMetadataWithDonorbatches <- metadata %>% 
    pivot_longer(batch_1:batch_3, names_to = "batch", values_to = "batch_number") %>% 
    filter(!is.na(batch_number), project == project_filter) %>% 
    mutate(donor_batch = paste0(donor, "_", batch_number))

  return(patientMetadataWithDonorbatches)
}

isolateDonorAndPatientMetadata <- function(metadata, metadata_with_donor_batches_used, project_filter){
  DonorAndPatientMetadata <- metadata %>% 
    mutate(donor_batch = paste0(id, "_", fecal_donation_number)) %>% 
    filter(project == project_filter | donor_batch %in% metadata_with_donor_batches_used$donor_batch) %>% 
    mutate(group = if_else(project == "donor_batch", "FMT", group))

    return(DonorAndPatientMetadata)
}

createXaxis <- function(df){
  dfWithXaxis <- df %>%
    mutate(x_axis = if_else(stage == "inclusion", "Pre",
                            if_else(stage == "followup_30d", "Post", stage))) %>% 
    mutate(x_axis = if_else(project == "donor_batch", "Donor", x_axis))
}



### Metaphlan functions

calculateSpeciesRichness <- function(metaphlan, filter_species = 0){
  species_count_for_each_sample <- metaphlan %>% 
    select(-NCBI_tax_id) %>% 
    pivot_longer(-clade_name, names_to = "LibID") %>% 
    group_by(LibID) %>% 
    summarise(richness = sum(value > filter_species))

  return(species_count_for_each_sample)
}

transposeMetaphlan <- function(metaphlan){
  t_metaphlan <- metaphlan %>% 
    select(-NCBI_tax_id) %>% 
    pivot_longer(-clade_name, names_to = "LibID") %>% 
    pivot_wider(names_from = clade_name, values_from = value) %>% 
    arrange(LibID) %>% 
    column_to_rownames(var = "LibID")


  return(t_metaphlan)
}



plot_theme = theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) 