


beta_diversity <- function(metadata, t_taxdata, transform = F, method = "bray", abundance_threshold = 10^-34, stages = c("Pre", "Post"), seed = 1, verbose = F){
  
  stopifnot("method must be either sorensen or bray (default: bray)" = method %in% c("bray", "sorensen"))
  
  transform <- ifelse(method == "sorensen", F, transform)
  
  set.seed(seed)
  
  t_metaphlan_tib <- t_taxdata %>% 
    rownames_to_column(var = "LibID")
  
  metadata_patients <- metadata %>% 
    filter(x_axis != "Donor") %>% 
    arrange(id)
  
  metadata_donor <- metadata %>% 
    filter(x_axis == "Donor")
  
  
  t_metaphlan_donors <- t_metaphlan_tib %>% 
    left_join(select(metadata_donor, LibID, id), by = "LibID") %>% 
    relocate(id) %>% 
    select(-LibID) %>% 
    group_by(id) %>% 
    summarise(across(where(~is.numeric(.)), ~mean(.)))
  
  
  
  
  
  beta_diversity <- tibble()
  
  for (pt_ID in unique(metadata_patients$id)) {
    
    patient <- metadata_patients %>% 
      filter(id == pt_ID)
    
    
    stage_filter = all(stages == c("Pre", "Post"))
    
    if (stage_filter) {
      hasPrePostTaxData <- patient %>% 
        filter(x_axis %in% c("Pre", "Post")) %>% 
        filter(LibID %in% t_metaphlan_tib$LibID)
      
      
      if (length(hasPrePostTaxData$id) != 2) {
        if (verbose) {
          print(paste0(pt_ID, " did not complete treatment or sample was lost"))
          
        }
        next
      } 
    }
    
    
    patient_stages <- patient %>% 
      filter(x_axis %in% stages)
    
    isPatientFMT = patient_stages$group[1] == "FMT"
    
    if (isPatientFMT) {
      patient_donors <- patient %>% 
        pivot_longer(batch_1:batch_3, names_to = "batch", values_to = "batch_number") %>% 
        filter(!is.na(batch_number)) %>% 
        mutate(donor_batch = paste0(donor, "_", batch_number))
      
      donors <- metadata_donor %>% 
        filter(donor_batch %in% patient_donors$donor_batch)
      
      excluded_donors <- metadata_donor %>% 
        filter(!id %in% donors$id) %>%
        filter(LibID %in% t_metaphlan_tib$LibID) %>% 
        mutate(id_number = parse_number(id))
      
    } else {
      excluded_donors <- metadata_donor %>% 
        filter(LibID %in% t_metaphlan_tib$LibID) %>% 
        mutate(id_number = parse_number(id))
    }
    
    
    
    
    metaphlan_patient <- t_metaphlan_tib %>%
      filter(LibID %in% patient_stages$LibID) %>% 
      column_to_rownames(var = "LibID")
    
    nstages <- length(rownames(metaphlan_patient))
    
    set.seed(1)
    
    isPatientMP <- patient$project[1] == "MP"
    if (isPatientMP) {
      random_donors <- excluded_donors %>% 
        #filter(id_number %in% sample(unique(excluded_donors$id_number), 4)) %>% 
        arrange(id) %>% 
        group_by(id_number) %>% 
        mutate(n_samples = row_number()) %>% 
        filter(n_samples == sample(1:n(),1)) %>% 
        ungroup()
      
      
      metaphlan_random_donor <- t_metaphlan_tib %>%
        filter(LibID %in% random_donors$LibID) %>%
        column_to_rownames("LibID")
      
      
      metaphlan_random_donor_patient <- bind_rows(metaphlan_patient, metaphlan_random_donor) %>% 
        abundance_filter(abundance_threshold = abundance_threshold) %>% 
        hellinger_transform(transform = transform) %>%
        as.matrix()
      
      
      beta_dist_random <- beta_method(metaphlan_random_donor_patient, method)
      
      
      beta_tibble_random <- beta_dist_random %>%
        as.data.frame() %>% 
        rownames_to_column("donor_comparison") %>%
        as_tibble() %>% 
        filter(row_number() > nstages) %>% 
        select(1:(nstages+1)) %>% 
        pivot_longer(-donor_comparison,names_to = "LibID", values_to = "beta_diversity") %>% 
        mutate(comparison = "random")
      
      
    }
    
    
    
    
    
    if (isPatientFMT) {
      
      
      areAllDonorsInMetaphlan <- all(donors$LibID %in% t_metaphlan_tib$LibID)
      if (areAllDonorsInMetaphlan == FALSE) {
        missingDonorLib <- donors %>% 
          filter(!LibID %in% t_metaphlan_tib$LibID)
        
        multipleDonorBatches <- donors %>% 
          group_by(id) %>% 
          filter(n_distinct(donor_batch) > 1)
        
        n_missingLibs <- length(missingDonorLib$id)
        n_Batches <- length(multipleDonorBatches$id)
        
        
        # Check if donor received another batch for comparison.
        if (n_Batches - n_missingLibs == 0) {
          missingDonorID <- donors %>% 
            filter(id %in% missingDonorLib$id)
          
          metaphlan_missingDonor <- t_metaphlan_donors %>% 
            filter(id %in% missingDonorLib$id) %>% 
            column_to_rownames("LibID")
          
          metaphlan_actualDonor <- t_metaphlan_tib %>%
            filter(LibID %in% donors$LibID) %>% 
            column_to_rownames("LibID")
          
          metaphlan_donor <- bind_rows(metaphlan_actualDonor, metaphlan_missingDonor)
        } else {
          metaphlan_donor <- t_metaphlan_tib %>%
            filter(LibID %in% donors$LibID) %>%
            column_to_rownames("LibID")
        }
        
        
      } else {
        metaphlan_donor <- t_metaphlan_tib %>%
          filter(LibID %in% donors$LibID) %>%
          column_to_rownames("LibID")
      }
      
      
      
      metaphlan_donor_patient <- bind_rows(metaphlan_patient, metaphlan_donor) %>% 
        abundance_filter(abundance_threshold = abundance_threshold) %>% 
        hellinger_transform(transform = transform) %>% 
        as.matrix()
      
      
      beta_dist_donor <- beta_method(metaphlan_donor_patient, method)
      
      beta_tibble_donor <- beta_dist_donor %>%
        as.data.frame() %>% 
        rownames_to_column("donor_comparison") %>%
        as_tibble() %>% 
        filter(row_number() > nstages) %>% 
        select(1:(nstages+1)) %>% 
        pivot_longer(-donor_comparison,names_to = "LibID", values_to = "beta_diversity") %>% 
        mutate(comparison = "real_donor")
      
      beta_diversity <- bind_rows(beta_diversity,beta_tibble_donor)
      
    }
    
    
    if (isPatientMP) {
      beta_diversity <- bind_rows(beta_diversity,beta_tibble_random)
    }
    
    
  }
  
  metadata_beta_diversity <- left_join(metadata_patients, beta_diversity, by = "LibID")
  
  return(metadata_beta_diversity)
  
}

beta_method <- function(df, method) {
  if (method == "bray") {
    df <- vegdist(df) %>% 
      as.matrix()
  } 
  if (method == "sorensen") {
    df <- vegdist(df, binary = T) %>% 
      as.matrix()
  }
  return(df)
}

beta_iteration <- function(n, method, transform = T){
  beta_tib <- tibble() 
  
  for (i in 1:n) {
    
    print(paste0("Iteration: ", i, " out of ", n))
    
    b <- beta_diversity(metadata_MP, t_metaphlan, method, transform = transform, seed = i)
    
    b_tib <- tibble(iteration = i, b)
    
    beta_tib <- bind_rows(beta_tib, b_tib)
  }
  print("Done!")
  return(beta_tib)
}


summariseBetaDiversityOutput <- function(beta_df, donor_id) {
  metadata_beta_diversity_mean_donor <- beta_df %>% 
    left_join(donor_id) %>% 
    filter(!is.na(donor_comparison)) 
  
  summarise_real_donor <- metadata_beta_diversity_mean_donor %>% 
    filter(comparison == "real_donor") %>% 
    group_by(id, group, x_axis, comparison) %>%  
    summarise(beta_diversity = median(beta_diversity))
  
  summarise_random_donor <- metadata_beta_diversity_mean_donor %>% 
    filter(comparison != "real_donor") %>%
    group_by(id, group, x_axis, comparison) %>%
    summarise(beta_diversity = median(beta_diversity))
  
  
  metadata_beta_diversity_summarised <- bind_rows(summarise_real_donor, summarise_random_donor) %>% 
    mutate(comparison = if_else(comparison == "random" & group == "FMT", "Donors not received", 
                                if_else(comparison == "random" & group == "placebo", "All Donors", "Actual Donors"))
    ) %>%
    arrange(id)
  
  return(metadata_beta_diversity_summarised)
}

