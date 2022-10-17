library(tidyverse)
source("functions/metadata_functions.R")


m <- load_metadata()



m_diff <- m %>% 
  mutate(diff = paired_reads - non_human_reads) %>% 
  filter(!is.na(diff), !is.na(id))

m_diff %>% 
  group_by(project) %>% 
  summarise(
    mean_reads = mean(paired_reads),
    mean_non_hum_reads = mean(non_human_reads),
    mean_diff = mean(diff, na.rm = T),
    sd_diff = sd(diff, na.rm = T),
    diff_pct = 100 - mean_non_hum_reads / mean_reads * 100
    )
