library(ggplot2)



gg_alpha <- function(df, alpha_metric = "richness", project_filter = "MP", stat.test = "wilcox.test") {
  
  stopifnot(project_filter %in% c("MP", "NP"))
  
  if (alpha_metric == "richness") {
    ylim = c(0,100)
    label.y = 75
  } else {
    ylim = c(0,4.5)
    label.y = 3
  }
  
  df_rename <- df %>%
    rename(alpha_metric = alpha_metric)
  
  p.label = "p.format"
  p.list = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p < 0.0001", "p < 0.001", " p < 0.01", "p < 0.05", "ns"))
  
  test <- stat.test
  
  if (project_filter == "MP") {
    gg_richness_data <- df_rename %>% 
      filter(!is.na(remission) | x_axis == "Donor", x_axis %in% c("Pre", "Post", "Donor")) %>% 
      mutate(x_axis = factor(x_axis, levels = c("Pre", "Post", "Donor"))) %>% 
      arrange(id)
    
    
    test <- stat.test
    
    gg_alpha_rich <- gg_richness_data %>% 
      ggplot(aes(x = x_axis, y = alpha_metric, color = x_axis)) + 
      geom_point(aes(group = id), position = position_dodge(0.2)) +
      geom_boxplot(outlier.shape = NA) +
      geom_line(data = filter(gg_richness_data, x_axis !="Donor"), aes(group = id), alpha = 0.6, color = "grey70", position = position_dodge(0.2)) +
      geom_point(aes(group = id), position = position_dodge(0.2), size = 2) +
      stat_compare_means(label = p.label, method = test, comparisons = list(c("Pre", "Post")), paired = T, label.y = label.y, symnum.arg = p.list) +
      stat_compare_means(data = filter(gg_richness_data, group != "placebo"), label = p.label, method = test, comparisons = list(c("Pre", "Donor"), c("Post", "Donor")), symnum.arg = p.list) +
      coord_cartesian(ylim = ylim) + 
      facet_grid(.~group, scales = "free_x", space = "free") +
      plot_theme
    
  } else {
    
    
    gg_richness_data <- df_rename %>% 
      mutate(x_axis = factor(x_axis, levels = stages_all, labels = stages_labels)) %>% 
      mutate(lty = ifelse(x_axis == "Donor", NA, "solid"),
             shape = ifelse(x_axis == "Donor", NA, "pt_sample")) %>% 
      arrange(id)
    
    gg_alpha_rich <- gg_richness_data %>% 
      ggplot(aes(x = x_axis, y = alpha_metric, color = id)) + 
      geom_line(aes(group = id, linetype = lty), alpha = 0.6, color = "grey70", show.legend = F) +
      geom_point(aes(group = id, shape = shape), size = 2, show.legend = F) +
      geom_point(data = filter(gg_richness_data, x_axis == "Donor"), aes(group = id), size = 2, alpha = 0.6, position = position_jitter()) +
      coord_cartesian(ylim = ylim) + 
      plot_theme + 
      theme(legend.position = "right",
            axis.text.x = element_text(angle = 90, vjust = 0.35, hjust = 1)) 
    
  
    }
  
  
  return(gg_alpha_rich)
  
  
}



gg_beta <- function(df, beta_metric, project_filter = "MP", stat.test = "wilcox.test") {
  
  
  stopifnot(project_filter %in% c("MP", "NP"))
  stopifnot(beta_metric %in% c("bray", "sorensen"))
  
  
  p.label = "p.format"
  tip = 0.02
  if (beta_metric == "bray") {
    y_lab = "Similarity to Donors (Bray-Curtis)"
    y.l = 0.45
    
    if (project_filter == "MP") {
      ylim = c(0,0.5)
    } else {
      ylim = c(0,0.6)
    }
    
    
  } else {
    y_lab = "Similarity to Donors (Sørensen Coefficient)"
    y.l = 0.55
    
    if (project_filter == "MP") {
      ylim = c(0,0.6)
    } else {
      ylim = c(0,0.8)
    }
    
  }
  
  
  
  if (project_filter == "MP") {
    gg <- df %>% 
      ggplot(aes(x = x_axis, y = 1 - beta_diversity, color = x_axis)) +
      geom_boxplot(outlier.shape = NA) +
      geom_line(aes(group = id), alpha = 0.6, color = "grey70", position = position_dodge(0.2)) + 
      geom_point(aes(group = id), position = position_dodge(0.2), size = 2) +
      stat_compare_means(label = p.label, method = stat.test, comparisons = list(c("Pre", "Post")), paired = T, label.y = y.l, tip.length = tip) +
      coord_cartesian(ylim = ylim) +
      labs(title = "FMT group", y = y_lab, x = "") +
      facet_grid(.~comparison) +
      plot_theme
  } else{
    gg <- df %>% 
      ggplot(aes(x = x_axis, y = 1 - beta_diversity, color = id)) +
      geom_line(aes(group = id), alpha = 0.6, color = "grey70", position = position_dodge(0.2)) + 
      geom_point(aes(group = id), position = position_dodge(0.2), size = 2) +
      coord_cartesian(ylim = ylim) +
      labs(y = y_lab, x = "") +
      plot_theme
  }
  
  return(gg)
}




