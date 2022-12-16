library(ggplot2)

stages_pre_post <- c("Pre", "Post", "Donor")
stages_all <- c("Pre", "treatment_5", "treatment_10", "treatment_15", "treatment_21", "Post", "followup_1m", "followup_3m", "followup_6m", "followup_12m", "Donor")
stages_labels <- c("Pre", paste0("Day ", c(5, 10, 16, 28)), "Post", paste0("Follow-up ", c(1, 3, 6, 12), "M"), "Donor")


plot_theme <- theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )


gg_alpha <- function(df, x = stages_all, alpha_metric = "richness") {
  df_rename <- df %>%
    rename(alpha_metric = alpha_metric)

  gg_richness_data <- df_rename %>%
    filter(x_axis %in% x) %>%
    mutate(x_axis = factor(x_axis, levels = stages_all, labels = stages_labels)) %>%
    mutate(
      lty = ifelse(x_axis == "Donor", NA, "solid"),
      shape = ifelse(x_axis == "Donor", NA, "pt_sample")
    ) %>%
    arrange(id)

  gg_alpha_rich <- gg_richness_data %>%
    ggplot(aes(x = x_axis, y = alpha_metric, color = id)) +
    geom_line(aes(group = id, linetype = lty), alpha = 0.6, color = "grey70", show.legend = F) +
    geom_point(aes(group = id, shape = shape), size = 2, show.legend = F) +
    geom_point(data = filter(gg_richness_data, x_axis == "Donor"), aes(group = id), size = 2, alpha = 0.6, position = position_jitter()) +
    plot_theme +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 90, vjust = 0.35, hjust = 1)
    )

  return(gg_alpha_rich)
}

gg_beta <- function(df, method) {
  if (method == "bray") {
    y_lab <- "Similarity to Donors (Bray-Curtis)"
  } else if (method == "sorensen") {
    y_lab <- "Similarity to Donors (Sørensen Coefficient)"
  } else {
    stop("method is not valid.")
  }

  df_factor <- df %>%
    mutate(x_axis = factor(x_axis, levels = stages_all, labels = stages_labels))

  gg <- df_factor %>%
    ggplot(aes(x = x_axis, y = 1 - median_dissimilarity, color = id)) +
    geom_line(aes(group = id), alpha = 0.6, color = "grey70", position = position_dodge(0.2)) +
    geom_point(aes(group = id), position = position_dodge(0.2), size = 2) +
    labs(y = y_lab, x = "") +
    plot_theme +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.075, hjust = 1)
    )

  return(gg)
}
