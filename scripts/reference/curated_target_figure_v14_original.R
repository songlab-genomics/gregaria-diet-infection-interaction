# ==============================================================================
# COMPOSITE PUBLICATION FIGURE GENERATOR (v41 - Dynamic Chi-Square & Excel v14)
# Target Dimensions: 190 mm x 200 mm (Vector PDF Output)
# Font Constraints:  >= 6 pt across all elements
# ==============================================================================

# 1. DEPENDENCY SETUP
required_packages <- c("ggplot2", "dplyr", "patchwork", "grid", "gtable", "cowplot", "readxl")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)
library(gtable)
library(cowplot)
library(readxl)

# 2. PROMPT USER FOR EXCEL FILE LOCATION (Defaults to v14)
message("\n--- Select '2026-07-26_curated-hypoth_DEGs_v14.xlsx' File ---")
excel_file_path <- file.choose()

if (!file.exists(excel_file_path)) {
  stop("Selected Excel file does not exist!")
}

# 3. COLOR DEFINITIONS & CONSTANTS
COLOR_33   <- "#C8641A"  # Burnt Orange / Brown
COLOR_50   <- "#1B6B3A"  # Forest Green
COLOR_83   <- "#5A2A82"  # Royal Purple
COLOR_UP   <- "#9E1B22"  # Upregulated (Red)
COLOR_NS   <- "#CCCCCC"  # N.S. (Gray)
COLOR_DOWN <- "#3B59B6"  # Downregulated (Blue)

X_POS_NUMERIC <- c(0.78, 1.00, 1.22)
BAR_WIDTH     <- 0.18
SCALE_FACTOR  <- 0.82    # Scales 100% total height down to 82% for top p-value cushion

# Post-processor for continuous horizontal x-axis line
add_continuous_xaxis <- function(p) {
  g <- ggplotGrob(p)
  panel_indices <- grepl("^panel", g$layout$name)
  panel_cols    <- g$layout$l[panel_indices]
  panel_rows    <- g$layout$b[panel_indices]
  
  g <- gtable::gtable_add_grob(
    g,
    grobs = linesGrob(
      x = unit(c(0, 1), "npc"), 
      y = unit(c(0, 0), "npc"), 
      gp = gpar(col = "black", lwd = 0.6)
    ),
    t = max(panel_rows), 
    l = min(panel_cols), 
    r = max(panel_cols), 
    z = Inf, clip = "off", name = "continuous-x-axis"
  )
  return(cowplot::ggdraw(g))
}

# 4. READ & PROCESS DATA FROM EXCEL (v14: 8 Categories across rows 3-10)
sheet_data <- read_excel(excel_file_path, sheet = "Sum-by-umbrella_num", col_names = FALSE)
data_rows  <- sheet_data[3:10, ]

# Helper function to perform dynamic Chi-Square test on raw count matrices
calc_chisq_pvalue <- function(mat_counts) {
  # Remove columns that are entirely zero to prevent chisq.test errors
  mat_clean <- mat_counts[, colSums(mat_counts) > 0, drop = FALSE]
  
  if (ncol(mat_clean) < 2 || nrow(mat_clean) < 2) {
    return("N.S.D.")
  }
  
  # Compute Chi-Square Test (with simulated p-value for small sample sizes)
  res <- suppressWarnings(chisq.test(mat_clean, simulate.p.value = (any(mat_clean < 5))))
  p_val <- res$p.value
  
  if (is.na(p_val) || p_val >= 0.05) {
    return("N.S.D.")
  } else if (p_val < 0.0001) {
    return("p < 0.0001")
  } else if (p_val < 0.001) {
    return("p < 0.001")
  } else {
    return(sprintf("p = %.3f", p_val))
  }
}

parse_umbrella_data <- function(rows_df) {
  within_list  <- list()
  across_list  <- list()
  pval_w_list  <- list()
  pval_a_list  <- list()
  
  for (i in 1:nrow(rows_df)) {
    cat_name   <- as.character(rows_df[i, 2])
    total_degs <- as.numeric(rows_df[i, 3])
    
    # Within-diet raw counts
    w_33_up <- as.numeric(rows_df[i, 4]);  w_33_ns <- as.numeric(rows_df[i, 5]);  w_33_dn <- as.numeric(rows_df[i, 6])
    w_50_up <- as.numeric(rows_df[i, 7]);  w_50_ns <- as.numeric(rows_df[i, 8]);  w_50_dn <- as.numeric(rows_df[i, 9])
    w_83_up <- as.numeric(rows_df[i, 10]); w_83_ns <- as.numeric(rows_df[i, 11]); w_83_dn <- as.numeric(rows_df[i, 12])
    
    # Across-diet raw counts
    a_50_33_up <- as.numeric(rows_df[i, 13]); a_50_33_ns <- as.numeric(rows_df[i, 14]); a_50_33_dn <- as.numeric(rows_df[i, 15])
    a_83_33_up <- as.numeric(rows_df[i, 16]); a_83_33_ns <- as.numeric(rows_df[i, 17]); a_83_33_dn <- as.numeric(rows_df[i, 18])
    a_83_50_up <- as.numeric(rows_df[i, 19]); a_83_50_ns <- as.numeric(rows_df[i, 20]); a_83_50_dn <- as.numeric(rows_df[i, 21])
    
    # Dynamic Chi-Squared Calculations
    mat_within <- matrix(c(
      w_33_up, w_33_ns, w_33_dn,
      w_50_up, w_50_ns, w_50_dn,
      w_83_up, w_83_ns, w_83_dn
    ), nrow = 3, byrow = TRUE)
    
    mat_across <- matrix(c(
      a_50_33_up, a_50_33_ns, a_50_33_dn,
      a_83_33_up, a_83_33_ns, a_83_33_dn,
      a_83_50_up, a_83_50_ns, a_83_50_dn
    ), nrow = 3, byrow = TRUE)
    
    pval_w_list[[i]] <- data.frame(Category = cat_name, p_label = calc_chisq_pvalue(mat_within), stringsAsFactors = FALSE)
    pval_a_list[[i]] <- data.frame(Category = cat_name, p_label = calc_chisq_pvalue(mat_across), stringsAsFactors = FALSE)
    
    calc_pct <- function(u, n, d) {
      tot <- u + n + d
      if (is.na(tot) || tot == 0) return(c(0, 100, 0))
      return(c(u, n, d) / tot * 100)
    }
    
    p_w33 <- calc_pct(w_33_up, w_33_ns, w_33_dn)
    p_w50 <- calc_pct(w_50_up, w_50_ns, w_50_dn)
    p_w83 <- calc_pct(w_83_up, w_83_ns, w_83_dn)
    
    p_a50_33 <- calc_pct(a_50_33_up, a_50_33_ns, a_50_33_dn)
    p_a83_33 <- calc_pct(a_83_33_up, a_83_33_ns, a_83_33_dn)
    p_a83_50 <- calc_pct(a_83_50_up, a_83_50_ns, a_83_50_dn)
    
    w_df <- data.frame(
      Category   = cat_name,
      Diet       = rep(c("33", "50", "83"), each = 3),
      Diet_Color = rep(c(COLOR_33, COLOR_50, COLOR_83), each = 3),
      Status     = rep(c("Upregulated", "N.S.", "Downregulated"), times = 3),
      Percentage = c(p_w33, p_w50, p_w83),
      Total_DEGs = total_degs,
      stringsAsFactors = FALSE
    )
    
    a_df <- data.frame(
      Category   = cat_name,
      Comparison = rep(c("50 vs 33", "83 vs 33", "83 vs 50"), each = 3),
      Status     = rep(c("Upregulated", "N.S.", "Downregulated"), times = 3),
      Percentage = c(p_a50_33, p_a83_33, p_a83_50),
      Total_DEGs = total_degs,
      stringsAsFactors = FALSE
    )
    
    within_list[[i]] <- w_df
    across_list[[i]] <- a_df
  }
  
  return(list(
    within = do.call(rbind, within_list),
    across = do.call(rbind, across_list),
    pval_w = do.call(rbind, pval_w_list),
    pval_a = do.call(rbind, pval_a_list)
  ))
}

parsed_data   <- parse_umbrella_data(data_rows)
df_within_all <- parsed_data$within
df_across_all <- parsed_data$across
df_pval_w_all <- parsed_data$pval_w
df_pval_a_all <- parsed_data$pval_a

# 5. HELPER PLOTTING FUNCTIONS
build_within_plot <- function(df_sub, cat_levels, n_labels, p_vals) {
  df_sub$Category <- factor(df_sub$Category, levels = names(cat_levels))
  levels(df_sub$Category) <- unname(cat_levels)
  
  p_vals$Category <- factor(p_vals$Category, levels = names(cat_levels))
  levels(p_vals$Category) <- unname(cat_levels)
  
  df_sub$Status <- factor(df_sub$Status, levels = c("Downregulated", "N.S.", "Upregulated"))
  df_sub$Diet   <- factor(df_sub$Diet, levels = c("33", "50", "83"))
  df_sub$Scaled_Percentage <- df_sub$Percentage * SCALE_FACTOR
  
  diet_map <- data.frame(
    Diet = factor(c("33", "50", "83"), levels = c("33", "50", "83")),
    x_num = X_POS_NUMERIC
  )
  df_sub <- left_join(df_sub, diet_map, by = "Diet")
  
  labels_df <- data.frame(
    Category = factor(rep(unname(cat_levels), each = 3), levels = unname(cat_levels)),
    x_pos    = rep(X_POS_NUMERIC, times = length(cat_levels)),
    diet_num = c("33", "50", "83"),
    diet_col = c(COLOR_33, COLOR_50, COLOR_83)
  )
  
  n_df <- data.frame(
    Category = factor(unname(cat_levels), levels = unname(cat_levels)),
    n_text   = n_labels
  )
  
  p <- ggplot(df_sub, aes(x = x_num, y = Scaled_Percentage, fill = Status, color = Diet_Color)) +
    geom_bar(stat = "identity", width = BAR_WIDTH, linewidth = 0.4) +
    facet_wrap(~ Category, nrow = 1) +
    scale_x_continuous(limits = c(0.55, 1.45), expand = c(0, 0)) +
    scale_fill_manual(values = c("Upregulated" = COLOR_UP, "N.S." = COLOR_NS, "Downregulated" = COLOR_DOWN)) +
    scale_color_identity() +
    scale_y_continuous(
      breaks = c(0, 50 * SCALE_FACTOR, 100 * SCALE_FACTOR),
      labels = c("0", "50", "100"),
      limits = c(-35, 100),
      expand = c(0, 0)
    ) +
    geom_text(data = p_vals, aes(x = 1.0, y = 91, label = p_label), inherit.aes = FALSE, size = 2.1, fontface = "bold") +
    geom_text(data = labels_df, aes(x = x_pos, y = -6, label = diet_num, color = diet_col),
              inherit.aes = FALSE, fontface = "bold", size = 2.2) +
    geom_text(data = n_df, aes(x = 1.0, y = -22.5, label = n_text),
              inherit.aes = FALSE, size = 2.0, color = "black") +
    labs(y = "Percentage of DEGs (%)") +
    coord_cartesian(ylim = c(0, 98), clip = "off") +
    theme_classic(base_size = 7) +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_text(size = 7, color = "black", margin = margin(r = 3)),
      axis.title.x = element_blank(),
      axis.text.y  = element_text(size = 6, color = "black"),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank(),
      axis.line.y  = element_line(color = "black", linewidth = 0.4),
      strip.background = element_blank(),
      strip.text = element_text(size = 6.2, face = "bold", color = "black", lineheight = 0.85, margin = margin(b = 2.5, t = 1), vjust = 0.25),
      panel.spacing = unit(0.02, "lines"),
      plot.margin   = margin(t = 2, r = 6, b = 32, l = 2),
      legend.position = "none"
    )
  
  return(add_continuous_xaxis(p))
}

build_across_plot <- function(df_sub, cat_levels, n_labels, p_vals) {
  df_sub$Category <- factor(df_sub$Category, levels = names(cat_levels))
  levels(df_sub$Category) <- unname(cat_levels)
  
  p_vals$Category <- factor(p_vals$Category, levels = names(cat_levels))
  levels(p_vals$Category) <- unname(cat_levels)
  
  df_sub$Status     <- factor(df_sub$Status, levels = c("Downregulated", "N.S.", "Upregulated"))
  df_sub$Comparison <- factor(df_sub$Comparison, levels = c("50 vs 33", "83 vs 33", "83 vs 50"))
  df_sub$Scaled_Percentage <- df_sub$Percentage * SCALE_FACTOR
  
  comp_map <- data.frame(
    Comparison = factor(c("50 vs 33", "83 vs 33", "83 vs 50"), levels = c("50 vs 33", "83 vs 33", "83 vs 50")),
    x_num = X_POS_NUMERIC
  )
  df_sub <- left_join(df_sub, comp_map, by = "Comparison")
  
  labels_df <- data.frame(
    Category    = factor(rep(unname(cat_levels), each = 9), levels = unname(cat_levels)),
    x_pos       = rep(rep(X_POS_NUMERIC, each = 3), times = length(cat_levels)),
    y_pos       = rep(c(-4.0, -9.0, -14.0), times = 3 * length(cat_levels)),
    label_text  = rep(c("50", "vs", "33", "83", "vs", "33", "83", "vs", "50"), times = length(cat_levels)),
    label_col   = rep(c(COLOR_50, "black", COLOR_33, COLOR_83, "black", COLOR_33, COLOR_83, "black", COLOR_50), times = length(cat_levels)),
    font_weight = rep(c("bold", "plain", "bold"), times = 3 * length(cat_levels))
  )
  
  n_df <- data.frame(
    Category = factor(unname(cat_levels), levels = unname(cat_levels)),
    n_text   = n_labels
  )
  
  p <- ggplot(df_sub, aes(x = x_num, y = Scaled_Percentage, fill = Status)) +
    geom_bar(stat = "identity", width = BAR_WIDTH, color = "#222222", linewidth = 0.3) +
    facet_wrap(~ Category, nrow = 1) +
    scale_x_continuous(limits = c(0.55, 1.45), expand = c(0, 0)) +
    scale_fill_manual(values = c("Upregulated" = COLOR_UP, "N.S." = COLOR_NS, "Downregulated" = COLOR_DOWN)) +
    scale_y_continuous(
      breaks = c(0, 50 * SCALE_FACTOR, 100 * SCALE_FACTOR),
      labels = c("0", "50", "100"),
      limits = c(-35, 100),
      expand = c(0, 0)
    ) +
    geom_text(data = p_vals, aes(x = 1.0, y = 91, label = p_label), inherit.aes = FALSE, size = 2.1, fontface = "bold") +
    geom_text(data = labels_df, aes(x = x_pos, y = y_pos, label = label_text, color = label_col, fontface = font_weight),
              inherit.aes = FALSE, size = 1.9) +
    geom_text(data = n_df, aes(x = 1.0, y = -22.5, label = n_text),
              inherit.aes = FALSE, size = 2.0, color = "black") +
    scale_color_identity() +
    labs(y = "Percentage of DEGs (%)") +
    coord_cartesian(ylim = c(0, 98), clip = "off") +
    theme_classic(base_size = 7) +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_text(size = 7, color = "transparent", margin = margin(r = 3)),
      axis.title.x = element_blank(),
      axis.text.y  = element_text(size = 6, color = "transparent"),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x  = element_blank(),
      axis.line.y  = element_line(color = "black", linewidth = 0.4),
      strip.background = element_blank(),
      strip.text = element_text(size = 6.2, face = "bold", color = "black", lineheight = 0.85, margin = margin(b = 2.5, t = 1), vjust = 0.25),
      panel.spacing = unit(0.02, "lines"),
      plot.margin   = margin(t = 2, r = 6, b = 32, l = 2),
      legend.position = "none"
    )
  
  return(add_continuous_xaxis(p))
}

# 6. ROW DATA MAPS & SUBPLOT GENERATION (DYNAMIC CHI-SQUARED)
# Row 1: Immune-related
r1_map <- c(
  "Innate Immune Recognition & Signaling" = "Innate Immune\nRecognition &\nSignaling", 
  "Humoral Immunity"                      = " \nHumoral\nImmunity", 
  "Cellular Immunity"                     = " \nCellular\nImmunity", 
  "Detoxification"                        = " \n \nDetoxification"
)
n_r1   <- c("(n=45)", "(n=11)", "(n=8)", "(n=8)")
pval_r1_w <- filter(df_pval_w_all, Category %in% names(r1_map))
pval_r1_a <- filter(df_pval_a_all, Category %in% names(r1_map))

p_r1_w <- build_within_plot(filter(df_within_all, Category %in% names(r1_map)), r1_map, n_r1, pval_r1_w)
p_r1_a <- build_across_plot(filter(df_across_all, Category %in% names(r1_map)), r1_map, n_r1, pval_r1_a)

# Row 2: Protein Anabolism
r2_map <- c(
  "Protein Synthesis & Translation"          = " \nProtein\nSynthesis &\nTranslation", 
  "Nutrient Sensing & TOR/Insulin Signaling" = "Nutrient\nSensing &\nTOR/Insulin\nSignaling", 
  "Nitrogen & Amino Acid Metabolism"         = " \nNitrogen &\nAmino Acid\nMetabolism", 
  "Storage Hexamerins & Reserves"            = " \nStorage\nHexamerins &\nReserves"
)
n_r2   <- c("(n=51)", "(n=17)", "(n=18)", "(n=5)")
pval_r2_w <- filter(df_pval_w_all, Category %in% names(r2_map))
pval_r2_a <- filter(df_pval_a_all, Category %in% names(r2_map))

p_r2_w <- build_within_plot(filter(df_within_all, Category %in% names(r2_map)), r2_map, n_r2, pval_r2_w)
p_r2_a <- build_across_plot(filter(df_across_all, Category %in% names(r2_map)), r2_map, n_r2, pval_r2_a)

# 7. BANNERS, HEADERS, AND ROW LABELS
make_row_label <- function(text) {
  ggdraw() + draw_label(text, fontface = "bold", size = 7.5, angle = 90, x = 0.5, y = 0.5)
}

lbl_row1 <- make_row_label("Immune-related")
lbl_row2 <- make_row_label("Protein Anabolism")

banner_title  <- ggdraw() + draw_label("Curated DEGs relating to immune and metabolic pathways", fontface = "bold", size = 9.5, y = 0.5)
banner_within <- ggdraw() + draw_label("Within-diet comparisons\n(Infected vs Control)", fontface = "bold", size = 7.5, y = 0.5)
banner_across <- ggdraw() + draw_label("Across-diet comparisons\n(Higher protein vs Lower protein)", fontface = "bold", size = 7.5, y = 0.5)

# Shared Legend
dummy_legend_plot <- ggplot(df_within_all[1:3, ], aes(x = Diet, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = NULL, values = c("Upregulated" = COLOR_UP, "N.S." = COLOR_NS, "Downregulated" = COLOR_DOWN)) +
  theme_classic(base_size = 7) +
  theme(legend.position = "top", legend.key.size = unit(3, "mm"), legend.text = element_text(size = 6.5))

legend_grob <- cowplot::get_legend(dummy_legend_plot)

# 8. COMPOSITE ASSEMBLY USING PATCHWORK
grid_title  <- (plot_spacer() | plot_spacer() | banner_title | plot_spacer()) + plot_layout(widths = c(0.06, 0.08, 0.80, 0.06))
grid_legend <- (plot_spacer() | plot_spacer() | wrap_elements(legend_grob) | plot_spacer()) + plot_layout(widths = c(0.06, 0.08, 0.80, 0.06))

col_banners <- (plot_spacer() | plot_spacer() | banner_within | banner_across | plot_spacer()) + 
  plot_layout(widths = c(0.06, 0.08, 0.40, 0.40, 0.06))

grid_r1 <- (plot_spacer() | lbl_row1 | p_r1_w | p_r1_a | plot_spacer()) + plot_layout(widths = c(0.06, 0.08, 0.40, 0.40, 0.06))
grid_r2 <- (plot_spacer() | lbl_row2 | p_r2_w | p_r2_a | plot_spacer()) + plot_layout(widths = c(0.06, 0.08, 0.40, 0.40, 0.06))

grid_r3 <- plot_spacer()

body_grid <- (grid_r1 / grid_r2 / grid_r3) + plot_layout(heights = c(1, 1, 1))

v41_final_figure <- (
  grid_title /
  grid_legend /
  col_banners /
  body_grid
) + plot_layout(heights = c(0.04, 0.03, 0.05, 0.88))

# 9. SAVE PROMPT
message("\n--- Select Directory and Filename to Save Output PDF ---")
save_path <- file.choose(new = TRUE)

if (!grepl("\\.pdf$", save_path, ignore.case = TRUE)) {
  save_path <- paste0(save_path, ".pdf")
}

ggsave(
  filename = save_path,
  plot = v41_final_figure,
  width = 190,
  height = 200,
  units = "mm",
  dpi = 600,
  device = "pdf"
)

message("\nSuccess! Publication-ready composite figure saved to:\n", save_path)