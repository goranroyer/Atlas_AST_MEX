library(Seurat)
library(dplyr)
library(ggplot2)

#' Calculate Tau index for an expression vector
#' @param x : Numeric vector of expression across different tissues
#' @return Tau value between 0 and 1 (0 = ubiquitous, 1 = highly specific)
calculate_tau <- function(x) {
  if(all(is.na(x))) return(NA)
  if(min(x, na.rm = TRUE) < 0) return(NA)
  
  max_val <- max(x, na.rm = TRUE)
  if(max_val == 0) return(0)
  
  normalized <- 1 - (x / max_val)
  tau <- sum(normalized, na.rm = TRUE) / (length(x) - 1)
  
  return(tau)
}


#' Calculate Tau index for all genes in a Seurat object
#' @param seurat_obj Seurat object
#' @param group_by Variable to group cells (e.g., "germ_layer", ...)
#' @param assay Assay to use (default: "RNA")
#' @param slot Slot to use (default: "data" for normalized data)
#' @param min_pct Minimum percentage of cells expressing the gene (default: 0.1)
#' @param min_expr Minimum average expression to consider a gene (default: 0)
#' @return Data frame with results
calculate_tau_seurat <- function(seurat_obj, 
                                 group_by,
                                 assay = "RNA",
                                 slot = "data",
                                 min_pct = 0.1,
                                 min_expr = 0) {
  
  # Check that grouping variable exists
  if(!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Variable", group_by, "not found in metadata"))
  }
  
  # Extract expression matrix
  DefaultAssay(seurat_obj) <- assay
  expr_matrix <- GetAssayData(seurat_obj, slot = slot)
  
  # Get groups
  cell_groups <- seurat_obj@meta.data[[group_by]]
  unique_groups <- unique(cell_groups)
  
  cat("Calculating average expression per", group_by, "...\n")
  
  # Calculate average expression per group
  avg_expr <- sapply(unique_groups, function(group) {
    cells_in_group <- which(cell_groups == group)
    if(length(cells_in_group) > 0) {
      rowMeans(as.matrix(expr_matrix[, cells_in_group, drop = FALSE]))
    } else {
      rep(0, nrow(expr_matrix))
    }
  })
  
  rownames(avg_expr) <- rownames(expr_matrix)
  colnames(avg_expr) <- unique_groups
  
  cat("Calculating Tau indices...\n")
  
  # Calculate Tau for each gene
  tau_values <- apply(avg_expr, 1, calculate_tau)
  
  # Find group with maximum expression
  max_group <- apply(avg_expr, 1, function(x) {
    if(all(is.na(x)) || max(x, na.rm = TRUE) == 0) return(NA)
    colnames(avg_expr)[which.max(x)]
  })
  
  # Maximum expression
  max_expr <- apply(avg_expr, 1, max, na.rm = TRUE)
  
  # Overall average expression
  mean_expr <- rowMeans(avg_expr, na.rm = TRUE)
  
  # Calculate percentage of cells expressing the gene per group
  pct_expr <- sapply(unique_groups, function(group) {
    cells_in_group <- which(cell_groups == group)
    if(length(cells_in_group) > 0) {
      rowSums(as.matrix(expr_matrix[, cells_in_group, drop = FALSE]) > 0) / length(cells_in_group) * 100
    } else {
      rep(0, nrow(expr_matrix))
    }
  })
  colnames(pct_expr) <- paste0("pct_", unique_groups)
  
  # Create results data frame
  results <- data.frame(
    gene = rownames(avg_expr),
    tau = tau_values,
    max_group = max_group,
    max_expr = max_expr,
    mean_expr = mean_expr,
    avg_expr,
    pct_expr,
    stringsAsFactors = FALSE
  )
  
  # Filter lowly expressed genes
  results <- results %>%
    filter(max_expr > min_expr) %>%
    arrange(desc(tau))
  
  cat("Done! ", nrow(results), "genes analyzed\n")
  
  return(results)
}

#' Identify specific genes for each group
#' @param tau_results Results from calculate_tau_seurat
#' @param tau_threshold Tau threshold to consider a gene as specific (default: 0.8)
#' @param top_n Number of top genes to return per group (default: NULL = all)
#' @return List of specific genes per group
get_specific_genes <- function(tau_results, tau_threshold = 0.8, top_n = NULL) {
  
  # Filter by Tau threshold
  specific <- tau_results %>%
    filter(tau > tau_threshold) %>%
    filter(!is.na(max_group))
  
  # Group by max_group
  specific_by_group <- split(specific, specific$max_group)
  
  # Optional: keep only top N per group
  if(!is.null(top_n)) {
    specific_by_group <- lapply(specific_by_group, function(df) {
      df %>% arrange(desc(tau)) %>% head(top_n)
    })
  }
  
  return(specific_by_group)
}

#' #' Visualize Tau value distribution
#' #' @param tau_results Results from calculate_tau_seurat
#' #' @param tau_threshold Tau threshold to display (default: 0.8)
#' plot_tau_distribution <- function(tau_results, tau_threshold = 0.8) {
#'   
#'   p <- ggplot(tau_results, aes(x = tau)) +
#'     geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
#'     geom_vline(xintercept = tau_threshold, color = "red", linetype = "dashed", size = 1) +
#'     annotate("text", x = tau_threshold + 0.05, y = Inf, 
#'              label = paste("Threshold =", tau_threshold), 
#'              vjust = 2, color = "red", size = 4) +
#'     labs(title = "Tau Index Distribution",
#'          x = "Tau (0 = ubiquitous, 1 = specific)",
#'          y = "Number of genes") +
#'     theme_minimal() +
#'     theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
#'   
#'   print(p)
#'   return(p)
#' }
#' 
#' #' Count specific genes per group
#' #' @param tau_results Results from calculate_tau_seurat
#' #' @param tau_threshold Tau threshold (default: 0.8)
#' plot_specific_genes_count <- function(tau_results, tau_threshold = 0.8) {
#'   
#'   specific <- tau_results %>%
#'     filter(tau > tau_threshold) %>%
#'     filter(!is.na(max_group)) %>%
#'     count(max_group, name = "n_genes")
#'   
#'   p <- ggplot(specific, aes(x = reorder(max_group, n_genes), y = n_genes)) +
#'     geom_bar(stat = "identity", fill = "coral", color = "black") +
#'     coord_flip() +
#'     labs(title = paste("Number of specific genes (Tau >", tau_threshold, ")"),
#'          x = "Group",
#'          y = "Number of genes") +
#'     theme_minimal() +
#'     theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
#'   
#'   print(p)
#'   return(p)
#' }
#' 
#' 
#' #' Density plot comparing Tau distributions across groups
#' #' @param tau_results Results from calculate_tau_seurat
#' plot_tau_by_group <- function(tau_results) {
#'   
#'   specific_genes <- tau_results %>%
#'     filter(tau > 0.5)
#'   
#'   p <- ggplot(specific_genes, aes(x = tau, color = max_group)) +
#'     geom_density(linewidth = 1.2, fill = NA) +
#'     labs(title = "Tau distribution for genes by max expression group",
#'          x = "Tau",
#'          y = "Density",
#'          color = "Group") +
#'     theme_minimal() +
#'     theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
#'   
#'   print(p)
#'   p
#' }

