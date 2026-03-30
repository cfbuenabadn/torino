library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(RColorBrewer)

# Main wrapper function for plotting factors
plot_factors_wrapper <- function(EF, gene='NULL', strand = 'plus', 
                                colores = NULL, figsize = NULL, lm = TRUE, 
                                factor_list = NULL, ax = NULL, xlim = NULL, 
                                plot_coords = FALSE, lwidth = 3, mute_y_axis = FALSE) {
  
  # Extract EF data (excluding first and last rows)
#   EF <- as.data.frame(rds[[ebpmf]][['EF_smooth']][2:(nrow(rds[[ebpmf]][['EF_smooth']])-1), ])
#   rownames(EF) <- rds[[ebpmf]][['coords']][2:(length(rds[[ebpmf]][['coords']])-1)]
  coords <- rownames(EF)
  
  # Set xlim if not provided
  if (is.null(xlim)) {
    coord_parts1 <- strsplit(coords[1], ":")[[1]]
    coord_parts2 <- strsplit(coords[length(coords)], ":")[[1]]
    xlim1 <- as.numeric(coord_parts1[2]) - 1000
    xlim2 <- as.numeric(coord_parts2[2]) + 1000
    xlim <- c(xlim1, xlim2)
  }
  
  # Set column names
  colnames(EF) <- paste0('factor', 1:ncol(EF))
  
  # Filter factors if factor_list provided
  if (!is.null(factor_list)) {
    EF <- EF[, factor_list, drop = FALSE]
    colnames(EF) <- paste0('factor', 1:ncol(EF))
  }
  
  # Apply linear model transformation if requested
  for (factor in colnames(EF)) {
    if (lm) {
      EF[, factor] <- factor_lm(EF[, factor], strand)
    }
  }
  
  # Normalize and plot
  EF_norm <- sweep(EF, 2, apply(EF, 2, quantile, 0.99), FUN = "/")
  plot_factor_tracks(EF_norm, gene, ncol(EF), colores = colores, 
                    figsize = figsize, ax = ax, xlim = xlim, 
                    plot_coords = plot_coords, mute_y_axis = mute_y_axis, 
                    lwidth = lwidth)
}







# # Main wrapper function for plotting factors
# plot_factors_wrapper <- function(rds, gene, ebpmf = 'ebpmf_10', strand = 'plus', 
#                                 colores = NULL, figsize = NULL, lm = TRUE, 
#                                 factor_list = NULL, ax = NULL, xlim = NULL, 
#                                 plot_coords = FALSE, lwidth = 3, mute_y_axis = FALSE) {
  
#   # Extract EF data (excluding first and last rows)
#   EF <- as.data.frame(rds[[ebpmf]][['EF_smooth']][2:(nrow(rds[[ebpmf]][['EF_smooth']])-1), ])
#   rownames(EF) <- rds[[ebpmf]][['coords']][2:(length(rds[[ebpmf]][['coords']])-1)]
  
#   # Set xlim if not provided
#   if (is.null(xlim)) {
#     coord_parts1 <- strsplit(rds[[ebpmf]][['coords']][2], ":")[[1]]
#     coord_parts2 <- strsplit(rds[[ebpmf]][['coords']][length(rds[[ebpmf]][['coords']])-1], ":")[[1]]
#     xlim1 <- as.numeric(coord_parts1[2]) - 1000
#     xlim2 <- as.numeric(coord_parts2[2]) + 1000
#     xlim <- c(xlim1, xlim2)
#   }
  
#   # Set column names
#   colnames(EF) <- paste0('factor', 1:ncol(EF))
  
#   # Filter factors if factor_list provided
#   if (!is.null(factor_list)) {
#     EF <- EF[, factor_list, drop = FALSE]
#     colnames(EF) <- paste0('factor', 1:ncol(EF))
#   }
  
#   # Apply linear model transformation if requested
#   for (factor in colnames(EF)) {
#     if (lm) {
#       EF[, factor] <- factor_lm(EF[, factor], strand)
#     }
#   }
  
#   # Normalize and plot
#   EF_norm <- sweep(EF, 2, apply(EF, 2, quantile, 0.99), FUN = "/")
#   plot_factor_tracks(EF_norm, gene, ncol(EF), colores = colores, 
#                     figsize = figsize, ax = ax, xlim = xlim, 
#                     plot_coords = plot_coords, mute_y_axis = mute_y_axis, 
#                     lwidth = lwidth)
# }











# Helper function for linear model (placeholder - needs implementation based on your specific needs)
factor_lm <- function(y, strand) {
  # This is a placeholder - implement based on your specific linear model requirements
  return(as.numeric(y))
}

# Main function for plotting factor tracks
plot_factor_tracks <- function(EF, gene, K, title = NULL, fill = TRUE, smooth = FALSE, 
                              colores = NULL, q = 0.99, figsize = NULL, ax = NULL, 
                              xlim = NULL, plot_coords = NULL, lwidth = NULL, 
                              mute_y_axis = FALSE) {
  
  # Normalize data
  EF_norm <- pmin(sweep(EF, 2, apply(EF, 2, quantile, q), FUN = "/"), 1)
  
  # Set colors
  if (is.null(colores)) {
    if (K <= 5) {
      colores <- c('blue', 'orange', 'green', 'goldenrod', 'red')
    } else {
      colores <- rainbow(K)
    }
  }
  
  # Set figure size
  if (K <= 3) {
    S <- K
  } else {
    S <- K * 0.8
  }
  
  if (is.null(figsize)) {
    figsize <- c(15, S)
  }
  
  # Extract coordinates
  coord_parts_start <- strsplit(rownames(EF)[1], ":")[[1]]
  coord_parts_end <- strsplit(rownames(EF)[nrow(EF)], ":")[[1]]
  start <- as.numeric(coord_parts_start[2])
  end <- as.numeric(coord_parts_end[2])
  length_coords <- nrow(EF)
  
  coords <- seq(start, end, length.out = length_coords)
  
  # Create plots
  plot_list <- list()
  
  for (i in 1:K) {
    factor_name <- paste0('factor', i)
    scaled_y <- EF_norm[, factor_name]
    
    # Create data frame for ggplot
    plot_data <- data.frame(
      x = coords,
      y = scaled_y,
      factor = factor_name
    )
    
    p <- ggplot(plot_data, aes(x = x, y = y)) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")
      )
    
    if (fill) {
      p <- p + geom_ribbon(aes(ymin = 0, ymax = y), 
                          fill = colores[i], alpha = 0.3)
    }
    
    p <- p + geom_line(color = colores[i], alpha = 0.9, 
                      size = ifelse(is.null(lwidth), 1, lwidth/2))
    
    if (i < K || !plot_coords) {
      p <- p + theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
    }
    
    if (mute_y_axis) {
      p <- p + theme(axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.line.y = element_blank())
    }
    
    if (!is.null(xlim)) {
      p <- p + xlim(xlim[1], xlim[2])
    }
    
    plot_list[[i]] <- p
  }
  
  # Arrange plots
  if (!is.null(title)) {
    final_plot <- grid.arrange(grobs = plot_list, ncol = 1, 
                              top = textGrob(title, gp = gpar(fontsize = 12)))
  } else {
    final_plot <- grid.arrange(grobs = plot_list, ncol = 1)
  }
  
  return(final_plot)
}

# Function for plotting isoform annotations
plot_isoform_annotations <- function(annotation_exons, gene, colores = NULL, 
                                   start = NULL, end = NULL, figsize = NULL, 
                                   lwidth = 5, iso_order = NULL, axes = NULL, 
                                   xlim = NULL) {
  
  gene_exons <- annotation_exons[annotation_exons$gene_id == gene, ]
  
  if (is.null(iso_order)) {
    isoforms <- sort(unique(gene_exons$transcript_id))
  } else {
    isoforms <- paste(gene, iso_order, sep = ".")
  }
  
  isoform_dict <- list()
  for (i in seq_along(isoforms)) {
    isoform_name <- paste0('isoform_', i)
    df <- gene_exons[gene_exons$transcript_id == isoforms[i], ]
    df$transcript_id <- paste(gene, isoform_name, sep = ".")
    isoform_dict[[isoform_name]] <- list(df = df)
  }
  
  chrom <- annotation_exons$chrom[1]
  
  if (is.null(start)) {
    start <- as.character(min(gene_exons$start) - 1000)
  }
  if (is.null(end)) {
    end <- as.character(max(gene_exons$end) + 1000)
  }
  
  coords <- c(paste(chrom, start, sep = ":"), paste(chrom, end, sep = ":"))
  
  plot_gene_isoforms(isoform_dict, coords, color_list = colores, 
                    figsize = figsize, lwidth = lwidth, axes = axes, xlim = xlim)
}

# Function for plotting gene isoforms
plot_gene_isoforms <- function(isoforms_dict, coordinates, color_list = NULL, 
                              axes = NULL, figsize = NULL, lwidth = 5, xlim = NULL) {
  
  if (is.null(xlim)) {
    xlim1 <- as.numeric(strsplit(coordinates[1], ":")[[1]][2])
    xlim2 <- as.numeric(strsplit(coordinates[length(coordinates)], ":")[[1]][2])
  } else {
    xlim1 <- xlim[1]
    xlim2 <- xlim[2]
  }
  
  if (is.null(color_list)) {
    color_list <- rainbow(length(isoforms_dict))
  }
  
  K <- length(isoforms_dict)
  
  if (is.null(figsize)) {
    figsize <- c(20, 3)
  }
  
  plot_list <- list()
  
  for (i in 1:K) {
    isoform_df <- isoforms_dict[[paste0('isoform_', i)]][['df']]
    color <- color_list[i]
    
    p <- plot_isoform(isoform_df, color, xlim1, xlim2, lwidth)
    plot_list[[i]] <- p
  }
  
  final_plot <- grid.arrange(grobs = plot_list, ncol = 1)
  return(final_plot)
}

# Function for plotting individual isoform
plot_isoform <- function(isoform_df, color, xlim1, xlim2, lwidth) {
  
  plot_data <- data.frame()
  connection_data <- data.frame()
  
  first_end <- NULL
  last_start <- NULL
  
  for (i in 1:nrow(isoform_df)) {
    row <- isoform_df[i, ]
    start <- as.numeric(row$start)
    end <- as.numeric(row$end)
    
    if (i == 1) {
      first_end <- end
    }
    if (i == nrow(isoform_df)) {
      last_start <- start
    }
    
    # Add exon rectangle data
    exon_data <- data.frame(
      xmin = start,
      xmax = end,
      ymin = 0,
      ymax = 1,
      type = "exon"
    )
    plot_data <- rbind(plot_data, exon_data)
  }
  
  # Add connection line data
  if (!is.null(first_end) && !is.null(last_start)) {
    connection_data <- data.frame(
      x = c(first_end, last_start),
      y = c(0.5, 0.5)
    )
  }
  
  p <- ggplot() +
    geom_rect(data = plot_data, 
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = color, alpha = 0.7) +
    theme_void() +
    xlim(xlim1, xlim2) +
    ylim(0, 1)
  
  if (nrow(connection_data) > 0) {
    p <- p + geom_line(data = connection_data, 
                      aes(x = x, y = y), 
                      color = color, 
                      size = lwidth/2)
  }
  
  return(p)
}


plot_sorted_bar <- function(numbers, color, ylim = NULL, axes = TRUE) {
  # Sort the vector
  sorted_numbers <- sort(numbers)
  
  # Compute the median
  median_value <- median(sorted_numbers, na.rm = TRUE)
  
  # Plot the barplot
  barplot(
    sorted_numbers,
    col = color,
    border = color,
    space = 0,     # Thin bars
    xaxs = "i",    # Flush bars to edges
    ylim = ylim,
    axes = axes,
    names.arg = ""
  )
  
  # Add median line and marker
  abline(h = median_value, col = "black", lty = 2)  # Dotted line
  points(x = length(sorted_numbers) / 2, y = median_value, col = "black", pch = 16)  # Dot
  usr <- par("usr")
  segments(x0 = usr[1], y0 = usr[3], x1 = usr[2], y1 = usr[3])  # bottom border
  segments(x0 = usr[1], y0 = usr[3], x1 = usr[1], y1 = usr[4])  # left border
}


plot_expression_grid <- function(EL, meta, plot_fn = plot_sorted_bar) {
  row_labels <- colnames(EL)
  col_labels <- sort(unique(meta$Tissue))
  
  nrows <- length(row_labels)
  ncols <- length(col_labels)
  
  # Color palette (length 10 assumed)
  tab10 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
             "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
  
  par(mfrow = c(nrows, ncols),
      mar = c(0.5, 0.5, 0.2, 0.2),
      oma = c(6, 6, 0.5, 0.5))
  
  for (row_idx in seq_along(row_labels)) {
    row_x <- row_labels[row_idx]
    
    # Get data for all columns in this row to compute the shared ylim
    y_list <- lapply(col_labels, function(col_y) {
      tissue_samples <- meta[meta$Tissue == col_y, , drop = FALSE] %>% pull(Sample)
      y <- EL[tissue_samples, row_x]
      y
    })
    
    # Find max y value across this row
    row_ymax <- max(sapply(y_list, max, na.rm = TRUE), na.rm = TRUE)
    ylim <- c(0, row_ymax)
    
    for (col_idx in seq_along(col_labels)) {
      col_y <- col_labels[col_idx]
      y <- y_list[[col_idx]]
      
      # Whether to show y-axis ticks (only in 1st column)
      show_yticks <- (col_idx == 1)
      
      # Use corresponding color
      color <- tab10[(row_idx - 1) %% length(tab10) + 1]
      
      # Call user-defined plotting function with fixed ylim and color
      plot_fn(y, color = color, ylim = ylim, axes = show_yticks)
      
#       box()
      
      if (col_idx == 1) {
        mtext(row_x, side = 2, line = 1.5, las = 1, cex = 0.8)
      }
      
      if (row_idx == nrows) {
        mtext(col_y, side = 1, line = 1.5, las = 2, cex = 0.7)
      }
    }
  }
}
