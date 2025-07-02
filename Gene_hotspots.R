###############################################################################
### IDENTIFICATION OF HOTSPOTS OF GENES IN THE GENOME #########################
###############################################################################

########################
# Set up and data prep #
########################

# Set input files
gtfdir <- "~/Desktop/genespots/Inputs/dmel-all-r6.62.gtf"
clustersdir <- "~/Desktop/genespots/Inputs/defs_degs.txt"

# Load packages
library(GenomicRanges)
library(ggplot2)
library(KernSmooth)
library(dplyr)

# Read gtf file and transfrom into df
gtf <- rtracklayer::import(gtfdir)
genes_gr <- gtf[gtf$type == "gene"]
genes_gtf <- data.frame(
  gene_id = genes_gr$gene_id,
  seqnames = as.character(seqnames(genes_gr)),
  start = start(genes_gr),
  end = end(genes_gr),
  stringsAsFactors = FALSE
)

# Read the interested genes and combine with the gtf
df_clusters <- read.table(clustersdir, header = TRUE, sep = "\t", quote = "")
df_clusters <- df_clusters[, c(1,2,3)]   # Use to filter for cols containing: gene_id, gene_symbol and clusters. Remove or change colnums if needed
colnames(df_clusters) <- c("gene_id", "gene_symbol", "cluster")
merged_data <- merge(df_clusters, genes_gtf, by = "gene_id", all.x = TRUE)

# Transform chromosomes into factors
merged_data$seqnames <- factor(merged_data$seqnames)

# Define list of chromosomes and clusters to test
chromosomes <- unique(merged_data$seqnames)
#chromosomes <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")    # If interested only in particular chromosomes
clusters <- unique(merged_data$cluster)
#clusters <- c("upregulated", "downregulated")    # If interested only in up and down genes


########################
## Hotspot detection ###
########################

# Define bandwidth
bandwidth = 250000    # Sets up the size of the windows

# Define function to detect hotspots
analyze_hotspots <- function(data, chromosome, cluster_id, bw) {
  # Filter data
  chr_data <- data %>% 
    filter(seqnames == chromosome, cluster == cluster_id) %>%
    arrange(start)
  if (nrow(chr_data) < 10) return(NULL)
  # Get chr position range
  chr_range <- c(1, max(data$end[data$seqnames == chromosome], na.rm = TRUE))
  # kernel density calculation
  dens <- bkde(
    x = chr_data$start,
    bandwidth = bw,
    range.x = chr_range
  )
  # Identify hotspots
  d1 <- diff(dens$y)
  peaks <- which(d1[-length(d1)] > 0 & d1[-1] < 0) + 1
  if (length(peaks) == 0) return(NULL)
  # Filter hotspots
  min_density_threshold <- max(dens$y) * 0.90    # Keep only peaks with >90% of max
  peaks <- peaks[dens$y[peaks] > min_density_threshold]
  if (length(peaks) == 0) return(NULL)
  hotspot_positions <- dens$x[peaks]
  hotspot_densities <- dens$y[peaks]
  # Define gene ranges
  gene_ranges <- GRanges(
    seqnames = chr_data$seqnames,
    ranges = IRanges(start = chr_data$start, end = chr_data$end),
    gene_id = chr_data$gene_symbol
  )
  # Define hotspot ranges
  hotspot_ranges <- GRanges(
    seqnames = chromosome,
    ranges = IRanges(
      start = hotspot_positions - bw/2,
      end = hotspot_positions + bw/2
    )
  )
  # Define gene-hotspot overlaps
  overlaps <- findOverlaps(hotspot_ranges, gene_ranges)
  # Initialize results
  n_genes <- integer(length(hotspot_positions))
  genes <- character(length(hotspot_positions))
  # Loop to associate genes and hotspots
  for (i in seq_along(hotspot_positions)) {
    hits <- subjectHits(overlaps[queryHits(overlaps) == i])
    ids <- gene_ranges$gene_id[hits]
    n_genes[i] <- length(ids)
    genes[i] <- paste(ids, collapse = ", ")
  }
  # Get results
  results_df <- data.frame(
    chromosome = chromosome,
    cluster = cluster_id,
    position = hotspot_positions,
    density = hotspot_densities,
    n_genes = n_genes,
    genes = genes,
    stringsAsFactors = FALSE
  )
  # Remove hotspots with <3 genes
  results_df <- results_df[results_df$n_genes >= 3, ]    # Increase the number for stringency
  if (nrow(results_df) == 0) return(NULL)
  return(results_df)
}

# Test all combinations of chr and cluster
results <- list()
for (chr in chromosomes) {
  for (cl in clusters) {
    res <- analyze_hotspots(merged_data, chr, cl, bandwidth)
    if (!is.null(res)) results[[paste(chr, cl)]] <- res
  }
}

# Get results df
hotspots_df <- do.call(rbind, results)


########################
## Hotspot validation ##
########################

# Define function to validate identified hostpots
validate_hotspots <- function(hotspot_df, full_data, n_perm = 100, bw) {
  hotspot_df$p_value <- sapply(1:nrow(hotspot_df), function(i) {
    chr <- hotspot_df$chromosome[i]
    pos <- hotspot_df$position[i]
    cl  <- hotspot_df$cluster[i]
    # Obs density around hotspot
    window <- c(pos - bw, pos + bw)
    obs_data <- full_data %>% 
      filter(seqnames == chr, cluster == cl, start >= window[1], start <= window[2])
    obs_count <- nrow(obs_data)
    # Null density
    null_counts <- replicate(n_perm, {
      perm_data <- full_data %>%
        filter(seqnames == chr) %>%
        mutate(cluster_perm = sample(cluster))
      perm_cluster <- perm_data %>% 
        filter(cluster_perm == cl, start >= window[1], start <= window[2])
      nrow(perm_cluster)
    })
    # P-value
    sum(null_counts >= obs_count) / n_perm
  })
  hotspot_df$fdr <- p.adjust(hotspot_df$p_value, method = "BH")
  return(hotspot_df)
}

# Validate and sig hotspots
validated_hotspots <- validate_hotspots(hotspots_df, merged_data, n_perm = 100, bandwidth)
significant_hotspots <- validated_hotspots %>% filter(fdr < 0.05)


########################
## Density curves ######
########################

# Define function to get density curves
get_density_curve <- function(data, chromosome, cluster_id, bw) {
  chr_data <- data %>% 
    filter(seqnames == chromosome, cluster == cluster_id) %>%
    arrange(start)
  if (nrow(chr_data) < 1) return(NULL)
  # Get chrom range
  chr_range <- c(1, max(data$end[data$seqnames == chromosome], na.rm = TRUE))
  # Get density
  dens <- bkde(chr_data$start, bandwidth = bw, range.x = chr_range)
  data.frame(
    chromosome = chromosome,
    cluster = cluster_id,
    x = dens$x,
    y = dens$y
  )
}

# Get densities for all chr + cluster combinations
density_curves <- list()
for (chr in unique(merged_data$seqnames)) {
  for (cl in clusters) {
    curve <- get_density_curve(merged_data, chr, cl, bandwidth)
    if (!is.null(curve)) density_curves[[paste(chr, cl)]] <- curve
  }
}

# Get results
density_curve_df <- bind_rows(density_curves)


########################
## Null density curves #
########################

# Define function for null curve
get_null_density_curve <- function(data, chromosome, cluster_id, bw, n_perm = 100) {
  chr_data <- data %>% filter(seqnames == chromosome)
  chr_range <- c(1, max(chr_data$end, na.rm = TRUE))
  null_curves <- replicate(n_perm, {
    perm_data <- chr_data %>%
      mutate(cluster_perm = sample(cluster)) %>%
      filter(cluster_perm == cluster_id)
    if (nrow(perm_data) < 1) return(rep(NA, length(seq(chr_range[1], chr_range[2], length.out = 401))))
    bkde(perm_data$start, bandwidth = bw, range.x = chr_range)$y
  }, simplify = "array")
  # Mean densities
  mean_y <- rowMeans(null_curves, na.rm = TRUE)
  data.frame(
    chromosome = chromosome,
    cluster = cluster_id,
    x = seq(chr_range[1], chr_range[2], length.out = length(mean_y)),
    y = mean_y
  )
}

# Get null curves for each cluster and chr
null_density_curves <- list()
for (chr in unique(merged_data$seqnames)) {
  for (cl in clusters) {
    curve <- get_null_density_curve(merged_data, chr, cl, bandwidth)
    if (!is.null(curve)) null_density_curves[[paste(chr, cl)]] <- curve
  }
}
null_density_df <- bind_rows(null_density_curves)


########################
## Visualization #######
########################

## Filter chromosomes for significant hotspots only?
filter_chromosomes <- FALSE    # Change to TRUE to plot only chr with sig. hotspots

# Filtering
if (filter_chromosomes) {
  chromosome_filter <- unique(hotspots_df$chromosome)
  filtered_density_curve_df <- density_curve_df %>% 
    filter(chromosome %in% chromosome_filter)
  filtered_null_density_df <- null_density_df %>% 
    filter(chromosome %in% chromosome_filter)
  filtered_significant_hotspots <- significant_hotspots %>% 
    filter(chromosome %in% chromosome_filter)
} else {
  filtered_density_curve_df <- density_curve_df
  filtered_null_density_df <- null_density_df
  filtered_significant_hotspots <- significant_hotspots
}


# Plotting
ggplot() +
  # Null density
  geom_ribbon(
    data = filtered_null_density_df, 
    aes(x = x, ymin = 0, ymax = y), 
    fill = "gray70", alpha = 0.8
  ) +
  # Obs density
  geom_line(data = filtered_density_curve_df, aes(x = x, y = y), color = "gray30") +
  # Sig. hotspots
  geom_vline(data = filtered_significant_hotspots, 
             aes(xintercept = position), color = "red", linetype = "dashed") +
  geom_label(data = filtered_significant_hotspots,
             aes(x = position, y = density, label = paste0("FDR=", round(fdr, 3))),
             vjust = -0.5, size = 3) +
  facet_grid(cluster ~ chromosome, scales = "free_x") +
  labs(title = "Significant genomic hotspots with null densities",
       x = "Genomic position (Mb)",
       y = "Gene density") +
  #coord_cartesian(ylim = c(0, 4e-07)) +    # Set up to fix the Y axis
  scale_x_continuous(
    breaks = seq(0, max(filtered_null_density_df$x, na.rm = TRUE), by = 5000000),
    labels = function(x) x / 1e6
  )+
  theme_bw()

ggsave("Hostpots_plots.pdf", plot = last_plot())

write.table(significant_hotspots, file = "Sig_hotspots_table.txt", quote = FALSE, sep = "\t", row.names = FALSE, dec = )
