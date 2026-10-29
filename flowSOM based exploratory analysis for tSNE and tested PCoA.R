
# =============================================================================
# FlowSOM-based FACS immunophenotyping and bioinformatics analysis
# =============================================================================
#
# DESCRIPTION
#   Reproducible workflow for processing flow cytometry FCS files, performing
#   quality control, applying compensation and logicle transformation, training
#   a FlowSOM model, quantifying immune-cell populations, and visualising
#   high-dimensional marker-expression profiles.
#
# WORKFLOW
#   1. Configure paths, parameters, and experimental metadata.
#   2. Import and aggregate FCS files.
#   3. Apply compensation and logicle transformation.
#   4. Perform quality-control visualisation.
#   5. Train the FlowSOM self-organising map and identify metaclusters.
#   6. Quantify cluster abundance at the sample level.
#   7. Compare experimental groups using exploratory statistical analyses.
#   8. Generate PCoA and t-SNE visualisations.
#   9. Export tables, figures, the trained model, and session information.
#
# INPUT FILES
#   data/fcs/                    FCS files (not committed to GitHub)
#   metadata/sample_metadata.csv Sample-level experimental metadata
#
# METADATA FORMAT
#   File,Group
#   1,Vehicle
#   2,Treated
#
#   IMPORTANT: The File column must match the identifiers stored in the
#   aggregated flowFrame's "File" parameter. These identifiers may not be
#   identical to the original FCS filenames. Verify the mapping before
#   interpreting any group comparisons.
#
# IMPORTANT SCIENTIFIC CONSIDERATIONS
#   - Verify the compensation matrix and marker-channel identities.
#   - Confirm that the selected channels are fluorescence markers, not
#     scatter parameters, time, or file identifiers.
#   - Biological samples, not individual cells, are the units of inference
#     for sample-level group comparisons.
#   - PCoA and t-SNE are visualisation methods, not tests of significance.
#   - PERMANOVA can be affected by differences in group dispersion.
#   - Statistical models must reflect repeated measures, batches, and other
#     relevant experimental-design constraints.
#
# REPRODUCIBILITY
#   The random seed, FlowSOM parameters, and R session information are saved
#   or documented to facilitate reproducibility.
#
# =============================================================================


# =============================================================================
# 1. PACKAGES
# =============================================================================

# Install any missing CRAN packages before running this script.
# FlowSOM and flowCore are available through Bioconductor.

cran_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "Rtsne",
  "patchwork",
  "vegan",
  "ape",
  "pheatmap",
  "viridis"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("flowCore", quietly = TRUE)) {
  BiocManager::install("flowCore")
}

if (!requireNamespace("FlowSOM", quietly = TRUE)) {
  BiocManager::install("FlowSOM")
}

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_cran) > 0L) {
  install.packages(missing_cran)
}

# Load packages explicitly used in the workflow.
library(flowCore)
library(FlowSOM)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(Rtsne)
library(patchwork)
library(vegan)
library(ape)
library(pheatmap)
library(viridis)


# =============================================================================
# 2. CONFIGURATION
# =============================================================================

# Use paths relative to the GitHub repository root.
# Do not use personal computer paths or setwd() in the published script.

data_dir <- file.path("data", "fcs")
metadata_file <- file.path("metadata", "sample_metadata.csv")
results_dir <- "results"

# Create output directories if they do not already exist.
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(
  file.path(results_dir, "figures"),
  showWarnings = FALSE,
  recursive = TRUE
)
dir.create(
  file.path(results_dir, "tables"),
  showWarnings = FALSE,
  recursive = TRUE
)
dir.create(
  file.path(results_dir, "models"),
  showWarnings = FALSE,
  recursive = TRUE
)

# Set random seeds for reproducible sampling and dimensionality reduction.
random_seed <- 2026L
set.seed(random_seed)

# FlowSOM parameters.
som_xdim <- 5L
som_ydim <- 5L
n_metaclusters <- 8L

# Event aggregation target.
# This specifies the overall target passed to AggregateFlowFrames().
# It does not guarantee equal numbers of events per sample.
events_target_per_file <- 2800L

# Maximum number of events sampled per sample for t-SNE.
# Adjust according to available computational resources.
cells_per_sample <- 10000L

# t-SNE parameters.
tsne_perplexity <- 30
tsne_theta <- 0.5

# Marker selection.
# IMPORTANT:
# Replace this example with verified fluorescence-marker channel names
# whenever possible. These must match the column names in the FCS data.
#
# Example:
# marker_channels <- c("CD45", "CD3", "CD4", "CD8", "CD11b")
#
# If NULL, the script uses the numeric channel indices below.
marker_channels <- NULL
marker_indices <- 7:23

# Group variable used for group-level statistical comparisons.
group_variable <- "Group"

# Number of permutations for PERMANOVA and dispersion testing.
n_permutations <- 999L


# =============================================================================
# 3. INPUT VALIDATION AND EXPERIMENTAL METADATA
# =============================================================================

# Identify FCS files in the input directory.
fcs_files <- list.files(
  path = data_dir,
  pattern = "\\.fcs$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(fcs_files) == 0L) {
  stop(
    "No FCS files found. Check the data_dir configuration: ",
    normalizePath(data_dir, mustWork = FALSE)
  )
}

message("Number of FCS files detected: ", length(fcs_files))

# Save the input filenames for traceability.
writeLines(
  basename(fcs_files),
  con = file.path(results_dir, "input_files.txt")
)

# Import sample-level metadata.
if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Ensure the metadata contains the fields required for group comparisons.
required_columns <- c("File", group_variable)

if (!all(required_columns %in% colnames(metadata))) {
  stop(
    "Metadata must contain the following columns: ",
    paste(required_columns, collapse = ", ")
  )
}

metadata$File <- as.character(metadata$File)
metadata[[group_variable]] <- as.character(metadata[[group_variable]])

if (anyDuplicated(metadata$File)) {
  stop("The metadata contains duplicated sample identifiers.")
}

if (anyNA(metadata[[group_variable]]) ||
    any(metadata[[group_variable]] == "")) {
  stop("The experimental group variable contains missing values.")
}

# Additional metadata columns such as MouseID, Status, Timepoint, or Batch
# can be included in the CSV and used in downstream analyses.


# =============================================================================
# 4. IMPORT AND AGGREGATE FCS FILES
# =============================================================================

# Aggregate events from the input files.
# Aggregation facilitates a common FlowSOM model across samples.
#
# Confirm that this sampling strategy is suitable for the experiment.
# For strongly unequal sample sizes, consider explicit per-file sampling
# before aggregation instead of relying only on a global aggregation target.

ff <- AggregateFlowFrames(
  fcs_files,
  cTotal = events_target_per_file * length(fcs_files)
)

exprs_raw <- exprs(ff)

if (nrow(exprs_raw) == 0L) {
  stop("No events were retained after FCS aggregation.")
}

# The File parameter is used to trace events back to their source samples.
if (!"File" %in% colnames(exprs_raw)) {
  stop(
    "The aggregated flowFrame does not contain a 'File' parameter. ",
    "Inspect the aggregation output and establish sample provenance."
  )
}

# Save the observed file identifiers.
observed_file_ids <- unique(as.character(exprs_raw[, "File"]))

writeLines(
  observed_file_ids,
  con = file.path(results_dir, "observed_file_ids.txt")
)

# Count the events contributed by each source identifier.
events_per_file <- as.data.frame(
  table(File = as.character(exprs_raw[, "File"])),
  stringsAsFactors = FALSE
)

colnames(events_per_file) <- c("File", "n_events")

write.csv(
  events_per_file,
  file.path(results_dir, "tables", "events_per_file.csv"),
  row.names = FALSE
)

# Visualise event counts across source files.
p_events <- ggplot(
  events_per_file,
  aes(x = File, y = n_events)
) +
  geom_col(fill = "steelblue") +
  theme_bw() +
  labs(
    title = "Events retained per source file",
    x = "Source file identifier",
    y = "Number of events"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(
  file.path(results_dir, "figures", "events_per_file.png"),
  plot = p_events,
  width = 10,
  height = 5,
  dpi = 300
)


# =============================================================================
# 5. COMPENSATION AND LOGICLE TRANSFORMATION
# =============================================================================

# Retrieve the spillover matrix embedded in the aggregated flowFrame.
# Confirm that this matrix is appropriate for the panel and acquisition setup.
spillover <- ff@description[["$SPILLOVER"]]

if (is.null(spillover)) {
  stop(
    "No $SPILLOVER matrix was found. ",
    "Supply and validate the appropriate compensation matrix."
  )
}

# Apply fluorescence compensation.
ff_comp <- compensate(ff, spillover)

# Identify the channels to use for transformation and clustering.
all_channels <- colnames(exprs(ff_comp))

if (is.null(marker_channels)) {
  
  # Numeric selection is retained for compatibility with the original script.
  # Verify these indices against the actual panel before running the analysis.
  if (max(marker_indices) > length(all_channels) ||
      min(marker_indices) < 1L) {
    stop("The selected marker indices are outside the available channels.")
  }
  
  selected_channels <- all_channels[marker_indices]
  
} else {
  
  # Prefer explicit marker names for a reusable GitHub workflow.
  missing_channels <- setdiff(marker_channels, all_channels)
  
  if (length(missing_channels) > 0L) {
    stop(
      "The following marker channels are missing: ",
      paste(missing_channels, collapse = ", ")
    )
  }
  
  selected_channels <- marker_channels
}

if (anyDuplicated(selected_channels)) {
  stop("The selected marker channels contain duplicates.")
}

message("Channels selected for analysis:")
print(selected_channels)

# Estimate a logicle transformation from compensated fluorescence channels.
# Logicle transformation is commonly used for flow cytometry data because
# it accommodates positive and negative compensated fluorescence values.
#
# Do not apply this transformation blindly to scatter, time, or file channels.

transformation <- estimateLogicle(
  ff_comp,
  channels = selected_channels
)

ff_transformed <- transform(ff_comp, transformation)

exprs_transformed <- exprs(ff_transformed)

if (!all(selected_channels %in% colnames(exprs_transformed))) {
  stop("One or more selected channels are missing after transformation.")
}

if (any(!is.finite(exprs_transformed[, selected_channels, drop = FALSE]))) {
  stop("Non-finite values detected in transformed marker expression.")
}

message("Total events available after transformation: ", nrow(exprs_transformed))


# =============================================================================
# 6. QUALITY CONTROL: MARKER-EXPRESSION DISTRIBUTIONS
# =============================================================================

# Inspect marker-expression distributions after compensation and transformation.
# These plots can reveal unusual distributions, unstable channels, and
# potential differences between samples.

qc_df <- as.data.frame(
  exprs_transformed[, selected_channels, drop = FALSE]
)

qc_long <- tidyr::pivot_longer(
  tibble::as_tibble(qc_df),
  cols = everything(),
  names_to = "Marker",
  values_to = "Expression"
)

p_qc <- ggplot(
  qc_long,
  aes(x = Expression)
) +
  geom_histogram(
    bins = 80,
    fill = "steelblue",
    colour = NA
  ) +
  facet_wrap(
    ~ Marker,
    scales = "free",
    ncol = 4
  ) +
  theme_bw() +
  labs(
    title = "Transformed marker-expression distributions",
    x = "Transformed expression",
    y = "Number of events"
  )

ggsave(
  file.path(results_dir, "figures", "marker_expression_qc.png"),
  plot = p_qc,
  width = 14,
  height = 10,
  dpi = 300
)

# Compare marker-expression distributions between source files.
# This is a descriptive QC plot, not a formal batch-effect test.

qc_by_file <- as.data.frame(
  exprs_transformed[, selected_channels, drop = FALSE]
)

qc_by_file$File <- as.character(exprs_transformed[, "File"])

qc_by_file_long <- tidyr::pivot_longer(
  qc_by_file,
  cols = all_of(selected_channels),
  names_to = "Marker",
  values_to = "Expression"
)

p_qc_file <- ggplot(
  qc_by_file_long,
  aes(x = File, y = Expression)
) +
  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.2
  ) +
  facet_wrap(
    ~ Marker,
    scales = "free_y",
    ncol = 4
  ) +
  theme_bw() +
  labs(
    title = "Marker expression across source files",
    x = "Source file identifier",
    y = "Transformed expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(
  file.path(results_dir, "figures", "marker_expression_by_file.png"),
  plot = p_qc_file,
  width = 14,
  height = 10,
  dpi = 300
)


# =============================================================================
# 7. FLOW SOM CLUSTERING
# =============================================================================

# FlowSOM combines a self-organising map (SOM) with metaclustering to
# identify groups of events with similar multivariate marker-expression profiles.
#
# xdim and ydim control the dimensions of the SOM grid.
# nClus controls the number of metaclusters.
#
# These parameters should be selected based on the marker panel, data structure,
# biological interpretability, and stability of the resulting populations.

# Convert selected channel names to their corresponding column indices.
som_marker_indices <- match(
  selected_channels,
  colnames(exprs_transformed)
)

if (anyNA(som_marker_indices)) {
  stop("Unable to resolve marker channel indices for FlowSOM.")
}

set.seed(random_seed)

fsom <- FlowSOM(
  input = ff_transformed,
  colsToUse = som_marker_indices,
  seed = random_seed,
  xdim = som_xdim,
  ydim = som_ydim,
  nClus = n_metaclusters
)

# Remove appended annotation text from display labels, where present.
if (!is.null(fsom$prettyColnames)) {
  fsom$prettyColnames <- gsub(
    " <.*",
    "",
    fsom$prettyColnames
  )
}

# Save the trained model for reproducibility and downstream analysis.
saveRDS(
  fsom,
  file.path(results_dir, "models", "flowsom_model.rds")
)

# Generate the standard FlowSOM summary.
# Verify the function signature if using a different FlowSOM package version.
if ("FlowSOMmary" %in% getNamespaceExports("FlowSOM")) {
  FlowSOM::FlowSOMmary(
    fsom = fsom,
    plotFile = file.path(
      results_dir,
      "figures",
      "FlowSOM_summary.pdf"
    )
  )
} else {
  message(
    "FlowSOMmary() is not exported by the installed FlowSOM version. ",
    "Generate the summary using the plotting functions available in that version."
  )
}


# =============================================================================
# 8. ASSIGN EVENTS TO SOM NODES AND METACLUSTERS
# =============================================================================

# Extract the SOM-node assignment for each event.
# The metaclustering vector maps SOM nodes to broader cell populations.

if (is.null(fsom$map$mapping) ||
    is.null(fsom$metaclustering)) {
  stop("The FlowSOM object does not contain the expected mapping information.")
}

node_id <- fsom$map$mapping[, 1]

if (length(node_id) != nrow(exprs_transformed)) {
  stop(
    "The number of FlowSOM event assignments does not match the input events. ",
    "Check filtering, sampling, and event order."
  )
}

metacluster_id <- fsom$metaclustering[node_id]

if (anyNA(metacluster_id)) {
  stop("Some events could not be assigned to a metacluster.")
}

# Extract the original sample identifier for each event.
event_file_ids <- as.character(exprs_transformed[, "File"])

# Build a cell-level table containing marker expression, sample ID,
# SOM-node assignment, and metacluster assignment.
cell_df <- as.data.frame(
  exprs_transformed[, selected_channels, drop = FALSE]
)

cell_df$File <- event_file_ids
cell_df$SOM_node <- node_id
cell_df$Cluster <- factor(metacluster_id)

# Match sample-level metadata to event-level data.
# This requires metadata$File to match the actual flowFrame File identifiers.
unmatched_ids <- setdiff(unique(cell_df$File), metadata$File)

if (length(unmatched_ids) > 0L) {
  stop(
    "Some flowFrame file identifiers do not match metadata$File: ",
    paste(unmatched_ids, collapse = ", "),
    ". Inspect observed_file_ids.txt and correct the sample mapping."
  )
}

# Attach experimental metadata to every event.
cell_df <- cell_df %>%
  left_join(
    metadata,
    by = "File"
  )

if (anyNA(cell_df[[group_variable]])) {
  stop("Some events have missing experimental group assignments.")
}

# Export a cell-level cluster assignment table.
# This may be large; consider exporting only a sample or aggregate table
# if the repository does not need event-level data.
write.csv(
  cell_df %>%
    select(File, all_of(group_variable), SOM_node, Cluster),
  file.path(results_dir, "tables", "cell_cluster_assignments.csv"),
  row.names = FALSE
)


# =============================================================================
# 9. CLUSTER ABUNDANCE BY SAMPLE
# =============================================================================

# Quantify the number of cells assigned to each metacluster in each sample.
# Biological samples, rather than individual cells, are used for subsequent
# group-level statistical comparisons.

cluster_counts <- cell_df %>%
  count(
    File,
    .data[[group_variable]],
    Cluster,
    name = "n_cells"
  )

# Calculate the percentage of events assigned to each cluster within each
# sample. Percentages sum to 100% across the observed clusters for each sample.
cluster_percentages <- cluster_counts %>%
  group_by(File) %>%
  mutate(
    percentage = 100 * n_cells / sum(n_cells)
  ) %>%
  ungroup()

write.csv(
  cluster_counts,
  file.path(results_dir, "tables", "cluster_counts.csv"),
  row.names = FALSE
)

write.csv(
  cluster_percentages,
  file.path(results_dir, "tables", "cluster_percentages.csv"),
  row.names = FALSE
)

# Plot cluster abundance across samples.
p_cluster_abundance <- ggplot(
  cluster_percentages,
  aes(
    x = File,
    y = percentage,
    fill = Cluster
  )
) +
  geom_col() +
  theme_bw() +
  labs(
    title = "FlowSOM metacluster composition by sample",
    x = "Sample",
    y = "Cells per metacluster (%)",
    fill = "Metacluster"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(
  file.path(results_dir, "figures", "cluster_abundance_by_sample.png"),
  plot = p_cluster_abundance,
  width = 12,
  height = 6,
  dpi = 300
)

# Plot the abundance of each metacluster by experimental group.
p_cluster_group <- ggplot(
  cluster_percentages,
  aes(
    x = .data[[group_variable]],
    y = percentage,
    fill = .data[[group_variable]]
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_point(
    position = position_jitter(width = 0.12),
    alpha = 0.8,
    size = 2
  ) +
  facet_wrap(
    ~ Cluster,
    scales = "free_y"
  ) +
  theme_bw() +
  labs(
    title = "Metacluster abundance across experimental groups",
    x = "Experimental group",
    y = "Cells per metacluster (%)"
  ) +
  guides(fill = "none")

ggsave(
  file.path(results_dir, "figures", "cluster_abundance_by_group.png"),
  plot = p_cluster_group,
  width = 12,
  height = 8,
  dpi = 300
)


# =============================================================================
# 10. CLUSTER-WISE GROUP COMPARISONS
# =============================================================================

# This section performs exploratory cluster-wise comparisons using biological
# samples as independent observations.
#
# For exactly two independent groups, a Wilcoxon rank-sum test is used.
# P-values are adjusted using the Benjamini-Hochberg procedure.
#
# IMPORTANT:
#   - This is not a substitute for a model that accounts for repeated measures,
#     paired samples, covariates, or batch effects.
#   - Small sample sizes can provide limited statistical power.
#   - Cluster percentages are compositional and are not independent of one
#     another because their sum within each sample is fixed at 100%.
#   - Consider a design-appropriate compositional or count-based model when
#     the study requires formal differential-abundance inference.

group_levels <- unique(
  as.character(metadata[[group_variable]])
)

if (length(group_levels) == 2L) {
  
  cluster_stats <- cluster_percentages %>%
    group_by(Cluster) %>%
    summarise(
      n_group_1 = sum(
        .data[[group_variable]] == group_levels[1]
      ),
      n_group_2 = sum(
        .data[[group_variable]] == group_levels[2]
      ),
      median_group_1 = median(
        percentage[
          .data[[group_variable]] == group_levels[1]
        ],
        na.rm = TRUE
      ),
      median_group_2 = median(
        percentage[
          .data[[group_variable]] == group_levels[2]
        ],
        na.rm = TRUE
      ),
      p_value = {
        values_1 <- percentage[
          .data[[group_variable]] == group_levels[1]
        ]
        values_2 <- percentage[
          .data[[group_variable]] == group_levels[2]
        ]
        
        if (length(values_1) > 0L &&
            length(values_2) > 0L) {
          wilcox.test(
            values_1,
            values_2,
            exact = FALSE
          )$p.value
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    ) %>%
    mutate(
      group_1 = group_levels[1],
      group_2 = group_levels[2],
      median_difference = median_group_2 - median_group_1,
      adjusted_p_value = p.adjust(
        p_value,
        method = "BH"
      )
    )
  
  write.csv(
    cluster_stats,
    file.path(results_dir, "tables", "cluster_group_statistics.csv"),
    row.names = FALSE
  )
  
} else {
  
  message(
    "Cluster-wise Wilcoxon comparisons were skipped: ",
    "this section requires exactly two experimental groups. ",
    "Use a suitable multi-group model for the study design."
  )
}


# =============================================================================
# 11. MULTIVARIATE ANALYSIS OF SAMPLE-LEVEL CLUSTER ABUNDANCES
# =============================================================================

# Create a sample-by-cluster matrix of relative abundances.
# Each row represents one sample, and each column represents one metacluster.

cluster_wide <- cluster_percentages %>%
  select(
    File,
    all_of(group_variable),
    Cluster,
    percentage
  ) %>%
  pivot_wider(
    names_from = Cluster,
    values_from = percentage,
    values_fill = 0,
    names_prefix = "Cluster_"
  )

# Store one row of metadata per sample.
sample_metadata <- cluster_wide %>%
  select(
    File,
    all_of(group_variable)
  ) %>%
  distinct()

if (anyDuplicated(sample_metadata$File)) {
  stop("Sample metadata contains duplicated sample identifiers.")
}

# Build a numerical sample-by-cluster abundance matrix.
abundance_matrix <- cluster_wide %>%
  select(starts_with("Cluster_")) %>%
  as.data.frame()

rownames(abundance_matrix) <- cluster_wide$File

# Check sample order and numerical validity.
if (!identical(
  rownames(abundance_matrix),
  sample_metadata$File
)) {
  stop("Sample order differs between the abundance matrix and metadata.")
}

if (any(!is.finite(as.matrix(abundance_matrix)))) {
  stop("The cluster-abundance matrix contains non-finite values.")
}

# Bray-Curtis dissimilarity is commonly used for comparing community-like
# abundance profiles. Here it describes differences in metacluster composition.
cluster_distance <- vegan::vegdist(
  abundance_matrix,
  method = "bray"
)


# =============================================================================
# 12. PERMANOVA
# =============================================================================

# PERMANOVA tests whether sample-level cluster-abundance profiles differ
# between groups. Permutations must be compatible with the experimental design.
#
# The model below assumes independent samples and a single grouping variable.
# Modify the model and permutation scheme for paired, repeated-measures,
# blocked, or longitudinal experiments.

set.seed(random_seed)

permanova_result <- vegan::adonis2(
  cluster_distance ~ .data[[group_variable]],
  data = sample_metadata,
  permutations = n_permutations
)

print(permanova_result)

write.csv(
  as.data.frame(permanova_result),
  file.path(results_dir, "tables", "permanova_results.csv")
)


# =============================================================================
# 13. MULTIVARIATE DISPERSION TEST
# =============================================================================

# PERMANOVA can be influenced by differences in within-group dispersion.
# betadisper() assesses whether groups differ in their distance to their
# respective multivariate centroid.

dispersion_result <- vegan::betadisper(
  cluster_distance,
  group = factor(sample_metadata[[group_variable]])
)

set.seed(random_seed)

dispersion_test <- vegan::permutest(
  dispersion_result,
  permutations = n_permutations
)

print(dispersion_test)

writeLines(
  capture.output(dispersion_test),
  file.path(results_dir, "tables", "dispersion_test.txt")
)

# Export sample-level distances to group centroids for further inspection.
dispersion_table <- data.frame(
  File = names(dispersion_result$distances),
  distance_to_centroid = as.numeric(dispersion_result$distances)
)

write.csv(
  dispersion_table,
  file.path(results_dir, "tables", "dispersion_distances.csv"),
  row.names = FALSE
)


# =============================================================================
# 14. PRINCIPAL COORDINATES ANALYSIS (PCoA)
# =============================================================================

# PCoA visualises sample-level dissimilarities calculated from the
# cluster-abundance profiles.
#
# Axis percentages are derived from the eigenvalues. For distance matrices
# that yield negative eigenvalues, interpretation of these percentages
# requires additional care.

pcoa_result <- ape::pcoa(cluster_distance)

if (ncol(pcoa_result$vectors) < 2L) {
  stop("PCoA did not return enough axes for a two-dimensional plot.")
}

pcoa_scores <- as.data.frame(
  pcoa_result$vectors[, 1:2, drop = FALSE]
)

colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
pcoa_scores$File <- rownames(pcoa_scores)

pcoa_scores <- pcoa_scores %>%
  left_join(
    sample_metadata,
    by = "File"
  )

relative_eigenvalues <- pcoa_result$values$Relative_eig

axis1_pct <- 100 * relative_eigenvalues[1]
axis2_pct <- 100 * relative_eigenvalues[2]

# Plot samples coloured by experimental group.
p_pcoa <- ggplot(
  pcoa_scores,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = .data[[group_variable]]
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  coord_equal() +
  theme_bw() +
  labs(
    title = "PCoA of FlowSOM metacluster abundances",
    x = sprintf("PCoA1 (%.1f%%)", axis1_pct),
    y = sprintf("PCoA2 (%.1f%%)", axis2_pct),
    colour = "Experimental group"
  )

ggsave(
  file.path(results_dir, "figures", "pcoa_cluster_abundance.pdf"),
  plot = p_pcoa,
  width = 8,
  height = 6
)

ggsave(
  file.path(results_dir, "figures", "pcoa_cluster_abundance.png"),
  plot = p_pcoa,
  width = 8,
  height = 6,
  dpi = 300
)

write.csv(
  pcoa_scores,
  file.path(results_dir, "tables", "pcoa_coordinates.csv"),
  row.names = FALSE
)


# =============================================================================
# 15. CLUSTER-LEVEL MARKER-EXPRESSION HEATMAP
# =============================================================================

# Summarise marker expression by FlowSOM metacluster.
# Median expression is used to provide a robust descriptive summary.
#
# This helps assess whether clusters have distinct marker-expression profiles.
# Biological interpretation should be based on the marker panel and controls.

cluster_marker_medians <- cell_df %>%
  group_by(Cluster) %>%
  summarise(
    across(
      all_of(selected_channels),
      ~ median(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

# Convert to a matrix with clusters as rows and markers as columns.
heatmap_matrix <- as.data.frame(cluster_marker_medians)
rownames(heatmap_matrix) <- as.character(heatmap_matrix$Cluster)

heatmap_matrix$Cluster <- NULL
heatmap_matrix <- as.matrix(heatmap_matrix)

# Scale markers across clusters to highlight relative patterns.
# Constant markers cannot be scaled and are excluded from the scaled heatmap.
marker_sd <- apply(
  heatmap_matrix,
  2,
  sd,
  na.rm = TRUE
)

heatmap_matrix <- heatmap_matrix[
  ,
  is.finite(marker_sd) & marker_sd > 0,
  drop = FALSE
]

if (ncol(heatmap_matrix) > 0L) {
  
  heatmap_scaled <- scale(heatmap_matrix)
  
  pdf(
    file.path(results_dir, "figures", "cluster_marker_heatmap.pdf"),
    width = 10,
    height = 7
  )
  
  pheatmap::pheatmap(
    heatmap_scaled,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "Marker-expression profiles by FlowSOM metacluster",
    border_color = NA
  )
  
  dev.off()
  
  write.csv(
    heatmap_matrix,
    file.path(results_dir, "tables", "cluster_marker_medians.csv")
  )
  
} else {
  
  message(
    "Heatmap skipped because no markers with non-zero variance were available."
  )
}


# =============================================================================
# 16. t-SNE VISUALISATION OF FLOW SOM EVENTS
# =============================================================================

# t-SNE provides an exploratory two-dimensional visualisation of high-dimensional
# marker-expression data. It does not directly test statistical differences
# between groups or establish biological separation.
#
# Events are sampled separately within each source sample to reduce the
# influence of samples with disproportionately large event counts.

set.seed(random_seed)

tsne_cells <- cell_df %>%
  group_by(File) %>%
  slice_sample(
    n = min(cells_per_sample, n()),
    replace = FALSE
  ) %>%
  ungroup()

# Build a numerical event-by-marker matrix.
tsne_input <- as.matrix(
  tsne_cells[, selected_channels, drop = FALSE]
)

if (any(!is.finite(tsne_input))) {
  stop("The t-SNE input contains non-finite values.")
}

# Rtsne requires perplexity to be smaller than (number of cells - 1) / 3.
if (nrow(tsne_input) <= 3 * tsne_perplexity + 1) {
  stop(
    "Too few sampled events for the chosen t-SNE perplexity. ",
    "Reduce tsne_perplexity or increase the number of sampled events."
  )
}

# Run t-SNE with a fixed seed.
set.seed(random_seed)

tsne_result <- Rtsne(
  tsne_input,
  dims = 2,
  perplexity = tsne_perplexity,
  theta = tsne_theta,
  check_duplicates = FALSE,
  pca = TRUE,
  verbose = TRUE
)

# Add the t-SNE coordinates to the sampled cell table.
tsne_cells$tSNE1 <- tsne_result$Y[, 1]
tsne_cells$tSNE2 <- tsne_result$Y[, 2]

# Visualise metaclusters, faceted by experimental group.
p_tsne_cluster <- ggplot(
  tsne_cells,
  aes(
    x = tSNE1,
    y = tSNE2,
    colour = Cluster
  )
) +
  geom_point(
    size = 0.25,
    alpha = 0.8
  ) +
  facet_wrap(
    vars(.data[[group_variable]])
  ) +
  coord_equal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(
      fill = "grey90",
      colour = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.position = "right"
  ) +
  labs(
    title = "t-SNE of FlowSOM-assigned events",
    subtitle = "Events sampled within each source sample",
    x = "t-SNE 1",
    y = "t-SNE 2",
    colour = "Metacluster"
  )

ggsave(
  file.path(results_dir, "figures", "tsne_by_group.pdf"),
  plot = p_tsne_cluster,
  width = 12,
  height = 7
)

ggsave(
  file.path(results_dir, "figures", "tsne_by_group.tiff"),
  plot = p_tsne_cluster,
  width = 12,
  height = 7,
  dpi = 600,
  compression = "lzw"
)


# =============================================================================
# 17. t-SNE COLOURED BY INDIVIDUAL MARKER EXPRESSION
# =============================================================================

# Generate a separate t-SNE plot for each marker.
# The .data pronoun ensures that ggplot maps the actual values of the selected
# marker column, rather than interpreting the marker name as a literal string.

marker_plots <- lapply(
  selected_channels,
  function(marker_of_interest) {
    
    ggplot(
      tsne_cells,
      aes(
        x = tSNE1,
        y = tSNE2,
        colour = .data[[marker_of_interest]]
      )
    ) +
      geom_point(
        size = 0.25,
        alpha = 0.8
      ) +
      scale_colour_viridis_c() +
      coord_equal() +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        legend.position = "right"
      ) +
      labs(
        title = paste("Marker expression:", marker_of_interest),
        x = "t-SNE 1",
        y = "t-SNE 2",
        colour = marker_of_interest
      )
  }
)

names(marker_plots) <- selected_channels

# Combine marker plots into a multi-panel figure.
combined_marker_plot <- patchwork::wrap_plots(
  marker_plots,
  ncol = 3
)

ggsave(
  file.path(results_dir, "figures", "tsne_marker_expression.tiff"),
  plot = combined_marker_plot,
  width = 18,
  height = 18,
  dpi = 600,
  compression = "lzw"
)

# Export the sampled cell coordinates and marker expression for reproducibility.
write.csv(
  tsne_cells,
  file.path(results_dir, "tables", "tsne_cell_coordinates.csv"),
  row.names = FALSE
)


# =============================================================================
# 18. SAVE THE ANALYSIS CONFIGURATION AND SESSION INFORMATION
# =============================================================================

# Save key parameters used in the analysis.
analysis_configuration <- list(
  random_seed = random_seed,
  events_target_per_file = events_target_per_file,
  som_xdim = som_xdim,
  som_ydim = som_ydim,
  n_metaclusters = n_metaclusters,
  selected_channels = selected_channels,
  cells_per_sample = cells_per_sample,
  tsne_perplexity = tsne_perplexity,
  tsne_theta = tsne_theta,
  group_variable = group_variable,
  n_permutations = n_permutations
)

saveRDS(
  analysis_configuration,
  file.path(results_dir, "models", "analysis_configuration.rds")
)

# Record the R version, operating system, and package versions.
writeLines(
  capture.output(sessionInfo()),
  con = file.path(results_dir, "sessionInfo.txt")
)

# Print the final output directory.
message("Analysis workflow completed.")
message("Results directory: ", normalizePath(results_dir))

# =============================================================================
# END OF SCRIPT
# =============================================================================
