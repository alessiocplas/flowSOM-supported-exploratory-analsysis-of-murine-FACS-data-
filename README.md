# FlowSOM-Based Exploratory Analysis of Murine FACS Data

Exploratory analysis of flow cytometry (FACS) data from single-cell suspensions. This repository documents a reproducible R workflow for FlowSOM-based unsupervised clustering, immune cell population exploration, dimensionality reduction, and visualization to characterize tissue isolated immunophenotypes.

### Data preprocessing and compensation

FCS files were initially inspected and preprocessed in [FlowJo](https://www.flowjo.com/). Events were exported following a gating strategy on live cells using the viable cell staining, according to the experimental gating strategy. Compensation settings were optimized in independent experiments, and the corresponding compensation matrix was embedded in the exported FCS files.

The R workflow imports the FCS files and applies compensation and downstream transformations before FlowSOM-based clustering and exploratory analysis.

**Analysis workflow**

The pipeline includes:

- **FACS data preprocessing**: Import and aggregation of FCS files, using compensation matrices embedded in the original files.

- **Compensation and transformation**: Application of compensation and logicle transformation to prepare fluorescence intensity data for downstream analysis.

- **Quality control**: Cells were previously assessed in FlowJo, and cells with low Live staining were exported for downstream computational analysis. The compensation matrix was optimized in independent experiments before acquisition and is included in the FCS files.

- **Unsupervised clustering**: Application of FlowSOM self-organizing maps (SOMs) to identify distinct immune cell populations based on marker expression profiles.

- **Cluster characterization**: Visualization of marker expression patterns and cluster-level phenotypic profiles using heatmaps.

- **Dimensionality reduction**: t-distributed stochastic neighbor embedding (t-SNE) for exploratory visualization of cellular heterogeneity and population structure.

- **Population-level analysis**: Assessment of cluster abundance and distribution across experimental samples, complemented by multivariate analyses, including principal coordinates analysis (PCoA), PERMANOVA, and tests for multivariate dispersion.

**Objectives**

The main objective is to provide an exploratory and reproducible computational framework for investigating immune cell heterogeneity in the central nervous system during EAE. By combining unsupervised clustering with dimensionality reduction and statistical analysis, this workflow facilitates the identification and characterization of immune cell populations potentially involved in neuroinflammatory processes.

Dependencies

- R

- flowCore

- FlowSOM

- ggplot2

- pheatmap

- Rtsne

- dplyr

- phyloseq-compatible data handling where applicable

The scripts and associated outputs are intended to facilitate reproducibility, transparent data exploration, and further development of computational approaches for high-dimensional flow cytometry analysis.
