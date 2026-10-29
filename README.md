# FlowSOM-Based Exploratory Analysis of Murine FACS Data

Exploratory analysis of flow cytometry (FACS) data from single-cell suspensions. This repository documents a reproducible R workflow for FlowSOM-based unsupervised clustering, immune cell population exploration, dimensionality reduction, and visualization to characterize spinal cord immunophenotypes.

### Data preprocessing and compensation

FCS files were initially inspected and preprocessed in FlowJo. Events were exported following a gating strategy on live cells using the viable cell staining, according to the experimental gating strategy. Compensation settings were optimized in independent experiments, and the corresponding compensation matrix was embedded in the exported FCS files.

The R workflow imports the FCS files and applies compensation and downstream transformations before FlowSOM-based clustering and exploratory analysis.
