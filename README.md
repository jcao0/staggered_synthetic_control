# Staggered Synthetic Control

This repository contains two related resources for staggered synthetic control.

## Replication Files

The folder [`replication package`](replication%20package/) contains the replication files for Cao, Lu, and Wu, "Synthetic Control Inference for Staggered Adoption."

The replication package has two main components:

- [`Simulation`](replication%20package/Simulation/): simulation files for the results in Section 3, including Figure 1.
- [`Intergovernmental coordination and criminality`](replication%20package/Intergovernmental%20coordination%20and%20criminality/): empirical replication files for Section 4, including Table 1 and Figures 2 and 3.

See [`replication package/README.md`](replication%20package/README.md) for software requirements, running times, data availability, and step-by-step replication instructions.

## R Package

The folder [`stagsynth`](stagsynth/) contains an R package for staggered synthetic control estimation and inference. The package is written by Zhanchao Fu.

To install the package locally from this repository, run:

```r
install.packages("stagsynth", repos = NULL, type = "source")
```

or from the terminal:

```sh
R CMD INSTALL stagsynth
```

The main function is 'ssc()', which estimates heterogeneous, event-time, and overall average treatment effects for staggered adoption designs, and constructs placebo-in-time confidence intervals and p-values. The package also includes helper functions for synthetic control weights (synthetic_control(), synthetic_control_batch()), reshaping panel data from long format (panel_to_matrices()), and checking the smallest eigenvalue of the SSC design matrix (ssc_min_eigenvalue()). S3 methods (print, summary, plot) are provided for "ssc" objects.
