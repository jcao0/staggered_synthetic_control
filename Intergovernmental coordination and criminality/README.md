## Application study: Intergovernmental Coordination and Criminality

This repository contains replication files for the application studies in Cao, Lu, and Wu (2020), “Synthetic Control Inference for Staggered Adoption,” . It replicates the empirical analysis in Alcocer M (2025), “Increasing Intergovernmental Coordination to Fight Crime: Evidence from Mexico.”

It includes the data and code needed to reproduce the main tables and figures, specifically:

-   Table 1: Smallest eigenvalues of the sample analogue of the design matrix.

-   Figure 2: Treatment scheme for different outcomes.

-   Figure 3: Treatment effects estimated using staggered synthetic control (SSC) and generalized synthetic control (GSC) methods.

### Data

-   The analysis relies on municipal-level crime and cartel data originally provided by the authors of "Increasing Intergovernmental Coordination to Fight Crime: Evidence from Mexico" and hosted at the Dataverse archive ([doi:10.7910/DVN/PBTTDM](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi%3A10.7910%2FDVN%2FPBTTDM)).

### Directory Structure

-   [`raw_data/`](raw_data) — Contains the raw data for all analyses.

    -   [`psrm_cartel_data.RData`](raw_data/psrm_cartel_data.RData) — Annual cartel-related data, obtained from the Dataverse replication archive (doi:10.7910/DVN/PBTTDM).

    -   [`psrm_crime_data.RData`](raw_data/psrm_crime_data.RData) — Monthly crime data, obtained from the same Dataverse replication archive.

-   [`data_cleaning.R`](data_cleaning.R) — R script that processes the raw data in `raw_data/` and saves cleaned data to analysis folders.

-   [`smallest_eigenvalues/`](smallest_eigenvalues) — Contains data and code to compute the smallest eigenvalues of the sample analogue of the design matrices and export Table 1.

    -   [`cleaned_data/`](smallest_eigenvalues/cleaned_data) — Cleaned outcome data produced by [`data_cleaning.R`](data_cleaning.R).

        -   [`psrm_cartel_data.csv`](smallest_eigenvalues/cleaned_data/psrm_cartel_data.csv) — Cleaned dataset for cartel outcomes.
        -   [`psrm_crime_data.csv`](smallest_eigenvalues/cleaned_data/psrm_crime_data.csv) — Cleaned dataset for homicide and theft outcomes.

    -   [`calculate_eigenvalue_for_design_matrix.m`](smallest_eigenvalues/calculate_eigenvalue_for_design_matrix.m) — MATLAB code that builds the sample analogue of the design matrices, computes their smallest eigenvalues, and writes the summary table.

    -   [`functions/`](smallest_eigenvalues/functions) — MATLAB functions used for implementing SSC method.

        -   [`att_ci.m`](smallest_eigenvalues/functions/att_ci.m), [`att_event_ci.m`](smallest_eigenvalues/functions/att_event_ci.m) — Implement SSC for overall ATT and event-time ATT effects with inference.
        -   [`ssc.m`](smallest_eigenvalues/functions/ssc.m) — Implements SSC for individual treatment effects and event-time ATT without inference.
        -   [`synthetic_control.m`](smallest_eigenvalues/functions/synthetic_control.m), [`synthetic_control_batch.m`](smallest_eigenvalues/functions/synthetic_control_batch.m) — Estimate synthetic control weights for one unit and for a batch of units.
        -   [`ssc_inference.m`](smallest_eigenvalues/functions/ssc_inference.m) — Implement SSC inference
        -   [`vline.m`](smallest_eigenvalues/functions/vline.m), [`vline_license.txt`](smallest_eigenvalues/functions/vline_license.txt) — Plotting utility and its license.

    -   [`output/`](smallest_eigenvalues/output) — Output folder for the eigenvalue calculations.

        -   [`Table1_smallest_eigenvalue.csv`](smallest_eigenvalues/output/Table1_smallest_eigenvalue.csv) — Table 1: smallest eigenvalues of the sample analogue of the design matrix for each outcome.

-   [`treatment_scheme/`](treatment_scheme) — Contains code and data that produce Figure 2 that visualize the timing of treatment adoption.

    -   [`cleaned_data/`](smallest_eigenvalues/cleaned_data) — Cleaned outcome data produced by [`data_cleaning.R`](data_cleaning.R).
        -   [`psrm_cartel_data.csv`](smallest_eigenvalues/cleaned_data/psrm_cartel_data.csv) — Cleaned dataset for cartel outcomes.
        -   [`psrm_crime_data.csv`](smallest_eigenvalues/cleaned_data/psrm_crime_data.csv) — Cleaned dataset for homicide and theft outcomes.
    -   [`plot_treatment_scheme.R`](treatment_scheme/plot_treatment_scheme.R) — Plot the treatment scheme for different outcome groups and save it as Figure 3.
    -   [`output/`](treatment_scheme/output) — Output folder for the treatment scheme figure.
        -   [`Figure2_treatment_scheme.png`](treatment_scheme/output/Figure2_treatment_scheme.png) — Figure that visualizes the treatment scheme for different outcome groups.

-   [`treatment_effect/`](treatment_effect) — Contains data and code to estimate treatment effects and to plot Figure 2.

    -   [`cleaned_data/`](smallest_eigenvalues/cleaned_data) — Cleaned outcome data produced by [`data_cleaning.R`](data_cleaning.R).

        -   [`psrm_cartel_data.csv`](smallest_eigenvalues/cleaned_data/psrm_cartel_data.csv) — Cleaned dataset for cartel outcomes.
        -   [`psrm_crime_data.csv`](smallest_eigenvalues/cleaned_data/psrm_crime_data.csv) — Cleaned dataset for homicide and theft outcomes.

    -   [`estimation_ssc.m`](treatment_effect/estimation_ssc.m) — MATLAB script that estimate treatment effects by SSC method and saves results.

    -   [`estimation_gsc.R`](treatment_effect/estimation_gsc.R) — R script that estimate treatment effects by GSC method and saves results.

    -   [`plot_results.m`](treatment_effect/plot_results.m) — R script that reads SSC and GSC results and produces Figure 3 that compares the methods.

    -   [`functions/`](smallest_eigenvalues/functions) — MATLAB functions used for implementing SSC method.

        -   [`att_ci.m`](smallest_eigenvalues/functions/att_ci.m), [`att_event_ci.m`](smallest_eigenvalues/functions/att_event_ci.m) — Implement SSC for overall ATT and event-time ATT effects with inference.
        -   [`ssc.m`](smallest_eigenvalues/functions/ssc.m) — Implements SSC for individual treatment effects and event-time ATT without inference.
        -   [`synthetic_control.m`](smallest_eigenvalues/functions/synthetic_control.m), [`synthetic_control_batch.m`](smallest_eigenvalues/functions/synthetic_control_batch.m) — Estimate synthetic control weights for one unit and for a batch of units.
        -   [`ssc_inference.m`](smallest_eigenvalues/functions/ssc_inference.m) — Implement SSC inference
        -   [`vline.m`](smallest_eigenvalues/functions/vline.m), [`vline_license.txt`](smallest_eigenvalues/functions/vline_license.txt) — Plotting utility and its license.

    -   [`output/`](treatment_effect/output) — Output folder for the SSC and GSC estimation.

        -   [`results_gsc.csv`](treatment_effect/output/results_gsc.csv) — Estimated treatment effects by [`estimation_gsc.R`](treatment_effect/estimation_gsc.R) and related statistics from the GSC approach.
        -   [`results_ssc.csv`](treatment_effect/output/results_ssc.csv) — Estimated treatment effects by [`estimation_ssc.m`](treatment_effect/estimation_ssc.m) and related statistics from the SSC approach.
        -   [`Figure3_application_results.png`](treatment_effect/output/Figure3_application_results.png) — Figure that compares the SSC and GSC approach, produced by [`plot_results.m`](treatment_effect/plot_results.m) .

### Requirements

-   R (4.5.1)
    -   required packages: dplyr, gsynth, panelView, ggplot2, cowplot
-   MATLAB (R2018b)
    -   required toolbox: Optimization Toolbox

### Replication Steps

1.  **Prepare cleaned data (R).**
    -   Run [`data_cleaning.R`](data_cleaning.R) in the top-level project directory.
        -   This script loads the raw data from [`raw_data/`](raw_data), applies the same sample restrictions as in the paper, and writes cleaned CSV files to [`smallest_eigenvalues/cleaned_data`](smallest_eigenvalues/cleaned_data), [`treatment_effect/cleaned_data`](treatment_effect/cleaned_data), and [`treatment_scheme/cleaned_data`](treatment_scheme/cleaned_data).
2.  **Compute smallest eigenvalues of the sample analogue of the design matrices (MATLAB).**
    -   Open directory [`smallest_eigenvalues/`](smallest_eigenvalues).
    -   Run [`calculate_eigenvalue_for_design_matrix.m`](smallest_eigenvalues/calculate_eigenvalue_for_design_matrix.m) by MATLAB
        -   This script reads the cleaned data and writes [`output/`](treatment_effect/output)[`Table1_smallest_eigenvalue.csv`](smallest_eigenvalues/output/Table1_smallest_eigenvalue.csv), summarizing the smallest eigenvalues for each outcome.
    -   Output CSV file will be saved in [`output/`](smallest_eigenvalues/output) directory.
3.  **Visualize treatment scheme (R):**
    -   Open directory [`treatment_scheme/`](treatment_scheme).
    -   Run [`plot_treatment_scheme.R`](treatment_scheme/plot_treatment_scheme.R) by R
        -   This script creates [`Figure2_treatment_scheme.png`](treatment_scheme/output/Figure2_treatment_scheme.png), which displays the timing of treatment adoption for homicide, theft, and cartel outcomes.
    -   Output figure will be saved in [`output/`](treatment_scheme/output) directory.
4.  **Estimate treatment effects (R and MATLAB).**
    -   Open directory [`treatment_effect/`](treatment_effect).
    -   Run [`estimation_gsc.R`](treatment_effect/estimation_gsc.R) by R.
        -   The script computes the results by GSC method that uses the "gsynth" package and the cleaned data in [`treatment_effect/cleaned_data`](treatment_effect/cleaned_data) to compute event-time ATT estimates and confidence intervals for all outcomes, and saves [`results_gsc.csv`](treatment_effect/output/results_gsc.csv).
    -   Run [`estimation_ssc.m`](treatment_effect/estimation_ssc.m) by MATLAB.
        -   The script computes the results by SSC method that uses the SSC functions in [`treatment_effect/functions`](treatment_effect/functions) and the cleaned data in [`treatment_effect/cleaned_data`](treatment_effect/cleaned_data) to compute event-time ATT estimates and confidence intervals for all outcomes, and saves [`results_ssc.csv`](treatment_effect/output/results_ssc.csv).
    -   Run [`plot_results.m`](treatment_effect/plot_results.m) by MATLAB.
        -   This script produces a figure that compares the two methods: it reads [`results_ssc.csv`](treatment_effect/output/results_ssc.csv) and [`results_gsc.csv`](treatment_effect/output/results_gsc.csv), plots figures that compare SSC and GSC event-time ATT estimates across all outcomes, and saves [`Figure3_application_results.png`](treatment_effect/output/Figure3_application_results.png).
    -   All output files will be saved in [`output/`](treatment_effect/output) directory.

Following these steps will reproduce the main empirical tables and figures implemented in this code package.