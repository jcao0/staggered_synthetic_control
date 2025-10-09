## Overview

This repository contains replication files for Cao, Lu and Wu (2020), "Synthetic Control Inference for Staggered Adoption: Estimating the Dynamic Effects of Board Gender Diversity Policies."

There are two main parts to the replication:

1.  Data cleaning
2.  Replication of tables and figures, including Table 1, Table 3 and all the figures in the main text and appendix.

## Date availability

All data sources are publicly available. Data used are listed below:

1.  Data on corporate board’s female representation across EU countries from European Institute for Gender Equality: this data is stored in *wmidm_bus_bus\_\_wmid_comp_compbm.xlsx* and can be downloaded from <https://eige.europa.eu/gender-statistics/dgs/indicator/eustrat_ges_lead_bus__wmid_comp__ind53_top/datatable>.
2.  Labor outcome variables across EU countries from the EU Labor Force Survey (“LFS”):  (1) full-time employment data is stored in \*lfsq_epgais.xlsx\* and can be downloaded from <https://ec.europa.eu/eurostat/databrowser/view/LFSQ_EPGAIS__custom_6036320/default/table>. (2) working hours data is stored in *lfsq_ewhuis.xlsx* and can be downloaded from <https://ec.europa.eu/eurostat/databrowser/view/lfsq_ewhuis/default/table>.
3.  Announcement dates of the board gender policies: this data is stored \*policy.dta\* and is hand-collected from each EU country’s official website, with additional guidance from Seierstad et al. (2017).


## Folder content

This replication package is organized according to the following structure:

```         
.
├── README.md                               
└── replication_files/
    |── data cleaning/
    |   |── EIGE board gender data/                
    |   |   |── policy.dta                  Raw dataset of policy announcement data
    |   |   |── wmidm_bus_bus__wmid_comp_compbm.xlsx  Raw datasets for corporate board’s female ratio
    |   |   |── jasa board gender clean.do  Stata script for cleaning board’s female ratio
    |   |   |── jasa_boardgendereige.dta    Cleaned Stata dataset.
    |   |   └── data_boardgendereige.csv    Cleaned stata dataset stored as csv file
    |   |── LFS full time employment data/
    |   |   |── policy.dta                  Raw dataset of policy announcement data
    |   |   |── lfsq_epgais.xlsx            Raw dataset of full-time employment
    |   |   |── jasa ft employment clean.do  Stata script for cleaning full-time employment data.
    |   |   |── jasa_ft.dta                 Cleaned Stata dataset.
    |   |   └── data_ft.csv                 Cleaned dataset stored as csv file.
    |   └── LFS work hours data/
    |       |── policy.dta                  Raw dataset of policy announcement data
    |       |── lfsq_ewhuis.xlsx            Raw dataset of working hours
    |       |── jasa work hours clean.do    Stata script for cleaning working hours data.
    |       |── jasa_hr.dta                 Cleaned Stata dataset.
    |       └── data_hr.csv                 Cleaned dataset stored as csv file
    ├── Figure 1/                           
    │   ├── data_boardgendereige.csv        Board gender ratio as input data for Figure 1.
    │   ├── event_time_att.m                MATLAB script for event-time ATT estimation with our method
    │   ├── functions/                      Core MATLAB functions for our synthetic control method (SSC).
    │   │   ├── synthetic_control.m         Classic synthetic control weights algorithms for one treated unit.
    │   │   ├── synthetic_control_batch.m   Our synthetic control weights algorithms for multiple treated unit.
    │   │   ├── ssc.m                       Estimation
    │   │   ├── ssc_inference.m             Inference
    │   │   ├── att_event_ci.m              ATT confidence intervals.
    │   │   ├── vline.m, vline_license.txt  Plot helper and license.
    │   ├── graph1a.png, graph1b.png        Output plots for Figure 1.
    ├── Figure 2–6/                         Same structure as Figure 1, using different datasets and scripts.
    ├── Figure A1/                          
    │   ├── jasa graphA1.do                 Stata script for Figure A1.
    │   ├── jasa_boardgendereige.dta        Board gender ratio as as input stata data file.
    │   └── figureA1.png                    Output plot for Figure A1.
    ├── Figure A2/                          
    │   ├── results_staggered_did.R, results_staggered_did.csv R code and results for ATT estimation with staggered DID
    │   ├── results_staggered_synthetic_control.m MATLAB script for ATT estimation with our method
    │   ├── data_boardgendereige.csv        Board gender ratio as input dataset.
    │   ├── functions/                      Core MATLAB functions for our method (same as ./replication_files/Figure 1/functions).
    │   └── FigureA2.png                    Output plot for Figure A.2.
    ├── Table 1/                         
    │   ├── jasa_boardgendereige.dta        Data on Board gender ratio.
    │   ├── jasa_ft.dta                     Data on full-time employment rate
    │   ├── jasa_hr.dta                     Data on working hours
    │   ├── jasa table1.do                  Stata scripts for Table 1.
    │   ├── *.tex, *.txt                    Exported LaTeX and text tables.
    └── Table 3/                         
        ├── data_boardgendereige.csv        Board gender ratio as input data
        ├── event_time_att.m                MATLAB scripts for estimation.
        ├── functions/                      Core MATLAB functions for our method (same as./replication_files/Figure 1/functions).
        └── table3.csv                      Output data for Table 3

```

## Replication instructions

To replicate the cleaned data , one can navigate to the  *replication_files/data cleaning/* folder and run the appropriate cleaning files (see above). The origianl process was conducted using Stata 14.2. 

Tables and figures can be immediately replicated by navigating to the corresponding folder in *replication_files/* and running the relevant code as shown below:

| Figure/Table \# | Original study software   | Relevant code |
|-----------------|---------------------------|----------------------------------------|
| Figure 1        | MATLAB R2018b             | Figure 1/event_time_att.m |
| Figure 2        | MATLAB R2018b             | Figure 2/event_time_att_comparison.m |
| Figure 3        | MATLAB R2018b             | Figure 3/event_time_att_comparison_other.m |
| Figure 4        | MATLAB R2018b             | Figure 4/event_time_att_comparison_other.m |
| Figure 5        | MATLAB R2018b             | Figure 5/event_time_att_comparison_other.m |
| Figure 6        | MATLAB R2018b             | Figure 6/event_time_att_comparison_other.m |
| Figure A.1      | Stata 14.2                | Figure A1/jasa graphA1.do |
| Figure A.2[^1]  | R 4.5.1 and MATLAB R2018b | Figure A2/results_staggered_did.R <br> Figure A2/results_staggered_synthetic_control.m |
| Table 1         | Stata 14.2                | Table 1/jasa table1.do |
| Table 3         | MATLAB R2018b             | Table 3/event_time_att.m |

[^1]: To replicate Figure A.2, first run the R code *results_staggered_did.R* to produce results using the staggered DID method. Then, run the Matlab code *results_staggered_synthetic_control.m* to produce results using our method and perform the comparison.

## Data Reference

Seierstad, Cathrine, Patricia Gabaldón, and Heike Mensi-Klarbach. 2017. Gender diversity in the boardroom. Volume 1, The use of different quota regulations. Palgrave Macmillan.
