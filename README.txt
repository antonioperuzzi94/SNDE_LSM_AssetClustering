TITLE: "A Dynamic Latent-Space Model for Assets Clustering"

AUTHORS:      Roberto Casarin & Antonio Peruzzi

AVAILABLE AT:    Studies in Nonlinear Dynamics and Econometrics


DATE:        October 2023


Tested on R version 4.2.2 (2022-10-31) -- "Innocent and Trusting"

-------------------------------------------------------------------------------
The following Repository contains the files (scripts and data) used to reproduce 
the results of the paper "A Dynamic Latent-Space Model for Assets Clustering".


-------------------------------------------------------------------------------

        %%%%%%%%%%%%%% PRELIMINARIES  %%%%%%%%%%%%%%

-------------------------------------------------------------------------------

- The folder "Code" contains the R scripts;

- The folder "Data" contains the dataset used in this work;

- The folder "Figures" contains the figures in the manuscript.

-------------------------------------------------------------------------------

        %%%%%%%%%%%%%%  R SCRIPTS  %%%%%%%%%%%%%%

-------------------------------------------------------------------------------

* SNDE_01_motivation_plot.R

This script is used to reproduce Figure 1

* SNDE_02_properties_plot.R

This script is used to reproduce Figure 2

* SNDE_03_synthetic_data.R

This script is used to generate the simulation dataset in Section 3.3 and reproduce Figure 4

* SNDE_04_simulation.R

This script is used to reproduce Figure 5 and 6 and 10 and Table 1

* SNDE_05_SP100_model.R

This script is used to run the mcmc algorithm on the SP100 data

* SNDE_06_DAX40_model.R

This script is used to run the mcmc algorithm on the DAX40 data

* SNDE_07_Results_Plots.R

This script is used to reproduce Figure 8 and 9 and Table 2 and 3

* SNDE_08_sequential_Robustness_SP100.R

This script is used to run the mcmc algorithm on the SP100 data for the sequential robustness check

* SNDE_09_sequential_Robustness_DAX40.R

This script is used to run the mcmc algorithm on the DAX40 data for the sequential robustness check

* SNDE_10_Robustness_Plots.R

This script is used to reproduce Figure 12


-------------------------------------------------------------------------------

        %%%%%%%%%%%%%%  C++ Code  %%%%%%%%%%%%%%

-------------------------------------------------------------------------------

* Utils.cpp

This script contains some useful c++ functions such as "procrustes_preprocessing" used to carry out the procrustes transformation pre-processing.
