
Data and code for: Structural asymmetry in biotic interactions as a tool to understand and predict ecological persistence

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8087366.svg)](https://doi.org/10.5281/zenodo.8087366)

### Summary

The raw data is available in the /Data folder, in .csv files that contain the information needed to create the interaction matrices for our empirical examples. Every other manipulation is stored in the /Results folder.

All results and figures from the study are already generated in the relevant folders. To re-run the analyses, run each script in order. Scripts with prefix 01 are used to generate data from the empirical community matrices; scripts with prefix 02 are used to process and derive different metrics and results from 01_scripts; scripts with prefix 03 analyse these metrics and provide the inputs for the tables; and, finally, scripts with prefix 04 are use to plot figures.

In particular, figures and tables of the main text can be generated via the following scripts:

- Figure 1: "04_plot_fig1_fig2.ipynb"
- Figure 2: "04_plot_fig1_fig2.ipynb"
- Figure 3: "04_plot_fig3.R"
- Figure 4: "04_plot_fig4a.ipynb" and "04_plot_fig4b.R"
- Figure 5: "04_plot_fig5.R"

##### Notes:

- there are different auxiliary scripts both in the main /Scripts folder and in one subfolder.
- some of the scripts are very computationally demanding, especially those scripts with the prefix 01. 
- the Box2_02 scripts create files of such size (~3GB) that they could not be stored in this public repository (the outputs of Box2_02 files are the inputs of Box2_03 and Box2_04 files).
- the scripts used to represent three-dimensional images (Figs. 1, 2, 4a, S3, S4, S5 and S8) are in python (jupyter notebook), but all data processing and analysis is in R files.
