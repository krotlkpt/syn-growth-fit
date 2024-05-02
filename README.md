# Synechocystis growth fit

## Overview

This repository contains Jupyter notebooks for fitting the _Synechpcystis_ sp. PCC 6803 model from the publication 'A quantitative description of light-limited cyanobacterial growth using flux balance analysis.' to experimental data and generating corresponding plots.
Use 'n_Data' to adjust the amount of data points to simulate (this can lead to long ciputing times)

## Usage
### Download data or use local

By setting the parameter'download' to "**True**", the notebook will download the model sbml file and the experimental data dircetly from source. By dafault, local versions of the model and data are used.

### Adjusting Data Points

You can control the number of data points to simulate by adjusting the n_Data parameter in the notebook. Note that a higher number of data points may result in longer computation times.

### Customizing Plot Style

Feel free to customize the style of the plots to your preference. You can modify plot styles within the notebooks to achieve the desired visualizations.

### Saving Plots

If you wish to save the generated plots, simply change the savefigs parameter to **True** in the notebook. Additionally, define the image format using the img_type parameter.

It's also possible to change parameter values to inspect the effect on the simulations after fitting.

## Getting Started

To get started, simply open the desired notebook and follow the instructions provided within. Experiment with different parameter values and plot styles to gain insights into Synechocystis growth dynamics.

Run in Binder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fitbgit.biologie.hu-berlin.de%2Fhoeper%2Fsyn-growth-fit/binder)
