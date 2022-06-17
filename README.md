# Feature Selection using e-values

This repository contains code and data for the experiments in the paper

> [Feature Selection using e-values](https://arxiv.org/abs/2206.05391)
>
> Subhabrata Majumdar and Snigdhansu Chatterjee
>
> International Conference on Machine Learning (ICML), 2022.

## Setup

You should have the [R statistical software](https://www.r-project.org) installed in your system. To install the dependencies, run the following in bash:
```
Rscript install_dependencies.R
```

Running the code for the Mixed Integer Optimization (MIO) method requires installation of the [Gurobi Optimizer](https://www.gurobi.com/documentation/9.5/quickstart_windows/software_installation_guid.html), available under free academic licenses, and the R package by provided by the authors:
```
Rscript -e "install.packages('devtools')"
Rscript -e "devtools::install_github(repo='ryantibs/best-subset', subdir='bestsubset')"
```

## Contents
The different subdirectories of this repository contain implementation codes for specific experiments. Details are available in README files inside the respective subdirectory.

| Subdirectory | Method |
|---|---|
| `lm` | Linear models |
| `lmhd` | High-dimensional linear models |
| `lmm` | Linear mixed models |
| `monsoon` | Indian monsoon data analysis |
| `fmri` | fMRI data analysis |

## Contact
Please feel free to [email](mailto:zoom.subha@gmail.com) with questions of comments.
