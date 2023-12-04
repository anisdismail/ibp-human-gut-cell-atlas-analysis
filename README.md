# Human Gut Cell Atlas Analysis Pipeline

## Introduction

This repository contains the analysis pipeline for the Human Gut Cell Atlas project. The Human Gut Cell Atlas is a comprehensive resource that aims to profile the cell types and states in the human gut.

## Table of Contents

- [Features](#features)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Features

- High-throughput analysis of single-cell RNA sequencing (scRNA-seq) data from the human gut.
- Identification and annotation of cell types and states.
- Visualization of the cell atlas and analysis results.

## Getting Started

### Prerequisites

- List any prerequisites here, such as software dependencies, data sources, and tools.

### Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/human-gut-cell-atlas-analysis.git
   ```

2. Change into the project directory:

   ```bash
   cd human-gut-cell-atlas-analysis
   ```

3. Install the required dependencies:

   ```bash
   pip install -r requirements.txt
   ```

## Usage

Provide detailed instructions on how to use the analysis pipeline. Include examples and usage scenarios, if applicable. 
# Creating R kernel
Run this code in your terminal to set up the R kernel for running any jupyter notebooks that require R.
```bash
conda create -n r_kernel
conda activate r_kernel
conda config --append channels conda-forge
conda install r-recommended r-irkernel
conda install Jupyter
conda install -c conda-forge nb_conda_kernels
conda install -c conda-forge r-dplyr
conda run -n r_kernel Rscript -e 'IRkernel::installspec(name="r_kernel", displayname="R 3.3")'
```
Alternative more concise code to set up the same kernel.
```bash
# Create a Conda environment with R 4.3.2
conda create -n r_kernel r-base=4.3.2 r-irkernel jupyter nb_conda_kernels r-dplyr

# Activate the environment
conda activate r_kernel

# Install IRkernel in the R environment
Rscript -e 'IRkernel::installspec(name="r_kernel", displayname="R 4.3.2")'
```
### Alternatives to installing R packages via github
Method 1. Use Bioconductor instead of github. The newest version (3.19) of Bioconductor asks for Rv4.4, which has not been released yet. Instead, Bioconductor 3.18 will work with Rv4.3.
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version=3.18)

BiocManager::install("STdeconvolve")
```
Method 2. Subsistute your own username/password where appropriate.

```R
install.packages("usethis")

usethis::use_git_config(user.name = "name", user.email = "name@gmail.com")

  

# The above command should send you to a URL to generate

# Personal Access Token (PAT)

# Make sure to copy down the PAT for use below

usethis::create_github_token()

usethis::use_git_config(helper="cache --timeout=2600000") #> cache timeout ~30 days

# Insert PAT into the command below

pat<-readline(prompt="please copy PAT here")

credentials::set_github_pat("pat")

usethis::edit_r_environ()

  

usethis::git_sitrep()

# The above command should contain something along the lines of

# Personal access token: <found in env var>

  

# Now for the seurat stuff
install.packages("remotes")

remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)

# If the above command fails, run this one and try running it again

#remotes::install_github("mojaveazure/seurat-object", "seurat5")

remotes::install_github("mojaveazure/seurat-disk")

  

devtools::install_github('satijalab/seurat-data')
```

## Contributing

If you would like to contribute to this project, please follow the guidelines in [CONTRIBUTING.md](CONTRIBUTING.md).

## License

This project is licensed under the [License Name] - see the [LICENSE.md](LICENSE.md) file for details.
```
