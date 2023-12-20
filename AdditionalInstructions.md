# Additional Instructions

## Installing SPOTlight

Before installing SPOTlight, in your conda environment (where your R kernel is): please install Libxml2 in your environment, like so:

```
conda activate (my_r_environment)
conda install -c conda-forge libxml2
```

In order to read h5 objects in R, the `rhdf5` package is required. It can be installed in R using:

```
BiocManager::install("rhdf5")
```

The `anndata` package allows for the reading in of `h5ad` files. This is installed in R simply using:
```
install.packages("anndata")
```