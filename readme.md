# Readme
This repository contains a simplified dataset as a practice of data analysis. The data is based on a Visium slide of Ampulla but for simplicity we only provide the highly variable genes and the spatial coordinates of the spots. The data files are located in the data folder. Genes are ranked by the variance of residuals from SCTransform.

- `data/ampulla_coordinates.csv`: Spatial coordinates for each spot.
- `data/ampulla_raw_counts.csv`: Raw gene count data for each spot.
- `data/ampulla_normalized_counts.csv`:  Normalized gene count data for each spot.
- `ampulla_13_41_S1_b.png`: a figure showing the basic information of the Visium slide. For left to right, the first figure is the H&E of the slide colored by the counts of all gene in each spot. Note that we only provide a subset of genes. The second figure is UMAP of the Visium spots. We use the `ampulla_raw_counts.csv` to calculate the top50 PCs and then calculate the UMAP. The spots are clustered and colored by groups. The thrid figure is the H&E image but colored by the groups.
- `dataset.R`: A script explaining the data generation process.
