library(Seurat)
H5_FILENAME <- "filtered_feature_bc_matrix.h5"
ampulla_path <- "/workspaces/cell-cell-interaction/kim-data/E.794/analyzed/Sample_Visium_13_41_S1_b/outs"
# isthmus_path <- "/workspaces/cell-cell-interaction/kim-data/E.794/analyzed/Sample_Visium_13_38_S1_b/outs"
OUT_DIR <- "/workspaces/cell-cell-interaction/compbio-hw/data/"
n_gene_vec <- c(3000, 200)
####################################################
# Ampulla
####################################################

for (n_gene in n_gene_vec){
  out_dir <- paste0(OUT_DIR, n_gene)
  # Load the Visium data
  # Replace 'path_to_data' with the actual path to your Visium dataset
  ampulla_obj <- Load10X_Spatial(
    ampulla_path,
    filename = H5_FILENAME,
    assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE
  )
  
  ampulla_obj <- SCTransform(ampulla_obj, assay = "Spatial", variable.features.n = n_gene)
  # An object of class Seurat 
  # 36601 features across 1870 samples within 1 assay 
  # Active assay: Spatial (36601 features, 0 variable features)
  #  1 layer present: counts
  #  1 image present: slice1
  
  # VariableFeatures(ampulla_obj)
  
  # File 1: Spatial Coordinates of barcodes ----------
  coordinates_df <- ampulla_obj@images[["slice1"]]@coordinates[, c("row", "col")]
  # head(coordinates_df)
  #                     row col
  # AAACCGGGTAGGTACC-1  42  28
  # AAACCGTTCGTCCAGG-1  52  42
  # AAACCTCATGAAGTTG-1  37  19
  # AAACGAAGAACATACC-1   6  64
  # AAACGAGACGGTTGAT-1  35  79
  # AAACTGCTGGCTCCAA-1  45  67
  
  file_name <- paste0(out_dir, "/ampulla_coordinates.csv")
  write.csv(coordinates_df, file_name)
  
  # File 2: Raw Transcriptomics Counts ----------
  raw_assay <- subset(
    ampulla_obj@assays[["Spatial"]],
    features = VariableFeatures(ampulla_obj)
  )
  # Extracting the count matrix
  raw_count_mtx <- raw_assay@layers[["counts"]]
  # Converting count matrix to data frame
  raw_count_df <- as.data.frame(as.matrix(raw_count_mtx))
  rownames(raw_count_df) <- rownames(raw_assay)
  colnames(raw_count_df) <- colnames(raw_assay)
  # Reorder the rows by ariableFeatures(ampulla_obj)
  raw_count_df <- raw_count_df[VariableFeatures(ampulla_obj), ]
  # raw_count_df[1:5, 1:5]
  #       AAACCGGGTAGGTACC-1 AAACCGTTCGTCCAGG-1 AAACCTCATGAAGTTG-1 AAACGAAGAACATACC-1 AAACGAGACGGTTGAT-1
  # IGLC2                 47                 11                  4                  9                 26
  # IGHG3                 17                  5                  2                  9                 16
  # IGKC                  14                 10                  3                  6                 16
  # MT2A                   8                  0                  0                  2                  2
  # IGHG1                 15                  3                  0                  3                 10
  
  file_name <- paste0(out_dir, "/ampulla_raw_counts.csv")
  write.csv(raw_count_df, file_name)
  
  # File 3: Normalized Transcriptomics Counts ----------
  sct_count_mtx <- ampulla_obj@assays[["SCT"]]@scale.data
  # Reorder the rows by ariableFeatures(ampulla_obj)
  sct_count_mtx <- sct_count_mtx[VariableFeatures(ampulla_obj), ]
  # sct_count_mtx[1:5, 1:5]
  #       AAACCGGGTAGGTACC-1 AAACCGTTCGTCCAGG-1 AAACCTCATGAAGTTG-1 AAACGAAGAACATACC-1 AAACGAGACGGTTGAT-1
  # IGLC2           6.603713         -1.0478664          0.2023871         -1.3896850          -2.322335
  # IGHG3           2.157385         -1.0083370          0.1320593          0.3247651          -1.501051
  # IGKC            1.855709          0.9357456          1.3652242         -0.2063727          -1.063739
  # MT2A            2.452951         -2.0077930         -1.4606833         -0.7388417          -2.017789
  # IGHG1           2.967789         -0.8634551         -1.1905125         -0.7896924          -1.213939
  
  file_name <- paste0(out_dir, "/ampulla_normalized_counts.csv")
  write.csv(sct_count_mtx, file_name)
}

ampulla_obj_test <- Load10X_Spatial(
  ampulla_path,
  filename = H5_FILENAME,
  assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE
)
ampulla_obj_test <- SCTransform(ampulla_obj_test, assay = "Spatial", variable.features.n = n_gene)

ampulla_obj_test <- FindVariableFeatures(ampulla_obj_test, nfeatures=200)

rownames(ampulla_obj_test)
