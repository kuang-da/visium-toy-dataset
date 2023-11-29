library(Seurat)
H5_FILENAME <- "filtered_feature_bc_matrix.h5"
ampulla_path <- "/workspaces/cell-cell-interaction/kim-data/E.794/analyzed/Sample_Visium_13_41_S1_b/outs"
# isthmus_path <- "/workspaces/cell-cell-interaction/kim-data/E.794/analyzed/Sample_Visium_13_38_S1_b/outs"
OUT_DIR <- "/workspaces/cell-cell-interaction/compbio-hw/data"
####################################################
# Ampulla
####################################################

# Load the Visium data
# Replace 'path_to_data' with the actual path to your Visium dataset
ampulla_obj <- Load10X_Spatial(
  ampulla_path,
  filename = H5_FILENAME,
  assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE
)
ampulla_obj <- SCTransform(ampulla_obj, assay = "Spatial")
VariableFeatures(ampulla_obj)

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

file_name <- paste0(OUT_DIR, "/ampulla_coordinates.csv")
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
# raw_count_df[1:5, 1:5]
#           AAACCGGGTAGGTACC-1 AAACCGTTCGTCCAGG-1 AAACCTCATGAAGTTG-1 AAACGAAGAACATACC-1 AAACGAGACGGTTGAT-1
# LINC02593                  0                  0                  0                  0                  0
# SAMD11                     0                  0                  0                  0                  0
# HES4                       1                  0                  0                  1                  2
# ISG15                      2                  0                  0                  0                  0
# B3GALT6                    0                  0                  0                  0                  1

file_name <- paste0(OUT_DIR, "/ampulla_raw_counts.csv")
write.csv(raw_count_df, file_name)

# File 3: Normalized Transcriptomics Counts ----------
sct_count_mtx <- ampulla_obj@assays[["SCT"]]@scale.data
# sct_count_mtx[1:5, 1:5]
#           AAACCGGGTAGGTACC-1 AAACCGTTCGTCCAGG-1 AAACCTCATGAAGTTG-1 AAACGAAGAACATACC-1 AAACGAGACGGTTGAT-1
# LINC02593         -0.1400006         -0.1278043        -0.02688107         -0.1205156         -0.2921160
# SAMD11            -0.3872973         -0.3724581        -0.14972631         -0.3631890         -0.6482071
# HES4               1.0965559         -0.5304079        -0.23695647          1.2786578          0.6489449
# ISG15              6.3440350         -0.2967609        -0.09539067         -0.2889317         -0.5301926
# B3GALT6           -0.2464190         -0.2357010        -0.05483691         -0.2290093          1.5370411

file_name <- paste0(OUT_DIR, "/ampulla_normalized_counts.csv")
write.csv(sct_count_mtx, file_name)