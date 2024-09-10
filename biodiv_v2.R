################################################################################
##                              run biodivMapR_v2                             ##
################################################################################
# clean environment
rm(list=ls(all=TRUE));gc()
# Define working directory as the directory where the script is located
if (rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# ################################################################################
# # Run following code once to install bigRaster package
#
# devtools::install_gitlab('jbferet/bigRaster')
################################################################################
library(biodivMapR)
ls("package:biodivMapR")


################################################################################
# path for the image to process
input_raster_path <- '/Users/Oleksandr/biodivMapR_v2/est_0524/est0824.tif'

################################################################################
# parameters for biodivMapR
# Define path for master output directory where files produced during the process are saved
output_dir <- '/Users/Oleksandr/biodivMapR_v2/RESULTS_biodivMapR_v2/est0524'
dir.create(output_dir,showWarnings = F, recursive = T)
# window size for computation of spectral diversity
window_size <- 10
# computational parameters
nbCPU <- 8
maxRows <- 500
################################################################################
##                    produce a mask based on NDVI thresholding               ##
################################################################################
# get template for S2 data
S2_hdr_file <- system.file("extdata", "HDR/SENTINEL_2.hdr", package="biodivMapR")
S2_hdr <- read_ENVI_header(HDRpath = S2_hdr_file)
input_mask_path <- NULL
# produce a radiometric filter which will be stored in 'output_dir/mask_update'
updated_mask_path <- radiometric_filtering(input_raster_path = input_raster_path,
                                           output_dir = output_dir,
                                           input_rast_wl = S2_hdr$wavelength,
                                           input_mask_path = NULL,
                                           NDVI_Thresh = 0.6, Blue_Thresh = 500, NIR_Thresh = 1500)

PCA_info <- perform_PCA(input_raster_path = input_raster_path,
                        output_dir = output_dir,
                        input_rast_wl = S2_hdr$wavelength,
                        input_mask_path = updated_mask_path,
                        Continuum_Removal = TRUE,
                        TypePCA = 'SPCA', maxRows = maxRows)

# path of the raster resulting from dimensionality reduction
PCA_Files <- PCA_info$PCA_Files
# path for the updated mask
Input_Mask_File <- PCA_info$MaskPath

Selected_PC <- select_PCA_components(output_dir = output_dir,
                                     pca_rast_path = PCA_info$PCA_Files$PCA)

################################################################################
##            run biodivMapR on theselection of spectral indices              ##
################################################################################
ab_info <- biodivMapR_full(input_raster_path = PCA_info$PCA_Files$PCA,
                           SelectBands = Selected_PC,
                           output_dir = output_dir,
                           window_size = window_size,
                           input_mask_path = updated_mask_path,
                           nbCPU = nbCPU, maxRows = maxRows)

