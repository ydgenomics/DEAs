# /software/conda/Anaconda/bin/R 250220
# cloud-image: convert--08
library(optparse)
library(sceasy)
library(reticulate)
option_list <- list(
  make_option(c("-i", "--input_file"),
    type = "character", default = NULL,
    help = "Path to input file for convrting"
  ),
  make_option(c("-o", "--output_file"),
    type = "character", default = NULL,
    help = "Path to output file for convrting"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

use_python("/software/conda/Anaconda/bin/python")
loompy <- reticulate::import('loompy')
# opt$input_file <- "/data/work/bbknn/Cer_test_convert_BBKNNR_integrated.rds"
temp0 <- readRDS(opt$input_file)
print(temp0)

colnames(temp0@meta.data)

temp0[["RNA"]] <- as(temp0[["RNA"]], "Assay")

#
sceasy::convertFormat(temp0, from="seurat", to="anndata",assay = "RNA",outFile = opt$output_file)