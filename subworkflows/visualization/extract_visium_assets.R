#!/usr/bin/env Rscript
# Extract everything the interactive visualization needs directly from a Visium Seurat .rds,
# so the visualization can run automatically (no manually-prepared image / coordinates):
#   - tissue_image.png      : the tissue image stored in the object
#   - tissue_positions.csv   : barcode, in_tissue, array_row/col, pxl_row/col_in_fullres
#   - scale_factor.txt       : factor mapping full-res pixel coords onto that image
#
# Usage:
#   Rscript extract_visium_assets.R --sp_input sample.rds --out_dir vis_assets
#
# If the object has no image (@images empty), this exits with an error: in that case run the
# pipeline without the visualization step (do_visu=false).

suppressMessages({
  library(Seurat)
  library(png)
})

par <- list(sp_input = NA, out_dir = "vis_assets")
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
par[names(args)] <- args
dir.create(par$out_dir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(par$sp_input)
if (length(obj@images) == 0)
  stop("The spatial object has no image (@images is empty); cannot build the visualization. ",
       "Re-run with do_visu=false.")
im <- obj@images[[1]]

# --- Spot positions (full-res pixel coordinates) ---
co <- im@coordinates    # columns: tissue, row, col, imagerow, imagecol
coords <- data.frame(
  barcode            = rownames(co),
  in_tissue          = co$tissue,
  array_row          = co$row,
  array_col          = co$col,
  pxl_row_in_fullres = co$imagerow,
  pxl_col_in_fullres = co$imagecol
)
# No header / no row names: the visualizer reads this with header=None and sets columns by position.
write.table(coords, file.path(par$out_dir, "tissue_positions.csv"),
            sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

# --- Tissue image (@image is an H x W x 3 array in [0,1]) ---
img <- im@image
png::writePNG(img, file.path(par$out_dir, "tissue_image.png"))
H <- nrow(img); W <- ncol(img)

# --- Scale factor ---
# Pick the scalefactor (hires or lowres) under which the full-res coordinates actually fit the
# stored image -> that is the resolution of the image embedded in the object. Robust whether the
# object kept the hires or the lowres image.
sf <- im@scale.factors
cand <- c(hires = sf$hires, lowres = sf$lowres)
fits <- sapply(cand, function(s) (max(co$imagerow) * s <= H * 1.02) && (max(co$imagecol) * s <= W * 1.02))
scale_factor <- if (any(fits)) unname(cand[which(fits)[1]]) else (W / max(co$imagecol))
writeLines(format(scale_factor, digits = 10), file.path(par$out_dir, "scale_factor.txt"))

cat("Wrote:", par$out_dir, "| image", H, "x", W,
    "| spots", nrow(coords), "| scale_factor", scale_factor, "\n")
