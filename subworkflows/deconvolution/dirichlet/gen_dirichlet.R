Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(Matrix)
library(Seurat)

par <- list(
    cell_min = 5
)

library(gtools)
library(data.table)


# Replace default values by user input
args <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)
par[names(args)] <- args
print(par)


generate_dirichlet_dataframe <- function(alpha, rows, columns) {
  # Check that the length of alpha matches the number of columns
  if (length(alpha) != length(columns)) {
    stop("The length of alpha must match the number of columns.")
  }
  
  # Generate values from the Dirichlet distribution
  dirichlet_values <- rdirichlet(rows, alpha)
  
  # Convert to a data.frame and name the columns
  df <- as.data.frame(dirichlet_values)
  colnames(df) <- columns
  
  return(df)
}

save_dataframe_to_tsv <- function(df, file_path) {
  # Save the data.frame to a TSV file
  fwrite(df, file_path, sep = "\t")
}

sp_input <- readRDS(par$sp_input)


# Parameters
print(length(sp_input$gold_standard_priorregion$celltype))
alpha_value = 0.5
alpha <- rep(alpha_value, length(sp_input$gold_standard_priorregion$celltype))
rows <- sp_input$counts@Dim[2]
columns <- sp_input$gold_standard_priorregion$celltype  


# Generate the data.frame and save it to a TSV file
df <- generate_dirichlet_dataframe(alpha, rows, columns)
file_path <- par$output  # path to the output file

# Remove all spaces and dots from cell names, sort them
colnames(df) <- stringr::str_replace_all(colnames(df), "[/ .&-]", "")
# deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

# if (nrow(deconv_matrix) != ncol(spatial_data$counts)){
#     cat ("okkkkkk\n")
#   message("The following rows were removed, possibly due to low number of genes: ",
#           paste0("'", colnames(spatial_data$counts)[!colnames(spatial_data$counts) %in% rownames(deconv_matrix)], "'", collapse=", "))
# }


save_dataframe_to_tsv(df, file_path)

print("TSV file generated successfully.")
