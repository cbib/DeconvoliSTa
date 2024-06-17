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
  # Vérifier que la longueur de alpha correspond au nombre de colonnes
  if (length(alpha) != length(columns)) {
    stop("La longueur de alpha doit correspondre au nombre de colonnes.")
  }
  
  # Générer les valeurs selon la distribution de Dirichlet
  dirichlet_values <- rdirichlet(rows, alpha)
  
  # Convertir en data.frame et nommer les colonnes
  df <- as.data.frame(dirichlet_values)
  colnames(df) <- columns
  
  return(df)
}

save_dataframe_to_tsv <- function(df, file_path) {
  # Sauvegarder le DataFrame dans un fichier TSV
  fwrite(df, file_path, sep = "\t")
}

sp_input <- readRDS(par$sp_input)


# Paramètres
print(length(sp_input$gold_standard_priorregion$celltype))
alpha_value = 0.5
alpha <- rep(alpha_value, length(sp_input$gold_standard_priorregion$celltype))
rows <- sp_input$counts@Dim[2]
columns <- sp_input$gold_standard_priorregion$celltype  


# Générer le DataFrame et sauvegarder dans un fichier TSV
df <- generate_dirichlet_dataframe(alpha, rows, columns)
file_path <- par$output  # Spécifiez le chemin de votre fichier de sortie

# Remove all spaces and dots from cell names, sort them
colnames(df) <- stringr::str_replace_all(colnames(df), "[/ .&-]", "")
# deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

# if (nrow(deconv_matrix) != ncol(spatial_data$counts)){
#     cat ("okkkkkk\n")
#   message("The following rows were removed, possibly due to low number of genes: ",
#           paste0("'", colnames(spatial_data$counts)[!colnames(spatial_data$counts) %in% rownames(deconv_matrix)], "'", collapse=", "))
# }


save_dataframe_to_tsv(df, file_path)

print("Fichier TSV généré avec succès.")
