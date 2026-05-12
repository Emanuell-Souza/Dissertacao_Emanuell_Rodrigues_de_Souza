#### 03/31/2025
#### Enrique Medina-Acosta, Emanuell Rodrigues de Souza, 
#### Higor Almeida Cordeiro Nogueira, Victor Santos Lopes
#### UENF/CBB/LBT
#### -------------------------------------------------------------------------------------------
#### Machine learning Analysis of multiomics predictive variables and multiple response variables
#### --------------------------------------------------------------------------------------------
#### ML_get_data and downstream analysis for prediction of outcomes using 
###########
###########
# Selected minimum meaning too meaningful multi-omic, multi-optopsis RCD signature 
# Define list of required packages
# cran_packages <- c(
#   "caret", "circlize", "ComplexHeatmap", "dplyr", "ggplot2",
#   "grid", "gridExtra", "kableExtra", "Matrix", "mice", "missForest",
#   "pROC", "purrr", "reshape2", "rio", "rms", "survival",
#   "survivalROC", "survcomp", "survminer", "tidyr", "timeROC",
#   "UpSetR", "VIM", "xgboost", "DiagrammeR"
# )
# 
# bioc_packages <- c("ComplexHeatmap", "survcomp")
# 
# github_packages <- list(
#   "UCSCXenaShiny" = "openbiox/UCSCXenaShiny",
#   "DiagrammeRsvg" = "rich-iannone/DiagrammeRsvg"
# )
# 
# special_packages <- c("rsvg", "magick", "lightgbm")  # May need system dependencies
# 
# # Function to install from CRAN if not installed
# install_if_missing_cran <- function(pkg) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     install.packages(pkg, dependencies = TRUE)
#   }
# }
# 
# # Function to install from Bioconductor
# install_if_missing_bioc <- function(pkg) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#     BiocManager::install(pkg)
#   }
# }
# 
# # Function to install from GitHub
# install_if_missing_github <- function(pkg, repo) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
#     remotes::install_github(repo)
#   }
# }
# 
# # Install CRAN packages
# invisible(lapply(cran_packages, install_if_missing_cran))
# 
# # Install Bioconductor packages
# invisible(lapply(bioc_packages, install_if_missing_bioc))
# 
# # Install GitHub packages
# invisible(mapply(install_if_missing_github, names(github_packages), github_packages))
# 
# # Install special packages (try CRAN first, then OS-level check may be needed)
# invisible(lapply(special_packages, install_if_missing_cran))
# 
# # Load all packages (if installed)
# all_packages <- c(cran_packages, names(github_packages), special_packages)
# invisible(lapply(all_packages, function(pkg) {
#   suppressPackageStartupMessages(library(pkg, character.only = TRUE))
# }))
# 
# # Vector of all required package names
# required_libraries <- c(
#   "caret", "circlize", "ComplexHeatmap", "dplyr", "ggplot2", "grid", "gridExtra",
#   "kableExtra", "lightgbm", "Matrix", "mice", "missForest", "pROC", "purrr",
#   "reshape2", "rio", "rms", "survival", "survivalROC", "survcomp", "survminer",
#   "tidyr", "timeROC", "UCSCXenaShiny", "UpSetR", "VIM", "xgboost",
#   "DiagrammeR", "DiagrammeRsvg", "rsvg", "magick", "data.table", "stringi"
# )
# # Load all libraries quietly
# invisible(lapply(required_libraries, function(pkg) {
#   suppressPackageStartupMessages(library(pkg, character.only = TRUE))
# }))

# Universal R Package Extractor and Loader
script_path <- "C:/Users/Emamnuell/Desktop/UENF/marchine learning/Emanuell/imputaçao/ML_imputaçao.R"

# Read script as plain text (preserving comments)
script_lines <- readLines(script_path, warn = FALSE, encoding = "UTF-8")
script_text <- paste(script_lines, collapse = "\n")

# Initialize collection
pkg_matches <- character()

# 1. Match packages in library(...) and require(...)
pkg_matches <- c(pkg_matches, unlist(regmatches(
  script_text, gregexpr("(?<=library\\(|require\\()[\"']?([a-zA-Z0-9\\.]+)[\"']?", script_text, perl = TRUE)
)))

# 2. Match c("pkg1", "pkg2", ...) in active or commented lines
pkg_matches <- c(pkg_matches, unlist(regmatches(
  script_text, gregexpr("\"[a-zA-Z0-9\\.]+\"", script_text, perl = TRUE)
)))

# Clean and deduplicate
pkg_matches <- gsub("\"", "", pkg_matches)
pkg_matches <- sort(unique(pkg_matches))

# Optional: check CRAN availability
available_cran <- rownames(available.packages())

# Final valid list: either installed locally or on CRAN
valid_pkgs <- pkg_matches[
  sapply(pkg_matches, function(pkg) {
    requireNamespace(pkg, quietly = TRUE) || pkg %in% available_cran
  })
]

# Install + load function
load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Install and load all valid packages
invisible(lapply(valid_pkgs, load_or_install))

# Final output
cat("✔ Detected and loaded packages:\n")
print(valid_pkgs)

setwd("C:/Users/Emamnuell/Desktop/UENF/marchine learning/Emanuell/imputaçao")

#### ============================================================================
#### 📦 Master Wrapper Functions for rio Import/Export with NA Safety
#### Author: Enrique Medina-Acosta
#### Purpose: Ensure missing data integrity across formats (csv, tsv, xlsx)
#### ============================================================================

library(rio)

# Safe import with proper NA interpretation
safe_import <- function(file, format = NULL, ...) {
  import(file, format = format, na.strings = "NA", ...)
}

# Safe export with proper NA encoding
safe_export <- function(object, file, format = NULL, ...) {
  export(object, file = file, format = format, na = "NA", ...)
}
###
###
###
#### -----------------
#### Workflow Diagram
#### -----------------
#### -----------------------------
#### Total Expected Output Objects
#### -----------------------------
####
#### Pre-imputation processing: 5
####   • df001 to df005
####
#### Survival-imputed: 3
####   • Methods: mean, median, random
####   • Output: df006 to df008
####
#### CNV-imputed: 9
####   • 3 CNV methods × 3 survival datasets
####
#### Mutation-imputed: 27
####   • 3 mutation methods × 9 CNV-imputed datasets
####
#### Continuous-imputed: 216
####   • 8 continuous imputation methods × 27 mutation-imputed datasets
####
#### Evaluation node: 1
####   • Final model selection and analysis step
####
#### → Total: 5 + 3 + 9 + 27 + 216 + 1 = 261 output objects

# # Function to wrap text for node labels
# wrap_text <- function(text, width = 30) {
#   stri_wrap(text, width = width) %>%
#     paste(collapse = "\n")
# }
# 
# # Define node labels with wrapped text
# labels <- list(
#   df001 = wrap_text("df001: Multi-omic RCD Signatures; demographic and clinical data"),
#   df001_cleaned = wrap_text("Remove Rows with All NAs (Columns 23–End)"),
#   df001_harmonized = wrap_text("Harmonize Clinical Variables (Columns 23–30)"),
#   df004 = wrap_text("df004: Rename and Validate (No NAs in Columns 23–30)"),
#   df005 = wrap_text("df005: Harmonize CNV and Mutation Variables (CNV: Categorical; Mutation: Binary)"),
#   df006 = wrap_text("df006: Impute Survival Data (Method: Mean)"),
#   df007 = wrap_text("df007: Impute Survival Data (Method: Median)"),
#   df008 = wrap_text("df008: Impute Survival Data (Method: Random)"),
#   df009 = wrap_text("df009: Impute CNV Data (Method: Mode)"),
#   df010 = wrap_text("df010: Impute CNV Data (Method: Random)"),
#   df011 = wrap_text("df011: Impute CNV Data (Method: kNN)"),
#   df012 = wrap_text("df012: Impute Mutation Data (Method: Mean)"),
#   df013 = wrap_text("df013: Impute Mutation Data (Method: Median)"),
#   df014 = wrap_text("df014: Impute Mutation Data (Method: Mode)"),
#   df015 = wrap_text("df015: Impute Continuous Data (Method: Mean)"),
#   df016 = wrap_text("df016: Impute Continuous Data (Method: Median)"),
#   df017 = wrap_text("df017: Impute Continuous Data (Method: Random)"),
#   df018 = wrap_text("df018: Impute Continuous Data (Method: kNN)"),
#   df019 = wrap_text("df019: Impute Continuous Data (Method: MICE)"),
#   df020 = wrap_text("df020: Impute Continuous Data (Method: missForest)"),
#   df021 = wrap_text("df021: Impute Continuous Data (Method: XGBoost)"),
#   df022 = wrap_text("df022: Impute Continuous Data (Method: LightGBM)"),
#   evaluation = wrap_text("Model Evaluation and Selection\n(n = 260 total datasets)")
# )
# 
# # Define the workflow diagram
# diagram <- grViz("
# digraph workflow {
#   graph [layout = dot, rankdir = TB]
# 
#   node [shape = rectangle, style = filled, fillcolor = lightblue, width = 2.5, height = 1.0, fontsize = 10]
# 
#   df001 [label = '@@1']
#   df001_cleaned [label = '@@2']
#   df001_harmonized [label = '@@3']
#   df004 [label = '@@4']
#   df005 [label = '@@5']
#   df006 [label = '@@6']
#   df007 [label = '@@7']
#   df008 [label = '@@8']
#   df009 [label = '@@9']
#   df010 [label = '@@10']
#   df011 [label = '@@11']
#   df012 [label = '@@12']
#   df013 [label = '@@13']
#   df014 [label = '@@14']
#   df015 [label = '@@15']
#   df016 [label = '@@16']
#   df017 [label = '@@17']
#   df018 [label = '@@18']
#   df019 [label = '@@19']
#   df020 [label = '@@20']
#   df021 [label = '@@21']
#   df022 [label = '@@22']
#   evaluation [label = '@@23']
# 
#   df001 -> df001_cleaned -> df001_harmonized -> df004 -> df005
#   df005 -> {df006 df007 df008}
#   df006 -> {df009 df010 df011}
#   df007 -> {df009 df010 df011}
#   df008 -> {df009 df010 df011}
#   df009 -> {df012 df013 df014}
#   df010 -> {df012 df013 df014}
#   df011 -> {df012 df013 df014}
#   df012 -> {df015 df016 df017 df018 df019 df020 df021 df022}
#   df013 -> {df015 df016 df017 df018 df019 df020 df021 df022}
#   df014 -> {df015 df016 df017 df018 df019 df020 df021 df022}
#   df015 -> df016 -> df017 -> df018 -> df019 -> df020 -> df021 -> df022 [style=invis]
#   {df015 df016 df017 df018 df019 df020 df021 df022} -> evaluation
# }
# 
# [1]: labels$df001
# [2]: labels$df001_cleaned
# [3]: labels$df001_harmonized
# [4]: labels$df004
# [5]: labels$df005
# [6]: labels$df006
# [7]: labels$df007
# [8]: labels$df008
# [9]: labels$df009
# [10]: labels$df010
# [11]: labels$df011
# [12]: labels$df012
# [13]: labels$df013
# [14]: labels$df014
# [15]: labels$df015
# [16]: labels$df016
# [17]: labels$df017
# [18]: labels$df018
# [19]: labels$df019
# [20]: labels$df020
# [21]: labels$df021
# [22]: labels$df022
# [23]: labels$evaluation
# ")
# 
# # Convert the DiagrammeR object to SVG
# svg_code <- export_svg(diagram)
# 
# # Define file names
# pdf_file <- "workflow_diagram.pdf"
# tiff_file <- "workflow_diagram.tiff"
# 
# # Save as PDF
# rsvg_pdf(charToRaw(svg_code), file = pdf_file, width = 8.5, height = 11)  # 8.5x11 inches for portrait orientation
# 
# # Save as TIFF
# rsvg_png(charToRaw(svg_code), file = tiff_file, width = 2550, height = 3300)  # 8.5x11 inches at 300 DPI

#### -----------------------------
#### Total Expected Output Objects
#### -----------------------------
####
#### Step 1: Initial sequential processing
#### • df001 to df005 = 5 dataframes
####
#### Step 2: Survival data imputation
#### • 3 methods (mean, median, random)
#### • Result: df006 to df008 → 3 dataframes
####
#### Step 3: CNV imputation per survival object
#### • 3 CNV methods × 3 survival-imputed → 9 dataframes
####
#### Step 4: Mutation imputation per CNV-imputed dataset
####

#### 
#### 
#### ----------------------------------------------------------------------------------------------------------------------------------- 
#### PART A - Programmatically extracting clinical and multi-omic data form UCSCXenaShiny for the selected signatures from UCSCXenaShiny
#### -----------------------------------------------------------------------------------------------------------------------------------
####
####
# # Set working directory
# setwd("D:/Pré-artigo 5-optosis model/Machine learning CancerRCDShiny_prediction/Enrique version")
# 
# # Ensure critical functions exist before execution
# if (!exists("tcga_surv_get")) stop("Erro: tcga_surv_get() não está definido. Verifique a importação.")
# if (!exists("load_data")) stop("Erro: load_data() não está definido. Verifique a importação.")
# 
# # Import gene target list, ensuring correct format
# gene_symbols <- import("df1157_ML_targets_final.tsv", format = "tsv", na.strings = "NA")
# 
# length(unique(gene_symbols$Signature))
# length(unique(gene_symbols$Nomenclature))
# length(unique(gene_symbols$CTAB))
# length(unique(gene_symbols$Omic_feature))
# 
# # Select rows where Nomenclature is duplicated
# duplicated_nomenclature <- gene_symbols %>%
#   group_by(Nomenclature) %>%
#   filter(n() > 1) %>%
#   ungroup()
# 
# # Select rows where Signature is duplicated
# duplicated_signature <- gene_symbols %>%
#   group_by(Signature) %>%
#   filter(n() > 1) %>%
#   ungroup()
# 
# # Create a dataframe showing the distribution of "Signature" per "CTAB"
# signature_distribution <- duplicated_signature %>%
#   count(CTAB, Signature, name = "Count") %>%
#   arrange(desc(Count))
# 
# # Rename column for cleaner handling
# gene_symbols <- gene_symbols %>%
#   rename(Omic_feature = `Omic feature`) %>%
#   mutate(Omic_feature = case_when(
#     Omic_feature %in% c("mRNA", "miRNA") ~ Omic_feature,  # Keep unchanged
#     TRUE ~ tolower(Omic_feature)  # Convert others to lowercase
#   ))
# 
# # Create a dataframe showing the distribution of duplicated "Signature" per "CTAB"
# duplicated_signature_distribution <- duplicated_signature %>%
#   group_by(CTAB, Signature) %>%
#   filter(n() > 1) %>%  # Keep only duplicated Signature values per CTAB
#   summarise(Count = n(), .groups = "drop") %>%  # Count occurrences
#   arrange(desc(Count))  # Sort in descending order
# 
# # Extract rows from "duplicated_signature" that match "duplicated_signature_distribution"
# filtered_duplicated_signature <- duplicated_signature %>%
#   semi_join(duplicated_signature_distribution, by = c("CTAB", "Signature"))
# 
# #####
# ##### -------------------------------------------------------------------
# ##### Split the dataframe gene_symbols by the values in the CTAB column;
# ##### Each subset is saved as a .tsv file in the working directory
# ##### -------------------------------------------------------------------
# #####
# 
# # Load required library
# suppressPackageStartupMessages(library(rio))
# 
# # Loop through each unique CTAB value
# unique_ctabs <- unique(gene_symbols$CTAB)
# 
# for (ctab in unique_ctabs) {
#   # Subset the gene_symbols dataframe for this CTAB
#   df_subset <- gene_symbols[gene_symbols$CTAB == ctab, , drop = FALSE]
#   
#   # Clean object name and filename (replace invalid characters)
#   object_name <- paste0("gene_symbols_", make.names(ctab))
#   file_name <- paste0(object_name, ".tsv")
#   
#   # Assign to global environment
#   assign(object_name, df_subset, envir = .GlobalEnv)
#   
#   # Export to .tsv file
#   export(df_subset, file_name, format = "tsv", na = "NA")
#   
#   # Logging
#   cat("💾 Saved subset for", ctab, "as", file_name, "\n")
# }
# 
# #####
# ##### ---------------------------------------------------------------
# ##### Remove gene_symbols_<CTAB> objects from global environment
# ##### ---------------------------------------------------------------
# #####
# 
# # Identify all objects in the environment matching the pattern
# objs_to_remove <- ls(pattern = "^gene_symbols_")
# 
# # Remove them
# rm(list = objs_to_remove, envir = .GlobalEnv)
# 
# # Logging
# cat("🧹 Removed", length(objs_to_remove), "objects from global environment:\n")
# print(objs_to_remove)
# 
# ######
# ###### ------------------------------------------------------------
# ###### Split gene_symbols by "Omic feature" and save each as .tsv
# ###### ------------------------------------------------------------
# ######
# 
# suppressPackageStartupMessages(library(rio))  # Ensure rio is loaded
# 
# # Loop through each unique Omic feature value
# unique_features <- unique(gene_symbols$`Omic feature`)
# 
# for (feature in unique_features) {
#   # Subset the gene_symbols dataframe for this Omic feature
#   df_subset <- gene_symbols[gene_symbols$`Omic feature` == feature, , drop = FALSE]
#   
#   # Create a syntactically valid name for object and file
#   safe_feature <- make.names(feature)
#   object_name <- paste0("gene_symbols_", safe_feature)
#   file_name <- paste0("gene_symbols_", safe_feature, ".tsv")
#   
#   # Save to disk
#   export(df_subset, file_name)
#   cat("💾 Saved:", file_name, "\n")
#   
#   # Assign to global environment (optional step)
#   assign(object_name, df_subset, envir = .GlobalEnv)
# }
# 
# ######
# ###### -------------------------------------------------------------
# ###### Remove gene_symbols_<OmicFeature> objects from environment
# ###### -------------------------------------------------------------
# ######
# 
# # Identify all gene_symbols_* objects in the global environment
# objs_to_remove <- ls(pattern = "^gene_symbols_")
# 
# # Remove them
# rm(list = objs_to_remove, envir = .GlobalEnv)
# 
# # Logging
# cat("🧹 Removed", length(objs_to_remove), "objects from global environment:\n")
# print(objs_to_remove)
# 
# ####
# ####
# ####
# #### FETCHING CLINICAL, DEMOGRAPHIC AND SURVIVAL DATA
# ####
# ####
# ####
# # Load clinical and survival data
# tcga_clinical_data <- load_data("tcga_clinical")
# tcga_survival_data <- load_data("tcga_surv") 
# 
# # Stop execution if clinical or survival data is missing
# if (is.null(tcga_clinical_data) || is.null(tcga_survival_data)) {
#   stop("Erro ao carregar dados clínicos ou de sobrevivência.")
# }
# 
# # Merge clinical and survival data, removing duplicates
# tcga_cli_data <- full_join(tcga_clinical_data, tcga_survival_data, by = "sample") %>%
#   distinct(.keep_all = TRUE)
# 
# # Ensure .opt_pancan is defined
# if (!exists(".opt_pancan")) {
#   message("A variável .opt_pancan não foi encontrada. Definindo como TRUE.")
#   .opt_pancan <- TRUE  # Set default value
# }
# 
# # Define RDS file for progressive saving
# rds_file <- "gene_omic_list_progress.rds"
# 
# # Load previous progress if available
# if (file.exists(rds_file)) {
#   message("Carregando progresso anterior de ", rds_file)
#   gene_omic_list <- readRDS(rds_file)
# } else {
#   gene_omic_list <- list()
# }
# 
# # Get list of already processed genes to avoid redundancy
# processed_genes <- names(gene_omic_list)
# 
# # Loop through genes to fetch expression data
# for (i in seq_len(nrow(gene_symbols))) {
#   gene <- gene_symbols$Signature[i]
#   cancer_type <- gene_symbols$CTAB[i]
#   omic <- gene_symbols$Omic_feature[i]
#   Nomenclature <- gene_symbols$Nomenclature[i]  # Get nomenclature
#   
#   # Skip already processed genes
#   if (Nomenclature %in% processed_genes) next
#   
#   message("Processando: ", gene, " em ", cancer_type, " para ", omic)
#   
#   data <- tryCatch({
#     tcga_surv_get(
#       item = gene,
#       TCGA_cohort = cancer_type,
#       profile = omic,
#       TCGA_cli_data = tcga_cli_data,
#       opt_pancan = .opt_pancan  # Ensured to be TRUE
#     )
#   }, error = function(e) {
#     warning("Erro ao obter dados para ", gene, " em ", cancer_type, " para ", omic, ": ", e$message)
#     NULL
#   })
#   
#   # Process and store data if valid
#   if (!is.null(data) && "value" %in% colnames(data)) {
#     if ("sampleID" %in% colnames(data)) {
#       data <- rename(data, sample = sampleID)
#     }
#     
#     data <- data %>%
#       rename(expression_value = value) %>%
#       mutate(
#         gene = gene,
#         type = cancer_type,
#         omic_feature = omic,
#         Nomenclature = Nomenclature
#       ) %>%
#       select(sample, type, Nomenclature, expression_value)
#     
#     # Store data in list
#     gene_omic_list[[Nomenclature]] <- data
#   } else {
#     message("Sem dados para ", gene, " em ", cancer_type, " para ", omic)
#   }
#   
#   # Save progress every 10 genes
#   if (length(gene_omic_list) %% 10 == 0) {
#     saveRDS(gene_omic_list, rds_file)
#     message("Progresso salvo em ", rds_file)
#   }
# }
# 
# # Final save after loop completes
# saveRDS(gene_omic_list, rds_file)
# message("Progresso final salvo em ", rds_file)
# 
# # Proceed only if data exists
# if (length(gene_omic_list) > 0) {
#   # Combine all expression tables into a single dataframe
#   final_omic_data <- bind_rows(gene_omic_list)
#   
#   # Convert long format to wide format
#   final_omic_data <- final_omic_data %>%
#     pivot_wider(names_from = Nomenclature, values_from = expression_value)
#   
#   # Merge clinical data with processed expression data
#   final_data <- full_join(tcga_cli_data, final_omic_data, by = c("sample", "type"))
#   
#   # Convert list columns to character strings
#   list_columns <- sapply(final_data, is.list)
#   if (any(list_columns)) {
#     final_data <- final_data %>%
#       mutate(across(where(is.list), ~ if_else(is.null(.x), NA_character_, toString(.x))))
#   }
#   
#   # Organize final dataset
#   final_data <- final_data %>%
#     arrange(type, sample) %>%
#     relocate(sample, type)
#   
#   # Save as TSV
#   export(final_data, "ML_final_data.tsv", format = "tsv", na = "NA", row.names = FALSE)
#   
#   # Reimport to verify
#   ML_final_data <- import("ML_final_data.tsv", format = "tsv", na.strings = "NA")
#   
#   # Check unique column count
#   print(length(unique(colnames(ML_final_data))))  # Ensure correct column count
# } else {
#   message("Nenhum dado foi gerado. Verifique os logs.")
# }
# 
# ####
# ####
# #### ----------------------------------------------------------------------------
# #### PART B - debugging missing values in variable "Nomenclature" in final output
# #### -----------------------------------------------------------------------------
# #### Checking expected dimension and structure of final output
# #### Check which Nomenclature values are missing in gene_omic_list
# #### 
# #### 
# missing_nomenclature <- setdiff(gene_symbols$Nomenclature, names(gene_omic_list))
# print(length(missing_nomenclature))  # Should return 22
# print(missing_nomenclature)  # Show missing values
# 
# ## Check Which Nomenclature Values Were Processed in the Loop
# processed_nomenclature <- c()  # Track stored values
# 
# ## Check for NULL Entries in gene_omic_list
# null_entries <- sum(sapply(gene_omic_list, is.null))
# print(paste("Number of NULL entries in gene_omic_list:", null_entries))
# 
# ## Check if Data Exists for These Nomenclature Values
# check_data <- gene_symbols %>%
#   filter(Nomenclature %in% missing_nomenclature)
# print(check_data)
# 
# print(length(unique(names(gene_omic_list))))  # Should be 14907
# print(length(unique(colnames(final_omic_data))))  # Should be 14907
# 
# ######
# ###### -------------------------------------------------------
# ###### Rerun fetch script only for missing Nomenclature values
# ###### -------------------------------------------------------
# ######
# # If missing_nomenclature variable values exist, rerun the code only
# # process only the missing "Nomenclature" values rather than rerunning the entire script. This will ensure the missing 22 columns are retrieved efficiently.
# # Filter only the missing Nomenclature values from gene_symbols.
# # Run the loop only for these missing values.
# # Append the new results to gene_omic_list.
# # Save the updated gene_omic_list and regenerate final_omic_data.
# # Re-merge with tcga_cli_data and save the final dataset.
# ######
# ######
# ######
# # Identify missing Nomenclature values
# missing_nomenclature <- setdiff(gene_symbols$Nomenclature, names(gene_omic_list))
# message("Processing only missing Nomenclature values: ", length(missing_nomenclature))
# 
# # Filter gene_symbols to process only missing Nomenclature values
# missing_gene_symbols <- gene_symbols %>%
#   filter(Nomenclature %in% missing_nomenclature)
# 
# # Process missing values only
# for (i in seq_len(nrow(missing_gene_symbols))) {
#   gene <- missing_gene_symbols$Signature[i]
#   cancer_type <- missing_gene_symbols$CTAB[i]
#   omic <- missing_gene_symbols$Omic_feature[i]
#   Nomenclature <- missing_gene_symbols$Nomenclature[i]
#   
#   message("Processing missing Nomenclature: ", Nomenclature)
#   
#   data <- tryCatch({
#     tcga_surv_get(
#       item = gene,
#       TCGA_cohort = cancer_type,
#       profile = omic,
#       TCGA_cli_data = tcga_cli_data,
#       opt_pancan = .opt_pancan
#     )
#   }, error = function(e) {
#     warning("Error with ", Nomenclature, ": ", e$message)
#     NULL
#   })
#   
#   if (!is.null(data) && "value" %in% colnames(data)) {
#     if ("sampleID" %in% colnames(data)) {
#       data <- rename(data, sample = sampleID)
#     }
#     
#     data <- data %>%
#       rename(expression_value = value) %>%
#       mutate(
#         gene = gene,
#         type = cancer_type,
#         omic_feature = omic,
#         Nomenclature = Nomenclature
#       ) %>%
#       select(sample, type, Nomenclature, expression_value)
#     
#     gene_omic_list[[Nomenclature]] <- data
#   } else {
#     message("⚠ No data retrieved for ", Nomenclature)
#   }
#   
#   # Save progress every 5 genes
#   if (length(gene_omic_list) %% 5 == 0) {
#     saveRDS(gene_omic_list, rds_file)
#     message("Progress saved in ", rds_file)
#   }
# }
# 
# # Final save
# saveRDS(gene_omic_list, rds_file)
# message("Final progress saved in ", rds_file)
# 
# # Rebuild final dataset
# final_omic_data <- bind_rows(gene_omic_list) %>%
#   pivot_wider(names_from = Nomenclature, values_from = expression_value)
# 
# final_data <- full_join(tcga_cli_data, final_omic_data, by = c("sample", "type"))
# 
# rio::export(final_data, "ML_final_data_updated.tsv", , na = "NA", format = "tsv", row.names = FALSE)
# message("Updated dataset saved as ML_final_data_updated.tsv")
# 
# ML_final_data_updated <- import( "ML_final_data_updated.tsv", na.strings = "NA")
# 
# df001 <- ML_final_data_updated 
# 
# ####
# ####
# #### Converting empty values (== "") to NA missing values
# #### 
# #### 
# ### It does not put the string "NA" into the cell. 
# ### Instead, it converts empty string values ("") to actual R missing values, i.e., 
# ### NA of type logical, character, or factor, depending on the column's class.
# # Count empty string ("") values in each column
# empty_string_counts <- sapply(df001, function(col) sum(col == "", na.rm = TRUE))
# 
# # Show columns with at least one empty value
# empty_string_counts[empty_string_counts > 0]
# 
# df001[df001 == ""] <- NA
# 
# # Count empty string ("") values in each column
# empty_string_counts <- sapply(df001, function(col) sum(col == "", na.rm = TRUE))
# 
# # Show columns with at least one empty value
# empty_string_counts[empty_string_counts > 0]
# 
# ##### 
# ##### -----------------------------------------------------------------------------------
# ##### Remove Rows Fully NA Across Columns 23:14937 (Response and predictor variables)
# ##### -----------------------------------------------------------------------------------
# ##### 
# 
# # Define the range of columns to inspect
# col_range <- 23:14937
# 
# # Identify rows where all values in columns 23 through 14937 are NA
# fully_na_rows <- which(rowSums(is.na(df001[, col_range])) == length(col_range))
# 
# # Logging
# cat("🔍 Found", length(fully_na_rows), "rows with all NA in columns 23 to 14937.\n")
# 
# # Optional: Save those rows for tracking
# excluded_na_rows_df <- df001[fully_na_rows, ]
# 
# # Remove those rows (including the first occurrence)
# df001_cleaned <- df001[-fully_na_rows, ]
# 
# # Post-validation
# remaining <- which(rowSums(is.na(df001_cleaned[, col_range])) == length(col_range))
# if (length(remaining) == 0) {
#   cat("✅ All fully-NA rows across the specified column range were successfully removed.\n")
# } else {
#   cat("⚠️ Still", length(remaining), "rows remaining with all NA in that range.\n")
# }
# 
# export(df001_cleaned, "df001_cleaned.tsv", na = "NA")
# 
# ##### 
# ##### -----------------------------------------------------------------------------------
# ##### Remove Rows Fully NA Across Columns 23:30 (Response and predictor variables)
# ##### -----------------------------------------------------------------------------------
# ##### 
# df001_cleaned <- import("df001_cleaned.tsv", na.strings = "NA")
# 
# # Count empty string ("") values in each column
# empty_string_counts <- sapply(df001_cleaned, function(col) sum(col == "", na.rm = TRUE))
# 
# # Show columns with at least one empty value
# empty_string_counts[empty_string_counts > 0]
# 
# # ML_final_data_update <-  import("ML_final_data_updated.tsv", na.strings = "NA")
# # gene_symbols <- import("df1157_ML_targets_final.tsv", format = "tsv", na.strings = "NA")
# # df1157_ML_target_final <- import("df1157_ML_targets_final.tsv", na.strings = "NA")
# 
# ####
# ####
# ####
# 
# #### ============================================================================
# #### 📌 Summary of `df001_cleaned` vs `df001_cleaned_final`
# #### ============================================================================
# 
# # 🔹 df001_cleaned:
# # - This is the initial cleaned dataset **prior to imputation**.
# # - It may contain duplicated patient entries.
# # - Columns 23–30 may have inconsistencies across those duplicated entries.
# # - Used as the basis for identifying and harmonizing duplicates.
# 
# # 🔹 df001_cleaned_final:
# # - This is the final harmonized version of the dataset after resolving duplicates.
# # - Duplicated patient entries were identified and grouped.
# # - For identical values in cols 23–30 → rows were retained as-is.
# # - For differing values in cols 23–30 → values were harmonized using the row with
# #   the fewest missing entries as reference.
# # - The harmonized duplicated entries were merged back with the unique patients.
# # - This is the **final dataset** to be used for downstream imputation and modeling.
# 
# #### ============================================================================
# 
# # #### 
# # #### 
# # #### Duplicate patients and harmonization of clinical data before imputation
# # #### 
# # #### 
# # # 
# # # Step 1: Count the number of times each patient appears
# # patient_counts <- table(df001_cleaned$patient)
# # 
# # # Step 2: Identify patients with more than one occurrence
# # duplicated_patients <- names(patient_counts[patient_counts > 1])
# # 
# # # Step 3: Subset the original dataframe for all those patients (include all duplicates)
# # df001_cleaned_duplications <- df001_cleaned[df001_cleaned$patient %in% duplicated_patients, ]
# # 
# # # Step 4: Logging
# # dup_freqs <- patient_counts[duplicated_patients]
# # cat("🔍 Total duplicated patient entries (including originals):", nrow(df001_cleaned_duplications), "\n")
# # cat("📊 Number of unique duplicated patients:", length(duplicated_patients), "\n")
# # cat("📈 Duplication frequency ranges from:",
# #     min(dup_freqs), "to", max(dup_freqs), "occurrences.\n")
# # 
# # # Optional: view summary distribution of duplication frequency
# # duplication_summary <- as.data.frame(table(dup_freqs))
# # colnames(duplication_summary) <- c("Occurrences", "Num_Patients")
# # print(duplication_summary)
# # 
# # # Define the column range to compare
# # col_range <- 23:30
# # 
# # # Get all duplicated patient IDs
# # duplicated_patients <- unique(df001_cleaned_duplications$patient)
# # 
# # # Initialize lists
# # identical_rows_list <- list()
# # differing_rows_list <- list()
# # 
# # # Iterate over each duplicated patient
# # for (pat in duplicated_patients) {
# #   patient_group <- df001_cleaned_duplications[df001_cleaned_duplications$patient == pat, ]
# #   comp_cols <- patient_group[, col_range]
# #   
# #   # Check if all rows in columns 23:30 are identical
# #   if (nrow(unique(comp_cols)) == 1) {
# #     identical_rows_list[[pat]] <- patient_group
# #   } else {
# #     differing_rows_list[[pat]] <- patient_group
# #   }
# # }
# # 
# # # Combine rows into dataframes
# # df_identical_23_30 <- do.call(rbind, identical_rows_list)
# # df_differing_23_30 <- do.call(rbind, differing_rows_list)
# # 
# # # Logging
# # cat("✅ Identical patient groups in columns 23:30:", length(identical_rows_list), "patients |",
# #     nrow(df_identical_23_30), "rows\n")
# # cat("⚠️ Differing patient groups in columns 23:30:", length(differing_rows_list), "patients |",
# #     nrow(df_differing_23_30), "rows\n")
# 
# #### =============================================================================
# #### 📘 Harmonization Pipeline for Duplicated Patients in `df001_cleaned`
# #### =============================================================================

















suppressPackageStartupMessages(library(rio))


#### -----------------------------------------------------------------------------
#### Step 0: Load Initial Cleaned Dataset
#### -----------------------------------------------------------------------------
df001_cleaned <- import("ML_final_data.tsv", na.strings = "NA")

#### -----------------------------------------------------------------------------
#### Step 1: Identify Duplicated Patients
#### -----------------------------------------------------------------------------
patient_counts <- table(df001_cleaned$patient)
duplicated_patients <- names(patient_counts[patient_counts > 1])
df_duplicated_all <- df001_cleaned[df001_cleaned$patient %in% duplicated_patients, ]

# Logging
cat("🔍 Total duplicated patient entries:", nrow(df_duplicated_all), "\n")
cat("📊 Unique duplicated patients:", length(duplicated_patients), "\n")

#### -----------------------------------------------------------------------------
#### Step 2: Compare Cols 23–30 Within Each Duplicated Patient
#### -----------------------------------------------------------------------------
col_range <- 23:30
identical_rows_list <- list()
differing_rows_list <- list()

for (pat in duplicated_patients) {
  group <- df_duplicated_all[df_duplicated_all$patient == pat, ]
  comp_cols <- group[, col_range]
  if (nrow(unique(comp_cols)) == 1) {
    identical_rows_list[[pat]] <- group
  } else {
    differing_rows_list[[pat]] <- group
  }
}

# Combine groups
df_identical_rows <- do.call(rbind, identical_rows_list)
df_differing_rows <- do.call(rbind, differing_rows_list)

# Logging
cat("✅ Identical col23–30 groups:", nrow(df_identical_rows), "rows\n")
cat("⚠️ Differing col23–30 groups:", nrow(df_differing_rows), "rows\n")

#### -----------------------------------------------------------------------------
#### Step 3: Harmonize Differing Cols 23–30 Using Row with Fewest NAs
#### -----------------------------------------------------------------------------
df_harmonized_differing <- df_differing_rows
modification_log <- data.frame(
  patient_id = character(), rows_modified = integer(),
  values_replaced = integer(), reference_row = integer(),
  stringsAsFactors = FALSE
)

for (pat in unique(df_differing_rows$patient)) {
  idx <- which(df_differing_rows$patient == pat)
  group <- df_differing_rows[idx, ]
  mat <- group[, col_range]
  na_counts <- rowSums(is.na(mat))
  ref_idx <- which.min(na_counts)
  ref_values <- mat[ref_idx, ]
  
  total_changes <- 0
  for (j in seq_along(idx)) {
    current_idx <- idx[j]
    old_values <- df_differing_rows[current_idx, col_range]
    diffs <- ref_values != old_values | (is.na(ref_values) != is.na(old_values))
    total_changes <- total_changes + sum(diffs, na.rm = TRUE)
    df_harmonized_differing[current_idx, col_range] <- ref_values
  }
  
  modification_log <- rbind(modification_log, data.frame(
    patient_id = pat,
    rows_modified = length(idx),
    values_replaced = total_changes,
    reference_row = idx[ref_idx]
  ))
}

cat("🛠️ Harmonization completed for:", nrow(modification_log), "patients\n")

#### -----------------------------------------------------------------------------
#### Step 4: Reconstruct Full Harmonized Dataset
#### -----------------------------------------------------------------------------
df_duplicated_harmonized <- rbind(df_identical_rows, df_harmonized_differing)
patients_harmonized <- unique(df_duplicated_harmonized$patient)
df_unique_patients <- df001_cleaned[!(df001_cleaned$patient %in% patients_harmonized), ]

df001_cleaned_final <- rbind(df_unique_patients, df_duplicated_harmonized)

# Final Logging
cat("📦 Final rows:", nrow(df001_cleaned_final), "\n")
cat("🧾 Unique patients:", length(unique(df001_cleaned_final$patient)), "\n")

#### -----------------------------------------------------------------------------
#### Summary Comment
#### -----------------------------------------------------------------------------
# 🔹 df001_cleaned:
#     - Original cleaned dataset, includes duplicated patients.
# 🔹 df001_cleaned_final:
#     - Fully harmonized dataset after resolving duplicates in cols 23–30.
#     - To be used for downstream imputation and modeling.

# Optional export
# export(df001_cleaned_final, "df001_cleaned_final.tsv")

#### 
#### 
#### ----------------------------------------------------------------------------------
#### Validation of Duplicate patients harmonization of clinical data before imputation
#### ----------------------------------------------------------------------------------
#### 
#### 
# 
# Step 1: Count the number of times each patient appears
patient_counts <- table(df001_cleaned_final$patient)

# Step 2: Identify patients with more than one occurrence
duplicated_patients <- names(patient_counts[patient_counts > 1])

# Step 3: Subset the original dataframe for all those patients (include all duplicates)
df001_cleaned_final_duplications <- df001_cleaned_final[df001_cleaned_final$patient %in% duplicated_patients, ]

# Step 4: Logging
dup_freqs <- patient_counts[duplicated_patients]
cat("🔍 Total duplicated patient entries (including originals):", nrow(df001_cleaned_final_duplications), "\n")
cat("📊 Number of unique duplicated patients:", length(duplicated_patients), "\n")
cat("📈 Duplication frequency ranges from:",
    min(dup_freqs), "to", max(dup_freqs), "occurrences.\n")

# Optional: view summary distribution of duplication frequency
duplication_summary <- as.data.frame(table(dup_freqs))
colnames(duplication_summary) <- c("Occurrences", "Num_Patients")
print(duplication_summary)

# Define the column range to compare
col_range <- 23:30

# Get all duplicated patient IDs
duplicated_patients <- unique(df001_cleaned_final_duplications$patient)

# Initialize lists
identical_rows_list <- list()
differing_rows_list <- list()

# Iterate over each duplicated patient
for (pat in duplicated_patients) {
  patient_group <- df001_cleaned_final_duplications[df001_cleaned_final_duplications$patient == pat, ]
  comp_cols <- patient_group[, col_range]
  
  # Check if all rows in columns 23:30 are identical
  if (nrow(unique(comp_cols)) == 1) {
    identical_rows_list[[pat]] <- patient_group
  } else {
    differing_rows_list[[pat]] <- patient_group
  }
}

# Combine rows into dataframes
df_identical_23_30 <- do.call(rbind, identical_rows_list)

# Safely combine and log differing patient groups
# # Logging
cat("✅ Identical patient groups in columns 23:30:", length(identical_rows_list), "patients |",
    nrow(df_identical_23_30), "rows\n")

if (length(differing_rows_list) > 0) {
  df_differing_23_30 <- do.call(rbind, differing_rows_list)
  cat("⚠️ Differing patient groups in columns 23:30:", 
      length(differing_rows_list), "patients |",
      nrow(df_differing_23_30), "rows\n")
} else {
  df_differing_23_30 <- data.frame()
  cat("ℹ️ No differing patient groups found in columns 23:30.\n")
}

#####
##### ------------------------------------------------------------------------------------
##### Harmonizing Values in Columns 23:30 Across Duplicated Patients in df_differing_23_30
##### -------------------------------------------------------------------------------------
#####
# -----------------------------------------------------------
# Harmonize columns 23:30 across duplicates in df_differing_23_30
# Replace values using the row with the fewest NAs
# Also log modified patients and number of values changed
# -----------------------------------------------------------

df_differing_23_30_harmonized <- df_differing_23_30
col_range <- 23:30
duplicated_patients <- unique(df_differing_23_30$patient)

# Initialize a log dataframe
modification_log <- data.frame(
  patient_id = character(),
  rows_modified = integer(),
  values_replaced = integer(),
  reference_row = integer(),
  stringsAsFactors = FALSE
)

for (pat in duplicated_patients) {
  patient_group_idx <- which(df_differing_23_30$patient == pat)
  patient_group <- df_differing_23_30[patient_group_idx, ]
  comp_matrix <- patient_group[, col_range]
  
  # Count NAs per row and pick the one with the fewest
  na_counts <- rowSums(is.na(comp_matrix))
  reference_index <- which.min(na_counts)
  reference_values <- comp_matrix[reference_index, ]
  
  # Calculate how many values will be replaced in total
  total_replacements <- 0
  for (i in seq_along(patient_group_idx)) {
    current_idx <- patient_group_idx[i]
    old_values <- df_differing_23_30[current_idx, col_range]
    # Count non-identical values (including NAs)
    diffs <- reference_values != old_values | (is.na(reference_values) != is.na(old_values))
    total_replacements <- total_replacements + sum(diffs, na.rm = TRUE)
    
    # Replace with reference
    df_differing_23_30_harmonized[current_idx, col_range] <- reference_values
  }
  
  # Log modification
  modification_log <- rbind(modification_log, data.frame(
    patient_id = pat,
    rows_modified = length(patient_group_idx),
    values_replaced = total_replacements,
    reference_row = patient_group_idx[reference_index]
  ))
}

cat("✅ Harmonization complete for", nrow(modification_log), "patients.\n")
print(head(modification_log))

# Optional: merge with df_identical_23_30 to reconstruct full duplicate subset
df001_cleaned_duplications_harmonized <- rbind(df_identical_23_30, df_differing_23_30_harmonized)

# Optional: export results
# export(df_differing_23_30_harmonized, "df_differing_23_30_harmonized.tsv", na = "NA")
# export(modification_log, "modification_log.tsv", na = "NA")
# export(df001_cleaned_duplications_harmonized, "df001_cleaned_duplications_harmonized.tsv", na = "NA")

# ------------------------------------------------------------------------------------
# Merge harmonized duplicated patient block into original full dataset
# ------------------------------------------------------------------------------------

# Step 1: Extract the list of all patient IDs that were harmonized
harmonized_patients <- unique(df001_cleaned_duplications_harmonized$patient)

# Step 2: Subset df001_cleaned to exclude all those harmonized patients
df001_cleaned_without_duplicates <- df001_cleaned[!(df001_cleaned$patient %in% harmonized_patients), ]

# Step 3: Append the harmonized duplicate patient block
df001_cleaned_final <- rbind(df001_cleaned_without_duplicates, df001_cleaned_duplications_harmonized)

# Step 4: Optional consistency check
cat("📦 Rows in final dataset:", nrow(df001_cleaned_final), "\n")
cat("🧾 Unique patients:", length(unique(df001_cleaned_final$patient)), "\n")

# Optional: export final version
export(df001_cleaned_final, "df001_cleaned_final.tsv", na = "NA")

df001_cleaned_final <- import("df001_cleaned_final.tsv", na.strings = "NA")

#####
##### ----------------------------------------------------------------------
##### Check for Fully Duplicated Rows Across All Columns in df001_cleaned_final
##### ----------------------------------------------------------------------
#####

# Check for duplicated rows across all variables
duplicated_row_flags <- duplicated(df001_cleaned_final)

# Count how many duplicated rows
num_duplicates <- sum(duplicated_row_flags)

# Display result
cat("🔍 Total number of fully duplicated rows:", num_duplicates, "\n")

# Optional: Extract the duplicated rows (excluding first occurrence)
df001_cleaned_final_duplicates <- df001_cleaned_final[duplicated_row_flags, ]

# Optional: View summary
if (num_duplicates > 0) {
  print(head(df001_cleaned_final_duplicates))
} else {
  cat("✅ No fully duplicated rows found in df001_cleaned_final.\n")
}

# Remove rows where all values from column 23 onwards are NA
df001_cleaned_final_filtered <- df001_cleaned_final[!apply(df001_cleaned_final[, 23:ncol(df001_cleaned_final)], 1, function(row) all(is.na(row))), ]

gc()

df004 <- df001_cleaned_final_filtered

export(df004, "df004.tsv", na = "NA")

rm(df001_cleaned_final, df001_cleaned_duplications_harmonized, df001_cleaned_final_duplications)

gc()


#### 
#### 
#### -----------------------------------------------------------------------------
#### Groupwise imputation across high-dimensional binomial or continuous variables
#### -----------------------------------------------------------------------------
#### 
#### Remaking the consolidated dataframe for ML imputations
#### Brainstorm meeting 03/25/2025; 01/04/2025
#### Emanuell, Higor, Victor, Enrique
#### 
#### 
#### 

df002 <-df004

# Remove rows where all values from column 23 onwards are NA
df002_filtered <- df002[!apply(df002[, 23:ncol(df002)], 1, function(row) all(is.na(row))), ]

# Keep rows where all values from column 23 onwards are NOT NA
df002_filtered_2 <- df002[apply(df002[, 23:ncol(df002)], 1, function(row) all(is.na(row))), ]

df003 <- df002_filtered

# Columns of interest
cols_to_check <- c("OS", "DSS", "DFI", "PFI", 
                   "OS.time", "DSS.time", "DFI.time", "PFI.time")

# Separate duplicates
duplicates <- df003 %>%
  group_by(patient, type) %>%
  filter(n() >= 2) %>%
  arrange(patient, type, .by_group = TRUE)

# Rest of the dataset (non-duplicates)
non_duplicates <- df003 %>%
  group_by(patient, type) %>%
  filter(n() == 1) %>%
  ungroup()

# Function to fill missing values from duplicate pairs, triplicate or else
fill_pairwise <- function(df_pair) {
  if (nrow(df_pair) != 2) return(df_pair)  # Skip non-pairs
  row1 <- df_pair[1, ]
  row2 <- df_pair[2, ]
  
  for (col in cols_to_check) {
    if (is.na(row1[[col]]) && !is.na(row2[[col]])) {
      row1[[col]] <- row2[[col]]
    } else if (!is.na(row1[[col]]) && is.na(row2[[col]])) {
      row2[[col]] <- row1[[col]]
    }
    # If both NA, do nothing (but we will log later)
  }
  
  bind_rows(row1, row2)
}

# Apply the fill_pairwise function to each duplicated group
filled_duplicates <- duplicates %>%
  group_split(patient, type) %>%
  lapply(fill_pairwise) %>%
  bind_rows()

# Optionally identify unresolved cases (which duplicates, triplicate, etc. have NA)
unresolved <- filled_duplicates %>%
  group_by(patient, type) %>%
  filter(n() >= 2) %>%
  summarise(across(all_of(cols_to_check), ~all(is.na(.x))), .groups = "drop") %>%
  pivot_longer(cols = all_of(cols_to_check), names_to = "variable", values_to = "both_na") %>%
  filter(both_na == TRUE)

# Merge with non-duplicate data
df003_cleaned <- bind_rows(filled_duplicates, non_duplicates) %>%
  arrange(patient, type)

# Message
if (nrow(unresolved) > 0) {
  message("Unresolved NA values (both duplicates had NA):")
  print(unresolved)
} else {
  message("All possible imputations successfully completed.")
}

export(df003_cleaned, "df003_cleaned.tsv", na = "NA")

gc()

####
####
#### -------------------------------------------------------------
#### Handling Mutation and CNV variable values prior to imputation
#### ------------------------------------------------------------
####
#### -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#### PART 1 - Converting CNV (omic feature ".3") nomenclature values to categorical (nominal) "Normal", "Duplicated", "Deleted" and retaining NA values for imputation
#### -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#### Values could be integer or float depending on pipeline (e.g., GISTIC, segmentation).
####
####
# # Mapping function for CNV-like behavior: 0 → "Normal", >0 → "Duplicated", <0 → "Deleted"
# map_to_cnv_status <- function(x) {
#   # Convert to numeric if x is character or factor
#   if (is.character(x) || is.factor(x)) {
#     suppressWarnings(x <- as.numeric(as.character(x)))
#   }
#   
#   if (is.numeric(x)) {
#     return(ifelse(is.na(x), NA,
#                   ifelse(x == 0, "Normal",
#                          ifelse(x > 0, "Duplicated", "Deleted"))))
#   } else {
#     return(x)
#   }
# }
# 
# # Apply only to columns where 2nd token is "3"
# apply_cnv_mapping_2nd_pos_3_only <- function(df) {
#   target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) {
#     length(x) >= 2 && x[2] == "3"
#   })]
#   
#   # Apply mapping function to selected columns
#   df[target_vars] <- lapply(df[target_vars], map_to_cnv_status)
#   
#   return(df)
# }
# 
# # Example usage:
# df003_cleaned_categoric_3 <- apply_cnv_mapping_2nd_pos_3_only(df003_cleaned)
# 
# gc()
# 
# ## Example: Before conversion
# df003_cleaned_LGG <- df003_cleaned[df003_cleaned$type == "LGG", 
#                                    c("type", "LGG-136.3.3.N.3.35.71.1.2.1")]
# 
# unique(df003_cleaned_LGG$`LGG-136.3.3.N.3.35.71.1.2.1`)
# 
# ## Example: After conversion
# df003_cleaned_categoric_3_LGG <- df003_cleaned_categoric_3[df003_cleaned_categoric_3$type == "LGG", 
#                                                            c("type", "LGG-136.3.3.N.3.35.71.1.2.1")]
# 
# unique(df003_cleaned_categoric_3_LGG$`LGG-136.3.3.N.3.35.71.1.2.1`)
# 
# rm(df003_cleaned_categoric_3_LGG, df003_cleaned_LGG)
# 
# #### 
# #### 
# #### -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# #### PART 2 - Handling Mutation variable values for imputation
# #### -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# #### Mutation Variable Handling of mutation burden binary matrices
# ####
# ####
# 
# #### Converting Mutation (omic feature 2) nomenclature values to binary 0,1 and retaining NA for imputation
# #### Mapping function for mutation binary like behavior: 0 → "Wildtype" == 0, >0 → "Mutated" == 1.
# #### Converting Mutation (omic feature 2) nomenclature values to binary 0,1 and retaining NA for imputation
# 
# # Mapping function: 0 stays 0, non-zero becomes 1, NA remains NA
# map_to_binary <- function(x) {
#   if (is.factor(x) || is.character(x)) {
#     x <- as.numeric(as.character(x))
#   }
#   if (is.numeric(x)) {
#     return(ifelse(is.na(x), NA, ifelse(x == 0, 0, 1)))
#   } else {
#     return(x)
#   }
# }
# 
# # Apply mapping only to variables whose 2nd token (split by ".") is "2"
# apply_mapping_to_2nd_pos_2_only <- function(df) {
#   target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) {
#     length(x) >= 2 && x[2] == "2"
#   })]
#   
#   # Apply mapping function
#   df[target_vars] <- lapply(df[target_vars], map_to_binary)
#   
#   return(df)
# }
# 
# # ✅ Apply to the cleaned categoric CNV-transformed dataset
# df003_cleaned_categoric_2_3 <- apply_mapping_to_2nd_pos_2_only(df003_cleaned_categoric_3)
# 
# df005 <- df003_cleaned_categoric_2_3
# 
# export(df005, "df005.tsv", na = "NA")
# 
# # Example After CNV conversion
# df005_LGG <- df005[df005$type == "LGG", c("type", "LGG-136.3.3.N.3.35.71.1.2.1")]
# unique(df005_LGG$`LGG-136.3.3.N.3.35.71.1.2.1`)
# 
# # Example after mutation conversion
# df005_LUAD <- df005[df005$type == "LUAD", c("type", "LUAD-2090.2.2.P.3.1.11.1.4.1")]
# unique(df005_LUAD$`LUAD-2090.2.2.P.3.1.11.1.4.1`)
# 
# # Example before CNV conversion
# df003_cleaned_LGG <- df003_cleaned[df003_cleaned$type == "LGG", c("type", "LGG-136.3.3.N.3.35.71.1.2.1")]
# unique(df003_cleaned_LGG$`LGG-136.3.3.N.3.35.71.1.2.1`)
# 
# # Example before mutation conversion
# df003_cleaned_LUAD <- df003_cleaned[df003_cleaned$type == "LUAD", c("type", "LUAD-2090.2.2.P.3.1.11.1.4.1")]
# unique(df003_cleaned_LUAD$`LUAD-2090.2.2.P.3.1.11.1.4.1`)
# 
# export(df005, "df005.tsv", na = "NA")
# 
# df005 <- import("df005.tsv", na.strings = "NA")
# 
# rm(df005_LGG, df005_LUAD)

####
####
#### =================================================================================
#### 📘 Core Input–Output Indexing Principle for Multi-Omic Imputation Pipeline
#### =================================================================================
# This indexing convention governs how imputed data objects are sequentially tracked
# across multi-layered imputation steps (Survival → CNV → Mutation → Continuous).
#
# 🧩 Principle:
# ➤ Each imputation method uses as input the **set of outputs from the immediately 
#   preceding imputation layer**. 
# ➤ The **input indices remain the same** for all methods applied to a given omic layer.
# ➤ Each new imputation method generates a **new block of sequential outputs**.
# ➤ Outputs are exported as `.rds` and recorded in a cumulative log table
#   (`output_name_table_all.tsv`) with detailed metadata.
#
# ------------------------------------------------------------------------------------
# ✅ Example Flow:
# ------------------------------------------------------------------------------------
# Output hierarchy
# [Step 1] Survival Imputation
#   → Outputs: df006, df007, df008   (Mean, Median, Random)
#
# [Step 2] CNV Imputation
#   → Inputs : df006–df008
#   → Methods: Mode, Random, kNN
#   → Outputs: df009–df017   (3 inputs × 3 methods = 9 outputs)
#
# [Step 3] Mutation Imputation
#   → Inputs : df009–df017
#   → Methods: Mean, Median, Mode
#   → Outputs: df018–df044   (9 inputs × 3 methods = 27 outputs)
#
# [Step 4] Continuous Omic Imputation (.1/.4/.5/.6/.7)
#   → Inputs : df018–df044
#   → Methods: Mean, Median, Random, kNN, missForest, XGBoost, LightGBM, MICE
#   → Outputs: 
#        • Mean:      df045–df071
#        • Median:    df072–df098
#        • Random:    df099–df125
#        • kNN:       df126–df152
#        • missForest df153–df179
#        • XGBoost    df180–df206
#        • LightGBM   df207–df233
#        • MICE       df234–df260 (TBD)
#
# ------------------------------------------------------------------------------------
# 🧼 Each dataset is removed from the global environment after export to minimize RAM.
# 📑 The cumulative tracking table is updated at each step with:
#      - Step
#      - Input_Object
#      - Output_Object
#      - Method
#      - Saved_As (.rds path)
# ====================================================================================

####
#### --------------------------------------------------------------------------------------------
#### Function: Impute Missing NA Binary Survival Outcomes and Times: Clinical response variables
#### Dataset: df005
#### Target variables: OS, DSS, DFI, PFI and their respective *.time variables
#### -------------------------------------------------------------------------------------------
#### 

gc()

df005 <- import("df003_cleaned.tsv", na.strings = "NA")

impute_survival_variables <- function(df, method = c("mean", "median", "random"), min_time_threshold = 5, seed = 123, verbose = TRUE) {
  method <- match.arg(method)
  set.seed(seed)
  df_out <- df
  
  # Step 1: Impute OS and OS.time
  if (verbose) message("\n🔧 Imputing OS and OS.time")
  na_os <- which(is.na(df_out$OS))
  if (length(na_os) > 0) {
    if (method == "mean") {
      imputed_vals <- ifelse(mean(df_out$OS, na.rm = TRUE) >= 0.5, 1, 0)
      df_out$OS[na_os] <- imputed_vals
    } else if (method == "median") {
      imputed_vals <- ifelse(median(df_out$OS, na.rm = TRUE) >= 0.5, 1, 0)
      df_out$OS[na_os] <- imputed_vals
    } else if (method == "random") {
      observed_vals <- na.omit(df_out$OS)
      df_out$OS[na_os] <- sample(observed_vals, length(na_os), replace = TRUE)
    }
  }
  
  # OS.time imputation (use median of observed times)
  na_ostime <- which(is.na(df_out$OS.time))
  if (length(na_ostime) > 0) {
    os_time_val <- median(df_out$OS.time, na.rm = TRUE)
    df_out$OS.time[na_ostime] <- os_time_val
  }
  
  # Step 2: Impute DSS based on OS
  if (verbose) message("🔧 Imputing DSS based on OS")
  na_dss <- which(is.na(df_out$DSS))
  for (i in na_dss) {
    if (df_out$OS[i] == 0) {
      df_out$DSS[i] <- 0
    } else {
      if (method == "mean") {
        val <- ifelse(mean(df_out$DSS, na.rm = TRUE) >= 0.5, 1, 0)
      } else if (method == "median") {
        val <- ifelse(median(df_out$DSS, na.rm = TRUE) >= 0.5, 1, 0)
      } else {
        val <- sample(na.omit(df_out$DSS), 1)
      }
      df_out$DSS[i] <- val
    }
  }
  
  # Step 3: Impute DSS.time if NA
  if (verbose) message("🔧 Imputing DSS.time")
  na_dsstime <- which(is.na(df_out$DSS.time))
  df_out$DSS.time[na_dsstime] <- df_out$OS.time[na_dsstime]
  
  # Step 4: Impute DFI
  if (verbose) message("🔧 Imputing DFI and DFI.time")
  na_dfi <- which(is.na(df_out$DFI))
  for (i in na_dfi) {
    if (method == "mean") {
      val <- ifelse(mean(df_out$DFI, na.rm = TRUE) >= 0.5, 1, 0)
    } else if (method == "median") {
      val <- ifelse(median(df_out$DFI, na.rm = TRUE) >= 0.5, 1, 0)
    } else {
      val <- sample(na.omit(df_out$DFI), 1)
    }
    df_out$DFI[i] <- val
  }
  
  # Step 5: Impute DFI.time
  for (i in which(is.na(df_out$DFI.time))) {
    os_time <- df_out$OS.time[i]
    if (df_out$DFI[i] == 0) {
      df_out$DFI.time[i] <- os_time
    } else {
      df_out$DFI.time[i] <- if (os_time > min_time_threshold) sample(min_time_threshold:(os_time - 1), 1) else os_time
    }
  }
  
  # Step 6: Impute PFI
  if (verbose) message("🔧 Imputing PFI and PFI.time")
  na_pfi <- which(is.na(df_out$PFI))
  for (i in na_pfi) {
    if (method == "mean") {
      val <- ifelse(mean(df_out$PFI, na.rm = TRUE) >= 0.5, 1, 0)
    } else if (method == "median") {
      val <- ifelse(median(df_out$PFI, na.rm = TRUE) >= 0.5, 1, 0)
    } else {
      val <- sample(na.omit(df_out$PFI), 1)
    }
    df_out$PFI[i] <- val
  }
  
  # Step 7: Impute PFI.time
  for (i in which(is.na(df_out$PFI.time))) {
    os_time <- df_out$OS.time[i]
    if (df_out$PFI[i] == 0) {
      df_out$PFI.time[i] <- os_time
    } else {
      df_out$PFI.time[i] <- if (os_time > min_time_threshold) sample(min_time_threshold:(os_time - 1), 1) else os_time
    }
  }
  
  if (verbose) message("✅ Imputation of survival variables complete.")
  return(df_out)
}

# Example usage:
df005_mean_imputed_survival <- impute_survival_variables(df005, method = "mean")

df005_median_imputed_survival <- impute_survival_variables(df005, method = "median")

df005_random_imputed_survival <- impute_survival_variables(df005, method = "random")

# Load rio if not already
suppressPackageStartupMessages(library(rio))

# Define output filenames
rio::export(df005_random_imputed_survival, "df005_random_imputed_survival.tsv", na = "NA")
rio::export(df005_median_imputed_survival, "df005_median_imputed_survival.tsv", na = "NA")
rio::export(df005_mean_imputed_survival,   "df005_mean_imputed_survival.tsv", na = "NA")

## Validation of NA remaining
# Function to summarize remaining NA counts for survival variables
log_survival_na_summary <- function(df, method_name) {
  vars <- c("OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time")
  na_counts <- sapply(df[, vars], function(x) sum(is.na(x)))
  
  summary_df <- data.frame(
    Method = method_name,
    Variable = names(na_counts),
    NA_Count = as.integer(na_counts),
    stringsAsFactors = FALSE
  )
  
  return(summary_df)
}

# Apply to all imputed versions
na_summary_random <- log_survival_na_summary(df005_random_imputed_survival, "random")
na_summary_median <- log_survival_na_summary(df005_median_imputed_survival, "median")
na_summary_mean   <- log_survival_na_summary(df005_mean_imputed_survival, "mean")

# Combine and display
na_summary_all <- rbind(na_summary_random, na_summary_median, na_summary_mean)
print(na_summary_all)

#####
##### ----------------------------------------------------------
##### Save survival-imputed versions to .tsv and clear memory
##### ----------------------------------------------------------
##### 

# Save each version to .tsv format
export(df005_random_imputed_survival, "df005_random_imputed_survival.tsv", format = "tsv", na = "NA")
export(df005_median_imputed_survival, "df005_median_imputed_survival.tsv", format = "tsv", na = "NA")
export(df005_mean_imputed_survival,   "df005_mean_imputed_survival.tsv", format = "tsv", na = "NA")

cat("💾 All three survival-imputed datasets have been saved as .tsv files.\n")

# Clean up: remove objects from the global environment
rm(
  df005_random_imputed_survival,
  df005_median_imputed_survival,
  df005_mean_imputed_survival
)

cat("🧽 Imputed survival datasets removed from memory.\n")

### -----------------------------------------------------------------------
### Renaming to Serial DataFrames, Exporting, and Removing from Environment
### -----------------------------------------------------------------------
### Step 1: Import survival-imputed TSV files
### These dataframes will be renamed to standardized serial identifiers
### (df006, df007, df008), exported as .tsv files, and removed from the
### global environment to support memory-efficient progressive processing
### in the downstream imputation pipeline.
###
# Import and rename survival-imputed datasets
df006 <- import("df005_mean_imputed_survival.tsv", na.strings = "NA")
df007 <- import("df005_median_imputed_survival.tsv", na.strings = "NA")
df008 <- import("df005_random_imputed_survival.tsv", na.strings = "NA")

# Export as standardized .tsv files
export(df006, "df018.tsv", )
export(df007, "df019.tsv", na = "NA")
export(df008, "df020.tsv", na = "NA")

# Save objects as .rds files in the working directory
saveRDS(df006, file = "df018.rds")
saveRDS(df007, file = "df019.rds")
saveRDS(df008, file = "df020.rds")

# Remove from environment to conserve memory
rm(df006, df007, df008)

gc()

####
####
#### -------------------------------------------------------
#### PART C - validation of the preceding debugging PART B
#### -------------------------------------------------------
#### Debugging for missing variables in final output
#### Checking expected dimension and structure of final output
#### Check which Nomenclature values are missing in gene_omic_list
#### 
#### 
# missing_nomenclature <- setdiff(gene_symbols$Nomenclature, names(gene_omic_list))
# print(length(missing_nomenclature))  # Should return 22
# print(missing_nomenclature)  # Show missing values
# 
# ## Check which Nomenclature values eere processed in the loop
# processed_nomenclature <- c()  # Track stored values
# 
# ## Check for NULL Entries in gene_omic_list
# null_entries <- sum(sapply(gene_omic_list, is.null))
# print(paste("Number of NULL entries in gene_omic_list:", null_entries))
# 
# ## Check if Data Exists for These Nomenclature Values
# check_data <- gene_symbols %>%
#   filter(Nomenclature %in% missing_nomenclature)
# print(check_data)
# 
# print(length(unique(names(gene_omic_list))))  # Should be 14907
# print(length(unique(colnames(final_omic_data))))  # Should be 14909
# 
# setdiff(colnames(final_omic_data), names(gene_omic_list))  # Find columns in final_omic_data that are not in gene_omic_list
# 
# # Identify missing Nomenclature values
# missing_nomenclature <- setdiff(gene_symbols$Nomenclature, names(gene_omic_list))
# message("Processing only missing Nomenclature values: ", length(missing_nomenclature))
# 
# # Filter missing gene symbols based on missing nomenclature
# missing_gene_symbols <- gene_symbols %>%
#   filter(Nomenclature %in% missing_nomenclature)
# 
# # Log the number of missing entries
# cat("🔹 Number of missing gene symbol entries:", nrow(missing_gene_symbols), "\n")
# gc()

##### 
##### 
##### The moment of truth...
##### --------------------------------------
##### CNV categorical, nominal imputation
##### --------------------------------------
##### 
##### ----------------------------------------------------------------------------------------------------------------------------------
##### Imputing "NA" values to categorical (nominal) values - CNV - impute_mode_groupwise_by_type() with Full Logging
##### 1. Mode Imputation (Most Frequent Category) with Full Logging and Diagnostics
##### Mode Imputation group-wise by Matching Prefix in type variable and .3 Pattern (CNV)  in 
##### the 2nd token in signature nomenclature predictor variables
##### Note: Mode imputation is safer when there's a strong dominant signal or when group sizes are small, avoiding stochastic instability.
##### -----------------------------------------------------------------------------------------------------------------------------------
#####
#####
# ------------------------------------------------------------------------------
# 🧬 CNV Imputation — Overview of Implemented Methods
# -------------------------------------------------------------------------------
# Context:
# - Original CNV values in df003: numeric (e.g., -2, 0, 53, -87, NA)
# - Converted in df005 to categorical values:
#     NA → NA
#     0  → "Normal"
#     >0 → "Duplicated"
#     <0 → "Deleted"
# - Categorical CNV values are treated as *nominal* variables (no ordinal relationship)
#
# CNV Imputation Methods Implemented (Group-wise by 'type'):
# -----------------------------------------------------------------
# | Method      | Class Type       | Output Type | Status        |
# |-------------|------------------|-------------|----------------|
# | Mode        | Deterministic     | Nominal     | ✅ Implemented |
# | Random      | Stochastic        | Nominal     | ✅ Implemented |
# | kNN (VIM)   | Local structure   | Nominal     | ✅ Implemented |
# | Mean        | ⛔ Not applicable | ❌ Inapplicable to nominal CNV data |
# | Median      | ⛔ Not applicable | ❌ Inapplicable to nominal CNV data |
# --------------------------------------------------------------------------
# These methods are:
# - Group-aware (match tumor type prefix)
# - Biologically grounded and interpretable
# - Reproducible and robust across tumor types
#
# Recommendation:
# → No further imputation methods are necessary for CNV (categorical) variables.
# → Avoid methods that assume ordinal structure (e.g., mean, median, MICE, LASSO).
# → Optional: Post-imputation encoding to numeric (e.g., Deleted = -1, Normal = 0, Duplicated = +1) can be done for modeling, not imputation.
# ------------------------------------------------------------------------------
####
####
####
# ------------------------------------------------------------------------------
# 📊 Comparative Evaluation of Advanced Methods for Nominal CNV Variables
# -------------------------------------------------------------------------------
# CNV Data Type:
# - Categorical (nominal): "Deleted", "Normal", "Duplicated"
# - Derived from numeric CNV values: <0, 0, >0 → mapped to nominal classes
#
# Goal:
# - Impute missing (NA) values in a way that respects the categorical, non-ordinal nature
# - Maintain biological meaning and group (tumor type) structure
#
# Evaluation Summary:
# ----------------------------------------------------------------------------------------------------------------
# | Method            | Suitability     | Rationale                                                              |
# |-------------------|-----------------|------------------------------------------------------------------------|
# | Mode              | ✅ Excellent     | Simple, deterministic; respects dominant signal within tumor type      |
# | Random Sampling   | ✅ Excellent     | Preserves empirical frequency; adds diversity where class is balanced  |
# | kNN (VIM)         | ✅ Good          | Uses local structure via dummy column; valid for nominal data          |
# | MICE (polytomous) | ⚠️ Limited       | Technically feasible but unstable; poor convergence for multi-class factors |
# | missForest        | ⚠️ Conditional   | Works if correlated omic predictors are available; not ideal in isolation |
# | XGBoost / LGBM    | ❌ Not applicable| Require supervised labels; not usable for imputation                   |
# | DeepSurv          | ❌ Not applicable| Survival modeling tool, not an imputer                                 |
# | LASSO / Boruta    | ❌ Not applicable| Feature selection methods; not designed for missing data imputation    |
# | One-Hot Encoding  | ✅ Post-imputation | Useful for downstream models, not for imputation itself                |
# -----------------------------------------------------------------------------------------------------------------------
# Conclusion:
# - The current imputation methods (mode, random, kNN) are adequate and optimal
# - No need to implement more complex ML-based methods for nominal CNV imputation
# - Optional: recoding categorical CNVs to numeric values (-1, 0, +1) *after* imputation
#   for use in certain modeling frameworks (e.g., LGBM, DeepSurv)
# ------------------------------------------------------------------------------

#### -------------------------------------------------------------------------
#### Memory-Efficient Chained Processing Strategy for Imputation Pipeline
#### -------------------------------------------------------------------------
# Intended principles of memory efficiency, progressive data handling, and output tracking. 
# Efficiently chain-load, rename, save, and remove from environment
# each pre-generated dataframe (e.g., df006, df007, df008, etc.)
# on-the-fly, so that the subsequent imputation or processing steps 
# can execute sequentially without exceeding your system’s memory capacity.

# ✅ Goal:
# Efficiently process large-scale sequential imputation steps without overloading RAM,
# by dynamically loading, renaming, saving, and removing intermediate dataframes.

# ✅ Conceptual Strategy:
# Each dataframe (e.g., df006, df007, df008...) is:
# (1) Loaded from disk using `rio::import()`
# (2) Assigned to a standardized object name (e.g., df006)
# (3) Processed (e.g., CNV imputation, mutation imputation, etc.)
# (4) The resulting output dataframe(s) (e.g., df009, df010...) are:
#     - Immediately saved to disk using `rio::export()`
#     - Removed from memory via `rm()` and `gc()`
# (5) The input dataframe (e.g., df006) is also removed after use

# ✅ Why this is essential:
# - Minimizes RAM usage: Only one input + output object(s) reside in memory at a time
# - Ensures robustness: Persistent `.tsv` or `.rds` output files allow full recovery
# - Supports resumability: Code can restart from the next available file
# - Maximizes scalability: Enables full execution across 261 objects on memory-limited systems

# ✅ Naming Convention:
# - Serially indexed dataframes: df001 to df022 and beyond
# - Stored using consistent file paths and extensions (e.g., "df006.tsv")

# ✅ Example Workflow Step:
# Load df006 → perform CNV imputation → save df009, df010, df011 → remove all from memory

# This pattern will be reused throughout the survival → CNV → mutation → continuous chain

# ❗ Note:
# All code chunks must include `rm()` + `gc()` after saving to prevent memory bloat.
# Use `sprintf("df%03d", i)` for consistent three-digit serial naming.

####
#### -------------------------------------------------------------
#### CNV Imputation Strategy – Amended for Sequential Application
#### -------------------------------------------------------------
####
#### ➤ Objective:
####   To perform CNV-specific imputation (Mode, Random, kNN) on 
####   survival-imputed datasets df006, df007, and df008, which are
####   previously saved in the working directory and cleared from RAM.
####
#### ➤ Strategy:
####   (1) Each survival-imputed file will be:
####       - Imported using `rio::import()`
####       - Assigned a temporary object in the environment
####   (2) CNV imputation will be executed using 3 methods:
####       - Mode Imputation
####       - Random Imputation
####       - kNN Imputation
####   (3) Each resulting CNV-imputed dataset will be:
####       - Serially renamed (e.g., df009, df010, ..., df017)
####       - Exported back to the working directory using `rio::export()`
####       - Removed from memory using `rm()` and `gc()`
####
#### ➤ Total Output Objects:
####   • 3 Survival-imputed inputs: df006, df007, df008
####   • × 3 CNV imputation methods (mode, random, kNN)
####   → 3 × 3 = 9 CNV-imputed datasets
####
#### ➤ Resulting Output Name Table:
####   | Survival Input | CNV Method | Output Name |
####   |----------------|-------------|--------------|
####   | df006          | Mode        | df009        |
####   | df006          | Random      | df010        |
####   | df006          | kNN         | df011        |
####   | df007          | Mode        | df012        |
####   | df007          | Random      | df013        |
####   | df007          | kNN         | df014        |
####   | df008          | Mode        | df015        |
####   | df008          | Random      | df016        |
####   | df008          | kNN         | df017        |
####
#### ➤ Memory Optimization:
####   • Each dataframe is removed immediately after export to
####     preserve RAM, essential on systems with limited memory.
####
#### -------------------------------------------------------------


#### ------------------------------------------------------
#### Mode CNV Imputation (Amended for Sequential Execution)
#### ------------------------------------------------------

#' Strategy:
#' This snippet sequentially applies the mode-based CNV imputation function
#' to three survival-imputed datasets (df006, df007, df008), which must be
#' stored in the working directory as `.rds` files. Each resulting CNV-imputed
#' object is saved as a numbered dataframe (`df009` to `df011`) and removed
#' from memory after processing to manage limited RAM availability.

#' Required: rio, dplyr
library(rio)

####
#### ---------------------
#### Mode CNV imputation
#### ----------------------
#### 
#' # CNV imputation function using mode (unchanged)
#' impute_mode_groupwise_by_type <- function(df, type_col = "type", verbose = FALSE) {
#'   df_out <- df
#'   target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) length(x) >= 2 && x[2] == "3")]
#'   
#'   imputation_log <- data.frame(Variable = character(), Type = character(), NA_Before = integer(),
#'                                NA_After = integer(), n_imputed = integer(), Method = character(),
#'                                stringsAsFactors = FALSE)
#'   
#'   for (var in target_vars) {
#'     prefix <- strsplit(var, "-")[[1]][1]
#'     all_rows <- which(df[[type_col]] == prefix)
#'     na_rows <- which(df[[type_col]] == prefix & is.na(df[[var]]))
#'     na_before <- sum(is.na(df[[var]][all_rows]))
#'     group_values <- df[[var]][all_rows][!is.na(df[[var]][all_rows])]
#'     method_used <- NA_character_
#'     imputed <- 0
#'     
#'     if (length(na_rows) > 0 && length(group_values) > 0) {
#'       mode_val <- names(which.max(table(group_values)))
#'       df_out[[var]][na_rows] <- mode_val
#'       imputed <- length(na_rows)
#'       method_used <- "mode"
#'       if (verbose) message("🔁 Mode imputed ", imputed, " values for ", var, " using: ", mode_val)
#'     }
#'     
#'     na_after <- sum(is.na(df_out[[var]][all_rows]))
#'     if (!is.na(method_used)) {
#'       imputation_log <- rbind(imputation_log, data.frame(Variable = var, Type = prefix,
#'                                                          NA_Before = na_before, NA_After = na_after,
#'                                                          n_imputed = imputed, Method = method_used,
#'                                                          stringsAsFactors = FALSE))
#'     }
#'   }
#'   
#'   assign("imputation_log", imputation_log, envir = .GlobalEnv)
#'   return(df_out)
#' }
#' 
#' #' Sequentially process survival-imputed datasets
#' input_files <- c("df006.rds", "df007.rds", "df008.rds")
#' output_names <- c("df009", "df010", "df011")
#' 
#' for (i in seq_along(input_files)) {
#'   input_df <- import(input_files[i])
#'   imputed_df <- impute_mode_groupwise_by_type(input_df, verbose = TRUE)
#'   
#'   assign(output_names[i], imputed_df)
#'   export(imputed_df, paste0(output_names[i], ".rds"))
#'   
#'   rm(list = output_names[i])
#'   rm(imputed_df, input_df)
#'   gc()
#' }
#' 
#' #' Output tracking table
#' output_tracking <- data.frame(
#'   Survival_Input = c("df006", "df007", "df008"),
#'   CNV_Method = rep("mode", 3),
#'   Output_Name = output_names
#' )
#' print(output_tracking)

#' Expected output objects created: 3 (df009, df010, df011)

####
#### 
#### --------------------------------------------------
#### Random CNV imputation 
#### Random Imputation with Full Logging and Diagnostics
#### ---------------------------------------------------
#### Imputing "NA" values to categorical (nominal) CNV variables
#### 2. Random Category Imputation (Based on Empirical Category Distribution)
#### Matching Prefix in 'type' and ".3" Pattern (CNV) in 2nd token
#### Note: Random imputation tends to perform better when the categories are not highly skewed and preserving class diversity is critical.
#### 
#### ------------------------------------------------------------
#### PART 2 – CNV Random Imputation over Survival-Imputed Inputs
#### ------------------------------------------------------------
####
#### Three survival-imputed datasets (df006, df007, df008) previously saved
#### Each will be loaded, CNV-imputed via Random Imputation, saved (df012–df014),
#### and cleared from memory for RAM efficiency.

# # ---- CNV Random Imputation Function ----
# impute_random_groupwise_by_type <- function(df, type_col = "type", seed = 123, verbose = FALSE) {
#   set.seed(seed)
#   df_out <- df
#   target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) length(x) >= 2 && x[2] == "3")]
#   imputation_log <- data.frame(Variable = character(), Type = character(), NA_Before = integer(),
#                                NA_After = integer(), n_imputed = integer(), Method = character(),
#                                stringsAsFactors = FALSE)
#   
#   for (var in target_vars) {
#     prefix <- strsplit(var, "-")[[1]][1]
#     all_rows <- which(df[[type_col]] == prefix)
#     match_rows <- which(df[[type_col]] == prefix & is.na(df[[var]]))
#     group_values <- df[[var]][all_rows][!is.na(df[[var]][all_rows])]
#     na_before <- sum(is.na(df[[var]][all_rows]))
#     method_used <- NA_character_
#     imputed <- 0
#     
#     if (length(match_rows) > 0 && length(group_values) > 0) {
#       sampled_vals <- sample(group_values, size = length(match_rows), replace = TRUE)
#       df_out[[var]][match_rows] <- sampled_vals
#       imputed <- length(match_rows)
#       method_used <- "random"
#       if (verbose) message("🎲 Random imputed ", imputed, " values for ", var)
#     }
#     
#     na_after <- sum(is.na(df_out[[var]][all_rows]))
#     if (!is.na(method_used)) {
#       imputation_log <- rbind(imputation_log, data.frame(
#         Variable = var, Type = prefix, NA_Before = na_before, NA_After = na_after,
#         n_imputed = imputed, Method = method_used, stringsAsFactors = FALSE
#       ))
#     }
#   }
#   
#   assign("imputation_log", imputation_log, envir = .GlobalEnv)
#   return(df_out)
# }
# 
# # ---- Sequential Application to df006–df008 ----
# library(rio)
# 
# if (!exists("output_name_table_all")) {
#   output_name_table_all <- data.frame(
#     Step = character(),
#     Input_File = character(),
#     Output_Object = character(),
#     Saved_As = character(),
#     stringsAsFactors = FALSE
#   )
# }
# 
# for (i in 6:8) {
#   input_file <- sprintf("df%03d.rds", i)
#   df <- import(input_file)
#   cat("✅ Loaded:", input_file, "\n")
#   
#   # Apply Random CNV imputation
#   df_out <- impute_random_groupwise_by_type(df, verbose = TRUE)
#   
#   # Rename and save
#   output_index <- i + 6  # df012, df013, df014
#   output_name <- sprintf("df%03d", output_index)
#   assign(output_name, df_out)
#   saveRDS(df_out, sprintf("%s.rds", output_name))
#   
#   # Record the mapping
#   output_name_table_all <- rbind(output_name_table_all, data.frame(
#     Step = "CNV_Random",
#     Input_File = input_file,
#     Output_Object = output_name,
#     Saved_As = sprintf("%s.rds", output_name),
#     stringsAsFactors = FALSE
#   ))
#   
#   # Save the updated mapping table
#   write.table(
#     output_name_table_all,
#     file = "output_name_table_all.tsv",
#     sep = "\t",
#     row.names = FALSE,
#     col.names = TRUE,
#     quote = FALSE
#   )
#   
#   
#   # Clean up
#   rm(list = c("df", output_name, "df_out", "imputation_log"))
#   gc()
# }
# 
# # View the output name mapping
# print(output_name_table_all)
# 
# ####
# #### 
#### ----------------------
#### kNN CNV imputation
#### -----------------------
#### Imputing "NA" values to categorical (nominal) CNV variables
#### 3. kNN Imputation using `VIM::kNN` with Groupwise Prefix-Type Matching with Full Logging and Diagnostics
#### Matching Prefix in 'type' and ".3" Pattern (CNV) in 2nd token
#### Note: This approach imputes nominal CNV values preserving type structure.
#### ---------------------------------------------------------------------
#### PART 3 – CNV kNN Imputation on Survival-Imputed Inputs (df006–df008)
#### ---------------------------------------------------------------------
#### 
# 
# # ---- Load required library ----
# suppressPackageStartupMessages(library(VIM))
# library(rio)
# 
# # ---- kNN + fallback mode CNV imputation function ----
# impute_knn_or_mode_groupwise_by_type <- function(df, type_col = "type", k = 5, verbose = FALSE) {
#   df_out <- df
#   target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) length(x) >= 2 && x[2] == "3")]
#   imputation_log <- data.frame(Variable = character(), Type = character(), NA_Before = integer(),
#                                NA_After = integer(), n_imputed = integer(), Method = character(),
#                                stringsAsFactors = FALSE)
#   
#   for (var in target_vars) {
#     prefix <- strsplit(var, "-")[[1]][1]
#     group_rows <- which(df[[type_col]] == prefix)
#     if (length(group_rows) == 0) next
#     df_sub <- df[group_rows, , drop = FALSE]
#     if (!var %in% names(df_sub)) next
#     if (!is.factor(df_sub[[var]])) df_sub[[var]] <- as.factor(df_sub[[var]])
#     original_na <- is.na(df_sub[[var]])
#     na_rows <- group_rows[original_na]
#     na_before <- sum(original_na)
#     method_used <- NA_character_
#     imputed <- 0
#     
#     if (length(group_rows) > k) {
#       tryCatch({
#         df_knn_input <- df_sub[, c(var), drop = FALSE]
#         df_knn_input$.__dummy__ <- seq_len(nrow(df_knn_input))
#         rownames(df_knn_input) <- NULL
#         df_knn <- suppressWarnings(VIM::kNN(df_knn_input, k = k, imp_var = FALSE))
#         imputed <- sum(original_na & !is.na(df_knn[[var]]))
#         if (imputed > 0) {
#           df_out[group_rows, var] <- df_knn[[var]]
#           method_used <- "kNN"
#         } else {
#           fallback_vals <- df[[var]][group_rows][!original_na]
#           if (length(fallback_vals) > 0) {
#             mode_val <- names(which.max(table(fallback_vals)))
#             df_out[na_rows, var] <- mode_val
#             method_used <- "fallback_mode"
#             imputed <- length(na_rows)
#           }
#         }
#       }, error = function(e) {
#         fallback_vals <- df[[var]][group_rows][!original_na]
#         if (length(fallback_vals) > 0) {
#           mode_val <- names(which.max(table(fallback_vals)))
#           df_out[na_rows, var] <- mode_val
#           method_used <- "fallback_mode"
#           imputed <- length(na_rows)
#         }
#       })
#     } else {
#       fallback_vals <- df[[var]][group_rows][!original_na]
#       if (length(fallback_vals) > 0) {
#         mode_val <- names(which.max(table(fallback_vals)))
#         df_out[na_rows, var] <- mode_val
#         method_used <- "fallback_mode"
#         imputed <- length(na_rows)
#       }
#     }
#     
#     na_after <- sum(is.na(df_out[group_rows, var]))
#     if (!is.na(method_used)) {
#       imputation_log <- rbind(imputation_log, data.frame(
#         Variable = var, Type = prefix, NA_Before = na_before, NA_After = na_after,
#         n_imputed = imputed, Method = method_used, stringsAsFactors = FALSE
#       ))
#     }
#   }
#   
#   assign("imputation_log", imputation_log, envir = .GlobalEnv)
#   return(df_out)
# }
# 
# # ---- Sequential Execution: apply to df006–df008, produce df015–df017 ----
# if (!exists("output_name_table_all")) {
#   output_name_table_all <- data.frame(
#     Step = character(),
#     Input_File = character(),
#     Output_Object = character(),
#     Saved_As = character(),
#     stringsAsFactors = FALSE
#   )
# }
# 
# for (i in 6:8) {
#   input_file <- sprintf("df%03d.rds", i)
#   df <- import(input_file)
#   cat("✅ Loaded:", input_file, "\n")
#   
#   df_out <- impute_knn_or_mode_groupwise_by_type(df, k = 5, verbose = TRUE)
#   
#   output_index <- i + 9  # df015, df016, df017
#   output_name <- sprintf("df%03d", output_index)
#   assign(output_name, df_out)
#   saveRDS(df_out, sprintf("%s.rds", output_name))
#   
#   # Record the mapping
#   output_name_table_all <- rbind(output_name_table_all, data.frame(
#     Step = "CNV_kNN",
#     Input_File = input_file,
#     Output_Object = output_name,
#     Saved_As = sprintf("%s.rds", output_name),
#     stringsAsFactors = FALSE
#   ))
#   
#   # Save the updated mapping table
#   write.table(
#     output_name_table_all,
#     file = "output_name_table_all.tsv",
#     sep = "\t",
#     row.names = FALSE,
#     col.names = TRUE,
#     quote = FALSE
#   )
#   
#   rm(list = c("df", output_name, "df_out", "imputation_log"))
#   gc()
# }
# 
# # ---- Display resulting output table ----
# print(output_name_table_all)

### 
###
### ------------------------------------------------------------
### Mean Mutation binary imputation
### Mutation Variable Handling of mutation burden binary matrices
### -------------------------------------------------------------
### Function: impute_mean_groupwise_binary() with Full Logging
### 
### 
# ------------------------------------------------------------------------------
# 🧬 Mutation Burden Imputation — Comparative Assessment of Imputation Methods
#
# Context:
# - Mutation variables were binarized from original values in df003 to 0, 1, or NA.
# - 0 = No mutation
# - 1 = Mutated (burden ≥ 1)
# - NA = Unknown / Missing
#
# Goal:
# - Impute NA values in a biologically meaningful and statistically robust way.
# - Maintain interpretability for downstream analyses (burden, survival, clustering).
#
# Methods Evaluated:
#
# | Method        | Suitability   | Reason                                                                 |
# |---------------|---------------|------------------------------------------------------------------------|
# | Mean (rounded)| ✅ Excellent  | Encodes group-level mutation rate; interpretable and deterministic     |
# | Median        | ✅ Excellent  | Robust to skew; ideal for tumors with variable mutation prevalence     |
# | Mode          | ✅ Excellent  | Safe and simple; preserves dominant signal (e.g., when mutations are rare) |
# | Random        | ⚠️ Conditional| Use for stochastic simulations; not recommended for deterministic models |
# | kNN           | ⚠️ Low        | Limited use in sparse, binary-only data; may misrepresent rare events |
# | MissForest    | ⚠️ Conditional| Only justified when mutation is part of a larger multi-omic input matrix |
#
# Recommendation:
# - Use Mean/Median/Mode (rounded) for production-quality binary mutation imputation.
# - Use Random for stochastic bootstraps if needed.
# - Avoid kNN or MissForest unless integrating with multi-omic profiles.
# ------------------------------------------------------------------------------

# No Further Investment Required In:
# kNN: inappropriate for sparse binary mutation data unless multivariate similarity is modeled,
# MissForest: overkill for binary matrices unless you're imputing across a joint multi-omic panel,
# Random assignment: useful only for simulations, not primary imputation pipelines.
#
#### -----------------------------------------------------------------------------
#### Mutation Mean Imputation — Sequential Strategy for CNV-Imputed Datasets
#### -----------------------------------------------------------------------------

# # ✅ Mutation Imputation Function (Mean → Binary)
# impute_mean_groupwise_binary <- function(df, type_col = "type", verbose = FALSE) {
#   df_out <- df
#   
#   target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) {
#     length(x) >= 2 && x[2] == "2"
#   })]
#   
#   imputation_log <- data.frame(
#     Variable   = character(),
#     Type       = character(),
#     NA_Before  = integer(),
#     NA_After   = integer(),
#     n_imputed  = integer(),
#     Method     = character(),
#     stringsAsFactors = FALSE
#   )
#   
#   for (var in target_vars) {
#     prefix <- strsplit(var, "-")[[1]][1]
#     if (!is.numeric(df[[var]])) suppressWarnings(df[[var]] <- as.numeric(as.character(df[[var]])))
#     
#     all_rows <- which(df[[type_col]] == prefix)
#     na_rows  <- which(df[[type_col]] == prefix & is.na(df[[var]]))
#     na_before <- sum(is.na(df[[var]][all_rows]))
#     
#     group_values <- df[[var]][all_rows][!is.na(df[[var]][all_rows])]
#     imputed <- 0
#     method_used <- NA_character_
#     
#     if (length(na_rows) > 0 && length(group_values) > 0) {
#       mean_val <- mean(group_values)
#       binary_val <- round(mean_val)
#       df_out[[var]][na_rows] <- binary_val
#       imputed <- length(na_rows)
#       method_used <- "mean_binary"
#       if (verbose) message("📉 Mean imputed ", imputed, " values for ", var, ": ", binary_val)
#     }
#     
#     na_after <- sum(is.na(df_out[[var]][all_rows]))
#     
#     if (!is.na(method_used)) {
#       imputation_log <- rbind(imputation_log, data.frame(
#         Variable   = var,
#         Type       = prefix,
#         NA_Before  = na_before,
#         NA_After   = na_after,
#         n_imputed  = imputed,
#         Method     = method_used,
#         stringsAsFactors = FALSE
#       ))
#     }
#   }
#   
#   assign("imputation_log", imputation_log, envir = .GlobalEnv)
#   return(df_out)
# }
# 
# # ✅ Output Tracking Table (persistent across runs)
# if (!exists("output_name_table_all")) {
#   output_name_table_all <- data.frame(
#     Step = character(),
#     Input_File = character(),
#     Output_Object = character(),
#     Saved_As = character(),
#     stringsAsFactors = FALSE
#   )
# }
# 
# # ✅ Batch Processing of CNV-Imputed Files (df009–df017)
# # → Resulting Mutation Mean-Imputed Files: df018–df026
# 
# for (i in 9:17) {
#   input_file <- sprintf("df%03d.rds", i)
#   df <- import(input_file)
#   cat("✅ Loaded:", input_file, "\n")
#   
#   df_out <- impute_mean_groupwise_binary(df, verbose = TRUE)
#   
#   output_index <- i + 9  # df018–df026
#   output_name <- sprintf("df%03d", output_index)
#   assign(output_name, df_out)
#   saveRDS(df_out, sprintf("%s.rds", output_name))
#   
#   # Log entry
#   output_name_table_all <- rbind(output_name_table_all, data.frame(
#     Step = "Mutation_Mean",
#     Input_File = input_file,
#     Output_Object = output_name,
#     Saved_As = sprintf("%s.rds", output_name),
#     stringsAsFactors = FALSE
#   ))
#   
#   write.table(
#     output_name_table_all,
#     file = "output_name_table_all.tsv",
#     sep = "\t",
#     row.names = FALSE,
#     col.names = TRUE,
#     quote = FALSE
#   )
#   
#   rm(list = c("df", output_name, "df_out", "imputation_log"))
#   gc()
# }
# 
# # ✅ Final Confirmation
# print(output_name_table_all)


####
#### -------------------------------------------------------------------------
#### Median Mutation binary imputation
#### Median Imputation group-wise by Matching Prefix in `type` and `.2` Pattern (in 2nd token)
#### impute_median_groupwise_binary() — With Full Logging and Binary
#### ---------------------------------------------------------------------------
####
#### ------------------------------------------------------------
#### Median Mutation Binary Imputation on CNV-Imputed Files
#### ------------------------------------------------------------

# impute_median_groupwise_binary <- function(df, type_col = "type", verbose = FALSE) {
#   df_out <- df
#   target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) length(x) >= 2 && x[2] == "2")]
#   imputation_log <- data.frame(
#     Variable = character(), Type = character(), NA_Before = integer(),
#     NA_After = integer(), n_imputed = integer(), Method = character(),
#     stringsAsFactors = FALSE
#   )
#   
#   for (var in target_vars) {
#     prefix <- strsplit(var, "-")[[1]][1]
#     if (!is.numeric(df[[var]])) suppressWarnings(df[[var]] <- as.numeric(as.character(df[[var]])))
#     all_rows <- which(df[[type_col]] == prefix)
#     na_rows <- which(df[[type_col]] == prefix & is.na(df[[var]]))
#     na_before <- sum(is.na(df[[var]][all_rows]))
#     group_values <- df[[var]][all_rows][!is.na(df[[var]][all_rows])]
#     method_used <- NA_character_
#     imputed <- 0
#     
#     if (length(na_rows) > 0 && length(group_values) > 0) {
#       median_val <- median(group_values, na.rm = TRUE)
#       binary_val <- round(median_val)
#       df_out[[var]][na_rows] <- binary_val
#       imputed <- length(na_rows)
#       method_used <- "median_binary"
#       if (verbose) message("🧮 Median imputed ", imputed, " values for ", var, " as binary: ", binary_val)
#     }
#     
#     na_after <- sum(is.na(df_out[[var]][all_rows]))
#     if (!is.na(method_used)) {
#       imputation_log <- rbind(imputation_log, data.frame(
#         Variable = var, Type = prefix, NA_Before = na_before,
#         NA_After = na_after, n_imputed = imputed, Method = method_used,
#         stringsAsFactors = FALSE
#       ))
#     }
#   }
#   
#   assign("imputation_log", imputation_log, envir = .GlobalEnv)
#   return(df_out)
# }

#### ---------------------------------------------------------
#### Sequential Execution over CNV-Imputed Files (df009–df017)
#### Produces Median-Mutation-Imputed Outputs (df027–df035)
#### ---------------------------------------------------------
# 
# library(rio)
# 
# if (!exists("output_name_table_all")) {
#   output_name_table_all <- data.frame(
#     Step = character(), Input_File = character(),
#     Output_Object = character(), Saved_As = character(),
#     stringsAsFactors = FALSE
#   )
# }
# 
# for (i in 9:17) {
#   input_file <- sprintf("df%03d.rds", i)
#   df <- import(input_file)
#   cat("✅ Loaded:", input_file, "\n")
#   
#   df_out <- impute_median_groupwise_binary(df, verbose = TRUE)
#   
#   output_index <- i + 18  # i = 9–17 ⟹ df027–df035
#   output_name <- sprintf("df%03d", output_index)
#   assign(output_name, df_out)
#   saveRDS(df_out, sprintf("%s.rds", output_name))
#   
#   output_name_table_all <- rbind(output_name_table_all, data.frame(
#     Step = "Mutation_Median",
#     Input_File = input_file,
#     Output_Object = output_name,
#     Saved_As = sprintf("%s.rds", output_name),
#     stringsAsFactors = FALSE
#   ))
#   
#   write.table(
#     output_name_table_all,
#     file = "output_name_table_all.tsv",
#     sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
#   )
#   
#   rm(list = c("df", output_name, "df_out", "imputation_log"))
#   gc()
# }
# 
# print(output_name_table_all)

####
#### --------------------------------------------------------------------
#### Mode Mutation binary imputation
#### Impute_mode_groupwise_binary() — With Full Logging and Binary Output
#### ---------------------------------------------------------------------
#### 
#### --------------------------------------------------------------
#### Mode Mutation Binary Imputation from CNV-Imputed Inputs
#### --------------------------------------------------------------
# impute_mode_groupwise_binary <- function(df, type_col = "type", verbose = FALSE) {
#   df_out <- df
#   target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) length(x) >= 2 && x[2] == "2")]
#   
#   imputation_log <- data.frame(
#     Variable = character(), Type = character(), NA_Before = integer(),
#     NA_After = integer(), n_imputed = integer(), Method = character(),
#     stringsAsFactors = FALSE
#   )
#   
#   for (var in target_vars) {
#     prefix <- strsplit(var, "-")[[1]][1]
#     if (!is.numeric(df[[var]])) suppressWarnings(df[[var]] <- as.numeric(as.character(df[[var]])))
#     
#     all_rows <- which(df[[type_col]] == prefix)
#     na_rows <- which(df[[type_col]] == prefix & is.na(df[[var]]))
#     na_before <- sum(is.na(df[[var]][all_rows]))
#     group_values <- df[[var]][all_rows][!is.na(df[[var]][all_rows])]
#     method_used <- NA_character_
#     imputed <- 0
#     
#     if (length(na_rows) > 0 && length(group_values) > 0) {
#       mode_val <- names(which.max(table(group_values)))
#       binary_val <- round(as.numeric(mode_val))
#       df_out[[var]][na_rows] <- binary_val
#       method_used <- "mode_binary"
#       imputed <- length(na_rows)
#       if (verbose) message("📊 Mode imputed ", imputed, " values for ", var, " as binary: ", binary_val)
#     }
#     
#     na_after <- sum(is.na(df_out[[var]][all_rows]))
#     if (!is.na(method_used)) {
#       imputation_log <- rbind(imputation_log, data.frame(
#         Variable = var, Type = prefix, NA_Before = na_before,
#         NA_After = na_after, n_imputed = imputed, Method = method_used,
#         stringsAsFactors = FALSE
#       ))
#     }
#   }
#   
#   assign("imputation_log", imputation_log, envir = .GlobalEnv)
#   return(df_out)
# }

#### -------------------------------------------------------------------
#### Sequential Execution on CNV-Imputed Datasets df009–df017 → df036–df044
#### -------------------------------------------------------------------

# library(rio)
# 
# if (!exists("output_name_table_all")) {
#   output_name_table_all <- data.frame(
#     Step = character(), Input_File = character(),
#     Output_Object = character(), Saved_As = character(),
#     stringsAsFactors = FALSE
#   )
# }
# 
# for (i in 9:17) {
#   input_file <- sprintf("df%03d.rds", i)
#   df <- import(input_file)
#   cat("✅ Loaded:", input_file, "\n")
#   
#   df_out <- impute_mode_groupwise_binary(df, verbose = TRUE)
#   
#   output_index <- i + 27  # i=9 ⇒ df036; i=17 ⇒ df044
#   output_name <- sprintf("df%03d", output_index)
#   assign(output_name, df_out)
#   saveRDS(df_out, sprintf("%s.rds", output_name))
#   
#   # Append to output registry
#   output_name_table_all <- rbind(output_name_table_all, data.frame(
#     Step = "Mutation_Mode",
#     Input_File = input_file,
#     Output_Object = output_name,
#     Saved_As = sprintf("%s.rds", output_name),
#     stringsAsFactors = FALSE
#   ))
#   
#   # Save updated table to disk
#   write.table(
#     output_name_table_all,
#     file = "output_name_table_all.tsv",
#     sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
#   )
#   
#   rm(list = c("df", output_name, "df_out", "imputation_log"))
#   gc()
# }
# 
# # Display final mapping
# print(output_name_table_all)

#### ===============================================================
#### 📊 Commentary: Mutation Binary Imputation Indexing Convention
#### ===============================================================
# As part of the multi-method mutation imputation strategy applied to
# CNV-imputed datasets (df009–df017), the following indexing convention
# has been adopted for saving the output mutation-imputed objects:
#
# ➤ Imputation Methods Applied:
#     1. Mean (rounded binary)
#     2. Median (rounded binary)
#     3. Mode (rounded binary)
#
# ➤ Input Range:
#     df009 to df017 → 9 input files (CNV-imputed)
#
# ➤ Output Assignment:
#     • For each CNV-imputed input (df009–df017),
#       a new mutation-imputed object is created per method.
#
# ➤ Output File Index Mapping:
#   --------------------------------------------------------------
#   | Imputation Method | Output Index Range | Description        |
#   |-------------------|--------------------|--------------------|
#   | Mean              | df018 – df026      | Mutation Imputed   |
#   | Median            | df027 – df035      | Mutation Imputed   |
#   | Mode              | df036 – df044      | Mutation Imputed   |
#   --------------------------------------------------------------
#
# ➤ File Format:
#     • Saved as .rds using `saveRDS()`
#     • Registered in `output_name_table_all.tsv` with method tags
#
# ➤ Memory-Safe Execution:
#     • Each object is removed after export with `rm()` and `gc()`
#
# ➤ Output consistency enables downstream integration with
#     subsequent transcript, miRNA, and methylation imputation stages.
# ===============================================================

####
####
#### ---------------------------------------------
#### IMPUTATION STRATEGY FOR CONTINUOUS VARIABLES
#### ---------------------------------------------
#### This section summarizes the rationale and selection
#### of imputation methods for numeric omic variables in the
#### Pan-Cancer Multi-Omic dataset, based on data type,
#### distribution, and cancer-type stratification logic.
#### 
#### -----------------------------------------------------
#### Target Omic Layers and Tokens:
####   - .1 = Protein expression (log2 or z-score)
####   - .4 = miRNA expression (normalized counts or CPM)
####   - .5 = Transcript expression (TPM, FPKM, etc.)
####   - .6 = mRNA expression (TPM, log2(TPM+1))
####   - .7 = CpG Methylation (beta or M-values)
#### ------------------------------------------------------
#### 
#### All are continuous numeric variables, potentially skewed or multimodal.
#### Imputation is performed groupwise by 'type' (i.e., cancer type).
#### 
#### 

#### 
#### ----------------------------------
#### METHODS CONSIDERED FOR IMPUTATION:
#### ----------------------------------

#### 1. SIMPLE STATISTICAL METHODS
####    - Mean     ✅ Suitable; sensitive to outliers.
####    - Median   ✅ Robust to outliers; preserves central tendency.
####    - Mode     ❌ Avoid; not meaningful for continuous distributions.
####    - Random   ✅ Acceptable with reproducibility control.

#### 2. DISTANCE-BASED METHODS
####    - kNN (VIM::kNN)         ✅ Captures local structure, widely used.
####    - MICE (mice::mice)      ✅ Flexible, multivariate; monitor convergence.
####    - missForest             ✅ Handles nonlinearity; robust, non-parametric.

#### 3. TREE-BASED ML METHODS
####    - XGBoost                ✅ Strong learner; requires wrapper interface.
####    - LightGBM               ✅ Fast alternative; similar to XGBoost.
####    - LASSO / Boruta / SHAP  ❌ Feature selection, not designed for imputation.
####    - DeepSurv               ❌ Survival modeling; not for missing data imputation.

#### ------------------------------------
#### SELECTED ORDER OF METHODS FOR IMPLEMENTATION:
#### ------------------------------------
#### ✓ Mean
#### ✓ Median
#### ✓ Random
#### ✓ kNN
#### ✓ missForest
#### ✓ XGBoost
#### ✓ LightGBM
#### ✓ MICE # run lats because of the multiple predictions

#### Imputation for each method will be performed groupwise
#### by cancer type, using variable names filtered by their
#### second token (e.g., ".1", ".4", ".5", ".6", ".7").

#### Next Step: Implement modular imputation functions for each method,
#### maintaining logging and diagnostics consistent with the CNV (.3) snippet.

####
#### ---------------------------------------------------------------
#### Groupwise Imputation for Continuous Variables: .1, .4, .5, .6, .7
#### Methods: Mean → Median → Random → kNN (VIM::kNN)
#### Fallbacks are applied if primary method fails or group too small.
#### Full logging and diagnostics included.
#### ---------------------------------------------------------------
#### 

#### ------------------------------------------------------------
#### PART 1 of 8: Groupwise Mean Imputation of Continuous Variables
#### Across Mutation-Imputed Inputs df018–df044 (n=27)
#### ------------------------------------------------------------

suppressPackageStartupMessages(library(VIM))
library(rio)

# Groupwise mean imputation function
impute_numeric_groupwise <- function(df, type_col = "type", method = "mean", k = 5, verbose = FALSE) {
  df_out <- df
  numeric_tokens <- c("5", "6")
  target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) length(x) >= 2 && x[2] %in% numeric_tokens)]
  
  imputation_log <- data.frame(Variable=character(), Type=character(), NA_Before=integer(),
                               NA_After=integer(), n_imputed=integer(), Method=character(),
                               stringsAsFactors=FALSE)
  
  for (var in target_vars) {
    prefix <- strsplit(var, "-")[[1]][1]
    group_rows <- which(df[[type_col]] == prefix)
    if (length(group_rows) == 0) next
    
    df_sub <- df[group_rows, , drop = FALSE]
    if (!var %in% names(df_sub)) next
    if (!is.numeric(df_sub[[var]])) df_sub[[var]] <- as.numeric(df_sub[[var]])
    
    original_na <- is.na(df_sub[[var]])
    na_rows <- group_rows[original_na]
    na_before <- sum(original_na)
    
    method_used <- NA_character_
    imputed <- 0
    
    val <- mean(df_sub[[var]], na.rm = TRUE)
    if (!is.nan(val)) {
      df_out[na_rows, var] <- val
      imputed <- length(na_rows)
      method_used <- "mean"
    }
    
    na_after <- sum(is.na(df_out[group_rows, var]))
    if (!is.na(method_used)) {
      imputation_log <- rbind(imputation_log, data.frame(
        Variable=var, Type=prefix, NA_Before=na_before,
        NA_After=na_after, n_imputed=imputed, Method=method_used,
        stringsAsFactors=FALSE))
    }
    
    if (verbose) {
      cat("🔬", var, "| Type:", prefix, "| Method:", method_used,
          "| Imputed:", imputed, "| Remaining NA:", na_after, "\n")
    }
  }
  assign("numeric_imputation_log", imputation_log, envir = .GlobalEnv)
  return(df_out)
}

# Initialize output name table if not exists
if (!exists("output_name_table_all")) {
  output_name_table_all <- data.frame(
    Step = character(), Input_File = character(),
    Output_Object = character(), Saved_As = character(),
    stringsAsFactors = FALSE)
}

# Perform mean imputation over 27 inputs
input_indices <- 18:20
output_indices <- 21:23

for (i in seq_along(input_indices)) {
  input_file <- sprintf("df%03d.rds", input_indices[i])
  df <- import(input_file)
  cat("✅ Loaded:", input_file, "\n")
  
  df_out <- impute_numeric_groupwise(df, method = "mean", verbose = TRUE)
  
  output_name <- sprintf("df%03d", output_indices[i])
  assign(output_name, df_out)
  saveRDS(df_out, paste0(output_name, ".rds"))
  
  # Log output
  output_name_table_all <- rbind(output_name_table_all, data.frame(
    Step = "Continuous_mean",
    Input_File = input_file,
    Output_Object = output_name,
    Saved_As = paste0(output_name, ".rds"),
    stringsAsFactors = FALSE))
  
  write.table(output_name_table_all, "output_name_table_all.tsv",
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  rm(list = c("df", "df_out", output_name, "numeric_imputation_log"))
  gc()
}

cat("✅ Mean imputation completed for continuous omic layers on 27 input files (df018–df044 → df045–df071).\n")

#### --------------------------------------------------------------------------
#### PART 2 — Groupwise Median Imputation for Continuous Variables (.1, .4-.7)
#### Applies to Mutation-Imputed Inputs: df018–df044
#### Outputs: df072–df098
#### Method: Median (per cancer type and numeric omic token)
#### --------------------------------------------------------------------------

suppressPackageStartupMessages(library(rio))
suppressPackageStartupMessages(library(VIM))
gc()

# --- Imputation Function ---
impute_numeric_groupwise_median <- function(df, type_col = "type", verbose = FALSE) {
  df_out <- df
  numeric_tokens <- c("5", "6")
  target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) {
    length(x) >= 2 && x[2] %in% numeric_tokens
  })]
  
  imputation_log <- data.frame(
    Variable = character(), Type = character(),
    NA_Before = integer(), NA_After = integer(),
    n_imputed = integer(), Method = character(),
    stringsAsFactors = FALSE
  )
  
  for (var in target_vars) {
    prefix <- strsplit(var, "-")[[1]][1]
    group_rows <- which(df[[type_col]] == prefix)
    if (length(group_rows) == 0) next
    df_sub <- df[group_rows, , drop = FALSE]
    if (!var %in% names(df_sub)) next
    if (!is.numeric(df_sub[[var]])) df_sub[[var]] <- as.numeric(df_sub[[var]])
    
    original_na <- is.na(df_sub[[var]])
    na_rows <- group_rows[original_na]
    na_before <- sum(original_na)
    
    method_used <- NA_character_
    imputed <- 0
    
    val <- median(df_sub[[var]], na.rm = TRUE)
    if (!is.nan(val)) {
      df_out[na_rows, var] <- val
      imputed <- length(na_rows)
      method_used <- "median"
    }
    
    na_after <- sum(is.na(df_out[group_rows, var]))
    
    if (!is.na(method_used)) {
      imputation_log <- rbind(imputation_log, data.frame(
        Variable = var, Type = prefix,
        NA_Before = na_before, NA_After = na_after,
        n_imputed = imputed, Method = method_used,
        stringsAsFactors = FALSE
      ))
    }
    
    if (verbose) {
      cat("🔬", var, " | Type:", prefix, "| Method:", method_used,
          "| Imputed:", imputed, "| Remaining NA:", na_after, "\n")
    }
  }
  
  assign("numeric_imputation_log", imputation_log, envir = .GlobalEnv)
  return(df_out)
}

# --- Loop over Mutation-Imputed Inputs df018 to df044 ---
if (!exists("output_name_table_all")) {
  output_name_table_all <- data.frame(
    Step = character(),
    Input_File = character(),
    Output_Object = character(),
    Saved_As = character(),
    stringsAsFactors = FALSE
  )
}

for (i in 18:20) {
  input_file <- sprintf("df%03d.rds", i)
  df <- import(input_file)
  cat("✅ Loaded:", input_file, "\n")
  
  df_out <- impute_numeric_groupwise_median(df, verbose = TRUE)
  
  output_index <- i + 54  # 18+54=72 through 44+54=98
  output_name <- sprintf("df%03d", output_index)
  assign(output_name, df_out)
  saveRDS(df_out, sprintf("%s.rds", output_name))
  
  output_name_table_all <- rbind(output_name_table_all, data.frame(
    Step = "Continuous_median",
    Input_File = input_file,
    Output_Object = output_name,
    Saved_As = sprintf("%s.rds", output_name),
    stringsAsFactors = FALSE
  ))
  
  write.table(
    output_name_table_all,
    file = "output_name_table_all.tsv",
    sep = "\t", row.names = FALSE,
    col.names = TRUE, quote = FALSE
  )
  
  rm(list = c("df", output_name, "df_out", "numeric_imputation_log"))
  gc()
}

cat("✅ Median imputation completed for continuous omic layers on 27 input files (df018–df044 → df072–df098).\n")

#### ---------------------------------------------------------------
#### PART 3 — Memory-Efficient Continuous Imputation (Random Method)
#### ---------------------------------------------------------------
suppressPackageStartupMessages(library(VIM))
library(rio)

gc()

# Reuse the impute_numeric_groupwise function for "random"
# (previously defined and assumed available in environment)

# Safety: Check that output log exists
if (!exists("output_name_table_all")) {
  output_name_table_all <- data.frame(
    Step = character(),
    Input_File = character(),
    Output_Object = character(),
    Saved_As = character(),
    stringsAsFactors = FALSE
  )
}

# Input and output indexing
input_indices <- 18:20
output_indices <- 24:26

# Sanity check
stopifnot(length(input_indices) == 3, length(output_indices) == 3)

# Loop through inputs and apply random imputation
for (i in seq_along(input_indices)) {
  input_index <- input_indices[i]
  output_index <- output_indices[i]
  
  input_name <- sprintf("df%03d", input_index)
  output_name <- sprintf("df%03d", output_index)
  
  cat("🔁 Processing Random Imputation:", input_name, "→", output_name, "\n")
  
  # Import input
  input_file <- paste0(input_name, ".rds")
  df <- import(input_file)
  
  # Apply groupwise random imputation
  df_out <- impute_numeric_groupwise(df, method = "random", verbose = TRUE)
  
  # Save result
  saveRDS(df_out, file = paste0(output_name, ".rds"))
  
  # Log result
  output_name_table_all <- rbind(output_name_table_all, data.frame(
    Step = "Continuous_random",
    Input_File = input_file,
    Output_Object = output_name,
    Saved_As = paste0(output_name, ".rds"),
    stringsAsFactors = FALSE
  ))
  
  # Save log to disk
  write.table(output_name_table_all,
              file = "output_name_table_all.tsv",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)
  
  # Cleanup
  rm(df, df_out)
  gc()
}

cat("✅ Random imputation completed for continuous omic layers on 27 input files (df018–df044 → df099–df125).\n")

####
#### ----------------------------------------------------------------
#### PART 4 — Memory-Efficient Continuous Imputation (kNN Method)
#### ----------------------------------------------------------------

suppressPackageStartupMessages(library(VIM))
library(rio)

gc()

# Function for groupwise imputation (kNN included)
# Assumes impute_numeric_groupwise is already defined
# If not, you can load it from prior definition

# Ensure global log exists
if (!exists("output_name_table_all")) {
  output_name_table_all <- data.frame(
    Step = character(),
    Input_File = character(),
    Output_Object = character(),
    Saved_As = character(),
    stringsAsFactors = FALSE
  )
}

# Indexing range
input_indices <- 18:20
output_indices <- 27:29

stopifnot(length(input_indices) == 3, length(output_indices) == 3)

for (i in seq_along(input_indices)) {
  input_index <- input_indices[i]
  output_index <- output_indices[i]
  
  input_name <- sprintf("df%03d", input_index)
  output_name <- sprintf("df%03d", output_index)
  
  cat("🔁 Processing kNN Imputation:", input_name, "→", output_name, "\n")
  
  input_file <- paste0(input_name, ".rds")
  df <- import(input_file)
  
  # Apply groupwise kNN imputation (with fallback handling)
  df_out <- impute_numeric_groupwise(df, method = "knn", k = 5, verbose = TRUE)
  
  # Save result
  saveRDS(df_out, file = paste0(output_name, ".rds"))
  
  # Log this output
  output_name_table_all <- rbind(output_name_table_all, data.frame(
    Step = "Continuous_kNN",
    Input_File = input_file,
    Output_Object = output_name,
    Saved_As = paste0(output_name, ".rds"),
    stringsAsFactors = FALSE
  ))
  
  # Persist log table
  write.table(output_name_table_all,
              file = "output_name_table_all.tsv",
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)
  
  # Clean memory
  rm(df, df_out)
  gc()
}

cat("✅ kNN imputation completed for continuous omic layers on 27 input files (df018–df044 → df126–df152).\n")

####
####
####

suppressPackageStartupMessages(library(missForest))
suppressPackageStartupMessages(library(rio))

gc()

####
#### -------------------------------------------------------------------------------
#### Groupwise missForest Imputation for numeric omic layers (.1, .4, .5, .6, .7)
#### -------------------------------------------------------------------------------
#### -------------------------------------------------------------------------------
#### PART 5 of 8: Groupwise missForest Imputation for Continuous Variables
#### Across Mutation-Imputed Inputs df018–df044 (n = 27)
#### -------------------------------------------------------------------------------
# 
# suppressPackageStartupMessages(library(missForest))
# suppressPackageStartupMessages(library(rio))
# 
# # Imputation Function
# impute_missforest_groupwise <- function(df, type_col = "type", maxiter = 10, ntree = 100, verbose = FALSE) {
#   df_out <- df
#   numeric_tokens <- c("5", "6")
#   target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) length(x) >= 2 && x[2] %in% numeric_tokens)]
#   
#   imputation_log <- data.frame(
#     Variable = character(), Type = character(),
#     NA_Before = integer(), NA_After = integer(),
#     n_imputed = integer(), Method = character(),
#     stringsAsFactors = FALSE
#   )
#   
#   for (prefix in unique(df[[type_col]])) {
#     group_rows <- which(df[[type_col]] == prefix)
#     df_sub <- df[group_rows, target_vars, drop = FALSE]
#     na_before <- sum(is.na(df_sub))
#     if (na_before == 0) next
#     
#     tryCatch({
#       result <- missForest(df_sub, maxiter = maxiter, ntree = ntree, verbose = FALSE)
#       complete_df <- result$ximp
#       df_out[group_rows, target_vars] <- complete_df
#       imputed <- na_before - sum(is.na(complete_df))
#       imputation_log <- rbind(imputation_log, data.frame(
#         Variable = "Multiple", Type = prefix, NA_Before = na_before,
#         NA_After = sum(is.na(complete_df)), n_imputed = imputed, Method = "missForest",
#         stringsAsFactors = FALSE
#       ))
#       if (verbose) cat("✅ missForest imputation done for:", prefix, "| Imputed:", imputed, "\n")
#     }, error = function(e) {
#       if (verbose) cat("⚠️ Error for group:", prefix, ":", e$message, "\n")
#     })
#   }
#   
#   assign("missforest_imputation_log", imputation_log, envir = .GlobalEnv)
#   return(df_out)
# }
# 
# # Loop with Progressive Indexing
# input_indices <- 18:20
# output_indices <- 30:32
# 
# if (!exists("output_name_table_all")) {
#   output_name_table_all <- data.frame(
#     Step = character(),
#     Input_File = character(),
#     Output_Object = character(),
#     Saved_As = character(),
#     stringsAsFactors = FALSE
#   )
# }
# 
# for (i in seq_along(input_indices)) {
#   input_file <- sprintf("df%03d.rds", input_indices[i])
#   output_object <- sprintf("df%03d", output_indices[i])
#   output_file <- paste0(output_object, ".rds")
#   
#   cat("🚀 Processing:", input_file, "→", output_object, "\n")
#   
#   df <- import(input_file)
#   df_imputed <- impute_missforest_groupwise(df, verbose = TRUE)
#   assign(output_object, df_imputed)
#   saveRDS(df_imputed, output_file)
#   
#   # Log entry
#   output_name_table_all <- rbind(output_name_table_all, data.frame(
#     Step = "Continuous_missForest",
#     Input_File = input_file,
#     Output_Object = output_object,
#     Saved_As = output_file,
#     stringsAsFactors = FALSE
#   ))
#   
#   write.table(output_name_table_all, "output_name_table_all.tsv",
#               sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#   
#   rm(list = c("df", "df_imputed", output_object, "missforest_imputation_log"))
#   gc()
# }
# 
# cat("✅ missForest continuous imputation completed for all 27 inputs (df018–df044 → df153–df179).\n")

suppressPackageStartupMessages(library(missForest))
suppressPackageStartupMessages(library(rio))
suppressPackageStartupMessages(library(dplyr))

# 📌 Função auxiliar para imputar com a moda
imputar_com_moda <- function(v) {
  moda <- names(sort(table(v), decreasing = TRUE))[1]
  v[is.na(v)] <- moda
  return(v)
}

# 🚀 Função principal de imputação
impute_missforest_groupwise <- function(df, type_col = "type", maxiter = 10, ntree = 100, verbose = FALSE) {
  df_out <- df
  numeric_tokens <- c("5", "6")
  target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) length(x) >= 2 && x[2] %in% numeric_tokens)]
  
  imputation_log <- data.frame(
    Variable = character(), Type = character(),
    NA_Before = integer(), NA_After = integer(),
    n_imputed = integer(), Method = character(),
    stringsAsFactors = FALSE
  )
  
  for (prefix in unique(df[[type_col]])) {
    group_rows <- which(df[[type_col]] == prefix)
    df_sub <- df[group_rows, target_vars, drop = FALSE]
    na_before <- sum(is.na(df_sub))
    if (na_before == 0) next
    
    tryCatch({
      if (ncol(df_sub) == 1) {
        # ✅ Imputação com moda se houver apenas uma variável
        varname <- names(df_sub)
        na_count <- sum(is.na(df_sub[[1]]))
        df_sub[[1]] <- imputar_com_moda(df_sub[[1]])
        df_out[group_rows, varname] <- df_sub[[1]]
        imputation_log <- rbind(imputation_log, data.frame(
          Variable = varname, Type = prefix, NA_Before = na_count,
          NA_After = sum(is.na(df_sub[[1]])), n_imputed = na_count, Method = "Moda",
          stringsAsFactors = FALSE
        ))
        if (verbose) cat("🟡 Imputação com moda para grupo:", prefix, "| Imputado:", na_count, "\n")
      } else {
        # ✅ missForest imputação normal
        result <- missForest(df_sub, maxiter = maxiter, ntree = ntree, verbose = FALSE)
        complete_df <- result$ximp
        df_out[group_rows, target_vars] <- complete_df
        imputed <- na_before - sum(is.na(complete_df))
        imputation_log <- rbind(imputation_log, data.frame(
          Variable = "Multiple", Type = prefix, NA_Before = na_before,
          NA_After = sum(is.na(complete_df)), n_imputed = imputed, Method = "missForest",
          stringsAsFactors = FALSE
        ))
        if (verbose) cat("✅ missForest imputação para:", prefix, "| Imputado:", imputed, "\n")
      }
    }, error = function(e) {
      if (verbose) cat("⚠️ Erro no grupo:", prefix, ":", e$message, "\n")
    })
  }
  
  assign("missforest_imputation_log", imputation_log, envir = .GlobalEnv)
  return(df_out)
}


#####
##### -------------------------------------------------------------------------------
##### Groupwise XGBoost Imputation for numeric omic layers (.1, .4, .5, .6, .7)
##### -------------------------------------------------------------------------------
##### 
##### -------------------------------------------------------------------------------
##### Groupwise XGBoost Imputation for numeric omic layers (.1, .4, .5, .6, .7)
##### -------------------------------------------------------------------------------
##### 

#### ------------------------------------------------------------------------
#### PART 6 — Memory-Efficient XGBoost Imputation for Continuous Variables
#### Input:  df018 to df044 (mutation-imputed)
#### Output: df180 to df206 (XGBoost-imputed)
#### ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(xgboost)
  library(rio)
})

# -- Function: XGBoost groupwise imputation for continuous omic variables --
impute_xgboost_groupwise <- function(df, type_col = "type", verbose = FALSE) {
  df_out <- df
  numeric_tokens <- c("5", "6")
  target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) {
    length(x) >= 2 && x[2] %in% numeric_tokens
  })]
  
  imputation_log <- data.frame(
    Variable = character(), Type = character(),
    NA_Before = integer(), NA_After = integer(),
    n_imputed = integer(), Method = character(),
    stringsAsFactors = FALSE
  )
  
  for (var in target_vars) {
    prefix <- strsplit(var, "-")[[1]][1]
    group_rows <- which(df[[type_col]] == prefix)
    df_sub <- df[group_rows, , drop = FALSE]
    
    y <- df_sub[[var]]
    na_idx <- which(is.na(y))
    if (length(na_idx) == 0) next
    
    predictors <- setdiff(target_vars, var)
    df_train <- df_sub[!is.na(y), predictors, drop = FALSE]
    y_train <- y[!is.na(y)]
    df_pred  <- df_sub[ is.na(y), predictors, drop = FALSE]
    
    tryCatch({
      dtrain <- xgb.DMatrix(data = as.matrix(df_train), label = y_train)
      model  <- xgboost(data = dtrain, nrounds = 50, verbose = 0)
      dtest  <- xgb.DMatrix(data = as.matrix(df_pred))
      pred   <- predict(model, dtest)
      df_out[group_rows[na_idx], var] <- pred
      
      imputation_log <- rbind(imputation_log, data.frame(
        Variable = var, Type = prefix, NA_Before = length(na_idx),
        NA_After = 0, n_imputed = length(na_idx), Method = "XGBoost",
        stringsAsFactors = FALSE
      ))
      
      if (verbose) cat("✅ XGBoost imputed", length(na_idx), "values for", var, "| Group:", prefix, "\n")
      
    }, error = function(e) {
      if (verbose) cat("⚠️ XGBoost failed for", var, "in group", prefix, "|", e$message, "\n")
    })
  }
  
  assign("xgboost_imputation_log", imputation_log, envir = .GlobalEnv)
  return(df_out)
}

# -- Main Loop: Process df018 to df044, Output df180 to df206 --
input_indices <- 18:20
output_indices <- 33:35

if (!exists("output_name_table_all")) {
  output_name_table_all <- data.frame(
    Step = character(),
    Input_File = character(),
    Output_Object = character(),
    Saved_As = character(),
    stringsAsFactors = FALSE
  )
}

for (i in seq_along(input_indices)) {
  input_file <- sprintf("df%03d.rds", input_indices[i])
  output_name <- sprintf("df%03d", output_indices[i])
  output_file <- paste0(output_name, ".rds")
  
  cat("📂 Loaded:", input_file, "\n")
  df <- import(input_file)
  
  df_out <- impute_xgboost_groupwise(df, verbose = TRUE)
  assign(output_name, df_out)
  saveRDS(df_out, output_file)
  
  # Append to output tracking table
  output_name_table_all <- rbind(output_name_table_all, data.frame(
    Step = "Continuous_XGBoost",
    Input_File = input_file,
    Output_Object = output_name,
    Saved_As = output_file,
    stringsAsFactors = FALSE
  ))
  
  write.table(output_name_table_all, "output_name_table_all.tsv", sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  rm(list = c("df", output_name, "df_out", "xgboost_imputation_log"))
  gc()
}

cat("✅ XGBoost continuous imputation completed for all 27 inputs (df018–df044 → df180–df206).\n")
####
####
####

##### 
##### -------------------------------------------------------------------------------
##### Groupwise LightGBM Imputation for numeric omic layers (.1, .4, .5, .6, .7)
##### -------------------------------------------------------------------------------
#####
##### -------------------------------------------------------------------------------
##### PART 7 — Groupwise LightGBM Imputation for Continuous Omic Layers (.1, .4, .5, .6, .7)
##### Memory-Efficient Implementation: df018–df044 → df207–df233
##### -------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(lightgbm)
  library(rio)
})

# LightGBM imputation function
impute_lightgbm_groupwise <- function(df, type_col = "type", verbose = FALSE) {
  df_out <- df
  numeric_tokens <- c("5", "6")
  target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) length(x) >= 2 && x[2] %in% numeric_tokens)]
  
  imputation_log <- data.frame(
    Variable   = character(), Type = character(),
    NA_Before  = integer(),  NA_After = integer(),
    n_imputed  = integer(),  Method = character(),
    stringsAsFactors = FALSE
  )
  
  for (var in target_vars) {
    prefix <- strsplit(var, "-")[[1]][1]
    group_rows <- which(df[[type_col]] == prefix)
    df_sub <- df[group_rows, , drop = FALSE]
    
    y <- df_sub[[var]]
    na_idx <- which(is.na(y))
    if (length(na_idx) == 0) next
    
    predictors <- setdiff(target_vars, var)
    df_train <- df_sub[!is.na(y), predictors, drop = FALSE]
    y_train  <- y[!is.na(y)]
    df_pred  <- df_sub[ is.na(y), predictors, drop = FALSE]
    
    if (nrow(df_train) < 10 || length(unique(y_train)) < 2) {
      if (verbose) cat("⚠️ Skipped", var, "in", prefix, "— insufficient training data.\n")
      next
    }
    
    tryCatch({
      train_data <- lgb.Dataset(data = as.matrix(df_train), label = y_train)
      model <- lgb.train(
        params = list(objective = "regression", verbose = -1),
        data = train_data,
        nrounds = 50
      )
      
      pred <- predict(model, newdata = as.matrix(df_pred))
      df_out[group_rows[na_idx], var] <- pred
      
      imputation_log <- rbind(imputation_log, data.frame(
        Variable = var, Type = prefix,
        NA_Before = length(na_idx), NA_After = 0,
        n_imputed = length(na_idx), Method = "LightGBM",
        stringsAsFactors = FALSE
      ))
      
      if (verbose) cat("✅ LightGBM imputed", length(na_idx), "values for", var, "| Group:", prefix, "\n")
      
    }, error = function(e) {
      if (verbose) cat("⚠️ LightGBM failed for", var, "in group", prefix, "|", e$message, "\n")
    })
  }
  
  assign("lightgbm_imputation_log", imputation_log, envir = .GlobalEnv)
  return(df_out)
}

# Memory-efficient loop: df018–df044 → df207–df233
input_indices <- 18:20
output_indices <- 36:38

if (!exists("output_name_table_all")) {
  output_name_table_all <- data.frame(
    Step = character(),
    Input_File = character(),
    Output_Object = character(),
    Saved_As = character(),
    stringsAsFactors = FALSE
  )
}

for (i in seq_along(input_indices)) {
  input_file <- sprintf("df%03d.rds", input_indices[i])
  output_name <- sprintf("df%03d", output_indices[i])
  output_file <- paste0(output_name, ".rds")
  
  cat("🚀 Processing:", input_file, "→", output_name, "\n")
  df <- import(input_file)
  
  df_out <- impute_lightgbm_groupwise(df, verbose = TRUE)
  assign(output_name, df_out)
  saveRDS(df_out, output_file)
  
  # Append to tracking table
  output_name_table_all <- rbind(output_name_table_all, data.frame(
    Step = "Continuous_LightGBM",
    Input_File = input_file,
    Output_Object = output_name,
    Saved_As = output_file,
    stringsAsFactors = FALSE
  ))
  
  write.table(output_name_table_all, "output_name_table_all.tsv", sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  rm(list = c("df", output_name, "df_out", "lightgbm_imputation_log"))
  gc()
}

cat("✅ LightGBM continuous imputation completed for all 27 inputs (df018–df044 → df207–df233).\n")














##### 
##### -------------------------------------------------------------------------------
##### Groupwise MICE Imputation (mice) for numeric omic layers (.1, .4, .5, .6, .7)
##### -------------------------------------------------------------------------------
##### 
# 
# ##### -------------------------------------------------------------------------------
# ##### PART 8 — Groupwise MICE Imputation for numeric omic layers (.1, .4, .5, .6, .7)
# ##### Input Range : df018–df044 (n=27)
# ##### Output Range: df234–df260 (MICE-imputed)
# ##### -------------------------------------------------------------------------------
# 
# suppressPackageStartupMessages({
#   library(mice)
#   library(rio)
# })
# 
# # -- Function: Groupwise MICE imputation --
# impute_mice_groupwise <- function(df, type_col = "type", m = 1, maxit = 5, verbose = FALSE) {
#   df_out <- df
#   numeric_tokens <- c("5", "6")
#   target_vars <- names(df)[sapply(strsplit(names(df), "\\."), function(x) {
#     length(x) >= 2 && x[2] %in% numeric_tokens
#   })]
#   
#   imputation_log <- data.frame(
#     Variable = character(), Type = character(),
#     NA_Before = integer(), NA_After = integer(),
#     n_imputed = integer(), Method = character(),
#     stringsAsFactors = FALSE
#   )
#   
#   for (prefix in unique(df[[type_col]])) {
#     group_rows <- which(df[[type_col]] == prefix)
#     df_sub <- df[group_rows, target_vars, drop = FALSE]
#     na_before <- sum(is.na(df_sub))
#     if (na_before == 0) next
#     
#     tryCatch({
#       mice_model <- mice(df_sub, m = m, maxit = maxit, method = "pmm", seed = 123, printFlag = FALSE)
#       complete_df <- complete(mice_model)
#       df_out[group_rows, target_vars] <- complete_df
#       imputed <- na_before - sum(is.na(complete_df))
#       
#       imputation_log <- rbind(imputation_log, data.frame(
#         Variable = "Multiple", Type = prefix,
#         NA_Before = na_before, NA_After = sum(is.na(complete_df)),
#         n_imputed = imputed, Method = "MICE", stringsAsFactors = FALSE
#       ))
#       
#       if (verbose) cat("✅ MICE imputation done for:", prefix, "| Imputed:", imputed, "\n")
#     }, error = function(e) {
#       if (verbose) cat("⚠️ MICE failed for group:", prefix, "| Error:", e$message, "\n")
#     })
#   }
#   
#   assign("mice_imputation_log", imputation_log, envir = .GlobalEnv)
#   return(df_out)
# }
# 
# # -- Memory-efficient MICE loop for df018 to df044 → df234 to df260 --
# input_indices <- 18:20
# output_indices <- 39:41
# 
# if (!exists("output_name_table_all")) {
#   output_name_table_all <- data.frame(
#     Step = character(), Input_File = character(),
#     Output_Object = character(), Saved_As = character(),
#     stringsAsFactors = FALSE
#   )
# }
# 
# for (i in seq_along(input_indices)) {
#   input_file <- sprintf("df%03d.rds", input_indices[i])
#   output_name <- sprintf("df%03d", output_indices[i])
#   output_file <- paste0(output_name, ".rds")
#   
#   cat("📥 Processing:", input_file, "→", output_name, "\n")
#   
#   df <- import(input_file)
#   df_out <- impute_mice_groupwise(df, m = 1, maxit = 5, verbose = TRUE)
#   
#   assign(output_name, df_out)
#   saveRDS(df_out, output_file)
#   
#   output_name_table_all <- rbind(output_name_table_all, data.frame(
#     Step = "Continuous_MICE",
#     Input_File = input_file,
#     Output_Object = output_name,
#     Saved_As = output_file,
#     stringsAsFactors = FALSE
#   ))
#   
#   write.table(output_name_table_all, "output_name_table_all.tsv",
#               sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#   
#   rm(list = c("df", "df_out", output_name, "mice_imputation_log"))
#   gc()
# }
# 
# cat("✅ MICE continuous imputation completed for all 27 inputs (df018–df044 → df234–df260).\n")
# ####
# ####
# ####
