setwd("D:/ML_multioptosis/Enrique/Part 1")
# Define required packages (fallback in case dynamic detection fails)
cran_packages <- c(
  "caret", "circlize", "ComplexHeatmap", "dplyr", "ggplot2", "grid", "gridExtra",
  "kableExtra", "Matrix", "mice", "missForest", "pROC", "purrr", "reshape2",
  "rio", "rms", "survival", "survivalROC", "survcomp", "survminer",
  "tidyr", "timeROC", "UpSetR", "VIM", "xgboost", "DiagrammeR"
)

bioc_packages <- c("ComplexHeatmap", "survcomp")

github_packages <- list(
  "UCSCXenaShiny" = "openbiox/UCSCXenaShiny",
  "DiagrammeRsvg" = "rich-iannone/DiagrammeRsvg"
)

special_packages <- c("rsvg", "magick", "lightgbm")

# ---- Universal Installation Helpers ----
install_if_missing_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

install_if_missing_bioc <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg)
  }
}

install_if_missing_github <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github(repo)
  }
}

# ---- Install Declared Packages ----
invisible(lapply(cran_packages, install_if_missing_cran))
invisible(lapply(bioc_packages, install_if_missing_bioc))
invisible(mapply(install_if_missing_github, names(github_packages), github_packages))
invisible(lapply(special_packages, install_if_missing_cran))

# ---- Dynamic Package Detection from Script ----
script_path <- "harmonization_code.R"  # OK

script_lines <- readLines(script_path, warn = FALSE, encoding = "UTF-8")
script_text <- paste(script_lines, collapse = "\n")

pkg_matches <- character()

# Detect 'library()' or 'require()' calls
pkg_matches <- c(pkg_matches, unlist(regmatches(
  script_text, gregexpr("(?<=library\\(|require\\()[\"']?([a-zA-Z0-9\\.]+)[\"']?", script_text, perl = TRUE)
)))

# Also collect package names in vectors like c("pkg1", "pkg2", ...)
pkg_matches <- c(pkg_matches, unlist(regmatches(
  script_text, gregexpr("\"[a-zA-Z0-9\\.]+\"", script_text, perl = TRUE)
)))

# Clean and deduplicate
pkg_matches <- gsub("\"", "", pkg_matches)
pkg_matches <- sort(unique(pkg_matches))

# 🛡️ Exclude known problematic or non-relevant packages
excluded_pkgs <- c("plasma", "red", "text", "topics")
pkg_matches <- setdiff(pkg_matches, excluded_pkgs)
message("🛡️ Excluded problematic or non-relevant packages: ",
        paste(excluded_pkgs, collapse = ", "))


# Validate availability on system or CRAN
available_cran <- rownames(available.packages())
valid_pkgs <- pkg_matches[
  sapply(pkg_matches, function(pkg) {
    requireNamespace(pkg, quietly = TRUE) || pkg %in% available_cran
  })
]

# Install + load detected packages
load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
invisible(lapply(valid_pkgs, load_or_install))

message("✔ Detected and loaded packages:\n")
print(valid_pkgs)

# Set working directory
grep("setwd", readLines("harmonization_code.R"), value = TRUE)

if (file.exists("harmonization_code.R")) {
  grep("setwd", readLines("harmonization_code.R"), value = TRUE)
}

setwd("D:/ML_multioptosis/Enrique/Part 1")

#### ============================================================================
#### 📦 Master Wrapper Functions for rio Import/Export with NA Safety
#### ============================================================================

library(rio)

safe_import_tsv <- function(file, format = NULL, ...) {
  import(file, format = format, na.strings = "NA", ...)
}

# ---- Safe TSV writer (define once) ----
if (!exists("safe_export_tsv", mode = "function")) {
  safe_export_tsv <- function(x, file, na = "NA") {
    stopifnot(is.character(file), length(file) == 1L)
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(
      data.table::as.data.table(x),
      file = file,
      sep = "\t",
      na = na,
      quote = FALSE
    )
  }
}



#### ============================================================================
#### 📦 Installation Script for PDF/Image Processing Packages
#### ============================================================================

required_pkgs <- c("pdftools", "magick", "fs", "tidyverse")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]

if (length(new_pkgs) > 0) {
  message("Installing missing packages: ", paste(new_pkgs, collapse = ", "))
  tryCatch({
    install.packages(new_pkgs, repos = "https://cloud.r-project.org", dependencies = TRUE)
    message("Successfully installed packages")
  }, error = function(e) {
    warning("Installation failed: ", e$message)
    tryCatch({
      install.packages(new_pkgs, repos = "https://cran.rstudio.com", dependencies = TRUE)
      message("Successfully installed using backup mirror")
    }, error = function(e) {
      stop("Critical installation failure. Please check internet connection or try:\n",
           'install.packages(c("', paste(new_pkgs, collapse = '", "'), '"))')
    })
  })
} else {
  message("All required packages already installed")
}

# Verify successful loading
success <- sapply(required_pkgs, require, character.only = TRUE)
if (all(success)) {
  message("\nAll packages loaded successfully:")
  print(sessionInfo()[c("otherPkgs", "loadedOnly")])
} else {
  warning("Failed to load: ", paste(names(success)[!success], collapse = ", "))
}

# Linux-specific system requirement hints
if (Sys.info()["sysname"] == "Linux") {
  message("\nLinux system detected. You may also need to run:")
  message('sudo apt-get install libpoppler-cpp-dev libmagick++-dev')
}

message("✔ Detected and loaded packages:\n")
print(valid_pkgs)

# 👇 Add this block to show what's actually loaded
message("\n📦 Final list of packages successfully loaded in this session:\n")
loaded_pkgs <- valid_pkgs[sapply(valid_pkgs, function(pkg) pkg %in% loadedNamespaces())]
print(loaded_pkgs)



#### -----------------------------------------------------------------------------
#### Step 5: harmonize demographic data
#### -----------------------------------------------------------------------------
setwd("D:/ML_multioptosis/Enrique/Part 1")

## ==== Helpers & Operators (definir no topo do script) ====
`%||%` <- function(a, b) if (is.null(a)) b else a

norm_status <- function(x) {
  trimws(toupper(as.character(x)))
}

# ======================= Part 1 - load the data =======================
dir.create("harmonization_results")
dir.create("data_table_results")

# df005 <- import("df003_cleaned.tsv", na.strings = "NA")
df005 <- import("df005.rds")

# 📌 Define range of columns to evaluate
colunas_para_avaliar <- 1:26

# 🧼 Replace "" with NA in categorical variables
df005[colunas_para_avaliar] <- lapply(df005[colunas_para_avaliar], function(x) {
  x[x == ""] <- NA
  return(x)
})

# =======================
# Audit "" → NA per column
# =======================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)   # if using readr::write_excel_csv2 to export, for example
  library(tidyr)
  library(purrr)
})

# ---- 0) Parameters
arquivo <- "df005.rds"
colunas_para_avaliar <- 1:26

# ---- 1) Read "before" data
df_before <- import(arquivo, na.strings = "NA")

# Sanity check: ensure that the range exists
stopifnot(max(colunas_para_avaliar) <= ncol(df_before))

# Helper: summarize per column
summarize_empty_na <- function(df, cols) {
  tibble(
    col = colnames(df)[cols],
    n_empty = map_int(cols, ~{
      x <- df[[.x]]
      # Count empty strings (only for character vectors)
      if (is.character(x)) sum(x == "", na.rm = TRUE) else 0L
    }),
    n_na = map_int(cols, ~ sum(is.na(df[[.x]])))
  )
}

# ---- 2) Metrics "before"
audit_before <- summarize_empty_na(df_before, colunas_para_avaliar) %>%
  rename(n_empty_before = n_empty, n_na_before = n_na)

# ---- 3) Harmonization ("" → NA) in target columns
df_after <- df_before
df_after[colunas_para_avaliar] <- lapply(df_after[colunas_para_avaliar], function(x) {
  if (is.character(x)) {
    x[x == ""] <- NA
  }
  x
})

# ---- 4) Metrics "after"
audit_after <- summarize_empty_na(df_after, colunas_para_avaliar) %>%
  rename(n_empty_after = n_empty, n_na_after = n_na)

# ---- 5) Join and calculate deltas
audit <- audit_before %>%
  inner_join(audit_after, by = "col") %>%
  mutate(
    delta_empty = n_empty_after - n_empty_before,  # should be 0 or negative
    delta_na    = n_na_after    - n_na_before,     # should increase when there were ""
    ok_sem_vazio = n_empty_after == 0
  ) %>%
  arrange(desc(!ok_sem_vazio), desc(delta_na))

# ---- 6) Assertions/guarantees
if (any(!audit$ok_sem_vazio)) {
  problem_cols <- paste(audit$col[!audit$ok_sem_vazio], collapse = ", ")
  stop(paste0("Audit failure: there are still empty strings in the columns: ", problem_cols))
}

# Optional: also prevent delta_empty from being positive (should never increase)
if (any(audit$delta_empty > 0)) {
  stop("Audit failure: 'n_empty' increased after harmonization (unexpected behavior).")
}

# ---- 7) Result
print(audit, n = nrow(audit))

# Export audit as TSV
readr::write_tsv(audit, "harmonization_results/audit_Parte_1.tsv")

## ============================================================
## PART 2 — BATCHID IMPUTATION BY WEIGHTED RANDOM DRAW + AUDIT
##   - Stratified by "type"
##   - Draws missing BatchId from the empirical distribution of
##     observed batches within each type (proportional to counts)
##   - Preserves non-missing values; audits validity and expectations
##   - Exports audits to:
##     - harmonization_results/audit_Parte_2_BatchId_random_draw_por_type.tsv
##     - harmonization_results/audit_Parte_2_BatchId_random_draw_rowlevel.tsv
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

## -----------------------------
## 0) Output directory + seed
## -----------------------------
out_dir <- "harmonization_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Reprodutibilidade do sorteio
set.seed(123)  # ajuste/remova conforme sua política

## -----------------------------
## 1) "Before" state
## -----------------------------
stopifnot(exists("df005"))
stopifnot(all(c("type", "BatchId") %in% names(df005)))

df006_before <- df005 %>%
  mutate(
    BatchId = as.character(BatchId),
    .rowid  = dplyr::row_number()
  )

## -----------------------------
## 2) Distribuição empírica p(BatchId | type)
## -----------------------------
counts <- df006_before %>%
  filter(!is.na(BatchId)) %>%
  count(type, BatchId, name = "n_val")

probs_by_type <- counts %>%
  group_by(type) %>%
  mutate(
    total_type = sum(n_val),
    p          = n_val / total_type
  ) %>%
  ungroup()

# lookup por type para sorteio
draw_tables <- probs_by_type %>%
  group_split(type) %>%
  setNames(unique(probs_by_type$type))

draw_one_batch <- function(type_value) {
  if (!type_value %in% names(draw_tables)) return(NA_character_)
  tab <- draw_tables[[type_value]]
  p <- tab$p / sum(tab$p)              # segurança numérica
  sample(tab$BatchId, size = 1L, replace = TRUE, prob = p)
}

## -----------------------------
## 3) Imputação por random draw ponderado
## -----------------------------
df006_after <- df006_before %>%
  mutate(
    BatchId_after = if_else(
      is.na(BatchId),
      vapply(type, draw_one_batch, FUN.VALUE = character(1)),
      BatchId
    ),
    imputed_flag = is.na(BatchId) & !is.na(BatchId_after)
  )

## -----------------------------
## 4) Auditorias — integridade e expectativas
## -----------------------------
# (a) Nenhum valor não-NA pode ter sido alterado
changed_non_na <- df006_after %>%
  filter(!is.na(BatchId) & BatchId_after != BatchId)
stopifnot(nrow(changed_non_na) == 0)

# (b) Toda imputação deve pertencer ao suporte observado do 'type'
valid_support <- df006_after %>%
  filter(imputed_flag) %>%
  select(.rowid, type, BatchId_after) %>%
  left_join(
    probs_by_type %>% distinct(type, BatchId) %>% mutate(in_support = TRUE),
    by = c("type", "BatchId_after" = "BatchId")
  ) %>%
  mutate(in_support = !is.na(in_support))

if (nrow(valid_support) > 0 && any(!valid_support$in_support)) {
  bad <- valid_support %>% filter(!in_support)
  stop("Audit failure: imputations include BatchId not observed in the respective 'type'.")
}

# (c) Expectativas por estrato
before_stats <- df006_before %>%
  group_by(type) %>%
  summarise(
    n                  = n(),
    n_na_before        = sum(is.na(BatchId)),
    n_non_na_before    = n - n_na_before,
    n_unique_non_na    = n_distinct(BatchId, na.rm = TRUE),
    .groups = "drop"
  )

after_stats <- df006_after %>%
  group_by(type) %>%
  summarise(
    n_after    = n(),
    n_na_after = sum(is.na(BatchId_after)),
    n_imputed  = sum(imputed_flag),
    .groups = "drop"
  )

# (d) Distância de variação total (opcional)
dist_by_type <- {
  p_before <- df006_before %>%
    filter(!is.na(BatchId)) %>%
    count(type, BatchId, name = "n_before") %>%
    group_by(type) %>%
    mutate(p_before = n_before / sum(n_before)) %>%
    ungroup()
  
  p_after <- df006_after %>%
    filter(!is.na(BatchId_after)) %>%
    count(type, BatchId_after, name = "n_after_by_batch") %>%
    group_by(type) %>%
    mutate(p_after = n_after_by_batch / sum(n_after_by_batch)) %>%
    ungroup() %>%
    rename(BatchId = BatchId_after)
  
  full_join(p_before %>% select(type, BatchId, p_before),
            p_after  %>% select(type, BatchId, p_after),
            by = c("type", "BatchId")) %>%
    mutate(
      p_before = replace_na(p_before, 0),
      p_after  = replace_na(p_after,  0),
      abs_diff = abs(p_after - p_before)
    ) %>%
    group_by(type) %>%
    summarise(TVD = 0.5 * sum(abs_diff), .groups = "drop")
}

## -----------------------------
## 5) Auditoria — NÍVEL DE LINHA (qual BatchId foi imputado)
## -----------------------------
# Construir suporte e probabilidade por 'type' para join
support_by_type <- probs_by_type %>%
  arrange(type, desc(p), BatchId) %>%
  group_by(type) %>%
  summarise(
    support      = paste(BatchId, collapse = " | "),
    .groups = "drop"
  )

# Auditoria row-level: apenas linhas imputadas
audit_rowlevel <- df006_after %>%
  filter(imputed_flag) %>%
  transmute(
    .rowid,
    type,
    BatchId_before = NA_character_,          # explicitamos que era NA
    BatchId_after,
    # probabilidade empírica do valor sorteado dentro do 'type'
    # (juntamos com probs_by_type para obter 'p')
    # e anexamos o suporte do estrato
  ) %>%
  left_join(
    probs_by_type %>% select(type, BatchId, p),
    by = c("type", "BatchId_after" = "BatchId")
  ) %>%
  left_join(support_by_type, by = "type") %>%
  arrange(type, desc(p), .rowid)

# Opcional: incluir um identificador de amostra, se existir na base
id_cols <- c("sample", "Sample", "patient_id", "case_id", "barcode")
id_col_present <- id_cols[id_cols %in% names(df006_before)]
if (length(id_col_present) > 0) {
  id_col <- id_col_present[1]
  audit_rowlevel <- df006_after %>%
    filter(imputed_flag) %>%
    select(.rowid, type, !!sym(id_col)) %>%
    right_join(audit_rowlevel, by = c(".rowid","type")) %>%
    relocate(!!sym(id_col), .after = .rowid)
}

## -----------------------------
## 6) Consolidação e exportação das auditorias
## -----------------------------
audit_summary <- before_stats %>%
  left_join(after_stats,  by = "type") %>%
  left_join(
    probs_by_type %>%
      group_by(type) %>%
      summarise(
        mode_like = BatchId[which.max(p)],
        max_p     = max(p),
        support   = paste(sort(unique(BatchId)), collapse = " | "),
        .groups = "drop"
      ),
    by = "type"
  ) %>%
  left_join(dist_by_type, by = "type") %>%
  mutate(
    expectation_met =
      case_when(
        n_non_na_before > 0 ~ n_na_after == 0L,
        n_non_na_before == 0 ~ n_na_after == n
      )
  ) %>%
  select(type, n, n_unique_non_na,
         n_na_before, n_non_na_before,
         n_na_after, n_imputed,
         mode_like, max_p, support, TVD, expectation_met) %>%
  arrange(desc(!expectation_met), desc(n_imputed), type)

# Exceção se expectativas falharem
if (any(!audit_summary$expectation_met, na.rm = TRUE)) {
  prob_types <- paste(audit_summary$type[isFALSE(audit_summary$expectation_met)], collapse = ", ")
  stop(paste0(
    "Audit failure: expectation not met for types: ", prob_types,
    ". Check if there were enough non-NA values to define the empirical distribution or inconsistencies in 'type'."
  ))
}

# Exportações
audit_path_summary <- file.path(out_dir, "audit_Parte_2_BatchId_random_draw_por_type.tsv")
audit_path_rowlvl  <- file.path(out_dir, "audit_Parte_2_BatchId_random_draw_rowlevel.tsv")

readr::write_tsv(audit_summary, audit_path_summary)
readr::write_tsv(audit_rowlevel, audit_path_rowlvl)

## -----------------------------
## 7) Objeto final
## -----------------------------
df006 <- df006_after %>%
  mutate(BatchId = BatchId_after) %>%
  select(-BatchId_after, -imputed_flag, -.rowid)

rio::export(df006,"data_table_results/df006.tsv")

cat("✓ Weighted random draw imputation audit completed.\n")
cat("→ Summary report:  ", normalizePath(audit_path_summary), "\n", sep = "")
cat("→ Row-level report:", normalizePath(audit_path_rowlvl),  "\n", sep = "")

## ============================================================
## PART 3 — HARMONIZATION OF AJCC PATHOLOGIC TUMOR STAGE + AUDIT
##   Rules (in the order applied by your code):
##   R1) If stage is NA, OS==0, DSS==0, DFI==0, PFI==0 and tumor_status=="TUMOR FREE"  -> "T0"
##   R2) If stage is NA, OS==0, DSS==0, DFI==0, PFI==0 and tumor_status=="WITH TUMOR" -> "Missing"
##   R3) If stage is NA and new_tumor_event_type contains "Distant Metastasis|Metastatic" -> "Stage IV"
##   (R1 and R2 precede R3 because you fill first, and R3 only runs on remaining NAs)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(stringr)
})

out_dir <- "harmonization_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## -----------------------------
## 0) Preconditions
## -----------------------------
stopifnot(exists("df006"))
req <- c("ajcc_pathologic_tumor_stage","OS","DSS","DFI","PFI","tumor_status","new_tumor_event_type")
stopifnot(all(req %in% names(df006)))

## -----------------------------
## 1) "Before" state
## -----------------------------
df007_before <- df006 %>%
  mutate(
    ajcc_pathologic_tumor_stage = as.character(ajcc_pathologic_tumor_stage),
    tumor_status                = as.character(tumor_status),
    new_tumor_event_type        = as.character(new_tumor_event_type),
    .rowid = dplyr::row_number()
  )

stage_before <- df007_before$ajcc_pathologic_tumor_stage

## -----------------------------
## 2) Rules — flags computed on the "before" state
##     (without normalizations: literal equality for tumor_status,
##      regex case-insensitive for metastasis)
## -----------------------------
r1_flag <- is.na(df007_before$ajcc_pathologic_tumor_stage) &
  df007_before$OS  == 0 &
  df007_before$DSS == 0 &
  df007_before$DFI == 0 &
  df007_before$PFI == 0 &
  df007_before$tumor_status == "TUMOR FREE"

r2_flag <- is.na(df007_before$ajcc_pathologic_tumor_stage) &
  df007_before$OS  == 0 &
  df007_before$DSS == 0 &
  df007_before$DFI == 0 &
  df007_before$PFI == 0 &
  df007_before$tumor_status == "WITH TUMOR"

r3_flag <- is.na(df007_before$ajcc_pathologic_tumor_stage) &
  !is.na(df007_before$new_tumor_event_type) &
  str_detect(df007_before$new_tumor_event_type, regex("Distant Metastasis|Metastatic", ignore_case = TRUE))

## Flag for multiple concurrent triggers
multi_hit <- (r1_flag | r2_flag) & r3_flag

## -----------------------------
## 3) EXPECTED result (applying rules in the same order as your code)
## -----------------------------
stage_expected <- case_when(
  r1_flag ~ "T0",
  r2_flag ~ "Missing",
  r3_flag ~ "Stage IV",
  TRUE    ~ stage_before
)

## -----------------------------
## 4) Application of your BLOCK (OBSERVED result)
## -----------------------------
df007_after <- df007_before %>%
  mutate(
    ajcc_pathologic_tumor_stage = if_else(
      is.na(ajcc_pathologic_tumor_stage) &
        OS == 0 & DSS == 0 & DFI == 0 & PFI == 0 &
        tumor_status == "TUMOR FREE",
      "T0",
      ajcc_pathologic_tumor_stage
    )
  ) %>%
  mutate(
    ajcc_pathologic_tumor_stage = if_else(
      is.na(ajcc_pathologic_tumor_stage) &
        OS == 0 & DSS == 0 & DFI == 0 & PFI == 0 &
        tumor_status == "WITH TUMOR",
      "Missing",
      ajcc_pathologic_tumor_stage
    )
  ) %>%
  mutate(
    ajcc_pathologic_tumor_stage = case_when(
      is.na(ajcc_pathologic_tumor_stage) &
        !is.na(new_tumor_event_type) &
        str_detect(new_tumor_event_type, regex("Distant Metastasis|Metastatic", ignore_case = TRUE)) ~ "Stage IV",
      TRUE ~ ajcc_pathologic_tumor_stage
    )
  )

stage_after <- df007_after$ajcc_pathologic_tumor_stage

## -----------------------------
## 5) Audit — formal checks
## -----------------------------

# (a) No non-NA value should have been altered
changed_non_na <- which(!is.na(stage_before) & stage_before != stage_after)
if (length(changed_non_na) > 0) {
  stop("Audit failure: non-NA values of 'ajcc_pathologic_tumor_stage' were altered in ",
       length(changed_non_na), " row(s). Examples (rowid): ",
       paste(head(df007_before$.rowid[changed_non_na], 10), collapse = ", "))
}

# (b) Values introduced by the block must belong to the allowed set
introduced <- unique(stage_after[is.na(stage_before) & !is.na(stage_after)])
allowed   <- c("T0","Missing","Stage IV")
if (length(setdiff(introduced, allowed)) > 0) {
  stop("Audit failure: labels outside the allowed set were introduced: ",
       paste(setdiff(introduced, allowed), collapse = ", "))
}

# (c) Observed result == expected result (row by row)
mismatch_idx <- which(stage_after != stage_expected | xor(is.na(stage_after), is.na(stage_expected)))
if (length(mismatch_idx) > 0) {
  stop("Audit failure: ", length(mismatch_idx),
       " row(s) do not match between expected and observed. Examples (rowid): ",
       paste(head(df007_before$.rowid[mismatch_idx], 10), collapse = ", "))
}

## -----------------------------
## 6) Summaries and reports
## -----------------------------
# Effective changes
changed_idx <- which(stage_before != stage_after | xor(is.na(stage_before), is.na(stage_after)))

audit_rowlevel <- tibble(
  rowid                 = df007_before$.rowid[changed_idx],
  type                  = df007_before$type[changed_idx] %||% NA_character_,
  stage_before          = stage_before[changed_idx],
  stage_after           = stage_after[changed_idx],
  R1_T0                 = r1_flag[changed_idx],
  R2_Missing            = r2_flag[changed_idx],
  R3_StageIV            = r3_flag[changed_idx],
  multi_hit_R1R2_vs_R3  = multi_hit[changed_idx],
  new_tumor_event_type  = df007_before$new_tumor_event_type[changed_idx],
  tumor_status          = df007_before$tumor_status[changed_idx]
)

# Global and per-type summaries
audit_summary <- tibble(
  total_rows               = nrow(df007_before),
  n_stage_before_NA        = sum(is.na(stage_before)),
  n_changed                = length(changed_idx),
  n_R1_T0_applied          = sum(r1_flag),
  n_R2_Missing_applied     = sum(r2_flag),
  n_R3_StageIV_applied     = sum(r3_flag),
  n_multi_hit              = sum(multi_hit)
)

audit_by_type <- audit_rowlevel %>%
  mutate(rule = case_when(
    R1_T0 ~ "R1_T0",
    R2_Missing ~ "R2_Missing",
    R3_StageIV ~ "R3_StageIV",
    TRUE ~ "Other"   # should not occur
  )) %>%
  count(type, rule, name = "n_applied") %>%
  arrange(type, desc(n_applied))

# Exports
readr::write_tsv(audit_rowlevel, file.path(out_dir, "audit_Parte_3_ajcc_stage_rowlevel.tsv"))
readr::write_tsv(audit_by_type,  file.path(out_dir, "audit_Parte_3_ajcc_stage_by_type.tsv"))
readr::write_tsv(audit_summary,  file.path(out_dir, "audit_Parte_3_ajcc_stage_summary.tsv"))


## -----------------------------
## 7) Update working object and log
## -----------------------------
df007 <- df007_after %>% select(-.rowid)
cat("✓ Part 3 audit completed.\n")
cat("→ Row-level: ", normalizePath(file.path(out_dir, "audit_Parte_3_ajcc_stage_rowlevel.tsv")), "\n", sep = "")
cat("→ By type:   ", normalizePath(file.path(out_dir, "audit_Parte_3_ajcc_stage_by_type.tsv")),  "\n", sep = "")
cat("→ Summary:   ", normalizePath(file.path(out_dir, "audit_Parte_3_ajcc_stage_summary.tsv")),  "\n", sep = "")

rio::export(df007,"data_table_results/df007.tsv")

## ============================================================
## PART 4 — HARMONIZATION OF OS / vital_status + AUDIT
##   Rules replicated exactly from your code:
##   R1) If vital_status=="Alive" and OS is NA → OS := 0
##   R2) If vital_status=="Dead"  and OS is NA → OS := 1
##   R3) Then: vital_status := "Dead"  if OS==1;
##             vital_status := "Alive" if OS==0;
##             otherwise keep previous value.
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

out_dir <- "harmonization_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## -----------------------------
## 0) Preconditions
## -----------------------------
stopifnot(exists("df007"))
req <- c("OS","vital_status")
stopifnot(all(req %in% names(df007)))

## -----------------------------
## 1) "Before" state
## -----------------------------
df008_before <- df007 %>%
  mutate(
    # do not force type here to avoid altering your logic;
    # just ensure strings in vital_status
    vital_status = as.character(vital_status),
    .rowid = dplyr::row_number()
  )

OS_before          <- df008_before$OS
vital_before       <- df008_before$vital_status
has_type_column    <- "type" %in% names(df008_before)

## -----------------------------
## 2) EXPECTED — apply rules exactly as in your code
## -----------------------------
OS_expected <- ifelse(df008_before$vital_status == "Alive" & is.na(df008_before$OS), 0,
                      ifelse(df008_before$vital_status == "Dead"  & is.na(df008_before$OS), 1,
                             df008_before$OS))

vital_expected <- ifelse(OS_expected == 1, "Dead",
                         ifelse(OS_expected == 0, "Alive",
                                df008_before$vital_status))

## -----------------------------
## 3) OBSERVED — apply your BLOCK for comparison
## -----------------------------
df008_after <- df008_before

# R1/R2: imputation of OS only when NA based on vital_status
df008_after$OS <- ifelse(df008_after$vital_status == "Alive" & is.na(df008_after$OS), 0,
                         ifelse(df008_after$vital_status == "Dead"  & is.na(df008_after$OS), 1,
                                df008_after$OS))

# R3: map vital_status from OS
df008_after$vital_status <- ifelse(df008_after$OS == 1, "Dead",
                                   ifelse(df008_after$OS == 0, "Alive",
                                          df008_after$vital_status))

OS_after    <- df008_after$OS
vital_after <- df008_after$vital_status

## -----------------------------
## 4) Audit — formal checks
## -----------------------------

# (a) Non-NA OS values should not have been altered
changed_non_na_OS_idx <- which(!is.na(OS_before) & !is.na(OS_after) & (OS_before != OS_after))
if (length(changed_non_na_OS_idx) > 0) {
  stop("Audit failure: non-NA values of 'OS' were altered in ",
       length(changed_non_na_OS_idx), " row(s). Examples (rowid): ",
       paste(head(df008_before$.rowid[changed_non_na_OS_idx], 10), collapse = ", "))
}

# (b) Final consistency: OS==1 → vital_status=='Dead'; OS==0 → 'Alive'
inconsistent_idx <- which( (OS_after == 1 & vital_after != "Dead") |
                             (OS_after == 0 & vital_after != "Alive") )
if (length(inconsistent_idx) > 0) {
  stop("Audit failure: final inconsistency between OS and vital_status in ",
       length(inconsistent_idx), " row(s). Examples (rowid): ",
       paste(head(df008_before$.rowid[inconsistent_idx], 10), collapse = ", "))
}

# (c) Observed == Expected (row by row)
mismatch_idx <- which( ( !(is.na(OS_after)    & is.na(OS_expected))    & OS_after    != OS_expected ) |
                         ( !(is.na(vital_after) & is.na(vital_expected)) & vital_after != vital_expected ) )
if (length(mismatch_idx) > 0) {
  stop("Audit failure: ", length(mismatch_idx),
       " row(s) differ between observed and expected. Examples (rowid): ",
       paste(head(df008_before$.rowid[mismatch_idx], 10), collapse = ", "))
}

# (d) Domain of OS after harmonization (optional but recommended)
vals_OS <- unique(OS_after[!is.na(OS_after)])
invalid_OS_vals <- setdiff(vals_OS, c(0, 1))
if (length(invalid_OS_vals) > 0) {
  warning("Warning: values of OS outside {0,1,NA} were found after harmonization: ",
          paste(invalid_OS_vals, collapse = ", "),
          ". Your block does not change original non-NA values; check their source.")
}

## -----------------------------
## 5) Report construction
## -----------------------------
# Rule flags applied (computed on 'before')
r1_flag <- is.na(OS_before) & vital_before == "Alive"   # OS := 0
r2_flag <- is.na(OS_before) & vital_before == "Dead"    # OS := 1
# "Map-only": OS already 0/1 and vital was adjusted (if needed)
map_only_flag <- !r1_flag & !r2_flag & (
  (OS_before == 0 & vital_before != "Alive") |
    (OS_before == 1 & vital_before != "Dead")
)

# Which rows changed (OS and/or vital_status)
os_changed_flag    <- ( (is.na(OS_before) != is.na(OS_after)) | (!is.na(OS_before) & !is.na(OS_after) & OS_before != OS_after) )
vital_changed_flag <- ( (is.na(vital_before) != is.na(vital_after)) | (!is.na(vital_before) & !is.na(vital_after) & vital_before != vital_after) )
changed_idx <- which(os_changed_flag | vital_changed_flag)

audit_rowlevel <- tibble(
  rowid            = df008_before$.rowid[changed_idx],
  type             = if (has_type_column) df008_before$type[changed_idx] else NA_character_,
  OS_before        = OS_before[changed_idx],
  vital_before     = vital_before[changed_idx],
  OS_after         = OS_after[changed_idx],
  vital_after      = vital_after[changed_idx],
  R1_OS_from_Alive = r1_flag[changed_idx],
  R2_OS_from_Dead  = r2_flag[changed_idx],
  MapVital_only    = map_only_flag[changed_idx]
)

# Global summary
audit_summary <- tibble(
  total_rows           = nrow(df008_before),
  n_OS_NA_before       = sum(is.na(OS_before)),
  n_OS_NA_after        = sum(is.na(OS_after)),
  n_OS_imputed_R1      = sum(r1_flag),  # OS := 0
  n_OS_imputed_R2      = sum(r2_flag),  # OS := 1
  n_vital_mapped_only  = sum(map_only_flag),
  n_rows_changed       = length(changed_idx)
)

# By type (if exists)
if (has_type_column) {
  audit_by_type <- audit_rowlevel %>%
    mutate(rule = case_when(
      R1_OS_from_Alive ~ "R1_OS:=0_from_Alive",
      R2_OS_from_Dead  ~ "R2_OS:=1_from_Dead",
      MapVital_only    ~ "R3_MapVital_only",
      TRUE             ~ "Other" # should not occur
    )) %>%
    count(type, rule, name = "n_applied") %>%
    arrange(type, desc(n_applied))
} else {
  audit_by_type <- tibble(
    info = "Column 'type' does not exist; report by type not generated."
  )
}

## -----------------------------
## 6) Exports and object update
## -----------------------------
readr::write_tsv(audit_rowlevel, file.path(out_dir, "audit_Parte_4_OS_vital_rowlevel.tsv"))
readr::write_tsv(audit_summary,  file.path(out_dir, "audit_Parte_4_OS_vital_summary.tsv"))
readr::write_tsv(audit_by_type,  file.path(out_dir, "audit_Parte_4_OS_vital_by_type.tsv"))

# Update df008 according to your transformation
df008 <- df008_after %>% select(-.rowid)

cat("✓ Part 4 audit completed.\n")
cat("→ Row-level: ", normalizePath(file.path(out_dir, "audit_Parte_4_OS_vital_rowlevel.tsv")), "\n", sep = "")
cat("→ Summary:   ", normalizePath(file.path(out_dir, "audit_Parte_4_OS_vital_summary.tsv")),  "\n", sep = "")
cat("→ By type:   ", normalizePath(file.path(out_dir, "audit_Parte_4_OS_vital_by_type.tsv")),  "\n", sep = "")

rio::export(df008,"data_table_results/df008.tsv")

## ============================================================
## PART 5 — HARMONIZATION OF DFI & DFI.cr + AUDIT
##   Regra aplicada (exata, case-sensitive):
##   Se is.na(Evento) & new_tumor_event_type %in% valid_events → Evento := 1
##   Caso contrário, manter valor original.
##   Observação: DFI ∈ {0,1,NA}; DFI.cr ∈ {0,1,2,NA}. Aqui só imputamos "1".
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(stringr)
})

out_dir <- "harmonization_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## -----------------------------
## 0) Preconditions
## -----------------------------
stopifnot(exists("df008"))
req_base <- c("new_tumor_event_type","tumor_status","vital_status")
stopifnot(all(req_base %in% names(df008)))

has_type_column <- "type" %in% names(df008)
has_DFI         <- "DFI"     %in% names(df008)
has_DFIcr       <- "DFI.cr"  %in% names(df008)
has_time_DFI    <- "DFI.time"     %in% names(df008)
has_time_DFIcr  <- "DFI.time.cr"  %in% names(df008)

if (!has_DFI & !has_DFIcr) {
  stop("Nenhuma coluna de endpoint encontrada: é necessário ao menos 'DFI' ou 'DFI.cr'.")
}

## -----------------------------
## 1) "Before" state + vector of valid events
## -----------------------------
df009_before <- df008 %>%
  mutate(
    new_tumor_event_type = as.character(new_tumor_event_type),
    .rowid = dplyr::row_number()
  )

valid_events <- c(
  "Recurrence",
  "Distant Metastasis",
  "Locoregional Recurrence",
  "New Primary Tumor",
  "Biochemical evidence of disease",
  "Progression of Disease"
)

## Utilitários de auditoria e export
write_if_rows <- function(tbl, path) {
  if (!is.null(tbl) && nrow(tbl) > 0) readr::write_tsv(tbl, path)
}

## ============================================================
## 2) BLOCO — DFI (binário) — membership exato NA→1
## ============================================================
if (has_DFI) {
  DFI_before   <- df009_before$DFI
  event_before <- df009_before$new_tumor_event_type
  time_before  <- if (has_time_DFI) df009_before$DFI.time else rep(NA_real_, nrow(df009_before))
  
  ## 2.1 EXPECTED
  candidate_flag_DFI <- is.na(DFI_before) & !is.na(event_before) & (event_before %in% valid_events)
  DFI_expected       <- ifelse(candidate_flag_DFI, 1, DFI_before)
  
  ## 2.2 OBSERVED
  df009_after_DFI <- df009_before %>%
    mutate(
      DFI = case_when(
        is.na(DFI) & new_tumor_event_type %in% valid_events ~ 1,
        TRUE ~ DFI
      )
    )
  
  DFI_after  <- df009_after_DFI$DFI
  time_after <- if (has_time_DFI) df009_after_DFI$DFI.time else rep(NA_real_, nrow(df009_after_DFI))
  
  ## 2.3 Auditorias
  # (a) não alterar não-NA
  changed_non_na_idx <- which(!is.na(DFI_before) & !is.na(DFI_after) & (DFI_before != DFI_after))
  if (length(changed_non_na_idx) > 0) {
    stop("Audit failure (DFI): non-NA values were altered in ",
         length(changed_non_na_idx), " row(s). Examples (rowid): ",
         paste(head(df009_before$.rowid[changed_non_na_idx], 10), collapse = ", "))
  }
  
  # (b) Observed == Expected
  mismatch_idx <- which( ( !(is.na(DFI_after) & is.na(DFI_expected)) ) & (DFI_after != DFI_expected) )
  if (length(mismatch_idx) > 0) {
    stop("Audit failure (DFI): ", length(mismatch_idx),
         " row(s) differ between observed and expected. Examples (rowid): ",
         paste(head(df009_before$.rowid[mismatch_idx], 10), collapse = ", "))
  }
  
  # (c) domínio
  vals_DFI <- unique(DFI_after[!is.na(DFI_after)])
  invalid_vals <- setdiff(vals_DFI, c(0,1))
  if (length(invalid_vals) > 0) {
    warning("Domain warning (DFI): outside {0,1,NA}: ", paste(invalid_vals, collapse = ", "))
  }
  
  # (d) origem das imputações
  imputed_idx <- which(is.na(DFI_before) & !is.na(DFI_after))
  if (length(imputed_idx) > 0) {
    bad_origin <- imputed_idx[ !(event_before[imputed_idx] %in% valid_events) | (DFI_after[imputed_idx] != 1) ]
    if (length(bad_origin) > 0) {
      stop("Audit failure (DFI): imputations not from valid_events OR imputed value != 1. Examples (rowid): ",
           paste(head(df009_before$.rowid[bad_origin], 10), collapse = ", "))
    }
  }
  
  ## 2.4 Relatórios DFI
  imputed_flag <- is.na(DFI_before) & !is.na(DFI_after) & (DFI_after == 1)
  changed_idx  <- which(imputed_flag)
  
  audit_rowlevel_DFI <- tibble(
    rowid                = df009_before$.rowid[changed_idx],
    type                 = if (has_type_column) df009_before$type[changed_idx] else NA_character_,
    DFI_before           = DFI_before[changed_idx],
    DFI_after            = DFI_after[changed_idx],
    new_tumor_event_type = event_before[changed_idx],
    candidate_rule       = candidate_flag_DFI[changed_idx],
    imputed              = imputed_flag[changed_idx]
  )
  
  audit_summary_DFI <- tibble(
    total_rows        = nrow(df009_before),
    n_DFI_NA_before   = sum(is.na(DFI_before)),
    n_candidates_rule = sum(candidate_flag_DFI),
    n_imputed         = sum(imputed_flag),
    n_DFI_NA_after    = sum(is.na(DFI_after)),
    n_rows_changed    = length(changed_idx)
  ) %>% mutate(check_balance = (n_DFI_NA_before - n_imputed) == n_DFI_NA_after)
  
  if (has_type_column) {
    audit_by_type_DFI <- tibble(
      type            = df009_before$type,
      candidate_rule  = candidate_flag_DFI,
      imputed         = imputed_flag
    ) %>%
      group_by(type) %>%
      summarise(
        n_candidates_rule = sum(candidate_rule),
        n_imputed         = sum(imputed),
        .groups = "drop"
      ) %>%
      arrange(desc(n_imputed), desc(n_candidates_rule), type)
  } else {
    audit_by_type_DFI <- tibble(info = "Column 'type' does not exist; report by type not generated.")
  }
  
  lower_valid <- tolower(trimws(valid_events))
  diag_idx <- which(is.na(DFI_before) & !is.na(event_before) &
                      tolower(trimws(event_before)) %in% lower_valid &
                      !(event_before %in% valid_events))
  diag_case_mismatch_DFI <- tibble(
    rowid = df009_before$.rowid[diag_idx],
    type  = if (has_type_column) df009_before$type[diag_idx] else NA_character_,
    event_original   = event_before[diag_idx],
    event_normalized = tolower(trimws(event_before[diag_idx]))
  )
  
  freq_DFI_event_type <- df009_after_DFI %>%
    group_by(DFI, new_tumor_event_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(DFI, desc(n))
  
  freq_DFI_tumor_vital <- df009_after_DFI %>%
    group_by(DFI, tumor_status, vital_status) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n))
  
  freq_pivot_DFI <- df009_after_DFI %>%
    filter(!is.na(DFI)) %>%
    count(DFI, tumor_status, new_tumor_event_type) %>%
    pivot_wider(names_from = DFI, values_from = n, values_fill = 0)
  
  # Evento × Tempo — binário
  if (has_time_DFI) {
    time_check_DFI <- tibble(
      event_col = "DFI",
      time_col  = "DFI.time",
      n_event_1_time_NA = sum(!is.na(DFI_after) & DFI_after %in% c(1L,"1") & is.na(time_after)),
      n_event_0_time_pos= sum(!is.na(DFI_after) & DFI_after %in% c(0L,"0") & !is.na(time_after) & time_after > 0),
      n_time_neg        = sum(!is.na(time_after) & time_after < 0)
    )
    readr::write_tsv(time_check_DFI, file.path(out_dir, "audit_Parte_5_DFI_time_check.tsv"))
  }
  
  ## 2.5 Atualiza objeto de trabalho com DFI harmonizado
  df009_mid <- df009_after_DFI %>% select(-.rowid)
  
  ## 2.6 Exports DFI
  write_if_rows(audit_rowlevel_DFI,   file.path(out_dir, "audit_Parte_5_DFI_rowlevel.tsv"))
  write_if_rows(audit_summary_DFI,    file.path(out_dir, "audit_Parte_5_DFI_summary.tsv"))
  write_if_rows(audit_by_type_DFI,    file.path(out_dir, "audit_Parte_5_DFI_by_type.tsv"))
  write_if_rows(diag_case_mismatch_DFI, file.path(out_dir, "audit_Parte_5_DFI_case_mismatch.tsv"))
  readr::write_tsv(freq_DFI_event_type, file.path(out_dir, "freq_Parte_5_DFI_event_type.tsv"))
  readr::write_tsv(freq_DFI_tumor_vital, file.path(out_dir, "freq_Parte_5_DFI_tumor_vital.tsv"))
  readr::write_tsv(freq_pivot_DFI,      file.path(out_dir, "freq_Parte_5_DFI_pivot.tsv"))
  
  cat("✓ Part 5 (DFI) audit completed.\n")
} else {
  df009_mid <- df009_before %>% select(-.rowid)
}

## ============================================================
## 3) BLOCO — DFI.cr (competitivo) — membership exato NA→1
## ============================================================
if (has_DFIcr) {
  # Recarrega .rowid (pode ter sido removido)
  df009_before2 <- df009_mid %>%
    mutate(.rowid = dplyr::row_number())
  
  DFIcr_before   <- df009_before2$`DFI.cr`
  event_before2  <- df009_before2$new_tumor_event_type
  timecr_before  <- if (has_time_DFIcr) df009_before2$`DFI.time.cr` else rep(NA_real_, nrow(df009_before2))
  
  ## 3.1 EXPECTED
  candidate_flag_DFIcr <- is.na(DFIcr_before) & !is.na(event_before2) & (event_before2 %in% valid_events)
  DFIcr_expected       <- ifelse(candidate_flag_DFIcr, 1, DFIcr_before)
  
  ## 3.2 OBSERVED
  df009_after_DFIcr <- df009_before2 %>%
    mutate(
      `DFI.cr` = case_when(
        is.na(`DFI.cr`) & new_tumor_event_type %in% valid_events ~ 1,
        TRUE ~ `DFI.cr`
      )
    )
  
  DFIcr_after  <- df009_after_DFIcr$`DFI.cr`
  timecr_after <- if (has_time_DFIcr) df009_after_DFIcr$`DFI.time.cr` else rep(NA_real_, nrow(df009_after_DFIcr))
  
  ## 3.3 Auditorias
  # (a) não alterar não-NA
  changed_non_na_idx <- which(!is.na(DFIcr_before) & !is.na(DFIcr_after) & (DFIcr_before != DFIcr_after))
  if (length(changed_non_na_idx) > 0) {
    stop("Audit failure (DFI.cr): non-NA values were altered in ",
         length(changed_non_na_idx), " row(s). Examples (rowid): ",
         paste(head(df009_before2$.rowid[changed_non_na_idx], 10), collapse = ", "))
  }
  
  # (b) Observed == Expected
  mismatch_idx <- which( ( !(is.na(DFIcr_after) & is.na(DFIcr_expected)) ) & (DFIcr_after != DFIcr_expected) )
  if (length(mismatch_idx) > 0) {
    stop("Audit failure (DFI.cr): ", length(mismatch_idx),
         " row(s) differ between observed and expected. Examples (rowid): ",
         paste(head(df009_before2$.rowid[mismatch_idx], 10), collapse = ", "))
  }
  
  # (c) domínio {0,1,2}
  vals_DFIcr <- unique(DFIcr_after[!is.na(DFIcr_after)])
  invalid_vals <- setdiff(vals_DFIcr, c(0,1,2))
  if (length(invalid_vals) > 0) {
    warning("Domain warning (DFI.cr): outside {0,1,2,NA}: ", paste(invalid_vals, collapse = ", "))
  }
  
  # (d) origem das imputações
  imputed_idx <- which(is.na(DFIcr_before) & !is.na(DFIcr_after))
  if (length(imputed_idx) > 0) {
    bad_origin <- imputed_idx[ !(event_before2[imputed_idx] %in% valid_events) | (DFIcr_after[imputed_idx] != 1) ]
    if (length(bad_origin) > 0) {
      stop("Audit failure (DFI.cr): imputations not from valid_events OR imputed value != 1. Examples (rowid): ",
           paste(head(df009_before2$.rowid[bad_origin], 10), collapse = ", "))
    }
  }
  
  ## 3.4 Relatórios DFI.cr
  imputed_flag_cr <- is.na(DFIcr_before) & !is.na(DFIcr_after) & (DFIcr_after == 1)
  changed_idx_cr  <- which(imputed_flag_cr)
  
  audit_rowlevel_DFIcr <- tibble(
    rowid                = df009_before2$.rowid[changed_idx_cr],
    type                 = if (has_type_column) df009_before2$type[changed_idx_cr] else NA_character_,
    DFIcr_before         = DFIcr_before[changed_idx_cr],
    DFIcr_after          = DFIcr_after[changed_idx_cr],
    new_tumor_event_type = event_before2[changed_idx_cr],
    candidate_rule       = candidate_flag_DFIcr[changed_idx_cr],
    imputed              = imputed_flag_cr[changed_idx_cr]
  )
  
  audit_summary_DFIcr <- tibble(
    total_rows          = nrow(df009_before2),
    n_DFIcr_NA_before   = sum(is.na(DFIcr_before)),
    n_candidates_rule   = sum(candidate_flag_DFIcr),
    n_imputed           = sum(imputed_flag_cr),
    n_DFIcr_NA_after    = sum(is.na(DFIcr_after)),
    n_rows_changed      = length(changed_idx_cr)
  ) %>% mutate(check_balance = (n_DFIcr_NA_before - n_imputed) == n_DFIcr_NA_after)
  
  if (has_type_column) {
    audit_by_type_DFIcr <- tibble(
      type            = df009_before2$type,
      candidate_rule  = candidate_flag_DFIcr,
      imputed         = imputed_flag_cr
    ) %>%
      group_by(type) %>%
      summarise(
        n_candidates_rule = sum(candidate_rule),
        n_imputed         = sum(imputed),
        .groups = "drop"
      ) %>%
      arrange(desc(n_imputed), desc(n_candidates_rule), type)
  } else {
    audit_by_type_DFIcr <- tibble(info = "Column 'type' does not exist; report by type not generated.")
  }
  
  lower_valid <- tolower(trimws(valid_events))
  diag_idx_cr <- which(is.na(DFIcr_before) & !is.na(event_before2) &
                         tolower(trimws(event_before2)) %in% lower_valid &
                         !(event_before2 %in% valid_events))
  diag_case_mismatch_DFIcr <- tibble(
    rowid = df009_before2$.rowid[diag_idx_cr],
    type  = if (has_type_column) df009_before2$type[diag_idx_cr] else NA_character_,
    event_original   = event_before2[diag_idx_cr],
    event_normalized = tolower(trimws(event_before2[diag_idx_cr]))
  )
  
  freq_DFIcr_event_type <- df009_after_DFIcr %>%
    group_by(`DFI.cr`, new_tumor_event_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(`DFI.cr`, desc(n))
  
  freq_DFIcr_tumor_vital <- df009_after_DFIcr %>%
    group_by(`DFI.cr`, tumor_status, vital_status) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n))
  
  freq_pivot_DFIcr <- df009_after_DFIcr %>%
    filter(!is.na(`DFI.cr`)) %>%
    count(`DFI.cr`, tumor_status, new_tumor_event_type) %>%
    pivot_wider(names_from = `DFI.cr`, values_from = n, values_fill = 0)
  
  # Evento × Tempo — competitivo
  if (has_time_DFIcr) {
    time_check_DFIcr <- tibble(
      event_col = "DFI.cr",
      time_col  = "DFI.time.cr",
      n_event_12_time_NA = sum(!is.na(DFIcr_after) & DFIcr_after %in% c(1L,2L) & is.na(timecr_after)),
      n_event_0_time_pos = sum(!is.na(DFIcr_after) & DFIcr_after %in% c(0L)   & !is.na(timecr_after) & timecr_after > 0),
      n_time_neg         = sum(!is.na(timecr_after) & timecr_after < 0)
    )
    readr::write_tsv(time_check_DFIcr, file.path(out_dir, "audit_Parte_5_DFIcr_time_check.tsv"))
  }
  
  ## 3.5 Atualiza objeto final com DFI.cr harmonizado
  df009 <- df009_after_DFIcr %>% select(-.rowid)
  
  ## 3.6 Exports DFI.cr
  write_if_rows(audit_rowlevel_DFIcr,   file.path(out_dir, "audit_Parte_5_DFIcr_rowlevel.tsv"))
  write_if_rows(audit_summary_DFIcr,    file.path(out_dir, "audit_Parte_5_DFIcr_summary.tsv"))
  write_if_rows(audit_by_type_DFIcr,    file.path(out_dir, "audit_Parte_5_DFIcr_by_type.tsv"))
  write_if_rows(diag_case_mismatch_DFIcr, file.path(out_dir, "audit_Parte_5_DFIcr_case_mismatch.tsv"))
  readr::write_tsv(freq_DFIcr_event_type,  file.path(out_dir, "freq_Parte_5_DFIcr_event_type.tsv"))
  readr::write_tsv(freq_DFIcr_tumor_vital, file.path(out_dir, "freq_Parte_5_DFIcr_tumor_vital.tsv"))
  readr::write_tsv(freq_pivot_DFIcr,       file.path(out_dir, "freq_Parte_5_DFIcr_pivot.tsv"))
  
  cat("✓ Part 5 (DFI.cr) audit completed.\n")
} else {
  df009 <- df009_mid
}

## -----------------------------
## 4) Export final object
## -----------------------------
rio::export(df009, "data_table_results/df009.tsv")

## Mensagens finais
if (has_DFI) {
  cat("→ DFI Row-level:   ", normalizePath(file.path(out_dir, "audit_Parte_5_DFI_rowlevel.tsv")), "\n", sep = "")
  cat("→ DFI Summary:     ", normalizePath(file.path(out_dir, "audit_Parte_5_DFI_summary.tsv")),  "\n", sep = "")
  cat("→ DFI By-type:     ", normalizePath(file.path(out_dir, "audit_Parte_5_DFI_by_type.tsv")),  "\n", sep = "")
  cat("→ DFI Frequencies: ", normalizePath(file.path(out_dir, "freq_Parte_5_DFI_event_type.tsv")), "\n", sep = "")
  cat("                    ", normalizePath(file.path(out_dir, "freq_Parte_5_DFI_tumor_vital.tsv")), "\n", sep = "")
  cat("                    ", normalizePath(file.path(out_dir, "freq_Parte_5_DFI_pivot.tsv")),      "\n", sep = "")
  if (exists("time_check_DFI")) {
    cat("→ DFI Time-check:  ", normalizePath(file.path(out_dir, "audit_Parte_5_DFI_time_check.tsv")), "\n", sep = "")
  }
  if (exists("diag_case_mismatch_DFI") && nrow(diag_case_mismatch_DFI) > 0) {
    cat("→ DFI Case-mismatch:", normalizePath(file.path(out_dir, "audit_Parte_5_DFI_case_mismatch.tsv")), "\n", sep = "")
  }
}
if (has_DFIcr) {
  cat("→ DFI.cr Row-level:   ", normalizePath(file.path(out_dir, "audit_Parte_5_DFIcr_rowlevel.tsv")), "\n", sep = "")
  cat("→ DFI.cr Summary:     ", normalizePath(file.path(out_dir, "audit_Parte_5_DFIcr_summary.tsv")),  "\n", sep = "")
  cat("→ DFI.cr By-type:     ", normalizePath(file.path(out_dir, "audit_Parte_5_DFIcr_by_type.tsv")),  "\n", sep = "")
  cat("→ DFI.cr Frequencies: ", normalizePath(file.path(out_dir, "freq_Parte_5_DFIcr_event_type.tsv")), "\n", sep = "")
  cat("                      ", normalizePath(file.path(out_dir, "freq_Parte_5_DFIcr_tumor_vital.tsv")), "\n", sep = "")
  cat("                      ", normalizePath(file.path(out_dir, "freq_Parte_5_DFIcr_pivot.tsv")),      "\n", sep = "")
  if (exists("time_check_DFIcr")) {
    cat("→ DFI.cr Time-check:  ", normalizePath(file.path(out_dir, "audit_Parte_5_DFIcr_time_check.tsv")), "\n", sep = "")
  }
  if (exists("diag_case_mismatch_DFIcr") && nrow(diag_case_mismatch_DFIcr) > 0) {
    cat("→ DFI.cr Case-mismatch:", normalizePath(file.path(out_dir, "audit_Parte_5_DFIcr_case_mismatch.tsv")), "\n", sep = "")
  }
}


## ============================================================
## PART 6 — Harmonize new_tumor_event_dx_days_to ↔ PFI / PFI.time + AUDIT
##   Rules (evaluated on the "before" state):
##   R1) If PFI==1 & !is.na(PFI.time) & is.na(dx)           → dx := PFI.time
##   R2) If is.na(PFI) & is.na(PFI.time) & !is.na(dx)       → PFI := 1; PFI.time := dx
##   (where dx = new_tumor_event_dx_days_to)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

out_dir <- "harmonization_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## -----------------------------
## 0) Preconditions
## -----------------------------
stopifnot(exists("df009"))
req <- c("PFI","PFI.time","new_tumor_event_dx_days_to")
stopifnot(all(req %in% names(df009)))
has_type_column <- "type" %in% names(df009)

## -----------------------------
## 1) "Before" snapshot
## -----------------------------
df010_before <- df009 %>%
  mutate(
    dx0   = suppressWarnings(as.numeric(new_tumor_event_dx_days_to)),
    pfi0  = suppressWarnings(as.numeric(PFI)),
    pfit0 = suppressWarnings(as.numeric(PFI.time)),
    .rowid = dplyr::row_number()
  )

dx_before    <- df010_before$dx0
pfi_before   <- df010_before$pfi0
pfit_before  <- df010_before$pfit0

## -----------------------------
## 2) Flags (evaluated on "before")
## -----------------------------
# R1: pull dx from PFI.time when PFI==1 and dx is NA
flag_R1_dx_from_PFI <- (!is.na(pfi_before) & pfi_before == 1) &
  (!is.na(pfit_before)) &
  (is.na(dx_before))

# R2: infer PFI and PFI.time from dx when both were NA
flag_R2_fill_PFI_from_dx <- (is.na(pfi_before) & is.na(pfit_before) & !is.na(dx_before))

## -----------------------------
## 3) EXPECTED — apply rules R1 and R2 on the snapshot
## -----------------------------
dx_expected   <- ifelse(flag_R1_dx_from_PFI, pfit_before, dx_before)
pfi_expected  <- ifelse(flag_R2_fill_PFI_from_dx, 1, pfi_before)
pfit_expected <- ifelse(flag_R2_fill_PFI_from_dx, dx_before, pfit_before)

## -----------------------------
## 4) OBSERVED — apply harmonization using "before" values
## -----------------------------
df010_after <- df010_before %>%
  mutate(
    # R1: dx from PFI.time (only when eligible)
    new_dx   = if_else(flag_R1_dx_from_PFI, pfit0, dx0),
    # R2: PFI and PFI.time from dx (only when eligible)
    new_PFI  = if_else(flag_R2_fill_PFI_from_dx, 1, pfi0),
    new_PFIT = if_else(flag_R2_fill_PFI_from_dx, dx0, pfit0)
  ) %>%
  # write into final columns
  mutate(
    new_tumor_event_dx_days_to = new_dx,
    PFI      = new_PFI,
    PFI.time = new_PFIT
  )

dx_after   <- df010_after$new_tumor_event_dx_days_to
pfi_after  <- df010_after$PFI
pfit_after <- df010_after$PFI.time

## -----------------------------
## 5) Audit — formal checks
## -----------------------------

# (a) Observed == Expected (row by row)
mismatch_dx   <- which( !(is.na(dx_after)   & is.na(dx_expected))   & (dx_after   != dx_expected) )
mismatch_pfi  <- which( !(is.na(pfi_after)  & is.na(pfi_expected))  & (pfi_after  != pfi_expected) )
mismatch_pfit <- which( !(is.na(pfit_after) & is.na(pfit_expected)) & (pfit_after != pfit_expected) )
if (length(mismatch_dx) + length(mismatch_pfi) + length(mismatch_pfit) > 0) {
  bad_rows <- unique(c(mismatch_dx, mismatch_pfi, mismatch_pfit))
  stop("Audit failure: differences between observed and expected in ",
       length(bad_rows), " row(s). Examples (rowid): ",
       paste(head(df010_before$.rowid[bad_rows], 10), collapse = ", "))
}

# (b) Immutability outside the rules
changed_dx_outside  <- which((!flag_R1_dx_from_PFI) &
                               ((is.na(dx_before) != is.na(dx_after)) |
                                  (!is.na(dx_before) & !is.na(dx_after) & dx_before != dx_after)))
changed_pfi_outside <- which((!flag_R2_fill_PFI_from_dx) &
                               ((is.na(pfi_before) != is.na(pfi_after)) |
                                  (!is.na(pfi_before) & !is.na(pfi_after) & pfi_before != pfi_after)))
changed_pfit_outside <- which((!flag_R2_fill_PFI_from_dx) &
                                ((is.na(pfit_before) != is.na(pfit_after)) |
                                   (!is.na(pfit_before) & !is.na(pfit_after) & pfit_before != pfit_after)))
if (length(changed_dx_outside) > 0 || length(changed_pfi_outside) > 0 || length(changed_pfit_outside) > 0) {
  stop("Audit failure: changes detected outside R1/R2 conditions.")
}

# (c) Domains and sanity (warnings)
invalid_pfi_vals <- setdiff(unique(pfi_after[!is.na(pfi_after)]), c(0,1))
if (length(invalid_pfi_vals) > 0) {
  warning("Warning: PFI values outside domain {0,1,NA}: ",
          paste(invalid_pfi_vals, collapse = ", "),
          ". Harmonization does not correct pre-existing values.")
}
neg_times_idx <- which(!is.na(dx_after) & dx_after < 0 | !is.na(pfit_after) & pfit_after < 0)
if (length(neg_times_idx) > 0) {
  warning("Warning: negative times detected in ", length(neg_times_idx), " row(s). Check data source.")
}

## -----------------------------
## 6) Reports
## -----------------------------
changed_idx <- which(
  (is.na(dx_before)   != is.na(dx_after))   | (!is.na(dx_before)   & !is.na(dx_after)   & dx_before   != dx_after)   |
    (is.na(pfi_before)  != is.na(pfi_after))  | (!is.na(pfi_before)  & !is.na(pfi_after)  & pfi_before  != pfi_after)  |
    (is.na(pfit_before) != is.na(pfit_after)) | (!is.na(pfit_before) & !is.na(pfit_after) & pfit_before != pfit_after)
)

audit_rowlevel <- tibble(
  rowid  = df010_before$.rowid[changed_idx],
  type   = if (has_type_column) df010_before$type[changed_idx] else NA_character_,
  dx_before   = dx_before[changed_idx],
  dx_after    = dx_after[changed_idx],
  pfi_before  = pfi_before[changed_idx],
  pfi_after   = pfi_after[changed_idx],
  pfit_before = pfit_before[changed_idx],
  pfit_after  = pfit_after[changed_idx],
  R1_dx_from_PFI   = flag_R1_dx_from_PFI[changed_idx],
  R2_fill_from_dx  = flag_R2_fill_PFI_from_dx[changed_idx]
)

audit_summary <- tibble(
  total_rows            = nrow(df010_before),
  n_R1_dx_from_PFI      = sum(flag_R1_dx_from_PFI),
  n_R2_fill_from_dx     = sum(flag_R2_fill_PFI_from_dx),
  n_dx_changed          = sum((is.na(dx_before) != is.na(dx_after)) |
                                (!is.na(dx_before) & !is.na(dx_after) & dx_before != dx_after)),
  n_pfi_changed         = sum((is.na(pfi_before) != is.na(pfi_after)) |
                                (!is.na(pfi_before) & !is.na(pfi_after) & pfi_before != pfi_after)),
  n_pfit_changed        = sum((is.na(pfit_before) != is.na(pfit_after)) |
                                (!is.na(pfit_before) & !is.na(pfit_after) & pfit_before != pfit_after)),
  check_R1_impacts_dx   = n_R1_dx_from_PFI == n_dx_changed - n_R2_fill_from_dx,  # dx changes by R1; R2 may also affect dx? (no)
  check_R2_impacts_pfi  = n_R2_fill_from_dx == n_pfi_changed,
  check_R2_impacts_pfit = n_R2_fill_from_dx == n_pfit_changed
)

if (has_type_column) {
  audit_by_type <- tibble(
    type = df010_before$type,
    R1 = flag_R1_dx_from_PFI,
    R2 = flag_R2_fill_PFI_from_dx
  ) %>%
    group_by(type) %>%
    summarise(n_R1 = sum(R1), n_R2 = sum(R2), .groups = "drop") %>%
    arrange(desc(n_R1 + n_R2), type)
} else {
  audit_by_type <- tibble(info = "Column 'type' does not exist; report by type not generated.")
}

# Export reports
readr::write_tsv(audit_rowlevel, file.path(out_dir, "audit_Parte_6_dxPFI_rowlevel.tsv"))
readr::write_tsv(audit_summary,  file.path(out_dir, "audit_Parte_6_dxPFI_summary.tsv"))
readr::write_tsv(audit_by_type,  file.path(out_dir, "audit_Parte_6_dxPFI_by_type.tsv"))

## -----------------------------
## 7) Update working object and log
## -----------------------------
df010 <- df010_after %>%
  select(-dx0, -pfi0, -pfit0, -new_dx, -new_PFI, -new_PFIT, -.rowid)

rio::export(df010,"data_table_results/df010.tsv")

cat("✓ Part 6 harmonization + audit completed.\n")
cat("→ Row-level: ", normalizePath(file.path(out_dir, "audit_Parte_6_dxPFI_rowlevel.tsv")), "\n", sep = "")
cat("→ Summary:   ", normalizePath(file.path(out_dir, "audit_Parte_6_dxPFI_summary.tsv")),  "\n", sep = "")
cat("→ By-type:   ", normalizePath(file.path(out_dir, "audit_Parte_6_dxPFI_by_type.tsv")),  "\n", sep = "")

## ============================================================
## PART 7 — "Meaningful NA" + AUDIT
##   R1) cause_of_death := "No event" if (OS==0 & is.na(cause_of_death))
##   R2) new_tumor_event_type := "No event"      if (PFI==0 & is.na(new_tumor_event_type))
##   R3) new_tumor_event_site := "No event"      if (PFI==0 & is.na(new_tumor_event_site))
##   R4) new_tumor_event_site_other := "No event" if (PFI==0 & is.na(new_tumor_event_site_other))
##   R5) new_tumor_event_dx_days_to := PFI.time   if (is.na(dx) & PFI==0)
##   R6) new_tumor_event_dx_days_to := NA         if (dx < 0)   [after R5]
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

out_dir <- "harmonization_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## -----------------------------
## 0) Preconditions
## -----------------------------
stopifnot(exists("df010"))
req <- c("OS","PFI","PFI.time","cause_of_death",
         "new_tumor_event_type","new_tumor_event_site","new_tumor_event_site_other",
         "new_tumor_event_dx_days_to")
stopifnot(all(req %in% names(df010)))
has_type_column <- "type" %in% names(df010)

## -----------------------------
## 1) "Before" snapshot (explicit types)
## -----------------------------
df011_before <- df010 %>%
  mutate(
    OS   = suppressWarnings(as.numeric(OS)),
    PFI  = suppressWarnings(as.numeric(PFI)),
    PFI.time = suppressWarnings(as.numeric(PFI.time)),
    cause_of_death            = as.character(cause_of_death),
    new_tumor_event_type      = as.character(new_tumor_event_type),
    new_tumor_event_site      = as.character(new_tumor_event_site),
    new_tumor_event_site_other= as.character(new_tumor_event_site_other),
    dx0  = suppressWarnings(as.numeric(new_tumor_event_dx_days_to)),
    .rowid = dplyr::row_number()
  )

# "Before" vectors
OS0   <- df011_before$OS
PFI0  <- df011_before$PFI
PFIT0 <- df011_before$PFI.time
C0    <- df011_before$cause_of_death
T0    <- df011_before$new_tumor_event_type
S0    <- df011_before$new_tumor_event_site
SO0   <- df011_before$new_tumor_event_site_other
DX0   <- df011_before$dx0

## -----------------------------
## 2) Flags (evaluated on "before")
## -----------------------------
f_R1_cause <- (OS0 == 0) & is.na(C0)
f_R2_type  <- (PFI0 == 0) & is.na(T0)
f_R3_site  <- (PFI0 == 0) & is.na(S0)
f_R4_siteo <- (PFI0 == 0) & is.na(SO0)
f_R5_dx    <- is.na(DX0) & (PFI0 == 0)              # if PFIT0 is NA, result remains NA (consistent with your code)

## -----------------------------
## 3) EXPECTED — apply R1..R6 on the snapshot
## -----------------------------
C_exp  <- ifelse(f_R1_cause, "No event", C0)
T_exp  <- ifelse(f_R2_type,  "No event", T0)
S_exp  <- ifelse(f_R3_site,  "No event", S0)
SO_exp <- ifelse(f_R4_siteo, "No event", SO0)

DX1_exp <- ifelse(f_R5_dx, PFIT0, DX0)              # R5
DX2_exp <- ifelse(DX1_exp < 0, NA_real_, DX1_exp)   # R6 (after R5)

## -----------------------------
## 4) OBSERVED — apply your BLOCK for comparison
## -----------------------------
df011_after <- df010 %>%
  mutate(
    cause_of_death = ifelse(OS == 0 & is.na(cause_of_death), "No event", cause_of_death),
    new_tumor_event_type = ifelse(PFI == 0 & is.na(new_tumor_event_type), "No event", new_tumor_event_type),
    new_tumor_event_site = ifelse(PFI == 0 & is.na(new_tumor_event_site), "No event", new_tumor_event_site),
    new_tumor_event_site_other = ifelse(PFI == 0 & is.na(new_tumor_event_site_other), "No event", new_tumor_event_site_other)
  ) %>%
  mutate(
    new_tumor_event_dx_days_to = ifelse(
      is.na(new_tumor_event_dx_days_to) & PFI == 0,
      PFI.time,
      new_tumor_event_dx_days_to
    ),
    new_tumor_event_dx_days_to = ifelse(
      new_tumor_event_dx_days_to < 0,
      NA_real_,
      new_tumor_event_dx_days_to
    )
  ) %>%
  mutate(
    # Types for comparison
    cause_of_death            = as.character(cause_of_death),
    new_tumor_event_type      = as.character(new_tumor_event_type),
    new_tumor_event_site      = as.character(new_tumor_event_site),
    new_tumor_event_site_other= as.character(new_tumor_event_site_other),
    dx_obs = suppressWarnings(as.numeric(new_tumor_event_dx_days_to))
  )

C_obs  <- df011_after$cause_of_death
T_obs  <- df011_after$new_tumor_event_type
S_obs  <- df011_after$new_tumor_event_site
SO_obs <- df011_after$new_tumor_event_site_other
DX_obs <- df011_after$dx_obs

## -----------------------------
## 5) Audit — formal checks
## -----------------------------

# (a) Observed == Expected (all target columns)
m_mismatch <- which( (C_obs  != C_exp  | xor(is.na(C_obs),  is.na(C_exp)))  |
                       (T_obs  != T_exp  | xor(is.na(T_obs),  is.na(T_exp)))  |
                       (S_obs  != S_exp  | xor(is.na(S_obs),  is.na(S_exp)))  |
                       (SO_obs != SO_exp | xor(is.na(SO_obs), is.na(SO_exp))) |
                       (DX_obs != DX2_exp| xor(is.na(DX_obs), is.na(DX2_exp))) )
if (length(m_mismatch) > 0) {
  stop("Audit failure: ", length(m_mismatch),
       " row(s) differ between observed and expected. Examples (rowid): ",
       paste(head(df011_before$.rowid[m_mismatch], 10), collapse = ", "))
}

# (b) Do not alter non-NA values in 'event*' and 'cause_of_death'
changed_non_na <- function(bef, aft) {
  which(!is.na(bef) & !is.na(aft) & bef != aft)
}
bad_cause <- changed_non_na(C0,  C_obs)
bad_type  <- changed_non_na(T0,  T_obs)
bad_site  <- changed_non_na(S0,  S_obs)
bad_siteo <- changed_non_na(SO0, SO_obs)
if (length(bad_cause)+length(bad_type)+length(bad_site)+length(bad_siteo) > 0) {
  stop("Audit failure: non-NA values were changed in event/cause fields.")
}

# (c) Only the allowed label was introduced
introduced_vals <- function(bef, aft) {
  unique(aft[ is.na(bef) & !is.na(aft) ])
}
intro_cause <- introduced_vals(C0,  C_obs)
intro_type  <- introduced_vals(T0,  T_obs)
intro_site  <- introduced_vals(S0,  S_obs)
intro_siteo <- introduced_vals(SO0, SO_obs)
allowed <- "No event"
if (any(setdiff(na.omit(intro_cause), allowed) %>% length() > 0) ||
    any(setdiff(na.omit(intro_type),  allowed) %>% length() > 0) ||
    any(setdiff(na.omit(intro_site),  allowed) %>% length() > 0) ||
    any(setdiff(na.omit(intro_siteo), allowed) %>% length() > 0)) {
  stop("Audit failure: labels different from 'No event' were introduced.")
}

# (d) Domain of DX after treatment (non-negative or NA)
left_neg <- which(!is.na(DX_obs) & DX_obs < 0)
if (length(left_neg) > 0) {
  stop("Audit failure: negative values remain in 'new_tumor_event_dx_days_to' after harmonization.")
}

## -----------------------------
## 6) Reports
## -----------------------------
changed_idx <- which(
  (C0  != C_obs  | xor(is.na(C0),  is.na(C_obs)))  |
    (T0  != T_obs  | xor(is.na(T0),  is.na(T_obs)))  |
    (S0  != S_obs  | xor(is.na(S0),  is.na(S_obs)))  |
    (SO0 != SO_obs | xor(is.na(SO0), is.na(SO_obs))) |
    (DX0 != DX_obs | xor(is.na(DX0), is.na(DX_obs)))
)

audit_rowlevel <- tibble(
  rowid  = df011_before$.rowid[changed_idx],
  type   = if (has_type_column) df011_before$type[changed_idx] else NA_character_,
  # Before/After
  cause_before = C0[changed_idx],  cause_after = C_obs[changed_idx],
  type_before  = T0[changed_idx],  type_after  = T_obs[changed_idx],
  site_before  = S0[changed_idx],  site_after  = S_obs[changed_idx],
  siteo_before = SO0[changed_idx], siteo_after = SO_obs[changed_idx],
  dx_before = DX0[changed_idx],    dx_after    = DX_obs[changed_idx],
  # Flags
  R1_cause = f_R1_cause[changed_idx],
  R2_type  = f_R2_type[changed_idx],
  R3_site  = f_R3_site[changed_idx],
  R4_siteo = f_R4_siteo[changed_idx],
  R5_dx    = f_R5_dx[changed_idx],
  R6_neg2NA= (!is.na(DX1_exp) & DX1_exp < 0)[changed_idx]
)

audit_summary <- tibble(
  total_rows              = nrow(df011_before),
  n_R1_cause              = sum(f_R1_cause),
  n_R2_type               = sum(f_R2_type),
  n_R3_site               = sum(f_R3_site),
  n_R4_siteo              = sum(f_R4_siteo),
  n_R5_dx_from_PFItime    = sum(f_R5_dx),
  n_R6_neg2NA_applied     = sum(!is.na(DX1_exp) & DX1_exp < 0),
  # how many "No event" were actually introduced
  n_intro_cause           = sum(is.na(C0)  & !is.na(C_obs)  & C_obs  == "No event"),
  n_intro_type            = sum(is.na(T0)  & !is.na(T_obs)  & T_obs  == "No event"),
  n_intro_site            = sum(is.na(S0)  & !is.na(S_obs)  & S_obs  == "No event"),
  n_intro_siteo           = sum(is.na(SO0) & !is.na(SO_obs) & SO_obs == "No event"),
  # dx changed
  n_dx_changed            = sum((is.na(DX0) != is.na(DX_obs)) | (!is.na(DX0) & !is.na(DX_obs) & DX0 != DX_obs)),
  # sanity checks balance
  check_cause_intro       = n_intro_cause == n_R1_cause,
  check_type_intro        = n_intro_type  == n_R2_type,
  check_site_intro        = n_intro_site  == n_R3_site,
  check_siteo_intro       = n_intro_siteo == n_R4_siteo
)

if (has_type_column) {
  audit_by_type <- audit_rowlevel %>%
    mutate(rule = case_when(
      R1_cause ~ "R1_cause:=No event",
      R2_type  ~ "R2_type:=No event",
      R3_site  ~ "R3_site:=No event",
      R4_siteo ~ "R4_siteo:=No event",
      R5_dx    ~ "R5_dx:=PFI.time",
      R6_neg2NA~ "R6_neg2NA",
      TRUE     ~ "Other"
    )) %>%
    tidyr::separate_rows(rule, sep = ";") %>%
    count(type, rule, name = "n_applied") %>%
    arrange(type, desc(n_applied))
} else {
  audit_by_type <- tibble(info = "Column 'type' does not exist; report by type not generated.")
}

## -----------------------------
## 7) Exports and update working object
## -----------------------------
readr::write_tsv(audit_rowlevel, file.path(out_dir, "audit_Parte_7_NoEvent_rowlevel.tsv"))
readr::write_tsv(audit_summary,  file.path(out_dir, "audit_Parte_7_NoEvent_summary.tsv"))
readr::write_tsv(audit_by_type,  file.path(out_dir, "audit_Parte_7_NoEvent_by_type.tsv"))

# Update working object
df011 <- df011_after %>%
  dplyr::select(-dplyr::any_of(c("dx_obs", ".rowid")))

rio::export(df011,"data_table_results/df011.tsv")

cat("✓ Part 7 audit completed.\n")
cat("→ Row-level: ", normalizePath(file.path(out_dir, "audit_Parte_7_NoEvent_rowlevel.tsv")), "\n", sep = "")
cat("→ Summary:   ", normalizePath(file.path(out_dir, "audit_Parte_7_NoEvent_summary.tsv")),  "\n", sep = "")
cat("→ By-type:   ", normalizePath(file.path(out_dir, "audit_Parte_7_NoEvent_by_type.tsv")),  "\n", sep = "")

## ============================================================
## PART 8 (v3) — Harmonization of cause_of_death + AUDIT (multi-endpoint aware)
## Regras resumidas:
##   Classe derivada (prioriza DSS_cr; senão DSS; exige OS==1):
##     death_cancer := OS==1 & DSS_cr==1  (ou DSS==1)
##     death_other  := OS==1 & DSS_cr==2  (ou DSS==0)
##   Cond0: cause missing & sem info de óbito → keep
##   Cond1: cause missing & death_other                       → "Other, non-malignant disease"
##   Cond2: cause missing & death_cancer & event type missing → use type
##   Cond3: cause missing & death_cancer & Metastasis & site missing → use type
##   Cond4: cause missing & death_cancer & Metastasis & site present → randomize {type, site}
##   Cond5: cause missing & death_cancer & não-metastático    → use type
##   Pós:  mapear "type" → "tumor_name" (se cause == type)
## Auditoria:
##   - Row-level: antes/esperado/depois + Cond1–Cond5
##   - Summary: totais, #mudanças, #por regra, #eq esperado=observado
##   - By-type: distribuição por tumor e regra
##   - Diagnósticos extras (cross-tabs entre faltantes)
## ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(purrr); library(tibble); library(stringr)
})

## -----------------------------
## 0) Parâmetros e pré-condições
## -----------------------------
out_dir <- "harmonization_results"
out_res <- "data_table_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_res, recursive = TRUE, showWarnings = FALSE)

stopifnot(exists("df011"))
stopifnot(exists("df005"))   # usaremos a ordem original do df005

req  <- c("type", "OS", "cause_of_death",
          "new_tumor_event_type", "new_tumor_event_site")
opt  <- c("DSS", "DSS_cr")
stopifnot(all(req %in% names(df011)))

if (!any(opt %in% names(df011))) {
  warning("PART 8: Nenhuma coluna de DSS/DSS_cr encontrada; apenas Cond0 (keep) será possível quando cause_of_death estiver ausente.")
}

## Helper
is_true <- function(x) !is.na(x) & x

## -----------------------------
## 1) Snapshot 'before' (tipagem + trims)
## -----------------------------
df012_before <- df011 %>%
  mutate(
    OS     = suppressWarnings(as.numeric(OS)),
    DSS    = if ("DSS"    %in% names(.)) suppressWarnings(as.numeric(DSS))    else NA_real_,
    DSS_cr = if ("DSS_cr" %in% names(.)) suppressWarnings(as.numeric(DSS_cr)) else NA_real_,
    cause_of_death       = as.character(cause_of_death),
    new_tumor_event_type = as.character(new_tumor_event_type),
    new_tumor_event_site = as.character(new_tumor_event_site),
    cause0  = cause_of_death,
    cause0t = ifelse(is.na(cause_of_death), NA_character_, trimws(cause_of_death)),
    type0t  = ifelse(is.na(new_tumor_event_type), NA_character_, trimws(new_tumor_event_type)),
    site0t  = ifelse(is.na(new_tumor_event_site), NA_character_, trimws(new_tumor_event_site)),
    .rowid  = dplyr::row_number()
  )

cause0 <- df012_before$cause0

## -----------------------------
## 2) Diagnóstico de conflitos DSS vs DSS_cr
## -----------------------------
diag_conflict <- tibble()
if (all(c("DSS","DSS_cr") %in% names(df012_before))) {
  conflict_idx <- which(
    !is.na(df012_before$DSS) & !is.na(df012_before$DSS_cr) &
      ( (df012_before$DSS == 1 & df012_before$DSS_cr != 1) |
          (df012_before$DSS == 0 & df012_before$DSS_cr == 1) )
  )
  diag_conflict <- df012_before %>%
    filter(.rowid %in% conflict_idx) %>%
    select(.rowid, type, OS, DSS, DSS_cr,
           new_tumor_event_type, new_tumor_event_site, cause_of_death)
  if (nrow(diag_conflict) > 0) {
    readr::write_tsv(diag_conflict, file.path(out_dir, "audit_Parte_8_conflict_DSS_vs_DSScr.tsv"))
  }
}

## -----------------------------
## 3) Mapeamento tumor (Study Abbrev → Nome)
## -----------------------------
tumor_map <- tribble(
  ~type, ~tumor_name,
  "LAML","Acute Myeloid Leukemia",               "ACC","Adrenocortical carcinoma",
  "BLCA","Bladder Urothelial Carcinoma",         "LGG","Brain Lower Grade Glioma",
  "BRCA","Breast invasive carcinoma",            "CESC","Cervical squamous cell carcinoma and endocervical adenocarcinoma",
  "CHOL","Cholangiocarcinoma",                   "LCML","Chronic Myelogenous Leukemia",
  "COAD","Colon adenocarcinoma",                 "CNTL","Controls",
  "ESCA","Esophageal carcinoma",                 "FPPP","FFPE Pilot Phase II",
  "GBM","Glioblastoma multiforme",               "HNSC","Head and Neck squamous cell carcinoma",
  "KICH","Kidney Chromophobe",                   "KIRC","Kidney renal clear cell carcinoma",
  "KIRP","Kidney renal papillary cell carcinoma","LIHC","Liver hepatocellular carcinoma",
  "LUAD","Lung adenocarcinoma",                  "LUSC","Lung squamous cell carcinoma",
  "DLBC","Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
  "MESO","Mesothelioma",                         "MISC","Miscellaneous",
  "OV","Ovarian serous cystadenocarcinoma",      "PAAD","Pancreatic adenocarcinoma",
  "PCPG","Pheochromocytoma and Paraganglioma",   "PRAD","Prostate adenocarcinoma",
  "READ","Rectum adenocarcinoma",                "SARC","Sarcoma",
  "SKCM","Skin Cutaneous Melanoma",              "STAD","Stomach adenocarcinoma",
  "TGCT","Testicular Germ Cell Tumors",          "THYM","Thymoma",
  "THCA","Thyroid carcinoma",                    "UCS","Uterine Carcinosarcoma",
  "UCEC","Uterine Corpus Endometrial Carcinoma", "UVM","Uveal Melanoma"
)

## -----------------------------
## 4) Classe derivada de óbito e flags
## -----------------------------
metastasis_re <- regex("distant\\s+metastasis|metastatic", ignore_case = TRUE)

death_cancer <- with(df012_before, ifelse(
  !is.na(DSS_cr), OS == 1 & DSS_cr == 1,
  ifelse(!is.na(DSS),   OS == 1 & DSS   == 1, NA)
))
death_other <- with(df012_before, ifelse(
  !is.na(DSS_cr), OS == 1 & DSS_cr == 2,
  ifelse(!is.na(DSS),   OS == 1 & DSS   == 0, NA)
))
has_death_info <- !is.na(death_cancer) | !is.na(death_other)

## -----------------------------
## 5) Aplicação das regras (EXPECTED e OBSERVED)
## -----------------------------
set.seed(123)

apply_rules_by_type <- function(sub_df) {
  n <- nrow(sub_df)
  r <- runif(n)  # usado apenas em Cond4
  
  miss_cause <- is.na(sub_df$cause0t) | sub_df$cause0t == ""
  meta_ev    <- str_detect(sub_df$type0t, metastasis_re)
  site_abs   <- is.na(sub_df$site0t) | sub_df$site0t == ""
  site_pres  <- !site_abs
  type_miss  <- is.na(sub_df$type0t) | sub_df$type0t == ""
  
  idx  <- match(sub_df$.rowid, df012_before$.rowid)
  dc   <- death_cancer[idx]
  do   <- death_other[idx]
  info <- has_death_info[idx]
  
  expected <- sub_df %>%
    mutate(
      cause_expected = case_when(
        miss_cause & !info                             ~ cause0,                              # Cond0
        miss_cause & is_true(do)                       ~ "Other, non-malignant disease",      # Cond1
        miss_cause & is_true(dc) & type_miss           ~ type,                                # Cond2
        miss_cause & is_true(dc) & meta_ev & site_abs  ~ type,                                # Cond3
        miss_cause & is_true(dc) & meta_ev & site_pres ~ ifelse(r > 0.5, type, sub_df$site0t),# Cond4
        miss_cause & is_true(dc) & !meta_ev            ~ type,                                # Cond5
        TRUE                                           ~ cause0
      ),
      cond1 = miss_cause & is_true(do),
      cond2 = miss_cause & is_true(dc) & type_miss,
      cond3 = miss_cause & is_true(dc) & meta_ev & site_abs,
      cond4 = miss_cause & is_true(dc) & meta_ev & site_pres,
      cond5 = miss_cause & is_true(dc) & !meta_ev
    ) %>%
    select(.rowid, type, cause_expected, cond1:cond5)
  
  observed <- sub_df %>%
    mutate(
      cause_of_death = case_when(
        miss_cause & !info                             ~ cause0,                              # Cond0
        miss_cause & is_true(do)                       ~ "Other, non-malignant disease",      # Cond1
        miss_cause & is_true(dc) & type_miss           ~ type,                                # Cond2
        miss_cause & is_true(dc) & meta_ev & site_abs  ~ type,                                # Cond3
        miss_cause & is_true(dc) & meta_ev & site_pres ~ ifelse(r > 0.5, type, sub_df$site0t),# Cond4
        miss_cause & is_true(dc) & !meta_ev            ~ type,                                # Cond5
        TRUE                                           ~ cause0
      )
    ) %>%
    select(.rowid, type, cause_of_death)
  
  list(expected = expected, observed = observed)
}

split_list     <- df012_before %>% split(.$type) %>% map(apply_rules_by_type)
expected_stage1 <- map_dfr(split_list, "expected")
observed_stage1 <- map_dfr(split_list, "observed")

## Pós-processamento: mapear type → tumor_name (para expected; para observed faremos in-place abaixo)
expected_final <- expected_stage1 %>%
  left_join(tumor_map, by = "type") %>%
  mutate(cause_expected_final = ifelse(cause_expected == type, tumor_name, cause_expected)) %>%
  select(.rowid, type, cause_expected, cause_expected_final, cond1:cond5)

## ============================================================
## >>> ATUALIZAÇÃO IN-PLACE (NÃO ALTERA ORDEM DE COLUNAS) <<<
## ============================================================
df012_after <- df012_before
idx_obs <- match(df012_after$.rowid, observed_stage1$.rowid)
df012_after$cause_of_death <- observed_stage1$cause_of_death[idx_obs]

## Substituir rótulo igual a 'type' por 'tumor_name' sem mover colunas:
df012_after <- df012_after %>%
  left_join(tumor_map, by = "type") %>%
  mutate(cause_of_death = ifelse(cause_of_death == type, tumor_name, cause_of_death)) %>%
  select(-tumor_name)

## -----------------------------
## 6) AUDIT — comparações e relatórios
## -----------------------------
obs_tbl <- df012_after %>%
  transmute(.rowid, type, cause_obs = trimws(cause_of_death))

cmp <- expected_final %>%
  mutate(cause_expected_final = trimws(cause_expected_final)) %>%
  left_join(obs_tbl, by = c(".rowid","type")) %>%
  mutate(eq = ( (is.na(cause_obs) & is.na(cause_expected_final)) |
                  (!is.na(cause_obs) & !is.na(cause_expected_final) &
                     cause_obs == cause_expected_final) ))

before_trim <- trimws(df012_before$cause0)[match(cmp$.rowid, df012_before$.rowid)]
after_trim  <- cmp$cause_obs
changed_idx <- which( (is.na(before_trim) != is.na(after_trim)) |
                        (!is.na(before_trim) & !is.na(after_trim) & before_trim != after_trim) )

audit_rowlevel <- cmp %>%
  transmute(
    rowid = .rowid, type,
    cause_before   = before_trim,
    cause_expected = cause_expected_final,
    cause_after    = cause_obs,
    cond1, cond2, cond3, cond4, cond5,
    expected_equals_observed = eq
  ) %>%
  filter(rowid %in% cmp$.rowid[changed_idx])

audit_summary <- tibble(
  total_rows                 = nrow(df012_before),
  n_changed                  = nrow(audit_rowlevel),
  n_expected_equals_observed = sum(cmp$eq, na.rm = TRUE),
  n_expected_not_equal       = sum(!cmp$eq, na.rm = TRUE),
  n_cond1_applied            = sum(cmp$cond1, na.rm = TRUE),
  n_cond2_applied            = sum(cmp$cond2, na.rm = TRUE),
  n_cond3_applied            = sum(cmp$cond3, na.rm = TRUE),
  n_cond4_applied            = sum(cmp$cond4, na.rm = TRUE),
  n_cond5_applied            = sum(cmp$cond5, na.rm = TRUE)
)

audit_by_type <- audit_rowlevel %>%
  mutate(rule = case_when(
    cond1 ~ "Cond1_non_malignant",
    cond2 ~ "Cond2_type_missing_event",
    cond3 ~ "Cond3_metastasis_site_abs",
    cond4 ~ "Cond4_metastasis_site_pres_rand",
    cond5 ~ "Cond5_type_non_metastatic",
    TRUE  ~ "Other_or_NoRule"
  )) %>%
  count(type, rule, name = "n_applied") %>%
  arrange(type, desc(n_applied))

readr::write_tsv(audit_rowlevel, file.path(out_dir, "audit_Parte_8v3_cause_rowlevel.tsv"))
readr::write_tsv(audit_summary,  file.path(out_dir, "audit_Parte_8v3_cause_summary.tsv"))
readr::write_tsv(audit_by_type,  file.path(out_dir, "audit_Parte_8v3_cause_by_type.tsv"))

## -----------------------------
## 7) Diagnósticos extras (não alteram regras)
## -----------------------------
diag_tbl <- df012_before %>%
  mutate(
    cause_miss = is.na(cause_of_death) | cause_of_death == "",
    DSS_miss   = is.na(DSS),
    OS1        = !is.na(OS)  & OS  == 1,
    DSS0       = !is.na(DSS) & DSS == 0,
    DSS1       = !is.na(DSS) & DSS == 1,
    metastasis = str_detect(new_tumor_event_type,
                            regex("Distant Metastasis|Metastatic", ignore_case = TRUE)),
    site_abs   = is.na(new_tumor_event_site) | new_tumor_event_site == "",
    site_pres  = !site_abs
  )

n_cond0_keep <- sum(diag_tbl$cause_miss & diag_tbl$DSS_miss, na.rm = TRUE)

xtab_OS_DSS <- diag_tbl %>%
  filter(cause_miss) %>%
  mutate(
    OS_cat  = case_when(is.na(OS)  ~ "OS=NA",  OS==0 ~ "OS=0",  OS==1 ~ "OS=1",  TRUE ~ "OS=?"),
    DSS_cat = case_when(is.na(DSS) ~ "DSS=NA", DSS==0 ~ "DSS=0", DSS==1 ~ "DSS=1", TRUE ~ "DSS=?")
  ) %>%
  count(OS_cat, DSS_cat, name = "n")

meta_diag <- diag_tbl %>%
  filter(cause_miss) %>%
  summarise(
    n_miss_cause            = n(),
    n_OS1                   = sum(OS1,  na.rm = TRUE),
    n_DSS0                  = sum(DSS0, na.rm = TRUE),
    n_DSS1                  = sum(DSS1, na.rm = TRUE),
    n_metastasis            = sum(metastasis, na.rm = TRUE),
    n_meta_site_abs         = sum(metastasis & site_abs,  na.rm = TRUE),
    n_meta_site_pres        = sum(metastasis & site_pres, na.rm = TRUE),
    .groups = "drop"
  )

n_if_DSScr <- NA_integer_
if ("DSS_cr" %in% names(df012_before)) {
  n_if_DSScr <- df012_before %>%
    mutate(cause_miss = is.na(cause_of_death) | cause_of_death == "") %>%
    summarise(n = sum(cause_miss & !is.na(DSS_cr) & DSS_cr == 1, na.rm = TRUE)) %>%
    pull(n)
}

audit_summary_extra <- audit_summary %>%
  mutate(
    n_cond0_keep_DSS_missing       = n_cond0_keep,
    n_cause_missing_total          = sum(diag_tbl$cause_miss, na.rm = TRUE),
    n_cause_missing_OS1            = sum(diag_tbl$cause_miss & diag_tbl$OS1, na.rm = TRUE),
    n_cause_missing_DSS0           = sum(diag_tbl$cause_miss & diag_tbl$DSS0, na.rm = TRUE),
    n_cause_missing_DSS1           = sum(diag_tbl$cause_miss & diag_tbl$DSS1, na.rm = TRUE),
    n_cause_missing_metastasis     = meta_diag$n_metastasis,
    n_cause_missing_metastasis_site_abs  = meta_diag$n_meta_site_abs,
    n_cause_missing_metastasis_site_pres = meta_diag$n_meta_site_pres,
    n_diag_if_use_DSScr_as_fallback     = n_if_DSScr
  )

readr::write_tsv(xtab_OS_DSS,         file.path(out_dir, "audit_Parte_8_diag_xtab_OSxDSS_among_missing_cause.tsv"))
readr::write_tsv(meta_diag,           file.path(out_dir, "audit_Parte_8_diag_missing_cause_metastasis.tsv"))
readr::write_tsv(audit_summary_extra, file.path(out_dir, "audit_Parte_8v3_cause_summary_with_diagnostics.tsv"))

## -----------------------------
## 8) Remover auxiliares e PRESERVAR ORDEM DO df005
## -----------------------------
aux_cols <- c("cause0","cause0t","type0t","site0t",".rowid")
df012 <- df012_after %>% dplyr::select(-dplyr::any_of(aux_cols))

## Alinhar a ordem final à do df005, mantendo colunas extras de df012 ao final
baseline_order     <- names(df005)
present_in_df005   <- intersect(baseline_order, names(df012))
extras_df012       <- setdiff(names(df012), baseline_order)
df012 <- df012[, c(present_in_df005, extras_df012)]

## Auditoria da ordem
order_audit <- tibble::tibble(
  coluna   = present_in_df005,
  pos_df005= match(present_in_df005, baseline_order),
  pos_df012= match(present_in_df005, names(df012)),
  ok       = pos_df012 == pos_df005
)
readr::write_tsv(order_audit, file.path(out_dir, "audit_Parte_8v3_order_df012_vs_df005.tsv"))
readr::write_lines(
  paste0("Order matches df005 for common columns? ",
         all(order_audit$ok, na.rm = TRUE)),
  file.path(out_dir, "audit_Parte_8v3_order_summary.txt")
)

## -----------------------------
## 9) Persistir
## -----------------------------
readr::write_tsv(df012, file.path(out_res, "df012.tsv"))

## ============================================================
## PART 9 — AGE IMPUTATION (mice::rf) + ROBUST AUDIT (refinado)
##   - Imputa por RF quando há preditores válidos dentro de 'type'
##   - Fallback: mediana dentro de 'type' (se houver observados)
##   - Auditoria por .rowid (não altera observados)
##   - Imputados SEM casas decimais (ceil)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(purrr); library(tibble); library(mice)
  library(forcats)
})

# Dependência explícita do RF do mice::rf
if (!requireNamespace("randomForest", quietly = TRUE)) {
  stop("O método mice::rf requer o pacote 'randomForest'. Instale-o: install.packages('randomForest')")
}

out_dir <- "harmonization_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## 0) Preconditions
## -----------------------------
stopifnot(exists("df012"))
req <- c("age_at_initial_pathologic_diagnosis","gender","type","tumor_status","vital_status")
missing_req <- setdiff(req, names(df012))
if (length(missing_req) > 0) {
  stop("Missing required columns: ", paste(missing_req, collapse = ", "))
}

## -----------------------------
## 1) Before snapshot + coerções
## -----------------------------
df013_before <- df012 %>%
  mutate(
    age_at_initial_pathologic_diagnosis = suppressWarnings(as.numeric(age_at_initial_pathologic_diagnosis)),
    gender       = as.factor(gender),
    tumor_status = as.factor(tumor_status),
    vital_status = as.factor(vital_status),
    .rowid = dplyr::row_number()
  )

## -----------------------------
## 2) Imputation by `type` with guarded RF
## -----------------------------
set.seed(123)
vars_rf <- c("age_at_initial_pathologic_diagnosis","gender","tumor_status","vital_status")

meta_list <- list()
method_per_row <- rep(NA_character_, nrow(df013_before))

# mínimo de casos completos por type para tentar RF
min_complete <- 10L  # ajuste conforme necessário

df013_after <- df013_before %>%
  split(.$type) %>%
  imap_dfr(function(sub_df, this_type) {
    
    idx_global   <- sub_df$.rowid
    n_before_na  <- sum(is.na(sub_df$age_at_initial_pathologic_diagnosis))
    sub_rf       <- sub_df %>% dplyr::select(all_of(vars_rf))
    
    # Sem NA → nada a fazer
    if (n_before_na == 0) {
      meta_list[[this_type]] <<- tibble(
        type = this_type, n_rows = nrow(sub_df),
        n_age_NA_before = 0L, n_age_imputed = 0L, n_age_NA_after = 0L,
        method_used = "none_needed",
        used_predictors = NA_character_,
        dropped_columns = NA_character_,
        error_message = NA_character_
      )
      return(sub_df)
    }
    
    # Remover colunas sem variação (após NA removidos)
    distinct_counts <- sub_rf %>% summarise(across(everything(), ~ dplyr::n_distinct(., na.rm = TRUE)))
    valid_cols <- names(distinct_counts)[as.integer(distinct_counts[1, ]) > 1]
    dropped_cols <- setdiff(names(sub_rf), valid_cols)
    sub_rf <- sub_rf %>% dplyr::select(all_of(valid_cols))
    
    preds <- setdiff(names(sub_rf), "age_at_initial_pathologic_diagnosis")
    
    # Fallback mediana (fechado)
    fallback_median <- function(err_msg = NA_character_, tag = "median_within_type") {
      age_obs <- sub_df$age_at_initial_pathologic_diagnosis[!is.na(sub_df$age_at_initial_pathologic_diagnosis)]
      if (length(age_obs) == 0) {
        imputed <- sub_df$age_at_initial_pathologic_diagnosis
        n_imp <- 0L
        list(values = imputed, n_imputed = n_imp,
             method = "fallback_none_no_observed_age",
             err = err_msg, used_preds = NA_character_)
      } else {
        med <- stats::median(age_obs, na.rm = TRUE)
        imputed <- sub_df$age_at_initial_pathologic_diagnosis
        tgt <- is.na(imputed)
        # >>> arredondar para cima apenas nos imputados
        imputed[tgt] <- ceiling(med)
        method_per_row[idx_global[tgt]] <<- "median_within_type"
        list(values = imputed, n_imputed = sum(tgt),
             method = tag, err = err_msg, used_preds = NA_character_)
      }
    }
    
    # Guardas antes de tentar RF
    if (length(preds) < 1) {
      res <- fallback_median(err_msg = "No valid predictors after zero-variance drop.",
                             tag = "rf_guard_fallback_median_within_type")
    } else {
      cc <- stats::complete.cases(sub_rf[, c("age_at_initial_pathologic_diagnosis", preds), drop = FALSE])
      n_cc <- sum(cc)
      if (n_cc < min_complete) {
        res <- fallback_median(err_msg = paste0("Too few complete cases for RF within type (n_cc=", n_cc, ")."),
                               tag = "rf_guard_fallback_median_within_type")
      } else {
        # Tentar mice::rf
        meth <- make.method(sub_rf)
        meth["age_at_initial_pathologic_diagnosis"] <- "rf"
        
        predM <- make.predictorMatrix(sub_rf)
        predM["age_at_initial_pathologic_diagnosis","age_at_initial_pathologic_diagnosis"] <- 0
        
        res <- tryCatch({
          imp  <- mice(data = sub_rf, method = meth, predictorMatrix = predM,
                       m = 1, seed = 123, print = FALSE)
          comp <- complete(imp, 1)
          age_imp <- comp$age_at_initial_pathologic_diagnosis
          
          # >>> arredondar para cima apenas nas posições originalmente NA
          tgt <- is.na(sub_df$age_at_initial_pathologic_diagnosis) & !is.na(age_imp)
          age_imp[tgt] <- ceiling(age_imp[tgt])
          
          method_per_row[idx_global[tgt]] <<- "rf"
          
          used_pred_names <- names(which(imp$predictorMatrix["age_at_initial_pathologic_diagnosis", ] == 1))
          used_pred_names <- if (length(used_pred_names) == 0) NA_character_ else paste(used_pred_names, collapse = " | ")
          
          list(values = age_imp,
               n_imputed = sum(tgt),
               method = "rf", err = NA_character_,
               used_preds = used_pred_names)
        }, error = function(e) {
          fb <- fallback_median(err_msg = conditionMessage(e),
                                tag = "rf_fallback_median_within_type")
          fb
        })
      }
    }
    
    meta_list[[this_type]] <<- tibble(
      type = this_type,
      n_rows = nrow(sub_df),
      n_age_NA_before = n_before_na,
      n_age_imputed   = res$n_imputed,
      n_age_NA_after  = sum(is.na(res$values)),
      method_used     = res$method,
      used_predictors = res$used_preds,
      dropped_columns = if (length(dropped_cols) == 0) NA_character_ else paste(dropped_cols, collapse = " | "),
      error_message   = res$err
    )
    
    sub_df$age_at_initial_pathologic_diagnosis <- res$values
    sub_df
  })

# rastreio do método por linha
df013_after$.__age_impute_method <- method_per_row

## -----------------------------
## 3) Audit — by .rowid
## -----------------------------
cmp_age <- df013_after %>%
  select(.rowid, type,
         age_after  = age_at_initial_pathologic_diagnosis,
         method_line = .__age_impute_method) %>%
  left_join(df013_before %>% select(.rowid, age_before = age_at_initial_pathologic_diagnosis),
            by = ".rowid") %>%
  mutate(
    age_before_is_na = is.na(age_before),
    age_after_is_na  = is.na(age_after),
    age_imputed      = age_before_is_na & !age_after_is_na,
    changed_obs      = !age_before_is_na & !age_after_is_na & (age_before != age_after)
  )

# (a) Observados não podem mudar
if (any(cmp_age$changed_obs, na.rm = TRUE)) {
  bad <- cmp_age %>% filter(changed_obs) %>% slice_head(n = 10)
  stop("Audit failure: observed (non-NA) ages were changed.\n",
       "Examples (rowid | type | before | after):\n",
       paste(apply(bad[, c(".rowid","type","age_before","age_after")], 1, paste, collapse = " | "),
             collapse = "\n"))
}

# (b) Não pode haver perda de informação
if (any(!cmp_age$age_before_is_na & cmp_age$age_after_is_na, na.rm = TRUE)) {
  bad <- cmp_age %>% filter(!age_before_is_na & age_after_is_na) %>% slice_head(n = 10)
  stop("Audit failure: information loss occurred (age became NA where it was not NA before).\n",
       "Examples (rowid | type | before | after):\n",
       paste(apply(bad[, c(".rowid","type","age_before","age_after")], 1, paste, collapse = " \n"))
  )
}

# (c) Plausibilidade
implaus_idx <- which(!is.na(cmp_age$age_after) & (cmp_age$age_after < 0 | cmp_age$age_after > 120))
if (length(implaus_idx) > 0) {
  warning("Ages outside plausible range [0,120]: ", length(implaus_idx),
          " row(s). Examples (rowid): ",
          paste(head(cmp_age$.rowid[implaus_idx], 10), collapse = ", "))
}

## -----------------------------
## 4) Reports
## -----------------------------
rowlevel <- cmp_age %>%
  filter(age_imputed) %>%
  transmute(
    rowid = .rowid, type, method_line,
    age_before, age_after
  )

meta_by_type <- bind_rows(meta_list)

by_type_counts <- cmp_age %>%
  group_by(type) %>%
  summarise(
    n_rows_calc          = dplyr::n(),
    n_age_NA_before_calc = sum(age_before_is_na),
    n_age_imputed_calc   = sum(age_imputed),
    n_age_NA_after_calc  = sum(age_after_is_na),
    .groups = "drop"
  )

by_type <- by_type_counts %>%
  left_join(meta_by_type, by = "type") %>%
  mutate(
    check_counts = (n_age_NA_before_calc == n_age_imputed_calc + n_age_NA_after_calc)
  ) %>%
  relocate(
    type,
    n_rows_calc, n_age_NA_before_calc, n_age_imputed_calc, n_age_NA_after_calc,
    check_counts, method_used, used_predictors, dropped_columns, error_message
  )

summary_tbl <- tibble(
  total_rows      = nrow(cmp_age),
  n_age_NA_before = sum(cmp_age$age_before_is_na),
  n_age_imputed   = sum(cmp_age$age_imputed),
  n_age_NA_after  = sum(cmp_age$age_after_is_na),
  check_balance   = (sum(cmp_age$age_before_is_na) ==
                       sum(cmp_age$age_imputed) + sum(cmp_age$age_after_is_na)),
  n_implausible_0_120 = length(implaus_idx)
)

readr::write_tsv(rowlevel,    file.path(out_dir, "audit_Parte_9_age_rowlevel.tsv"))
readr::write_tsv(by_type,     file.path(out_dir, "audit_Parte_9_age_by_type.tsv"))
readr::write_tsv(summary_tbl, file.path(out_dir, "audit_Parte_9_age_summary.tsv"))
readr::write_tsv(meta_by_type,file.path(out_dir, "audit_Parte_9_age_meta_by_type.tsv"))

## -----------------------------
## 5) Update working object
## -----------------------------
df013 <- df013_after %>%
  arrange(.rowid) %>%
  select(-.rowid, -.__age_impute_method)

cat("✓ Part 9 audit (mice::rf + guarded fallback; ceil nos imputados) completed.\n")
cat("→ Row-level: ", normalizePath(file.path(out_dir, "audit_Parte_9_age_rowlevel.tsv")),  "\n", sep = "")
cat("→ By-type:   ", normalizePath(file.path(out_dir, "audit_Parte_9_age_by_type.tsv")),   "\n", sep = "")
cat("→ Summary:   ", normalizePath(file.path(out_dir, "audit_Parte_9_age_summary.tsv")),   "\n", sep = "")
cat("→ Meta-type: ", normalizePath(file.path(out_dir, "audit_Parte_9_age_meta_by_type.tsv")),"\n", sep = "")

rio::export(df013,"data_table_results/df013.tsv")
df014 <- df013

## ============================================================
## PART 11 — IMPUTATION OF tumor_status (ranger) + AUDIT (robust if df015 missing)
## ============================================================
# 📌 The "permutation" metric was chosen to estimate predictor variable importance
# in the Random Forest model because it is less prone to bias, especially in binary 
# classification tasks such as imputing the tumor_status_bin ("WITH TUMOR" vs "TUMOR FREE").
#
# 👉 Alternatives:
# - importance = "impurity": uses purity gain (e.g., Gini index) at each split,
#   but tends to overestimate variables with many levels or imbalanced distributions. 
#   Not recommended when the goal is impartial evaluation of predictive relevance, 
#   especially for categorical variables.
#
# - importance = "impurity_corrected": applies statistical corrections to reduce 
#   the bias of the "impurity" metric. While computationally efficient and more 
#   reliable than the previous one, it still does not outperform the empirical 
#   robustness of "permutation".
#
# ✅ Therefore, "permutation" is the most appropriate choice in scenarios where 
# imputation accuracy and fair variable selection are essential, such as imputing 
# binary variables with clinical or diagnostic impact.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(tibble); library(ranger); library(purrr)
})

out_dir <- "harmonization_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## 0) Preconditions
## -----------------------------
stopifnot(exists("df014"))
req_pred <- c("age_at_initial_pathologic_diagnosis", "gender", "vital_status", "OS")
req_all  <- c(req_pred, "tumor_status")
miss_14  <- setdiff(req_all, names(df014))
if (length(miss_14) > 0) {
  stop("Missing columns in df014: ", paste(miss_14, collapse = ", "))
}

## -----------------------------
## 0.a) If df015 does not exist, generate it from the original pipeline
## -----------------------------
if (!exists("df015")) {
  message("df015 missing; generating from df014 for audit...")
  df015 <- df014 %>%
    mutate(
      tumor_status_bin = case_when(
        tumor_status == "WITH TUMOR" ~ 1,
        tumor_status == "TUMOR FREE" ~ 0,
        TRUE ~ NA_real_
      ),
      gender       = as.factor(gender),
      vital_status = as.factor(vital_status)
    )
  
  cols_preditoras <- c("age_at_initial_pathologic_diagnosis", "gender", "vital_status", "OS")
  
  train <- df015 %>%
    filter(!is.na(tumor_status_bin)) %>%
    select(all_of(cols_preditoras), tumor_status_bin) %>%
    filter(if_all(everything(), ~ !is.na(.)))
  
  test <- df015 %>%
    filter(is.na(tumor_status_bin)) %>%
    select(all_of(cols_preditoras))
  
  set.seed(123)
  modelo_rf_tmp <- ranger(
    formula        = tumor_status_bin ~ .,
    data           = train,
    num.trees      = 500,
    importance     = "permutation",
    probability    = FALSE,
    classification = TRUE
  )
  
  # Prediction (same as your block; if there are NA predictors in the test, attempt and warn)
  predictions <- tryCatch(
    predict(modelo_rf_tmp, data = test)$predictions,
    error = function(e) {
      warning("Prediction for test rows failed (possible NA in predictors). ",
              "Affected rows will remain NA in the imputation. Message: ", conditionMessage(e))
      rep(NA_real_, nrow(test))
    }
  )
  predictions <- ifelse(predictions >= 0.5, 1, 0)
  
  # Impute only where target was NA
  na_idx <- which(is.na(df015$tumor_status_bin))
  if (length(na_idx) > 0) {
    df015$tumor_status_bin[na_idx] <- predictions
  }
  
  df015 <- df015 %>%
    mutate(
      tumor_status = case_when(
        tumor_status_bin == 1 ~ "WITH TUMOR",
        tumor_status_bin == 0 ~ "TUMOR FREE",
        TRUE ~ NA_character_
      )
    ) %>%
    select(-tumor_status_bin)
}

## -----------------------------
## 1) Snapshots: BEFORE (df014) and AFTER (df015)
## -----------------------------
df14 <- df014 %>%
  mutate(
    .rowid = dplyr::row_number(),
    age_at_initial_pathologic_diagnosis = suppressWarnings(as.numeric(age_at_initial_pathologic_diagnosis)),
    gender       = as.factor(as.character(gender)),
    vital_status = as.factor(as.character(vital_status)),
    OS           = suppressWarnings(as.numeric(OS)),
    tumor_status_bin_before = dplyr::case_when(
      tumor_status == "WITH TUMOR" ~ 1,
      tumor_status == "TUMOR FREE" ~ 0,
      TRUE ~ NA_real_
    )
  )

df15 <- df015 %>%
  mutate(
    .rowid = dplyr::row_number(),
    tumor_status_bin_after = dplyr::case_when(
      tumor_status == "WITH TUMOR" ~ 1,
      tumor_status == "TUMOR FREE" ~ 0,
      TRUE ~ NA_real_
    )
  )

stopifnot(nrow(df14) == nrow(df15))  # same number of rows

## -----------------------------
## 2) EXPECTED reconstruction (model and predictions from df14)
## -----------------------------
cols_preditoras <- c("age_at_initial_pathologic_diagnosis", "gender", "vital_status", "OS")

base_exp <- df14 %>%
  select(.rowid, all_of(cols_preditoras), tumor_status_bin_before)

train <- base_exp %>%
  filter(!is.na(tumor_status_bin_before)) %>%
  select(all_of(cols_preditoras), tumor_status_bin = tumor_status_bin_before) %>%
  filter(if_all(everything(), ~ !is.na(.)))

test <- base_exp %>%
  filter(is.na(tumor_status_bin_before)) %>%
  select(.rowid, all_of(cols_preditoras))

test_incomplete <- test %>%
  mutate(n_na = rowSums(is.na(across(all_of(cols_preditoras))))) %>%
  filter(n_na > 0)
n_test_incomplete <- nrow(test_incomplete)

set.seed(123)
modelo_rf_exp <- ranger(
  formula        = tumor_status_bin ~ .,
  data           = train,
  num.trees      = 500,
  importance     = "permutation",
  probability    = FALSE,
  classification = TRUE
)

if (nrow(test) > 0) {
  pred_raw <- tryCatch(
    predict(modelo_rf_exp, data = select(test, all_of(cols_preditoras)))$predictions,
    error = function(e) {
      warning("EXPECTED prediction failed (likely NA in test predictors). ",
              "expected_bin will remain NA for these rows. Message: ", conditionMessage(e))
      rep(NA_real_, nrow(test))
    }
  )
  pred_bin <- ifelse(pred_raw >= 0.5, 1, 0)
  
  expected_pred <- tibble(
    .rowid                 = test$.rowid,
    expected_bin           = as.numeric(pred_bin),
    expected_tumor_status  = ifelse(expected_bin == 1, "WITH TUMOR", "TUMOR FREE")
  )
} else {
  expected_pred <- tibble(.rowid = integer(), expected_bin = numeric(), expected_tumor_status = character())
}

## -----------------------------
## 3) Comparison OBSERVED vs EXPECTED (only imputed rows)
## -----------------------------
cmp <- df14 %>%
  select(.rowid, type, all_of(cols_preditoras),
         tumor_status_before = tumor_status, tumor_status_bin_before) %>%
  left_join(df15 %>% select(.rowid, tumor_status_after = tumor_status, tumor_status_bin_after), by = ".rowid") %>%
  left_join(expected_pred, by = ".rowid") %>%
  mutate(
    was_imputed_target = is.na(tumor_status_bin_before),
    eq_imputed_bin     = ifelse(was_imputed_target, tumor_status_bin_after == expected_bin, NA),
    eq_imputed_label   = ifelse(was_imputed_target, tumor_status_after     == expected_tumor_status, NA),
    changed_observed   = !was_imputed_target & !is.na(tumor_status_before) &
      (tumor_status_before != tumor_status_after)
  )

## -----------------------------
## 4) Formal checks
## -----------------------------
# (a) Observed (non-NA before) must not change
if (any(cmp$changed_observed, na.rm = TRUE)) {
  bad <- cmp %>% filter(changed_observed) %>% slice_head(n = 10)
  stop("Audit failure: observed 'tumor_status' (non-NA) was altered.\n",
       "Examples (rowid | type | before | after):\n",
       paste(apply(bad[, c(".rowid","type","tumor_status_before","tumor_status_after")], 1, paste, collapse = " | "),
             collapse = "\n"))
}

# (b) Observed vs. expected agreement in imputed rows (where expected is defined)
mismatch <- cmp %>%
  filter(was_imputed_target & !is.na(expected_tumor_status) &
           (is.na(eq_imputed_label) | !eq_imputed_label))
if (nrow(mismatch) > 0) {
  diag <- mismatch %>% slice_head(n = 10)
  stop("Audit failure: ", nrow(mismatch),
       " imputed row(s) diverge from expected (same model/seed).\n",
       "Examples (rowid | type | expected | observed):\n",
       paste(apply(diag[, c(".rowid","type","expected_tumor_status","tumor_status_after")], 1, paste, collapse = " | "),
             collapse = "\n"))
}

# (c) Label domain
valid_labels <- c("WITH TUMOR","TUMOR FREE")
invalid_after <- setdiff(unique(na.omit(cmp$tumor_status_after)), valid_labels)
if (length(invalid_after) > 0) {
  stop("Audit failure: labels outside expected domain: ",
       paste(invalid_after, collapse = " | "))
}

## -----------------------------
## 5) Model metrics (training)
## -----------------------------
oob_err <- modelo_rf_exp$prediction.error
conf_mat <- if (!is.null(modelo_rf_exp$confusion.matrix)) {
  as.data.frame.matrix(modelo_rf_exp$confusion.matrix)
} else NULL

varimp <- tibble(
  variable = names(modelo_rf_exp$variable.importance),
  importancePermutation = as.numeric(modelo_rf_exp$variable.importance)
) %>% arrange(desc(importancePermutation))

## -----------------------------
## 6) Reports
## -----------------------------
audit_rowlevel <- cmp %>%
  filter(was_imputed_target) %>%
  transmute(
    rowid = .rowid, type,
    age_at_initial_pathologic_diagnosis, gender, vital_status, OS,
    tumor_status_before,
    expected_tumor_status, tumor_status_after,
    match_expected = eq_imputed_label
  )

audit_by_type <- audit_rowlevel %>%
  count(type, tumor_status_after, name = "n") %>%
  group_by(type) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(type, desc(n))

audit_summary <- tibble(
  total_rows             = nrow(cmp),
  n_imputed              = sum(cmp$was_imputed_target),
  n_imputed_match        = sum(cmp$was_imputed_target & cmp$eq_imputed_label, na.rm = TRUE),
  n_imputed_mismatch     = sum(cmp$was_imputed_target & !cmp$eq_imputed_label, na.rm = TRUE),
  n_changed_observed     = sum(cmp$changed_observed, na.rm = TRUE),
  test_incomplete_rows   = n_test_incomplete,
  oob_error              = oob_err
)

if (!is.null(conf_mat)) {
  conf_tbl <- conf_mat %>%
    tibble::rownames_to_column(var = "predicted") %>%
    as_tibble() %>%
    pivot_longer(cols = -predicted, names_to = "true", values_to = "n")
  readr::write_tsv(conf_tbl, file.path(out_dir, "audit_Part_11_tumorStatus_confusion.tsv"))
}

readr::write_tsv(audit_rowlevel, file.path(out_dir, "audit_Part_11_tumorStatus_rowlevel.tsv"))
readr::write_tsv(audit_by_type,  file.path(out_dir, "audit_Part_11_tumorStatus_by_type.tsv"))
readr::write_tsv(audit_summary,  file.path(out_dir, "audit_Part_11_tumorStatus_summary.tsv"))
readr::write_tsv(varimp,         file.path(out_dir, "audit_Part_11_tumorStatus_varimp.tsv"))

## -----------------------------
## 7) Observation about test completeness
## -----------------------------
if (n_test_incomplete > 0) {
  warning("There are ", n_test_incomplete, " row(s) in the test set with NA in some predictor. ",
          "Your original block does not remove NAs in the test; in different datasets this may cause predict() to fail. ",
          "Options: filter incomplete test rows or impute predictors beforehand.")
}

cat("✓ Part 11 audit completed.\n")
cat("→ Row-level: ", normalizePath(file.path(out_dir, "audit_Part_11_tumorStatus_rowlevel.tsv")), "\n", sep = "")
cat("→ Summary:   ", normalizePath(file.path(out_dir, "audit_Part_11_tumorStatus_summary.tsv")),  "\n", sep = "")
cat("→ By-type:   ", normalizePath(file.path(out_dir, "audit_Part_11_tumorStatus_by_type.tsv")),  "\n", sep = "")
cat("→ VarImp:    ", normalizePath(file.path(out_dir, "audit_Part_11_tumorStatus_varimp.tsv")),   "\n", sep = "")
if (!is.null(conf_mat)) {
  cat("→ Confusion: ", normalizePath(file.path(out_dir, "audit_Part_11_tumorStatus_confusion.tsv")), "\n", sep = "")
}
rio::export(df015,"data_table_results/df015.tsv")

## ============================================================
## PART 12 — IMPUTATION of initial_pathologic_dx_year (ranger)
##            sem round-trip / sem limpeza global
##            + auditoria "freeze" das colunas 31+
## ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
  library(tibble); library(ranger); library(purrr)
})

out_dir <- "harmonization_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## 0) Preconditions
## -----------------------------
stopifnot(exists("df015"))
req_pred <- c("age_at_initial_pathologic_diagnosis","gender","vital_status","tumor_status","type")
req_all  <- c(req_pred, "initial_pathologic_dx_year")
miss_15  <- setdiff(req_all, names(df015))
if (length(miss_15) > 0) stop("Missing columns in df015: ", paste(miss_15, collapse = ", "))

## -----------------------------
## 0.a) Construção de df016 sem alterar outras colunas
## -----------------------------
# Cópia direta: não re-infere tipos, não limpa globalmente
df016 <- df015

# --- Congelamento das colunas 31+ (exceto a coluna-alvo) para auditoria:
lock_from <- 31L
n_cols    <- ncol(df016)
locked_idx <- if (n_cols >= lock_from) lock_from:n_cols else integer(0)
locked_cols <- names(df016)[locked_idx]
locked_cols_excl_target <- setdiff(locked_cols, "initial_pathologic_dx_year")

snapshot_vals_31p_before  <- if (length(locked_cols_excl_target) > 0) df016[, locked_cols_excl_target, drop = FALSE] else NULL
snapshot_class_31p_before <- if (length(locked_cols_excl_target) > 0) vapply(df016[locked_cols_excl_target], function(x) paste(class(x), collapse="|"), character(1)) else NULL

## -----------------------------
## 1) Data frame auxiliar de modelagem (sem tocar df016)
## -----------------------------
# As coerções/fatores são feitas APENAS no frame auxiliar 'mf', não em df016
mf <- df016 %>%
  mutate(
    .rowid = dplyr::row_number(),
    y_year = suppressWarnings(as.numeric(initial_pathologic_dx_year)),
    age_num = suppressWarnings(as.numeric(age_at_initial_pathologic_diagnosis)),
    gender_f = factor(as.character(gender)),
    vital_f  = factor(as.character(vital_status)),
    tumor_f  = factor(as.character(tumor_status)),
    type_f   = factor(as.character(type))
  )

pred_cols <- c("age_num","gender_f","vital_f","tumor_f","type_f")

# Conjunto de treino (alvo observado + preditores completos)
train_idx <- which(!is.na(mf$y_year) & stats::complete.cases(mf[, c("y_year", pred_cols), drop = FALSE]))
# Conjunto de teste (alvo NA + preditores completos)
test_idx  <- which(is.na(mf$y_year)  & stats::complete.cases(mf[, pred_cols, drop = FALSE]))

## -----------------------------
## 2) Treino e predição
## -----------------------------
if (length(test_idx) > 0 && length(train_idx) > 0) {
  set.seed(123)
  modelo_reg <- ranger(
    formula        = y_year ~ .,
    data           = mf[train_idx, c("y_year", pred_cols), drop = FALSE],
    num.trees      = 500,
    importance     = "permutation",
    classification = FALSE
  )
  preds <- predict(modelo_reg, data = mf[test_idx, pred_cols, drop = FALSE])$predictions
  preds_int <- as.integer(round(preds))
  
  # (opcional) sanidade de faixa plausível
  plaus_min <- 1900L; plaus_max <- 2025L
  preds_int[preds_int < plaus_min | preds_int > plaus_max] <- preds_int[preds_int < plaus_min | preds_int > plaus_max]
  
  # Mapear pelos .rowid (evita desalinhamento)
  target_rows <- mf$.rowid[test_idx]
  df016$initial_pathologic_dx_year[target_rows] <- preds_int
} else {
  message("Sem linhas elegíveis para imputação (ou sem treino ou sem teste completos).")
}

## -----------------------------
## 3) Auditoria — não alterar colunas 31+ (exceto alvo)
## -----------------------------
if (length(locked_cols_excl_target) > 0) {
  snapshot_vals_31p_after  <- df016[, locked_cols_excl_target, drop = FALSE]
  snapshot_class_31p_after <- vapply(df016[locked_cols_excl_target], function(x) paste(class(x), collapse="|"), character(1))
  
  # (a) Conteúdo idêntico
  if (!identical(snapshot_vals_31p_before, snapshot_vals_31p_after)) {
    # Diagnóstico rápido
    diff_cols <- locked_cols_excl_target[vapply(locked_cols_excl_target, function(cn) !identical(snapshot_vals_31p_before[[cn]], snapshot_vals_31p_after[[cn]]), logical(1))]
    stop("Audit failure (freeze 31+): valores alterados em colunas fora do alvo: ",
         paste(head(diff_cols, 10), collapse = ", "),
         if (length(diff_cols) > 10) "... (mais colunas).")
  }
  # (b) Classe idêntica
  if (!identical(snapshot_class_31p_before, snapshot_class_31p_after)) {
    diff_cls <- names(which(snapshot_class_31p_before != snapshot_class_31p_after))
    stop("Audit failure (freeze 31+): classes alteradas em colunas fora do alvo: ",
         paste(head(diff_cls, 10), collapse = ", "),
         if (length(diff_cls) > 10) "... (mais colunas).")
  }
}

## -----------------------------
## 4) Relatórios mínimos (opcional)
## -----------------------------
df016_before <- df015 %>% transmute(.rowid = dplyr::row_number(),
                                    initial_pathologic_dx_year_before = suppressWarnings(as.numeric(initial_pathologic_dx_year)))
df016_after  <- df016 %>% transmute(.rowid = dplyr::row_number(),
                                    initial_pathologic_dx_year_after  = suppressWarnings(as.numeric(initial_pathologic_dx_year)))

cmp <- df016_before %>%
  left_join(df016_after, by = ".rowid") %>%
  mutate(
    was_na_before = is.na(initial_pathologic_dx_year_before),
    was_imputed   = was_na_before & !is.na(initial_pathologic_dx_year_after),
    changed_obs   = !was_na_before & (initial_pathologic_dx_year_before != initial_pathologic_dx_year_after)
  )

if (any(cmp$changed_obs, na.rm = TRUE)) {
  bad <- cmp %>% filter(changed_obs) %>% slice_head(n = 10)
  stop("Audit failure: anos observados foram alterados.\n",
       paste(apply(bad[, c(".rowid","initial_pathologic_dx_year_before","initial_pathologic_dx_year_after")], 1, paste, collapse = " | "),
             collapse = "\n"))
}

readr::write_tsv(
  tibble(
    total_rows    = nrow(cmp),
    n_na_before   = sum(cmp$was_na_before),
    n_imputed     = sum(cmp$was_imputed),
    n_changed_obs = sum(cmp$changed_obs, na.rm = TRUE)
  ),
  file.path(out_dir, "audit_Parte_12_dxyear_summary.tsv")
)

cat("✓ Part 12 (sem round-trip; freeze 31+) concluída.\n")

rio::export(df016,"data_table_results/df016.tsv")



## ============================================================
## PART 13 — HARMONIZATION of new_tumor_event_type/_site + AUDIT
##   Extensões: integra PFI, PFI.1, PFI.2 e PFS com prioridade
##   Conservador: apenas preenche NA; não altera valores não-NA
## ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(tibble); library(readr); library(readxl)
})

out_dir <- "harmonization_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## 0) Preconditions
## -----------------------------
stopifnot(exists("df016"))
req_cols <- c("type","new_tumor_event_type",
              "new_tumor_event_site","new_tumor_event_site_other")
miss_16  <- setdiff(req_cols, names(df016))
if (length(miss_16) > 0) stop("Missing columns in df016: ", paste(miss_16, collapse=", "))

## colunas de endpoint tratadas como opcionais (se ausentes, criamos NA)
opt_endpoints <- c("DFI","PFI","PFI.1","PFI.2","PFS")
for (nm in opt_endpoints) {
  if (!nm %in% names(df016)) df016[[nm]] <- NA_real_
}

## -----------------------------
## 1) HARMONIZATION (work on tmp017; promote to df017 only after guard)
## -----------------------------
tmp017 <- df016

## Normalizador de strings “vazia”
is_empty <- function(x) is.na(x) | x == ""

## 1A) DFI==1 & NA(new_tumor_event_type) -> "Recurrence" (recorrência verdadeira)
tmp017 <- tmp017 %>%
  mutate(
    new_tumor_event_type = case_when(
      is_empty(new_tumor_event_type) & !is.na(DFI) & DFI == 1 ~ "Recurrence",
      TRUE ~ new_tumor_event_type
    )
  )

## 1B) PFI-family (PFI, PFI.1, PFI.2) & NA(event_type) & DFI != 1
##     → imputação no domínio permitido (reprodutível; semente fixa)
eventos_possiveis <- c(
  "Distant Metastasis","Locoregional Recurrence","New Primary Tumor",
  "Metastatic","Locoregional Disease","Progression of Disease","Recurrence"
)

## semente fixa p/ reprodutibilidade (evita set.seed(Sys.time()))
set.seed(777)

## flag: qualquer PFI* == 1
pfifam_signal <- with(tmp017, (PFI == 1) | (`PFI.1` == 1) | (`PFI.2` == 1))

idx_pfi_impute <- which(
  is_empty(tmp017$new_tumor_event_type) &
    !isTRUE(tmp017$DFI == 1) &
    isTRUE(pfifam_signal)
)

if (length(idx_pfi_impute) > 0) {
  tmp017$new_tumor_event_type[idx_pfi_impute] <-
    sample(eventos_possiveis, length(idx_pfi_impute), replace = TRUE)
}

## 1C) PFS==1 & NA(event_type) & (nenhum PFI* acionado) & DFI != 1
##     → “Progression of Disease” (determinístico; não aleatório)
idx_pfs_impute <- which(
  is_empty(tmp017$new_tumor_event_type) &
    !isTRUE(tmp017$DFI == 1) &
    isTRUE(tmp017$PFS == 1) &
    !isTRUE(pfifam_signal)        # prioridade: só PFS se ninguém da família PFI acionou
)

if (length(idx_pfs_impute) > 0) {
  tmp017$new_tumor_event_type[idx_pfs_impute] <- "Progression of Disease"
}

## 1D) SITE rules (locoregional e mapeamentos específicos), apenas quando site vazio
eventos_locorregionais <- c(
  "Locoregional Recurrence","Locoregional Disease","Recurrence","Biochemical evidence of disease"
)

tmp017 <- tmp017 %>%
  mutate(
    new_tumor_event_site = case_when(
      new_tumor_event_type %in% eventos_locorregionais &
        (is_empty(new_tumor_event_site)) ~ type,
      
      new_tumor_event_type == "New primary melanoma" &
        (is_empty(new_tumor_event_site)) ~ "SKIM",
      
      new_tumor_event_type == "Regional lymph node" &
        (is_empty(new_tumor_event_site)) ~ "Lymph Node",
      
      new_tumor_event_type == "Intrahepatic Recurrence" &
        (is_empty(new_tumor_event_site)) ~ "Liver",
      
      new_tumor_event_type == "Locoregional (Urothelial tumor event)" &
        (is_empty(new_tumor_event_site)) ~ "Bladder",
      
      new_tumor_event_type == "Intrapleural Progression" &
        (is_empty(new_tumor_event_site)) ~ "Ipsilateral Pleura",
      
      new_tumor_event_type == "Progression of Disease" &
        (is_empty(new_tumor_event_site)) ~ "No event",
      
      TRUE ~ new_tumor_event_site
    )
  )

## 1E) Extrahepatic Recurrence -> usar moda de sites já observados
moda_site_extrahepatic <- tmp017 %>%
  filter(new_tumor_event_type == "Extrahepatic Recurrence",
         !is_empty(new_tumor_event_site)) %>%
  count(new_tumor_event_site, sort = TRUE) %>%
  head(1) %>% pull(new_tumor_event_site)

tmp017 <- tmp017 %>%
  mutate(
    new_tumor_event_site = case_when(
      new_tumor_event_type == "Extrahepatic Recurrence" &
        (is_empty(new_tumor_event_site)) ~ moda_site_extrahepatic,
      TRUE ~ new_tumor_event_site
    )
  )

## 1F) Converter códigos TCGA -> órgãos (SITE)
tcga_to_organ <- c(
  LAML="Bone marrow", ACC="Adrenal gland", BLCA="Bladder", LGG="Brain",
  BRCA="Breast", CESC="Cervix", CHOL="Intrahepatic bile ducts",
  LCML="Bone marrow and peripheral blood", COAD="Colon", CNTL="Various normal tissues",
  ESCA="Esophagus", GBM="Brain", HNSC="Head and neck region",
  KICH="Kidney", KIRC="Kidney", KIRP="Kidney",
  LIHC="Liver", LUAD="Lung", LUSC="Lung", DLBC="Lymphatic system",
  MESO="Pleura", MISC="Various", OV="Ovary", PAAD="Pancreas",
  PCPG="Adrenal medulla and paraganglia", PRAD="Prostate", READ="Rectum",
  SARC="Connective tissue", SKCM="Skin", STAD="Stomach", TGCT="Testis",
  THYM="Thymus", THCA="Thyroid gland", UCS="Uterus", UCEC="Endometrium", UVM="Uvea"
)

tmp017 <- tmp017 %>%
  mutate(
    new_tumor_event_site = if_else(
      !is_empty(new_tumor_event_site) & new_tumor_event_site %in% names(tcga_to_organ),
      unname(tcga_to_organ[new_tumor_event_site]),
      new_tumor_event_site
    )
  )

## 1G) Padronização de *_site_other
tmp017 <- tmp017 %>%
  mutate(
    new_tumor_event_site_other = case_when(
      is_empty(new_tumor_event_site) ~ "No event",
      !str_detect(str_to_lower(str_squish(new_tumor_event_site)), "^other, specify") ~ "No event",
      TRUE ~ new_tumor_event_site_other
    )
  )

## 1H) (Opcional) Imputação de *_site_other quando site == "Other, specify"
if (file.exists("TCGA_Metastasis_Complete_Table.xlsx")) {
  cancer_type_ontology <- list(
    Epithelial     = c("LUAD","LUSC","BRCA","BLCA","UCEC","UCS","CESC","ESCA","STAD","COAD","READ","PRAD","OV","THCA","HNSC"),
    HepatoBiliary  = c("LIHC","CHOL","PAAD"),
    Renal          = c("KIRP","KIRC","KICH"),
    Neuroendocrine = c("ACC","PCPG"),
    CNS_Glia       = c("GBM","LGG"),
    Mesenchymal    = c("SARC","MESO"),
    Hematologic    = c("DLBC","THYM","LAML"),
    GermCell       = c("TGCT"),
    Melanocytic    = c("SKCM","UVM")
  )
  code_to_group <- map_dfr(names(cancer_type_ontology), function(gr) {
    tibble(code = cancer_type_ontology[[gr]], group = gr)
  }) %>% deframe()
  
  met_raw <- read_excel("TCGA_Metastasis_Complete_Table.xlsx")
  met_use <- met_raw %>%
    filter(`Metastatic Site` != "no metastasis data found") %>%
    mutate(
      prop_raw = suppressWarnings(as.numeric(`Estimated Proportion (%)`)),
      prop_raw = ifelse(is.na(prop_raw),
                        suppressWarnings(as.numeric(str_replace_all(`Estimated Proportion (%)`, "[^0-9\\.]", ""))),
                        prop_raw),
      prop = ifelse(prop_raw > 1, prop_raw/100, prop_raw)
    ) %>%
    filter(!is.na(prop), prop > 0) %>%
    transmute(code = str_trim(`TCGA Code`),
              site = str_trim(`Metastatic Site`),
              prop)
  
  met_by_code <- met_use %>%
    group_by(code, site) %>% summarise(p = sum(prop), .groups="drop") %>%
    group_by(code) %>%
    summarise(destinos = list(site), pesos_raw = list(p), soma = sum(p), .groups="drop") %>%
    mutate(
      residual = pmax(0, 1 - soma),
      destinos = map2(destinos, residual, ~ if (.y > 0) c(.x,"__FALLBACK_GROUP__") else .x),
      pesos    = map2(pesos_raw, residual, ~ if (.y > 0) c(.x, .y) else (.x / sum(.x)))
    )
  
  metastasis_table <- met_by_code %>%
    transmute(type = code, organ = destinos, weight = pesos) %>%
    tidyr::unnest(c(organ, weight)) %>%
    filter(organ != "__FALLBACK_GROUP__", is.finite(weight), weight > 0) %>%
    group_by(type) %>% mutate(weight = weight / sum(weight)) %>% ungroup()
  
  ont_groups <- cancer_type_ontology
  
  impute_random_site <- function(tcga_code) {
    dir_tbl <- metastasis_table %>% filter(type == tcga_code)
    if (nrow(dir_tbl) > 0) {
      return(sample(dir_tbl$organ, size = 1, prob = dir_tbl$weight))
    }
    grp <- names(ont_groups)[vapply(ont_groups, function(g) tcga_code %in% g, logical(1))]
    if (length(grp) > 0) {
      grp_types <- ont_groups[[grp]]
      grp_tbl <- metastasis_table %>%
        filter(type %in% grp_types) %>%
        group_by(organ) %>% summarise(w = sum(weight), .groups="drop") %>%
        mutate(w = w / sum(w)) %>% filter(is.finite(w), w > 0)
      if (nrow(grp_tbl) > 0) return(sample(grp_tbl$organ, size = 1, prob = grp_tbl$w))
    }
    return(NA_character_)
  }
  
  set.seed(777)
  idx_pred <- which(
    str_to_lower(str_squish(replace_na(tmp017$new_tumor_event_site, ""))) == "other, specify" &
      (is_empty(tmp017$new_tumor_event_site_other))
  )
  if (length(idx_pred) > 0) {
    tmp017$new_tumor_event_site_other[idx_pred] <- vapply(
      tmp017$type[idx_pred], impute_random_site, FUN.VALUE = character(1)
    )
  }
}

## 1I) Safe cleanup: trocar "" por NA em colunas character
tmp017[] <- lapply(tmp017, function(x) {
  if (is.character(x)) { x[x == ""] <- NA; return(x) } else x
})

## Integrity guard
stopifnot(nrow(tmp017) == nrow(df016))
df017 <- tmp017
df018 <- df017  # se você usa df018 como “dataset after part 13”


## -----------------------------
## 2) AUDIT (df016 vs df017), incluindo gatilhos DFI/PFI*/PFS
## -----------------------------
bef <- df016 %>%
  transmute(
    .rowid = dplyr::row_number(),
    type = as.character(type),
    DFI = suppressWarnings(as.numeric(DFI)),
    PFI = suppressWarnings(as.numeric(PFI)),
    PFI.1 = suppressWarnings(as.numeric(`PFI.1`)),
    PFI.2 = suppressWarnings(as.numeric(`PFI.2`)),
    PFS = suppressWarnings(as.numeric(PFS)),
    evt_type_before = as.character(new_tumor_event_type),
    site_before     = as.character(new_tumor_event_site),
    site_other_before = as.character(new_tumor_event_site_other)
  )

aft <- df017 %>%
  transmute(
    .rowid = dplyr::row_number(),
    type = as.character(type),
    DFI = suppressWarnings(as.numeric(DFI)),
    PFI = suppressWarnings(as.numeric(PFI)),
    PFI.1 = suppressWarnings(as.numeric(`PFI.1`)),
    PFI.2 = suppressWarnings(as.numeric(`PFI.2`)),
    PFS = suppressWarnings(as.numeric(PFS)),
    evt_type_after = as.character(new_tumor_event_type),
    site_after     = as.character(new_tumor_event_site),
    site_other_after = as.character(new_tumor_event_site_other)
  )

stopifnot(nrow(bef) == nrow(aft))

cmp <- bef %>%
  left_join(aft, by = c(".rowid","type","DFI","PFI","PFI.1","PFI.2","PFS")) %>%
  mutate(
    was_na_evt      = is_empty(evt_type_before),
    was_na_site     = is_empty(site_before),
    is_na_site_aft  = is_empty(site_after),
    changed_type    = evt_type_before != evt_type_after & !(is.na(evt_type_before) & is.na(evt_type_after)),
    changed_site    = site_before     != site_after     & !(is.na(site_before)     & is.na(site_after)),
    
    pfi_family_1    = (PFI == 1) | (PFI.1 == 1) | (PFI.2 == 1),
    trigger_DFI     = was_na_evt & (DFI == 1) & (evt_type_after == "Recurrence"),
    trigger_PFIfam  = was_na_evt & !isTRUE(DFI == 1) & pfi_family_1 & (evt_type_after %in% !!eventos_possiveis),
    trigger_PFS     = was_na_evt & !isTRUE(DFI == 1) & !pfi_family_1 & (PFS == 1) & (evt_type_after == "Progression of Disease")
  )

## Domínio e invariantes
eventos_locorregionais <- c(
  "Locoregional Recurrence","Locoregional Disease","Recurrence","Biochemical evidence of disease"
)
tcga_codes <- names(tcga_to_organ)

# C1: DFI==1 & NA(evt_before) -> 'Recurrence'
c1_idx <- which(cmp$trigger_DFI)
c1_ok  <- length(c1_idx) == sum(cmp$DFI == 1 & cmp$was_na_evt, na.rm = TRUE)

# C2: PFI-family imputou apenas dentro do domínio permitido
c2_idx <- which(cmp$trigger_PFIfam)
c2_ok  <- length(c2_idx) == sum(cmp$was_na_evt & !isTRUE(cmp$DFI == 1) & cmp$pfi_family_1, na.rm = TRUE) &&
  all(cmp$evt_type_after[c2_idx] %in% eventos_possiveis)

# C3: PFS imputou exatamente “Progression of Disease”
c3_idx <- which(cmp$trigger_PFS)
c3_ok  <- length(c3_idx) == sum(cmp$was_na_evt & !isTRUE(cmp$DFI == 1) & !cmp$pfi_family_1 & cmp$PFS == 1, na.rm = TRUE) &&
  all(cmp$evt_type_after[c3_idx] == "Progression of Disease")

# C4: Não alterar tipo quando não era alvo (não-NA before)
c4_idx <- which(!cmp$was_na_evt)
c4_ok  <- length(c4_idx) == 0 || all(
  cmp$evt_type_after[c4_idx] == cmp$evt_type_before[c4_idx] |
    (is.na(cmp$evt_type_before[c4_idx]) & is.na(cmp$evt_type_after[c4_idx]))
)

# C5: Locorregional com site vazio antes -> site preenchido e sem código TCGA cru
locoreg_mask <- cmp$evt_type_after %in% eventos_locorregionais & cmp$was_na_site
c5_idx <- which(locoreg_mask)
c5_ok_nonempty <- length(c5_idx) == 0 || all(!cmp$is_na_site_aft[c5_idx])
c5_no_tcga_left <- length(c5_idx) == 0 || all(!(cmp$site_after[c5_idx] %in% tcga_codes))

# C6: Regras específicas (event->site)
rule_expect <- tribble(
  ~event,                          ~site_expected,
  "New primary melanoma",          "SKIM",
  "Regional lymph node",           "Lymph Node",
  "Intrahepatic Recurrence",       "Liver",
  "Locoregional (Urothelial tumor event)", "Bladder",
  "Intrapleural Progression",      "Ipsilateral Pleura",
  "Progression of Disease",        "No event"
)
c6_checks <- rule_expect %>%
  mutate(
    n_applicable = map2_int(event, site_expected, ~{
      idx <- which(cmp$evt_type_after == .x & cmp$was_na_site)
      length(idx)
    }),
    n_ok = map2_int(event, site_expected, ~{
      idx <- which(cmp$evt_type_after == .x & cmp$was_na_site)
      if (length(idx) == 0) return(0L)
      sum(cmp$site_after[idx] == .y, na.rm = TRUE)
    }),
    ok = (n_applicable == n_ok)
  )

# C7: Extrahepatic Recurrence com site vazio antes -> site preenchido
c7_idx <- which(cmp$evt_type_after == "Extrahepatic Recurrence" & cmp$was_na_site)
c7_ok_nonempty <- length(c7_idx) == 0 || all(!cmp$is_na_site_aft[c7_idx])

# C8/C9: *_site_other padronizado
low_site <- str_to_lower(str_squish(replace_na(cmp$site_after, "")))
c8_idx_empty <- which(is_empty(cmp$site_after))
c8_ok_empty  <- length(c8_idx_empty) == 0 ||
  all(replace_na(cmp$site_other_after[c8_idx_empty], "") == "No event")
c8_idx_not_other <- which(!str_detect(low_site, "^other, specify") & !is_empty(cmp$site_after))
c8_ok_not_other  <- length(c8_idx_not_other) == 0 ||
  all(replace_na(cmp$site_other_after[c8_idx_not_other], "") == "No event")

## Reports
rowlevel <- cmp %>%
  mutate(
    trigger = case_when(
      trigger_DFI    ~ "DFI",
      trigger_PFIfam ~ "PFI_family",
      trigger_PFS    ~ "PFS",
      TRUE ~ NA_character_
    ),
    rule_applied = case_when(
      trigger_DFI    ~ "DFI==1→type='Recurrence'",
      trigger_PFIfam ~ "PFI*==1 & NA→type∈domain",
      trigger_PFS    ~ "PFS==1 & NA→'Progression of Disease'",
      TRUE ~ NA_character_
    ),
    site_rule_expected = case_when(
      evt_type_after == "New primary melanoma" & was_na_site ~ "SKIM",
      evt_type_after == "Regional lymph node"  & was_na_site ~ "Lymph Node",
      evt_type_after == "Intrahepatic Recurrence" & was_na_site ~ "Liver",
      evt_type_after == "Locoregional (Urothelial tumor event)" & was_na_site ~ "Bladder",
      evt_type_after == "Intrapleural Progression" & was_na_site ~ "Ipsilateral Pleura",
      evt_type_after == "Progression of Disease" & was_na_site ~ "No event",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(changed_type | changed_site | !is.na(rule_applied) | !is.na(site_rule_expected)) %>%
  select(.rowid, type, DFI, PFI, `PFI.1`, `PFI.2`, PFS,
         evt_type_before, evt_type_after,
         site_before, site_after,
         site_other_before, site_other_after,
         trigger, rule_applied, site_rule_expected)

by_type <- cmp %>%
  transmute(type,
            imputed_type_DFI = trigger_DFI,
            imputed_type_PFI = trigger_PFIfam,
            imputed_type_PFS = trigger_PFS,
            site_locoreg     = (evt_type_after %in% eventos_locorregionais & was_na_site & !is_na_site_aft),
            site_specific    = (was_na_site & evt_type_after %in% rule_expect$event & !is_na_site_aft)) %>%
  group_by(type) %>%
  summarise(
    n = n(),
    n_imputed_type_DFI = sum(imputed_type_DFI, na.rm = TRUE),
    n_imputed_type_PFI = sum(imputed_type_PFI, na.rm = TRUE),
    n_imputed_type_PFS = sum(imputed_type_PFS, na.rm = TRUE),
    n_site_locoreg     = sum(site_locoreg, na.rm = TRUE),
    n_site_specific    = sum(site_specific, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(type)

audit_summary <- tibble(
  total = nrow(cmp),
  # triggers
  n_trigger_DFI      = sum(cmp$trigger_DFI,    na.rm = TRUE),
  n_trigger_PFI_fam  = sum(cmp$trigger_PFIfam, na.rm = TRUE),
  n_trigger_PFS      = sum(cmp$trigger_PFS,    na.rm = TRUE),
  # checks
  c1_DFI_all_targeted_became_Recurrence = c1_ok,
  c2_PFI_domain_ok                      = c2_ok,
  c3_PFS_strict_prog_label_ok           = c3_ok,
  c4_preserved_nonNA_type               = c4_ok,
  n_locoreg_rules_applied               = length(c5_idx),
  c5_locoreg_nonempty                   = c5_ok_nonempty,
  c5_locoreg_no_tcga_left               = c5_no_tcga_left,
  c6_all_specific_rules_ok              = all((rule_expect %>% left_join(c6_checks, by=c("event","site_expected")))$ok),
  c7_exhep_nonempty                     = c7_ok_nonempty,
  c8_empty_site_other_is_NoEvent        = c8_ok_empty,
  c8_non_other_site_other_is_NoEvent    = c8_ok_not_other
)

notes <- c()
if (!c1_ok) notes <- c(notes, "Nem todos DFI==1 & NA(type) viraram 'Recurrence'.")
if (!c2_ok) notes <- c(notes, "PFI-family imputou fora do domínio permitido.")
if (!c3_ok) notes <- c(notes, "PFS imputou rótulo diferente de 'Progression of Disease' ou fora do alvo.")
if (!c4_ok) notes <- c(notes, "Houve alteração de tipo onde não havia NA (fora do alvo).")
if (!c5_ok_nonempty) notes <- c(notes, "Locorregional permaneceu sem site após a regra.")
if (!c5_no_tcga_left) notes <- c(notes, "Permaneceram códigos TCGA crus no site após conversão.")
if (!all((rule_expect %>% left_join(c6_checks, by=c("event","site_expected")))$ok))
  notes <- c(notes, "Alguma regra específica event→site não foi satisfeita.")
if (!c7_ok_nonempty) notes <- c(notes, "Extrahepatic Recurrence seguiu sem site após imputação de moda.")
if (!c8_ok_empty) notes <- c(notes, "site_other ≠ 'No event' quando site é vazio.")
if (!c8_ok_not_other) notes <- c(notes, "site_other ≠ 'No event' quando site não é 'Other, specify'.")

audit_notes <- tibble(note = if (length(notes)==0) "OK: no inconsistencies detected." else notes)

## Export reports
readr::write_tsv(rowlevel,     file.path(out_dir, "audit_Parte_13_events_rowlevel.tsv"))
readr::write_tsv(by_type,      file.path(out_dir, "audit_Parte_13_events_by_type.tsv"))
readr::write_tsv(audit_summary,file.path(out_dir, "audit_Parte_13_events_summary.tsv"))
readr::write_tsv(c6_checks,    file.path(out_dir, "audit_Parte_13_events_specific_rules.tsv"))
readr::write_tsv(audit_notes,  file.path(out_dir, "audit_Parte_13_events_notes.tsv"))

## Final sanity
stopifnot(nrow(df017) == nrow(df016))
cat("✓ Part 13 completed. Rows preserved: ", nrow(df017), "\n", sep = "")
cat("→ Row-level: ", normalizePath(file.path(out_dir, "audit_Parte_13_events_rowlevel.tsv")), "\n", sep = "")
cat("→ Summary:   ", normalizePath(file.path(out_dir, "audit_Parte_13_events_summary.tsv")),  "\n", sep = "")
cat("→ By-type:   ", normalizePath(file.path(out_dir, "audit_Parte_13_events_by_type.tsv")),  "\n", sep = "")
cat("→ Specific:  ", normalizePath(file.path(out_dir, "audit_Parte_13_events_specific_rules.tsv")),  "\n", sep = "")
cat("→ Notes:     ", normalizePath(file.path(out_dir, "audit_Parte_13_events_notes.tsv")),  "\n", sep = "")

rio::export(df017,"data_table_results/df017.tsv")


## ============================================================
## PART 14 — Diagnostics + Tolerant Harmonization + Audit
##   Target: ajcc_pathologic_tumor_stage, clinical_stage
##   Robust rules: NA | "" | (tolower(str_squish(.)) == "no data") → "Missing"
## ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(purrr); library(tibble); library(stringr)
})

out_dir <- "harmonization_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## 0) Preconditions
## -----------------------------
stopifnot(exists("df017"))
stage_vars <- c("ajcc_pathologic_tumor_stage", "clinical_stage")
faltam_17 <- setdiff(stage_vars, names(df017))
if (length(faltam_17) > 0) stop("Faltam colunas em df017: ", paste(faltam_17, collapse=", "))

## -----------------------------
## 1) Diagnostics of current AFTER (if df018 already exists)
## -----------------------------
if (exists("df018")) {
  diag_res <- tibble(column = stage_vars) %>%
    rowwise() %>%
    mutate(
      n_na            = sum(is.na(df018[[column]])),
      n_empty         = sum(df018[[column]] == "", na.rm = TRUE),
      n_nodata_exact  = sum(df018[[column]] == "No data", na.rm = TRUE),
      n_nodata_folded = sum(tolower(str_squish(df018[[column]])) == "no data", na.rm = TRUE),
      examples_variants =
        paste(
          head(
            unique(
              df018[[column]][
                !is.na(df018[[column]]) &
                  df018[[column]] != "" &
                  tolower(str_squish(df018[[column]])) == "no data" &
                  df018[[column]] != "No data"
              ]), 5),
          collapse=" | "
        )
    ) %>% ungroup()
  
  readr::write_tsv(diag_res, file.path(out_dir, "audit_Parte_14_stage_diagnostics_before_fix.tsv"))
}

## -----------------------------
## 2) Tolerant harmonization
## -----------------------------
harmonize_stage <- function(x) {
  # Convert to char, trim, and map variations of "No data"
  x_chr <- as.character(x)
  x_trim <- str_squish(x_chr)
  
  to_missing <-
    is.na(x_chr) |
    x_trim == "" |
    tolower(x_trim) == "no data"
  
  x_trim[to_missing] <- "Missing"
  x_trim
}

# Rebuild df018 from df017 with tolerant harmonization
df018 <- df017 %>%
  mutate(across(all_of(stage_vars), harmonize_stage))

stopifnot(nrow(df018) == nrow(df017))

## -----------------------------
## 3) Robust audit (join by .rowid; direct AFTER checks)
## -----------------------------

# Snapshots
df_before <- df017 %>%
  transmute(
    .rowid = dplyr::row_number(),
    type_before = as.character(type),
    DFI_before  = suppressWarnings(as.numeric(DFI)),
    PFI_before  = suppressWarnings(as.numeric(PFI)),
    across(all_of(stage_vars), ~ as.character(.), .names = "before_{col}")
  )

df_after  <- df018 %>%
  transmute(
    .rowid = dplyr::row_number(),
    type_after = as.character(type),
    DFI_after  = suppressWarnings(as.numeric(DFI)),
    PFI_after  = suppressWarnings(as.numeric(PFI)),
    across(all_of(stage_vars), ~ as.character(.), .names = "after_{col}")
  )

stopifnot(nrow(df_before) == nrow(df_after))

# Expected (from df017, with the same tolerant rule)
expected <- df017 %>%
  transmute(
    .rowid = dplyr::row_number(),
    across(all_of(stage_vars), harmonize_stage, .names = "exp_{col}")
  )

# Observed vs Expected comparison (by .rowid only)
cmp <- df_before %>%
  select(.rowid, starts_with("before_"),
         type_before, DFI_before, PFI_before) %>%
  left_join(df_after %>% select(.rowid, starts_with("after_"),
                                type_after, DFI_after, PFI_after),
            by = ".rowid") %>%
  left_join(expected, by = ".rowid")

mk_col_cmp <- function(var) {
  tibble(
    .rowid = cmp$.rowid,
    column = var,
    before   = cmp[[paste0("before_", var)]],
    after    = cmp[[paste0("after_",  var)]],
    expected = cmp[[paste0("exp_",    var)]]
  ) %>%
    mutate(
      changed = !( (before == after) | (is.na(before) & is.na(after)) ),
      match_expected = ( (after == expected) | (is.na(after) & is.na(expected)) )
    )
}
cmp_long <- bind_rows(lapply(stage_vars, mk_col_cmp))

# Direct AFTER check (no NA/""/"No data" should remain under tolerant rule)
direct_after <- tibble(column = stage_vars) %>%
  rowwise() %>%
  mutate(
    n_na     = sum(is.na(df018[[column]])),
    n_empty  = sum(df018[[column]] == "", na.rm = TRUE),
    n_nodata = sum(tolower(str_squish(df018[[column]])) == "no data", na.rm = TRUE)
  ) %>% ungroup()

if (any(direct_after$n_na > 0 | direct_after$n_empty > 0 | direct_after$n_nodata > 0)) {
  readr::write_tsv(direct_after, file.path(out_dir, "audit_Parte_14_stage_residuals_after.tsv"))
  stop("Falha de auditoria: ainda existem resíduos após harmonização tolerante. ",
       "Veja 'audit_Parte_14_stage_residuals_after.tsv'.")
}

# Obs vs Exp divergences (should be none)
divs <- cmp_long %>% filter(!match_expected)
if (nrow(divs) > 0) {
  readr::write_tsv(
    divs %>% transmute(.rowid, column,
                       before = replace_na(before, "<NA>"),
                       after  = replace_na(after,  "<NA>"),
                       expected = replace_na(expected, "<NA>")),
    file.path(out_dir, "audit_Parte_14_stage_divergences.tsv")
  )
  stop("Falha de auditoria: Observado ≠ Esperado em ", nrow(divs), " célula(s). ",
       "Detalhes em 'audit_Parte_14_stage_divergences.tsv'.")
}

# Informative DRIFT (type/DFI/PFI)
drift_summary <- cmp %>%
  transmute(
    type_changed = !(type_before == type_after | (is.na(type_before) & is.na(type_after))),
    DFI_changed  = !(DFI_before  == DFI_after  | (is.na(DFI_before) & is.na(DFI_after))),
    PFI_changed  = !(PFI_before  == PFI_after  | (is.na(PFI_before) & is.na(PFI_after)))
  ) %>%
  summarise(
    n_type_changed = sum(type_changed, na.rm = TRUE),
    n_DFI_changed  = sum(DFI_changed,  na.rm = TRUE),
    n_PFI_changed  = sum(PFI_changed,  na.rm = TRUE)
  )

## -----------------------------
## 4) Final reports
## -----------------------------
audit_rowlevel <- cmp_long %>%
  filter(changed) %>%
  transmute(rowid = .rowid, column, value_before = before, value_after = after) %>%
  arrange(column, rowid)

summ_col <- function(var) {
  v_before <- df017[[var]]
  v_after  <- df018[[var]]
  tibble(
    column = var,
    class_before = paste(class(v_before), collapse = "|"),
    class_after  = paste(class(v_after),  collapse = "|"),
    n_rows = length(v_before),
    n_before_NA      = sum(is.na(v_before)),
    n_before_empty   = sum(v_before == "", na.rm = TRUE),
    n_before_NoDataF = sum(tolower(str_squish(v_before)) == "no data", na.rm = TRUE),
    n_before_Missing = sum(v_before == "Missing", na.rm = TRUE),
    n_after_NA       = sum(is.na(v_after)),
    n_after_empty    = sum(v_after == "", na.rm = TRUE),
    n_after_NoDataF  = sum(tolower(str_squish(v_after)) == "no data", na.rm = TRUE),
    n_after_Missing  = sum(v_after == "Missing", na.rm = TRUE),
    n_changed_total  = sum(!( (v_before == v_after) | (is.na(v_before) & is.na(v_after)) ))
  )
}
audit_summary <- bind_rows(lapply(stage_vars, summ_col))

by_type <- purrr::map_dfr(stage_vars, function(var) {
  tibble(
    type = as.character(df018$type),
    value = as.character(df018[[var]])
  ) %>%
    mutate(column = var, is_missing = (value == "Missing")) %>%
    group_by(column, type) %>%
    summarise(
      n = n(),
      n_missing = sum(is_missing, na.rm = TRUE),
      frac_missing = ifelse(n > 0, n_missing / n, NA_real_),
      .groups = "drop"
    )
})

# Export
readr::write_tsv(diag_res %||% tibble(), file.path(out_dir, "audit_Parte_14_stage_diagnostics_before_fix.tsv"))
readr::write_tsv(audit_rowlevel, file.path(out_dir, "audit_Parte_14_stage_rowlevel.tsv"))
readr::write_tsv(audit_summary,  file.path(out_dir, "audit_Parte_14_stage_summary.tsv"))
readr::write_tsv(by_type,        file.path(out_dir, "audit_Parte_14_stage_by_type.tsv"))
readr::write_tsv(drift_summary,  file.path(out_dir, "audit_Parte_14_stage_drift_summary.tsv"))

cat("✓ Auditoria Parte 14 (tolerante) concluída.\n")
cat("→ Diagnostics (antes do fix): ", normalizePath(file.path(out_dir, "audit_Parte_14_stage_diagnostics_before_fix.tsv")), "\n", sep = "")
cat("→ Row-level:                 ", normalizePath(file.path(out_dir, "audit_Parte_14_stage_rowlevel.tsv")), "\n", sep = "")
cat("→ Summary:                   ", normalizePath(file.path(out_dir, "audit_Parte_14_stage_summary.tsv")),  "\n", sep = "")
cat("→ By-type:                   ", normalizePath(file.path(out_dir, "audit_Parte_14_stage_by_type.tsv")),  "\n", sep = "")
cat("→ Drift:                     ", normalizePath(file.path(out_dir, "audit_Parte_14_stage_drift_summary.tsv")),  "\n", sep = "")

rio::export(df018,"data_table_results/df018.tsv")


## ============================================================
## PART 15 — Harmonize birth_days_to a partir de age_at_initial_pathologic_diagnosis + AUDIT
##   Regra (avaliada no estado "before"):
##   R1) Se is.na(birth_days_to) & !is.na(age_at_initial_pathologic_diagnosis)
##       → birth_days_to := -365 * age_at_initial_pathologic_diagnosis
##   Obs.: Não altera valores não-NA já existentes em birth_days_to.
## ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble)
})

out_dir <- "harmonization_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
has_type_column <- "type" %in% names(df018)

## -----------------------------
## 0) Precondições
## -----------------------------
stopifnot(exists("df018"))
req <- c("birth_days_to", "age_at_initial_pathologic_diagnosis")
miss <- setdiff(req, names(df018))
if (length(miss) > 0) stop("Colunas ausentes em df018: ", paste(miss, collapse = ", "))

## Parâmetro opcional: usar 365.25 (anos médios) em vez de 365 exatos.
USE_365_25 <- FALSE
K_DAYS_PER_YEAR <- if (USE_365_25) 365.25 else 365

## -----------------------------
## 1) Snapshot "before" com tipos explícitos
## -----------------------------
df018_before <- df018 %>%
  mutate(
    birth0 = suppressWarnings(as.numeric(birth_days_to)),
    age0   = suppressWarnings(as.numeric(age_at_initial_pathologic_diagnosis)),
    .rowid = dplyr::row_number()
  )

b0 <- df018_before$birth0
a0 <- df018_before$age0

## -----------------------------
## 2) Flags (avaliadas no "before")
## -----------------------------
flag_R1 <- is.na(b0) & !is.na(a0)

## -----------------------------
## 3) EXPECTED — aplicação formal da regra
## -----------------------------
b_expected <- ifelse(flag_R1, -K_DAYS_PER_YEAR * a0, b0)

## -----------------------------
## 4) OBSERVED — aplicação via mutate/case_when (somente NA → imputação)
## -----------------------------
df019_after <- df018_before %>%
  mutate(
    birth_days_to = case_when(
      is.na(birth0) & !is.na(a0) ~ -K_DAYS_PER_YEAR * a0,
      TRUE                       ~ birth0
    )
  )

b_after <- df019_after$birth_days_to

## -----------------------------
## 5) Auditorias
## -----------------------------

# (a) Observed == Expected (comparação linha a linha com manejo de NA==NA)
mismatch_idx <- which( !(is.na(b_after) & is.na(b_expected)) & (b_after != b_expected) )
if (length(mismatch_idx) > 0) {
  stop("Audit failure: ", length(mismatch_idx),
       " linha(s) diferem entre observed e expected. Exemplos (rowid): ",
       paste(head(df018_before$.rowid[mismatch_idx], 10), collapse = ", "))
}

# (b) Imutabilidade de valores não-NA preexistentes
changed_non_na_idx <- which(!is.na(b0) & !is.na(b_after) & (b0 != b_after))
if (length(changed_non_na_idx) > 0) {
  stop("Audit failure: valores não-NA de 'birth_days_to' foram alterados em ",
       length(changed_non_na_idx), " linha(s). Exemplos (rowid): ",
       paste(head(df018_before$.rowid[changed_non_na_idx], 10), collapse = ", "))
}

# (c) Sanidade de domínio (aviso): datas de nascimento em dias devem ser ≤ 0
pos_idx <- which(!is.na(b_after) & b_after > 0)
if (length(pos_idx) > 0) {
  warning("Warning: 'birth_days_to' positivos detectados em ", length(pos_idx),
          " linha(s) após harmonização (espera-se ≤ 0). Verifique a origem/escala.")
}

# (d) Sanidade de idade (aviso): idades fora de faixa usual (0–120)
age_oob_idx <- which(!is.na(a0) & (a0 < 0 | a0 > 120))
if (length(age_oob_idx) > 0) {
  warning("Warning: valores de 'age_at_initial_pathologic_diagnosis' fora de 0–120 anos em ",
          length(age_oob_idx), " linha(s). Não alterados; verifique a fonte.")
}

## -----------------------------
## 6) Relatórios
## -----------------------------
imputed_idx <- which(is.na(b0) & !is.na(b_after))
audit_rowlevel <- tibble(
  rowid         = df018_before$.rowid[imputed_idx],
  type          = if (has_type_column) df018_before$type[imputed_idx] else NA_character_,
  birth_before  = b0[imputed_idx],
  birth_after   = b_after[imputed_idx],
  age_years     = a0[imputed_idx],
  rule_R1       = flag_R1[imputed_idx],
  factor_days_per_year = K_DAYS_PER_YEAR
)

audit_summary <- tibble(
  total_rows           = nrow(df018_before),
  n_birth_NA_before    = sum(is.na(b0)),
  n_age_nonNA          = sum(!is.na(a0)),
  n_imputed_R1         = sum(flag_R1),
  n_birth_NA_after     = sum(is.na(b_after)),
  check_balance        = (n_birth_NA_before - n_imputed_R1) == n_birth_NA_after,
  n_birth_positive     = length(pos_idx),
  n_age_out_of_bounds  = length(age_oob_idx),
  days_per_year_used   = K_DAYS_PER_YEAR
)

audit_by_type <- if (has_type_column) {
  tibble(type = df018_before$type, R1 = flag_R1) %>%
    group_by(type) %>%
    summarise(n_R1 = sum(R1), .groups = "drop") %>%
    arrange(desc(n_R1), type)
} else {
  tibble(info = "Coluna 'type' ausente; relatório por tipo não gerado.")
}

readr::write_tsv(audit_rowlevel, file.path(out_dir, "audit_Parte_15_birth_rowlevel.tsv"))
readr::write_tsv(audit_summary,  file.path(out_dir, "audit_Parte_15_birth_summary.tsv"))
readr::write_tsv(audit_by_type,  file.path(out_dir, "audit_Parte_15_birth_by_type.tsv"))

## -----------------------------
## 7) Atualizar objeto de trabalho e exportar
## -----------------------------
df019 <- df019_after %>%
  select(-.rowid, -birth0, -age0)

rio::export(df019, "data_table_results/df019.tsv")

cat("✓ Part 18→19 (birth_days_to) concluída.\n")
cat("→ Row-level: ", normalizePath(file.path(out_dir, "audit_Parte_15_birth_rowlevel.tsv")), "\n", sep = "")
cat("→ Summary:   ", normalizePath(file.path(out_dir, "audit_Parte_15_birth_summary.tsv")),  "\n", sep = "")
cat("→ By-type:   ", normalizePath(file.path(out_dir, "audit_Parte_15_birth_by_type.tsv")),  "\n", sep = "")


## ============================================================
## PART 16 — death_days_to := OS.time quando OS==1 & OS.time!=NA & death_days_to==NA
##            + AUDITORIA com razões do "por que não mudou"
## (versão sem 'df16_*' para não confundir com df006)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(tibble)
})

out_dir <- "harmonization_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(exists("df019"))
req <- c("OS", "OS.time", "death_days_to")
miss <- setdiff(req, names(df019))
if (length(miss) > 0) stop("PART 16: faltam colunas em df019: ", paste(miss, collapse = ", "))

has_type <- "type" %in% names(df019)

## 1) Snapshot BEFORE
p16_before <- df019 %>%
  mutate(
    OS        = suppressWarnings(as.numeric(OS)),
    OS_time   = suppressWarnings(as.numeric(OS.time)),
    death0    = suppressWarnings(as.numeric(death_days_to)),
    .rowid    = dplyr::row_number()
  )

## 2) Regra R1
flag16_R1 <- (p16_before$OS == 1) & !is.na(p16_before$OS_time) & is.na(p16_before$death0)

## 3) EXPECTED
death_expected16 <- ifelse(flag16_R1, p16_before$OS_time, p16_before$death0)

## 4) OBSERVED (aplica no dataset completo e gera df020)
df020 <- df019 %>%
  mutate(
    death_days_to = ifelse(flag16_R1, p16_before$OS_time, death_days_to)
  )

## 5) AFTER (apenas visão para auditoria; NÃO substitui df020)
p16_after_view <- df020 %>%
  mutate(death_after = suppressWarnings(as.numeric(death_days_to))) %>%
  transmute(
    .rowid        = p16_before$.rowid,
    type          = if (has_type) p16_before$type else NA_character_,
    OS            = p16_before$OS,
    OS.time       = p16_before$OS_time,
    death_before  = p16_before$death0,
    death_after   = death_after,
    changed_by_R1 = flag16_R1
  )

## Checagens
changed_obs <- which(!is.na(p16_before$death0) &
                       !is.na(p16_after_view$death_after) &
                       (p16_before$death0 != p16_after_view$death_after))
if (length(changed_obs) > 0) {
  bad <- head(p16_after_view[changed_obs, c(".rowid","type","death_before","death_after")], 10)
  stop("Audit failure (Part 16): valores não-NA de 'death_days_to' foram alterados.\n",
       "Exemplos (rowid | type | before | after):\n",
       paste(apply(bad, 1, paste, collapse = " | "), collapse = "\n"))
}

neq <- which(!(is.na(p16_after_view$death_after) & is.na(death_expected16)) &
               (p16_after_view$death_after != death_expected16))
if (length(neq) > 0) {
  bad <- head(p16_after_view[neq, c(".rowid","type","death_before","death_after")], 10)
  stop("Audit failure (Part 16): observado != esperado em ", length(neq), " linha(s).\n",
       "Exemplos (rowid | type | before | after):\n",
       paste(apply(bad, 1, paste, collapse = " | "), collapse = "\n"))
}

## 6) Razões
reason16 <- dplyr::case_when(
  flag16_R1 ~ "Changed by R1: OS==1 & OS.time present & death_days_to was NA",
  !is.na(p16_before$death0) ~ "Unchanged: death_days_to already filled",
  is.na(p16_before$OS) ~ "Unchanged: OS is NA",
  p16_before$OS != 1 ~ "Unchanged: OS != 1",
  is.na(p16_before$OS_time) ~ "Unchanged: OS.time is NA",
  TRUE ~ "Unchanged: not eligible (other)"
)

audit_rowlevel <- tibble(
  rowid  = p16_before$.rowid,
  type   = if (has_type) p16_before$type else NA_character_,
  OS     = p16_before$OS,
  OS.time= p16_before$OS_time,
  death_days_to_before = p16_before$death0,
  death_days_to_after  = p16_after_view$death_after,
  changed_by_R1 = flag16_R1,
  reason = reason16
)

audit_summary <- tibble(
  total_rows       = nrow(p16_before),
  n_changed_R1     = sum(flag16_R1, na.rm = TRUE),
  n_already_filled = sum(!is.na(p16_before$death0)),
  n_OS_NA          = sum(is.na(p16_before$OS)),
  n_OS_not_1       = sum(!is.na(p16_before$OS) & p16_before$OS != 1),
  n_OS_time_NA     = sum(is.na(p16_before$OS_time))
)

audit_by_type <- if (has_type) {
  audit_rowlevel %>%
    dplyr::mutate(changed = ifelse(changed_by_R1, "changed_by_R1", "unchanged")) %>%
    dplyr::count(type, changed, reason, name = "n") %>%
    dplyr::arrange(type, dplyr::desc(changed), dplyr::desc(n))
} else {
  tibble(info = "Column 'type' not available; by-type report omitted.")
}

## 7) Exports
readr::write_tsv(audit_rowlevel, file.path(out_dir, "audit_Parte_16_deathFromOS_rowlevel.tsv"))
readr::write_tsv(audit_summary,  file.path(out_dir, "audit_Parte_16_deathFromOS_summary.tsv"))
if (has_type) readr::write_tsv(audit_by_type, file.path(out_dir, "audit_Parte_16_deathFromOS_by_type.tsv"))

cat("✓ Part 16 concluída.\n")
cat("→ Row-level: ", normalizePath(file.path(out_dir, "audit_Parte_16_deathFromOS_rowlevel.tsv")), "\n", sep = "")
cat("→ Summary:   ", normalizePath(file.path(out_dir, "audit_Parte_16_deathFromOS_summary.tsv")),  "\n", sep = "")
if (has_type) {
  cat("→ By-type:   ", normalizePath(file.path(out_dir, "audit_Parte_16_deathFromOS_by_type.tsv")),  "\n", sep = "")
}

rio::export(df020, "data_table_results/df020.tsv")
saveRDS(df020, "data_table_results/df020.rds")
