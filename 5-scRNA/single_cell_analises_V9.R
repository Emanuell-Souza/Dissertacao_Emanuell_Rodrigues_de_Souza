## ======================================================================
## Pipeline scRNA-seq — versão otimizada
## Integração RPCA + marcadores + assinaturas + anotação automática
## Suporte a genes e transcritos (ENST...)
## FILTRO ATIVO: apenas assinaturas de camada ômica 6
##   (2º token da nomenclatura após split por '.')
##   Ex: "BRCA-693.6.3.N..." → camada 6 ✓
##       "BRCA-693.5.3.N..." → camada 5 ✗
## ======================================================================
## MUDANÇAS PRINCIPAIS vs versão original:
##   1. map_queries_to_features   → lookup vetorizado (sem loop for em R)
##   2. annotate_clusters         → opera sobre matriz esparsa; colSums
##                                  por cluster sem coerção para denso
##   3. feature_map               → construído 1x por câncer, passado por
##                                  referência via environment
##   4. score_signature           → sem regex de limpeza desnecessária
##   5. save_marker_exploration   → plots combinados; sem loop de ggsave
##                                  por gene (salva 1 PDF multipágina)
##   6. checkpoint                → NÃO guarda integrated_obj no RDS;
##                                  guarda apenas metadados e sumários
##   7. future                    → DESATIVADO globalmente; o Seurat usa
##                                  objetos S4 grandes que estouram o
##                                  limite de 500 MiB dos workers ao ser
##                                  serializado. Rodar em sequential é
##                                  obrigatório para evitar o erro.
##   8. JoinLayers                → chamado apenas quando necessário
##   9. Rm() explícito            → libera memória após cada câncer
##  10. CORREÇÃO: .log_cancer     → função declarada e retorno atribuído
##  11. CORREÇÃO CRÍTICA: conflito de namespace dplyr::slice vs
##      Seurat/Matrix::slice → todas as funções dplyr agora são
##      chamadas com o prefixo dplyr:: explicitamente nas funções
##      que operam sobre grouped_df para evitar o erro:
##      "método não aplicável para 'slice' aplicado a um objeto de
##      classe grouped_df"
##  12. NOVO: get_omic_layer      → extrai o 2º token da nomenclatura
##      e filtra signatures_df mantendo apenas camada ômica == "6"
## ======================================================================

set.seed(123)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(rio)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(openxlsx)
  library(readr)
})

## -----------------------------------------------------------------------
## Paralelismo — DESATIVADO intencionalmente.
## -----------------------------------------------------------------------
.use_future <- FALSE

if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
  message("Plano future definido como 'sequential'.")
} else {
  message("Pacote {future} não instalado. Rodando em modo serial.")
}

## -----------------------------------------------------------------------
## 1) Diretórios
## -----------------------------------------------------------------------
setwd("~/Emanuell/scRNA_seq")

base_dir     <- "~/Emanuell/scRNA_seq"
samples_root <- file.path(base_dir, "amostras")
results_root <- file.path(base_dir, "resultados")

max_samples_per_cancer <- 3   # altere para Inf para usar todas

ensure_directory <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

ensure_directory(results_root)

## -----------------------------------------------------------------------
## 2) Assinaturas — FILTRO: apenas camada ômica 6
## -----------------------------------------------------------------------

## Extrai a camada ômica da nomenclatura (2º token após split por '.')
## Ex: "BRCA-693.6.3.N.2.95.95.1.1.1" → "6"
##     "BRCA-693.5.3.N.2.95.95.1.1.1" → "5"
get_omic_layer <- function(nomenclatura) {
  sapply(nomenclatura, function(x) {
    parts <- strsplit(x, "\\.")[[1]]
    if (length(parts) >= 2) parts[2] else NA_character_
  }, USE.NAMES = FALSE)
}

signature_file <- file.path(base_dir, "Assinaturas_Omicas_SuperLearner.tsv")

if (!file.exists(signature_file)) stop("Arquivo não encontrado: ", signature_file)

signatures_df <- readr::read_tsv(
  signature_file,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
) %>%
  dplyr::mutate(dplyr::across(c(Nomenclature, Signature, Tipo),
                              ~ trimws(as.character(.)))) %>%
  dplyr::filter(
    !is.na(Nomenclature), Nomenclature != "",
    !is.na(Signature),    Signature    != "",
    !is.na(Tipo),         Tipo         != ""
  ) %>%
  dplyr::rename(
    nomeclatura = Nomenclature,
    assinaturas = Signature,
    type        = Tipo
  ) %>%
  ## ── FILTRO: manter apenas assinaturas de camada ômica 6 ─────────────
  dplyr::filter(get_omic_layer(nomeclatura) == "6")

message("Assinaturas carregadas (camada ômica 6): ", nrow(signatures_df))
message("Tipos de câncer presentes: ",
        paste(unique(signatures_df$type), collapse = ", "))

if (nrow(signatures_df) == 0)
  stop("Nenhuma assinatura de camada ômica 6 encontrada no arquivo.")

safe_export <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  
  tryCatch(
    {
      openxlsx::write.xlsx(x, file = file, overwrite = TRUE)
    },
    error = function(e) {
      warning("Falha ao salvar XLSX: ", file, " | ", e$message)
      
      csv_file <- sub("\\.xlsx$", ".csv", file)
      readr::write_csv(x, csv_file)
      
      message("Arquivo salvo como CSV alternativo: ", csv_file)
    }
  )
}

## -----------------------------------------------------------------------
## 3) Genes marcadores
## -----------------------------------------------------------------------
marker_file <- file.path(base_dir, "Genes_markers.xlsx")

if (!file.exists(marker_file)) stop("Arquivo não encontrado: ", marker_file)

markers_ref <- rio::import(marker_file) %>%
  dplyr::mutate(dplyr::across(c(Gene, Cell_type, Type),
                              ~ trimws(as.character(.)))) %>%
  dplyr::filter(
    !is.na(Gene), Gene != "",
    !is.na(Cell_type), Cell_type != "",
    !is.na(Type), Type != ""
  ) %>%
  dplyr::distinct()

message("Marcadores carregados: ", nrow(markers_ref))

## -----------------------------------------------------------------------
## 4) Genes para exploração visual
## -----------------------------------------------------------------------
genes_to_test <- c(
  "GCSAML",
  "EPCAM", "KRT8", "KRT18", "KRT19",
  "CD3D", "CD3E", "CD4", "CD8A", "TRBC1",
  "MS4A1", "CD79A",
  "NKG7", "GNLY",
  "LYZ", "S100A8", "S100A9", "LST1",
  "CD68", "SPP1", "C1QC",
  "COL1A1", "COL1A2", "DCN", "LUM",
  "ACTA2", "TAGLN", "MYLK",
  "PECAM1", "VWF", "KDR",
  "MKI67", "TOP2A",
  "HBB", "HBA1"
)

## -----------------------------------------------------------------------
## 5) Funções utilitárias
## -----------------------------------------------------------------------
parse_signature_features <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(character(0))
  x <- gsub("[()]", "", trimws(x))
  parts <- trimws(unlist(strsplit(x, "\\+")))
  unique(parts[nzchar(parts)])
}

sanitize_filename <- function(x) {
  x <- trimws(x)
  x <- gsub("[/\\\\:*?\"<>|]", "_", x)
  gsub("\\s+", "_", x)
}

is_transcript_id <- function(x) grepl("^ENST", x, ignore.case = TRUE)

classify_signature_input <- function(x) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return("empty")
  p <- mean(is_transcript_id(x))
  if (p == 1) "transcript_only" else if (p > 0) "mixed" else "gene_only"
}

format_elapsed_time <- function(seconds) {
  h <- seconds %/% 3600
  m <- (seconds %% 3600) %/% 60
  s <- round(seconds %% 60, 2)
  sprintf("%02d:%02d:%05.2f", h, m, s)
}

## -----------------------------------------------------------------------
## 5A) .log_cancer
## -----------------------------------------------------------------------
.log_cancer <- function(cp, session_id, cancer_name, status, t0) {
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cp$cancer_timing_log <- dplyr::bind_rows(
    cp$cancer_timing_log,
    data.frame(
      session_id        = session_id,
      cancer_type       = cancer_name,
      status            = status,
      start_time        = format(t0, "%Y-%m-%d %H:%M:%S"),
      end_time          = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      elapsed_seconds   = elapsed,
      elapsed_formatted = format_elapsed_time(elapsed),
      stringsAsFactors  = FALSE
    )
  )
  cp
}

## -----------------------------------------------------------------------
## 5B) feature_map
## -----------------------------------------------------------------------
build_feature_map <- function(obj, assay_expr = "RNA") {
  fn <- rownames(obj[[assay_expr]])
  fm <- data.frame(
    seurat_feature = fn,
    feature_id     = fn,
    feature_name   = fn,
    feature_type   = ifelse(grepl("^ENST", fn, ignore.case = TRUE),
                            "transcript", "gene_or_other"),
    stringsAsFactors = FALSE
  )
  rownames(fm) <- fn
  fm
}

get_feature_map <- function(obj, assay_expr = "RNA") {
  fm <- obj@misc$feature_map
  if (!is.null(fm) && nrow(fm) > 0) {
    if (!"feature_type" %in% colnames(fm)) {
      fm$feature_type <- ifelse(grepl("^ENST", fm$feature_id, ignore.case = TRUE),
                                "transcript", "gene_or_other")
    }
    rownames(fm) <- fm$seurat_feature
    return(fm)
  }
  build_feature_map(obj, assay_expr)
}

## -----------------------------------------------------------------------
## 5C) map_queries_to_features — totalmente vetorizado
## -----------------------------------------------------------------------
map_queries_to_features <- function(obj, query_vec, assay_expr = "RNA",
                                    feature_map = NULL) {
  if (is.null(feature_map)) feature_map <- get_feature_map(obj, assay_expr)
  
  assay_features <- rownames(obj[[assay_expr]])
  qv <- unique(trimws(as.character(query_vec)))
  qv <- qv[!is.na(qv) & qv != ""]
  
  if (length(qv) == 0) {
    empty_report <- data.frame(
      query = character(0), matched_feature = character(0),
      matched_feature_id = character(0), matched_feature_name = character(0),
      matched_feature_type = character(0), match_source = character(0),
      stringsAsFactors = FALSE
    )
    return(list(features = character(0), mapping_report = empty_report))
  }
  
  n   <- length(qv)
  hit <- rep(NA_character_, n)
  src <- rep(NA_character_, n)
  names(hit) <- qv
  
  fm_sf   <- feature_map$seurat_feature
  fm_id   <- feature_map$feature_id
  fm_nm   <- feature_map$feature_name
  fm_sf_u <- toupper(fm_sf)
  fm_id_u <- toupper(fm_id)
  fm_nm_u <- toupper(fm_nm)
  af_u    <- toupper(assay_features)
  qv_u    <- toupper(qv)
  
  .apply_match <- function(hit, src, needle, haystack, sf_vec, label) {
    unresolved <- which(is.na(hit))
    if (length(unresolved) == 0) return(list(hit = hit, src = src))
    idx <- match(needle[unresolved], haystack)
    found <- !is.na(idx)
    pos   <- unresolved[found]
    hit[pos] <- sf_vec[idx[found]]
    src[pos] <- label
    list(hit = hit, src = src)
  }
  
  r <- list(hit = hit, src = src)
  r <- .apply_match(r$hit, r$src, qv,   fm_sf,          fm_sf,          "seurat_feature_exact")
  r <- .apply_match(r$hit, r$src, qv,   fm_id,          fm_sf,          "feature_id_exact")
  r <- .apply_match(r$hit, r$src, qv,   fm_nm,          fm_sf,          "feature_name_exact")
  r <- .apply_match(r$hit, r$src, qv_u, fm_sf_u,        fm_sf,          "seurat_feature_ci")
  r <- .apply_match(r$hit, r$src, qv_u, fm_id_u,        fm_sf,          "feature_id_ci")
  r <- .apply_match(r$hit, r$src, qv_u, fm_nm_u,        fm_sf,          "feature_name_ci")
  r <- .apply_match(r$hit, r$src, qv,   assay_features, assay_features, "assay_exact")
  r <- .apply_match(r$hit, r$src, qv_u, af_u,           assay_features, "assay_ci")
  
  hit <- r$hit
  src <- r$src
  
  resolved_sf   <- hit[!is.na(hit)]
  resolved_qv   <- qv[!is.na(hit)]
  unresolved_qv <- qv[is.na(hit)]
  
  if (length(resolved_sf) > 0) {
    fm_idx <- match(resolved_sf, feature_map$seurat_feature)
    fm_ok  <- !is.na(fm_idx)
    
    rep_found <- data.frame(
      query                = resolved_qv,
      matched_feature      = resolved_sf,
      matched_feature_id   = ifelse(fm_ok, feature_map$feature_id  [fm_idx], resolved_sf),
      matched_feature_name = ifelse(fm_ok, feature_map$feature_name[fm_idx], resolved_sf),
      matched_feature_type = ifelse(fm_ok, feature_map$feature_type[fm_idx],
                                    ifelse(grepl("^ENST", resolved_sf, TRUE),
                                           "transcript", "gene_or_other")),
      match_source         = src[!is.na(hit)],
      stringsAsFactors     = FALSE
    )
  } else {
    rep_found <- data.frame(
      query = character(0), matched_feature = character(0),
      matched_feature_id = character(0), matched_feature_name = character(0),
      matched_feature_type = character(0), match_source = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  rep_miss <- if (length(unresolved_qv) > 0) {
    data.frame(
      query                = unresolved_qv,
      matched_feature      = NA_character_,
      matched_feature_id   = NA_character_,
      matched_feature_name = NA_character_,
      matched_feature_type = NA_character_,
      match_source         = "not_found",
      stringsAsFactors     = FALSE
    )
  } else NULL
  
  mapping_report <- dplyr::bind_rows(rep_found, rep_miss)
  
  list(features = unique(resolved_sf), mapping_report = mapping_report)
}

## -----------------------------------------------------------------------
## 5D) Checkpoint leve e robusto
## -----------------------------------------------------------------------
initialize_checkpoint <- function(results_root) {
  d <- file.path(results_root, "_checkpoint")
  ensure_directory(d)
  
  list(
    checkpoint_dir     = d,
    checkpoint_file    = file.path(d, "scRNA_seq_checkpoint_LIGHT.rds"),
    session_log_file   = file.path(d, "session_history.xlsx"),
    cancer_timing_file = file.path(d, "cancer_timing_history.xlsx")
  )
}

load_or_create_checkpoint <- function(checkpoint_file) {
  empty_cp <- list(
    completed_cancers    = character(0),
    completed_signatures = list(),
    session_history      = data.frame(),
    cancer_timing_log    = data.frame(),
    last_finished_cancer = NA_character_
  )
  
  if (!file.exists(checkpoint_file)) {
    return(list(checkpoint = empty_cp, resumed = FALSE))
  }
  
  cp <- tryCatch(
    readRDS(checkpoint_file),
    error = function(e) {
      warning(
        "Checkpoint antigo/corrompido ou grande demais não foi carregado: ",
        e$message,
        "\nCriando checkpoint leve novo."
      )
      empty_cp
    }
  )
  
  for (fld in names(empty_cp)) {
    if (is.null(cp[[fld]])) cp[[fld]] <- empty_cp[[fld]]
  }
  
  ## Segurança: nunca carregar/salvar resultados pesados no checkpoint
  cp$all_results <- NULL
  
  list(checkpoint = cp, resumed = TRUE)
}

save_checkpoint_files <- function(cp, paths) {
  cp_save <- cp
  
  ## Segurança absoluta: remove qualquer objeto pesado acidental
  cp_save$all_results <- NULL
  
  tmp_file <- paste0(paths$checkpoint_file, ".tmp")
  saveRDS(cp_save, tmp_file, compress = TRUE)
  file.rename(tmp_file, paths$checkpoint_file)
  
  if (nrow(cp_save$session_history) > 0) {
    safe_export(cp_save$session_history, paths$session_log_file)
  }
  
  if (nrow(cp_save$cancer_timing_log) > 0) {
    safe_export(cp_save$cancer_timing_log, paths$cancer_timing_file)
  }
  
  invisible(TRUE)
}

## -----------------------------------------------------------------------
## 6) score_signature
## -----------------------------------------------------------------------
score_signature <- function(obj, up, down = NULL, sig_name = "Sig",
                            assay_expr = "RNA") {
  stopifnot(assay_expr %in% names(obj@assays))
  DefaultAssay(obj) <- assay_expr
  
  fm <- get_feature_map(obj, assay_expr)
  
  up_res <- map_queries_to_features(obj, up,
                                    assay_expr = assay_expr, feature_map = fm)
  dn_res <- map_queries_to_features(obj,
                                    if (is.null(down)) character(0) else down,
                                    assay_expr = assay_expr, feature_map = fm)
  
  genes_up <- up_res$features
  genes_dn <- dn_res$features
  
  message(sprintf("  [%s] UP: %d/%d encontrados  |  DOWN: %d/%d encontrados",
                  sig_name,
                  length(genes_up), length(up),
                  length(genes_dn),
                  length(if (is.null(down)) character(0) else down)))
  
  if (length(genes_up) == 0 && length(genes_dn) == 0)
    stop(sprintf("Nenhum item da assinatura '%s' encontrado no assay %s.",
                 sig_name, assay_expr))
  
  old <- grep(paste0("^", sig_name, "_(MS\\d+|UP_Score|DOWN_Score|Composite|UCell)$"),
              colnames(obj@meta.data), value = TRUE)
  if (length(old) > 0) obj@meta.data[, old] <- NULL
  
  feature_list <- list()
  if (length(genes_up) > 0) feature_list$UP   <- genes_up
  if (length(genes_dn) > 0) feature_list$DOWN <- genes_dn
  
  obj <- AddModuleScore(obj, features = feature_list,
                        name = paste0(sig_name, "_MS"), assay = assay_expr)
  
  n_fl     <- length(feature_list)
  new_cols <- tail(colnames(obj@meta.data), n_fl)
  colnames(obj@meta.data)[match(new_cols, colnames(obj@meta.data))] <-
    paste0(sig_name, "_", names(feature_list), "_Score")
  
  has_up   <- paste0(sig_name, "_UP_Score")   %in% colnames(obj@meta.data)
  has_down <- paste0(sig_name, "_DOWN_Score") %in% colnames(obj@meta.data)
  
  obj[[paste0(sig_name, "_Composite")]] <- {
    if (has_up && has_down) {
      obj[[paste0(sig_name, "_UP_Score")]] - obj[[paste0(sig_name, "_DOWN_Score")]]
    } else if (has_up) {
      obj[[paste0(sig_name, "_UP_Score")]]
    } else {
      -obj[[paste0(sig_name, "_DOWN_Score")]]
    }
  }
  
  if (requireNamespace("UCell", quietly = TRUE) && length(genes_up) > 0) {
    obj <- UCell::AddModuleScore_UCell(
      obj, features = list(tmp = genes_up),
      name = paste0(sig_name, "_UCell"), assay = assay_expr
    )
  }
  
  obj@misc[[paste0(sig_name, "_mapping_report_up")]]   <- up_res$mapping_report
  obj@misc[[paste0(sig_name, "_mapping_report_down")]] <- dn_res$mapping_report
  
  DefaultAssay(obj) <- if ("integrated" %in% names(obj@assays)) "integrated" else assay_expr
  obj
}

## -----------------------------------------------------------------------
## 7) annotate_clusters_by_marker_union — CORRIGIDO: dplyr:: explícito
##    para evitar conflito com Seurat/Matrix::slice em grouped_df
## -----------------------------------------------------------------------
annotate_clusters_by_marker_union <- function(
    integrated_obj,
    cancer_name,
    markers_ref,
    assay_expr       = "RNA",
    slot_expr        = "data",
    expr_threshold   = 0,
    union_threshold  = 0.50,
    min_detected_genes = 2
) {
  DefaultAssay(integrated_obj) <- assay_expr
  
  if (!"seurat_clusters" %in% colnames(integrated_obj@meta.data))
    stop("Coluna 'seurat_clusters' não encontrada.")
  
  ## CORREÇÃO: usar dplyr::filter explicitamente para evitar máscara do Seurat
  cancer_markers <- dplyr::filter(markers_ref, Type == cancer_name)
  cluster_ids    <- sort(unique(as.character(integrated_obj$seurat_clusters)))
  
  empty_return <- function(reason) {
    list(
      integrated_obj        = integrated_obj,
      coverage_table        = data.frame(),
      marker_mapping_report = data.frame(),
      panel_gene_status     = data.frame(),
      cluster_assignments   = data.frame(),
      audit_unassigned      = data.frame(
        cancer_type        = cancer_name,
        cluster            = cluster_ids,
        assigned_cell_type = NA_character_,
        reason             = reason,
        stringsAsFactors   = FALSE
      ),
      audit_conflicts = data.frame()
    )
  }
  
  if (nrow(cancer_markers) == 0) {
    warning("Nenhum marcador encontrado para o tipo ", cancer_name)
    return(empty_return("no_marker_reference_for_this_cancer_type"))
  }
  
  expr_mat <- tryCatch({
    if (slot_expr == "data")
      LayerData(integrated_obj, assay = assay_expr, layer = "data")
    else if (slot_expr == "counts")
      LayerData(integrated_obj, assay = assay_expr, layer = "counts")
    else
      GetAssayData(integrated_obj, assay = assay_expr, slot = slot_expr)
  }, error = function(e)
    stop("Falha ao obter matriz de expressão: ", e$message))
  
  cluster_vec <- setNames(
    as.character(integrated_obj$seurat_clusters),
    colnames(integrated_obj)
  )
  
  fm         <- get_feature_map(integrated_obj, assay_expr)
  cell_types <- unique(cancer_markers$Cell_type)
  
  coverage_rows <- vector("list", length(cell_types) * length(cluster_ids))
  mapping_rows  <- vector("list", length(cell_types))
  panel_rows    <- vector("list", length(cell_types))
  k <- 1L
  
  for (ct in cell_types) {
    ## CORREÇÃO: dplyr::filter em vez de filter (possível conflito com stats::filter)
    panel_df       <- dplyr::filter(cancer_markers, Cell_type == ct)
    marker_queries <- unique(panel_df$Gene)
    
    resolved <- map_queries_to_features(
      obj        = integrated_obj,
      query_vec  = marker_queries,
      assay_expr = assay_expr,
      feature_map = fm
    )
    resolved_features <- resolved$features
    
    rep <- resolved$mapping_report
    if (nrow(rep) > 0) {
      rep$cancer_type <- cancer_name
      rep$Cell_type   <- ct
    }
    mapping_rows[[ct]] <- rep
    
    panel_rows[[ct]] <- data.frame(
      cancer_type     = cancer_name,
      Cell_type       = ct,
      query_gene      = marker_queries,
      found_in_object = marker_queries %in%
        rep$query[!is.na(rep$matched_feature)],
      stringsAsFactors = FALSE
    )
    
    rf_use   <- intersect(resolved_features, rownames(expr_mat))
    sub_expr <- if (length(rf_use) > 0) expr_mat[rf_use, , drop = FALSE] else NULL
    
    for (cl in cluster_ids) {
      cluster_cells   <- names(cluster_vec)[cluster_vec == cl]
      n_cells_cluster <- length(cluster_cells)
      
      row_base <- data.frame(
        cancer_type                  = cancer_name,
        cluster                      = cl,
        Cell_type                    = ct,
        n_cells_cluster              = n_cells_cluster,
        n_panel_genes_input          = length(marker_queries),
        n_panel_genes_found          = 0L,
        n_detected_genes_in_cluster  = 0L,
        n_cells_union_positive       = 0L,
        union_coverage               = 0,
        pass_union_threshold         = FALSE,
        stringsAsFactors             = FALSE
      )
      
      if (n_cells_cluster == 0 || is.null(sub_expr) || length(rf_use) == 0) {
        coverage_rows[[k]] <- row_base
        k <- k + 1L
        next
      }
      
      cl_cells <- intersect(cluster_cells, colnames(sub_expr))
      if (length(cl_cells) == 0) {
        coverage_rows[[k]] <- row_base
        k <- k + 1L
        next
      }
      
      cl_sub        <- sub_expr[, cl_cells, drop = FALSE]
      gene_row_sums <- Matrix::rowSums(cl_sub > expr_threshold)
      n_detected    <- sum(gene_row_sums > 0)
      cell_col_sums <- Matrix::colSums(cl_sub > expr_threshold)
      n_union_pos   <- sum(cell_col_sums > 0)
      union_cov     <- n_union_pos / n_cells_cluster
      
      pass <- (union_cov >= union_threshold) &&
        (n_detected >= min_detected_genes)
      
      row_base$n_panel_genes_found         <- length(rf_use)
      row_base$n_detected_genes_in_cluster <- n_detected
      row_base$n_cells_union_positive      <- n_union_pos
      row_base$union_coverage              <- union_cov
      row_base$pass_union_threshold        <- pass
      
      coverage_rows[[k]] <- row_base
      k <- k + 1L
    }
  }
  
  coverage_table        <- dplyr::bind_rows(coverage_rows[seq_len(k - 1L)])
  marker_mapping_report <- dplyr::bind_rows(mapping_rows)
  panel_gene_status     <- dplyr::bind_rows(panel_rows)
  
  if (nrow(coverage_table) == 0) return(empty_return("coverage_table_empty"))
  
  ## CORREÇÃO PRINCIPAL: dplyr::filter e dplyr::slice(1) explícitos
  ## O Seurat exporta uma função slice() que mascara dplyr::slice()
  ## em objetos grouped_df, causando o erro reportado.
  pass_table <- dplyr::filter(coverage_table, pass_union_threshold)
  
  if (nrow(pass_table) == 0) {
    integrated_obj$cell_type_auto <- "Unassigned"
    res <- empty_return("no_cell_type_reached_union_threshold")
    res$integrated_obj        <- integrated_obj
    res$coverage_table        <- coverage_table
    res$marker_mapping_report <- marker_mapping_report
    res$panel_gene_status     <- panel_gene_status
    return(res)
  }
  
  n_passed_by_cluster <- pass_table %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(n_passed_cell_types = dplyr::n(), .groups = "drop")
  
  ## CORREÇÃO: dplyr::slice(1) em vez de slice(1)
  best_assignment <- pass_table %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(
      dplyr::desc(union_coverage),
      dplyr::desc(n_detected_genes_in_cluster),
      dplyr::desc(n_panel_genes_found),
      .by_group = TRUE
    ) %>%
    dplyr::slice(1L) %>%      ## ← CORREÇÃO CRÍTICA
    dplyr::ungroup() %>%
    dplyr::select(
      cluster,
      assigned_cell_type = Cell_type,
      union_coverage,
      n_detected_genes_in_cluster,
      n_panel_genes_found
    )
  
  cluster_assignments <- data.frame(cluster = cluster_ids,
                                    stringsAsFactors = FALSE) %>%
    dplyr::left_join(best_assignment,    by = "cluster") %>%
    dplyr::left_join(n_passed_by_cluster, by = "cluster") %>%
    dplyr::mutate(
      cancer_type         = cancer_name,
      n_passed_cell_types = dplyr::coalesce(n_passed_cell_types, 0L)
    ) %>%
    dplyr::select(
      cancer_type, cluster, assigned_cell_type, union_coverage,
      n_detected_genes_in_cluster, n_panel_genes_found, n_passed_cell_types
    )
  
  audit_unassigned <- dplyr::filter(
    cluster_assignments,
    is.na(assigned_cell_type) | assigned_cell_type == ""
  ) %>%
    dplyr::mutate(reason = "no_cell_type_reached_union_threshold")
  
  ## CORREÇÃO: dplyr::summarise com dplyr:: explícito
  audit_conflicts <- dplyr::filter(cluster_assignments, n_passed_cell_types > 1) %>%
    dplyr::left_join(
      pass_table %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          passed_cell_types = paste(Cell_type,           collapse = " | "),
          passed_coverages  = paste(round(union_coverage, 4), collapse = " | "),
          .groups = "drop"
        ),
      by = "cluster"
    )
  
  assign_vec <- setNames(cluster_assignments$assigned_cell_type,
                         cluster_assignments$cluster)
  ct_auto    <- assign_vec[as.character(integrated_obj$seurat_clusters)]
  ct_auto[is.na(ct_auto)] <- "Unassigned"
  integrated_obj$cell_type_auto <- unname(ct_auto)
  
  list(
    integrated_obj        = integrated_obj,
    coverage_table        = coverage_table,
    marker_mapping_report = marker_mapping_report,
    panel_gene_status     = panel_gene_status,
    cluster_assignments   = cluster_assignments,
    audit_unassigned      = audit_unassigned,
    audit_conflicts       = audit_conflicts
  )
}

## -----------------------------------------------------------------------
## 8) Exploração de marcadores — 1 PDF multipágina
## -----------------------------------------------------------------------
save_marker_exploration <- function(integrated_obj, outdir, cancer_name,
                                    genes_to_test) {
  message("  >>> Explorando genes marcadores customizados")
  DefaultAssay(integrated_obj) <- "RNA"
  
  mapped        <- map_queries_to_features(integrated_obj, genes_to_test, "RNA")
  genes_present <- mapped$features
  
  gene_report <- mapped$mapping_report %>%
    dplyr::mutate(
      cancer_type        = cancer_name,
      present_in_dataset = !is.na(matched_feature)
    )
  safe_export(gene_report, file.path(outdir, paste0(cancer_name, "_tested_genes_presence.xlsx")))
  
  if (length(genes_present) == 0) {
    message("  Nenhum gene de interesse encontrado.")
    return(invisible(NULL))
  }
  
  p_fp <- FeaturePlot(integrated_obj, features = genes_present, min.cutoff = 0,
                      ncol = min(4L, length(genes_present)))
  ggsave(file.path(outdir, paste0(cancer_name, "_FeaturePlot_markers.pdf")),
         p_fp,
         width  = max(10, 3 * min(4L, length(genes_present))),
         height = max(8,  3 * ceiling(length(genes_present) / min(4L, length(genes_present)))))
  
  p_dot <- DotPlot(integrated_obj, features = genes_present,
                   group.by = "seurat_clusters") + RotatedAxis()
  ggsave(file.path(outdir, paste0(cancer_name, "_DotPlot_markers.pdf")),
         p_dot,
         width  = max(10, 0.35 * length(genes_present) + 6),
         height = 6)
  
  pdf(file.path(outdir, paste0(cancer_name, "_VlnPlots_markers.pdf")),
      width = 9, height = 6)
  for (g in genes_present) {
    print(VlnPlot(integrated_obj, features = g,
                  group.by = "seurat_clusters", pt.size = 0) +
            ggtitle(paste(cancer_name, "-", g)))
  }
  dev.off()
  
  invisible(NULL)
}

## -----------------------------------------------------------------------
## 8B) Plots de clusters e tipos celulares
## -----------------------------------------------------------------------
save_cluster_annotation_plots <- function(integrated_obj, outdir, cancer_name) {
  ensure_directory(outdir)
  if (!"umap" %in% names(integrated_obj@reductions)) {
    message("  UMAP não disponível. Pulando plots de cluster.")
    return(invisible(NULL))
  }
  
  p_cl <- DimPlot(integrated_obj, reduction = "umap",
                  group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
    ggtitle(paste0(cancer_name, " - clusters"))
  ggsave(file.path(outdir, paste0(cancer_name, "_UMAP_clusters.pdf")),
         p_cl, width = 8, height = 6)
  
  if ("cell_type_auto" %in% colnames(integrated_obj@meta.data)) {
    p_ct <- DimPlot(integrated_obj, reduction = "umap",
                    group.by = "cell_type_auto", label = TRUE, repel = TRUE) +
      ggtitle(paste0(cancer_name, " - cell types"))
    ggsave(file.path(outdir, paste0(cancer_name, "_UMAP_cell_types.pdf")),
           p_ct, width = 10, height = 7)
  }
  invisible(NULL)
}

## -----------------------------------------------------------------------
## 8C) Sumário de distribuição de assinatura
## -----------------------------------------------------------------------
summarize_signature_distribution <- function(integrated_obj, sig_feat) {
  meta <- integrated_obj@meta.data
  if (!sig_feat %in% colnames(meta))
    stop("Coluna de assinatura não encontrada: ", sig_feat)
  
  meta$cluster        <- as.character(integrated_obj$seurat_clusters)
  meta$score          <- meta[[sig_feat]]
  meta$cell_type_auto <- if (!is.null(meta$cell_type_auto)) meta$cell_type_auto else "Unassigned"
  
  by_cluster <- meta %>%
    dplyr::group_by(cluster, cell_type_auto) %>%
    dplyr::summarise(
      n_cells      = dplyr::n(),
      score_mean   = mean(score, na.rm = TRUE),
      score_median = median(score, na.rm = TRUE),
      score_sd     = sd(score, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(score_mean))
  
  by_celltype <- meta %>%
    dplyr::group_by(cell_type_auto) %>%
    dplyr::summarise(
      n_cells      = dplyr::n(),
      score_mean   = mean(score, na.rm = TRUE),
      score_median = median(score, na.rm = TRUE),
      score_sd     = sd(score, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(score_mean))
  
  list(summary_by_cluster = by_cluster, summary_by_celltype = by_celltype)
}

## -----------------------------------------------------------------------
## 9) build_integrated_object
## -----------------------------------------------------------------------
build_integrated_object <- function(cancer_dir, cancer_name,
                                    samples_to_use = NULL) {
  message("\n======================================================")
  message("Construindo objeto integrado: ", cancer_name)
  message("======================================================")
  
  if (is.null(samples_to_use)) {
    matrix_files <- list.files(cancer_dir, "_matrix\\.mtx\\.gz$",
                               full.names = TRUE)
    if (length(matrix_files) == 0) {
      message("Nenhum *_matrix.mtx.gz encontrado em ", cancer_dir)
      return(NULL)
    }
    samples <- sort(unique(sub("_matrix\\.mtx\\.gz$", "",
                               basename(matrix_files))))
  } else {
    samples <- samples_to_use
  }
  
  message("Amostras a processar: ", length(samples))
  
  seurat_objects <- list()
  
  for (sample_id in samples) {
    mtx_path  <- file.path(cancer_dir, paste0(sample_id, "_matrix.mtx.gz"))
    feat_path <- file.path(cancer_dir, paste0(sample_id, "_features.tsv.gz"))
    bc_path   <- file.path(cancer_dir, paste0(sample_id, "_barcodes.tsv.gz"))
    
    if (!all(file.exists(mtx_path, feat_path, bc_path))) {
      message("Arquivos incompletos para ", sample_id, ". Pulando.")
      next
    }
    
    expr_matrix <- as(readMM(mtx_path), "CsparseMatrix")
    genes_tb    <- read.table(feat_path, header = FALSE, stringsAsFactors = FALSE)
    barcodes    <- read.table(bc_path,   header = FALSE, stringsAsFactors = FALSE)
    
    if (nrow(expr_matrix) != nrow(genes_tb) || ncol(expr_matrix) != nrow(barcodes))
      stop("Inconsistência de dimensões na amostra ", sample_id)
    
    feature_df <- data.frame(
      feature_id   = trimws(as.character(genes_tb[[1]])),
      feature_name = trimws(if (ncol(genes_tb) >= 2)
        as.character(genes_tb[[2]])
        else as.character(genes_tb[[1]])),
      stringsAsFactors = FALSE
    )
    feature_df$feature_name[feature_df$feature_name == ""] <-
      feature_df$feature_id[feature_df$feature_name == ""]
    feature_df$feature_type <- ifelse(
      grepl("^ENST", feature_df$feature_id, ignore.case = TRUE),
      "transcript", "gene_or_other"
    )
    
    rownames(expr_matrix) <- make.unique(feature_df$feature_id)
    colnames(expr_matrix) <- paste(sample_id, barcodes$V1, sep = "_")
    feature_df$seurat_feature <- rownames(expr_matrix)
    
    seu <- CreateSeuratObject(
      counts       = expr_matrix,
      project      = paste0(cancer_name, "_", sample_id),
      min.cells    = 2,
      min.features = 200
    )
    seu$sample_id   <- sample_id
    seu$cancer_type <- cancer_name
    
    fm <- feature_df[match(rownames(seu), feature_df$seurat_feature), , drop = FALSE]
    rownames(fm) <- fm$seurat_feature
    seu@misc$feature_map <- fm
    
    seu <- NormalizeData(seu, verbose = FALSE)
    seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)
    
    if (ncol(seu) > 35000) {
      set.seed(123)
      seu <- subset(seu, cells = sample(Cells(seu), 35000))
    }
    
    seurat_objects[[sample_id]] <- seu
    message("  ", sample_id, " | células=", ncol(seu), " | features=", nrow(seu),
            " | transcritos=",
            sum(fm$feature_type == "transcript", na.rm = TRUE))
  }
  
  if (length(seurat_objects) == 0) {
    message("Nenhuma amostra válida para ", cancer_name)
    return(NULL)
  }
  
  .post_process <- function(obj, is_integrated = FALSE) {
    DefaultAssay(obj) <- if (is_integrated) "integrated" else "RNA"
    
    obj <- tryCatch(
      ScaleData(obj, verbose = FALSE),
      error = function(e) {
        vf <- VariableFeatures(obj)
        if (length(vf) == 0)
          vf <- VariableFeatures(FindVariableFeatures(obj, verbose = FALSE))
        ScaleData(obj, features = vf, verbose = FALSE)
      }
    )
    
    obj <- tryCatch(
      RunPCA(obj, npcs = 50, verbose = FALSE),
      error = function(e) RunPCA(obj, npcs = 30, verbose = FALSE)
    )
    
    pca_emb <- tryCatch(Embeddings(obj, "pca"), error = function(e) NULL)
    if (is.null(pca_emb) || anyNA(pca_emb)) {
      warning("PCA com NA — pulando clustering/UMAP.")
      return(obj)
    }
    
    n_dims <- min(30L, ncol(pca_emb))
    obj <- tryCatch(
      FindNeighbors(obj, dims = seq_len(n_dims), verbose = FALSE),
      error = function(e)
        FindNeighbors(obj, dims = seq_len(min(20L, ncol(pca_emb))), verbose = FALSE)
    )
    obj <- tryCatch(
      FindClusters(obj, resolution = 0.1, verbose = FALSE),
      error = function(e) FindClusters(obj, resolution = 0.05, verbose = FALSE)
    )
    obj <- tryCatch(
      RunUMAP(obj, dims = seq_len(n_dims), umap.method = "uwot", verbose = FALSE),
      error = function(e) { message("UMAP omitido."); obj }
    )
    obj
  }
  
  if (length(seurat_objects) == 1) {
    message("Uma única amostra — sem integração.")
    integrated_obj <- .post_process(seurat_objects[[1]], is_integrated = FALSE)
    integrated_obj@misc$feature_map <- seurat_objects[[1]]@misc$feature_map
    
  } else {
    message("  Integração RPCA (", length(seurat_objects), " amostras)")
    features <- SelectIntegrationFeatures(seurat_objects, nfeatures = 3000)
    
    seurat_objects <- lapply(seurat_objects, function(o) {
      o <- ScaleData(o, features = features, verbose = FALSE)
      RunPCA(o, features = features, verbose = FALSE)
    })
    
    anchors <- tryCatch(
      FindIntegrationAnchors(seurat_objects, anchor.features = features,
                             reduction = "rpca"),
      error = function(e) { message("RPCA falhou: ", e$message); NULL }
    )
    
    if (is.null(anchors)) {
      message("Fallback: usando primeira amostra.")
      integrated_obj <- .post_process(seurat_objects[[1]], FALSE)
      integrated_obj@misc$feature_map <- seurat_objects[[1]]@misc$feature_map
    } else {
      integrated_obj <- IntegrateData(anchorset = anchors)
      
      if (!"integrated" %in% names(integrated_obj@assays)) {
        message("Assay 'integrated' não gerado. Fallback.")
        integrated_obj <- .post_process(seurat_objects[[1]], FALSE)
        integrated_obj@misc$feature_map <- seurat_objects[[1]]@misc$feature_map
      } else {
        integrated_obj <- .post_process(integrated_obj, is_integrated = TRUE)
        DefaultAssay(integrated_obj) <- "RNA"
        
        all_fm <- dplyr::bind_rows(lapply(seurat_objects, function(o) {
          fm <- o@misc$feature_map
          if (is.null(fm)) return(NULL)
          as.data.frame(fm, stringsAsFactors = FALSE)
        })) %>%
          dplyr::distinct(seurat_feature, .keep_all = TRUE) %>%
          dplyr::filter(seurat_feature %in% rownames(integrated_obj[["RNA"]]))
        
        if (nrow(all_fm) == 0)
          all_fm <- build_feature_map(integrated_obj, "RNA")
        rownames(all_fm) <- all_fm$seurat_feature
        integrated_obj@misc$feature_map <- all_fm
      }
    }
  }
  
  DefaultAssay(integrated_obj) <- "RNA"
  
  if (!"seurat_clusters" %in% colnames(integrated_obj@meta.data)) {
    message("Clustering ausente — criando cluster único.")
    integrated_obj$seurat_clusters <- factor(rep(1, ncol(integrated_obj)))
  }
  
  if (inherits(integrated_obj[["RNA"]], "Assay5")) {
    layers_now <- SeuratObject::Layers(integrated_obj, assay = "RNA")
    needs_join <- sum(grepl("^counts\\.", layers_now)) > 1 ||
      sum(grepl("^data\\.",   layers_now)) > 1
    if (needs_join) {
      integrated_obj <- tryCatch(JoinLayers(integrated_obj),
                                 error = function(e) {
                                   message("JoinLayers falhou: ", e$message)
                                   integrated_obj
                                 })
    }
  }
  
  markers <- tryCatch(
    FindAllMarkers(integrated_obj, only.pos = TRUE,
                   min.pct = 0.25, logfc.threshold = 0.25),
    error = function(e) {
      message("FindAllMarkers falhou: ", e$message)
      data.frame()
    }
  )
  
  ## CORREÇÃO: dplyr::slice_max e dplyr:: em vez de slice_max sem namespace
  top10_markers <- if (nrow(markers) > 0) {
    markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
      dplyr::ungroup()
  } else data.frame()
  
  list(
    integrated_obj = integrated_obj,
    valid_samples  = names(seurat_objects),
    markers        = markers,
    top10_markers  = top10_markers
  )
}

## -----------------------------------------------------------------------
## 10) process_one_signature
## -----------------------------------------------------------------------
process_one_signature <- function(
    integrated_obj_annotated,
    annotation_res_base,
    cancer_name,
    results_root,
    valid_samples,
    markers,
    top10_markers,
    sig_name,
    my_signature_up,
    my_signature_down = character(0),
    genes_to_test,
    rodada_name = NULL
) {
  sig_safe <- sanitize_filename(sig_name)
  
  if (!is.null(rodada_name)) {
    outdir <- file.path(results_root, cancer_name, rodada_name, sig_safe)
    cancer_name <- paste0(cancer_name, "_", rodada_name)
  } else {
    outdir <- file.path(results_root, cancer_name, sig_safe)
  }
  ensure_directory(outdir)
  
  message("\n--- Assinatura: ", sig_name, " (", cancer_name, ")")
  
  integrated_obj <- integrated_obj_annotated
  
  save_cluster_annotation_plots(integrated_obj, outdir, cancer_name)
  
  tabs <- list(
    cluster_marker_union_coverage = annotation_res_base$coverage_table,
    marker_mapping_report         = annotation_res_base$marker_mapping_report,
    panel_gene_status             = annotation_res_base$panel_gene_status,
    cluster_auto_annotations      = annotation_res_base$cluster_assignments,
    audit_unassigned_clusters     = annotation_res_base$audit_unassigned,
    audit_conflicting_clusters    = annotation_res_base$audit_conflicts
  )
  for (nm in names(tabs)) {
    if (!is.null(tabs[[nm]]) && nrow(tabs[[nm]]) > 0)
      safe_export(tabs[[nm]],
                  file.path(outdir, paste0(cancer_name, "_", nm, ".xlsx")))
  }
  
  if (nrow(markers) > 0)
    safe_export(markers,
                file.path(outdir, paste0(cancer_name, "_markers_all.xlsx")))
  
  if (nrow(top10_markers) > 0) {
    safe_export(top10_markers,
                file.path(outdir, paste0(cancer_name, "_top10_markers.xlsx")))
    
    top_feats <- map_queries_to_features(
      integrated_obj, unique(top10_markers$gene), "RNA"
    )$features
    if (length(top_feats) > 0) {
      pdf(file.path(outdir, paste0(cancer_name, "_Heatmap_top10.pdf")),
          width = 12, height = 10)
      print(DoHeatmap(integrated_obj, features = top_feats,
                      group.by = "seurat_clusters") + NoLegend())
      dev.off()
    }
  }
  
  save_marker_exploration(integrated_obj, outdir, cancer_name, genes_to_test)
  
  integrated_obj <- score_signature(
    obj        = integrated_obj,
    up         = my_signature_up,
    down       = my_signature_down,
    sig_name   = sig_safe,
    assay_expr = "RNA"
  )
  
  sig_feat <- paste0(sig_safe, "_Composite")
  
  for (direction in c("up", "down")) {
    mr <- integrated_obj@misc[[paste0(sig_safe, "_mapping_report_", direction)]]
    if (!is.null(mr) && nrow(mr) > 0)
      safe_export(mr, file.path(outdir,
                                paste0(cancer_name, "_", sig_safe,
                                       "_mapping_", direction, ".xlsx")))
  }
  
  if ("umap" %in% names(integrated_obj@reductions) &&
      sig_feat %in% colnames(integrated_obj@meta.data)) {
    p_fp <- FeaturePlot(integrated_obj, features = sig_feat,
                        reduction = "umap") + scale_color_viridis_c()
    ggsave(file.path(outdir, paste0(cancer_name, "_", sig_safe, "_FeaturePlot.pdf")),
           p_fp, width = 8, height = 6)
  }
  
  if (sig_feat %in% colnames(integrated_obj@meta.data)) {
    p_vln <- VlnPlot(integrated_obj, features = sig_feat,
                     group.by = "seurat_clusters", pt.size = 0)
    ggsave(file.path(outdir, paste0(cancer_name, "_", sig_safe, "_Vln_cluster.pdf")),
           p_vln, width = 10, height = 6)
    
    if ("cell_type_auto" %in% colnames(integrated_obj@meta.data)) {
      p_ct <- VlnPlot(integrated_obj, features = sig_feat,
                      group.by = "cell_type_auto", pt.size = 0)
      ggsave(file.path(outdir, paste0(cancer_name, "_", sig_safe, "_Vln_celltype.pdf")),
             p_ct, width = 12, height = 6)
    }
  }
  
  sig_sum <- summarize_signature_distribution(integrated_obj, sig_feat)
  for (nm in names(sig_sum)) {
    if (nrow(sig_sum[[nm]]) > 0)
      safe_export(sig_sum[[nm]],
                  file.path(outdir,
                            paste0(cancer_name, "_", sig_safe, "_", nm, ".xlsx")))
  }
  
  
  safe_export(integrated_obj@meta.data,
              file.path(outdir, paste0(cancer_name, "_", sig_safe, "_metadata.xlsx")))
  
  map_up    <- integrated_obj@misc[[paste0(sig_safe, "_mapping_report_up")]]
  n_found_up <- if (!is.null(map_up)) sum(!is.na(map_up$matched_feature)) else 0L
  
  report_df <- data.frame(
    cancer_type              = cancer_name,
    n_samples                = length(valid_samples),
    samples                  = paste(valid_samples, collapse = ", "),
    n_cells_total            = ncol(integrated_obj),
    n_features_total         = nrow(integrated_obj),
    signature_name           = sig_name,
    signature_name_safe      = sig_safe,
    signature_items_original = paste(my_signature_up, collapse = " + "),
    signature_input_type     = classify_signature_input(my_signature_up),
    n_signature_items_input  = length(my_signature_up),
    n_signature_items_found  = n_found_up,
    prop_found               = ifelse(length(my_signature_up) > 0,
                                      n_found_up / length(my_signature_up),
                                      NA_real_),
    n_clusters_total         = length(unique(as.character(integrated_obj$seurat_clusters))),
    n_clusters_unassigned    = nrow(annotation_res_base$audit_unassigned),
    n_clusters_conflicting   = nrow(annotation_res_base$audit_conflicts),
    stringsAsFactors         = FALSE
  )
  
  safe_export(report_df,
              file.path(outdir, paste0(cancer_name, "_", sig_safe, "_report.xlsx")))
  
  
  message("  OK: ", cancer_name, " / ", sig_name)
  
  rm(integrated_obj)
  gc(verbose = FALSE)
  
  invisible(report_df)
}

## -----------------------------------------------------------------------
## 11) Checkpoint + timer global
## -----------------------------------------------------------------------
cp_paths     <- initialize_checkpoint(results_root)
script_start <- Sys.time()
cp_loaded    <- load_or_create_checkpoint(cp_paths$checkpoint_file)
checkpoint   <- cp_loaded$checkpoint
resumed      <- cp_loaded$resumed

session_id <- paste0("session_", format(script_start, "%Y%m%d_%H%M%S"))

checkpoint$session_history <- dplyr::bind_rows(
  checkpoint$session_history,
  data.frame(
    session_id                         = session_id,
    session_start                      = format(script_start, "%Y-%m-%d %H:%M:%S"),
    resumed_from_checkpoint            = resumed,
    last_finished_cancer_before_resume =
      ifelse(is.na(checkpoint$last_finished_cancer), "",
             checkpoint$last_finished_cancer),
    session_end                        = NA_character_,
    session_elapsed_seconds            = NA_real_,
    session_elapsed_formatted          = NA_character_,
    stringsAsFactors                   = FALSE
  )
)
save_checkpoint_files(checkpoint, cp_paths)

message("\n======================================================")
message("INÍCIO: ", format(script_start, "%Y-%m-%d %H:%M:%S"),
        if (resumed) "  [RETOMADA DE CHECKPOINT]" else "  [EXECUÇÃO NOVA]")
message("Filtro ativo: apenas assinaturas de camada ômica 6")
message("======================================================\n")

## -----------------------------------------------------------------------
## 12) Descobre cânceres
## -----------------------------------------------------------------------
cancer_dirs  <- list.dirs(samples_root, full.names = TRUE, recursive = FALSE)
cancer_names <- basename(cancer_dirs)

if (length(cancer_dirs) == 0)
  stop("Nenhuma pasta de câncer em: ", samples_root)

## Restringe os cânceres processados apenas aos que têm assinaturas
## de camada ômica 6 — evita processar amostras sem utilidade
cancers_with_sigs <- unique(signatures_df$type)
cancer_mask       <- cancer_names %in% cancers_with_sigs
cancer_dirs       <- cancer_dirs[cancer_mask]
cancer_names      <- cancer_names[cancer_mask]

if (length(cancer_dirs) == 0)
  stop("Nenhuma pasta de câncer corresponde aos tipos com assinaturas de camada 6.")

message("Cânceres com assinaturas de camada 6: ",
        paste(cancer_names, collapse = ", "))
if (length(checkpoint$completed_cancers) > 0)
  message("Já concluídos: ", paste(checkpoint$completed_cancers, collapse = ", "))

## -----------------------------------------------------------------------
## 13) Loop principal por câncer
## -----------------------------------------------------------------------
for (i in seq_along(cancer_dirs)) {
  cancer_dir  <- cancer_dirs[i]
  cancer_name <- cancer_names[i]
  
  matrix_files <- list.files(cancer_dir, "_matrix\\.mtx\\.gz$", full.names = TRUE)
  if (length(matrix_files) == 0) {
    message("\n[SKIP] ", cancer_name, ": Nenhum *_matrix.mtx.gz encontrado.")
    next
  }
  
  all_samples <- sort(unique(sub("_matrix\\.mtx\\.gz$", "", basename(matrix_files))))
  
  batch_size <- max_samples_per_cancer
  if (is.infinite(batch_size)) batch_size <- length(all_samples)
  num_batches <- ceiling(length(all_samples) / batch_size)
  
  for (b in seq_len(num_batches)) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx   <- min(b * batch_size, length(all_samples))
    batch_samples <- all_samples[start_idx:end_idx]
    
    rodada_name <- paste0("rodada_", b, "_n", length(batch_samples))
    rodada_id   <- paste0(cancer_name, "_", rodada_name)
    
    if (rodada_id %in% checkpoint$completed_cancers) {
      message("\n[SKIP] ", rodada_id, " já totalmente concluído.")
      next
    }
    
    t0 <- Sys.time()
    message("\n########## ", cancer_name, " (", rodada_name, ") ##########")
    
    cancer_sigs <- dplyr::filter(signatures_df, type == cancer_name)
    
    if (nrow(cancer_sigs) == 0) {
      message("Sem assinaturas de camada 6 para ", cancer_name, ". Pulando.")
      checkpoint <- .log_cancer(checkpoint, session_id, rodada_id,
                                "skipped_no_signature", t0)
      checkpoint$completed_cancers <- unique(c(checkpoint$completed_cancers, rodada_id))
      save_checkpoint_files(checkpoint, cp_paths)
      next
    }
    
    message("  Assinaturas de camada 6 encontradas: ", nrow(cancer_sigs))
    
    if (is.null(checkpoint$completed_signatures[[rodada_id]])) {
      checkpoint$completed_signatures[[rodada_id]] <- character(0)
    }
    
    sigs_to_run <- cancer_sigs$nomeclatura[!cancer_sigs$nomeclatura %in% checkpoint$completed_signatures[[rodada_id]]]
    
    if (length(sigs_to_run) == 0) {
      message("  Todas as assinaturas desta rodada já foram processadas.")
      checkpoint$completed_cancers <- unique(c(checkpoint$completed_cancers, rodada_id))
      save_checkpoint_files(checkpoint, cp_paths)
      next
    }
    
    message("  Assinaturas restantes para processar: ", length(sigs_to_run))
    
    bundle <- tryCatch(
      build_integrated_object(cancer_dir, rodada_id, samples_to_use = batch_samples),
      error = function(e) {
        message("[ERRO] build_integrated_object: ", e$message)
        NULL
      }
    )
    
    if (is.null(bundle)) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      checkpoint$cancer_timing_log <- dplyr::bind_rows(
        checkpoint$cancer_timing_log,
        data.frame(
          session_id        = session_id, cancer_type = rodada_id,
          status            = "failed_build",
          start_time        = format(t0, "%Y-%m-%d %H:%M:%S"),
          end_time          = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          elapsed_seconds   = elapsed,
          elapsed_formatted = format_elapsed_time(elapsed),
          stringsAsFactors  = FALSE
        )
      )
      save_checkpoint_files(checkpoint, cp_paths)
      next
    }
    
    annot <- tryCatch(
      annotate_clusters_by_marker_union(
        bundle$integrated_obj, cancer_name, markers_ref,
        assay_expr = "RNA", slot_expr = "data",
        expr_threshold = 0, union_threshold = 0.50, min_detected_genes = 2
      ),
      error = function(e) {
        message("[ERRO] annotate_clusters: ", e$message)
        NULL
      }
    )
    
    if (is.null(annot)) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      checkpoint$cancer_timing_log <- dplyr::bind_rows(
        checkpoint$cancer_timing_log,
        data.frame(
          session_id        = session_id, cancer_type = rodada_id,
          status            = "failed_annotation",
          start_time        = format(t0, "%Y-%m-%d %H:%M:%S"),
          end_time          = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          elapsed_seconds   = elapsed,
          elapsed_formatted = format_elapsed_time(elapsed),
          stringsAsFactors  = FALSE
        )
      )
      save_checkpoint_files(checkpoint, cp_paths)
      next
    }
    
    integrated_obj_ann <- annot$integrated_obj
    
    for (j in seq_len(nrow(cancer_sigs))) {
      sig_name   <- cancer_sigs$nomeclatura[j]
      sig_string <- cancer_sigs$assinaturas[j]
      sig_up     <- parse_signature_features(sig_string)
      
      if (sig_name %in% checkpoint$completed_signatures[[rodada_id]]) {
        next
      }
      
      if (length(sig_up) == 0) {
        message("  Assinatura vazia: ", sig_name, ". Pulando.")
        checkpoint$completed_signatures[[rodada_id]] <- unique(c(checkpoint$completed_signatures[[rodada_id]], sig_name))
        save_checkpoint_files(checkpoint, cp_paths)
        next
      }
      
      res <- tryCatch(
        process_one_signature(
          integrated_obj_annotated = integrated_obj_ann,
          annotation_res_base      = annot,
          cancer_name              = cancer_name,
          results_root             = results_root,
          valid_samples            = bundle$valid_samples,
          markers                  = bundle$markers,
          top10_markers            = bundle$top10_markers,
          sig_name                 = sig_name,
          my_signature_up          = sig_up,
          my_signature_down        = character(0),
          genes_to_test            = genes_to_test,
          rodada_name              = rodada_name
        ),
        error = function(e) {
          message("[ERRO] ", rodada_id, "/", sig_name, ": ", e$message)
          NULL
        }
      )
      
      rm(res)
      gc(verbose = FALSE)
      
      # Salvar checkpoint após cada assinatura
      checkpoint$completed_signatures[[rodada_id]] <- unique(c(checkpoint$completed_signatures[[rodada_id]], sig_name))
      save_checkpoint_files(checkpoint, cp_paths)
    }
    
    checkpoint$completed_cancers <- unique(c(checkpoint$completed_cancers, rodada_id))
    checkpoint$last_finished_cancer <- rodada_id
    
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    checkpoint$cancer_timing_log <- dplyr::bind_rows(
      checkpoint$cancer_timing_log,
      data.frame(
        session_id        = session_id, cancer_type = rodada_id,
        status            = "completed",
        start_time        = format(t0, "%Y-%m-%d %H:%M:%S"),
        end_time          = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        elapsed_seconds   = elapsed,
        elapsed_formatted = format_elapsed_time(elapsed),
        stringsAsFactors  = FALSE
      )
    )
    
    save_checkpoint_files(checkpoint, cp_paths)
    
    rm(bundle, annot, integrated_obj_ann)
    gc(verbose = FALSE)
    
    message("CONCLUÍDO: ", rodada_id, " | Tempo: ", format_elapsed_time(elapsed))
  }
}

## -----------------------------------------------------------------------
## 14) Finalização
## -----------------------------------------------------------------------
script_end     <- Sys.time()
script_elapsed <- as.numeric(difftime(script_end, script_start, units = "secs"))

idx <- which(checkpoint$session_history$session_id == session_id)
if (length(idx) == 1) {
  checkpoint$session_history$session_end[idx]               <- format(script_end, "%Y-%m-%d %H:%M:%S")
  checkpoint$session_history$session_elapsed_seconds[idx]   <- script_elapsed
  checkpoint$session_history$session_elapsed_formatted[idx] <- format_elapsed_time(script_elapsed)
}

save_checkpoint_files(checkpoint, cp_paths)

message("\n======================================================")
message("FIM: ", format(script_end, "%Y-%m-%d %H:%M:%S"))
message("TEMPO TOTAL: ", format_elapsed_time(script_elapsed))
message("Cânceres concluídos: ", paste(checkpoint$completed_cancers, collapse = ", "))
message("======================================================")