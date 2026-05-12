# ======================================================================================
# RSF PIPELINE (RSF + Boruta + Clinicas opcionais) COM MAPA DE COBERTURA POR ENDPOINT
# + SOMENTE MULTIFENOTIPO (SEM MONOFENOTIPO)
# + ÊNFASE EM C-INDEX E AUC TEMPORAL
# + PAINEL DE IMPORTÂNCIA APENAS NOS MODOS SEM BORUTA
# + TABELA DE VARIÁVEIS USADAS EM CADA PREDIÇÃO
# + MODOS COM BORUTA: exporta Importance + RelativeImp
# + MODOS SEM BORUTA:
#       RSF PRELIMINAR -> CONSENSUS IMPORTANCE -> SELEÇÃO -> RSF FINAL
#       exporta painel preliminar + variáveis selecionadas + painel final
# + CORREÇÃO:
#       minimal_depth <- extraído de max.subtree()$order
#       split_frequency <- proxy via max.subtree()$count
# + NOVO:
#       inclusão automática de variáveis fenotípicas extras conforme os fenótipos
#       presentes nas assinaturas selecionadas
# + CORREÇÃO DE PROGRESSO:
#       progress key agora inclui modo_vars + omic_vars
# ======================================================================================

suppressPackageStartupMessages({
  library(rio)
  library(dplyr)
  library(randomForest)
  library(survival)
  library(survminer)
  library(caret)
  library(pROC)
  library(Boruta)
  library(randomForestSRC)
  library(timeROC)
  library(gridExtra)
  library(grid)
  library(fmsb)
  library(scales)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(corrplot)
  library(survex)
  library(doParallel)
  library(foreach)
})

# 📂 Diretório
setwd("D:/ML_Emanuell/ML_RSF/linear")

# ======================================================================================
# TEMPO GLOBAL DA EXECUÇÃO - INÍCIO
# ======================================================================================
execucao_inicio <- Sys.time()

formatar_duracao_global <- function(inicio, fim) {
  dur <- as.numeric(difftime(fim, inicio, units = "secs"))
  
  horas <- floor(dur / 3600)
  minutos <- floor((dur %% 3600) / 60)
  segundos <- dur %% 60
  
  paste0(
    sprintf("%02d", horas), "h ",
    sprintf("%02d", minutos), "min ",
    sprintf("%05.2f", segundos), "s"
  )
}

cat("\n")
cat("============================================================\n")
cat("🕒 INÍCIO DA EXECUÇÃO GLOBAL\n")
cat("Horário de início:", format(execucao_inicio, "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n\n")

# ======================================================================================
# 0) PASTA DE PROGRESSO FORA DE output
# ======================================================================================
progress_dir <- file.path(getwd(), "progress")
if (!dir.exists(progress_dir)) {
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
}

# NOVO NOME para não reaproveitar chaves antigas incorretas
progress_file <- file.path(progress_dir, "progress_global_rsf_consensus_v2.rds")

# ======================================================================================
# 1) ARQUIVO DE MAPEAMENTO: qual df usar para cada (cancer_type, metric)
# ======================================================================================
coverage_map <- readr::read_tsv("final_ML.tsv", show_col_types = FALSE)

coverage_map <- coverage_map %>%
  mutate(
    cancer_type = as.character(cancer_type),
    metric      = as.character(metric),
    assigned_df = as.character(assigned_df)
  )

required_cols <- c("cancer_type", "metric", "assigned_df")
missing_required <- setdiff(required_cols, colnames(coverage_map))
if (length(missing_required) > 0) {
  stop("❌ coverage_map.tsv não possui as colunas obrigatórias: ",
       paste(missing_required, collapse = ", "))
}

coverage_map <- coverage_map %>%
  filter(!is.na(assigned_df), assigned_df != "")

if (nrow(coverage_map) == 0) {
  stop("❌ coverage_map.tsv não possui linhas válidas com assigned_df.")
}

dup_pairs <- coverage_map %>%
  count(cancer_type, metric) %>%
  filter(n > 1)

if (nrow(dup_pairs) > 0) {
  cat("⚠️ Foram encontradas combinações duplicadas em coverage_map.tsv.\n")
  print(dup_pairs)
  stop("❌ Cada combinação (cancer_type, metric) deve apontar para apenas um assigned_df.")
}

# ======================================================================================
# 2) CONFIGURAÇÕES AUXILIARES
# ======================================================================================

# Fenótipo = 4º componente global da assinatura
# Após remover o prefixo do câncer (ex: "BLCA-"), isso vira o 3º componente restante
FENOTIPO_COMPONENT_POS <- 3

pheno_demog_map <- list(
  pheno1_TMB = c("Non_silent_per_Mb"),
  pheno2_MSI = c(
    "Total_nb_MSI_events", "MSI_3utr", "MSI_5utr", "MSI_mono",
    "MSI_di", "MSI_tri", "MSI_tetra",
    "MSI_category_nb_from_TCGA_consortium"
  ),
  pheno3_stemness = c("RNAss")
)

# ======================================================================================
# 3) FUNÇÕES AUXILIARES
# ======================================================================================
safe_select <- function(df, cols) {
  cols_ok <- intersect(cols, colnames(df))
  if (length(cols_ok) == 0) return(NULL)
  df %>% dplyr::select(all_of(cols_ok))
}

# CHAVE CORRIGIDA
make_progress_key <- function(modo_vars, omic_vars, tipo, metrica, omic_label, pheno_label) {
  paste(modo_vars, omic_vars, tipo, metrica, omic_label, pheno_label, sep = "||")
}

get_df_file_for_pair <- function(tipo_cancer, metrica, coverage_map, base_dir = getwd()) {
  hit <- coverage_map %>%
    filter(cancer_type == tipo_cancer, metric == metrica)
  
  if (nrow(hit) == 0) {
    return(NA_character_)
  }
  
  df_file <- hit$assigned_df[1]
  
  if (!grepl("\\.rds$", df_file, ignore.case = TRUE)) {
    df_file <- paste0(df_file, ".rds")
  }
  
  if (file.exists(df_file)) {
    return(normalizePath(df_file, winslash = "/", mustWork = FALSE))
  }
  
  full_path <- file.path(base_dir, df_file)
  if (file.exists(full_path)) {
    return(normalizePath(full_path, winslash = "/", mustWork = FALSE))
  }
  
  return(full_path)
}

df_cache <- new.env(parent = emptyenv())

load_df_cached <- function(path) {
  key <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (exists(key, envir = df_cache, inherits = FALSE)) {
    return(get(key, envir = df_cache, inherits = FALSE))
  }
  obj <- readRDS(path)
  assign(key, obj, envir = df_cache)
  obj
}

save_global_progress <- function(progress_file, resultados_geral, progress_keys) {
  progress_obj <- list(
    resultados_geral   = resultados_geral,
    progress_keys      = progress_keys,
    ultima_atualizacao = Sys.time()
  )
  saveRDS(progress_obj, progress_file)
}

# ------------------------------------------------------------------
# FUNÇÕES PARA FENÓTIPOS DAS ASSINATURAS
# ------------------------------------------------------------------
extrair_fenotipo_assinatura <- function(nome_assinatura, posicao_fenotipo = FENOTIPO_COMPONENT_POS) {
  x <- as.character(nome_assinatura)
  
  # remove apenas o prefixo do tipo tumoral, seja com "-" ou "_"
  x2 <- sub("^[^_-]+[_-]", "", x)
  
  partes <- strsplit(x2, "\\.")[[1]]
  
  if (length(partes) < posicao_fenotipo) return(NA_character_)
  
  partes[posicao_fenotipo]
}

obter_fenotipos_presentes <- function(selected_vars, posicao_fenotipo = FENOTIPO_COMPONENT_POS) {
  if (length(selected_vars) == 0) return(character(0))
  
  phenos <- vapply(
    selected_vars,
    FUN = extrair_fenotipo_assinatura,
    FUN.VALUE = character(1),
    posicao_fenotipo = posicao_fenotipo
  )
  
  phenos <- unique(phenos[!is.na(phenos) & phenos %in% c("1", "2", "3")])
  phenos
}

obter_variaveis_fenotipicas_extras <- function(selected_vars, posicao_fenotipo = FENOTIPO_COMPONENT_POS) {
  phenos_presentes <- obter_fenotipos_presentes(
    selected_vars = selected_vars,
    posicao_fenotipo = posicao_fenotipo
  )
  
  vars_extra <- character(0)
  
  if ("1" %in% phenos_presentes) {
    vars_extra <- c(vars_extra, pheno_demog_map$pheno1_TMB)
  }
  if ("2" %in% phenos_presentes) {
    vars_extra <- c(vars_extra, pheno_demog_map$pheno2_MSI)
  }
  if ("3" %in% phenos_presentes) {
    vars_extra <- c(vars_extra, pheno_demog_map$pheno3_stemness)
  }
  
  unique(vars_extra)
}

# ======================================================================================
# 3A) TABELA DE VARIÁVEIS USADAS
# ======================================================================================
montar_tabela_variaveis_usadas <- function(
    predictor_names,
    selected_vars,
    clinicas_relevantes,
    tipo_cancer,
    metrica,
    omic_label,
    pheno_label,
    modo_vars,
    omic_root,
    pheno_root,
    importance_boruta_df = NULL,
    importance_panel_df = NULL,
    thr_vimp = NA_real_,
    thr_depth = NA_real_,
    thr_split = NA_real_
) {
  
  tab <- data.frame(
    Variable = predictor_names,
    stringsAsFactors = FALSE
  )
  
  tab <- tab %>%
    dplyr::mutate(
      Tipo_variavel = ifelse(Variable %in% clinicas_relevantes, "clinica", "omica"),
      Origem_selecao = dplyr::case_when(
        Variable %in% clinicas_relevantes ~ "Cox_clinico",
        Variable %in% selected_vars & modo_vars %in% c(1, 2) ~ "Boruta",
        Variable %in% selected_vars & modo_vars %in% c(3, 4) ~ "Consensus_RSF",
        TRUE ~ "Outra"
      ),
      Usada_no_modelo_final = TRUE,
      Tipo = tipo_cancer,
      Metrica = metrica,
      Omic_layer = omic_label,
      Pheno_layer = pheno_label,
      Modo_vars = modo_vars,
      Omic_mode = omic_root,
      Pheno_mode = pheno_root
    )
  
  if (modo_vars %in% c(1, 2)) {
    
    if (is.null(importance_boruta_df) || nrow(importance_boruta_df) == 0) {
      importance_boruta_df <- data.frame(
        Variable = predictor_names,
        Importance = NA_real_,
        RelativeImp = NA_real_,
        stringsAsFactors = FALSE
      )
    }
    
    tab <- tab %>%
      dplyr::left_join(importance_boruta_df, by = "Variable")
    
    if (!"Importance" %in% names(tab))  tab$Importance  <- NA_real_
    if (!"RelativeImp" %in% names(tab)) tab$RelativeImp <- NA_real_
  }
  
  if (modo_vars %in% c(3, 4)) {
    
    if (is.null(importance_panel_df) || nrow(importance_panel_df) == 0) {
      importance_panel_df <- data.frame(
        Variable = predictor_names,
        VIMP = NA_real_,
        minimal_depth = NA_real_,
        split_frequency = NA_real_,
        crit_vimp = NA,
        crit_depth = NA,
        crit_split = NA,
        score_importancia = NA_integer_,
        importante_consenso = NA,
        rank_vimp = NA_real_,
        rank_depth = NA_real_,
        rank_split = NA_real_,
        rank_consenso = NA_real_,
        stringsAsFactors = FALSE
      )
    }
    
    tab <- tab %>%
      dplyr::left_join(importance_panel_df, by = "Variable") %>%
      dplyr::mutate(
        Regra_VIMP = "VIMP > 0",
        Regra_MinDepth = ifelse(is.na(thr_depth),
                                "minimal depth < média depth (NA)",
                                paste0("minimal depth < média depth = ", round(thr_depth, 4))),
        Regra_SplitFreq = ifelse(is.na(thr_split),
                                 "split frequency > percentil 75 (NA)",
                                 paste0("split frequency > percentil 75 = ", round(thr_split, 4))),
        Threshold_VIMP = thr_vimp,
        Threshold_MinDepth = thr_depth,
        Threshold_SplitFreq = thr_split
      )
  }
  
  tab
}

# ======================================================================================
# 3B) FUNÇÕES AUXILIARES PARA IMPORTÂNCIA
# ======================================================================================
extrair_importancia_boruta_modo <- function(rsf_model, predictor_names) {
  out <- data.frame(
    Variable = predictor_names,
    Importance = NA_real_,
    RelativeImp = NA_real_,
    stringsAsFactors = FALSE
  )
  
  imp <- tryCatch(rsf_model$importance, error = function(e) NULL)
  
  if (is.null(imp) || length(imp) == 0) {
    return(out)
  }
  
  tmp <- data.frame(
    Variable = names(imp),
    VIMP_raw = as.numeric(imp),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(VIMP_raw))
  
  tmp <- tmp %>%
    dplyr::mutate(
      Importance = pmax(VIMP_raw, 0)
    )
  
  mx <- suppressWarnings(max(tmp$Importance, na.rm = TRUE))
  if (!is.finite(mx) || is.na(mx) || mx <= 0) {
    tmp$RelativeImp <- NA_real_
  } else {
    tmp$RelativeImp <- tmp$Importance / mx
  }
  
  out <- out %>%
    dplyr::select(Variable) %>%
    dplyr::left_join(
      tmp %>% dplyr::select(Variable, Importance, RelativeImp),
      by = "Variable"
    )
  
  out
}

extrair_vimp_df <- function(rsf_model) {
  out <- data.frame(
    Variable = character(0),
    VIMP = numeric(0),
    stringsAsFactors = FALSE
  )
  
  imp <- tryCatch(rsf_model$importance, error = function(e) NULL)
  if (is.null(imp) || length(imp) == 0) return(out)
  
  out <- data.frame(
    Variable = names(imp),
    VIMP = as.numeric(imp),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(VIMP))
  
  out
}

extrair_min_depth_df <- function(rsf_model) {
  out <- data.frame(
    Variable = character(0),
    minimal_depth = numeric(0),
    stringsAsFactors = FALSE
  )
  
  md_obj <- tryCatch({
    randomForestSRC::var.select(rsf_model, method = "md", verbose = FALSE)
  }, error = function(e) {
    cat("⚠️ Não foi possível extrair minimal depth via var.select():", e$message, "\n")
    NULL
  })
  
  if (is.null(md_obj)) return(out)
  
  if (!is.null(md_obj$order) && is.matrix(md_obj$order)) {
    ord <- md_obj$order
    
    var_names <- rownames(ord)
    
    if (is.null(var_names) || length(var_names) == 0) {
      if (!is.null(md_obj$topvars) && length(md_obj$topvars) == nrow(ord)) {
        var_names <- md_obj$topvars
      } else if (!is.null(names(md_obj$count)) && length(md_obj$count) == nrow(ord)) {
        var_names <- names(md_obj$count)
      } else {
        var_names <- paste0("Var", seq_len(nrow(ord)))
      }
    }
    
    out <- data.frame(
      Variable = as.character(var_names),
      minimal_depth = as.numeric(ord[, 1]),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::filter(!is.na(Variable), !is.na(minimal_depth))
    
    return(out)
  }
  
  cand <- NULL
  
  if (is.data.frame(md_obj)) {
    cand <- md_obj
  } else if (is.list(md_obj)) {
    for (nm in c("depth", "varselect", "md.obj", "topvars")) {
      if (!is.null(md_obj[[nm]]) && is.data.frame(md_obj[[nm]])) {
        cand <- md_obj[[nm]]
        break
      }
    }
  }
  
  if (is.null(cand) || nrow(cand) == 0) return(out)
  
  cn <- colnames(cand)
  var_col <- cn[cn %in% c("variable", "Variable", "var", "xvar", "name")]
  depth_col <- cn[cn %in% c("depth", "Depth", "minimal.depth", "minimal_depth")]
  
  if (length(var_col) == 0 && !is.null(rownames(cand))) {
    cand$Variable <- rownames(cand)
    var_col <- "Variable"
  }
  
  if (length(var_col) == 0 || length(depth_col) == 0) return(out)
  
  out <- data.frame(
    Variable = as.character(cand[[var_col[1]]]),
    minimal_depth = as.numeric(cand[[depth_col[1]]]),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(Variable), !is.na(minimal_depth))
  
  out
}

extrair_splitfreq_df <- function(rsf_model, predictor_names) {
  out <- data.frame(
    Variable = predictor_names,
    split_frequency = NA_real_,
    stringsAsFactors = FALSE
  )
  
  md_obj <- tryCatch({
    randomForestSRC::max.subtree(rsf_model)
  }, error = function(e) NULL)
  
  if (!is.null(md_obj) && !is.null(md_obj$count)) {
    cnt <- md_obj$count
    
    tmp <- data.frame(
      Variable = names(cnt),
      split_frequency = as.numeric(cnt),
      stringsAsFactors = FALSE
    )
    
    out <- out %>%
      dplyr::left_join(tmp, by = "Variable", suffix = c("", ".tmp")) %>%
      dplyr::mutate(
        split_frequency = dplyr::coalesce(split_frequency.tmp, split_frequency)
      ) %>%
      dplyr::select(Variable, split_frequency)
    
    return(out)
  }
  
  cat("⚠️ Não foi possível extrair split_frequency; ficará NA.\n")
  out
}

construir_painel_importancia <- function(rsf_model, predictor_names) {
  vimp_df  <- extrair_vimp_df(rsf_model)
  depth_df <- extrair_min_depth_df(rsf_model)
  split_df <- extrair_splitfreq_df(rsf_model, predictor_names)
  
  painel <- data.frame(
    Variable = predictor_names,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(vimp_df,  by = "Variable") %>%
    dplyr::left_join(depth_df, by = "Variable") %>%
    dplyr::left_join(split_df, by = "Variable")
  
  if (!"VIMP" %in% names(painel)) painel$VIMP <- NA_real_
  if (!"minimal_depth" %in% names(painel)) painel$minimal_depth <- NA_real_
  if (!"split_frequency" %in% names(painel)) painel$split_frequency <- NA_real_
  
  thr_vimp  <- 0
  thr_depth <- if (all(is.na(painel$minimal_depth))) NA_real_ else mean(painel$minimal_depth, na.rm = TRUE)
  thr_split <- if (all(is.na(painel$split_frequency))) NA_real_ else as.numeric(stats::quantile(painel$split_frequency, 0.75, na.rm = TRUE, names = FALSE))
  
  painel <- painel %>%
    dplyr::mutate(
      crit_vimp  = !is.na(VIMP) & VIMP > thr_vimp,
      crit_depth = !is.na(minimal_depth) & !is.na(thr_depth) & minimal_depth < thr_depth,
      crit_split = !is.na(split_frequency) & !is.na(thr_split) & split_frequency > thr_split,
      score_importancia = as.integer(crit_vimp) + 
        ifelse(is.na(crit_depth), 0, as.integer(crit_depth)) + 
        ifelse(is.na(crit_split), 0, as.integer(crit_split)),
      importante_consenso = ifelse(all(is.na(minimal_depth)) & all(is.na(split_frequency)),
                                   score_importancia >= 1, 
                                   score_importancia >= 2)
    )
  
  painel <- painel %>%
    dplyr::mutate(
      rank_vimp  = ifelse(!is.na(VIMP), rank(-VIMP, ties.method = "min"), NA_real_),
      rank_depth = ifelse(!is.na(minimal_depth), rank(minimal_depth, ties.method = "min"), NA_real_),
      rank_split = ifelse(!is.na(split_frequency), rank(-split_frequency, ties.method = "min"), NA_real_)
    )
  
  painel <- painel %>%
    dplyr::mutate(
      rank_consenso = rowMeans(cbind(rank_vimp, rank_depth, rank_split), na.rm = TRUE)
    )
  
  painel$rank_consenso[is.nan(painel$rank_consenso)] <- NA_real_
  
  painel <- painel %>%
    dplyr::arrange(
      dplyr::desc(importante_consenso),
      dplyr::desc(score_importancia),
      rank_consenso,
      dplyr::desc(VIMP)
    )
  
  attr(painel, "threshold_vimp")  <- thr_vimp
  attr(painel, "threshold_depth") <- thr_depth
  attr(painel, "threshold_split") <- thr_split
  
  painel
}

selecionar_variaveis_por_consenso <- function(importance_panel, min_vars = 2) {
  if (is.null(importance_panel) || nrow(importance_panel) == 0) {
    return(character(0))
  }
  
  vars_consenso <- importance_panel %>%
    dplyr::filter(importante_consenso) %>%
    dplyr::pull(Variable)
  
  vars_consenso <- unique(as.character(vars_consenso))
  vars_consenso <- vars_consenso[!is.na(vars_consenso) & vars_consenso != ""]
  
  if (length(vars_consenso) < min_vars) {
    vars_consenso <- importance_panel %>%
      dplyr::arrange(
        dplyr::desc(score_importancia),
        rank_consenso,
        dplyr::desc(VIMP)
      ) %>%
      dplyr::slice_head(n = min_vars) %>%
      dplyr::pull(Variable)
    
    vars_consenso <- unique(as.character(vars_consenso))
    vars_consenso <- vars_consenso[!is.na(vars_consenso) & vars_consenso != ""]
  }
  
  vars_consenso
}

# ======================================================================================
# 4) FUNÇÕES PARA C-INDEX E AUC TEMPORAL
# ======================================================================================
calcular_cindex <- function(rsf_model, df_test_clean, pred_test) {
  
  cindex_oob <- tryCatch({
    err_vec <- rsf_model$err.rate
    if (is.null(err_vec) || length(err_vec) == 0) {
      NA_real_
    } else {
      1 - err_vec[length(err_vec)]
    }
  }, error = function(e) {
    cat("⚠️ Erro ao extrair C-index OOB:", e$message, "\n")
    NA_real_
  })
  
  if (nrow(df_test_clean) < 2 || sum(df_test_clean$status == 1, na.rm = TRUE) < 2) {
    return(list(
      cindex_oob = cindex_oob,
      cindex_test = NA_real_,
      cindex_test_sem_reverse = NA_real_,
      cindex_test_com_reverse = NA_real_,
      direcao_usada = NA_character_
    ))
  }
  
  marcador <- tryCatch({
    as.numeric(pred_test$predicted)
  }, error = function(e) {
    cat("⚠️ Erro ao extrair predição do teste:", e$message, "\n")
    rep(NA_real_, nrow(df_test_clean))
  })
  
  ok <- complete.cases(df_test_clean$time, df_test_clean$status, marcador)
  if (sum(ok) < 2) {
    return(list(
      cindex_oob = cindex_oob,
      cindex_test = NA_real_,
      cindex_test_sem_reverse = NA_real_,
      cindex_test_com_reverse = NA_real_,
      direcao_usada = NA_character_
    ))
  }
  
  time_ok   <- df_test_clean$time[ok]
  status_ok <- df_test_clean$status[ok]
  marcador  <- marcador[ok]
  
  if (length(unique(marcador)) < 2) {
    return(list(
      cindex_oob = cindex_oob,
      cindex_test = NA_real_,
      cindex_test_sem_reverse = NA_real_,
      cindex_test_com_reverse = NA_real_,
      direcao_usada = NA_character_
    ))
  }
  
  surv_obj <- survival::Surv(time_ok, status_ok)
  
  cindex_sem_reverse <- tryCatch({
    survival::concordance(surv_obj ~ marcador)$concordance
  }, error = function(e) {
    cat("⚠️ Erro no cálculo do C-index teste (sem reverse):", e$message, "\n")
    NA_real_
  })
  
  cindex_com_reverse <- tryCatch({
    survival::concordance(surv_obj ~ marcador, reverse = TRUE)$concordance
  }, error = function(e) {
    cat("⚠️ Erro no cálculo do C-index teste (com reverse):", e$message, "\n")
    NA_real_
  })
  
  candidatos <- c(cindex_sem_reverse, cindex_com_reverse)
  nomes_dir  <- c("sem_reverse", "com_reverse")
  
  if (all(is.na(candidatos))) {
    cindex_test <- NA_real_
    direcao_usada <- NA_character_
  } else {
    idx_best <- which.max(candidatos)
    cindex_test <- candidatos[idx_best]
    direcao_usada <- nomes_dir[idx_best]
  }
  
  return(list(
    cindex_oob = cindex_oob,
    cindex_test = cindex_test,
    cindex_test_sem_reverse = cindex_sem_reverse,
    cindex_test_com_reverse = cindex_com_reverse,
    direcao_usada = direcao_usada
  ))
}

calcular_auc_temporal <- function(df_test, marcador_risco_test, horizons_years = c(1, 3, 5, 10)) {
  horizons_days <- horizons_years * 365
  
  auc_results <- data.frame(
    Horizonte_anos = horizons_years,
    AUC = NA_real_,
    N_eventos = NA_integer_,
    N_risco = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  if (length(unique(na.omit(marcador_risco_test))) < 2) {
    cat("⚠️ Marcador de risco constante - não é possível calcular AUC\n")
    return(auc_results)
  }
  
  for (k in seq_along(horizons_days)) {
    t_cut <- horizons_days[k]
    
    eventos_ate_t <- sum(df_test$time <= t_cut & df_test$status == 1, na.rm = TRUE)
    auc_results$N_eventos[k] <- eventos_ate_t
    
    em_risco_t <- sum(df_test$time > t_cut, na.rm = TRUE)
    auc_results$N_risco[k] <- em_risco_t
    
    if (eventos_ate_t < 5) {
      # Usando cat diretamente em vez de log_skip com string problemática
      cat("⚠️ Eventos insuficientes para horizonte", horizons_years[k], "anos:", eventos_ate_t, "\n")
      next
    }
    
    roc_t <- tryCatch(
      timeROC(
        T = df_test$time,
        delta = df_test$status,
        marker = marcador_risco_test,
        cause = 1,
        weighting = "marginal",
        times = t_cut,
        iid = TRUE
      ),
      error = function(e) {
        cat("⚠️ Erro no timeROC para horizonte", horizons_years[k], "anos:", e$message, "\n")
        NULL
      }
    )
    
    if (!is.null(roc_t) && !is.null(roc_t$AUC) && length(roc_t$AUC) >= 1) {
      auc_results$AUC[k] <- tail(roc_t$AUC, 1)
    }
  }
  
  return(auc_results)
}

criar_tabela_metricas <- function(cindex_oob, cindex_test, auc_results) {
  metricas_df <- data.frame(
    Metrica = c("C-index OOB", "C-index Teste",
                paste0("AUC ", auc_results$Horizonte_anos, " ano(s)")),
    Valor = c(cindex_oob, cindex_test, auc_results$AUC),
    Interpretacao = c(
      ifelse(is.na(cindex_oob), "N/A",
             ifelse(cindex_oob >= 0.8, "Excelente",
                    ifelse(cindex_oob >= 0.7, "Bom",
                           ifelse(cindex_oob >= 0.6, "Aceitável", "Ruim")))),
      ifelse(is.na(cindex_test), "N/A",
             ifelse(cindex_test >= 0.8, "Excelente",
                    ifelse(cindex_test >= 0.7, "Bom",
                           ifelse(cindex_test >= 0.6, "Aceitável", "Ruim")))),
      sapply(auc_results$AUC, function(x) {
        ifelse(is.na(x), "N/A",
               ifelse(x >= 0.8, "Excelente",
                      ifelse(x >= 0.7, "Bom",
                             ifelse(x >= 0.6, "Aceitável", "Ruim"))))
      })
    ),
    N_eventos = c(NA, NA, auc_results$N_eventos),
    Em_risco = c(NA, NA, auc_results$N_risco),
    stringsAsFactors = FALSE
  )
  
  if (!is.na(cindex_oob) && !is.na(cindex_test)) {
    diff_overfitting <- abs(cindex_oob - cindex_test)
    overfitting_status <- ifelse(diff_overfitting < 0.05, "Baixo",
                                 ifelse(diff_overfitting < 0.10, "Moderado", "Alto"))
  } else {
    diff_overfitting <- NA_real_
    overfitting_status <- "N/A"
  }
  
  diagnostico_df <- data.frame(
    Metrica = "Diferença OOB-Teste",
    Valor = diff_overfitting,
    Interpretacao = overfitting_status,
    N_eventos = NA,
    Em_risco = NA,
    stringsAsFactors = FALSE
  )
  
  tabela_final <- rbind(metricas_df, diagnostico_df)
  return(tabela_final)
}

plotar_comparacao_metricas <- function(cindex_oob, cindex_test, auc_results,
                                       tipo_cancer, metrica, omic_label, pheno_label,
                                       output_dir) {
  df_plot <- data.frame(
    Metrica = c("C-index OOB", "C-index Teste",
                paste0("AUC ", auc_results$Horizonte_anos, "y")),
    Valor = c(cindex_oob, cindex_test, auc_results$AUC),
    Tipo = c(rep("C-index", 2), rep("AUC Temporal", nrow(auc_results))),
    stringsAsFactors = FALSE
  )
  
  df_plot <- df_plot[!is.na(df_plot$Valor), ]
  
  if (nrow(df_plot) == 0) {
    return(NULL)
  }
  
  p <- ggplot(df_plot, aes(x = Metrica, y = Valor, fill = Tipo)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "orange", alpha = 0.5) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "green", alpha = 0.5) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    labs(
      title = paste("Comparação de Métricas -", tipo_cancer, "-", metrica),
      subtitle = paste("Camada:", omic_label, "| Fenótipo:", pheno_label),
      y = "Valor",
      x = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    scale_fill_manual(values = c("C-index" = "#3498db", "AUC Temporal" = "#e74c3c")) +
    geom_text(aes(label = sprintf("%.3f", Valor)), vjust = -0.5, size = 3.5)
  
  ggsave(
    filename = file.path(output_dir,
                         paste0("Comparacao_Metricas_", tipo_cancer, "_", metrica, "_",
                                omic_label, "_", pheno_label, ".pdf")),
    plot = p,
    width = 10,
    height = 6
  )
  
  return(p)
}

plotar_auc_temporal <- function(auc_results, tipo_cancer, metrica,
                                omic_label, pheno_label, output_dir) {
  auc_valido <- auc_results[!is.na(auc_results$AUC), ]
  
  if (nrow(auc_valido) == 0) {
    return(NULL)
  }
  
  p <- ggplot(auc_valido, aes(x = Horizonte_anos, y = AUC)) +
    geom_line(color = "#e74c3c", linewidth = 1.2) +
    geom_point(aes(size = N_eventos), color = "#c0392b", alpha = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "orange", alpha = 0.5) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "green", alpha = 0.5) +
    scale_x_continuous(breaks = auc_valido$Horizonte_anos) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_size_continuous(name = "Nº Eventos", range = c(3, 8)) +
    labs(
      title = paste("AUC Temporal -", tipo_cancer, "-", metrica),
      subtitle = paste("Camada:", omic_label, "| Fenótipo:", pheno_label),
      x = "Horizonte (anos)",
      y = "AUC"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    ) +
    geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -1, size = 3.5)
  
  ggsave(
    filename = file.path(output_dir,
                         paste0("AUC_Temporal_", tipo_cancer, "_", metrica, "_",
                                omic_label, "_", pheno_label, ".pdf")),
    plot = p,
    width = 8,
    height = 6
  )
  
  return(p)
}

# ======================================================================================
# 5) CARREGAR PROGRESSO GLOBAL
# ======================================================================================
if (file.exists(progress_file)) {
  progress_obj <- tryCatch(
    readRDS(progress_file),
    error = function(e) {
      stop("❌ Não foi possível ler o arquivo global de progresso: ", e$message)
    }
  )
  
  resultados_geral <- progress_obj$resultados_geral
  progress_keys    <- progress_obj$progress_keys
  
  if (is.null(resultados_geral)) resultados_geral <- data.frame()
  if (is.null(progress_keys)) progress_keys <- character()
  
  if (nrow(resultados_geral) > 0) {
    needed_cols <- c("Tipo", "Metrica", "Omic_layer", "Pheno_layer")
    for (cc in needed_cols) {
      if (!cc %in% colnames(resultados_geral)) {
        resultados_geral[[cc]] <- NA_character_
      }
    }
  }
  
  cat("🔁 Retomando execução a partir do arquivo global:\n")
  cat("   ", progress_file, "\n")
  cat("✅ Combinações já processadas:", length(progress_keys), "\n")
  
} else {
  resultados_geral <- data.frame()
  progress_keys <- character()
  
  save_global_progress(
    progress_file = progress_file,
    resultados_geral = resultados_geral,
    progress_keys = progress_keys
  )
  
  cat("🆕 Arquivo global de progresso criado em:\n")
  cat("   ", progress_file, "\n")
}

# ======================================================================================
# 6) DEFINIÇÕES GERAIS
# ======================================================================================
audit_log_v4 <- list()
log_skip <- function(modo_vars, omic_vars, tipo_cancer, metrica, ...) {
  msg <- paste(...)
  cat("⚠️", msg, "\n")
  audit_log_v4[[length(audit_log_v4) + 1]] <<- data.frame(
    Modo = modo_vars, Omic_Layer = omic_vars, Cancer = tipo_cancer, Endpoint = metrica,
    Status = "SKIPPED", Razao = msg, stringsAsFactors = FALSE
  )
}
metricas <- c("OS", "DSS", "DFI", "PFI")
tipos_tumorais <- sort(unique(coverage_map$cancer_type))

# ======================================================================================
# 7) LOOP EXTERNO PARA TODOS OS MODELOS
# ======================================================================================
for (modo_vars in 1:4) {
  for (omic_vars in 1:2) {
    
    cat("\n======================================\n")
    cat("🔧 Rodando modelo com:\n")
    cat("   → modo_vars     =", modo_vars, "\n")
    cat("   → omic_vars     =", omic_vars, "\n")
    cat("   → fenotipo      = multifenotipo fixo\n")
    cat("======================================\n")
    
    omic_root  <- ifelse(omic_vars == 1, "multiomica", "mono_omica")
    pheno_root <- "multifenotipo"
    
    base_output_dir <- file.path(
      getwd(), "output", "modelo",
      paste0("modelo", modo_vars),
      omic_root,
      pheno_root
    )
    if (!dir.exists(base_output_dir)) dir.create(base_output_dir, recursive = TRUE)
    
    for (tipo_cancer in tipos_tumorais) {
      cat("\n\n🧪 Tipo tumoral:", tipo_cancer, "\n")
      
      for (metrica in metricas) {
        cat("\n======================\n🧪 Iniciando métrica (sobrevida):", metrica, "\n")
        
        df_file <- get_df_file_for_pair(
          tipo_cancer = tipo_cancer,
          metrica     = metrica,
          coverage_map = coverage_map,
          base_dir = getwd()
        )
        
        if (is.na(df_file)) {
          log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Nenhum assigned_df encontrado em coverage_map para", tipo_cancer, "+", metrica, " - pulando")
          next
        }
        
        if (!file.exists(df_file)) {
          log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Arquivo RDS não encontrado para", tipo_cancer, "+", metrica, ":", df_file, " - pulando")
          next
        }
        
        cat("📦 Usando dataframe:", basename(df_file), "para", tipo_cancer, "+", metrica, "\n")
        
        df12 <- tryCatch(
          load_df_cached(df_file),
          error = function(e) {
            cat("❌ Erro ao carregar", df_file, ":", e$message, "\n")
            NULL
          }
        )
        if (is.null(df12)) next
        
        if (!"type" %in% colnames(df12)) {
          log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "O dataframe", basename(df_file), "não possui coluna 'type' - pulando")
          next
        }
        
        df12_filter <- df12 %>% dplyr::filter(type == tipo_cancer)
        if (nrow(df12_filter) == 0) {
          log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "O dataframe", basename(df_file), "não possui amostras para", tipo_cancer, " - pulando")
          next
        }
        
        time_var <- paste0(metrica, ".time")
        padrao_regex_all <- paste0("^", tipo_cancer, ".*\\.")
        
        df12_filter_RF <- df12_filter %>%
          dplyr::select(
            1:37,
            dplyr::any_of(c(metrica, time_var)),
            dplyr::matches(padrao_regex_all)
          )
        
        if (!all(c(metrica, time_var) %in% colnames(df12_filter_RF))) {
          log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Endpoint", metrica, "ou", time_var, "não encontrado em", basename(df_file), " - pulando")
          next
        }
        
        n_omic_all <- sum(grepl(padrao_regex_all, names(df12_filter_RF)))
        if (n_omic_all == 0) {
          log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Nenhuma coluna ômica correspondente ao padrão", padrao_regex_all, "em", basename(df_file), "")
          next
        }
        
        names(df12_filter_RF) <- gsub("-", "_", names(df12_filter_RF))
        names(df12_filter_RF) <- make.names(names(df12_filter_RF), unique = TRUE)
        if (any(duplicated(names(df12_filter_RF)))) {
          df12_filter_RF <- df12_filter_RF[, !duplicated(names(df12_filter_RF))]
        }
        
        df_temp <- df12_filter_RF %>%
          dplyr::filter(!is.na(.data[[metrica]]) & !is.na(.data[[time_var]])) %>%
          dplyr::select(-any_of(c("status", "time")))
        
        df_temp <- df_temp %>%
          dplyr::mutate(
            status = ifelse(as.numeric(as.character(.data[[metrica]])) == 1, 1, 0),
            time   = as.numeric(.data[[time_var]])
          )
        
        if ("sample" %in% names(df_temp))  df_temp$sample  <- as.character(df_temp$sample)
        if ("patient" %in% names(df_temp)) df_temp$patient <- as.character(df_temp$patient)
        
        if (nrow(df_temp) < 20) {
          log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Poucas amostras com endpoint disponível para", tipo_cancer, metrica, " - pulando")
          next
        }
        
        if (length(unique(df_temp$status)) < 2) {
          log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Apenas uma classe no desfecho", metrica, " - pulando")
          next
        }
        
        sample_ids_full <- dplyr::case_when(
          "sample"  %in% names(df_temp) ~ as.character(df_temp$sample),
          "patient" %in% names(df_temp) ~ as.character(df_temp$patient),
          TRUE ~ as.character(seq_len(nrow(df_temp)))
        )
        stopifnot(length(sample_ids_full) == nrow(df_temp))
        
        set.seed(123)
        idx_train <- caret::createDataPartition(df_temp$status, p = 0.7, list = FALSE)
        df_train0  <- df_temp[idx_train, , drop = FALSE]
        df_test0   <- df_temp[-idx_train, , drop = FALSE]
        
        sample_ids_train0 <- sample_ids_full[idx_train]
        sample_ids_test0  <- sample_ids_full[-idx_train]
        
        if (length(unique(df_train0$status)) < 2) {
          log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Apenas uma classe no TREINO para", metrica, " - pulando")
          next
        }
        
        clinica_train0 <- df_train0 %>% dplyr::select(1:37)
        clinica_test0  <- df_test0  %>% dplyr::select(1:37)
        
        if (omic_vars == 1) {
          camada_list <- list(list(label = "all",  padrao = paste0("^", tipo_cancer, ".*\\.")))
        } else {
          camada_list <- list(
            list(label = "omic5", padrao = paste0("^", tipo_cancer, "_.*\\.5\\.")),
            list(label = "omic6", padrao = paste0("^", tipo_cancer, "_.*\\.6\\."))
          )
        }
        
        for (cam in camada_list) {
          omic_label  <- cam$label
          padrao_omic <- cam$padrao
          pheno_label <- "allPheno"
          
          dir_output_layer <- file.path(base_output_dir, tipo_cancer, metrica, omic_label)
          if (!dir.exists(dir_output_layer)) dir.create(dir_output_layer, recursive = TRUE)
          
          genes_train_layer <- df_train0 %>% dplyr::select(dplyr::matches(padrao_omic))
          genes_test_layer  <- df_test0  %>% dplyr::select(dplyr::matches(padrao_omic))
          
          if (ncol(genes_train_layer) < 2) {
            log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Camada", omic_label, "sem colunas suficientes - pulando")
            next
          }
          
          key <- make_progress_key(
            modo_vars = modo_vars,
            omic_vars = omic_vars,
            tipo = tipo_cancer,
            metrica = metrica,
            omic_label = omic_label,
            pheno_label = pheno_label
          )
          
          if (key %in% progress_keys) {
            cat("✅ Já processado:", key, " - pulando...\n")
            next
          }
          
          cat("➡️  Camada ômica:", omic_label, " | Fenótipo:", pheno_label, "\n")
          
          df_train_base <- df_train0
          df_test_base  <- df_test0
          clinica_train_base <- clinica_train0
          clinica_test_base  <- clinica_test0
          sample_ids_train_base <- sample_ids_train0
          sample_ids_test_base  <- sample_ids_test0
          
          # Profilaxia 1: Remover ômicas com MUITOS NAs (> 20%) antes do complete.cases
          na_prop_train <- colMeans(is.na(genes_train_layer))
          keep_cols_na <- na_prop_train <= 0.20
          if (sum(!keep_cols_na) > 0) {
            genes_train_layer <- genes_train_layer[, keep_cols_na, drop = FALSE]
            genes_test_layer  <- genes_test_layer[, keep_cols_na, drop = FALSE]
          }
          
          # Profilaxia 2: Remover NZV antes do modelo e Boruta
          nzv_pre <- caret::nearZeroVar(genes_train_layer)
          if (length(nzv_pre) > 0) {
            genes_train_layer <- genes_train_layer[, -nzv_pre, drop = FALSE]
            genes_test_layer  <- genes_test_layer[, -nzv_pre, drop = FALSE]
          }
          
          keep_tr <- complete.cases(genes_train_layer)
          genes_train <- genes_train_layer[keep_tr, , drop = FALSE]
          df_train    <- df_train_base[keep_tr, , drop = FALSE]
          clinica_train_full <- clinica_train_base[keep_tr, , drop = FALSE]
          sample_ids_train <- sample_ids_train_base[keep_tr]
          
          keep_te <- complete.cases(genes_test_layer)
          genes_test <- genes_test_layer[keep_te, , drop = FALSE]
          df_test    <- df_test_base[keep_te, , drop = FALSE]
          clinica_test_full <- clinica_test_base[keep_te, , drop = FALSE]
          sample_ids_test <- sample_ids_test_base[keep_te]
          
          if (length(sample_ids_train) != nrow(df_train)) {
            sample_ids_train <- if ("sample" %in% names(df_train)) as.character(df_train$sample) else
              if ("patient" %in% names(df_train)) as.character(df_train$patient) else
                as.character(seq_len(nrow(df_train)))
          }
          if (length(sample_ids_test) != nrow(df_test)) {
            sample_ids_test <- if ("sample" %in% names(df_test)) as.character(df_test$sample) else
              if ("patient" %in% names(df_test)) as.character(df_test$patient) else
                as.character(seq_len(nrow(df_test)))
          }
          
          stopifnot(length(sample_ids_train) == nrow(df_train))
          stopifnot(length(sample_ids_test)  == nrow(df_test))
          
          # === CORREÇÃO 1: Verificação de nrow(df_train) após remoção de NAs ===
          if (nrow(df_train) < 20) {
            cat("⚠️ Poucas amostras no TREINO após remoção NA ômico - pulando\n")
            next
          }
          if (length(unique(df_train$status)) < 2) {
            cat("⚠️ Apenas uma classe no TREINO após remoção NA ômico - pulando\n")
            next
          }
          if (nrow(df_test) < 10) {
            cat("⚠️ Poucas amostras no TESTE após remoção NA ômico - métricas podem falhar\n")
          }
          
          demog_vars_pheno <- character()
          
          # --------------------------------------------------------------------------
          # SELEÇÃO DE VARIÁVEIS
          # --------------------------------------------------------------------------
          importance_panel_pre <- data.frame()
          importance_consenso_pre <- data.frame()
          selected_vars_pre <- character()
          
          if (modo_vars %in% c(1, 2)) {
            genes_boruta <- genes_train[, colSums(is.na(genes_train)) == 0, drop = FALSE]
            variaveis_usadas_path <- file.path(dir_output_layer, paste0("Variaveis_Usadas_No_Modelo_", tipo_cancer, "_", metrica, "_", omic_label, "_", pheno_label, ".csv"))
            
            if (file.exists(variaveis_usadas_path)) {
              cat("  ♻️ Recuperando variáveis do arquivo existente (pulando Boruta!)...\n")
              df_vars <- read.csv(variaveis_usadas_path)
              selected_vars <- df_vars$Variable[df_vars$Origem_selecao == "Boruta"]
              if (length(selected_vars) == 0) selected_vars <- df_vars$Variable[df_vars$Tipo_variavel == "omica"]
            } else {
              if (ncol(genes_boruta) < 2) { cat("⚠️ Poucas colunas sem NA para Boruta - pulando\n"); next }
              
              cc_boruta <- complete.cases(genes_boruta, df_train$status)
              if (sum(cc_boruta) < 20) { cat("⚠️ Amostras completas insuficientes para Boruta - pulando\n"); next }
              
              set.seed(123)
              boruta_result <- tryCatch({
                Boruta(
                  x = genes_boruta[cc_boruta, , drop = FALSE],
                  y = factor(df_train$status[cc_boruta], labels = c("No", "Yes")),
                  doTrace = 0, maxRuns = 100
                )
              }, error = function(e) { cat("❌ Erro no Boruta:", e$message, "\n"); NULL })
              
              if (is.null(boruta_result)) next
              
              selected_vars <- getSelectedAttributes(boruta_result, withTentative = FALSE)
            }
            if (length(selected_vars) < 2) { cat("⚠️ Poucos genes selecionados por Boruta - pulando\n"); next }
            
            cat("Modulo", modo_vars, "- usando", length(selected_vars), "genes selecionados por Boruta\n")
            
          } else {
            genes_sem_na <- genes_train[, colSums(is.na(genes_train)) == 0, drop = FALSE]
            if (ncol(genes_sem_na) < 2) { cat("⚠️ Poucas colunas sem NA para modo", modo_vars, "- pulando\n"); next }
            
            df_pre_train <- data.frame(
              time   = df_train$time,
              status = df_train$status,
              genes_sem_na,
              check.names = FALSE
            )
            
            cc_pre <- complete.cases(df_pre_train)
            df_pre_train <- df_pre_train[cc_pre, , drop = FALSE]
            
            if (nrow(df_pre_train) < 20) {
              log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Amostras insuficientes para RSF preliminar nos modos", modo_vars, " - pulando")
              next
            }
            
            if (length(unique(df_pre_train$status)) < 2) {
              log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Apenas uma classe no RSF preliminar dos modos", modo_vars, " - pulando")
              next
            }
            
            set.seed(123)
            rsf_preliminar <- tryCatch({
              rfsrc(
                formula    = Surv(time, status) ~ .,
                data       = df_pre_train,
                ntree      = 2000,
                importance = "permute",
                na.action  = "na.omit"
              )
            }, error = function(e) {
              cat("❌ Erro ao treinar RSF preliminar:", e$message, "\n")
              NULL
            })
            
            if (is.null(rsf_preliminar)) next
            
            importance_panel_pre <- tryCatch({
              construir_painel_importancia(
                rsf_model = rsf_preliminar,
                predictor_names = colnames(df_pre_train)[!colnames(df_pre_train) %in% c("time", "status")]
              )
            }, error = function(e) {
              cat("⚠️ Erro ao construir painel preliminar de importância:", e$message, "\n")
              data.frame()
            })
            
            if (nrow(importance_panel_pre) == 0) {
              log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Painel preliminar de importância vazio - pulando")
              next
            }
            
            selected_vars_pre <- selecionar_variaveis_por_consenso(
              importance_panel = importance_panel_pre,
              min_vars = 2
            )
            
            if (length(selected_vars_pre) < 2) {
              log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Consenso selecionou menos de 2 variáveis - pulando")
              next
            }
            
            selected_vars <- intersect(selected_vars_pre, colnames(genes_sem_na))
            
            if (length(selected_vars) < 2) {
              log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "Variáveis selecionadas por consenso não estão disponíveis no treino final - pulando")
              next
            }
            
            importance_consenso_pre <- importance_panel_pre %>%
              dplyr::filter(Variable %in% selected_vars)
            
            cat("Modulo", modo_vars, "(sem Boruta): RSF preliminar selecionou",
                length(selected_vars), "variáveis por consenso para o RSF final\n")
          }
          
          # --------------------------------------------------------------------------
          # VARIÁVEIS FENOTÍPICAS EXTRAS DIRIGIDAS PELAS ASSINATURAS
          # --------------------------------------------------------------------------
          demog_vars_pheno <- obter_variaveis_fenotipicas_extras(
            selected_vars = selected_vars,
            posicao_fenotipo = FENOTIPO_COMPONENT_POS
          )
          
          phenos_presentes <- obter_fenotipos_presentes(
            selected_vars = selected_vars,
            posicao_fenotipo = FENOTIPO_COMPONENT_POS
          )
          
          demog_vars_pheno <- intersect(demog_vars_pheno, colnames(df_train))
          
          if (length(phenos_presentes) > 0) {
            cat("📌 Fenótipos detectados nas assinaturas selecionadas:",
                paste(phenos_presentes, collapse = ", "), "\n")
          } else {
            cat("📌 Nenhum fenótipo 1/2/3 detectado nas assinaturas selecionadas.\n")
          }
          
          if (length(demog_vars_pheno) > 0) {
            cat("📌 Variáveis fenotípicas extras incluídas:",
                paste(demog_vars_pheno, collapse = ", "), "\n")
          } else {
            cat("📌 Nenhuma variável fenotípica extra disponível para inclusão.\n")
          }
          
          # --------------------------------------------------------------------------
          # SELEÇÃO DE CLÍNICAS
          # --------------------------------------------------------------------------
          clinicas_relevantes <- c()
          pval_df <- data.frame()
          
          if (modo_vars %in% c(1, 3)) {
            
            clinicas_possiveis <- c(
              "age_at_initial_pathologic_diagnosis","gender","race",
              "ajcc_pathologic_tumor_stage","clinical_stage","histological_type",
              "histological_grade","menopause_status",
              "new_tumor_event_type","margin_status",
              "residual_tumor","initial_pathologic_dx_year",
              "new_tumor_event_site","new_tumor_event_site_other",
              "treatment_outcome_first_course","BestCall"
            ) #"tumor_status",
            
            clinicas_possiveis <- unique(c(clinicas_possiveis, demog_vars_pheno))
            clinicas_disponiveis <- intersect(clinicas_possiveis, colnames(df_train))
            
            for (var in clinicas_disponiveis) {
              vals <- df_train[[var]]
              if (sum(!is.na(vals)) < 10) next
              
              if (is.factor(vals) || is.character(vals)) {
                vals <- as.factor(vals)
                if (nlevels(vals) < 2) next
                tab_vals <- table(vals)
                if (any(tab_vals < 5)) next
              } else {
                if (sd(vals, na.rm = TRUE) == 0) next
              }
              
              df_tmp_var <- df_train[complete.cases(df_train[, c("time", "status", var)]), , drop = FALSE]
              if (nrow(df_tmp_var) < 10) next
              
              formula_uni <- as.formula(paste0("Surv(time, status) ~ `", var, "`"))
              cox_model <- tryCatch(
                withCallingHandlers(
                  coxph(formula_uni, data = df_tmp_var, control = coxph.control(iter.max = 50)),
                  warning = function(w) {
                    message("⚠️ Warning em coxph para variável ", var, ": ", conditionMessage(w))
                    invokeRestart("muffleWarning")
                  }
                ),
                error = function(e) {
                  message("⚠️ Erro em coxph para variável ", var, ": ", e$message)
                  NULL
                }
              )
              
              if (!is.null(cox_model)) {
                coef_sum <- summary(cox_model)$coefficients
                if (!is.null(coef_sum) && nrow(as.matrix(coef_sum)) > 0) {
                  coef_sum <- as.matrix(coef_sum)
                  pvals    <- coef_sum[, "Pr(>|z|)"]
                  termos   <- rownames(coef_sum)
                  
                  pval_df <- rbind(pval_df, data.frame(
                    Variavel = rep(var, length(pvals)),
                    Termo    = termos,
                    Pvalor   = pvals,
                    stringsAsFactors = FALSE
                  ))
                  
                  if (any(pvals < 0.05, na.rm = TRUE)) {
                    clinicas_relevantes <- c(clinicas_relevantes, var)
                  }
                }
              }
            }
            
            clinicas_relevantes <- unique(clinicas_relevantes)
            
            write.csv(
              pval_df,
              file = file.path(dir_output_layer,
                               paste0("Pvalores_Cox_", tipo_cancer, "_", metrica, "_", omic_label, "_", pheno_label, ".csv")),
              row.names = FALSE
            )
            
            if (length(clinicas_relevantes) == 0) {
              cat("⚠️ Nenhuma variável clínica/demográfica significativa (TREINO) - Continuando apenas com ômicas...\n")
            }
            
          } else {
            cat("Modulo", modo_vars, "(somente ômicas): sem clínicas/demográficas\n")
          }
          
          if (modo_vars %in% c(1, 3) && length(clinicas_relevantes) > 0) {
            clinicas_train_sel <- clinica_train_full %>% dplyr::select(all_of(clinicas_relevantes))
            clinicas_test_sel  <- clinica_test_full  %>% dplyr::select(all_of(clinicas_relevantes))
          } else {
            clinicas_train_sel <- NULL
            clinicas_test_sel  <- NULL
          }
          
          sel_train <- intersect(selected_vars, colnames(genes_train))
          sel_test  <- intersect(selected_vars, colnames(genes_test))
          
          if (length(sel_train) < 2) {
            cat("⚠️ selected_vars não bate com genes_train após filtros - pulando\n")
            next
          }
          if (length(sel_test)  < 2) {
            cat("⚠️ selected_vars não bate com genes_test após filtros - pulando\n")
            next
          }
          
          genes_train_sel <- as.data.frame(genes_train[, sel_train, drop = FALSE])
          genes_test_sel  <- as.data.frame(genes_test[,  sel_test,  drop = FALSE])
          
          stopifnot(nrow(df_train) == nrow(genes_train_sel))
          stopifnot(nrow(df_test)  == nrow(genes_test_sel))
          stopifnot(length(sample_ids_train) == nrow(df_train))
          stopifnot(length(sample_ids_test)  == nrow(df_test))
          
          df_model_train <- data.frame(
            sample_id = sample_ids_train,
            time      = df_train$time,
            status    = df_train$status,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
          df_model_test <- data.frame(
            sample_id = sample_ids_test,
            time      = df_test$time,
            status    = df_test$status,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
          
          if (!is.null(clinicas_train_sel)) {
            stopifnot(nrow(clinicas_train_sel) == nrow(df_train))
            df_model_train <- cbind(df_model_train, as.data.frame(clinicas_train_sel))
          }
          if (!is.null(clinicas_test_sel)) {
            stopifnot(nrow(clinicas_test_sel) == nrow(df_test))
            df_model_test <- cbind(df_model_test, as.data.frame(clinicas_test_sel))
          }
          
          df_model_train <- cbind(df_model_train, genes_train_sel)
          df_model_test  <- cbind(df_model_test,  genes_test_sel)
          
          df_model_train <- df_model_train %>%
            dplyr::mutate(across(where(is.logical), as.factor),
                          across(where(is.character), as.factor))
          df_model_test <- df_model_test %>%
            dplyr::mutate(across(where(is.logical), as.factor),
                          across(where(is.character), as.factor))
          
          X_train <- df_model_train %>% dplyr::select(-sample_id, -time, -status)
          X_train <- X_train[, colSums(!is.na(X_train)) > 0, drop = FALSE]
          if (ncol(X_train) == 0) { cat("⚠️ Sem preditores não-NA (TREINO) - pulando\n"); next }
          
          nzv <- caret::nearZeroVar(X_train)
          if (length(nzv) > 0 && length(nzv) < ncol(X_train)) {
            X_train <- X_train[, -nzv, drop = FALSE]
          }
          if (ncol(X_train) == 0) { cat("⚠️ Sem preditores após NZV (TREINO) - pulando\n"); next }
          
          common_cols <- colnames(X_train)
          
          X_test_raw <- df_model_test %>% dplyr::select(-sample_id, -time, -status)
          missing_in_test <- setdiff(common_cols, colnames(X_test_raw))
          if (length(missing_in_test) > 0) {
            for (mc in missing_in_test) X_test_raw[[mc]] <- NA
          }
          X_test <- X_test_raw %>% dplyr::select(all_of(common_cols))
          
          df_model_train <- cbind(
            sample_id = df_model_train$sample_id,
            time      = df_model_train$time,
            status    = df_model_train$status,
            X_train
          )
          df_model_test <- cbind(
            sample_id = df_model_test$sample_id,
            time      = df_model_test$time,
            status    = df_model_test$status,
            X_test
          )
          
          cc_train <- complete.cases(df_model_train[, -1])
          df_model_train <- df_model_train[cc_train, , drop = FALSE]
          
          # === CORREÇÃO 2: Verificação de nrow(df_model_train) ===
          if (nrow(df_model_train) < 20) {
            cat("⚠️ Amostras TREINO insuficientes para RSF - pulando\n")
            next
          }
          if (length(unique(df_model_train$status)) < 2) {
            cat("⚠️ Apenas uma classe no TREINO após filtros - pulando\n")
            next
          }
          
          cc_test <- complete.cases(df_model_test[, -1])
          df_model_test <- df_model_test[cc_test, , drop = FALSE]
          if (nrow(df_model_test) < 10) cat("⚠️ Amostras TESTE muito poucas para avaliação - métricas podem falhar\n")
          
          patient_ids_train <- df_model_train$sample_id
          patient_ids_test  <- df_model_test$sample_id
          
          df_train_clean <- df_model_train %>% dplyr::select(-sample_id)
          df_test_clean  <- df_model_test  %>% dplyr::select(-sample_id)
          
          # --------------------------------------------------------------------------
          # RSF FINAL
          # --------------------------------------------------------------------------
          importance_mode <- "permute"
          
          rsf_model <- tryCatch({
            rfsrc(
              formula    = Surv(time, status) ~ .,
              data       = df_train_clean,
              ntree      = 5000,
              importance = importance_mode,
              na.action  = "na.omit"
            )
          }, error = function(e) {
            cat("❌ Erro ao treinar RSF:", e$message, "\n")
            NULL
          })
          if (is.null(rsf_model)) next
          
          pred_test <- tryCatch({
            predict(rsf_model, newdata = df_test_clean, na.action = "na.omit")
          }, error = function(e) {
            cat("❌ Erro na predição do RSF:", e$message, "\n")
            NULL
          })
          if (is.null(pred_test)) next
          
          marcador_risco_test <- tryCatch({
            as.numeric(pred_test$predicted)
          }, error = function(e) {
            cat("⚠️ Não foi possível extrair marcador de risco:", e$message, "\n")
            rep(NA_real_, nrow(df_test_clean))
          })
          
          # ==========================================================================
          # SALVAR MATERIAL PARA XAI POST-HOC
          # ==========================================================================
          tryCatch({
            explainer_dir <- file.path(dir_output_layer, "Explicabilidade")
            if (!dir.exists(explainer_dir)) {
              dir.create(explainer_dir, recursive = TRUE)
            }
            
            x_test_xai <- df_test_clean %>%
              dplyr::select(-time, -status)
            
            if (nrow(x_test_xai) == 0 || ncol(x_test_xai) == 0) {
              cat("⚠️ Dados insuficientes para preparar bundle de explicabilidade.\n")
            } else {
              
              bg_size <- min(150, nrow(x_test_xai))
              set.seed(123)
              bg_idx <- sample(seq_len(nrow(x_test_xai)), bg_size)
              
              xai_bundle <- list(
                meta = list(
                  tipo_cancer = tipo_cancer,
                  metrica = metrica,
                  omic_label = omic_label,
                  pheno_label = pheno_label,
                  modo_vars = modo_vars,
                  omic_vars = omic_vars,
                  criado_em = Sys.time()
                ),
                model = rsf_model,
                patient_ids_test = patient_ids_test,
                df_test_clean = df_test_clean,
                x_test_xai = x_test_xai,
                pred_test = marcador_risco_test,
                bg_idx = bg_idx,
                bg_data = x_test_xai[bg_idx, , drop = FALSE],
                bg_surv = survival::Surv(
                  df_test_clean$time[bg_idx],
                  df_test_clean$status[bg_idx]
                )
              )
              
              saveRDS(
                xai_bundle,
                file = file.path(
                  explainer_dir,
                  paste0(
                    "XAI_bundle_",
                    tipo_cancer, "_", metrica, "_", omic_label, "_", pheno_label, ".rds"
                  )
                )
              )
              
              resumo_xai <- data.frame(
                idx = seq_len(nrow(x_test_xai)),
                patient_id = as.character(patient_ids_test),
                risk_score = as.numeric(marcador_risco_test),
                time = df_test_clean$time,
                status = df_test_clean$status,
                stringsAsFactors = FALSE
              )
              
              resumo_xai$patient_id <- ifelse(
                is.na(resumo_xai$patient_id) | resumo_xai$patient_id == "",
                paste0("Paciente_", resumo_xai$idx),
                resumo_xai$patient_id
              )
              
              write.csv(
                resumo_xai,
                file = file.path(
                  explainer_dir,
                  paste0(
                    "Resumo_Pacientes_XAI_",
                    tipo_cancer, "_", metrica, "_", omic_label, "_", pheno_label, ".csv"
                  )
                ),
                row.names = FALSE
              )
              
              cat("✅ Bundle para XAI post-hoc salvo com sucesso.\n")
            }
            
          }, error = function(e) {
            cat("⚠️ Falha ao salvar bundle de explicabilidade:", e$message, "\n")
          })
          # ==========================================================================
          
          cindex_results <- calcular_cindex(rsf_model, df_test_clean, pred_test)
          
          if (!is.list(cindex_results) || !all(c("cindex_oob", "cindex_test") %in% names(cindex_results))) {
            log_skip(modo_vars, omic_vars, tipo_cancer, metrica, "calcular_cindex() retornou objeto inválido - pulando combinação")
            next
          }
          
          cindex_oob  <- cindex_results$cindex_oob
          cindex_test <- cindex_results$cindex_test
          cindex_test_sem_reverse <- cindex_results$cindex_test_sem_reverse
          cindex_test_com_reverse <- cindex_results$cindex_test_com_reverse
          direcao_cindex <- cindex_results$direcao_usada
          
          cat("📊 C-index OOB:", round(cindex_oob, 4), "\n")
          cat("📊 C-index Teste (sem reverse):", round(cindex_test_sem_reverse, 4), "\n")
          cat("📊 C-index Teste (com reverse):", round(cindex_test_com_reverse, 4), "\n")
          cat("📊 C-index Teste final:", round(cindex_test, 4), " | direção usada:", direcao_cindex, "\n")
          
          horizons_years <- c(1, 3, 5, 10)
          auc_results <- calcular_auc_temporal(
            df_test = df_test_clean,
            marcador_risco_test = marcador_risco_test,
            horizons_years = horizons_years
          )
          
          cat("📊 AUC Temporal:\n")
          for (i in seq_len(nrow(auc_results))) {
            if (!is.na(auc_results$AUC[i])) {
              cat("   ", auc_results$Horizonte_anos[i], "ano(s):",
                  round(auc_results$AUC[i], 4),
                  "(eventos:", auc_results$N_eventos[i], ")\n")
            }
          }
          
          tabela_metricas <- criar_tabela_metricas(cindex_oob, cindex_test, auc_results)
          
          write.csv(
            tabela_metricas,
            file = file.path(dir_output_layer,
                             paste0("Tabela_Metricas_", tipo_cancer, "_", metrica, "_",
                                    omic_label, "_", pheno_label, ".csv")),
            row.names = FALSE
          )
          
          plotar_comparacao_metricas(
            cindex_oob = cindex_oob,
            cindex_test = cindex_test,
            auc_results = auc_results,
            tipo_cancer = tipo_cancer,
            metrica = metrica,
            omic_label = omic_label,
            pheno_label = pheno_label,
            output_dir = dir_output_layer
          )
          
          plotar_auc_temporal(
            auc_results = auc_results,
            tipo_cancer = tipo_cancer,
            metrica = metrica,
            omic_label = omic_label,
            pheno_label = pheno_label,
            output_dir = dir_output_layer
          )
          
          horizons_days <- horizons_years * 365
          t_int <- tryCatch(as.numeric(pred_test$time.interest), error = function(e) numeric(0))
          Smat  <- pred_test$survival
          
          if (length(t_int) > 1 && !is.null(Smat) && nrow(as.matrix(Smat)) == nrow(df_test_clean)) {
            
            Smat <- as.matrix(Smat)
            
            get_S_at_horizons <- function(i) {
              approx(
                x = t_int,
                y = Smat[i, ],
                xout = horizons_days,
                rule = 2,
                ties = "ordered"
              )$y
            }
            
            S_mat <- tryCatch({
              do.call(rbind, lapply(seq_len(nrow(Smat)), get_S_at_horizons))
            }, error = function(e) {
              cat("⚠️ Erro ao interpolar sobrevivência por horizonte:", e$message, "\n")
              NULL
            })
            
            if (!is.null(S_mat)) {
              P_event <- 1 - S_mat
              colnames(P_event) <- paste0("Prob_Event_", horizons_years, "y")
              
              patient_probs <- data.frame(
                ID = patient_ids_test,
                P_event_1y  = P_event[, 1],
                P_event_3y  = P_event[, 2],
                P_event_5y  = P_event[, 3],
                P_event_10y = P_event[, 4],
                stringsAsFactors = FALSE
              )
              
              write.table(
                patient_probs,
                file = file.path(dir_output_layer,
                                 paste0("Patient_Probabilities_", tipo_cancer, "_", metrica, "_",
                                        omic_label, "_", pheno_label, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE
              )
            } else {
              cat("⚠️ Probabilidades por paciente não puderam ser calculadas\n")
            }
            
          } else {
            cat("⚠️ time.interest/survival inválidos - pulando probabilidades por paciente\n")
          }
          
          write.csv(
            auc_results,
            file = file.path(dir_output_layer,
                             paste0("AUC_by_Horizon_", tipo_cancer, "_", metrica, "_",
                                    omic_label, "_", pheno_label, ".csv")),
            row.names = FALSE
          )
          
          # ==========================================================================
          # IMPORTÂNCIA E TABELA DE VARIÁVEIS USADAS
          # ==========================================================================
          predictor_names <- colnames(df_train_clean)[!colnames(df_train_clean) %in% c("time", "status")]
          
          importance_panel <- data.frame()
          importance_consenso <- data.frame()
          importance_boruta_df <- data.frame()
          n_importantes_consenso <- NA_integer_
          thr_vimp <- NA_real_
          thr_depth <- NA_real_
          thr_split <- NA_real_
          
          if (modo_vars %in% c(3, 4)) {
            
            importance_panel <- tryCatch({
              construir_painel_importancia(rsf_model, predictor_names)
            }, error = function(e) {
              cat("⚠️ Erro ao construir painel de importância do RSF final:", e$message, "\n")
              data.frame(
                Variable = predictor_names,
                VIMP = NA_real_,
                minimal_depth = NA_real_,
                split_frequency = NA_real_,
                crit_vimp = NA,
                crit_depth = NA,
                crit_split = NA,
                score_importancia = NA_integer_,
                importante_consenso = NA,
                rank_vimp = NA_real_,
                rank_depth = NA_real_,
                rank_split = NA_real_,
                rank_consenso = NA_real_,
                stringsAsFactors = FALSE
              )
            })
            
            thr_vimp  <- attr(importance_panel, "threshold_vimp")
            thr_depth <- attr(importance_panel, "threshold_depth")
            thr_split <- attr(importance_panel, "threshold_split")
            
            cat("📌 Thresholds de importância (RSF final):\n")
            cat("   VIMP >", thr_vimp, "\n")
            cat("   Minimal depth <", ifelse(is.na(thr_depth), "NA", round(thr_depth, 4)), "\n")
            cat("   Split frequency >", ifelse(is.na(thr_split), "NA", round(thr_split, 4)), "\n")
            
            n_importantes_consenso <- sum(importance_panel$importante_consenso, na.rm = TRUE)
            cat("📌 Variáveis importantes por consenso no RSF final:", n_importantes_consenso, "de", nrow(importance_panel), "\n")
            
            importance_panel_export <- importance_panel %>%
              dplyr::mutate(
                threshold_vimp = thr_vimp,
                threshold_depth = thr_depth,
                threshold_split = thr_split
              )
            
            importance_consenso <- importance_panel %>%
              dplyr::filter(importante_consenso)
            
            write.csv(
              importance_panel_export,
              file = file.path(dir_output_layer,
                               paste0("Importance_Panel_", tipo_cancer, "_", metrica, "_",
                                      omic_label, "_", pheno_label, ".csv")),
              row.names = FALSE
            )
            
            write.csv(
              importance_consenso,
              file = file.path(dir_output_layer,
                               paste0("Importance_ConsensusOnly_", tipo_cancer, "_", metrica, "_",
                                      omic_label, "_", pheno_label, ".csv")),
              row.names = FALSE
            )
            
            if (nrow(importance_panel_pre) > 0) {
              write.csv(
                importance_panel_pre,
                file = file.path(dir_output_layer,
                                 paste0("Importance_Panel_PRELIMINAR_", tipo_cancer, "_", metrica, "_",
                                        omic_label, "_", pheno_label, ".csv")),
                row.names = FALSE
              )
            }
            
            if (nrow(importance_consenso_pre) > 0) {
              write.csv(
                importance_consenso_pre,
                file = file.path(dir_output_layer,
                                 paste0("Variaveis_Selecionadas_PRELIMINAR_", tipo_cancer, "_", metrica, "_",
                                        omic_label, "_", pheno_label, ".csv")),
                row.names = FALSE
              )
            }
            
          } else {
            
            cat("📌 Modo com Boruta: será exportada apenas Importance + RelativeImp na tabela de variáveis.\n")
            
            importance_boruta_df <- tryCatch({
              extrair_importancia_boruta_modo(rsf_model, predictor_names)
            }, error = function(e) {
              cat("⚠️ Erro ao extrair importance para modo com Boruta:", e$message, "\n")
              data.frame(
                Variable = predictor_names,
                Importance = NA_real_,
                RelativeImp = NA_real_,
                stringsAsFactors = FALSE
              )
            })
          }
          
          tabela_variaveis_usadas <- montar_tabela_variaveis_usadas(
            predictor_names = predictor_names,
            selected_vars = selected_vars,
            clinicas_relevantes = clinicas_relevantes,
            tipo_cancer = tipo_cancer,
            metrica = metrica,
            omic_label = omic_label,
            pheno_label = pheno_label,
            modo_vars = modo_vars,
            omic_root = omic_root,
            pheno_root = pheno_root,
            importance_boruta_df = importance_boruta_df,
            importance_panel_df = importance_panel,
            thr_vimp = thr_vimp,
            thr_depth = thr_depth,
            thr_split = thr_split
          )
          
          write.csv(
            tabela_variaveis_usadas,
            file = file.path(dir_output_layer,
                             paste0("Variaveis_Usadas_No_Modelo_", tipo_cancer, "_", metrica, "_",
                                    omic_label, "_", pheno_label, ".csv")),
            row.names = FALSE
          )
          
          res_temp <- data.frame(
            Assinatura   = paste(selected_vars, collapse = ";"),
            Tipo         = tipo_cancer,
            Metrica      = metrica,
            Omic_layer   = omic_label,
            Pheno_layer  = pheno_label,
            Modo_vars    = modo_vars,
            Omic_mode    = omic_root,
            Pheno_mode   = pheno_root,
            NumGenes     = length(selected_vars),
            C_index_OOB  = cindex_oob,
            C_index_Test = cindex_test,
            C_index_Test_sem_reverse = cindex_test_sem_reverse,
            C_index_Test_com_reverse = cindex_test_com_reverse,
            Cindex_direcao_usada = direcao_cindex,
            AUC_1y = ifelse(1 %in% auc_results$Horizonte_anos,
                            auc_results$AUC[auc_results$Horizonte_anos == 1], NA_real_),
            AUC_3y = ifelse(3 %in% auc_results$Horizonte_anos,
                            auc_results$AUC[auc_results$Horizonte_anos == 3], NA_real_),
            AUC_5y = ifelse(5 %in% auc_results$Horizonte_anos,
                            auc_results$AUC[auc_results$Horizonte_anos == 5], NA_real_),
            AUC_10y = ifelse(10 %in% auc_results$Horizonte_anos,
                             auc_results$AUC[auc_results$Horizonte_anos == 10], NA_real_),
            N_eventos_1y = ifelse(1 %in% auc_results$Horizonte_anos,
                                  auc_results$N_eventos[auc_results$Horizonte_anos == 1], NA_integer_),
            N_eventos_3y = ifelse(3 %in% auc_results$Horizonte_anos,
                                  auc_results$N_eventos[auc_results$Horizonte_anos == 3], NA_integer_),
            N_eventos_5y = ifelse(5 %in% auc_results$Horizonte_anos,
                                  auc_results$N_eventos[auc_results$Horizonte_anos == 5], NA_integer_),
            N_eventos_10y = ifelse(10 %in% auc_results$Horizonte_anos,
                                   auc_results$N_eventos[auc_results$Horizonte_anos == 10], NA_integer_),
            NumVars_Modelo = length(predictor_names),
            NumVars_Important_Consensus = n_importantes_consenso,
            PropVars_Important_Consensus = ifelse(!is.na(n_importantes_consenso) && length(predictor_names) > 0,
                                                  n_importantes_consenso / length(predictor_names),
                                                  NA_real_),
            Top1_VIMP = ifelse(nrow(importance_panel) > 0 &&
                                 "VIMP" %in% colnames(importance_panel) &&
                                 !all(is.na(importance_panel$VIMP)),
                               suppressWarnings(max(importance_panel$VIMP, na.rm = TRUE)),
                               NA_real_),
            Mean_MinDepth = ifelse(nrow(importance_panel) > 0 &&
                                     "minimal_depth" %in% colnames(importance_panel) &&
                                     !all(is.na(importance_panel$minimal_depth)),
                                   mean(importance_panel$minimal_depth, na.rm = TRUE),
                                   NA_real_),
            Mean_SplitFreq = ifelse(nrow(importance_panel) > 0 &&
                                      "split_frequency" %in% colnames(importance_panel) &&
                                      !all(is.na(importance_panel$split_frequency)),
                                    mean(importance_panel$split_frequency, na.rm = TRUE),
                                    NA_real_),
            DF_Usado     = basename(df_file),
            stringsAsFactors = FALSE
          )
          
          verdadeiro_status <- df_test_clean$status
          
          if (all(is.na(marcador_risco_test)) || length(unique(na.omit(marcador_risco_test))) < 2) {
            cat("⚠️ Marcador de risco constante/NA — pulando matriz de confusão\n")
            goto_confusion <- FALSE
          } else {
            cutoff <- stats::median(marcador_risco_test, na.rm = TRUE)
            if (!is.finite(cutoff)) cutoff <- stats::quantile(marcador_risco_test, 0.5, na.rm = TRUE, type = 7)
            
            if (!is.finite(cutoff)) {
              cat("⚠️ Cutoff inválido — pulando matriz de confusão\n")
              goto_confusion <- FALSE
            } else {
              risco_binario <- ifelse(marcador_risco_test > cutoff, 1L, 0L)
              cc_cm <- complete.cases(risco_binario, verdadeiro_status)
              risco_binario     <- risco_binario[cc_cm]
              verdadeiro_status <- verdadeiro_status[cc_cm]
              
              goto_confusion <- (length(risco_binario) >= 5 &&
                                   length(unique(verdadeiro_status)) == 2 &&
                                   length(unique(risco_binario)) == 2)
              if (!goto_confusion) cat("⚠️ Amostras/classes insuficientes — pulando matriz de confusão\n")
            }
          }
          
          if (goto_confusion) {
            confusion <- caret::confusionMatrix(
              factor(risco_binario,     levels = c(0, 1)),
              factor(verdadeiro_status, levels = c(0, 1)),
              positive = "1"
            )
            
            cm <- confusion$table
            TN <- cm[1, 1]; FN <- cm[1, 2]
            FP <- cm[2, 1]; TP <- cm[2, 2]
            
            Precision   <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
            Recall      <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
            Specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
            Accuracy    <- as.numeric(confusion$overall["Accuracy"])
            F1_score    <- ifelse(!is.na(Precision) && !is.na(Recall) && (Precision + Recall) > 0,
                                  2 * (Precision * Recall) / (Precision + Recall), NA_real_)
            FPR         <- ifelse(!is.na(Specificity), 1 - Specificity, NA_real_)
            
            Prop_FP <- ifelse((FP + TN) > 0, FP / (FP + TN), NA_real_)
            Prop_FN <- ifelse((FN + TP) > 0, FN / (FN + TP), NA_real_)
            
            res_temp <- res_temp %>%
              dplyr::mutate(
                Predicted_Risk     = TP + FP,
                Predicted_no_Risk  = FN + TN,
                Actual_Risk        = TP + FN,
                Actual_no_Risk     = FP + TN,
                Specificity        = Specificity,
                Accuracy           = Accuracy,
                FPR                = FPR,
                Prop_FP            = Prop_FP,
                Prop_FN            = Prop_FN,
                Precision          = Precision,
                Recall             = Recall,
                F1_score           = F1_score
              )
            
            row_metrics <- data.frame(
              AUC               = res_temp$AUC_5y,
              F1_score          = F1_score,
              Recall            = Recall,
              Precision         = Precision,
              Prop_FN           = Prop_FN,
              Prop_FP           = Prop_FP,
              FPR               = FPR,
              Accuracy          = Accuracy,
              Specificity       = Specificity,
              Actual_no_Risk    = FP + TN,
              Actual_Risk       = TP + FN,
              Predicted_no_Risk = FN + TN,
              Predicted_Risk    = TP + FP,
              stringsAsFactors  = FALSE
            )
            
          } else {
            res_temp <- res_temp %>%
              dplyr::mutate(
                Predicted_Risk     = NA_real_,
                Predicted_no_Risk  = NA_real_,
                Actual_Risk        = NA_real_,
                Actual_no_Risk     = NA_real_,
                Specificity        = NA_real_,
                Accuracy           = NA_real_,
                FPR                = NA_real_,
                Prop_FP            = NA_real_,
                Prop_FN            = NA_real_,
                Precision          = NA_real_,
                Recall             = NA_real_,
                F1_score           = NA_real_
              )
            
            row_metrics <- data.frame(
              AUC               = res_temp$AUC_5y,
              F1_score          = NA_real_,
              Recall            = NA_real_,
              Precision         = NA_real_,
              Prop_FN           = NA_real_,
              Prop_FP           = NA_real_,
              FPR               = NA_real_,
              Accuracy          = NA_real_,
              Specificity       = NA_real_,
              Actual_no_Risk    = NA_real_,
              Actual_Risk       = NA_real_,
              Predicted_no_Risk = NA_real_,
              Predicted_Risk    = NA_real_,
              stringsAsFactors  = FALSE
            )
          }
          
          row_metrics_transposta <- as.data.frame(t(row_metrics))
          row_metrics_transposta <- tibble::rownames_to_column(row_metrics_transposta, var = "Metrica")
          colnames(row_metrics_transposta)[2] <- "Valor"
          
          write.csv(
            row_metrics_transposta,
            file = file.path(dir_output_layer,
                             paste0("Performance_Table_", tipo_cancer, "_", metrica, "_",
                                    omic_label, "_", pheno_label, ".csv")),
            row.names = FALSE
          )
          
          radar_data <- data.frame(
            AUC_5y      = res_temp$AUC_5y,
            F1_score    = res_temp$F1_score,
            Recall      = res_temp$Recall,
            Specificity = res_temp$Specificity,
            Precision   = res_temp$Precision,
            Prop_FP     = res_temp$Prop_FP,
            Prop_FN     = res_temp$Prop_FN
          )
          radar_data <- pmin(pmax(radar_data, 0), 1)
          radar_data <- rbind(rep(1, ncol(radar_data)), rep(0, ncol(radar_data)), radar_data)
          
          plot_title <- paste("Radar Performance Metrics -", tipo_cancer, "-", metrica,
                              "(", omic_label, "|", pheno_label, ")")
          
          pdf(file = file.path(dir_output_layer,
                               paste0("Radar_", tipo_cancer, "_", metrica, "_",
                                      omic_label, "_", pheno_label, ".pdf")),
              width = 9, height = 9)
          radarchart(
            radar_data,
            axistype = 1,
            pcol  = "darkgreen",
            pfcol = scales::alpha("darkgreen", 0.3),
            plwd  = 2,
            title = plot_title,
            cglcol = "grey", cglty = 1,
            axislabcol = "black",
            caxislabels = seq(0, 1, 0.25),
            cglwd = 0.8
          )
          dev.off()
          
          try({
            pdf(file = file.path(dir_output_layer,
                                 paste0("RSF_Combined_", tipo_cancer, "_", metrica, "_",
                                        omic_label, "_", pheno_label, ".pdf")),
                width = 12, height = 6)
            par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
            plot(rsf_model, main = paste("OOB Error -", tipo_cancer, "-", metrica,
                                         "(", omic_label, "|", pheno_label, ")"))
            dev.off()
          }, silent = TRUE)
          
          resultados_geral <- dplyr::bind_rows(resultados_geral, res_temp)
          progress_keys <- unique(c(progress_keys, key))
          
          save_global_progress(
            progress_file = progress_file,
            resultados_geral = resultados_geral,
            progress_keys = progress_keys
          )
          
          cat("💾 Progresso global salvo em:", progress_file, "\n")
          
        } # fim camada
      } # fim métrica
    } # fim tipo
  } # fim omic_vars
} # fim modo_vars

# ======================================================================================
# 8) RESUMO GERAL AO FINAL DA EXECUÇÃO
# ======================================================================================
cat("\n📊 Gerando resumos finais...\n")

resumo_dir <- file.path(getwd(), "output", "resumos")
if (!dir.exists(resumo_dir)) dir.create(resumo_dir, recursive = TRUE)

if (nrow(resultados_geral) > 0) {
  
  resumo_geral <- resultados_geral %>%
    group_by(Tipo, Metrica) %>%
    summarise(
      N_modelos = n(),
      Media_Cindex_OOB = mean(C_index_OOB, na.rm = TRUE),
      SD_Cindex_OOB = sd(C_index_OOB, na.rm = TRUE),
      Media_Cindex_Test = mean(C_index_Test, na.rm = TRUE),
      SD_Cindex_Test = sd(C_index_Test, na.rm = TRUE),
      Media_AUC_1y = mean(AUC_1y, na.rm = TRUE),
      SD_AUC_1y = sd(AUC_1y, na.rm = TRUE),
      Media_AUC_3y = mean(AUC_3y, na.rm = TRUE),
      SD_AUC_3y = sd(AUC_3y, na.rm = TRUE),
      Media_AUC_5y = mean(AUC_5y, na.rm = TRUE),
      SD_AUC_5y = sd(AUC_5y, na.rm = TRUE),
      Media_AUC_10y = mean(AUC_10y, na.rm = TRUE),
      SD_AUC_10y = sd(AUC_10y, na.rm = TRUE),
      Media_PropVars_Important_Consensus = mean(PropVars_Important_Consensus, na.rm = TRUE),
      Media_NumVars_Important_Consensus = mean(NumVars_Important_Consensus, na.rm = TRUE),
      .groups = 'drop'
    )
  
  write.csv(
    resumo_geral,
    file = file.path(resumo_dir, "Resumo_Geral_Metricas.csv"),
    row.names = FALSE
  )
  
  resumo_modo <- resultados_geral %>%
    group_by(Modo_vars, Omic_mode, Pheno_mode, Metrica) %>%
    summarise(
      N_modelos = n(),
      Media_Cindex_Test = mean(C_index_Test, na.rm = TRUE),
      Media_AUC_5y = mean(AUC_5y, na.rm = TRUE),
      Media_PropVars_Important_Consensus = mean(PropVars_Important_Consensus, na.rm = TRUE),
      .groups = 'drop'
    )
  
  write.csv(
    resumo_modo,
    file = file.path(resumo_dir, "Resumo_por_Modo.csv"),
    row.names = FALSE
  )
  
  if (nrow(resultados_geral) > 10) {
    cols_metricas <- c("C_index_OOB", "C_index_Test", "AUC_1y", "AUC_3y", "AUC_5y", "AUC_10y")
    cols_presentes <- intersect(cols_metricas, colnames(resultados_geral))
    
    if (length(cols_presentes) >= 3) {
      metricas_cor <- resultados_geral %>%
        dplyr::select(all_of(cols_presentes)) %>%
        cor(use = "pairwise.complete.obs")
      
      pdf(file = file.path(resumo_dir, "Correlacao_Metricas.pdf"),
          width = 8, height = 7)
      corrplot(metricas_cor, method = "color", type = "upper",
               tl.col = "black", tl.srt = 45,
               title = "Correlação entre Métricas de Desempenho",
               mar = c(0, 0, 2, 0))
      dev.off()
      
      write.csv(
        metricas_cor,
        file = file.path(resumo_dir, "Matriz_Correlacao.csv")
      )
    }
  }
  
  if (nrow(resultados_geral) > 5) {
    p_comp <- ggplot(resultados_geral, aes(x = C_index_Test, y = AUC_5y, color = Metrica)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = FALSE, alpha = 0.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      labs(
        title = "Comparação: C-index Teste vs AUC 5 anos",
        x = "C-index (Teste)",
        y = "AUC (5 anos)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )
    
    ggsave(
      filename = file.path(resumo_dir, "Comparacao_Cindex_vs_AUC5y.pdf"),
      plot = p_comp,
      width = 8,
      height = 7
    )
  }
  
} else {
  cat("⚠️ resultados_geral está vazio. Nenhum resumo final foi gerado.\n")
}

cat("\n✅ Execução finalizada completamente.\n")
cat("📁 Arquivo global de progresso:\n")
cat("   ", progress_file, "\n")
cat("📁 Resumos salvos em:\n")
cat("   ", resumo_dir, "\n")

# ======================================================================================
# TEMPO GLOBAL DA EXECUÇÃO - FIM
# ======================================================================================
execucao_fim <- Sys.time()
duracao_total_formatada <- formatar_duracao_global(execucao_inicio, execucao_fim)
duracao_total_segundos <- as.numeric(difftime(execucao_fim, execucao_inicio, units = "secs"))

cat("\n")
cat("============================================================\n")
cat("🕒 FIM DA EXECUÇÃO GLOBAL\n")
cat("Horário de término:", format(execucao_fim, "%Y-%m-%d %H:%M:%S"), "\n")
cat("⏱️ Tempo total de execução:", duracao_total_formatada, "\n")
cat("⏱️ Tempo total em segundos:", round(duracao_total_segundos, 2), "\n")
cat("============================================================\n\n")

tempo_execucao_global <- data.frame(
  inicio = format(execucao_inicio, "%Y-%m-%d %H:%M:%S"),
  fim = format(execucao_fim, "%Y-%m-%d %H:%M:%S"),
  duracao_segundos = round(duracao_total_segundos, 2),
  duracao_formatada = duracao_total_formatada,
  stringsAsFactors = FALSE
)

arquivo_tempo_global <- file.path(getwd(), "output", "resumos", "tempo_execucao_global.csv")

if (!dir.exists(dirname(arquivo_tempo_global))) {
  dir.create(dirname(arquivo_tempo_global), recursive = TRUE, showWarnings = FALSE)
}

write.csv(tempo_execucao_global, arquivo_tempo_global, row.names = FALSE)

cat("📁 Log de tempo salvo em:\n")
cat("   ", arquivo_tempo_global, "\n")

# ======================================================================================
# 8) SALVAR LOG DE AUDITORIA
# ======================================================================================
audit_df <- do.call(rbind, audit_log_v4)
if (!is.null(audit_df) && nrow(audit_df) > 0) {
  audit_path <- file.path(base_output_dir, "..", "Audit_Log_Modelos_V4_Pulos.csv")
  write.csv(audit_df, audit_path, row.names = FALSE)
  cat("\n✅ Log de Auditoria V4 salvo em:", audit_path, "\n")
}