# ==============================================================================
# XAI SUPERLEARNER: SurvSHAP & SurvLIME  (V8)
# Execução independente + paralela (Linux, mclapply)
# + CHECKPOINT: salva progresso em RDS e retoma de onde parou
# ==============================================================================
# CORREÇÕES V8 (sobre V7):
#
#   [FIX-D] Beeswarm: categóricas em cinza
#           Variáveis categóricas (factor/character em x_all_df) agora recebem
#           uma paleta discreta de até 12 cores via ggnewscale::new_scale_color(),
#           com legenda "variavel=categoria" no lado direito. Numéricas mantêm
#           o gradiente contínuo azul→vermelho em escala independente.
#           Requer: install.packages("ggnewscale")
#
#   [FIX-A] Can't combine <factor> and <double>
#           safe_contrib_df constrói todos os data.frames individuais já com
#           stringsAsFactors=FALSE + as.character() explícito. safe_bind_rows()
#           mantido como segunda linha de defesa.
#
#   [FIX-B] método não aplicável para 'slice'
#           randomForestSRC exporta slice() e sobrescreve dplyr::slice().
#           Corrigido com dplyr::slice() explícito em plot_heatmap().
#
#   [FIX-C] substituto tem 0 linha, dados têm 300
#           lm() substituído por glmnet ridge (alpha=0): numericamente estável
#           em qualquer cenário de colinearidade.
#
#   [FIX-E] In argument: patient_id == pid  (NOVO V8)
#           patient_id pode ser factor em all_shap_df/all_lime_df. Forçado
#           as.character() antes do loop de pacientes e uso de indexação
#           base R ([df$col == val]) em vez de dplyr::filter para evitar
#           conflitos de tipo.
#
#   [FIX-F] Nenhum modelo com peso != 0  (NOVO V8)
#           Detecção robusta do nome da coluna de peso independentemente de
#           case/nome. Garante conversão numérica. Mostra diagnóstico no log.
#           Detecta também colunas Modelo e Algoritmo de forma flexível.
#           Infere Algoritmo pelo nome do modelo se coluna não existir.
#
#   [FIX-G] Nenhum resultado SHAP/LIME gerado  (NOVO V8)
#           Mostra caminhos tentados quando bundle/XAI não encontrado, para
#           facilitar diagnóstico. Contadores de tentativas e falhas por loop.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(survival)
  library(glmnet)
  library(randomForestSRC)
  library(xgboost)
  library(Matrix)
  library(parallel)
  library(ggplot2)
  library(scales)
  library(ggbeeswarm)
  library(ggnewscale)
})

# [FIX-B] Forçar verbos dplyr sobre randomForestSRC
slice     <- dplyr::slice
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
arrange   <- dplyr::arrange
group_by  <- dplyr::group_by
ungroup   <- dplyr::ungroup
left_join <- dplyr::left_join
inner_join  <- dplyr::inner_join
full_join   <- dplyr::full_join
rename      <- dplyr::rename

# ==============================================================================
# CONFIGURAÇÃO
# ==============================================================================
RSF_ROOT     <- "D:/ML_Emanuell/ML_RSF/linear"
XGB_ROOT     <- "D:/ML_Emanuell/ML_XGBoost"
SL_ROOT      <- "D:/superLERANING/PostHoc_Octo_SuperLearner"
XAI_OUT_ROOT <- "D:/superLERANING/XAI_SurvSHAP_SurvLIME"

PROGRESS_PATH   <- file.path(SL_ROOT, "checkpoint_progress_superlearner_v3_overfitting.rds")
N_WORKERS       <- 8
N_LIME_PERT     <- 300
TOP_N_FEAT      <- 20
TOP_N_GLOBAL    <- 30
HEATMAP_MAX_PAT <- 80

XAI_CHECKPOINT_PATH <- file.path(XAI_OUT_ROOT, "XAI_checkpoint_progress.rds")

if (!dir.exists(XAI_OUT_ROOT)) dir.create(XAI_OUT_ROOT, recursive = TRUE)

# ==============================================================================
# HELPERS DE CHECKPOINT
# ==============================================================================

unit_key <- function(u) paste(u$Tipo, u$Metrica, u$Omic, sep = "|")

read_checkpoint <- function(path = XAI_CHECKPOINT_PATH) {
  if (file.exists(path)) {
    tryCatch(readRDS(path), error = function(e) list())
  } else {
    list()
  }
}

write_checkpoint_entry <- function(key, status, error_msg = NULL,
                                   path = XAI_CHECKPOINT_PATH) {
  lock_path <- paste0(path, ".lock")
  
  waited <- 0
  while (file.exists(lock_path) && waited < 30) {
    Sys.sleep(0.2); waited <- waited + 0.2
  }
  
  tryCatch(writeLines(as.character(Sys.getpid()), lock_path), error = function(e) NULL)
  on.exit({
    if (file.exists(lock_path)) file.remove(lock_path)
  }, add = TRUE)
  
  cp        <- read_checkpoint(path)
  cp[[key]] <- list(
    key       = key,
    status    = status,
    timestamp = Sys.time(),
    error_msg = error_msg
  )
  
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(cp, tmp)
  file.rename(tmp, path)
  
  invisible(NULL)
}

is_done <- function(key, cp) {
  !is.null(cp[[key]]) && identical(cp[[key]]$status, "done")
}

summarise_checkpoint <- function(cp) {
  if (length(cp) == 0) {
    cat("  Checkpoint vazio — nenhuma unidade processada anteriormente.\n")
    return(invisible(NULL))
  }
  statuses <- sapply(cp, `[[`, "status")
  cat(sprintf("  Checkpoint: %d done | %d error | %d total\n",
              sum(statuses == "done"),
              sum(statuses == "error"),
              length(statuses)))
}

# ==============================================================================
# HELPERS: localização de arquivos
# ==============================================================================
find_rsf_bundle <- function(cancer, endpoint, omic, modulo) {
  for (mode in c("multiomica", "mono_omica")) {
    p <- file.path(
      RSF_ROOT, "output", "modelo", paste0("modelo", modulo),
      mode, "multifenotipo", cancer, endpoint, omic,
      "Explicabilidade",
      paste0("XAI_bundle_", cancer, "_", endpoint, "_", omic, "_allPheno.rds")
    )
    if (file.exists(p)) return(p)
  }
  NULL
}

find_xgb_xai <- function(cancer, endpoint, omic, modulo) {
  for (mode in c("multiomica", "mono_omica")) {
    base  <- file.path(
      XGB_ROOT, "output", "modelo", paste0("modelo", modulo),
      mode, "multifenotipo", cancer, endpoint, omic, "XAI_ready"
    )
    m_p <- file.path(base, paste0("xgb_model_", cancer, "_", endpoint, "_", omic, ".rds"))
    d_p <- file.path(base, paste0("data_xai_",  cancer, "_", endpoint, "_", omic, ".rds"))
    if (file.exists(m_p) && file.exists(d_p))
      return(list(model_path = m_p, data_path = d_p))
  }
  NULL
}

# ==============================================================================
# [FIX-A] safe_contrib_df: tipos seguros desde a origem
# ==============================================================================
safe_contrib_df <- function(patient_id, variable, value, method) {
  data.frame(
    patient_id = as.character(patient_id),
    variable   = as.character(variable),
    value      = as.numeric(value),
    method     = as.character(method),
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# PALETA E TEMA
# ==============================================================================
CLR_POS  <- "#c0392b"
CLR_NEG  <- "#2980b9"
CLR_SHAP <- "#e67e22"
CLR_LIME <- "#8e44ad"

LBL_POS <- "Risco (+)"
LBL_NEG <- "Protecao (-)"

theme_xai <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title      = element_text(face = "bold", size = base_size + 1),
      plot.subtitle   = element_text(color = "grey40", size = base_size - 1),
      legend.position = "bottom",
      axis.text.y     = element_text(size = base_size - 3)
    )
}

format_var_label <- function(var_name, raw_value) {
  if (is.null(raw_value) || length(raw_value) == 0 || is.na(raw_value))
    return(var_name)
  val_str <- if (is.numeric(raw_value)) sprintf("%.2f", raw_value) else as.character(raw_value)
  paste0(var_name, " = ", val_str)
}

# ==============================================================================
# PLOT 1 & 2: WATERFALL individual
# ==============================================================================
plot_waterfall <- function(df_contrib, patient_id, cancer, endpoint,
                           method_label, x_patient = NULL,
                           top_n = TOP_N_FEAT) {
  
  df <- df_contrib %>%
    dplyr::filter(is.finite(value)) %>%
    dplyr::arrange(desc(abs(value))) %>%
    head(top_n)
  
  if (nrow(df) == 0) return(NULL)
  
  df <- df %>%
    rowwise() %>%
    dplyr::mutate(
      var_label = if (!is.null(x_patient) && variable %in% names(x_patient)) {
        format_var_label(variable, x_patient[[variable]])
      } else { variable },
      direction = if (value >= 0) LBL_POS else LBL_NEG
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(var_label = factor(var_label, levels = var_label[order(value)]))
  
  max_val <- max(abs(df$value), na.rm = TRUE)
  
  ggplot(df, aes(x = var_label, y = value, fill = direction)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.3) +
    geom_text(
      aes(label = sprintf("%.4f", value)),
      hjust = ifelse(df$value >= 0, -0.10, 1.10),
      size  = 3.0, color = "grey20"
    ) +
    coord_flip(clip = "off") +
    scale_y_continuous(
      limits = c(-max_val * 1.15, max_val * 1.30),
      expand = expansion(mult = c(0.03, 0.08))
    ) +
    scale_fill_manual(values = setNames(c(CLR_POS, CLR_NEG), c(LBL_POS, LBL_NEG))) +
    labs(
      title    = paste0(method_label, "  --  Paciente: ", patient_id),
      subtitle = paste0("Cancer: ", cancer, "  |  Endpoint: ", endpoint),
      x = NULL, y = "Contribuicao para o risco do ensemble", fill = NULL
    ) +
    theme_xai() +
    theme(plot.margin = margin(5, 110, 5, 5))
}

# ==============================================================================
# PLOT 3: SHAP vs LIME AGREEMENT individual
# ==============================================================================
plot_shap_lime_agreement <- function(df_shap, df_lime, patient_id, cancer, endpoint,
                                     top_n = TOP_N_FEAT) {
  
  merged <- dplyr::full_join(
    df_shap %>% dplyr::select(variable, shap = value),
    df_lime  %>% dplyr::select(variable, lime = value),
    by = "variable"
  ) %>%
    dplyr::mutate(
      shap       = tidyr::replace_na(shap, 0),
      lime       = tidyr::replace_na(lime, 0),
      importance = (abs(shap) + abs(lime)) / 2,
      concordant = sign(shap) == sign(lime)
    ) %>%
    dplyr::arrange(desc(importance)) %>%
    head(top_n) %>%
    dplyr::mutate(variable = factor(variable, levels = variable[order((shap + lime) / 2)]))
  
  if (nrow(merged) == 0) return(NULL)
  
  df_long <- merged %>%
    dplyr::select(variable, SHAP = shap, LIME = lime) %>%
    tidyr::pivot_longer(cols = c(SHAP, LIME), names_to = "metodo", values_to = "value")
  
  ggplot(df_long, aes(x = variable, y = value, fill = metodo)) +
    geom_col(position = position_dodge(width = 0.65), width = 0.60,
             color = "white", linewidth = 0.25) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50") +
    coord_flip() +
    scale_fill_manual(values = c(SHAP = CLR_SHAP, LIME = CLR_LIME)) +
    labs(
      title    = paste0("SHAP vs LIME Agreement  --  Paciente: ", patient_id),
      subtitle = paste0(
        "Cancer: ", cancer, "  |  Endpoint: ", endpoint,
        "  |  Concordantes: ", sum(merged$concordant), "/", nrow(merged)
      ),
      x = NULL, y = "Contribuicao para o risco", fill = "Metodo"
    ) +
    theme_xai()
}

# ==============================================================================
# PLOT 4 & 5: BEESWARM POPULACIONAL
# ==============================================================================
PALETTE_CAT <- c(
  "#e41a1c","#377eb8","#4daf4a","#984ea3",
  "#ff7f00","#a65628","#f781bf","#999999",
  "#1b9e77","#d95f02","#7570b3","#e7298a"
)

plot_beeswarm <- function(df_all, method_label, x_all_df = NULL,
                          top_n = TOP_N_GLOBAL) {
  
  top_vars <- df_all %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(imp = mean(abs(value), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(imp)) %>%
    head(top_n) %>%
    dplyr::pull(variable)
  
  df_plot <- df_all %>%
    dplyr::filter(variable %in% top_vars) %>%
    dplyr::mutate(variable = factor(variable, levels = rev(top_vars)))
  
  if (nrow(df_plot) == 0) return(NULL)
  
  is_cat_var <- function(vname, x_df) {
    if (is.null(x_df) || !vname %in% colnames(x_df)) return(FALSE)
    col <- x_df[[vname]]
    is.factor(col) || is.character(col) ||
      (!is.numeric(col) && length(unique(col[!is.na(col)])) <= 10)
  }
  
  cat_vars <- if (!is.null(x_all_df)) {
    top_vars[sapply(top_vars, is_cat_var, x_df = x_all_df)]
  } else character(0)
  
  num_vars <- setdiff(top_vars, cat_vars)
  
  df_plot <- df_plot %>%
    dplyr::mutate(feat_z    = NA_real_,
                  raw_value = NA_real_,
                  cat_lbl   = NA_character_)
  
  num_range_tbl <- data.frame(variable = character(), v_min = numeric(),
                              v_max = numeric(), stringsAsFactors = FALSE)
  
  if (!is.null(x_all_df) && length(num_vars) > 0) {
    common_num <- intersect(num_vars, colnames(x_all_df))
    if (length(common_num) > 0) {
      x_num_long <- x_all_df %>%
        dplyr::select(dplyr::all_of(common_num)) %>%
        dplyr::mutate(patient_id = rownames(.)) %>%
        tidyr::pivot_longer(-patient_id, names_to = "variable", values_to = "raw_value") %>%
        dplyr::mutate(raw_value = suppressWarnings(as.numeric(raw_value))) %>%
        dplyr::group_by(variable) %>%
        dplyr::mutate(
          v_min  = min(raw_value, na.rm = TRUE),
          v_max  = max(raw_value, na.rm = TRUE),
          feat_z = ifelse(v_max > v_min,
                          (raw_value - v_min) / (v_max - v_min),
                          0.5)
        ) %>%
        dplyr::ungroup()
      
      num_range_tbl <- x_num_long %>%
        dplyr::group_by(variable) %>%
        dplyr::summarise(v_min = first(v_min), v_max = first(v_max), .groups = "drop") %>%
        as.data.frame()
      
      x_num_join <- x_num_long %>% dplyr::select(patient_id, variable, feat_z, raw_value)
      
      df_plot <- df_plot %>%
        dplyr::left_join(x_num_join, by = c("patient_id", "variable"),
                         suffix = c("", ".new")) %>%
        dplyr::mutate(
          feat_z    = dplyr::coalesce(feat_z.new,    feat_z),
          raw_value = dplyr::coalesce(raw_value.new, raw_value)
        ) %>%
        dplyr::select(-dplyr::any_of(c("feat_z.new", "raw_value.new")))
    }
  }
  
  if (!is.null(x_all_df) && length(cat_vars) > 0) {
    common_cat <- intersect(cat_vars, colnames(x_all_df))
    if (length(common_cat) > 0) {
      x_cat_long <- x_all_df %>%
        dplyr::select(dplyr::all_of(common_cat)) %>%
        dplyr::mutate(patient_id = rownames(.)) %>%
        tidyr::pivot_longer(-patient_id, names_to = "variable", values_to = "cat_val") %>%
        dplyr::mutate(
          cat_val = as.character(cat_val),
          cat_lbl = paste0(variable, " = ", cat_val)
        ) %>%
        dplyr::select(patient_id, variable, cat_lbl)
      
      df_plot <- df_plot %>%
        dplyr::left_join(x_cat_long, by = c("patient_id", "variable"),
                         suffix = c("", ".new")) %>%
        dplyr::mutate(cat_lbl = dplyr::coalesce(cat_lbl.new, cat_lbl)) %>%
        dplyr::select(-dplyr::any_of("cat_lbl.new"))
    }
  }
  
  all_cat_lbls <- sort(unique(df_plot$cat_lbl[!is.na(df_plot$cat_lbl)]))
  n_cat        <- length(all_cat_lbls)
  cat_colors   <- if (n_cat > 0) {
    cols <- rep(PALETTE_CAT, length.out = n_cat)
    setNames(cols, all_cat_lbls)
  } else setNames(character(0), character(0))
  
  df_num <- df_plot %>% dplyr::filter(!is.na(feat_z))
  df_cat <- df_plot %>% dplyr::filter(!is.na(cat_lbl))
  df_unk <- df_plot %>% dplyr::filter(is.na(feat_z) & is.na(cat_lbl))
  
  has_num <- nrow(df_num) > 0
  has_cat <- nrow(df_cat) > 0
  
  if (has_num && nrow(num_range_tbl) > 0) {
    global_min <- min(num_range_tbl$v_min, na.rm = TRUE)
    global_max <- max(num_range_tbl$v_max, na.rm = TRUE)
    fmt_num <- function(x) {
      if (abs(x) >= 1000 || (abs(x) < 0.01 && x != 0))
        formatC(x, format = "e", digits = 2)
      else
        formatC(x, format = "f", digits = 2)
    }
    colorbar_lbl_low  <- fmt_num(global_min)
    colorbar_lbl_high <- fmt_num(global_max)
    colorbar_title    <- paste0("Expressao\n",
                                colorbar_lbl_low, " (azul)  \u2192  ",
                                colorbar_lbl_high, " (vermelho)")
  } else {
    colorbar_title <- "Expressao\n(baixo \u2192 alto)"
    fmt_num <- function(x) formatC(x, format = "f", digits = 2)
    global_min <- NA; global_max <- NA
  }
  
  p <- ggplot(df_plot, aes(x = variable, y = value)) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "grey70", linetype = "dashed") +
    coord_flip()
  
  if (has_num) {
    p <- p +
      ggnewscale::new_scale_color() +
      geom_quasirandom(
        data     = df_num,
        aes(color = feat_z),
        groupOnX = TRUE, size = 1.6, alpha = 0.78, bandwidth = 0.4, varwidth = TRUE
      ) +
      scale_color_gradient2(
        low      = "#2166ac",
        mid      = "#f7f7f7",
        high     = "#b2182b",
        midpoint = 0.5,
        limits   = c(0, 1),
        breaks   = c(0, 0.5, 1),
        labels   = if (!is.na(global_min)) {
          mid_val <- (global_min + global_max) / 2
          c(colorbar_lbl_low, fmt_num(mid_val), colorbar_lbl_high)
        } else c("Baixo", "Medio", "Alto"),
        name  = colorbar_title,
        guide = guide_colorbar(
          direction      = "horizontal",
          barwidth       = unit(8, "cm"),
          barheight      = unit(0.35, "cm"),
          title.position = "top",
          title.hjust    = 0.5,
          ticks.colour   = "grey40",
          frame.colour   = "grey60",
          label.theme    = element_text(size = 7),
          order          = 1
        )
      )
  }
  
  if (has_cat) {
    n_cols_legend <- min(3L, ceiling(n_cat / 4L))
    n_cols_legend <- max(1L, n_cols_legend)
    
    p <- p +
      ggnewscale::new_scale_color() +
      geom_quasirandom(
        data     = df_cat,
        aes(color = cat_lbl),
        groupOnX = TRUE, size = 1.9, alpha = 0.85, bandwidth = 0.4, varwidth = TRUE
      ) +
      scale_color_manual(
        values = cat_colors,
        name   = "Categoria",
        guide  = guide_legend(
          ncol           = n_cols_legend,
          byrow          = TRUE,
          override.aes   = list(size = 3.5, alpha = 1),
          title.position = "top",
          title.hjust    = 0,
          label.theme    = element_text(size = 7),
          order          = 2
        )
      )
  }
  
  if (nrow(df_unk) > 0) {
    p <- p +
      geom_quasirandom(
        data     = df_unk,
        color    = "grey55",
        groupOnX = TRUE, size = 1.6, alpha = 0.55, bandwidth = 0.4, varwidth = TRUE
      )
  }
  
  subtitle_txt <- paste0(
    "Top ", top_n, " variaveis  |  cada ponto = 1 paciente  |  ",
    if (has_num && has_cat) "numericas: gradiente de expressao  |  categoricas: cores por grupo"
    else if (has_num)       "cor: valor de expressao (azul = baixo, vermelho = alto)"
    else if (has_cat)       "cor por categoria"
    else                    ""
  )
  
  p +
    labs(
      title    = paste0(method_label, " - Distribuicao populacional das contribuicoes"),
      subtitle = subtitle_txt,
      x = NULL, y = paste0("Contribuicao ", method_label, " para o risco (ensemble)")
    ) +
    theme_xai(base_size = 11) +
    theme(
      legend.position    = "bottom",
      legend.box         = "vertical",
      legend.box.just    = "left",
      legend.margin      = margin(t = 6, b = 2),
      legend.spacing.y   = unit(0.3, "cm"),
      legend.text        = element_text(size = 7),
      legend.title       = element_text(size = 8, face = "bold"),
      legend.key.size    = unit(0.38, "cm"),
      axis.text.y        = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(t = 5, r = 10, b = 5, l = 5)
    )
}

# ==============================================================================
# PLOT 6 & 7: HEATMAP
# [FIX-B] dplyr::slice() explícito
# ==============================================================================
plot_heatmap <- function(df_all, method_label, top_n = TOP_N_GLOBAL,
                         max_pat = HEATMAP_MAX_PAT) {
  
  top_vars <- df_all %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(imp = mean(abs(value), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(imp)) %>%
    head(top_n) %>%
    dplyr::pull(variable)
  
  df_heat <- df_all %>%
    dplyr::filter(variable %in% top_vars) %>%
    dplyr::group_by(patient_id, variable) %>%
    dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")
  
  pat_risk <- df_heat %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(total_risk = sum(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(total_risk))
  
  if (nrow(pat_risk) > max_pat) {
    idx      <- seq(1, nrow(pat_risk), length.out = max_pat)
    pat_risk <- dplyr::slice(pat_risk, round(idx))
  }
  
  df_heat <- df_heat %>%
    dplyr::filter(patient_id %in% pat_risk$patient_id) %>%
    dplyr::mutate(
      patient_id = factor(patient_id, levels = rev(pat_risk$patient_id)),
      variable   = factor(variable,   levels = rev(top_vars))
    )
  
  if (nrow(df_heat) == 0) return(NULL)
  
  lim_val <- max(abs(df_heat$value), na.rm = TRUE)
  
  ggplot(df_heat, aes(x = variable, y = patient_id, fill = value)) +
    geom_tile(color = "white", linewidth = 0.18) +
    scale_fill_gradient2(
      low = CLR_NEG, mid = "white", high = CLR_POS,
      midpoint = 0, limits = c(-lim_val, lim_val), name = "Contribuicao"
    ) +
    labs(
      title    = paste0("Heatmap ", method_label, "  --  Pacientes x Variaveis"),
      subtitle = paste0(
        "Top ", top_n, " variaveis  |  ",
        dplyr::n_distinct(df_heat$patient_id), " pacientes (ordenados por risco total)"
      ),
      x = "Variavel", y = "Paciente"
    ) +
    theme_xai(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 5),
      panel.grid  = element_blank()
    )
}

# ==============================================================================
# PLOT 8: SHAP-LIME AGREEMENT SUMMARY global
# ==============================================================================
plot_agreement_summary <- function(df_shap, df_lime, top_n = TOP_N_GLOBAL) {
  
  merged <- dplyr::inner_join(
    df_shap %>% dplyr::select(patient_id, variable, shap = value),
    df_lime  %>% dplyr::select(patient_id, variable, lime = value),
    by = c("patient_id", "variable")
  ) %>%
    dplyr::mutate(concordant = sign(shap) == sign(lime))
  
  if (nrow(merged) == 0) return(NULL)
  
  summary_df <- merged %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      pct_concordant = mean(concordant, na.rm = TRUE) * 100,
      mean_shap_abs  = mean(abs(shap), na.rm = TRUE),
      n_pat          = dplyr::n(),
      .groups        = "drop"
    ) %>%
    dplyr::arrange(desc(mean_shap_abs)) %>%
    head(top_n) %>%
    dplyr::mutate(
      variable        = factor(variable, levels = variable[order(pct_concordant)]),
      concordance_lbl = sprintf("%.0f%%", pct_concordant)
    )
  
  ggplot(summary_df, aes(x = variable, y = pct_concordant, fill = pct_concordant)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    geom_text(aes(label = concordance_lbl, hjust = -0.08), size = 3.0, color = "grey20") +
    geom_hline(yintercept = 50, linetype = "dashed", color = "grey60", linewidth = 0.5) +
    coord_flip(clip = "off") +
    scale_fill_gradient2(
      low = "#e74c3c", mid = "#f39c12", high = "#27ae60",
      midpoint = 50, limits = c(0, 100), name = "% Concordante"
    ) +
    scale_y_continuous(limits = c(0, 115), labels = function(x) paste0(x, "%")) +
    labs(
      title    = "SHAP-LIME Agreement por Variavel (global)",
      subtitle = paste0(
        "Top ", top_n, " variaveis por |SHAP|  |  ",
        "linha tracejada = 50%  |  n = ", max(summary_df$n_pat), " pacientes"
      ),
      x = NULL, y = "% pacientes com SHAP e LIME concordantes"
    ) +
    theme_xai() +
    theme(legend.position = "right")
}

# ==============================================================================
# XGB TREE SHAP
# ==============================================================================
xgb_shap_per_patient <- function(model, x_mat, patient_ids, feature_names, ensemble_weight) {
  cat("    [XGB SHAP] Calculando tree SHAP...\n")
  contrib <- tryCatch(
    predict(model, newdata = x_mat, predcontrib = TRUE),
    error = function(e) { cat("    [!] predcontrib falhou:", e$message, "\n"); NULL }
  )
  if (is.null(contrib)) return(NULL)
  
  contrib_mat <- contrib[, seq_len(ncol(contrib) - 1), drop = FALSE]
  colnames(contrib_mat) <- feature_names[seq_len(ncol(contrib_mat))]
  
  lapply(seq_len(nrow(contrib_mat)), function(i) {
    safe_contrib_df(
      patient_id = patient_ids[i],
      variable   = colnames(contrib_mat),
      value      = as.numeric(contrib_mat[i, ]) * ensemble_weight,
      method     = "Tree SHAP (XGB)"
    )
  }) %>% setNames(patient_ids)
}

# ==============================================================================
# RSF SHAP VIA PERMUTAÇÃO
# ==============================================================================
rsf_shap_per_patient <- function(model, x_test, patient_ids, ensemble_weight) {
  cat("    [RSF SHAP] Calculando permutation SHAP...\n")
  n_pat      <- nrow(x_test)
  feat_names <- colnames(x_test)
  
  base_pred <- tryCatch(
    as.numeric(predict(model, newdata = x_test)$predicted),
    error = function(e) { cat("    [!] RSF predict base falhou:", e$message, "\n"); NULL }
  )
  if (is.null(base_pred)) return(NULL)
  
  shap_mat <- matrix(0, nrow = n_pat, ncol = length(feat_names),
                     dimnames = list(NULL, feat_names))
  
  for (j in seq_along(feat_names)) {
    x_perm      <- x_test
    x_perm[, j] <- sample(x_test[, j])
    perm_pred   <- tryCatch(
      as.numeric(predict(model, newdata = x_perm)$predicted),
      error = function(e) base_pred
    )
    shap_mat[, j] <- base_pred - perm_pred
  }
  
  lapply(seq_len(n_pat), function(i) {
    safe_contrib_df(
      patient_id = patient_ids[i],
      variable   = feat_names,
      value      = as.numeric(shap_mat[i, ]) * ensemble_weight,
      method     = "Permutation SHAP (RSF)"
    )
  }) %>% setNames(patient_ids)
}

# ==============================================================================
# LIME VIA PERTURBAÇÃO GAUSSIANA
# [FIX-C] glmnet ridge (alpha=0) — estável mesmo com colinearidade
# ==============================================================================
lime_per_patient <- function(predict_fn, x_test, patient_ids, ensemble_weight,
                             N = N_LIME_PERT, sigma = 1.0) {
  cat("    [LIME] Calculando LIME por perturbacao Gaussiana (glmnet ridge)...\n")
  
  if (!is.data.frame(x_test)) x_test <- as.data.frame(x_test)
  
  num_cols <- sapply(x_test, is.numeric)
  if (!all(num_cols)) {
    cat(sprintf("    [LIME] Removendo %d colunas nao-numericas\n", sum(!num_cols)))
    x_test <- x_test[, num_cols, drop = FALSE]
  }
  
  if (ncol(x_test) == 0) {
    cat("    [LIME] Nenhuma coluna numerica disponivel. Pulando.\n")
    return(setNames(vector("list", nrow(x_test)), patient_ids))
  }
  
  sd_x <- apply(x_test, 2, sd, na.rm = TRUE)
  sd_x[sd_x == 0 | is.na(sd_x)] <- 1
  
  lapply(seq_len(nrow(x_test)), function(i) {
    
    x0     <- as.numeric(x_test[i, ])
    noise  <- matrix(rnorm(N * ncol(x_test), 0, 0.1), nrow = N)
    X_pert <- sweep(noise, 2, x0, "+")
    colnames(X_pert) <- colnames(x_test)
    
    y_pert <- tryCatch(
      as.numeric(predict_fn(as.data.frame(X_pert))),
      error = function(e) rep(NA_real_, N)
    )
    
    dists  <- sqrt(rowSums(sweep(noise, 2, sd_x, "/")^2))
    kern_w <- exp(-dists^2 / (2 * sigma^2))
    
    ok <- is.finite(y_pert) & is.finite(kern_w)
    if (sum(ok) < 10) return(NULL)
    
    X_ok <- X_pert[ok, , drop = FALSE]
    y_ok <- y_pert[ok]
    w_ok <- kern_w[ok]
    
    coefs <- tryCatch({
      nfolds <- max(3L, min(5L, floor(sum(ok) / 10L)))
      cv_fit <- glmnet::cv.glmnet(
        x           = X_ok,
        y           = y_ok,
        weights     = w_ok,
        alpha       = 0,
        nfolds      = nfolds,
        standardize = TRUE
      )
      as.numeric(coef(cv_fit, s = "lambda.min"))[-1]
    }, error = function(e) {
      tryCatch({
        fit <- glmnet::glmnet(X_ok, y_ok, weights = w_ok, alpha = 0,
                              lambda = 1e-4, standardize = TRUE)
        as.numeric(coef(fit))[-1]
      }, error = function(e2) rep(NA_real_, ncol(X_ok)))
    })
    
    if (all(is.na(coefs))) return(NULL)
    
    if (is.null(names(coefs)) || length(names(coefs)) == 0)
      names(coefs) <- colnames(X_ok)
    
    coefs <- coefs[!is.na(coefs)]
    if (length(coefs) == 0) return(NULL)
    
    safe_contrib_df(
      patient_id = patient_ids[i],
      variable   = names(coefs),
      value      = as.numeric(coefs) * ensemble_weight,
      method     = "Local LIME surrogate"
    )
  }) %>% setNames(patient_ids)
}

# ==============================================================================
# SALVAR PDF via cairo_pdf
# ==============================================================================
save_pdf <- function(p, path, width = 11, height = 7) {
  if (!is.null(p)) {
    cairo_pdf(path, width = width, height = height)
    print(p)
    dev.off()
  }
}

# ==============================================================================
# AGREGAÇÃO SEGURA
# ==============================================================================
safe_bind_rows <- function(list_of_dfs) {
  cleaned <- lapply(Filter(Negate(is.null), list_of_dfs), function(df) {
    if (nrow(df) == 0) return(df)
    df %>% dplyr::mutate(dplyr::across(dplyr::everything(), function(x) {
      if (is.factor(x)) as.character(x) else x
    }))
  })
  if (length(cleaned) == 0) return(data.frame())
  bind_rows(cleaned)
}

# ==============================================================================
# [FIX-F] HELPER: detecção robusta de colunas em weights_df
# ==============================================================================
parse_weights_df <- function(weights_df) {
  
  cat("  [weights] colunas encontradas:", paste(names(weights_df), collapse = ", "), "\n")
  
  # --- Coluna de peso ---
  peso_col <- names(weights_df)[tolower(names(weights_df)) %in%
                                  c("peso", "weight", "coef", "coefficient", "lambda", "alpha")]
  if (length(peso_col) == 0) {
    num_cols <- names(weights_df)[sapply(weights_df, is.numeric)]
    if (length(num_cols) >= 1) {
      peso_col <- num_cols[1]
      cat("  [weights] coluna de peso nao identificada pelo nome; usando primeira numerica:", peso_col, "\n")
    } else {
      stop("Nao foi possivel identificar coluna de peso em weights_df. Colunas: ",
           paste(names(weights_df), collapse = ", "))
    }
  } else {
    peso_col <- peso_col[1]
  }
  weights_df[[peso_col]] <- suppressWarnings(as.numeric(weights_df[[peso_col]]))
  
  # --- Coluna de nome de modelo ---
  modelo_col <- names(weights_df)[tolower(names(weights_df)) %in%
                                    c("modelo", "model", "name", "nome", "id")]
  if (length(modelo_col) == 0) modelo_col <- names(weights_df)[1]
  modelo_col <- modelo_col[1]
  
  # --- Coluna de algoritmo ---
  alg_col <- names(weights_df)[tolower(names(weights_df)) %in%
                                 c("algoritmo", "algorithm", "alg", "type", "tipo", "method")]
  alg_col <- if (length(alg_col) > 0) alg_col[1] else NULL
  
  cat(sprintf("  [weights] mapeamento: modelo='%s', peso='%s', algoritmo='%s'\n",
              modelo_col, peso_col,
              if (is.null(alg_col)) "INFERIR DO NOME" else alg_col))
  
  # --- Filtrar peso != 0 ---
  selected <- weights_df[!is.na(weights_df[[peso_col]]) & weights_df[[peso_col]] != 0, , drop = FALSE]
  
  # --- Renomear para nomes padrao ---
  names(selected)[names(selected) == modelo_col] <- "Modelo"
  names(selected)[names(selected) == peso_col]   <- "Peso"
  if (!is.null(alg_col) && alg_col != "Modelo" && alg_col != "Peso") {
    names(selected)[names(selected) == alg_col] <- "Algoritmo"
  }
  
  # --- Inferir algoritmo se necessario ---
  if (!"Algoritmo" %in% names(selected)) {
    selected$Algoritmo <- dplyr::case_when(
      grepl("(?i)rsf|random.forest|rf_", selected$Modelo) ~ "RSF",
      grepl("(?i)xgb|xgboost|boost",    selected$Modelo) ~ "XGB",
      TRUE ~ NA_character_
    )
    n_na <- sum(is.na(selected$Algoritmo))
    if (n_na > 0)
      cat(sprintf("  [weights] %d modelo(s) sem Algoritmo identificavel pelo nome do modelo\n", n_na))
  }
  
  cat(sprintf("  [weights] %d modelos com peso != 0\n", nrow(selected)))
  
  if (nrow(selected) == 0) {
    cat("  [DIAGNOSTICO] Distribuicao de pesos encontrada:\n")
    print(table(weights_df[[peso_col]], useNA = "always"))
  }
  
  selected
}

# ==============================================================================
# PROCESS UNIT: CANCER x ENDPOINT x OMIC
# ==============================================================================
process_unit <- function(u) {
  cancer   <- u$Tipo
  endpoint <- u$Metrica
  omic     <- u$Omic
  key      <- unit_key(u)
  
  cp <- read_checkpoint()
  if (is_done(key, cp)) {
    cat(sprintf("  [SKIP - ja concluida] %s | %s | %s\n", cancer, endpoint, omic))
    return(invisible(NULL))
  }
  
  cat(sprintf("\n====== %s | %s | %s ======\n", cancer, endpoint, omic))
  
  result <- tryCatch({
    
    sl_dir <- file.path(SL_ROOT, cancer, endpoint, omic)
    
    unit_root     <- file.path(XAI_OUT_ROOT, cancer, endpoint, omic)
    ind_plots_dir <- file.path(unit_root, "individual",   "plots")
    ind_tabs_dir  <- file.path(unit_root, "individual",   "tables")
    pop_plots_dir <- file.path(unit_root, "populacional", "plots")
    pop_tabs_dir  <- file.path(unit_root, "populacional", "tables")
    for (d in c(ind_plots_dir, ind_tabs_dir, pop_plots_dir, pop_tabs_dir))
      if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    
    # ── 1. Carregar ensemble ───────────────────────────────────────────────
    sl_mod_path <- file.path(sl_dir, "ensemble_model_glmnet_train_only.rds")
    w_path      <- file.path(sl_dir, "ensemble_weights.csv")
    s_path      <- file.path(sl_dir, "Input_Model_Summary.csv")
    
    for (p in c(sl_mod_path, w_path, s_path)) {
      if (!file.exists(p)) stop("Arquivo nao encontrado: ", p)
    }
    
    sl_model   <- readRDS(sl_mod_path)
    weights_df <- readr::read_csv(w_path, show_col_types = FALSE)
    sl_summary <- readr::read_csv(s_path, show_col_types = FALSE)
    
    # [FIX-F] detecção robusta de colunas e filtragem
    selected <- parse_weights_df(weights_df)
    
    if (nrow(selected) == 0) {
      cat("  [AVISO] Nenhum modelo com peso != 0. Pulando...\n")
      stop("SKIP_ZERO")
    }
    
    # Remover modelos sem algoritmo identificado
    if (any(is.na(selected$Algoritmo))) {
      cat(sprintf("  [AVISO] Removendo %d modelos sem Algoritmo identificado\n",
                  sum(is.na(selected$Algoritmo))))
      selected <- selected[!is.na(selected$Algoritmo), , drop = FALSE]
    }
    if (nrow(selected) == 0)
      stop("Nenhum modelo com Algoritmo identificavel (RSF/XGB) apos filtragem.")
    
    cat(sprintf("  Modelos selecionados (%d):\n", nrow(selected)))
    print(selected[, c("Modelo", "Algoritmo", "Peso")])
    
    # ── 2. Loop por sub-modelo ─────────────────────────────────────────────
    shap_parts         <- list()
    lime_parts         <- list()
    global_patient_ids <- NULL
    x_all_collected    <- list()
    
    n_tentativas <- 0L
    n_sucesso    <- 0L
    
    for (i in seq_len(nrow(selected))) {
      m_name <- as.character(selected$Modelo[i])
      alg    <- as.character(selected$Algoritmo[i])
      peso   <- as.numeric(selected$Peso[i])
      mod_n  <- gsub(".*_M", "", m_name)
      
      n_tentativas <- n_tentativas + 1L
      cat(sprintf("  [mod %d/%d] %s (alg=%s, peso=%.4f)\n",
                  i, nrow(selected), m_name, alg, peso))
      
      if (alg == "RSF") {
        bp <- find_rsf_bundle(cancer, endpoint, omic, mod_n)
        
        # [FIX-G] diagnóstico quando bundle não encontrado
        if (is.null(bp)) {
          for (mode in c("multiomica", "mono_omica")) {
            tried <- file.path(
              RSF_ROOT, "output", "modelo", paste0("modelo", mod_n),
              mode, "multifenotipo", cancer, endpoint, omic, "Explicabilidade",
              paste0("XAI_bundle_", cancer, "_", endpoint, "_", omic, "_allPheno.rds")
            )
            cat("    [!] RSF tentado:", tried, "| existe:", file.exists(tried), "\n")
          }
          next
        }
        
        b           <- readRDS(bp)
        rsf_model   <- b$model
        x_xai       <- b$x_test_xai
        patient_ids <- as.character(b$patient_ids_test)
        if (is.null(global_patient_ids)) global_patient_ids <- patient_ids
        
        if (length(x_all_collected) == 0 && is.data.frame(x_xai)) {
          x_df_raw <- as.data.frame(x_xai)
          rownames(x_df_raw) <- patient_ids
          x_all_collected[["rsf"]] <- x_df_raw
        }
        
        shap_r <- rsf_shap_per_patient(rsf_model, x_xai, patient_ids, peso)
        if (!is.null(shap_r)) {
          shap_parts[[m_name]] <- shap_r
          n_sucesso <- n_sucesso + 1L
        }
        
        rsf_predict_fn <- function(nd) as.numeric(predict(rsf_model, newdata = nd)$predicted)
        lime_r <- lime_per_patient(rsf_predict_fn, x_xai, patient_ids, peso)
        if (!is.null(lime_r)) lime_parts[[m_name]] <- lime_r
        
      } else if (alg == "XGB") {
        xp <- find_xgb_xai(cancer, endpoint, omic, mod_n)
        
        # [FIX-G] diagnóstico quando XAI não encontrado
        if (is.null(xp)) {
          for (mode in c("multiomica", "mono_omica")) {
            base <- file.path(
              XGB_ROOT, "output", "modelo", paste0("modelo", mod_n),
              mode, "multifenotipo", cancer, endpoint, omic, "XAI_ready"
            )
            m_try <- file.path(base, paste0("xgb_model_", cancer, "_", endpoint, "_", omic, ".rds"))
            d_try <- file.path(base, paste0("data_xai_",  cancer, "_", endpoint, "_", omic, ".rds"))
            cat("    [!] XGB base:", base, "\n")
            cat("        model existe:", file.exists(m_try),
                "| data existe:", file.exists(d_try), "\n")
          }
          next
        }
        
        xgb_mod     <- readRDS(xp$model_path)
        xgb_data    <- readRDS(xp$data_path)
        x_test_sp   <- xgb_data$x_test
        feat_names  <- xgb_data$feature_names
        patient_ids <- as.character(xgb_data$patient_ids_test)
        df_test     <- xgb_data$df_test
        if (is.null(global_patient_ids)) global_patient_ids <- patient_ids
        
        if (length(x_all_collected) == 0) {
          x_df_raw <- df_test %>% dplyr::select(-dplyr::any_of(c("time", "status")))
          rownames(x_df_raw) <- patient_ids
          x_all_collected[["xgb"]] <- x_df_raw
        }
        
        shap_x <- xgb_shap_per_patient(xgb_mod, x_test_sp, patient_ids, feat_names, peso)
        if (!is.null(shap_x)) {
          shap_parts[[m_name]] <- shap_x
          n_sucesso <- n_sucesso + 1L
        }
        
        x_df <- df_test %>% dplyr::select(-dplyr::any_of(c("time", "status")))
        xgb_predict_fn <- function(nd) {
          mm    <- Matrix::sparse.model.matrix(~ . - 1, data = nd)
          X_out <- matrix(0, nrow(mm), length(feat_names),
                          dimnames = list(NULL, feat_names))
          cn    <- intersect(colnames(mm), feat_names)
          X_out[, cn] <- as.matrix(mm[, cn])
          predict(xgb_mod, newdata = as(X_out, "dgCMatrix"))
        }
        lime_x <- lime_per_patient(xgb_predict_fn, x_df, patient_ids, peso)
        if (!is.null(lime_x)) lime_parts[[m_name]] <- lime_x
        
      } else {
        cat(sprintf("    [!] Algoritmo desconhecido '%s' — pulando\n", alg))
      }
    }
    
    cat(sprintf("  Sub-modelos: %d tentados | %d com SHAP gerado\n",
                n_tentativas, n_sucesso))
    
    if (length(shap_parts) == 0 && length(lime_parts) == 0) {
      cat("  [AVISO] Nenhum resultado SHAP/LIME gerado (modelos/arquivos ausentes). Pulando...\n")
      stop("SKIP_NO_XAI")
    }
    
    # ── 3. Agregar contribuições ───────────────────────────────────────────
    agg_fn <- function(parts) {
      rows <- safe_bind_rows(lapply(parts, function(pl) {
        safe_bind_rows(Filter(Negate(is.null), pl))
      }))
      if (nrow(rows) == 0) return(rows)
      rows %>%
        dplyr::group_by(patient_id, variable) %>%
        dplyr::summarise(
          value  = sum(value, na.rm = TRUE),
          method = dplyr::first(method),
          .groups = "drop"
        )
    }
    
    all_shap_df <- agg_fn(shap_parts)
    all_lime_df <- agg_fn(lime_parts)
    x_all_df    <- if (length(x_all_collected) > 0) x_all_collected[[1]] else NULL
    
    # [FIX-E] Garantir patient_id como character ANTES de qualquer filtro/loop
    if (nrow(all_shap_df) > 0) all_shap_df$patient_id <- as.character(all_shap_df$patient_id)
    if (nrow(all_lime_df)  > 0) all_lime_df$patient_id  <- as.character(all_lime_df$patient_id)
    
    # ── 4. Tabela de agreement ─────────────────────────────────────────────
    agreement_df <- NULL
    if (nrow(all_shap_df) > 0 && nrow(all_lime_df) > 0) {
      agreement_df <- dplyr::inner_join(
        all_shap_df %>% dplyr::select(patient_id, variable, shap = value),
        all_lime_df  %>% dplyr::select(patient_id, variable, lime = value),
        by = c("patient_id", "variable")
      ) %>% dplyr::mutate(concordant = sign(shap) == sign(lime))
    }
    
    # ── 5. Tabelas populacionais ───────────────────────────────────────────
    if (nrow(all_shap_df) > 0) {
      readr::write_csv(all_shap_df, file.path(pop_tabs_dir, "SHAP_all_patients.csv"))
      all_shap_df %>%
        dplyr::group_by(variable) %>%
        dplyr::summarise(mean_abs_shap = mean(abs(value), na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(desc(mean_abs_shap)) %>%
        readr::write_csv(file.path(pop_tabs_dir, "Global_SHAP_importance.csv"))
    }
    if (nrow(all_lime_df) > 0) {
      readr::write_csv(all_lime_df, file.path(pop_tabs_dir, "LIME_all_patients.csv"))
      all_lime_df %>%
        dplyr::group_by(variable) %>%
        dplyr::summarise(mean_abs_lime = mean(abs(value), na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(desc(mean_abs_lime)) %>%
        readr::write_csv(file.path(pop_tabs_dir, "Global_LIME_importance.csv"))
    }
    if (!is.null(agreement_df))
      readr::write_csv(agreement_df,
                       file.path(pop_tabs_dir, "SHAP_LIME_agreement_all_patients.csv"))
    
    # ── 6. Plots POPULACIONAIS ─────────────────────────────────────────────
    cat("  Gerando plots populacionais...\n")
    
    if (nrow(all_shap_df) > 0) {
      save_pdf(plot_beeswarm(all_shap_df, "SHAP", x_all_df),
               file.path(pop_plots_dir, "SHAP_beeswarm.pdf"), 13, 11)
      save_pdf(plot_heatmap(all_shap_df, "SHAP"),
               file.path(pop_plots_dir, "SHAP_heatmap_all_patients.pdf"), 16, 10)
    }
    if (nrow(all_lime_df) > 0) {
      save_pdf(plot_beeswarm(all_lime_df, "LIME", x_all_df),
               file.path(pop_plots_dir, "LIME_beeswarm.pdf"), 13, 11)
      save_pdf(plot_heatmap(all_lime_df, "LIME"),
               file.path(pop_plots_dir, "LIME_heatmap_all_patients.pdf"), 16, 10)
    }
    if (!is.null(agreement_df) && nrow(agreement_df) > 0) {
      save_pdf(plot_agreement_summary(all_shap_df, all_lime_df),
               file.path(pop_plots_dir, "SHAP_LIME_agreement_summary.pdf"), 10, 8)
    }
    cat("  Plots populacionais concluidos.\n")
    
    # ── 7. Plots INDIVIDUAIS ───────────────────────────────────────────────
    pids <- if (!is.null(global_patient_ids)) global_patient_ids else
      unique(c(all_shap_df$patient_id, all_lime_df$patient_id))
    
    # [FIX-E] garantir character
    pids <- as.character(pids)
    
    cat(sprintf("  Gerando plots individuais (%d pacientes)...\n", length(pids)))
    
    for (pid in pids) {
      safe_pid <- gsub("[^A-Za-z0-9._-]", "_", pid)
      pid_chr  <- as.character(pid)   # [FIX-E] tipo seguro para comparação
      
      x_pat <- if (!is.null(x_all_df) && pid_chr %in% rownames(x_all_df)) {
        as.list(x_all_df[pid_chr, , drop = FALSE])
      } else NULL
      
      # [FIX-E] indexação base R evita conflito de tipo em dplyr::filter
      df_s <- all_shap_df[all_shap_df$patient_id == pid_chr, , drop = FALSE]
      df_l <- all_lime_df[all_lime_df$patient_id  == pid_chr, , drop = FALSE]
      
      if (nrow(df_s) > 0) {
        p_shap <- plot_waterfall(df_s, pid_chr, cancer, endpoint,
                                 "Tree/Permutation SHAP", x_patient = x_pat)
        save_pdf(p_shap, file.path(ind_plots_dir, paste0("SHAP_", safe_pid, ".pdf")))
        readr::write_csv(df_s, file.path(ind_tabs_dir, paste0("SHAP_", safe_pid, ".csv")))
      }
      if (nrow(df_l) > 0) {
        p_lime <- plot_waterfall(df_l, pid_chr, cancer, endpoint,
                                 "Local LIME surrogate", x_patient = x_pat)
        save_pdf(p_lime, file.path(ind_plots_dir, paste0("LIME_", safe_pid, ".pdf")))
        readr::write_csv(df_l, file.path(ind_tabs_dir, paste0("LIME_", safe_pid, ".csv")))
      }
      if (nrow(df_s) > 0 && nrow(df_l) > 0) {
        p_agree <- plot_shap_lime_agreement(df_s, df_l, pid_chr, cancer, endpoint)
        save_pdf(p_agree,
                 file.path(ind_plots_dir, paste0("SHAP_vs_LIME_", safe_pid, ".pdf")))
      }
    }
    
    cat(sprintf("  [ok] %d pacientes | %s/%s/%s\n", length(pids), cancer, endpoint, omic))
    "done"
    
  }, error = function(e) {
    if (e$message == "SKIP_ZERO") return("skipped_zero_weights")
    if (e$message == "SKIP_NO_XAI") return("skipped_no_xai")
    
    cat(sprintf("  [ERRO] %s/%s/%s: %s\n", cancer, endpoint, omic, e$message))
    list(status = "error", msg = e$message)
  })
  
  if (identical(result, "done") || identical(result, "skipped_zero_weights") || identical(result, "skipped_no_xai")) {
    write_checkpoint_entry(key, status = "done")
    
    if (identical(result, "done")) {
      cat(sprintf("  [checkpoint] done gravado: %s\n", key))
    } else {
      cat(sprintf("  [checkpoint] done (skipped) gravado: %s\n", key))
    }
  } else {
    err_msg <- if (is.list(result)) result$msg else as.character(result)
    write_checkpoint_entry(key, status = "error", error_msg = err_msg)
    cat(sprintf("  [checkpoint] error gravado: %s\n", key))
  }
  
  invisible(NULL)
}

# ==============================================================================
# EXECUÇÃO PRINCIPAL
# ==============================================================================
cat("Inicio:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Checkpoint XAI:", XAI_CHECKPOINT_PATH, "\n\n")

if (!file.exists(PROGRESS_PATH))
  stop("[ERRO] Checkpoint de entrada nao encontrado: ", PROGRESS_PATH)

master_list <- readRDS(PROGRESS_PATH)
master_df   <- do.call(rbind, master_list)

# Filtrar para executar apenas multiomica (Omic == "all")
master_df <- subset(master_df, Omic == "all")

cp_atual <- read_checkpoint()
summarise_checkpoint(cp_atual)

units_all  <- split(master_df, seq(nrow(master_df)))
units_todo <- Filter(function(u) !is_done(unit_key(u), cp_atual), units_all)

actual_workers <- if (.Platform$OS.type == "windows") 1 else N_WORKERS

cat(sprintf(
  "Unidades totais: %d | Ja concluidas: %d | A processar: %d | Workers: %d (solicitados: %d)\n\n",
  length(units_all),
  length(units_all) - length(units_todo),
  length(units_todo),
  actual_workers,
  N_WORKERS
))

if (length(units_todo) == 0) {
  cat("Todas as unidades ja foram concluidas. Nada a fazer.\n")
} else {
  results <- mclapply(units_todo, process_unit, mc.cores = actual_workers)
}

cp_final <- read_checkpoint()
cat("\n======== RESUMO FINAL ========\n")
summarise_checkpoint(cp_final)

if (length(cp_final) > 0) {
  erros <- Filter(function(x) x$status == "error", cp_final)
  if (length(erros) > 0) {
    cat("\nUnidades com erro:\n")
    for (e in erros)
      cat(sprintf("  - %s: %s\n", e$key, e$error_msg))
  }
}

cat("\nFim:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Resultados em:", XAI_OUT_ROOT, "\n")
cat("Checkpoint em:", XAI_CHECKPOINT_PATH, "\n")

# ==============================================================================
# RESUMO DOS FIXES V8
# ==============================================================================
# [FIX-E] In argument: patient_id == pid
#   CAUSA RAIZ: patient_id pode chegar como factor em all_shap_df/all_lime_df
#   após aggregação, causando falha silenciosa no dplyr::filter com pid character.
#   SOLUÇÃO: as.character() forçado nos data.frames após agg_fn(), antes do loop;
#   uso de indexação base R (df[df$col == val]) dentro do loop de pacientes,
#   que é imune a conflitos de tipo e independe do namespace dplyr.
#
# [FIX-F] Nenhum modelo com peso != 0
#   CAUSA RAIZ: o nome da coluna de peso no CSV pode variar (Peso/Weight/coef),
#   o case pode diferir, ou os valores podem ser NA/character ao ser lido.
#   SOLUÇÃO: função parse_weights_df() que detecta automaticamente as colunas
#   por lista de nomes alternativos (case-insensitive), converte para numérico,
#   e emite diagnóstico completo quando falha. Detecta também Modelo e Algoritmo
#   de forma robusta, e infere Algoritmo pelo padrão do nome quando coluna ausente.
#
# [FIX-G] Nenhum resultado SHAP/LIME gerado
#   CAUSA RAIZ: todos os sub-modelos caíam em `next` silenciosamente por bundle
#   ou XAI não encontrado, sem mostrar qual caminho foi tentado.
#   SOLUÇÃO: exibição de todos os caminhos tentados (com file.exists()) quando
#   bundle/XAI não encontrado; contador de tentativas vs sucessos por loop.
#
# Para reprocessar unidades com erro:
#   cp <- readRDS("~/Emanuell/super_learning/XAI_SurvSHAP_SurvLIME/XAI_checkpoint_progress.rds")
#   erros <- names(Filter(function(x) x$status == "error", cp))
#   for (e in erros) cp[[e]] <- NULL
#   saveRDS(cp, "~/Emanuell/super_learning/XAI_SurvSHAP_SurvLIME/XAI_checkpoint_progress.rds")
#
# Para reprocessar uma unidade específica:
#   cp[["BRCA|OS|all"]] <- NULL
#   saveRDS(cp, "~/Emanuell/super_learning/XAI_SurvSHAP_SurvLIME/XAI_checkpoint_progress.rds")
#
# Para resetar tudo:
#   file.remove("~/Emanuell/super_learning/XAI_SurvSHAP_SurvLIME/XAI_checkpoint_progress.rds")
# ==============================================================================