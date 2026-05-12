# ======================================================================================
# ANÁLISE DE SENSIBILIDADE A FÁRMACOS POR CONTEXTO TUMORAL COM DEPMAP
# + Integração com secondary-screen-dose-response-curve-parameters.csv
# + MATCH EXPLÍCITO ENTRE TARGET DO FÁRMACO E GENE DA ASSINATURA
#
# Objetivo:
# Ler assinaturas da coluna Signature, manter a auditoria de componentes, e analisar
# sensibilidade a drogas nas linhagens celulares correspondentes a cada tipo tumoral
# TCGA, usando:
#   1) depmap_drug_sensitivity() + depmap_metadata()
#   2) secondary-screen-dose-response-curve-parameters.csv
#
# Regra central desta versão:
# Uma droga só será associada a uma assinatura se pelo menos um alvo anotado
# do fármaco coincidir com pelo menos um gene da assinatura.
#
# Entrada:
# G:\DepMap\drug_sensitivity\Assinaturas_Omicas_Uma_Por_Linha_Anotadas.csv
# G:\DepMap\drug_sensitivity\secondary-screen-dose-response-curve-parameters.csv
#
# Saída:
# G:\DepMap\drug_sensitivity
#
# Organização das saídas:
# - Resultados gerais: output_dir
# - Resultados estratificados por tipo tumoral: output_dir/type/<TIPO>/
#
# Interpretação:
# dependency mais negativo -> maior sensibilidade relativa ao composto
# IC50/EC50/AUC/R2 -> parâmetros farmacológicos complementares
# ======================================================================================

suppressPackageStartupMessages({
  library(depmap)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
})

# ======================================================================================
# 1) CONFIGURAÇÃO
# ======================================================================================

arquivo_assinaturas  <- "D:/DepMap/drug_sensitivity/Assinaturas_Omicas_SuperLearner.tsv"
arquivo_curve_params <- "D:/DepMap/drug_sensitivity/secondary-screen-dose-response-curve-parameters.csv"

output_dir <- "D:/DepMap/drug_sensitivity"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

type_dir <- file.path(output_dir, "type")
dir.create(type_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(output_dir, "analysis_log.txt")

max_linhas_por_arquivo <- 2000000L
max_assinaturas_por_bloco <- 25L

log_message <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  write(msg, file = log_file, append = TRUE)
}

if (file.exists(log_file)) file.remove(log_file)
file.create(log_file)

log_message("============================================================")
log_message("INÍCIO DA ANÁLISE DE DRUG SENSITIVITY NO DEPMAP")
log_message("Arquivo de assinaturas: ", arquivo_assinaturas)
log_message("Arquivo de curvas dose-resposta: ", arquivo_curve_params)
log_message("Diretório de saída geral: ", output_dir)
log_message("Diretório de saída por tipo: ", type_dir)
log_message("max_linhas_por_arquivo = ", max_linhas_por_arquivo)
log_message("max_assinaturas_por_bloco = ", max_assinaturas_por_bloco)
log_message("============================================================")

# ======================================================================================
# 2) FUNÇÕES AUXILIARES
# ======================================================================================

limpar_assinatura <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("^\\s*\\(", "") %>%
    str_replace_all("\\)\\s*$", "") %>%
    str_squish()
}

quebrar_componentes_assinatura <- function(x) {
  x_limpo <- limpar_assinatura(x)
  
  if (is.na(x_limpo) || x_limpo == "") {
    return(character(0))
  }
  
  x_limpo %>%
    str_split("\\s*\\+\\s*", simplify = FALSE) %>%
    .[[1]] %>%
    str_trim() %>%
    .[. != "" & !is.na(.)]
}

classificar_componente <- function(x) {
  x2 <- toupper(as.character(x))
  
  case_when(
    is.na(x2) | x2 == "" ~ "vazio",
    str_detect(x2, "^ENST[0-9]+") ~ "transcrito_ensembl",
    str_detect(x2, "^ENSG[0-9]+") ~ "gene_ensembl",
    str_detect(x2, "^NM_[0-9]+") ~ "transcrito_refseq_mrna",
    str_detect(x2, "^NR_[0-9]+") ~ "transcrito_refseq_ncrna",
    str_detect(x2, "^XM_[0-9]+") ~ "transcrito_modelo_refseq",
    str_detect(x2, "^XR_[0-9]+") ~ "transcrito_modelo_refseq",
    str_detect(x2, "ISOFORM|TRANSCRIPT") ~ "transcrito_textual",
    str_detect(x2, "-20[0-9]$|-2[0-9][0-9]$") ~ "possivel_isoforma",
    str_detect(x2, "\\.[0-9]+$") & str_detect(x2, "^ENS[A-Z]*[0-9]+") ~ "identificador_ensembl_versionado",
    str_detect(x2, "^[A-Z0-9._-]+$") ~ "candidato_gene_symbol",
    TRUE ~ "ambiguo"
  )
}

normalizar_token_alvo <- function(x) {
  x %>%
    as.character() %>%
    toupper() %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

quebrar_targets_farmaco <- function(x) {
  x_limpo <- normalizar_token_alvo(x)
  
  if (is.na(x_limpo) || x_limpo == "") {
    return(character(0))
  }
  
  partes <- x_limpo %>%
    str_split("\\s*[,;|/+&]\\s*|\\s+AND\\s+|\\s+OR\\s+", simplify = FALSE) %>%
    .[[1]] %>%
    str_trim()
  
  partes <- partes[!is.na(partes) & partes != ""]
  
  unique(c(x_limpo, partes))
}

interpretar_dependency_farmaco <- function(x) {
  case_when(
    is.na(x) ~ "Sem dado disponível",
    x < -1.5 ~ "Alta sensibilidade média ao composto",
    x < -0.5 ~ "Sensibilidade moderada ao composto",
    x <= 0.5 ~ "Pouco efeito médio do composto",
    x > 0.5 ~ "Baixa sensibilidade relativa / possível resistência",
    TRUE ~ "Indeterminado"
  )
}

interpretar_ic50 <- function(x) {
  case_when(
    is.na(x) ~ "Sem IC50 disponível",
    x < 0.1 ~ "Potência muito alta",
    x < 1   ~ "Potência alta",
    x < 10  ~ "Potência moderada",
    x >= 10 ~ "Potência baixa",
    TRUE ~ "Indeterminado"
  )
}

interpretar_r2 <- function(x) {
  case_when(
    is.na(x) ~ "Sem informação de ajuste",
    x >= 0.90 ~ "Ajuste excelente",
    x >= 0.75 ~ "Ajuste bom",
    x >= 0.50 ~ "Ajuste moderado",
    x < 0.50  ~ "Ajuste fraco",
    TRUE ~ "Indeterminado"
  )
}

safe_mean <- function(x) {
  x_ok <- x[!is.na(x)]
  if (length(x_ok) == 0) return(NA_real_)
  mean(x_ok)
}

safe_median <- function(x) {
  x_ok <- x[!is.na(x)]
  if (length(x_ok) == 0) return(NA_real_)
  median(x_ok)
}

safe_sd <- function(x) {
  x_ok <- x[!is.na(x)]
  if (length(x_ok) <= 1) return(NA_real_)
  sd(x_ok)
}

safe_min <- function(x) {
  x_ok <- x[!is.na(x)]
  if (length(x_ok) == 0) return(NA_real_)
  min(x_ok)
}

safe_max <- function(x) {
  x_ok <- x[!is.na(x)]
  if (length(x_ok) == 0) return(NA_real_)
  max(x_ok)
}

safe_log10 <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(x_num) & x_num > 0, log10(x_num), NA_real_)
}

safe_prop_true <- function(x) {
  x2 <- toupper(trimws(as.character(x)))
  x2 <- x2[!is.na(x2) & x2 != ""]
  if (length(x2) == 0) return(NA_real_)
  mean(x2 %in% c("TRUE", "T", "PASS", "PASSED", "YES", "Y", "1"))
}

safe_max_int <- function(x) {
  x_ok <- x[!is.na(x)]
  if (length(x_ok) == 0) return(NA_integer_)
  as.integer(max(x_ok))
}

collapse_unique_plus <- function(x) {
  x2 <- unique(na.omit(as.character(x)))
  x2 <- x2[x2 != ""]
  if (length(x2) == 0) return(NA_character_)
  paste(sort(x2), collapse = " + ")
}

collapse_unique_bar <- function(x) {
  x2 <- unique(na.omit(as.character(x)))
  x2 <- x2[x2 != ""]
  if (length(x2) == 0) return(NA_character_)
  paste(sort(x2), collapse = " | ")
}

calcular_n_partes <- function(n_linhas_estimadas, max_linhas_por_arquivo = 2000000L) {
  if (is.na(n_linhas_estimadas) || n_linhas_estimadas <= 0) {
    return(1L)
  }
  as.integer(max(1, ceiling(n_linhas_estimadas / max_linhas_por_arquivo)))
}

dividir_indices_em_partes <- function(n, n_partes) {
  if (n <= 0 || n_partes <= 0) {
    return(list(integer(0)))
  }
  
  ids <- seq_len(n)
  split(ids, ceiling(seq_along(ids) / ceiling(n / n_partes)))
}

exportar_em_partes <- function(df, output_prefix, max_linhas_por_arquivo = 2000000L) {
  n_total <- nrow(df)
  
  if (n_total == 0) {
    return(
      tibble(
        arquivo = NA_character_,
        parte = NA_integer_,
        n_linhas = 0L
      )
    )
  }
  
  n_partes <- calcular_n_partes(n_total, max_linhas_por_arquivo)
  partes_idx <- dividir_indices_em_partes(n_total, n_partes)
  
  mapa_partes <- purrr::imap_dfr(
    partes_idx,
    function(idx, parte_i) {
      arquivo_i <- paste0(output_prefix, "_Parte_", parte_i, ".csv")
      readr::write_csv(df[idx, , drop = FALSE], arquivo_i)
      
      tibble(
        arquivo = arquivo_i,
        parte = as.integer(parte_i),
        n_linhas = length(idx)
      )
    }
  )
  
  mapa_partes
}

# ======================================================================================
# 3) MAPA EXPLÍCITO TCGA -> DEPMAP
# ======================================================================================

criar_mapa_tcga_depmap <- function(tipos_unicos) {
  
  mapa_manual <- tribble(
    ~Tipo,   ~nivel_filtro,          ~lineage_alvo,               ~lineage_subtype_alvo,              ~lineage_sub_subtype_alvo,          ~observacao,
    "BLCA",  "lineage_subtype",      "urinary_tract",            "bladder_carcinoma",                NA_character_,                      "Bexiga; usar bladder_carcinoma",
    "BRCA",  "lineage",              "breast",                   NA_character_,                      NA_character_,                      "Mama; manter abrangente por heterogeneidade",
    "CESC",  "custom",               "cervix",                   NA_character_,                      NA_character_,                      "Cervix; incluir adenocarcinoma e squamous",
    "COAD",  "lineage_subtype",      "colorectal",               "colorectal_adenocarcinoma",        NA_character_,                      "Cólon; adenocarcinoma colorretal",
    "GBM",   "lineage_sub_subtype",  "central_nervous_system",   "glioma",                           "glioblastoma",                     "GBM específico",
    "HNSC",  "lineage_subtype",      "upper_aerodigestive",      "upper_aerodigestive_squamous",     NA_character_,                      "Cabeça e pescoço escamoso",
    "KIRC",  "lineage_sub_subtype",  "kidney",                   "renal_cell_carcinoma",             "clear_cell",                       "RCC clear cell",
    "LGG",   "custom",               "central_nervous_system",   "glioma",                           NA_character_,                      "LGG = astrocytoma/oligodendroglioma; excluir glioblastoma",
    "LUAD",  "lineage_sub_subtype",  "lung",                     "NSCLC",                            "NSCLC_adenocarcinoma",             "Adenocarcinoma de pulmão",
    "LUSC",  "lineage_sub_subtype",  "lung",                     "NSCLC",                            "NSCLC_squamous",                   "Escamoso de pulmão",
    "PRAD",  "lineage",              "prostate",                 NA_character_,                      NA_character_,                      "Próstata",
    "SKCM",  "lineage_subtype",      "skin",                     "melanoma",                         NA_character_,                      "Melanoma; manter todos os subsubtipos",
    "STAD",  "lineage_subtype",      "gastric",                  "gastric_adenocarcinoma",           NA_character_,                      "Adenocarcinoma gástrico",
    "THCA",  "lineage_subtype",      "thyroid",                  "thyroid_carcinoma",                NA_character_,                      "Carcinoma de tireoide",
    "THYM",  "sem_match",            NA_character_,              NA_character_,                      NA_character_,                      "Sem match claro no recorte atual do DepMap"
  )
  
  tibble(Tipo = tipos_unicos) %>%
    left_join(mapa_manual, by = "Tipo") %>%
    mutate(
      nivel_filtro = coalesce(nivel_filtro, "sem_match"),
      observacao = coalesce(observacao, "Tipo não mapeado manualmente")
    )
}

# ======================================================================================
# 4) FUNÇÃO DE FILTRAGEM POR TIPO TCGA
# ======================================================================================

filtrar_depmap_por_tipo_tcga <- function(df, tipo_tcga, mapa_tcga_depmap) {
  
  info <- mapa_tcga_depmap %>%
    filter(Tipo == tipo_tcga)
  
  if (nrow(info) == 0) {
    return(tibble())
  }
  
  nivel <- info$nivel_filtro[[1]]
  lin <- info$lineage_alvo[[1]]
  sub <- info$lineage_subtype_alvo[[1]]
  subsub <- info$lineage_sub_subtype_alvo[[1]]
  
  if (is.na(nivel) || nivel == "sem_match") {
    return(tibble())
  }
  
  if (nivel == "lineage") {
    return(df %>% filter(lineage == lin))
  }
  
  if (nivel == "lineage_subtype") {
    return(df %>% filter(lineage == lin, lineage_subtype == sub))
  }
  
  if (nivel == "lineage_sub_subtype") {
    return(df %>% filter(lineage == lin, lineage_subtype == sub, lineage_sub_subtype == subsub))
  }
  
  if (nivel == "custom" && tipo_tcga == "LGG") {
    return(
      df %>%
        filter(
          lineage == "central_nervous_system",
          lineage_subtype == "glioma",
          lineage_sub_subtype %in% c("astrocytoma", "oligodendroglioma")
        )
    )
  }
  
  if (nivel == "custom" && tipo_tcga == "CESC") {
    return(
      df %>%
        filter(
          lineage == "cervix",
          lineage_subtype %in% c("cervical_adenocarcinoma", "cervical_carcinoma", "cervical_squamous")
        )
    )
  }
  
  tibble()
}

# ======================================================================================
# 5) FUNÇÃO DE AUDITORIA DO CONTEXTO TUMORAL
# ======================================================================================

auditar_contexto_tcga_depmap <- function(df_base, tipo_tcga, mapa_tcga_depmap) {
  
  info <- mapa_tcga_depmap %>%
    filter(Tipo == tipo_tcga)
  
  if (nrow(info) == 0) {
    return(tibble(
      Tipo = tipo_tcga,
      nivel_filtro = "nao_mapeado",
      lineage_alvo = NA_character_,
      lineage_subtype_alvo = NA_character_,
      lineage_sub_subtype_alvo = NA_character_,
      observacao = "Tipo não presente no mapa",
      n_linhas_filtradas = 0L,
      n_depmap_ids = 0L,
      n_cell_lines = 0L
    ))
  }
  
  filtrado <- filtrar_depmap_por_tipo_tcga(df_base, tipo_tcga, mapa_tcga_depmap)
  
  tibble(
    Tipo = tipo_tcga,
    nivel_filtro = info$nivel_filtro[[1]],
    lineage_alvo = info$lineage_alvo[[1]],
    lineage_subtype_alvo = info$lineage_subtype_alvo[[1]],
    lineage_sub_subtype_alvo = info$lineage_sub_subtype_alvo[[1]],
    observacao = info$observacao[[1]],
    n_linhas_filtradas = nrow(filtrado),
    n_depmap_ids = n_distinct(filtrado$depmap_id),
    n_cell_lines = n_distinct(filtrado$cell_line)
  )
}

# ======================================================================================
# 6) LEITURA DA TABELA DE ASSINATURAS
# ======================================================================================

log_message("Lendo arquivo de assinaturas...")
assinaturas_raw <- read_table(arquivo_assinaturas, show_col_types = FALSE)

colunas_obrigatorias <- c("Signature", "Tipo")
faltantes_entrada <- setdiff(colunas_obrigatorias, names(assinaturas_raw))

if (length(faltantes_entrada) > 0) {
  stop("Faltam colunas obrigatórias no arquivo de assinaturas: ",
       paste(faltantes_entrada, collapse = ", "))
}

log_message("Dimensão assinaturas_raw: ", nrow(assinaturas_raw), " x ", ncol(assinaturas_raw))

assinaturas_tbl <- assinaturas_raw %>%
  mutate(
    row_id = row_number(),
    Signature = as.character(Signature),
    Tipo = as.character(Tipo),
    Signature_limpa = limpar_assinatura(Signature)
  )

# ======================================================================================
# 6B) AUDITORIA DE UNICIDADE DAS ASSINATURAS POR TIPO
# ======================================================================================

auditoria_unicidade_assinaturas <- assinaturas_tbl %>%
  group_by(Tipo) %>%
  summarise(
    n_linhas = n(),
    n_signatures_unicas = n_distinct(Signature),
    n_signatures_limpas_unicas = n_distinct(Signature_limpa),
    ha_duplicatas_signature = n_linhas > n_signatures_unicas,
    ha_duplicatas_signature_limpa = n_linhas > n_signatures_limpas_unicas,
    .groups = "drop"
  ) %>%
  arrange(Tipo)

write_csv(
  auditoria_unicidade_assinaturas,
  file.path(output_dir, "auditoria_unicidade_assinaturas_por_tipo.csv")
)

log_message("Auditoria de unicidade das assinaturas por tipo exportada.")

# ======================================================================================
# 7) EXPANSÃO DAS ASSINATURAS EM COMPONENTES
# ======================================================================================

log_message("Expandindo assinaturas em componentes...")

assinaturas_componentes <- assinaturas_tbl %>%
  mutate(componentes = map(Signature_limpa, quebrar_componentes_assinatura)) %>%
  tidyr::unnest_longer(componentes, values_to = "componente", keep_empty = TRUE) %>%
  mutate(
    componente = str_trim(as.character(componente)),
    componente_upper = toupper(componente),
    tipo_componente = classificar_componente(componente_upper),
    analisavel_para_biomarcador = tipo_componente %in% c("candidato_gene_symbol", "gene_ensembl")
  )

log_message("Número de componentes expandidos: ", nrow(assinaturas_componentes))

# ======================================================================================
# 8) AUDITORIA DOS COMPONENTES
# ======================================================================================

auditoria_componentes <- assinaturas_componentes %>%
  count(tipo_componente, analisavel_para_biomarcador, name = "n_componentes") %>%
  arrange(desc(n_componentes))

write_csv(
  auditoria_componentes,
  file.path(output_dir, "auditoria_componentes_assinaturas.csv")
)

componentes_nao_analisaveis <- assinaturas_componentes %>%
  filter(!analisavel_para_biomarcador | is.na(analisavel_para_biomarcador)) %>%
  select(row_id, Tipo, Signature, componente, tipo_componente) %>%
  distinct()

write_csv(
  componentes_nao_analisaveis,
  file.path(output_dir, "componentes_nao_analisaveis_para_biomarcador.csv")
)

log_message("Auditoria de componentes exportada.")

# ======================================================================================
# 9) CARREGAMENTO DOS DADOS DEPMAP E CURVAS
# ======================================================================================

log_message("Carregando depmap_drug_sensitivity()...")
data_drug <- depmap::depmap_drug_sensitivity()

log_message("Carregando depmap_metadata()...")
data_meta <- depmap::depmap_metadata()

log_message("Lendo secondary-screen-dose-response-curve-parameters.csv...")
curve_raw <- read_csv(arquivo_curve_params, show_col_types = FALSE)

log_message("Dimensão data_drug: ", nrow(data_drug), " x ", ncol(data_drug))
log_message("Dimensão data_meta: ", nrow(data_meta), " x ", ncol(data_meta))
log_message("Dimensão curve_raw: ", nrow(curve_raw), " x ", ncol(curve_raw))

# ======================================================================================
# 10) AUDITORIA DE COLUNAS DO DEPMAP E CURVE
# ======================================================================================

log_message("Nomes das colunas em data_drug: ", paste(names(data_drug), collapse = ", "))
log_message("Nomes das colunas em curve_raw: ", paste(names(curve_raw), collapse = ", "))

colunas_meta_necessarias <- c(
  "depmap_id", "cell_line", "cell_line_name", "lineage",
  "lineage_subtype", "lineage_sub_subtype",
  "primary_disease", "subtype_disease", "primary_or_metastasis"
)

colunas_meta_faltantes <- setdiff(colunas_meta_necessarias, names(data_meta))

if (length(colunas_meta_faltantes) > 0) {
  stop("Faltam colunas em data_meta: ", paste(colunas_meta_faltantes, collapse = ", "))
}

# ======================================================================================
# 11) PADRONIZAÇÃO E JUNÇÃO DOS DADOS
# ======================================================================================

data_meta2 <- data_meta %>%
  mutate(
    depmap_id = as.character(depmap_id),
    cell_line = as.character(cell_line),
    cell_line_name = as.character(cell_line_name),
    lineage = as.character(lineage),
    lineage_subtype = as.character(lineage_subtype),
    lineage_sub_subtype = as.character(lineage_sub_subtype),
    primary_disease = as.character(primary_disease),
    subtype_disease = as.character(subtype_disease),
    primary_or_metastasis = as.character(primary_or_metastasis)
  )

data_drug2 <- data_drug %>%
  mutate(
    depmap_id = as.character(depmap_id),
    cell_line = as.character(cell_line),
    compound = as.character(if ("compound" %in% names(.)) compound else NA_character_),
    dependency = suppressWarnings(as.numeric(if ("dependency" %in% names(.)) dependency else NA_real_)),
    broad_id = as.character(if ("broad_id" %in% names(.)) broad_id else NA_character_),
    name = as.character(if ("name" %in% names(.)) name else NA_character_),
    dose = suppressWarnings(as.numeric(if ("dose" %in% names(.)) dose else NA_real_)),
    screen_id = as.character(if ("screen_id" %in% names(.)) screen_id else NA_character_),
    moa = as.character(if ("moa" %in% names(.)) moa else NA_character_),
    target = as.character(if ("target" %in% names(.)) target else NA_character_),
    disease_area = as.character(if ("disease_area" %in% names(.)) disease_area else NA_character_),
    indication = as.character(if ("indication" %in% names(.)) indication else NA_character_),
    smiles = as.character(if ("smiles" %in% names(.)) smiles else NA_character_),
    phase = as.character(if ("phase" %in% names(.)) phase else NA_character_)
  )

drug_meta <- data_drug2 %>%
  left_join(
    data_meta2 %>%
      select(
        depmap_id,
        cell_line,
        cell_line_name,
        lineage,
        lineage_subtype,
        lineage_sub_subtype,
        primary_disease,
        subtype_disease,
        primary_or_metastasis
      ),
    by = c("depmap_id", "cell_line")
  )

curve_tbl <- curve_raw %>%
  rename_with(~ gsub("^disease\\.area$", "disease_area", .x)) %>%
  mutate(
    broad_id = as.character(if ("broad_id" %in% names(.)) broad_id else NA_character_),
    depmap_id = as.character(if ("depmap_id" %in% names(.)) depmap_id else NA_character_),
    ccle_name = as.character(if ("ccle_name" %in% names(.)) ccle_name else NA_character_),
    screen_id = as.character(if ("screen_id" %in% names(.)) screen_id else NA_character_),
    upper_limit = suppressWarnings(as.numeric(if ("upper_limit" %in% names(.)) upper_limit else NA_real_)),
    lower_limit = suppressWarnings(as.numeric(if ("lower_limit" %in% names(.)) lower_limit else NA_real_)),
    slope = suppressWarnings(as.numeric(if ("slope" %in% names(.)) slope else NA_real_)),
    r2 = suppressWarnings(as.numeric(if ("r2" %in% names(.)) r2 else NA_real_)),
    auc = suppressWarnings(as.numeric(if ("auc" %in% names(.)) auc else NA_real_)),
    ec50 = suppressWarnings(as.numeric(if ("ec50" %in% names(.)) ec50 else NA_real_)),
    ic50 = suppressWarnings(as.numeric(if ("ic50" %in% names(.)) ic50 else NA_real_)),
    name = as.character(if ("name" %in% names(.)) name else NA_character_),
    moa = as.character(if ("moa" %in% names(.)) moa else NA_character_),
    target = as.character(if ("target" %in% names(.)) target else NA_character_),
    disease_area = as.character(if ("disease_area" %in% names(.)) disease_area else NA_character_),
    indication = as.character(if ("indication" %in% names(.)) indication else NA_character_),
    smiles = as.character(if ("smiles" %in% names(.)) smiles else NA_character_),
    phase = as.character(if ("phase" %in% names(.)) phase else NA_character_),
    passed_str_profiling = as.character(if ("passed_str_profiling" %in% names(.)) passed_str_profiling else NA_character_),
    row_name = as.character(if ("row_name" %in% names(.)) row_name else NA_character_),
    log10_ic50 = safe_log10(ic50),
    log10_ec50 = safe_log10(ec50),
    interpretacao_ic50 = interpretar_ic50(ic50),
    interpretacao_r2 = interpretar_r2(r2)
  )

curve_tbl <- curve_tbl %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  tibble::as_tibble()

curve_tbl <- curve_tbl %>%
  mutate(
    depmap_id = as.character(depmap_id),
    broad_id = as.character(broad_id),
    name = as.character(name),
    screen_id = as.character(screen_id),
    ccle_name = as.character(ccle_name),
    moa = as.character(moa),
    target = as.character(target),
    disease_area = as.character(disease_area),
    indication = as.character(indication),
    smiles = as.character(smiles),
    phase = as.character(phase),
    passed_str_profiling = as.character(passed_str_profiling),
    row_name = as.character(row_name),
    upper_limit = suppressWarnings(as.numeric(upper_limit)),
    lower_limit = suppressWarnings(as.numeric(lower_limit)),
    slope = suppressWarnings(as.numeric(slope)),
    r2 = suppressWarnings(as.numeric(r2)),
    auc = suppressWarnings(as.numeric(auc)),
    ec50 = suppressWarnings(as.numeric(ec50)),
    ic50 = suppressWarnings(as.numeric(ic50))
  )

curve_tbl_best <- curve_tbl %>%
  mutate(
    prioridade_screen = case_when(
      screen_id == "MTS010" ~ 1L,
      TRUE ~ 2L
    )
  ) %>%
  arrange(depmap_id, broad_id, name, prioridade_screen, desc(r2), ic50) %>%
  group_by(depmap_id, broad_id, name) %>%
  slice_head(n = 1) %>%
  ungroup()

curve_meta <- curve_tbl_best %>%
  left_join(
    data_meta2 %>%
      select(
        depmap_id,
        cell_line,
        cell_line_name,
        lineage,
        lineage_subtype,
        lineage_sub_subtype,
        primary_disease,
        subtype_disease,
        primary_or_metastasis
      ),
    by = "depmap_id"
  )

drug_curve_meta <- drug_meta %>%
  left_join(
    curve_meta %>%
      select(
        depmap_id,
        broad_id,
        name,
        ccle_name,
        screen_id,
        upper_limit,
        lower_limit,
        slope,
        r2,
        auc,
        ec50,
        ic50,
        passed_str_profiling,
        row_name,
        log10_ic50,
        log10_ec50,
        interpretacao_ic50,
        interpretacao_r2
      ),
    by = c("depmap_id", "broad_id", "name"),
    suffix = c("", "_curve")
  ) %>%
  mutate(
    screen_id = coalesce(screen_id_curve, screen_id)
  ) %>%
  select(-screen_id_curve)

log_message("Junção drug + metadata concluída.")
log_message("Dimensão drug_meta: ", nrow(drug_meta), " x ", ncol(drug_meta))
log_message("Junção curve + metadata concluída.")
log_message("Dimensão curve_meta: ", nrow(curve_meta), " x ", ncol(curve_meta))
log_message("Objeto final enriquecido drug_curve_meta criado.")
log_message("Dimensão drug_curve_meta: ", nrow(drug_curve_meta), " x ", ncol(drug_curve_meta))

# ======================================================================================
# 12) MAPA TCGA -> DEPMAP
# ======================================================================================

tipos_unicos <- sort(unique(na.omit(assinaturas_tbl$Tipo)))
mapa_tipos <- criar_mapa_tcga_depmap(tipos_unicos)

write_csv(
  mapa_tipos,
  file.path(output_dir, "mapa_tcga_para_depmap.csv")
)

log_message("Mapa TCGA -> DepMap criado para ", nrow(mapa_tipos), " tipos.")

# ======================================================================================
# 13) AUDITORIA DO MAPA E DOS CONTEXTOS CAPTURADOS
# ======================================================================================

auditoria_contexto_mapa <- map_dfr(
  tipos_unicos,
  ~ auditar_contexto_tcga_depmap(
    df_base = drug_curve_meta,
    tipo_tcga = .x,
    mapa_tcga_depmap = mapa_tipos
  )
)

write_csv(
  auditoria_contexto_mapa,
  file.path(output_dir, "auditoria_contexto_tcga_depmap.csv")
)

auditoria_contexto_detalhada <- map_dfr(
  tipos_unicos,
  function(tipo_i) {
    filtrado_i <- filtrar_depmap_por_tipo_tcga(drug_curve_meta, tipo_i, mapa_tipos)
    
    if (nrow(filtrado_i) == 0) {
      return(tibble(
        Tipo = tipo_i,
        lineage = NA_character_,
        lineage_subtype = NA_character_,
        lineage_sub_subtype = NA_character_,
        n_depmap_ids = 0L,
        n_cell_lines = 0L
      ))
    }
    
    filtrado_i %>%
      distinct(depmap_id, cell_line, lineage, lineage_subtype, lineage_sub_subtype) %>%
      group_by(lineage, lineage_subtype, lineage_sub_subtype) %>%
      summarise(
        n_depmap_ids = n_distinct(depmap_id),
        n_cell_lines = n_distinct(cell_line),
        .groups = "drop"
      ) %>%
      mutate(Tipo = tipo_i) %>%
      relocate(Tipo)
  }
)

write_csv(
  auditoria_contexto_detalhada,
  file.path(output_dir, "auditoria_contexto_tcga_depmap_detalhada.csv")
)

log_message("Auditoria do contexto tumoral exportada.")

# ======================================================================================
# 14) FILTRAGEM POR TIPO TUMORAL
# ======================================================================================

log_message("Aplicando filtro tumoral hierárquico por Tipo...")

drug_por_tipo <- map_dfr(
  tipos_unicos,
  function(tipo_i) {
    filtrado_i <- filtrar_depmap_por_tipo_tcga(drug_curve_meta, tipo_i, mapa_tipos)
    
    if (nrow(filtrado_i) == 0) {
      return(tibble())
    }
    
    filtrado_i %>%
      mutate(Tipo = tipo_i)
  }
)

log_message("Registros após filtro hierárquico por Tipo: ", nrow(drug_por_tipo))

auditoria_pos_filtro_por_tipo <- drug_por_tipo %>%
  count(Tipo, name = "n_depois_filtro") %>%
  right_join(
    mapa_tipos %>%
      select(
        Tipo,
        nivel_filtro,
        lineage_alvo,
        lineage_subtype_alvo,
        lineage_sub_subtype_alvo,
        observacao
      ),
    by = "Tipo"
  ) %>%
  mutate(
    n_depois_filtro = coalesce(n_depois_filtro, 0L)
  ) %>%
  arrange(Tipo)

write_csv(
  auditoria_pos_filtro_por_tipo,
  file.path(output_dir, "auditoria_pos_filtro_analitico.csv")
)

# ======================================================================================
# 15) TABELA DE CONTEXTO ASSINATURA -> LINHAGEM
# ======================================================================================

linhagens_por_tipo <- drug_por_tipo %>%
  distinct(
    Tipo,
    depmap_id,
    cell_line,
    cell_line_name,
    lineage,
    lineage_subtype,
    lineage_sub_subtype,
    primary_disease,
    subtype_disease,
    primary_or_metastasis
  )

assinatura_contexto_linhagem <- assinaturas_tbl %>%
  left_join(
    linhagens_por_tipo,
    by = "Tipo",
    relationship = "many-to-many"
  ) %>%
  left_join(
    mapa_tipos %>%
      select(
        Tipo,
        nivel_filtro,
        lineage_alvo,
        lineage_subtype_alvo,
        lineage_sub_subtype_alvo,
        observacao
      ),
    by = "Tipo"
  ) %>%
  arrange(Tipo, Signature, cell_line)

write_csv(
  assinatura_contexto_linhagem,
  file.path(output_dir, "assinatura_contexto_linhagem.csv")
)

mapa_partes_contexto <- map_dfr(
  tipos_unicos,
  function(tipo_i) {
    tipo_dir_i <- file.path(type_dir, tipo_i)
    dir.create(tipo_dir_i, recursive = TRUE, showWarnings = FALSE)
    
    df_dep_i <- assinatura_contexto_linhagem %>% filter(Tipo == tipo_i)
    
    exportar_em_partes(
      df = df_dep_i,
      output_prefix = file.path(tipo_dir_i, paste0("assinatura_contexto_linhagem_", tipo_i)),
      max_linhas_por_arquivo = max_linhas_por_arquivo
    ) %>%
      mutate(
        Tipo = tipo_i,
        fonte = "contexto_linhagem"
      )
  }
)

write_csv(
  mapa_partes_contexto,
  file.path(output_dir, "mapa_arquivos_particionados_contexto.csv")
)

log_message("Tabela assinatura -> contexto de linhagens exportada, inclusive em partes por tipo.")

# ======================================================================================
# 15B) MATCH EXPLÍCITO ENTRE COMPONENTES DA ASSINATURA E TARGET DO FÁRMACO
# ======================================================================================

log_message("Construindo match explícito entre assinatura e alvo do fármaco...")

assinaturas_componentes_alvo <- assinaturas_componentes %>%
  filter(
    analisavel_para_biomarcador,
    !is.na(componente_upper),
    componente_upper != ""
  ) %>%
  distinct(
    row_id,
    Tipo,
    Signature,
    Signature_limpa,
    componente,
    componente_upper
  )

drug_por_tipo_targets <- drug_por_tipo %>%
  mutate(
    target_tokens = map(target, quebrar_targets_farmaco)
  ) %>%
  tidyr::unnest_longer(target_tokens, values_to = "target_token", keep_empty = TRUE) %>%
  mutate(
    target_token = normalizar_token_alvo(target_token)
  )

match_assinatura_alvo_farmaco <- drug_por_tipo_targets %>%
  inner_join(
    assinaturas_componentes_alvo,
    by = c("Tipo", "target_token" = "componente_upper"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    componente_match = componente,
    target_match = target_token
  ) %>%
  select(
    row_id,
    Tipo,
    Signature,
    Signature_limpa,
    componente_match,
    target_match,
    depmap_id,
    cell_line,
    cell_line_name,
    ccle_name,
    lineage,
    lineage_subtype,
    lineage_sub_subtype,
    primary_disease,
    subtype_disease,
    primary_or_metastasis,
    compound,
    name,
    broad_id,
    dose,
    dependency,
    screen_id,
    upper_limit,
    lower_limit,
    slope,
    r2,
    auc,
    ec50,
    ic50,
    log10_ec50,
    log10_ic50,
    interpretacao_ic50,
    interpretacao_r2,
    moa,
    target,
    disease_area,
    indication,
    phase,
    passed_str_profiling,
    row_name,
    smiles
  )

match_assinatura_alvo_farmaco_detalhado <- match_assinatura_alvo_farmaco %>%
  group_by(
    row_id,
    Tipo,
    Signature,
    Signature_limpa,
    depmap_id,
    cell_line,
    cell_line_name,
    ccle_name,
    lineage,
    lineage_subtype,
    lineage_sub_subtype,
    primary_disease,
    subtype_disease,
    primary_or_metastasis,
    compound,
    name,
    broad_id,
    dose,
    dependency,
    screen_id,
    upper_limit,
    lower_limit,
    slope,
    r2,
    auc,
    ec50,
    ic50,
    log10_ec50,
    log10_ic50,
    interpretacao_ic50,
    interpretacao_r2,
    moa,
    target,
    disease_area,
    indication,
    phase,
    passed_str_profiling,
    row_name,
    smiles
  ) %>%
  summarise(
    n_genes_alvo_em_comum = n_distinct(componente_match),
    genes_alvo_em_comum = collapse_unique_plus(componente_match),
    targets_match = collapse_unique_plus(target_match),
    .groups = "drop"
  )

write_csv(
  match_assinatura_alvo_farmaco,
  file.path(output_dir, "match_assinatura_alvo_farmaco_raw.csv")
)

write_csv(
  match_assinatura_alvo_farmaco_detalhado,
  file.path(output_dir, "match_assinatura_alvo_farmaco_detalhado.csv")
)

auditoria_match_assinatura_alvo <- assinaturas_tbl %>%
  left_join(
    match_assinatura_alvo_farmaco_detalhado %>%
      group_by(row_id, Tipo, Signature, Signature_limpa) %>%
      summarise(
        n_drogas_match = n_distinct(paste(broad_id, name, dose, sep = "||")),
        n_linhagens_match = n_distinct(cell_line),
        genes_alvo_em_comum_total = {
          genes_split <- str_split(genes_alvo_em_comum, "\\s*\\+\\s*")
          collapse_unique_plus(unlist(genes_split))
        },
        .groups = "drop"
      ),
    by = c("row_id", "Tipo", "Signature", "Signature_limpa")
  ) %>%
  mutate(
    n_drogas_match = coalesce(n_drogas_match, 0L),
    n_linhagens_match = coalesce(n_linhagens_match, 0L),
    possui_match_alvo = n_drogas_match > 0
  )

write_csv(
  auditoria_match_assinatura_alvo,
  file.path(output_dir, "auditoria_match_assinatura_alvo.csv")
)

log_message("Match assinatura-alvo do fármaco concluído.")
log_message("Dimensão match_assinatura_alvo_farmaco: ", nrow(match_assinatura_alvo_farmaco))
log_message("Dimensão match_assinatura_alvo_farmaco_detalhado: ", nrow(match_assinatura_alvo_farmaco_detalhado))

# ======================================================================================
# 16) RESULTADOS ESTRATIFICADOS POR TIPO TUMORAL
# ======================================================================================

log_message("Iniciando geração dos resultados estratificados por tipo tumoral...")

resumo_expansao_por_tipo <- map_dfr(
  tipos_unicos,
  function(tipo_i) {
    
    log_message("------------------------------------------------------------")
    log_message("Processando tipo tumoral: ", tipo_i)
    
    tipo_dir_i <- file.path(type_dir, tipo_i)
    dir.create(tipo_dir_i, recursive = TRUE, showWarnings = FALSE)
    
    drug_i <- drug_por_tipo %>%
      filter(Tipo == tipo_i)
    
    ass_i <- assinaturas_tbl %>%
      filter(Tipo == tipo_i) %>%
      select(row_id, Tipo, Signature, Signature_limpa)
    
    mapa_i <- mapa_tipos %>%
      filter(Tipo == tipo_i)
    
    unicidade_i <- auditoria_unicidade_assinaturas %>%
      filter(Tipo == tipo_i)
    
    write_csv(
      unicidade_i,
      file.path(tipo_dir_i, paste0("auditoria_unicidade_assinaturas_", tipo_i, ".csv"))
    )
    
    if (nrow(drug_i) == 0 && nrow(ass_i) == 0) {
      log_message("Tipo ", tipo_i, " sem drogas e sem assinaturas.")
      return(
        tibble(
          Tipo = tipo_i,
          n_linhas_drug = 0L,
          n_assinaturas = 0L,
          n_linhas_estimadas = 0,
          n_blocos_assinatura = 0L,
          status = "sem_dados"
        )
      )
    }
    
    n_drug_i <- nrow(drug_i)
    n_ass_i <- nrow(ass_i)
    n_estimado <- n_drug_i * n_ass_i
    
    log_message(
      "Tipo ", tipo_i,
      " | linhas drug = ", n_drug_i,
      " | assinaturas = ", n_ass_i,
      " | linhas estimadas expansão = ", n_estimado
    )
    
    write_csv(
      mapa_i,
      file.path(tipo_dir_i, paste0("mapa_", tipo_i, ".csv"))
    )
    
    auditoria_contexto_i <- auditoria_contexto_mapa %>%
      filter(Tipo == tipo_i)
    
    write_csv(
      auditoria_contexto_i,
      file.path(tipo_dir_i, paste0("auditoria_contexto_", tipo_i, ".csv"))
    )
    
    auditoria_pos_filtro_i <- auditoria_pos_filtro_por_tipo %>%
      filter(Tipo == tipo_i)
    
    write_csv(
      auditoria_pos_filtro_i,
      file.path(tipo_dir_i, paste0("auditoria_pos_filtro_", tipo_i, ".csv"))
    )
    
    assinatura_contexto_i <- assinatura_contexto_linhagem %>%
      filter(Tipo == tipo_i)
    
    write_csv(
      assinatura_contexto_i,
      file.path(tipo_dir_i, paste0("assinatura_contexto_linhagem_", tipo_i, "_geral.csv"))
    )
    
    componentes_i <- assinaturas_componentes %>%
      filter(Tipo == tipo_i)
    
    write_csv(
      componentes_i,
      file.path(tipo_dir_i, paste0("componentes_assinaturas_", tipo_i, ".csv"))
    )
    
    auditoria_match_i <- auditoria_match_assinatura_alvo %>%
      filter(Tipo == tipo_i)
    
    write_csv(
      auditoria_match_i,
      file.path(tipo_dir_i, paste0("auditoria_match_assinatura_alvo_", tipo_i, ".csv"))
    )
    
    n_blocos_estimados <- if (n_drug_i == 0) 0L else max(1L, ceiling((n_drug_i * n_ass_i) / max_linhas_por_arquivo))
    n_blocos_assinatura <- if (n_ass_i == 0) 0L else max(1L, ceiling(n_ass_i / max_assinaturas_por_bloco), n_blocos_estimados)
    
    if (nrow(drug_i) > 0 && nrow(ass_i) > 0) {
      
      detalhado_i <- match_assinatura_alvo_farmaco_detalhado %>%
        filter(Tipo == tipo_i)
      
      mapa_partes_dependency <- exportar_em_partes(
        df = detalhado_i %>%
          left_join(
            mapa_i %>%
              select(
                Tipo,
                nivel_filtro,
                lineage_alvo,
                lineage_subtype_alvo,
                lineage_sub_subtype_alvo,
                observacao
              ),
            by = "Tipo"
          ) %>%
          select(
            row_id,
            Tipo,
            Signature,
            Signature_limpa,
            genes_alvo_em_comum,
            n_genes_alvo_em_comum,
            targets_match,
            nivel_filtro,
            lineage_alvo,
            lineage_subtype_alvo,
            lineage_sub_subtype_alvo,
            observacao,
            depmap_id,
            cell_line,
            cell_line_name,
            ccle_name,
            lineage,
            lineage_subtype,
            lineage_sub_subtype,
            primary_disease,
            subtype_disease,
            primary_or_metastasis,
            compound,
            name,
            broad_id,
            dose,
            dependency,
            screen_id,
            upper_limit,
            lower_limit,
            slope,
            r2,
            auc,
            ec50,
            ic50,
            log10_ec50,
            log10_ic50,
            interpretacao_ic50,
            interpretacao_r2,
            moa,
            target,
            disease_area,
            indication,
            phase,
            passed_str_profiling,
            row_name,
            smiles
          ) %>%
          arrange(Signature, name, dependency),
        output_prefix = file.path(
          tipo_dir_i,
          paste0("depmap_drug_sensitivity_tabela_linhagens_com_assinatura_", tipo_i, "_match_alvo")
        ),
        max_linhas_por_arquivo = max_linhas_por_arquivo
      ) %>%
        mutate(
          Tipo = tipo_i,
          fonte = "dependency_match_alvo"
        )
      
      write_csv(
        mapa_partes_dependency,
        file.path(tipo_dir_i, paste0("mapa_partes_dependency_", tipo_i, ".csv"))
      )
      
      resumo_droga_tipo_i <- drug_i %>%
        group_by(Tipo, compound, name, broad_id, moa, target, dose, phase) %>%
        summarise(
          n_registros = n(),
          n_linhagens = n_distinct(cell_line),
          media_dependency_farmaco = safe_mean(dependency),
          mediana_dependency_farmaco = safe_median(dependency),
          dp_dependency_farmaco = safe_sd(dependency),
          min_dependency_farmaco = safe_min(dependency),
          max_dependency_farmaco = safe_max(dependency),
          media_auc = safe_mean(auc),
          mediana_auc = safe_median(auc),
          media_ic50 = safe_mean(ic50),
          mediana_ic50 = safe_median(ic50),
          media_ec50 = safe_mean(ec50),
          mediana_ec50 = safe_median(ec50),
          media_r2 = safe_mean(r2),
          melhor_ic50 = safe_min(ic50),
          prop_str_ok = safe_prop_true(passed_str_profiling),
          interpretacao_global = interpretar_dependency_farmaco(media_dependency_farmaco),
          interpretacao_ic50_global = interpretar_ic50(mediana_ic50),
          interpretacao_r2_global = interpretar_r2(media_r2),
          .groups = "drop"
        ) %>%
        arrange(media_dependency_farmaco)
      
      write_csv(
        resumo_droga_tipo_i,
        file.path(tipo_dir_i, paste0("depmap_drug_sensitivity_resumo_por_droga_tipo_", tipo_i, ".csv"))
      )
      
      resumo_droga_assinatura_i <- match_assinatura_alvo_farmaco_detalhado %>%
        filter(Tipo == tipo_i) %>%
        group_by(
          row_id,
          Tipo,
          Signature,
          Signature_limpa,
          compound,
          name,
          broad_id,
          moa,
          target,
          dose,
          phase
        ) %>%
        summarise(
          n_registros = n(),
          n_linhagens = n_distinct(cell_line),
          media_dependency_farmaco = safe_mean(dependency),
          mediana_dependency_farmaco = safe_median(dependency),
          dp_dependency_farmaco = safe_sd(dependency),
          min_dependency_farmaco = safe_min(dependency),
          max_dependency_farmaco = safe_max(dependency),
          media_auc = safe_mean(auc),
          mediana_auc = safe_median(auc),
          media_ic50 = safe_mean(ic50),
          mediana_ic50 = safe_median(ic50),
          media_ec50 = safe_mean(ec50),
          mediana_ec50 = safe_median(ec50),
          media_r2 = safe_mean(r2),
          melhor_ic50 = safe_min(ic50),
          prop_str_ok = safe_prop_true(passed_str_profiling),
          n_genes_alvo_em_comum = safe_max_int(n_genes_alvo_em_comum),
          genes_alvo_em_comum = collapse_unique_bar(genes_alvo_em_comum),
          targets_match = collapse_unique_bar(targets_match),
          interpretacao_global = interpretar_dependency_farmaco(media_dependency_farmaco),
          interpretacao_ic50_global = interpretar_ic50(mediana_ic50),
          interpretacao_r2_global = interpretar_r2(media_r2),
          .groups = "drop"
        ) %>%
        arrange(Signature, media_dependency_farmaco)
      
      mapa_resumo_dep <- exportar_em_partes(
        df = resumo_droga_assinatura_i,
        output_prefix = file.path(tipo_dir_i, paste0("depmap_drug_sensitivity_resumo_por_droga_assinatura_", tipo_i)),
        max_linhas_por_arquivo = max_linhas_por_arquivo
      )
      
      write_csv(
        mapa_resumo_dep,
        file.path(tipo_dir_i, paste0("mapa_partes_resumo_dependency_", tipo_i, ".csv"))
      )
      
      top_mais_i <- resumo_droga_tipo_i %>%
        filter(!is.na(media_dependency_farmaco)) %>%
        slice_min(order_by = media_dependency_farmaco, n = 20, with_ties = FALSE)
      
      write_csv(
        top_mais_i,
        file.path(tipo_dir_i, paste0("depmap_top_drogas_mais_sensiveis_", tipo_i, ".csv"))
      )
      
      top_menos_i <- resumo_droga_tipo_i %>%
        filter(!is.na(media_dependency_farmaco)) %>%
        slice_max(order_by = media_dependency_farmaco, n = 20, with_ties = FALSE)
      
      write_csv(
        top_menos_i,
        file.path(tipo_dir_i, paste0("depmap_top_drogas_menos_sensiveis_", tipo_i, ".csv"))
      )
      
      top_mais_ass_i <- resumo_droga_assinatura_i %>%
        filter(!is.na(media_dependency_farmaco)) %>%
        group_by(row_id, Tipo, Signature, Signature_limpa) %>%
        slice_min(order_by = media_dependency_farmaco, n = 20, with_ties = FALSE) %>%
        ungroup() %>%
        arrange(Signature, media_dependency_farmaco)
      
      mapa_top_mais_dep <- exportar_em_partes(
        df = top_mais_ass_i,
        output_prefix = file.path(tipo_dir_i, paste0("depmap_top_drogas_mais_sensiveis_por_assinatura_", tipo_i)),
        max_linhas_por_arquivo = max_linhas_por_arquivo
      )
      
      write_csv(
        mapa_top_mais_dep,
        file.path(tipo_dir_i, paste0("mapa_partes_top_mais_dependency_", tipo_i, ".csv"))
      )
      
      top_menos_ass_i <- resumo_droga_assinatura_i %>%
        filter(!is.na(media_dependency_farmaco)) %>%
        group_by(row_id, Tipo, Signature, Signature_limpa) %>%
        slice_max(order_by = media_dependency_farmaco, n = 20, with_ties = FALSE) %>%
        ungroup() %>%
        arrange(Signature, desc(media_dependency_farmaco))
      
      mapa_top_menos_dep <- exportar_em_partes(
        df = top_menos_ass_i,
        output_prefix = file.path(tipo_dir_i, paste0("depmap_top_drogas_menos_sensiveis_por_assinatura_", tipo_i)),
        max_linhas_por_arquivo = max_linhas_por_arquivo
      )
      
      write_csv(
        mapa_top_menos_dep,
        file.path(tipo_dir_i, paste0("mapa_partes_top_menos_dependency_", tipo_i, ".csv"))
      )
      
      rm(
        resumo_droga_tipo_i,
        resumo_droga_assinatura_i,
        top_mais_i,
        top_menos_i,
        top_mais_ass_i,
        top_menos_ass_i,
        mapa_partes_dependency,
        mapa_resumo_dep,
        mapa_top_mais_dep,
        mapa_top_menos_dep
      )
      gc()
      
      status_i <- "ok"
      
    } else {
      status_i <- "sem_expansao"
      log_message("Tipo ", tipo_i, " sem dados suficientes para tabela grande.")
    }
    
    tibble(
      Tipo = tipo_i,
      n_linhas_drug = n_drug_i,
      n_assinaturas = n_ass_i,
      n_linhas_estimadas = n_estimado,
      n_blocos_assinatura = n_blocos_assinatura,
      status = status_i
    )
  }
)

write_csv(
  resumo_expansao_por_tipo,
  file.path(output_dir, "resumo_expansao_tabela_grande_por_tipo.csv")
)

log_message("Resultados estratificados por tipo tumoral exportados com particionamento.")

# ======================================================================================
# 17) RESULTADOS GERAIS LEVES
# ======================================================================================

resumo_droga_tipo <- drug_por_tipo %>%
  group_by(Tipo, compound, name, broad_id, moa, target, dose, phase) %>%
  summarise(
    n_registros = n(),
    n_linhagens = n_distinct(cell_line),
    media_dependency_farmaco = safe_mean(dependency),
    mediana_dependency_farmaco = safe_median(dependency),
    dp_dependency_farmaco = safe_sd(dependency),
    min_dependency_farmaco = safe_min(dependency),
    max_dependency_farmaco = safe_max(dependency),
    media_auc = safe_mean(auc),
    mediana_auc = safe_median(auc),
    media_ic50 = safe_mean(ic50),
    mediana_ic50 = safe_median(ic50),
    media_ec50 = safe_mean(ec50),
    mediana_ec50 = safe_median(ec50),
    media_r2 = safe_mean(r2),
    melhor_ic50 = safe_min(ic50),
    prop_str_ok = safe_prop_true(passed_str_profiling),
    interpretacao_global = interpretar_dependency_farmaco(media_dependency_farmaco),
    interpretacao_ic50_global = interpretar_ic50(mediana_ic50),
    interpretacao_r2_global = interpretar_r2(media_r2),
    .groups = "drop"
  ) %>%
  arrange(Tipo, media_dependency_farmaco)

write_csv(
  resumo_droga_tipo,
  file.path(output_dir, "depmap_drug_sensitivity_resumo_por_droga_tipo.csv")
)

resumo_droga_assinatura <- match_assinatura_alvo_farmaco_detalhado %>%
  group_by(
    row_id,
    Tipo,
    Signature,
    Signature_limpa,
    compound,
    name,
    broad_id,
    moa,
    target,
    dose,
    phase
  ) %>%
  summarise(
    n_registros = n(),
    n_linhagens = n_distinct(cell_line),
    media_dependency_farmaco = safe_mean(dependency),
    mediana_dependency_farmaco = safe_median(dependency),
    dp_dependency_farmaco = safe_sd(dependency),
    min_dependency_farmaco = safe_min(dependency),
    max_dependency_farmaco = safe_max(dependency),
    media_auc = safe_mean(auc),
    mediana_auc = safe_median(auc),
    media_ic50 = safe_mean(ic50),
    mediana_ic50 = safe_median(ic50),
    media_ec50 = safe_mean(ec50),
    mediana_ec50 = safe_median(ec50),
    media_r2 = safe_mean(r2),
    melhor_ic50 = safe_min(ic50),
    prop_str_ok = safe_prop_true(passed_str_profiling),
    n_genes_alvo_em_comum = safe_max_int(n_genes_alvo_em_comum),
    genes_alvo_em_comum = collapse_unique_bar(genes_alvo_em_comum),
    targets_match = collapse_unique_bar(targets_match),
    interpretacao_global = interpretar_dependency_farmaco(media_dependency_farmaco),
    interpretacao_ic50_global = interpretar_ic50(mediana_ic50),
    interpretacao_r2_global = interpretar_r2(media_r2),
    .groups = "drop"
  ) %>%
  arrange(Tipo, Signature, media_dependency_farmaco)

write_csv(
  resumo_droga_assinatura,
  file.path(output_dir, "depmap_drug_sensitivity_resumo_por_droga_assinatura.csv")
)

top_drogas_mais_sensiveis <- resumo_droga_tipo %>%
  filter(!is.na(media_dependency_farmaco)) %>%
  group_by(Tipo) %>%
  slice_min(order_by = media_dependency_farmaco, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Tipo, media_dependency_farmaco)

write_csv(
  top_drogas_mais_sensiveis,
  file.path(output_dir, "depmap_top_drogas_mais_sensiveis.csv")
)

top_drogas_menos_sensiveis <- resumo_droga_tipo %>%
  filter(!is.na(media_dependency_farmaco)) %>%
  group_by(Tipo) %>%
  slice_max(order_by = media_dependency_farmaco, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Tipo, desc(media_dependency_farmaco))

write_csv(
  top_drogas_menos_sensiveis,
  file.path(output_dir, "depmap_top_drogas_menos_sensiveis.csv")
)

top_drogas_mais_sensiveis_assinatura <- resumo_droga_assinatura %>%
  filter(!is.na(media_dependency_farmaco)) %>%
  group_by(row_id, Tipo, Signature, Signature_limpa) %>%
  slice_min(order_by = media_dependency_farmaco, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Tipo, Signature, media_dependency_farmaco)

write_csv(
  top_drogas_mais_sensiveis_assinatura,
  file.path(output_dir, "depmap_top_drogas_mais_sensiveis_por_assinatura.csv")
)

top_drogas_menos_sensiveis_assinatura <- resumo_droga_assinatura %>%
  filter(!is.na(media_dependency_farmaco)) %>%
  group_by(row_id, Tipo, Signature, Signature_limpa) %>%
  slice_max(order_by = media_dependency_farmaco, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Tipo, Signature, desc(media_dependency_farmaco))

write_csv(
  top_drogas_menos_sensiveis_assinatura,
  file.path(output_dir, "depmap_top_drogas_menos_sensiveis_por_assinatura.csv")
)

log_message("Resultados gerais leves exportados.")

# ======================================================================================
# 18) AUDITORIA GERAL POR ASSINATURA
# ======================================================================================

resumo_componentes_assinatura <- assinaturas_componentes %>%
  group_by(row_id, Tipo, Signature) %>%
  summarise(
    n_componentes_total = sum(!is.na(componente) & componente != ""),
    n_componentes_biomarcador = sum(analisavel_para_biomarcador, na.rm = TRUE),
    n_componentes_nao_analisaveis = sum((!analisavel_para_biomarcador) & !is.na(componente) & componente != "", na.rm = TRUE),
    componentes_concat = paste(unique(na.omit(componente)), collapse = " + "),
    .groups = "drop"
  )

resumo_linhagens_por_assinatura <- assinatura_contexto_linhagem %>%
  group_by(row_id, Tipo, Signature, Signature_limpa) %>%
  summarise(
    n_depmap_ids_contexto = n_distinct(depmap_id[!is.na(depmap_id)]),
    n_cell_lines_contexto = n_distinct(cell_line[!is.na(cell_line)]),
    .groups = "drop"
  )

auditoria_assinaturas <- assinaturas_tbl %>%
  left_join(
    resumo_componentes_assinatura,
    by = c("row_id", "Tipo", "Signature")
  ) %>%
  left_join(
    resumo_linhagens_por_assinatura,
    by = c("row_id", "Tipo", "Signature", "Signature_limpa")
  ) %>%
  left_join(
    auditoria_match_assinatura_alvo %>%
      select(
        row_id,
        Tipo,
        Signature,
        Signature_limpa,
        n_drogas_match,
        n_linhagens_match,
        genes_alvo_em_comum_total,
        possui_match_alvo
      ),
    by = c("row_id", "Tipo", "Signature", "Signature_limpa")
  ) %>%
  left_join(
    mapa_tipos %>%
      select(
        Tipo,
        nivel_filtro,
        lineage_alvo,
        lineage_subtype_alvo,
        lineage_sub_subtype_alvo,
        observacao
      ),
    by = "Tipo"
  ) %>%
  mutate(
    n_componentes_total = coalesce(n_componentes_total, 0L),
    n_componentes_biomarcador = coalesce(n_componentes_biomarcador, 0L),
    n_componentes_nao_analisaveis = coalesce(n_componentes_nao_analisaveis, 0L),
    n_depmap_ids_contexto = coalesce(n_depmap_ids_contexto, 0L),
    n_cell_lines_contexto = coalesce(n_cell_lines_contexto, 0L),
    n_drogas_match = coalesce(n_drogas_match, 0L),
    n_linhagens_match = coalesce(n_linhagens_match, 0L),
    possui_match_alvo = coalesce(possui_match_alvo, FALSE),
    interpretacao_metodologica = case_when(
      possui_match_alvo ~ "Assinatura associada apenas a drogas cujo target coincide com pelo menos um gene da assinatura.",
      n_cell_lines_contexto > 0 ~ "Assinatura possui contexto tumoral no DepMap, mas sem droga com target coincidente na assinatura.",
      TRUE ~ "Assinatura sem linhagens compatíveis no filtro tumoral atual."
    )
  )

write_csv(
  auditoria_assinaturas,
  file.path(output_dir, "auditoria_geral_assinaturas.csv")
)

log_message("Auditoria geral por assinatura exportada.")

# ======================================================================================
# 18B) ASSINATURAS COM TODOS OS COMPONENTES SENSÍVEIS (ALVO DA DROGA)
# ======================================================================================

log_message("Filtrando assinaturas onde todos os componentes são alvo da mesma droga...")

assinaturas_todos_componentes_alvo <- resumo_droga_assinatura %>%
  left_join(
    resumo_componentes_assinatura %>%
      select(row_id, Tipo, Signature, n_componentes_biomarcador),
    by = c("row_id", "Tipo", "Signature")
  ) %>%
  filter(n_componentes_biomarcador > 0) %>%
  filter(n_genes_alvo_em_comum == n_componentes_biomarcador) %>%
  arrange(Tipo, Signature, media_dependency_farmaco)

write_csv(
  assinaturas_todos_componentes_alvo,
  file.path(output_dir, "assinaturas_todos_componentes_como_alvo.csv")
)

log_message("Análise de assinaturas com todos os componentes como alvo exportada. Encontradas: ", nrow(assinaturas_todos_componentes_alvo))

# ======================================================================================
# 19) SAÍDA FINAL NO CONSOLE
# ======================================================================================

cat("\n================ RESULTADO FINAL ================\n")
cat("Diretório geral:", normalizePath(output_dir), "\n")
cat("Diretório por tipo:", normalizePath(type_dir), "\n\n")

cat("Resumo por droga e tipo:\n")
print(resumo_droga_tipo, n = min(nrow(resumo_droga_tipo), 30))

cat("\nResumo por droga e assinatura com match de alvo:\n")
print(resumo_droga_assinatura, n = min(nrow(resumo_droga_assinatura), 30))

cat("\nTop drogas mais sensíveis por tipo:\n")
print(top_drogas_mais_sensiveis, n = min(nrow(top_drogas_mais_sensiveis), 30))

cat("\nTop drogas menos sensíveis por tipo:\n")
print(top_drogas_menos_sensiveis, n = min(nrow(top_drogas_menos_sensiveis), 30))

cat("\nAuditoria geral das assinaturas:\n")
print(auditoria_assinaturas, n = min(nrow(auditoria_assinaturas), 30))

cat("\nAuditoria de unicidade das assinaturas por tipo:\n")
print(auditoria_unicidade_assinaturas, n = nrow(auditoria_unicidade_assinaturas))

cat("\nResumo de expansão por tipo:\n")
print(resumo_expansao_por_tipo, n = nrow(resumo_expansao_por_tipo))

cat("\nAssinaturas onde TODOS os componentes são alvo da droga:\n")
print(assinaturas_todos_componentes_alvo, n = min(nrow(assinaturas_todos_componentes_alvo), 30))

cat("\n=================================================\n")

log_message("FIM DA ANÁLISE")