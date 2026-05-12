# ======================================================================================
# ANÁLISE DE ESSENCIALIDADE DE COMPONENTES DE ASSINATURAS ÔMICAS NO DEPMAP (CRISPR)
# Objetivo:
# Ler assinaturas da coluna Signature, decompor seus componentes (genes/transcritos),
# identificar quais genes são essenciais para manutenção da viabilidade celular
# em contexto tumoral, usando depmap_crispr() + depmap_metadata().
#
# Entrada:
# G:\DepMap\dependencia_celular\Assinaturas_Omicas_Uma_Por_Linha_Anotadas.csv
#
# Observação metodológica:
# O DepMap CRISPR é gene-level. Transcritos/isoformas não são avaliados diretamente
# como dependência funcional no objeto depmap_crispr(). Assim, componentes com cara
# de transcrito serão auditados separadamente.
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

arquivo_assinaturas <- "D:/DepMap/dependencia_celular/Assinaturas_Omicas_SuperLearner.tsv"
output_dir <- "D:/DepMap/dependencia_celular/depmap_assinaturas_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(output_dir, "analysis_log.txt")

log_message <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  write(msg, file = log_file, append = TRUE)
}

if (file.exists(log_file)) file.remove(log_file)
file.create(log_file)

log_message("============================================================")
log_message("INÍCIO DA ANÁLISE DE ESSENCIALIDADE DE ASSINATURAS NO DEPMAP")
log_message("Arquivo de entrada: ", arquivo_assinaturas)
log_message("Diretório de saída: ", output_dir)
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
    is.na(x2) | x2 == ""                                                         ~ "vazio",
    str_detect(x2, "^ENST[0-9]+")                                                ~ "transcrito_ensembl",
    str_detect(x2, "^ENSG[0-9]+")                                                ~ "gene_ensembl",
    str_detect(x2, "^NM_[0-9]+")                                                 ~ "transcrito_refseq_mrna",
    str_detect(x2, "^NR_[0-9]+")                                                 ~ "transcrito_refseq_ncrna",
    str_detect(x2, "^XM_[0-9]+")                                                 ~ "transcrito_modelo_refseq",
    str_detect(x2, "^XR_[0-9]+")                                                 ~ "transcrito_modelo_refseq",
    str_detect(x2, "ISOFORM|TRANSCRIPT")                                         ~ "transcrito_textual",
    str_detect(x2, "-20[0-9]$|-2[0-9][0-9]$")                                   ~ "possivel_isoforma",
    str_detect(x2, "\\.[0-9]+$") & str_detect(x2, "^ENS[A-Z]*[0-9]+")          ~ "identificador_ensembl_versionado",
    str_detect(x2, "^[A-Z0-9._-]+$")                                             ~ "candidato_gene_symbol",
    TRUE                                                                          ~ "ambiguo"
  )
}

classificar_dependency <- function(x) {
  case_when(
    is.na(x)   ~ "Sem dado",
    x < -1.5   ~ "Forte dependência",
    x < -0.5   ~ "Dependência moderada",
    x <= 0.5   ~ "Linhagem mantida / pouco efeito",
    x > 0.5    ~ "Possível vantagem proliferativa após knockout",
    TRUE       ~ "Indeterminado"
  )
}

interpretar_dependency <- function(classe) {
  case_when(
    classe == "Forte dependência" ~
      "Gene fortemente essencial neste contexto; knockout reduz fortemente a viabilidade",
    classe == "Dependência moderada" ~
      "Gene com contribuição relevante para viabilidade; knockout reduz moderadamente a viabilidade",
    classe == "Linhagem mantida / pouco efeito" ~
      "Gene pouco essencial neste contexto; linhagem tende a ser mantida após knockout",
    classe == "Possível vantagem proliferativa após knockout" ~
      "Knockout pode favorecer crescimento ou refletir contexto biológico específico",
    TRUE ~
      "Sem interpretação disponível"
  )
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
      observacao   = coalesce(observacao, "Tipo não mapeado manualmente")
    )
}

# ======================================================================================
# 4) FUNÇÃO DE FILTRAGEM POR TIPO TCGA
# ======================================================================================

filtrar_depmap_por_tipo_tcga <- function(df, tipo_tcga, mapa_tcga_depmap) {
  
  info <- mapa_tcga_depmap %>% filter(Tipo == tipo_tcga)
  
  if (nrow(info) == 0) return(tibble())
  
  nivel  <- info$nivel_filtro[[1]]
  lin    <- info$lineage_alvo[[1]]
  sub    <- info$lineage_subtype_alvo[[1]]
  subsub <- info$lineage_sub_subtype_alvo[[1]]
  
  if (is.na(nivel) || nivel == "sem_match") return(tibble())
  
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
  
  info <- mapa_tcga_depmap %>% filter(Tipo == tipo_tcga)
  
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
    Tipo                     = tipo_tcga,
    nivel_filtro             = info$nivel_filtro[[1]],
    lineage_alvo             = info$lineage_alvo[[1]],
    lineage_subtype_alvo     = info$lineage_subtype_alvo[[1]],
    lineage_sub_subtype_alvo = info$lineage_sub_subtype_alvo[[1]],
    observacao               = info$observacao[[1]],
    n_linhas_filtradas       = nrow(filtrado),
    n_depmap_ids             = n_distinct(filtrado$depmap_id),
    n_cell_lines             = n_distinct(filtrado$cell_line)
  )
}

# ======================================================================================
# 6) LEITURA DA TABELA DE ASSINATURAS
# ======================================================================================

log_message("Lendo arquivo de assinaturas...")
assinaturas_raw <- read_tsv(arquivo_assinaturas, show_col_types = FALSE)

colunas_obrigatorias <- c("Signature", "Tipo")
faltantes_entrada    <- setdiff(colunas_obrigatorias, names(assinaturas_raw))

if (length(faltantes_entrada) > 0) {
  stop("Faltam colunas obrigatórias no arquivo de assinaturas: ",
       paste(faltantes_entrada, collapse = ", "))
}

log_message("Dimensão assinaturas_raw: ", nrow(assinaturas_raw), " x ", ncol(assinaturas_raw))

assinaturas_tbl <- assinaturas_raw %>%
  mutate(
    row_id             = row_number(),
    Signature_original = Signature,                   # valor bruto preservado do CSV
    Signature          = as.character(Signature),
    Tipo               = as.character(Tipo),
    Signature_limpa    = limpar_assinatura(Signature)
  )

# ======================================================================================
# 7) EXPANSÃO DAS ASSINATURAS EM COMPONENTES
# ======================================================================================

log_message("Expandindo assinaturas em componentes...")

assinaturas_componentes <- assinaturas_tbl %>%
  mutate(componentes = map(Signature_limpa, quebrar_componentes_assinatura)) %>%
  tidyr::unnest_longer(componentes, values_to = "componente") %>%
  mutate(
    componente              = str_trim(as.character(componente)),
    componente_upper        = toupper(componente),
    tipo_componente         = classificar_componente(componente_upper),
    analisavel_no_depmap_crispr = tipo_componente %in% c("candidato_gene_symbol", "gene_ensembl")
  )

log_message("Número de componentes expandidos: ", nrow(assinaturas_componentes))

# ======================================================================================
# 8) AUDITORIA DOS COMPONENTES
# ======================================================================================

auditoria_componentes <- assinaturas_componentes %>%
  count(tipo_componente, analisavel_no_depmap_crispr, name = "n_componentes") %>%
  arrange(desc(n_componentes))

write_csv(auditoria_componentes,
          file.path(output_dir, "auditoria_componentes_assinaturas.csv"))

componentes_nao_analisaveis <- assinaturas_componentes %>%
  filter(!analisavel_no_depmap_crispr) %>%
  select(row_id, Tipo, Signature_original, Signature, componente, tipo_componente) %>%
  distinct()

write_csv(componentes_nao_analisaveis,
          file.path(output_dir, "componentes_nao_analisaveis_no_depmap_crispr.csv"))

log_message("Auditoria de componentes exportada.")

# ======================================================================================
# 9) CARREGAMENTO DOS DADOS DEPMAP
# ======================================================================================

log_message("Carregando depmap_crispr()...")
data_crisp <- depmap::depmap_crispr()

log_message("Carregando depmap_metadata()...")
data_meta <- depmap::depmap_metadata()

# Verificação de carregamento (evitar erros se o ExperimentHub falhar)
if (is.null(data_crisp) || nrow(data_crisp) == 0 || ncol(data_crisp) == 0) {
  log_message("ERRO CRÍTICO: data_crisp está vazio ou não foi carregado.")
  log_message("Isso geralmente ocorre por falha de conexão com o ExperimentHub (Bioconductor).")
  stop("Falha ao carregar depmap_crispr(). Verifique sua conexão ou tente novamente.")
}

if (is.null(data_meta) || nrow(data_meta) == 0 || ncol(data_meta) == 0) {
  log_message("ERRO CRÍTICO: data_meta está vazio ou não foi carregado.")
  log_message("Isso geralmente ocorre por falha de conexão com o ExperimentHub (Bioconductor).")
  stop("Falha ao carregar depmap_metadata(). Verifique sua conexão ou tente novamente.")
}

log_message("Dimensão data_crisp: ", nrow(data_crisp), " x ", ncol(data_crisp))
log_message("Dimensão data_meta : ", nrow(data_meta),  " x ", ncol(data_meta))

# ======================================================================================
# 10) AUDITORIA E PADRONIZAÇÃO DE COLUNAS DO DEPMAP
# ======================================================================================

# Algumas versões do pacote ou releases do DepMap podem mudar nomes (ex: ModelID vs depmap_id)
# Esta função tenta padronizar nomes conhecidos para garantir compatibilidade
padronizar_colunas_depmap <- function(df) {
  nm <- names(df)
  # Mapeamento de nomes alternativos comuns (Portal vs Pacote R)
  if (!("depmap_id" %in% nm) && ("ModelID" %in% nm)) df <- rename(df, depmap_id = ModelID)
  if (!("cell_line" %in% nm) && ("CellLine" %in% nm)) df <- rename(df, cell_line = CellLine)
  if (!("cell_line_name" %in% nm) && ("CellLineName" %in% nm)) df <- rename(df, cell_line_name = CellLineName)
  if (!("lineage" %in% nm) && ("OncotreeLineage" %in% nm)) df <- rename(df, lineage = OncotreeLineage)
  if (!("primary_disease" %in% nm) && ("OncotreePrimaryDisease" %in% nm)) df <- rename(df, primary_disease = OncotreePrimaryDisease)
  if (!("gene_name" %in% nm) && ("symbol" %in% nm)) df <- rename(df, gene_name = symbol)
  return(df)
}

data_meta  <- padronizar_colunas_depmap(data_meta)
data_crisp <- padronizar_colunas_depmap(data_crisp)

colunas_crisp_necessarias <- c(
  "depmap_id", "gene", "dependency", "entrez_id", "gene_name", "cell_line"
)

colunas_meta_necessarias <- c(
  "depmap_id", "cell_line", "cell_line_name", "lineage",
  "lineage_subtype", "lineage_sub_subtype",
  "primary_disease", "subtype_disease", "primary_or_metastasis"
)

colunas_crisp_faltantes <- setdiff(colunas_crisp_necessarias, names(data_crisp))
colunas_meta_faltantes  <- setdiff(colunas_meta_necessarias,  names(data_meta))

if (length(colunas_crisp_faltantes) > 0) {
  log_message("Colunas encontradas em data_crisp: ", paste(names(data_crisp), collapse = ", "))
  stop("Faltam colunas em data_crisp: ", paste(colunas_crisp_faltantes, collapse = ", "))
}

if (length(colunas_meta_faltantes) > 0) {
  log_message("Colunas encontradas em data_meta: ", paste(names(data_meta), collapse = ", "))
  stop("Faltam colunas em data_meta: ", paste(colunas_meta_faltantes, collapse = ", "))
}

log_message("Auditoria de colunas do DepMap: OK")

# ======================================================================================
# 11) PADRONIZAÇÃO E JUNÇÃO DOS DADOS DEPMAP
# ======================================================================================

data_crisp2 <- data_crisp %>%
  mutate(
    gene      = toupper(as.character(gene)),
    gene_name = toupper(as.character(gene_name)),
    cell_line = as.character(cell_line)
  )

data_meta2 <- data_meta %>%
  mutate(
    cell_line             = as.character(cell_line),
    cell_line_name        = as.character(cell_line_name),
    lineage               = as.character(lineage),
    lineage_subtype       = as.character(lineage_subtype),
    lineage_sub_subtype   = as.character(lineage_sub_subtype),
    primary_disease       = as.character(primary_disease),
    subtype_disease       = as.character(subtype_disease),
    primary_or_metastasis = as.character(primary_or_metastasis)
  )

crisp_meta <- data_crisp2 %>%
  left_join(
    data_meta2 %>%
      select(
        depmap_id, cell_line, cell_line_name,
        lineage, lineage_subtype, lineage_sub_subtype,
        primary_disease, subtype_disease, primary_or_metastasis
      ),
    by = c("depmap_id", "cell_line")
  )

log_message("Junção do DepMap concluída.")
log_message("Dimensão crisp_meta: ", nrow(crisp_meta), " x ", ncol(crisp_meta))

# ======================================================================================
# 12) SELEÇÃO DOS COMPONENTES ANALISÁVEIS
# ======================================================================================

componentes_genicos <- assinaturas_componentes %>%
  filter(analisavel_no_depmap_crispr) %>%
  mutate(componente_busca = componente_upper) %>%
  distinct(
    row_id, Tipo, Signature, Signature_original, Signature_limpa,
    componente, componente_upper, componente_busca, tipo_componente
  )

if (nrow(componentes_genicos) == 0) {
  stop("Nenhum componente analisável como gene foi encontrado na coluna Signature.")
}

log_message("Componentes gênicos analisáveis: ", nrow(componentes_genicos))
log_message("Genes únicos candidatos: ", n_distinct(componentes_genicos$componente_busca))

# ======================================================================================
# 13) CRUZAMENTO ENTRE COMPONENTES GÊNICOS E DEPMAP
# ======================================================================================

res_genes <- componentes_genicos %>%
  left_join(
    crisp_meta %>%
      select(
        depmap_id, cell_line, cell_line_name,
        gene, gene_name, entrez_id, dependency,
        lineage, lineage_subtype, lineage_sub_subtype,
        primary_disease, subtype_disease, primary_or_metastasis
      ),
    by = c("componente_busca" = "gene_name")
  )

nao_bateram <- res_genes %>%
  filter(is.na(depmap_id)) %>%
  select(
    row_id, Tipo, Signature, Signature_original, Signature_limpa,
    componente, componente_upper, componente_busca, tipo_componente
  ) %>%
  distinct()

res_fallback <- nao_bateram %>%
  left_join(
    crisp_meta %>%
      select(
        depmap_id, cell_line, cell_line_name,
        gene, gene_name, entrez_id, dependency,
        lineage, lineage_subtype, lineage_sub_subtype,
        primary_disease, subtype_disease, primary_or_metastasis
      ),
    by = c("componente_busca" = "gene")
  )

res_genes_ok      <- res_genes    %>% filter(!is.na(depmap_id))
res_fallback_ok   <- res_fallback %>% filter(!is.na(depmap_id))

res_genes_total <- bind_rows(res_genes_ok, res_fallback_ok) %>% distinct()

log_message("Registros após cruzamento com DepMap: ", nrow(res_genes_total))

# ======================================================================================
# 14) AUDITORIA DE COMPONENTES GÊNICOS NÃO ENCONTRADOS NO DEPMAP
# ======================================================================================

genes_encontrados_depmap <- res_genes_total %>%
  distinct(componente_busca) %>%
  pull(componente_busca)

auditoria_genes_assinaturas <- componentes_genicos %>%
  distinct(row_id, Tipo, Signature, Signature_original, componente, componente_busca, tipo_componente) %>%
  mutate(encontrado_no_depmap = componente_busca %in% genes_encontrados_depmap)

write_csv(auditoria_genes_assinaturas,
          file.path(output_dir, "auditoria_genes_assinaturas_vs_depmap.csv"))

log_message("Auditoria genes x DepMap exportada.")

# ======================================================================================
# 15) MAPA TCGA -> DEPMAP
# ======================================================================================

tipos_unicos <- sort(unique(na.omit(assinaturas_tbl$Tipo)))
mapa_tipos   <- criar_mapa_tcga_depmap(tipos_unicos)

write_csv(mapa_tipos,
          file.path(output_dir, "mapa_tcga_para_depmap.csv"))

log_message("Mapa TCGA -> DepMap criado para ", nrow(mapa_tipos), " tipos.")

# ======================================================================================
# 16) AUDITORIA DO MAPA E DOS CONTEXTOS CAPTURADOS
# ======================================================================================

auditoria_contexto_mapa <- map_dfr(
  tipos_unicos,
  ~ auditar_contexto_tcga_depmap(
    df_base       = crisp_meta,
    tipo_tcga     = .x,
    mapa_tcga_depmap = mapa_tipos
  )
)

write_csv(auditoria_contexto_mapa,
          file.path(output_dir, "auditoria_contexto_tcga_depmap.csv"))

auditoria_contexto_detalhada <- map_dfr(
  tipos_unicos,
  function(tipo_i) {
    filtrado_i <- filtrar_depmap_por_tipo_tcga(crisp_meta, tipo_i, mapa_tipos)
    
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

write_csv(auditoria_contexto_detalhada,
          file.path(output_dir, "auditoria_contexto_tcga_depmap_detalhada.csv"))

log_message("Auditoria do contexto tumoral exportada.")

# ======================================================================================
# 17) FILTRAGEM POR TIPO TUMORAL USANDO MAPA HIERÁRQUICO
# ======================================================================================

log_message("Aplicando filtro tumoral hierárquico por Tipo...")

res_por_tipo <- res_genes_total %>%
  group_split(Tipo, .keep = TRUE) %>%
  map_dfr(function(df_tipo) {
    tipo_i <- unique(df_tipo$Tipo)[1]
    
    if (length(tipo_i) == 0 || is.na(tipo_i) || nrow(df_tipo) == 0) return(tibble())
    
    filtrar_depmap_por_tipo_tcga(
      df            = df_tipo,
      tipo_tcga     = tipo_i,
      mapa_tcga_depmap = mapa_tipos
    )
  })

log_message("Registros após filtro hierárquico por Tipo: ", nrow(res_por_tipo))

auditoria_pos_filtro_analitico <- res_genes_total %>%
  count(Tipo, name = "n_antes_filtro") %>%
  full_join(
    res_por_tipo %>% count(Tipo, name = "n_depois_filtro"),
    by = "Tipo"
  ) %>%
  mutate(
    n_antes_filtro  = coalesce(n_antes_filtro,  0L),
    n_depois_filtro = coalesce(n_depois_filtro, 0L)
  ) %>%
  left_join(
    mapa_tipos %>%
      select(Tipo, nivel_filtro, lineage_alvo, lineage_subtype_alvo,
             lineage_sub_subtype_alvo, observacao),
    by = "Tipo"
  ) %>%
  arrange(Tipo)

write_csv(auditoria_pos_filtro_analitico,
          file.path(output_dir, "auditoria_pos_filtro_analitico.csv"))

# ======================================================================================
# 18) CLASSIFICAÇÃO INTERPRETATIVA
# ======================================================================================

res_por_tipo_class <- res_por_tipo %>%
  mutate(
    dependencia_classe      = classificar_dependency(dependency),
    interpretacao_biologica = interpretar_dependency(dependencia_classe)
  )

# ======================================================================================
# 19) TABELA PRINCIPAL POR LINHAGEM CELULAR
# ======================================================================================

tabela_linhagens_assinaturas <- res_por_tipo_class %>%
  left_join(
    mapa_tipos %>%
      select(Tipo, nivel_filtro, lineage_alvo, lineage_subtype_alvo, lineage_sub_subtype_alvo),
    by = "Tipo"
  ) %>%
  select(
    row_id,
    Tipo,
    nivel_filtro,
    lineage_alvo,
    lineage_subtype_alvo,
    lineage_sub_subtype_alvo,
    Signature_original,
    Signature,
    componente,
    componente_busca,
    gene_name,
    gene,
    entrez_id,
    depmap_id,
    cell_line,
    cell_line_name,
    lineage,
    lineage_subtype,
    lineage_sub_subtype,
    primary_disease,
    subtype_disease,
    primary_or_metastasis,
    dependency,
    dependencia_classe,
    interpretacao_biologica
  ) %>%
  arrange(Tipo, Signature_original, componente, dependency)

write_csv(tabela_linhagens_assinaturas,
          file.path(output_dir, "depmap_assinaturas_tabela_linhagens.csv"))

log_message("Tabela principal por linhagem exportada.")

# ======================================================================================
# 20) RESUMO POR COMPONENTE GÊNICO DENTRO DE CADA ASSINATURA E TIPO
# ======================================================================================

resumo_componente <- res_por_tipo_class %>%
  group_by(Tipo, Signature_original, Signature, componente, gene_name) %>%
  summarise(
    n_linhagens                = n(),
    media_dependency           = mean(dependency, na.rm = TRUE),
    mediana_dependency         = median(dependency, na.rm = TRUE),
    dp_dependency              = sd(dependency, na.rm = TRUE),
    min_dependency             = min(dependency, na.rm = TRUE),
    max_dependency             = max(dependency, na.rm = TRUE),
    n_forte_dependencia        = sum(dependencia_classe == "Forte dependência", na.rm = TRUE),
    n_dependencia_moderada     = sum(dependencia_classe == "Dependência moderada", na.rm = TRUE),
    n_linhagem_mantida         = sum(dependencia_classe == "Linhagem mantida / pouco efeito", na.rm = TRUE),
    n_vantagem_knockout        = sum(dependencia_classe == "Possível vantagem proliferativa após knockout", na.rm = TRUE),
    prop_forte_dependencia     = n_forte_dependencia    / n_linhagens,
    prop_dependencia_moderada  = n_dependencia_moderada / n_linhagens,
    prop_linhagem_mantida      = n_linhagem_mantida     / n_linhagens,
    prop_vantagem_knockout     = n_vantagem_knockout    / n_linhagens,
    classificacao_global = case_when(
      media_dependency < -1.5  ~ "Gene fortemente essencial no tecido analisado",
      media_dependency < -0.5  ~ "Gene moderadamente essencial no tecido analisado",
      media_dependency <= 0.5  ~ "Gene pouco essencial no tecido analisado",
      media_dependency > 0.5   ~ "Knockout pode ser vantajoso no tecido analisado",
      TRUE                     ~ "Indeterminado"
    ),
    .groups = "drop"
  ) %>%
  arrange(Tipo, Signature_original, media_dependency)

write_csv(resumo_componente,
          file.path(output_dir, "depmap_assinaturas_resumo_por_componente.csv"))

log_message("Resumo por componente exportado.")

# ======================================================================================
# 21) RESUMO POR ASSINATURA
# ======================================================================================

resumo_assinatura <- resumo_componente %>%
  group_by(Tipo, Signature_original, Signature) %>%
  summarise(
    n_componentes_avaliados                      = n(),
    media_das_medias_dependency                  = mean(media_dependency, na.rm = TRUE),
    mediana_das_medias_dependency                = median(media_dependency, na.rm = TRUE),
    n_componentes_forte_dependencia_media        = sum(media_dependency < -1.5, na.rm = TRUE),
    n_componentes_dependencia_moderada_ou_maior  = sum(media_dependency < -0.5, na.rm = TRUE),
    prop_componentes_forte_dependencia_media     = n_componentes_forte_dependencia_media / n_componentes_avaliados,
    prop_componentes_dependencia_moderada_ou_maior = n_componentes_dependencia_moderada_ou_maior / n_componentes_avaliados,
    veredito_assinatura = case_when(
      prop_componentes_forte_dependencia_media >= 0.5 ~
        "Assinatura enriquecida em genes fortemente essenciais no contexto tumoral",
      prop_componentes_dependencia_moderada_ou_maior >= 0.5 ~
        "Assinatura contém múltiplos genes com evidência de essencialidade",
      TRUE ~
        "Assinatura com baixa evidência global de essencialidade no DepMap CRISPR"
    ),
    .groups = "drop"
  ) %>%
  arrange(Tipo, media_das_medias_dependency)

write_csv(resumo_assinatura,
          file.path(output_dir, "depmap_assinaturas_resumo_por_assinatura.csv"))

log_message("Resumo por assinatura exportado.")

# ======================================================================================
# 22) TOP EXEMPLOS DE MAIOR ESSENCIALIDADE
# ======================================================================================

top_linhagens_mais_dependentes <- tabela_linhagens_assinaturas %>%
  group_by(Tipo, Signature_original, componente) %>%
  slice_min(order_by = dependency, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Tipo, Signature_original, componente, dependency)

write_csv(top_linhagens_mais_dependentes,
          file.path(output_dir, "depmap_assinaturas_top_linhagens_mais_dependentes.csv"))

log_message("Top linhagens mais dependentes exportado.")

# ======================================================================================
# 23) TABELA-MESTRA DE AUDITORIA POR ASSINATURA
# ======================================================================================

auditoria_assinaturas <- assinaturas_tbl %>%
  left_join(
    assinaturas_componentes %>%
      group_by(row_id, Tipo, Signature, Signature_original) %>%
      summarise(
        n_componentes_total        = n(),
        n_componentes_analisaveis  = sum(analisavel_no_depmap_crispr, na.rm = TRUE),
        n_componentes_nao_analisaveis = sum(!analisavel_no_depmap_crispr, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("row_id", "Tipo", "Signature", "Signature_original")
  ) %>%
  left_join(
    resumo_assinatura %>%
      select(Tipo, Signature, Signature_original, n_componentes_avaliados, veredito_assinatura),
    by = c("Tipo", "Signature", "Signature_original")
  ) %>%
  mutate(
    n_componentes_total           = coalesce(n_componentes_total,           0L),
    n_componentes_analisaveis     = coalesce(n_componentes_analisaveis,     0L),
    n_componentes_nao_analisaveis = coalesce(n_componentes_nao_analisaveis, 0L),
    n_componentes_avaliados       = coalesce(n_componentes_avaliados,       0L),
    veredito_assinatura           = coalesce(veredito_assinatura,
                                             "Sem componentes avaliáveis no DepMap CRISPR")
  )

write_csv(auditoria_assinaturas,
          file.path(output_dir, "auditoria_geral_assinaturas.csv"))

log_message("Auditoria geral por assinatura exportada.")

# ======================================================================================
# 23B) ASSINATURAS COM TODOS OS COMPONENTES ESSENCIAIS (threshold: média < -0.5)
# Critério : 100% dos componentes avaliados têm média de dependency < -0.5
# Checagem : componentes ausentes no DepMap são sinalizados para evitar falsos positivos
# ======================================================================================

THRESHOLD_ESSENCIALIDADE <- -0.5

log_message("============================================================")
log_message("ANÁLISE DE ASSINATURAS COM ESSENCIALIDADE MÍNIMA UNIVERSAL")
log_message("Threshold  : média de dependency < ", THRESHOLD_ESSENCIALIDADE)
log_message("Critério   : 100% dos componentes devem atingir o threshold")
log_message("============================================================")

# --- Cobertura: quantos componentes da assinatura foram encontrados no DepMap ---
cobertura_assinaturas <- auditoria_genes_assinaturas %>%
  group_by(Tipo, Signature, Signature_original) %>%
  summarise(
    n_componentes_na_assinatura  = n(),
    n_encontrados_no_depmap      = sum(encontrado_no_depmap,  na.rm = TRUE),
    n_nao_encontrados_no_depmap  = sum(!encontrado_no_depmap, na.rm = TRUE),
    genes_nao_encontrados        = paste(componente[!encontrado_no_depmap], collapse = "; "),
    cobertura_completa           = n_nao_encontrados_no_depmap == 0,
    .groups = "drop"
  )

# --- Classificação componente a componente ---
essencialidade_componentes <- resumo_componente %>%
  mutate(atinge_threshold = media_dependency < THRESHOLD_ESSENCIALIDADE)

# --- Avaliação por assinatura ---
avaliacao_assinaturas_essencialidade <- essencialidade_componentes %>%
  group_by(Tipo, Signature_original, Signature) %>%
  summarise(
    n_componentes_avaliados        = n(),
    n_componentes_acima_threshold  = sum(atinge_threshold,  na.rm = TRUE),
    n_componentes_abaixo_threshold = sum(!atinge_threshold, na.rm = TRUE),
    prop_componentes_essenciais    = n_componentes_acima_threshold / n_componentes_avaliados,
    media_dependency_mais_fraca    = max(media_dependency, na.rm = TRUE),
    media_dependency_mais_forte    = min(media_dependency, na.rm = TRUE),
    media_dependency_geral         = mean(media_dependency, na.rm = TRUE),
    assinatura_essencial           = if_else(
      sum(atinge_threshold, na.rm = TRUE) > 0,
      paste0("(", paste(gene_name[atinge_threshold], collapse = " + "), ")"),
      NA_character_
    ),
    componentes_que_falham         = paste(gene_name[!atinge_threshold], collapse = "; "),
    todos_essenciais               = all(atinge_threshold),
    .groups = "drop"
  ) %>%
  # join com cobertura para detectar genes ausentes no DepMap
  left_join(
    cobertura_assinaturas %>%
      select(Tipo, Signature, Signature_original,
             n_componentes_na_assinatura, n_nao_encontrados_no_depmap,
             genes_nao_encontrados, cobertura_completa),
    by = c("Tipo", "Signature", "Signature_original")
  ) %>%
  mutate(
    status_assinatura = case_when(
      n_componentes_avaliados == 0 ~
        "SEM_DADOS — nenhum componente foi avaliado no DepMap CRISPR",
      
      # cobertura incompleta: aviso independente do resultado
      !coalesce(cobertura_completa, TRUE) & todos_essenciais ~
        paste0(
          "APROVADA COM RESSALVA — ",
          n_nao_encontrados_no_depmap,
          " componente(s) ausente(s) no DepMap: ",
          genes_nao_encontrados
        ),
      
      !coalesce(cobertura_completa, TRUE) & !todos_essenciais ~
        paste0(
          "REJEITADA — ",
          n_componentes_acima_threshold, "/", n_componentes_avaliados,
          " componentes atingem o threshold; além disso ",
          n_nao_encontrados_no_depmap,
          " componente(s) ausente(s) no DepMap: ",
          genes_nao_encontrados
        ),
      
      todos_essenciais ~
        "APROVADA — todos os componentes atingem essencialidade mínima",
      
      n_componentes_acima_threshold == 0 ~
        "REJEITADA — nenhum componente atinge essencialidade mínima",
      
      TRUE ~
        paste0(
          "REJEITADA — apenas ",
          n_componentes_acima_threshold, "/", n_componentes_avaliados,
          " componentes atingem o threshold"
        )
    )
  ) %>%
  arrange(Tipo, desc(todos_essenciais), desc(prop_componentes_essenciais), media_dependency_geral)

# --- CSV de todas as assinaturas (aprovadas e rejeitadas) ---
write_csv(
  avaliacao_assinaturas_essencialidade,
  file.path(output_dir, "assinaturas_essencialidade_minima_universal.csv")
)

# --- CSV apenas das aprovadas (sem ressalva e com ressalva separados) ---
assinaturas_aprovadas_plenas <- avaliacao_assinaturas_essencialidade %>%
  filter(todos_essenciais & coalesce(cobertura_completa, TRUE))

assinaturas_aprovadas_com_ressalva <- avaliacao_assinaturas_essencialidade %>%
  filter(todos_essenciais & !coalesce(cobertura_completa, TRUE))

write_csv(assinaturas_aprovadas_plenas,
          file.path(output_dir, "assinaturas_todos_componentes_essenciais.csv"))

write_csv(assinaturas_aprovadas_com_ressalva,
          file.path(output_dir, "assinaturas_todos_componentes_essenciais_com_ressalva.csv"))

# --- Log detalhado por assinatura ---
log_message("")
log_message("------ DETALHAMENTO POR ASSINATURA ------")
log_message(sprintf("Total de assinaturas avaliadas      : %d", nrow(avaliacao_assinaturas_essencialidade)))
log_message(sprintf("Assinaturas APROVADAS (plenas)      : %d", nrow(assinaturas_aprovadas_plenas)))
log_message(sprintf("Assinaturas APROVADAS (com ressalva): %d", nrow(assinaturas_aprovadas_com_ressalva)))
log_message(sprintf("Assinaturas REJEITADAS              : %d",
                    nrow(avaliacao_assinaturas_essencialidade) -
                      nrow(assinaturas_aprovadas_plenas) -
                      nrow(assinaturas_aprovadas_com_ressalva)))
log_message("")

walk(seq_len(nrow(avaliacao_assinaturas_essencialidade)), function(i) {
  row <- avaliacao_assinaturas_essencialidade[i, ]
  
  log_message(sprintf(
    "[%s] Assinatura original: %s",
    row$Tipo, row$Signature_original
  ))
  log_message(sprintf("  Status : %s", row$status_assinatura))
  log_message(sprintf(
    "  Componentes avaliados: %d | Threshold atingido: %d | Prop: %.1f%%",
    row$n_componentes_avaliados,
    row$n_componentes_acima_threshold,
    row$prop_componentes_essenciais * 100
  ))
  log_message(sprintf(
    "  Dependency média (geral): %.3f | mais fraca: %.3f | mais forte: %.3f",
    row$media_dependency_geral,
    row$media_dependency_mais_fraca,
    row$media_dependency_mais_forte
  ))
  if (!is.na(row$assinatura_essencial) && row$assinatura_essencial != "") {
    log_message(sprintf("  Assinatura essencial : %s", row$assinatura_essencial))
  }
  if (!is.na(row$componentes_que_falham) && row$componentes_que_falham != "") {
    log_message(sprintf("  Falham : %s", row$componentes_que_falham))
  }
  if (!is.na(row$genes_nao_encontrados) && row$genes_nao_encontrados != "") {
    log_message(sprintf("  Ausentes no DepMap: %s", row$genes_nao_encontrados))
  }
  log_message("")
})

log_message("------ FIM DO DETALHAMENTO ------")
log_message("")
log_message("Arquivos exportados (bloco 23B):")
log_message("  - assinaturas_essencialidade_minima_universal.csv      (todas)")
log_message("  - assinaturas_todos_componentes_essenciais.csv         (aprovadas plenas)")
log_message("  - assinaturas_todos_componentes_essenciais_com_ressalva.csv (aprovadas com genes ausentes no DepMap)")

# ======================================================================================
# 23C) ESTRATIFICAÇÃO DAS ASSINATURAS POR GRAU DE ESSENCIALIDADE
# ======================================================================================

log_message("============================================================")
log_message("ANÁLISE DE ESTRATIFICAÇÃO POR GRAU DE ESSENCIALIDADE")
log_message("============================================================")

# Usamos a assinaturas_tbl como base para garantir que TODAS as assinaturas apareçam,
# mesmo aquelas que não puderam ser avaliadas por falta de genes no DepMap.
assinaturas_estratificadas <- assinaturas_tbl %>%
  left_join(
    avaliacao_assinaturas_essencialidade %>%
      select(Tipo, Signature_original, Signature, 
             n_componentes_avaliados, n_componentes_acima_threshold, 
             prop_componentes_essenciais, status_assinatura),
    by = c("Tipo", "Signature_original", "Signature")
  ) %>%
  mutate(
    grau_essencialidade = case_when(
      # Caso não tenha sido avaliada (zero componentes encontrados no DepMap)
      is.na(n_componentes_avaliados) | n_componentes_avaliados == 0 ~ "Sem dados (não avaliada no DepMap)",
      
      # Classificação baseada na proporção de componentes essenciais
      prop_componentes_essenciais == 1 ~ "Completamente essencial",
      prop_componentes_essenciais >= 0.75 ~ "Fortemente essencial",
      prop_componentes_essenciais >= 0.50 ~ "Moderadamente essencial",
      prop_componentes_essenciais > 0 ~ "Fracamente essencial",
      prop_componentes_essenciais == 0 ~ "Não essencial",
      
      TRUE ~ "Indeterminado"
    )
  ) %>%
  # Limpeza final: trocar NAs por descrições claras para evitar "lacunas" no CSV
  mutate(
    status_assinatura = coalesce(status_assinatura, "Nenhum componente gene/símbolo encontrado para esta assinatura"),
    n_componentes_avaliados = coalesce(n_componentes_avaliados, 0L),
    n_componentes_acima_threshold = coalesce(n_componentes_acima_threshold, 0L),
    prop_componentes_essenciais = coalesce(prop_componentes_essenciais, 0)
  )

write_csv(assinaturas_estratificadas,
          file.path(output_dir, "assinaturas_estratificadas_por_essencialidade.csv"))

resumo_estratificacao <- assinaturas_estratificadas %>%
  group_by(grau_essencialidade) %>%
  summarise(n_assinaturas = n(), .groups = "drop") %>%
  arrange(match(grau_essencialidade, c(
    "Completamente essencial",
    "Fortemente essencial",
    "Moderadamente essencial",
    "Fracamente essencial",
    "Não essencial",
    "Sem dados (não avaliada no DepMap)",
    "Indeterminado"
  )))

log_message("")
log_message("------ RESUMO DA ESTRATIFICAÇÃO ------")
walk(seq_len(nrow(resumo_estratificacao)), function(i) {
  row <- resumo_estratificacao[i, ]
  log_message(sprintf("  %s: %d", row$grau_essencialidade, row$n_assinaturas))
})
log_message("--------------------------------------")
log_message("")
log_message("Arquivo exportado (bloco 23C):")
log_message("  - assinaturas_estratificadas_por_essencialidade.csv")
log_message("")

# ======================================================================================
# 24) SAÍDA FINAL NO CONSOLE
# ======================================================================================

cat("\n================ RESULTADO FINAL ================\n")
cat("Diretório de saída:", normalizePath(output_dir), "\n\n")

cat("Resumo por assinatura:\n")
print(resumo_assinatura, n = min(nrow(resumo_assinatura), 30))

cat("\nResumo por componente:\n")
print(resumo_componente, n = min(nrow(resumo_componente), 30))

cat("\nAuditoria geral das assinaturas:\n")
print(auditoria_assinaturas, n = min(nrow(auditoria_assinaturas), 30))

cat("\nAuditoria do contexto TCGA -> DepMap:\n")
print(auditoria_contexto_mapa, n = nrow(auditoria_contexto_mapa))

cat("\nAssinaturas com todos os componentes essenciais (aprovadas plenas):\n")
print(assinaturas_aprovadas_plenas, n = min(nrow(assinaturas_aprovadas_plenas), 30))

cat("\nAssinaturas aprovadas com ressalva (genes ausentes no DepMap):\n")
print(assinaturas_aprovadas_com_ressalva, n = min(nrow(assinaturas_aprovadas_com_ressalva), 30))

cat("\nEstratificação das Assinaturas por Grau de Essencialidade:\n")
print(resumo_estratificacao)

cat("\n=================================================\n")

log_message("FIM DA ANÁLISE")