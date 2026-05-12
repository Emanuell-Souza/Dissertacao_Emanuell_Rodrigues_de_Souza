## ==========================================================
## spacial_signature_ORIGINAL.R (Reconstruído)
## Análise Visium 10x — Loop Cânceres/Amostras/Assinaturas
## Output para: D:/RCDSurvXai/stRNA
## ==========================================================
req <- c("Seurat","SeuratObject","sctransform","ggplot2","patchwork",
         "R.utils","readxl","readr","dplyr","stringr","Matrix")
to_install <- req[!sapply(req, requireNamespace, quietly=TRUE)]
if (length(to_install)) install.packages(to_install, dependencies=TRUE)
library(Seurat); library(SeuratObject); library(ggplot2); library(patchwork)
library(R.utils); library(dplyr); library(stringr); library(readxl)
library(readr); library(Matrix)

## 0) Caminhos
base_amostras  <- "E:/ST_analises/amostras"
base_results   <- "E:/ST_analises/resultados"
dir.create(base_results, showWarnings=FALSE, recursive=TRUE)
sig_table_path <- "E:/ST_analises/Assinaturas_Omicas_SuperLearner.tsv"
markers_path   <- "E:/ST_analises/Genes_markers.xlsx"
if (!file.exists(sig_table_path)) stop("Tabela de assinaturas não encontrada.")
if (!file.exists(markers_path))   stop("Tabela de marcadores não encontrada.")

expression_sig_all <- read_tsv(sig_table_path, show_col_types=FALSE)
markers_all        <- read_excel(markers_path)
cancers_to_run     <- unique(expression_sig_all$Tipo)

# --- Helpers ---
split_genes <- function(x) {
  if (is.null(x) || is.na(x)) return(character(0))
  x <- as.character(x)
  x <- str_replace_all(x, '[\u201c\u201d\u201e\u201f\u2018\u2019""\'\'()]', "")
  g <- unlist(strsplit(x, "[+;,|\\s]+", perl=TRUE))
  g <- toupper(trimws(g)); g <- g[nzchar(g)]; unique(g)
}
scale01 <- function(x) {
  r <- range(x, na.rm=TRUE)
  if (!is.finite(r[1])||!is.finite(r[2])||r[1]==r[2]) return(rep(NA_real_,length(x)))
  (x - r[1])/(r[2]-r[1])
}
apply_diverging_scale <- function(p) {
  p + scale_fill_gradientn(colours=c("#0033FF","#FFFFC5","#FF0000"),
                           limits=c(0,1), na.value="transparent")
}

## ========== LOOP DE CÂNCERES ==========
for (cancer in cancers_to_run) {
  cancer_dir <- file.path(base_amostras, cancer)
  if (!dir.exists(cancer_dir)) { cat("Dir não encontrado:", cancer, "\n"); next }

  cancer_res_dir <- file.path(base_results, cancer)
  dir.create(cancer_res_dir, recursive=TRUE, showWarnings=FALSE)
  cat("\n=========================================\n")
  cat("Processando Câncer:", cancer, "\n")

  cancer_status_file <- file.path(cancer_res_dir, "status_cancer.rds")
  if (file.exists(cancer_status_file)) {
    if (readRDS(cancer_status_file) == "completed") {
      cat("Câncer", cancer, "já processado. Pulando.\n"); next
    }
  }

  expression_sig <- expression_sig_all %>% filter(Tipo == cancer)
  markers <- markers_all %>% filter(Type == cancer)

  # Marker sigs (Excel)
  marker_sigs <- list()
  if (nrow(markers) > 0) {
    marker_sigs <- markers %>% group_by(Cell_type) %>%
      summarise(genes=list(unique(toupper(trimws(Gene)))), .groups="drop") %>%
      { setNames(.$genes, .$Cell_type) }
  }

  # Custom sigs (TSV)
  sig_df <- expression_sig %>%
    transmute(sig_name_raw=Nomenclature, gene_string=Signature) %>%
    mutate(sig_name_raw=ifelse(is.na(sig_name_raw)|!nzchar(sig_name_raw),
                               sprintf("Signature_%03d",dplyr::row_number()),sig_name_raw),
           sig_name=sig_name_raw %>% str_trim("both") %>%
             str_replace_all("[^A-Za-z0-9_.-]+","_") %>% 
             str_replace_all("_+","_") %>% str_trim("both"))
  if (any(duplicated(sig_df$sig_name))) sig_df$sig_name <- make.unique(sig_df$sig_name, sep=".")

  custom_sigs <- setNames(vector("list", nrow(sig_df)), sig_df$sig_name)
  for (i in seq_len(nrow(sig_df))) {
    gi <- split_genes(sig_df$gene_string[i])
    if (length(gi)>0) custom_sigs[[sig_df$sig_name[i]]] <- gi
  }
  custom_sigs <- custom_sigs[lengths(custom_sigs)>0]

  signatures <- c(marker_sigs, custom_sigs)
  signature_groups <- c(setNames(rep("Marker",length(marker_sigs)),names(marker_sigs)),
                        setNames(rep("Custom",length(custom_sigs)),names(custom_sigs)))

  # Descobrir amostras
  dirs_amostra <- list.dirs(cancer_dir, recursive=FALSE)
  dirs_amostra <- dirs_amostra[!basename(dirs_amostra) %in% c("drafts","_backup","Plots")]
  dirs_amostra <- dirs_amostra[!grepl("_(filtered|raw)$", basename(dirs_amostra), ignore.case=TRUE)]

  # Organizar .gz soltos
  gz_files <- list.files(cancer_dir, "\\.gz$", full.names=TRUE)
  if (length(gz_files)>0) {
    prefixes <- unique(sub("_(matrix|features|barcodes|scalefactors|tissue|aligned|detected).*$","",basename(gz_files)))
    for (p in prefixes) {
      if (!nzchar(p)) next
      p_dir <- file.path(cancer_dir, p)
      dir.create(p_dir, showWarnings=FALSE)
      p_files <- list.files(cancer_dir, paste0("^",p,"_.*\\.gz$"), full.names=TRUE)
      file.rename(p_files, file.path(p_dir, basename(p_files)))
      dirs_amostra <- unique(c(dirs_amostra, p_dir))
    }
  }
  if (length(dirs_amostra)==0) { cat("Nenhuma amostra para", cancer, "\n"); next }

  map_file <- file.path(cancer_res_dir, "amostras_mapping.txt")
  if (file.exists(map_file)) file.remove(map_file)
  amostra_idx <- 1

  ## ========== LOOP DE AMOSTRAS ==========
  for (amostra_dir in dirs_amostra) {
    cat("\n--- Amostra:", basename(amostra_dir), "(Amostra_", amostra_idx, ") ---\n")
    nome_amostra_output <- paste0("Amostra_", amostra_idx)
    sample_res_dir <- file.path(cancer_res_dir, nome_amostra_output)
    dir.create(sample_res_dir, recursive=TRUE, showWarnings=FALSE)
    cat(sprintf("%s\t%s\n", nome_amostra_output, basename(amostra_dir)), file=map_file, append=TRUE)

    setwd(amostra_dir); base_dir <- amostra_dir
    dir10x  <- file.path(base_dir,"filtered_feature_bc_matrix")
    dirspat <- file.path(base_dir,"spatial")
    dir.create(dir10x, showWarnings=FALSE, recursive=TRUE)
    dir.create(dirspat, showWarnings=FALSE, recursive=TRUE)

    # Descompactar
    for (f in list.files(base_dir,"\\.gz$",full.names=TRUE,recursive=FALSE)) {
      out <- sub("\\.gz$","",f)
      if (!file.exists(out)) gunzip(f, remove=FALSE, overwrite=TRUE)
    }

    # Mover arquivos (pick_best + move_all_matches do V4)
    pick_best <- function(paths) { info<-file.info(paths); paths[order(info$size,info$mtime,decreasing=TRUE)[1]] }
    move_all_matches <- function(pattern, dest_dir, dest_name, backup_dir) {
      hits <- list.files(base_dir,pattern,full.names=TRUE,recursive=FALSE,ignore.case=TRUE)
      hits <- hits[!dir.exists(hits)]; hits <- hits[basename(hits)!=dest_name]
      if (!length(hits)) return(invisible(FALSE))
      dir.create(backup_dir, showWarnings=FALSE, recursive=TRUE)
      best <- pick_best(hits); dest <- file.path(dest_dir, dest_name)
      if (file.exists(dest)) file.remove(dest)
      ok <- file.rename(best, dest)
      if (!ok) { file.copy(best, dest, overwrite=TRUE); file.remove(best) }
      others <- setdiff(hits, best)
      if (length(others)) for (o in others) file.rename(o, file.path(backup_dir, basename(o)))
      invisible(TRUE)
    }

    bk10x <- file.path(dir10x,"_backup")
    move_all_matches("(?:^|_)matrix\\.mtx$",   dir10x,"matrix.mtx",  bk10x)
    move_all_matches("(?:^|_)features\\.tsv$", dir10x,"features.tsv",bk10x)
    move_all_matches("(?:^|_)barcodes\\.tsv$", dir10x,"barcodes.tsv",bk10x)
    bksp <- file.path(dirspat,"_backup")
    move_all_matches("scalefactors.*json$",          dirspat,"scalefactors_json.json",bksp)
    move_all_matches("tissue_positions_list\\.csv$", dirspat,"tissue_positions_list.csv",bksp)
    move_all_matches("tissue_positions\\.csv$",      dirspat,"tissue_positions.csv",bksp)
    move_all_matches("tissue_hires_image\\.png$",    dirspat,"tissue_hires_image.png",bksp)
    move_all_matches("tissue_lowres_image\\.png$",   dirspat,"tissue_lowres_image.png",bksp)

    pos_list <- file.path(dirspat,"tissue_positions_list.csv")
    pos_new  <- file.path(dirspat,"tissue_positions.csv")
    
    # --- NOVO: Conversor de JSON + TIFF para HNSC ---
    hnsc_json <- list.files(base_dir, ".*V11Y11.*\\.json$", full.names=TRUE)[1]
    if (!file.exists(pos_new) && !file.exists(pos_list) && !is.na(hnsc_json)) {
      cat("  HNSC detectado: Extraindo coordenadas e convertendo TIFF para PNG (isso pode demorar uns segundos)...\n")
      tryCatch({
        json_data <- jsonlite::fromJSON(hnsc_json)
        spots <- if(!is.null(json_data$oligo)) json_data$oligo else if(!is.null(json_data$spots)) json_data$spots else json_data
        if (!is.null(spots) && nrow(spots) > 0) {
          # Ler barcodes reais para que o Seurat reconheça as células
          barcodes_file <- list.files(base_dir, "barcodes\\.tsv", full.names=TRUE, recursive=TRUE)[1]
          real_bcs <- if (!is.na(barcodes_file)) readLines(barcodes_file) else paste0("spot_", 1:nrow(spots))
          real_bcs <- gsub("-[0-9]+$", "-1", real_bcs) # padronizar formato
          
          # TRUQUE: O Seurat não aceita que a imagem tenha mais spots que a matriz.
          # Vamos criar o CSV de coordenadas APENAS para o número de barcodes que realmente existem.
          n_cells <- length(real_bcs)
          spots_to_use <- spots[1:min(n_cells, nrow(spots)), , drop=FALSE]
          
          coords_hnsc <- data.frame(
            barcode = real_bcs[1:nrow(spots_to_use)],
            in_tissue = 1,
            row = spots_to_use$row,
            col = spots_to_use$col,
            pxl_row = spots_to_use$imageY,
            pxl_col = spots_to_use$imageX
          )


          dir.create(dirspat, showWarnings=FALSE, recursive=TRUE)
          write.table(coords_hnsc, pos_new, sep=",", col.names=FALSE, row.names=FALSE, quote=FALSE)
          
          # Carregar TIFF pesado e converter
          tif_path <- list.files(base_dir, "\\.tif$", ignore.case=TRUE, full.names=TRUE)[1]
          if (!is.na(tif_path) && requireNamespace("magick", quietly=TRUE)) {

            img <- magick::image_read(tif_path)
            orig_w <- magick::image_info(img)$width
            target_w <- 1500
            scale_f <- target_w / orig_w
            
            img_low <- magick::image_scale(img, paste0(target_w, "x"))
            magick::image_write(img_low, file.path(dirspat, "tissue_lowres_image.png"), format="png")
            
            dia <- if("dia" %in% colnames(spots)) mean(spots$dia, na.rm=TRUE) else 112.9
            sf_json <- sprintf('{"tissue_hires_scalef": %f, "tissue_lowres_scalef": %f, "fiducial_diameter_fullres": 144, "spot_diameter_fullres": %f}',
                               scale_f, scale_f, dia)
            writeLines(sf_json, file.path(dirspat, "scalefactors_json.json"))
          }
        }
      }, error = function(e) cat("    ERRO na conversão do HNSC:", conditionMessage(e), "\n"))
    }


    if (file.exists(pos_list) && !file.exists(pos_new)) file.copy(pos_list,pos_new,overwrite=TRUE)
    if (file.exists(pos_list) && file.exists(pos_new)) file.rename(pos_list, file.path(bksp,"tissue_positions_list.csv"))
    hires_png <- file.path(dirspat,"tissue_hires_image.png")
    lowres_png <- file.path(dirspat,"tissue_lowres_image.png")
    if (file.exists(hires_png) && file.exists(lowres_png)) file.rename(hires_png, file.path(bksp,"tissue_hires_image.png"))

    # Carregar/criar objeto Seurat
    processed_rds <- file.path(base_dir, "visium_seurat_processed.rds")
    if (file.exists(processed_rds)) {
      cat("  Objeto já processado. Carregando.\n")
      obj <- tryCatch(readRDS(processed_rds), error = function(e) {
        cat("  ERRO ao ler RDS (corrompido?):", conditionMessage(e), "\n  Deletando para reprocessar na próxima execução.\n")
        file.remove(processed_rds)
        NULL
      })
      if (is.null(obj)) { cat("  ⚠️ Pulando amostra — rode o script novamente para reprocessar.\n"); next }
      # Atualizar objeto se necessário (compatibilidade Seurat v4 -> v5)
      obj <- tryCatch(UpdateSeuratObject(obj), error = function(e) {
        cat("  Aviso: UpdateSeuratObject falhou:", conditionMessage(e), "\n"); obj
      })
    } else {
      req10x  <- file.path(dir10x,  c("matrix.mtx","features.tsv","barcodes.tsv"))
      reqspat <- file.path(dirspat, c("scalefactors_json.json","tissue_lowres_image.png"))
      if (!all(file.exists(req10x))||!all(file.exists(reqspat))||!file.exists(pos_new)) {
        cat("  Arquivos essenciais não encontrados. Pulando.\n"); next
      }
      feat_path <- file.path(dir10x,"features.tsv")
      first_line <- readLines(feat_path, n=1)
      feat_col <- if(length(strsplit(first_line,"\\t")[[1]])>=2) 2 else 1
      counts <- ReadMtx(mtx=file.path(dir10x,"matrix.mtx"), features=feat_path,
                        cells=file.path(dir10x,"barcodes.tsv"), feature.column=feat_col,
                        cell.column=1, unique.features=TRUE)
      obj <- CreateSeuratObject(counts=counts, assay="Spatial", project="Visium")
      img <- Read10X_Image(image.dir=dirspat, image.name="tissue_lowres_image.png")
      DefaultAssay(img) <- "Spatial"; obj[["slice1"]] <- img
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^MT-")
      if ("in_tissue" %in% colnames(obj@meta.data)) obj <- subset(obj, subset=in_tissue==1)
      obj <- subset(obj, subset=nFeature_Spatial>200 & nFeature_Spatial<7000 & percent.mt<20)
      if (ncol(obj)<10) { cat("  Spots insuficientes após QC. Pulando.\n"); next }
      use_glm <- requireNamespace("glmGamPoi", quietly=TRUE)
      obj <- SCTransform(obj, assay="Spatial", method=if(use_glm)"glmGamPoi"else"poisson", verbose=FALSE)
      DefaultAssay(obj) <- "SCT"
      obj <- RunPCA(obj, npcs=50, verbose=FALSE)
      dims_use <- 1:min(30, ncol(obj)-1)
      obj <- FindNeighbors(obj, dims=dims_use, verbose=FALSE)
      obj <- FindClusters(obj, resolution=0.6, verbose=FALSE)
      obj <- RunUMAP(obj, dims=dims_use, verbose=FALSE)
      saveRDS(obj, processed_rds)
    }

    ## Criar dirs de output
    plots_dir  <- file.path(sample_res_dir,"Plots"); dir.create(plots_dir,recursive=TRUE,showWarnings=FALSE)
    coords_dir <- file.path(plots_dir,"coords");     dir.create(coords_dir,recursive=TRUE,showWarnings=FALSE)
    scores_dir <- file.path(sample_res_dir,"signature_scores"); dir.create(scores_dir,recursive=TRUE,showWarnings=FALSE)
    overlap_dir<- file.path(sample_res_dir,"overlap_custom_vs_celltypes"); dir.create(overlap_dir,recursive=TRUE,showWarnings=FALSE)
    det_dir    <- file.path(overlap_dir,"detailed_coords"); dir.create(det_dir,recursive=TRUE,showWarnings=FALSE)
    am_dir     <- file.path(plots_dir,"All_Markers"); dir.create(am_dir,recursive=TRUE,showWarnings=FALSE)

    completed_sigs_rds <- file.path(sample_res_dir,"completed_signatures.rds")
    completed_sigs <- if(file.exists(completed_sigs_rds)) readRDS(completed_sigs_rds) else character(0)

    assay_detect <- if("Spatial" %in% Assays(obj)) "Spatial" else DefaultAssay(obj)
    assay_score  <- if("SCT" %in% Assays(obj)) "SCT" else DefaultAssay(obj)
    # Compatível com SeuratObject v4 (slot=) e v5 (layer=)
    get_assay_safe <- function(obj, assay, layer_name) {
      tryCatch(
        GetAssayData(obj, assay=assay, layer=layer_name),
        error = function(e) GetAssayData(obj, assay=assay, slot=layer_name)
      )
    }
    mat_counts <- get_assay_safe(obj, assay_detect, "counts")
    mat_data   <- get_assay_safe(obj, assay_score,  "data")

    # Coordenadas
    imgs <- Images(obj)
    img_name <- if("slice1" %in% imgs) "slice1" else imgs[1]
    xy <- GetTissueCoordinates(obj, image=img_name)
    if (all(c("x","y") %in% colnames(xy))) {
      xy_df <- xy[,c("x","y"),drop=FALSE]
    } else if (all(c("imagecol","imagerow") %in% colnames(xy))) {
      xy_df <- data.frame(x=xy[,"imagecol"], y=xy[,"imagerow"]); rownames(xy_df)<-rownames(xy)
    } else {
      xy_df <- data.frame(x=xy[,1], y=xy[,2], row.names=rownames(xy))
    }
    xy_df <- xy_df[colnames(obj),,drop=FALSE]

    # Armazenar scores para overlap e all_markers
    marker_detect_list <- list()  # celltype -> logical vector de detected
    custom_detect_list <- list()  # sig -> logical vector de detected

    ## ========== LOOP DE ASSINATURAS ==========
    for (sig in names(signatures)) {
      if (sig %in% completed_sigs) { cat("  Sig", sig, "já concluída.\n"); next }
      genes <- unique(signatures[[sig]])
      genes_in_counts <- intersect(genes, rownames(mat_counts))
      genes_in_data   <- intersect(genes, rownames(mat_data))
      if (length(genes_in_counts)==0||length(genes_in_data)==0) {
        completed_sigs <- c(completed_sigs, sig); saveRDS(completed_sigs, completed_sigs_rds); next
      }
      if (!setequal(genes_in_counts,genes_in_data)) {
        common <- intersect(genes_in_counts,genes_in_data)
        genes_in_counts <- common; genes_in_data <- common
      }

      g_name <- signature_groups[[sig]]

      # Regra de detecção diferente por grupo:
      # - Marker (tipos celulares): >50% dos genes presentes no spot
      # - Custom (assinaturas ômicas): TODOS os genes presentes no spot
      n_genes <- length(genes_in_counts)
      genes_detected_per_spot <- Matrix::colSums(mat_counts[genes_in_counts,,drop=FALSE] > 0)
      if (g_name == "Marker") {
        detect <- genes_detected_per_spot > (n_genes * 0.5)   # estrito: > 50% (não >=)
      } else {
        detect <- genes_detected_per_spot == n_genes           # todos os genes
      }

      score_raw    <- Matrix::colMeans(mat_data[genes_in_data,,drop=FALSE])
      score_scaled <- scale01(score_raw)
      score_masked <- score_scaled; score_masked[!detect] <- NA_real_
      # col_name <- paste0(sig,"Detected_scaled")
      col_name <- paste0("Detected_scaled")
      obj@meta.data[[col_name]] <- score_masked

      g_dir <- file.path(plots_dir, g_name); dir.create(g_dir,recursive=TRUE,showWarnings=FALSE)
      sc_dir <- file.path(scores_dir, g_name); dir.create(sc_dir,recursive=TRUE,showWarnings=FALSE)
      c_dir <- file.path(coords_dir, g_name); dir.create(c_dir,recursive=TRUE,showWarnings=FALSE)

      # Salvar signature_scores TSV
      df_scores <- data.frame(barcode=colnames(obj), signature=sig, group=g_name,
                              score_raw=as.numeric(score_raw),
                              score_masked=ifelse(detect, as.numeric(score_raw), NA),
                              detected=detect, stringsAsFactors=FALSE)
      write_tsv(df_scores, file.path(sc_dir, paste0("scores_",sig,".tsv")))

      # Guardar detecção para overlap e all_markers
      if (g_name=="Marker") marker_detect_list[[sig]] <- detect
      if (g_name=="Custom") custom_detect_list[[sig]] <- detect

      n_valid <- sum(!is.na(score_masked))
      if (n_valid > 0) {
        p <- SpatialFeaturePlot(obj, features=col_name, 
                                pt.size.factor=3,   # ← adicione isso arroz
                                combine=FALSE)[[1]]
        p <- apply_diverging_scale(p) + 
          ggtitle(paste0(sig, "  (n=", n_valid, ")")) +
          theme(
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.title = element_text(size=9),
            legend.text  = element_text(size=8)
          )
        ggsave(file.path(g_dir, paste0("SpatialFeature_",sig,".pdf")), p, width=6, height=5, dpi=300)

        df_sig <- cbind(data.frame(barcode=colnames(obj), score=score_masked, stringsAsFactors=FALSE), xy_df) %>%
          dplyr::filter(!is.na(score))
        write_tsv(df_sig, file.path(c_dir, paste0("coords_",sig,".tsv")))
      }

      completed_sigs <- c(completed_sigs, sig)
      saveRDS(completed_sigs, completed_sigs_rds)
      saveRDS(obj, processed_rds)
      cat("  ✅ Sig", sig, "concluída.\n")
    }

    ## ========== ALL_MARKERS (Combinações de cell types) ==========
    if (length(marker_detect_list) >= 2) {
      cat("  Gerando All_Markers...\n")
      detect_mat <- do.call(cbind, lapply(names(marker_detect_list), function(ct) marker_detect_list[[ct]]))
      colnames(detect_mat) <- names(marker_detect_list)

      marker_simple <- apply(detect_mat, 1, function(row) {
        pos <- names(which(row))
        if (length(pos)==0) return("Sem marcador detectado")
        if (length(pos)==1) return(pos)
        if (length(pos)==2) return(paste(sort(pos), collapse=" / "))
        return("Misto >= 3 tipos")
      })
      marker_combination <- apply(detect_mat, 1, function(row) {
        pos <- names(which(row))
        if (length(pos)==0) return("Sem marcador detectado")
        paste(sort(pos), collapse=" / ")
      })
      n_marker_types <- rowSums(detect_mat)

      am_table <- data.frame(barcode=colnames(obj), marker_simple=marker_simple,
                             marker_combination=marker_combination, n_marker_types=n_marker_types,
                             xy_df, stringsAsFactors=FALSE)
      write_tsv(am_table, file.path(am_dir,"Spatial_All_CellType_Markers_table.tsv"))

      # Combinations summary
      combo_tab <- as.data.frame(table(marker_combination), stringsAsFactors=FALSE)
      colnames(combo_tab) <- c("combination","count")
      combo_tab <- combo_tab[order(-combo_tab$count),]
      write_tsv(combo_tab, file.path(am_dir,"Spatial_All_CellType_Markers_combinations_summary.tsv"))

      # Plot CLEAN (marker_simple)
      obj@meta.data[["marker_simple"]] <- marker_simple
      tryCatch({
        n_cats <- length(unique(marker_simple))
        pal <- rep_len(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                         "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
                         "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
                         "#B3B3B3","#1B9E77","#D95F02","#7570B3","#E7298A"), n_cats)
        names(pal) <- unique(marker_simple)
        p_clean <- SpatialDimPlot(obj, group.by="marker_simple", pt.size.factor=1.5, cols=pal) +
          ggtitle("All CellType Markers — CLEAN") +
          theme(plot.title=element_text(hjust=0.5,face="bold",size=10), legend.text=element_text(size=7))
        ggsave(file.path(am_dir,"Spatial_All_CellType_Markers_CLEAN_with_tissue.pdf"), p_clean, width=10, height=7)
      }, error=function(e) cat("  Erro plot CLEAN:", conditionMessage(e),"\n"))

      # Plot TOP COMBINATIONS
      obj@meta.data[["marker_combination"]] <- marker_combination
      tryCatch({
        top_combos <- names(sort(table(marker_combination[marker_combination!="Sem marcador detectado"]),decreasing=TRUE))[1:min(15,length(unique(marker_combination)))]
        mc_top <- ifelse(marker_combination %in% top_combos, marker_combination, "Outras combinações")
        obj@meta.data[["marker_combo_top"]] <- mc_top
        n_cats2 <- length(unique(mc_top))
        pal2 <- rep_len(pal, n_cats2); names(pal2) <- unique(mc_top)
        p_top <- SpatialDimPlot(obj, group.by="marker_combo_top", pt.size.factor=1.5, cols=pal2) +
          ggtitle("All CellType Markers — TOP COMBINATIONS") +
          theme(plot.title=element_text(hjust=0.5,face="bold",size=10), legend.text=element_text(size=7))
        ggsave(file.path(am_dir,"Spatial_All_CellType_Markers_TOP_COMBINATIONS_with_tissue.pdf"), p_top, width=10, height=7)
      }, error=function(e) cat("  Erro plot TOP:", conditionMessage(e),"\n"))
    }

    ## ========== OVERLAP CUSTOM vs CELLTYPES ==========
    if (length(custom_detect_list)>0 && length(marker_detect_list)>0) {
      cat("  Gerando overlap custom vs celltypes...\n")
      ov_rows <- list()
      abs_mat <- matrix(0L, nrow=length(custom_detect_list), ncol=length(marker_detect_list),
                        dimnames=list(names(custom_detect_list), names(marker_detect_list)))
      pct_mat <- abs_mat * 0.0

      for (csig in names(custom_detect_list)) {
        cdet <- custom_detect_list[[csig]]
        sig_det_dir <- file.path(det_dir, csig); dir.create(sig_det_dir, showWarnings=FALSE, recursive=TRUE)
        n_sig <- sum(cdet)
        for (mct in names(marker_detect_list)) {
          mdet <- marker_detect_list[[mct]]
          overlap <- cdet & mdet
          n_overlap <- sum(overlap); n_marker <- sum(mdet)
          abs_mat[csig, mct] <- n_overlap
          pct_sig <- if(n_sig>0) n_overlap/n_sig*100 else 0
          pct_mrk <- if(n_marker>0) n_overlap/n_marker*100 else 0
          pct_mat[csig, mct] <- pct_sig
          ov_rows[[length(ov_rows)+1]] <- data.frame(Signature=csig, CellType=mct,
            n_signature_spots=n_sig, n_marker_spots=n_marker, n_overlap_spots=n_overlap,
            pct_signature_overlap=pct_sig, pct_marker_overlap=pct_mrk, stringsAsFactors=FALSE)
          # detailed coords
          if (n_overlap > 0) {
            ov_barcodes <- colnames(obj)[overlap]
            ov_df <- data.frame(barcode=ov_barcodes, xy_df[ov_barcodes,,drop=FALSE], stringsAsFactors=FALSE)
            write_tsv(ov_df, file.path(sig_det_dir, paste0("overlap_",csig,"_vs_",mct,".tsv")))
          }
        }
      }
      ov_summary <- do.call(rbind, ov_rows)
      write_tsv(ov_summary, file.path(overlap_dir,"overlap_summary_table.tsv"))
      abs_df <- data.frame(Signature=rownames(abs_mat), abs_mat, check.names=FALSE)
      write_tsv(abs_df, file.path(overlap_dir,"matrix_overlap_absolute_counts.tsv"))
      pct_df <- data.frame(Signature=rownames(pct_mat), pct_mat, check.names=FALSE)
      write_tsv(pct_df, file.path(overlap_dir,"matrix_overlap_pct_relative_to_signature.tsv"))
    }

    ## Salvar status da amostra
    status_info <- list(cancer=cancer, amostra=nome_amostra_output, original=basename(amostra_dir),
                        n_spots=ncol(obj), n_signatures=length(signatures),
                        completed=length(completed_sigs), timestamp=Sys.time())
    saveRDS(status_info, file.path(sample_res_dir,"status_sample.rds"))

    amostra_idx <- amostra_idx + 1
  }

  saveRDS("completed", cancer_status_file)
  cat("\n🎉 Câncer", cancer, "finalizado.\n")
}

cat("\n🚀 Processamento completo.\n")
