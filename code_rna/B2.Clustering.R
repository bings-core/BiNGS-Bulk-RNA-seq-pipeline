# ml R/4.2.0
# R

################################################
### BiNGS Bulk RNA - Pipeline - Transcriptomics
################################################

### Load Libraries ###
.libPaths(c("/hpc/packages/minerva-rocky9/rpackages/4.2.0/site-library", "/hpc/packages/minerva-rocky9/rpackages/bioconductor/3.15"))
library(DESeq2)
library(fst)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(rtracklayer)

### Detect script location ###
script_path <- dirname(normalizePath(sys.frame(1)$ofile))
base_dir <- dirname(script_path)

### Define Directories ###
raw_dir <- file.path(base_dir, "data_rna/raw")
preprocessed_dir <- file.path(base_dir, "data_rna/preprocessed")
processed_dir <- file.path(base_dir, "data_rna/processed")
compiled_data_dir <- file.path(processed_dir, "compiled_data")
figures_dir <- file.path(processed_dir, "figures")
analysis_dir <- file.path(processed_dir, "analysis")
anno_dir <- file.path(base_dir, "supporting_files/annotation")

### Create Required Directories ###
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "clustering"), recursive = TRUE, showWarnings = FALSE)

#################################################
### Detect organism automatically from folders
#################################################

detected_organisms = list.dirs(anno_dir, recursive = FALSE, full.names = FALSE)

if (length(detected_organisms) == 1) {
    organism <- detected_organisms[1]
    message(paste0("✅ Detected organism: ", organism))
} else if (length(detected_organisms) == 0) {
    stop("❌ No organism folders detected under annotation/. Please prepare annotation files first.")
} else {
    message("⚠️ Multiple organisms detected:")
    for (org in detected_organisms) message("- ", org)
    
    repeat {
        cat("\nPlease enter the organism you want to use exactly as listed above: ")
        organism <- tolower(trimws(readLines(con = stdin(), n = 1)))
        if (organism %in% detected_organisms) {
            message(paste0("✅ Using organism: ", organism))
            break
        } else {
            message("❌ Invalid input. Try again.")
        }
    }
}

#################################################
### Load annotation GTF
#################################################

message("📜 Loading annotation GTF...")

if (organism == "homo_sapiens") {
    subfolder = "grch38_gencode_36"
    gtf_file <- file.path(anno_dir, organism, subfolder, "data", "gencode.v36.annotation.gtf")
    gtf_url <- "https://ulukag01.u.hpc.mssm.edu/supporting_files/annotation/homo_sapiens/grch38_gencode_36/data/gencode.v36.annotation.gtf"
} else if (organism == "mus_musculus") {
    subfolder = "grcm38_gencode_M25"
    gtf_file <- file.path(anno_dir, organism, subfolder, "data", "gencode.vM25.annotation.gtf")
    gtf_url <- "https://ulukag01.u.hpc.mssm.edu/supporting_files/annotation/mus_musculus/grcm38_gencode_M25/data/gencode.vM25.annotation.gtf"
} else {
    cat("\n⚠️ You picked a custom organism.\nPlease provide full path to a GTF file (*.annotation.gtf): ")
    gtf_file <- trimws(readLines(con = stdin(), n = 1))
    gtf_url <- NULL
}

if (!file.exists(gtf_file) && !is.null(gtf_url)) {
    message("⬇️ Downloading GTF file from: ", gtf_url)
    dir.create(dirname(gtf_file), recursive = TRUE, showWarnings = FALSE)
    download.file(url = gtf_url, destfile = gtf_file, mode = "wb")
}

gtf_data <- import(gtf_file)
gtf_df <- as.data.frame(gtf_data) %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  distinct(gene_id, .keep_all = TRUE)

#################################################
### Load compiled data
#################################################

sample_metadata <- fst::read_fst(file.path(compiled_data_dir, "dtl_sample_metadata.fst"))
dtl_salmon_gene_counts <- fst::read_fst(file.path(compiled_data_dir, "dtl_salmon_gene_counts.fst"))
dtl_salmon_gene_counts_medianratios <- fst::read_fst(file.path(compiled_data_dir, "dtl_salmon_gene_counts_medianratios.fst"))

# Prepare Counts Matrix
if ("gene_id" %in% colnames(dtl_salmon_gene_counts_medianratios)) {
  rownames(dtl_salmon_gene_counts_medianratios) <- dtl_salmon_gene_counts_medianratios$gene_id
  dtl_salmon_gene_counts_medianratios$gene_id <- NULL
}

dtl_counts <- as.matrix(round(dtl_salmon_gene_counts_medianratios))
mode(dtl_counts) <- "integer"

# Prepare Sample Metadata
coldata <- sample_metadata[, c("condition", "sequencing_batch")]
rownames(coldata) <- sample_metadata$sample_id

if (!all(rownames(coldata) == colnames(dtl_counts))) {
    stop("❌ Error: rownames(coldata) must match colnames(counts)!")
} else {
    message("🎯 Sample metadata and counts matrix match.")
}

#################################################
### Create DESeq2 Object and PCA
#################################################

dds <- DESeqDataSetFromMatrix(countData = dtl_counts, colData = coldata, design = ~ condition)
vsd <- vst(dds, blind = TRUE)

# Prompt user for metadata columns to use in PCA coloring
meta_cols_available <- setdiff(colnames(sample_metadata), c("sample_id", "file_name_1", "file_name_2", "file_path_1", "file_path_2", "alignment_bam", "transcript_quantification_files_salmon"))

if (interactive()) {
  cat("\nMetadata columns available for PCA coloring:\n")
  print(meta_cols_available)

  cat("\nPlease enter the metadata column names (comma-separated) you'd like to use for PCA coloring (press Enter to use default: condition): ")
  user_input <- trimws(readLines(con = stdin(), n = 1))
  if (nchar(user_input) == 0) {
    pca_color_vars <- "condition"
  } else {
    pca_color_vars <- unlist(strsplit(user_input, ","))
    pca_color_vars <- trimws(pca_color_vars)
  }
} else {
  message("ℹ️ Not interactive: defaulting to PCA colored by 'condition'")
  pca_color_vars <- "condition"
}

# Run PCA and save plots per selected metadata column
if (length(pca_color_vars) > 0) {
  pcaData <- plotPCA(vsd, ntop = 1000, returnData = TRUE)
  pcaData$name <- rownames(pcaData)
  pcaData <- pcaData %>% dplyr::select(-any_of("condition"))
  pcaData <- dplyr::left_join(pcaData, sample_metadata, by = c("name" = "sample_id"))
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  for (var in pca_color_vars) {
    if (!var %in% colnames(pcaData)) {
      warning("❌ Skipping unknown column: ", var)
      next
    }

    plot_obj <- ggplot(pcaData, aes_string(x = "PC1", y = "PC2", color = var, label = "name")) +
      geom_point(size = 3) +
      geom_text(vjust = 2, size = 3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      theme_light()

    out_path <- file.path(figures_dir, "clustering", paste0("preqc_pca_colored_by_", var, ".pdf"))
    pdf(out_path, width = 8, height = 8)
    print(plot_obj)
    dev.off()
    message("✅ PCA plot colored by ", var, " saved at: ", out_path)
  }
}

# Sample distance heatmap using sample_id labels
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Use sample_id from metadata for labeling
sample_ids <- sample_metadata$sample_id[match(colnames(vsd), sample_metadata$sample_id)]
rownames(sampleDistMatrix) <- sample_ids
colnames(sampleDistMatrix) <- sample_ids

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf(file.path(figures_dir, "clustering", "preqc_distance_matrix.pdf"), width = 8, height = 8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()
message("✅ Distance matrix heatmap saved at: ", file.path(figures_dir, "clustering", "preqc_distance_matrix.pdf"))

### Ask user if they want to remove low quality samples and redo PCA ###
cat("\nWould you like to remove low-quality samples and re-make PCA/Distance matrix plots? (yes/no): ")
remove_samples_answer <- tolower(trimws(readLines(con = stdin(), n = 1)))

if (remove_samples_answer == "yes") {
  
  cat("\nHere are the available sample IDs:\n")
  print(colnames(vsd))
  
  cat("\nPlease enter the sample IDs to remove (comma-separated like IFE_1,HF_3): ")
  low_quality_samples_input <- trimws(readLines(con = stdin(), n = 1))
  low_quality_samples <- unlist(strsplit(low_quality_samples_input, ","))
  low_quality_samples <- trimws(low_quality_samples)
  
  # Filter vsd and coldata
  vsd_filtered <- vsd[, !(colnames(vsd) %in% low_quality_samples)]
  coldata_filtered <- coldata[!(rownames(coldata) %in% low_quality_samples), ]
  
  ### PCA Plots (Post-QC) ###
  pcaData_postqc <- plotPCA(vsd_filtered, ntop = 1000, returnData = TRUE)
  pcaData_postqc$name <- rownames(pcaData_postqc)
  pcaData_postqc <- pcaData_postqc %>% dplyr::select(-any_of("condition"))
  pcaData_postqc <- dplyr::left_join(pcaData_postqc, sample_metadata, by = c("name" = "sample_id"))
  percentVar_postqc <- round(100 * attr(pcaData_postqc, "percentVar"))

  # Prompt once more to re-specify coloring vars or reuse earlier
  if (interactive()) {
    cat("\nMetadata columns available for post-QC PCA coloring:\n")
    print(meta_cols_available)

    cat("\nEnter metadata columns for post-QC PCA coloring (press Enter to reuse previous: ", paste(pca_color_vars, collapse = ", "), "): ")
    postqc_input <- trimws(readLines(con = stdin(), n = 1))
    if (nchar(postqc_input) > 0) {
      pca_color_vars <- trimws(unlist(strsplit(postqc_input, ",")))
    }
  }

  for (var in pca_color_vars) {
    if (!var %in% colnames(pcaData_postqc)) {
      warning("❌ Skipping unknown column: ", var)
      next
    }
    pca_plot <- ggplot(pcaData_postqc, aes_string(x = "PC1", y = "PC2", color = var, label = "name")) +
      geom_point(size = 3) +
      geom_text(vjust = 2, size = 3) +
      xlab(paste0("PC1: ", percentVar_postqc[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar_postqc[2], "% variance")) +
      theme_light()

    out_path <- file.path(figures_dir, "clustering", paste0("postqc_pca_colored_by_", var, ".pdf"))
    pdf(out_path, width = 8, height = 8)
    print(pca_plot)
    dev.off()
    message("✅ Post-QC PCA plot colored by ", var, " saved at: ", out_path)
  }

  ### Sample Distance Matrix (Post-QC) ###
  sampleDists_postqc <- dist(t(assay(vsd_filtered)))
  sampleDistMatrix_postqc <- as.matrix(sampleDists_postqc)

  # Use sample_id from metadata for labeling
  sample_ids_postqc <- sample_metadata$sample_id[match(colnames(vsd_filtered), sample_metadata$sample_id)]
  rownames(sampleDistMatrix_postqc) <- sample_ids_postqc
  colnames(sampleDistMatrix_postqc) <- sample_ids_postqc

  pdf(file.path(figures_dir, "clustering", "postqc_distance_matrix.pdf"), width = 8, height = 8)
  pheatmap(sampleDistMatrix_postqc,
           clustering_distance_rows = sampleDists_postqc,
           clustering_distance_cols = sampleDists_postqc,
           col = colors,
           display_numbers = FALSE,
           number_color = "black",
           fontsize_number = 6)
  dev.off()
  message("✅ Post-QC sample distance matrix saved at: ", file.path(figures_dir, "clustering", "postqc_distance_matrix.pdf"))

  # Replace for further steps
  vsd <- vsd_filtered
  coldata <- coldata_filtered
}

### Ask user if they want batch correction ###
if ("sequencing_batch" %in% colnames(coldata) && length(unique(coldata$sequencing_batch)) > 1) {
  
  cat("\nWould you like to apply batch correction before PCA? (yes/no): ")
  batch_correct_answer <- tolower(trimws(readLines(con = stdin(), n = 1)))
  
  if (batch_correct_answer == "yes") {
    
    ### Perform Batch Correction using limma::removeBatchEffect ###
    library(limma)
    
    message("Performing Batch Correction using limma::removeBatchEffect")
    vsd_corrected_matrix <- removeBatchEffect(assay(vsd),
                                              batch = coldata$sequencing_batch)
    
    # Save batch corrected counts
    batch_corrected_fst <- file.path(compiled_data_dir, "dtl_salmon_gene_counts_vsd_batchcorrected.fst")
    batch_corrected_csv <- file.path(compiled_data_dir, "dtl_salmon_gene_counts_vsd_batchcorrected.csv")
    
    fst::write_fst(as.data.frame(vsd_corrected_matrix), batch_corrected_fst)
    write.csv(as.data.frame(vsd_corrected_matrix), batch_corrected_csv, row.names = TRUE)
    
    message("✅ Batch corrected counts saved as FST:")
    print(batch_corrected_fst)
    message("✅ Batch corrected counts saved as CSV:")
    print(batch_corrected_csv)
    
    ### PCA Plot after Batch Correction ###
    vsd_corrected_for_pca <- vsd
    assay(vsd_corrected_for_pca) <- vsd_corrected_matrix
    
    pcaData_batchcorrected <- plotPCA(vsd_corrected_for_pca, ntop = 1000, returnData = TRUE)
    pcaData_batchcorrected$name <- rownames(pcaData_batchcorrected)
    pcaData_batchcorrected <- pcaData_batchcorrected %>% dplyr::select(-any_of("condition"))
    pcaData_batchcorrected <- dplyr::left_join(pcaData_batchcorrected, sample_metadata, by = c("name" = "sample_id"))
    percentVar_batchcorrected <- round(100 * attr(pcaData_batchcorrected, "percentVar"))
    
    pca_batchcorrected_plot <- ggplot(pcaData_batchcorrected, aes(x = PC1, y = PC2, color = condition, label = name)) +
      geom_point(size = 3) +
      geom_text(vjust = 2, size = 3) +
      xlab(paste0("PC1: ", percentVar_batchcorrected[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar_batchcorrected[2], "% variance")) +
      theme_light()
    
    pdf(file.path(figures_dir, "clustering", "postqc_pca_batchcorrected_colored_by_condition.pdf"), width = 8, height = 8)
    print(pca_batchcorrected_plot)
    dev.off()
    
    message("✅ Batch-corrected PCA plot saved:")
    print(file.path(figures_dir, "clustering", "postqc_pca_batchcorrected_colored_by_condition.pdf"))
    
    pca_batchcorrected_plot_batch <- ggplot(pcaData_batchcorrected, aes(x = PC1, y = PC2, color = sequencing_batch, label = name)) +
    geom_point(size = 3) +
    geom_text(vjust = 2, size = 3) +
    xlab(paste0("PC1: ", percentVar_batchcorrected[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar_batchcorrected[2], "% variance")) +
    theme_light()
  
    pdf(file.path(figures_dir, "clustering", "postqc_pca_batchcorrected_colored_by_batch.pdf"), width = 8, height = 8)
    print(pca_batchcorrected_plot_batch)
    dev.off()
    message("✅ Batch-corrected PCA plot (batch) saved:")
    print(file.path(figures_dir, "clustering", "postqc_pca_batchcorrected_colored_by_batch.pdf"))
    cat("\n⚠️ This plot is just to see if you have batch effect that needs to be corrected. The C and D scripts are not capable of running batch corrected analysis.")

  } else {
    message("Skipping batch corrected PCA for now as per user request.")
  }
} 
