# ml R/4.2.0
# R

################################################
### BiNGS Bulk RNA - Pipeline - Transcriptomics
################################################

### Load Libraries ###
.libPaths(c("/hpc/packages/minerva-rocky9/rpackages/4.2.0/site-library", "/hpc/packages/minerva-rocky9/rpackages/bioconductor/3.15"))
library(DESeq2)
library(fst)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(rtracklayer)
library(dplyr)
library(reshape2)
library(patchwork)

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
dir.create(file.path(analysis_dir, "differential_expression"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "differential_expression"), recursive = TRUE, showWarnings = FALSE)

### Detect organism ###
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

### Load GTF file ###
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
    cat("\n⚠️ You picked a custom organism. Please provide full path to a GTF file (*.annotation.gtf): ")
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

### Load Data ###
sample_metadata <- fst::read_fst(file.path(compiled_data_dir, "dtl_sample_metadata.fst"))
dtl_salmon_gene_counts <- fst::read_fst(file.path(compiled_data_dir, "dtl_salmon_gene_counts.fst"))
dtl_salmon_gene_counts_medianratios <- fst::read_fst(file.path(compiled_data_dir, "dtl_salmon_gene_counts_medianratios.fst"))
rownames(dtl_salmon_gene_counts_medianratios) <- dtl_salmon_gene_counts_medianratios$gene_id
dtl_salmon_gene_counts_medianratios <- dtl_salmon_gene_counts_medianratios[, -1, drop = FALSE]

rownames(dtl_salmon_gene_counts) <- dtl_salmon_gene_counts$gene_id
dtl_counts <- as.matrix(round(dtl_salmon_gene_counts[, -1, drop = FALSE]))
mode(dtl_counts) <- "integer"

coldata <- data.frame(condition = sample_metadata$condition, sequencing_batch = sample_metadata$sequencing_batch)
rownames(coldata) <- sample_metadata$sample_id
coldata$condition <- as.factor(coldata$condition)

### Optionally remove low-quality samples ###
cat("\nWould you like to remove low-quality samples before DE analysis? (yes/no): ")
remove_samples_answer <- tolower(trimws(readLines(con = stdin(), n = 1)))

if (remove_samples_answer == "yes") {
    cat("\nAvailable sample IDs:")
    print(colnames(dtl_counts))
    cat("\nPlease enter samples to remove (comma-separated): ")
    low_quality_samples_input <- trimws(readLines(con = stdin(), n = 1))
    low_quality_samples <- trimws(unlist(strsplit(low_quality_samples_input, ",")))

    dtl_counts <- dtl_counts[, !(colnames(dtl_counts) %in% low_quality_samples)]
    coldata <- coldata[!(rownames(coldata) %in% low_quality_samples), ]
}

cat("\nRunning the DEG analysis with the following samples:\n")
print(colnames(dtl_counts))

if (!all(rownames(coldata) == colnames(dtl_counts))) stop("❌ Sample metadata and counts matrix do not match.")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = dtl_counts, colData = coldata, design = ~ condition)

### Check for existing DEG files ###
deg_csv_dir <- file.path(analysis_dir, "differential_expression")
existing_deg_files <- list.files(deg_csv_dir, pattern = "^condition_.*_degs\\.csv$", full.names = TRUE)

if (length(existing_deg_files) > 0) {
    cat("\n📁 Detected existing DEG results:\n")
    existing_comparisons <- gsub("^condition_(.*)_degs\\.csv$", "\\1", basename(existing_deg_files))
    print(existing_comparisons)

    cat("\nWould you like to run more differential expression comparisons? (yes/no): ")
    more_de <- tolower(trimws(readLines(con = stdin(), n = 1)))
} else {
    more_de <- "yes"
}

### Comparison loop ###
if (more_de == "yes") {
    repeat {
        cat("\nAvailable conditions:")
        print(unique(dds$condition))
        cat("\nEnter perturbation group (e.g., treated): ")
        perturbation <- trimws(readLines(con = stdin(), n = 1))
        cat("Enter reference group (e.g., control): ")
        reference <- trimws(readLines(con = stdin(), n = 1))

        comp_name <- paste0(perturbation, "_vs_", reference)
        comp <- c("condition", perturbation, reference)

        message("🔍 Running DESeq2: ", comp_name)
        dds <- DESeq(dds)
        res <- results(dds, contrast = comp)

        res <- tryCatch({ lfcShrink(dds, contrast = comp, type = "apeglm") }, error = function(e) { res })
        resOrdered <- res[order(res$log2FoldChange), ]

        de <- as.data.frame(resOrdered)
        de$deg_status <- "Not_DEG"
        de$deg_status[de$log2FoldChange > 1 & de$padj < 0.05] <- "Up_reg"
        de$deg_status[de$log2FoldChange < -1 & de$padj < 0.05] <- "Down_reg"
        de$gene_id <- rownames(de)
        de <- left_join(de, gtf_df, by = "gene_id")

        write.csv(de, file.path(analysis_dir, "differential_expression", paste0("condition_", comp_name, "_degs.csv")), row.names = FALSE)

        volcano <- ggplot(de, aes(x = log2FoldChange, y = -log10(padj), color = deg_status)) +
            geom_point(alpha = 0.7) +
            geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
            geom_text_repel(data = subset(de, deg_status != "Not_DEG" & padj < 0.01 & abs(log2FoldChange) > 2),
                            aes(label = gene_name), max.overlaps = 10) +
            scale_color_manual(values = c("Up_reg" = "red", "Down_reg" = "blue", "Not_DEG" = "grey")) +
            theme_minimal()

        pdf(file.path(figures_dir, "differential_expression", paste0("volcanoplot_", comp_name, ".pdf")), width = 8, height = 8)
        print(volcano)
        dev.off()

        cat("\nAnother comparison? (yes/no): ")
        if (tolower(trimws(readLines(con = stdin(), n = 1))) != "yes") break
    }
}

### Final DEG Heatmap Across All Comparisons ###

message("🎯 Now it's time to plot all the differentially expressed genes you identified in a heatmap!")

# Gather all DEG CSVs
deg_files <- list.files(file.path(analysis_dir, "differential_expression"), pattern = "^condition_.*_degs\\.csv$", full.names = TRUE)

if (length(deg_files) == 0) {
    message("⚠️ No DEG CSV files found. Skipping final combined heatmap.")
} else {
    deg_list <- lapply(deg_files, read.csv)

    # Extract comparison names
    comp_names <- gsub("^condition_(.*)_degs\\.csv$", "\\1", basename(deg_files))

    # Build a merged DEG status matrix
    all_genes <- unique(unlist(lapply(deg_list, function(x) x$gene_id)))
    deg_status_matrix <- data.frame(gene_id = all_genes)

    for (i in seq_along(deg_list)) {
        deg_status_matrix <- left_join(deg_status_matrix, deg_list[[i]][, c("gene_id", "deg_status")], by = "gene_id")
        colnames(deg_status_matrix)[ncol(deg_status_matrix)] <- comp_names[i]
    }

    rownames(deg_status_matrix) <- deg_status_matrix$gene_id
    deg_status_matrix$gene_id <- NULL

    # Map gene_id to gene_name
    gene_name_map <- setNames(gtf_df$gene_name, gtf_df$gene_id)
    gene_names <- ifelse(is.na(gene_name_map[rownames(deg_status_matrix)]),
                         rownames(deg_status_matrix),
                         gene_name_map[rownames(deg_status_matrix)])
    gene_names <- make.unique(gene_names)

    rownames(deg_status_matrix) <- gene_names

    # Save the combined DEG overlaps table
    deg_overlap_file <- file.path(analysis_dir, "differential_expression", "deg_overlaps.csv")
    write.csv(deg_status_matrix, file = deg_overlap_file, row.names = TRUE)
    message("📋 DEG overlaps table saved: ", deg_overlap_file)

    # Final genes to plot: union of all up/down genes
    final_genes_to_plot <- rownames(deg_status_matrix)[apply(deg_status_matrix, 1, function(x) any(x != "Not_DEG", na.rm = TRUE))]

    final_norm_counts <- dtl_salmon_gene_counts_medianratios[, colnames(dtl_counts), drop = FALSE]
    rownames(final_norm_counts) <- make.unique(ifelse(is.na(gene_name_map[rownames(final_norm_counts)]),
                                                    rownames(final_norm_counts),
                                                    gene_name_map[rownames(final_norm_counts)]))

    final_norm_counts <- as.data.frame(final_norm_counts)
    final_norm_counts <- final_norm_counts[intersect(rownames(final_norm_counts), final_genes_to_plot), , drop = FALSE]
    final_norm_counts <- as.matrix(final_norm_counts)
    mode(final_norm_counts) <- "numeric"

    # Filter problematic rows
    final_norm_counts <- final_norm_counts[rowSums(is.na(final_norm_counts) | is.infinite(final_norm_counts)) == 0, , drop = FALSE]
    final_norm_counts <- final_norm_counts[rowSums(final_norm_counts) != 0, , drop = FALSE]

    # Adjust row annotations too
    row_anno <- deg_status_matrix[final_genes_to_plot, , drop = FALSE]
    row_anno <- row_anno[rownames(final_norm_counts), , drop = FALSE]

    anno_colors <- list()
    for (comp in colnames(row_anno)) {
        anno_colors[[comp]] <- c(Up_reg = "red", Down_reg = "blue", Not_DEG = "grey")
    }

    pdf(file.path(figures_dir, "differential_expression", "combined_deg_heatmap.pdf"), width = 8, height = 8)
    pheatmap(final_norm_counts,
             scale = "row",
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             annotation_row = row_anno,
             annotation_colors = anno_colors,
             show_rownames = FALSE)
    dev.off()
    message("🎨 Final combined DEG heatmap saved: ", file.path(figures_dir, "differential_expression", "combined_deg_heatmap.pdf"))
}

### Gene-level Boxplots ###
cat("\nWould you like to plot expression for specific genes? Enter comma-separated gene names or press Enter to skip: ")
gene_input <- trimws(readLines(con = stdin(), n = 1))

if (nchar(gene_input) > 0) {
    gene_list <- trimws(unlist(strsplit(gene_input, ",")))
    gene_name_map <- setNames(gtf_df$gene_name, gtf_df$gene_id)
    gene_id_map <- setNames(gtf_df$gene_id, gtf_df$gene_name)

    # Map gene names to normalized count matrix rownames
    norm_counts <- dtl_salmon_gene_counts_medianratios[, colnames(dtl_counts), drop = FALSE]
    rownames(norm_counts) <- make.unique(ifelse(is.na(gene_name_map[rownames(norm_counts)]),
                                                rownames(norm_counts),
                                                gene_name_map[rownames(norm_counts)]))

    missing_genes <- setdiff(gene_list, rownames(norm_counts))
    if (length(missing_genes) > 0) {
        message("⚠️ These genes were not found and will be skipped: ", paste(missing_genes, collapse = ", "))
    }

    gene_list_found <- intersect(gene_list, rownames(norm_counts))

    if (length(gene_list_found) == 0) {
        message("❌ No valid genes found. Skipping boxplot generation.")
    } else {
        plot_list <- list()
        plot_data_long <- as.data.frame(norm_counts[gene_list_found, , drop = FALSE])
        plot_data_long$gene <- rownames(plot_data_long)
        plot_data_long <- reshape2::melt(plot_data_long, id.vars = "gene", variable.name = "sample", value.name = "value")
        plot_data_long$condition <- sample_metadata$condition[match(plot_data_long$sample, sample_metadata$sample_id)]
        plot_data_long$log2_value <- log2(plot_data_long$value + 1)

        for (g in gene_list_found) {
            p <- ggplot(subset(plot_data_long, gene == g), aes(x = condition, y = log2_value, fill = condition)) +
                geom_boxplot(outlier.shape = NA, alpha = 0.8) +
                geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
                theme_minimal() +
                labs(title = g, y = "log2(Norm Count + 1)", x = "") +
                theme(legend.position = "none",
                      plot.title = element_text(hjust = 0.5, face = "bold"))
            plot_list[[g]] <- p
        }

        # Determine layout
        n_genes <- length(plot_list)
        ncol <- ifelse(n_genes >= 4, 2, 1)
        nrow <- ceiling(n_genes / ncol)
        pdf_width <- 6 * ncol
        pdf_height <- 4 * nrow

        pdf(file.path(figures_dir, "differential_expression", "gene_boxplots.pdf"), width = pdf_width, height = pdf_height)
        print(wrap_plots(plot_list, ncol = ncol))
        dev.off()

        message("📦 Gene expression boxplots saved: ", file.path(figures_dir, "differential_expression", "gene_boxplots.pdf"))
    }
} else {
    message("👋 Skipping gene-level boxplots.")
}
