# ml R/4.2.0
# R

################################################
### BiNGS Bulk RNA - Pipeline - Transcriptomics
################################################

### Load Libraries ###
.libPaths(c("/hpc/packages/minerva-rocky9/rpackages/4.2.0/site-library", "/hpc/packages/minerva-rocky9/rpackages/bioconductor/3.15"))
library(DESeq2)
library(msigdbr)
library(clusterProfiler)
library(fst)
library(dplyr)
library(ggplot2)
library(forcats)
library(pheatmap)
library(ggrepel)
library(enrichplot)

### Set Seed ###
set.seed(42)

### Detect Script Location ###
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
dir.create(file.path(analysis_dir, "functional_analysis"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "functional_analysis"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "functional_analysis", "gsea_enrichment_curves"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "functional_analysis", "gsea_barplots"), recursive = TRUE, showWarnings = FALSE)

### Detect Organism ###
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

### Load Annotation ###
if (organism == "homo_sapiens") {
    ensembl_transcript_annotation_path <- file.path(anno_dir, "homo_sapiens/grch38_gencode_36/data/ensembl_transcript_annotations.csv")
    msigdbr_species = "Homo sapiens"
} else if (organism == "mus_musculus") {
    ensembl_transcript_annotation_path <- file.path(anno_dir, "mus_musculus/grcm38_gencode_M25/data/ensembl_transcript_annotations.csv")
    msigdbr_species = "Mus musculus"
} else {
    stop("❌ Invalid organism.")
}

ensembl_transcript_annotation = read.csv(ensembl_transcript_annotation_path, stringsAsFactors = FALSE)

### Define Functional Analysis Gene Sets ###
functional_analysis_gene_sets = c(
    "C2:CP:KEGG: KEGG pathway database",
    "C2:CP:REACTOME: Reactome pathway database",
    "C2:CP:WIKIPATHWAYS: WikiPathways pathway database",
    "H: hallmark gene sets"
)

### Load MSigDB ###
msigdbr_all <- msigdbr(species = msigdbr_species)

### Select Gene Sets ###
functional_gene_sets <- msigdbr_all %>%
  filter(
    (gs_cat == "C2" & gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "CP:WIKIPATHWAYS")) |
    (gs_cat == "H")
  )

### Helper Functions ###
trim_gene_version = function(gene_id, split = ".") {
    sapply(strsplit(gene_id, split = split, fixed = TRUE), "[[", 1)
}

parse_comparison_name <- function(comp) {
  comp <- gsub("^condition_", "Comparison: ", comp)
  comp <- gsub("_vs_", " vs ", comp)
  return(comp)
}

### Load DEG Results ###
deg_dir <- file.path(analysis_dir, "differential_expression")
deg_files <- list.files(deg_dir, pattern = "_degs.csv$", full.names = TRUE)

pvalue_threshold <- 0.01

for (deg_file in deg_files) {
    comp_name <- gsub("_degs.csv", "", basename(deg_file))
    comp_clean <- parse_comparison_name(comp_name)
    message("📊 Processing DEG comparison: ", comp_name)

    res.df <- read.csv(deg_file)

    results_functional <- res.df
    results_functional$entrezgene_id <- ensembl_transcript_annotation$entrezgene_id[
        match(trim_gene_version(results_functional$gene_id), ensembl_transcript_annotation$ensembl_gene_id)
    ]

    results_functional <- results_functional %>%
        filter(is.finite(log2FoldChange), is.finite(pvalue), !is.na(entrezgene_id))

    results_functional$pvalue[results_functional$pvalue == 0] <- .Machine$double.eps
    results_functional$neg_log_10_p_value <- -log10(results_functional$pvalue)
    results_functional$combined_score <- results_functional$log2FoldChange * results_functional$neg_log_10_p_value

    res_gene_list <- results_functional$combined_score
    names(res_gene_list) <- results_functional$entrezgene_id
    res_gene_list <- sort(res_gene_list, decreasing = TRUE)

    ### Run GSEA ###
    gsea_results <- GSEA(
        geneList = res_gene_list,
        TERM2GENE = functional_gene_sets %>% select(gs_name, entrez_gene),
        verbose = TRUE,
        by = "fgsea",
        seed = TRUE,
        pvalueCutoff = 1,
        BPPARAM = BiocParallel::MulticoreParam(workers = 1)
    )

    gsea_result_path <- file.path(analysis_dir, "functional_analysis", paste0(comp_name, "_gsea_results.csv"))
    write.csv(gsea_results@result, gsea_result_path, row.names = FALSE)

    res_plotting <- gsea_results@result %>%
    filter(!is.na(NES)) %>%
    mutate(
        direction = ifelse(NES > 0, "Up", "Down"),
        neg_log10_padj = -log10(p.adjust),
        label = ifelse(p.adjust < pvalue_threshold, paste0("*", Description), Description)
    )

    ### Enrichment Curves for Significant Pathways ###
    gsea_sig_terms <- res_plotting %>% filter(p.adjust < pvalue_threshold)

    if (nrow(gsea_sig_terms) > 0) {
        for (term in gsea_sig_terms$ID) {
            pdf(file.path(figures_dir, "functional_analysis", "gsea_enrichment_curves", paste0(comp_name, "_", term, "_enrichment_curve.pdf")), width = 12, height = 6)
            print(gseaplot2(gsea_results, geneSetID = term, title = paste0(term, "\n", comp_clean, "\nNES=", round(gsea_sig_terms$NES[gsea_sig_terms$ID == term], 2), ", padj=", signif(gsea_sig_terms$p.adjust[gsea_sig_terms$ID == term], 2))))
            dev.off()
        }
    }

    ### Barplot: Top Up and Down GSEA ###
    top_up <- res_plotting %>% filter(NES > 0) %>% arrange(p.adjust) %>% head(10)
    top_down <- res_plotting %>% filter(NES < 0) %>% arrange(p.adjust) %>% head(10)

    combined <- bind_rows(top_down, top_up)

    if (nrow(combined) > 0) {

    combined <- combined %>%
      mutate(
        label = as.character(label),
        short_label = ifelse(
          nchar(label) > 60,
          paste0(substr(label, 1, 60), "..."),
          label
        ),
        short_label = ifelse(grepl("^\\*", label), paste0("*", short_label), short_label),
        NES_clamped = ifelse(NES > 1, 1, ifelse(NES < -1, -1, NES))  # Clamp NES
      )

    combined <- combined %>%
        arrange(NES)  # Negative NES first, positive NES later

    combined$short_label <- factor(combined$short_label, levels = combined$short_label)

    pdf(file.path(figures_dir, "functional_analysis", "gsea_barplots", paste0(comp_name, "_gsea_top10_barplot.pdf")), width = 15, height = 9)

    p <- ggplot(combined, aes(x = NES, y = short_label, fill = NES_clamped)) +
        geom_col() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        scale_fill_gradient2(
            low = "blue", mid = "white", high = "red", midpoint = 0,
            limits = c(-1, 1),
            name = "NES"
        ) +
        labs(
            x = "Normalized Enrichment Score (NES)",
            y = "Pathway",
            title = paste("Top Pathways -", comp_clean),
            caption = "Asterisked pathways are padj < 0.01 and statistically significant*"
        ) +
        theme_minimal(base_size = 14) +
        theme(
            axis.text.y = element_text(
            face = "bold"),
            plot.title = element_text(hjust = 0.5),
            plot.caption = element_text(size = 10, hjust = 0)
        ) +
        coord_cartesian(xlim = c(-2, 2))  # Small margin for bar edges

    print(p)
    dev.off()
}

}

message("✅ Functional analysis complete!")
message("📋 GSEA tables saved under: ", file.path(analysis_dir, "functional_analysis"))
message("🎨 Figures saved under: ", file.path(figures_dir, "functional_analysis"))
