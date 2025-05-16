# ml R/4.2.0
# R

###############################################
### BiNGS Bulk RNA - Pipeline - Transcriptomics
###############################################

### Load necessary libraries ###
.libPaths(c("/hpc/packages/minerva-rocky9/rpackages/4.2.0/site-library", 
            "/hpc/packages/minerva-rocky9/rpackages/bioconductor/3.15"))
library(tximport)
library(fst)
library(DESeq2)
library(AnnotationDbi)
library(dplyr)

### Detect script location ###
script_path <- dirname(normalizePath(sys.frame(1)$ofile))
base_dir <- dirname(script_path)

### Define Directories ###
raw_dir <- file.path(base_dir, "data_rna/raw")
preprocessed_dir <- file.path(base_dir, "data_rna/preprocessed")
processed_dir <- file.path(base_dir, "data_rna/processed")
compiled_data_dir <- file.path(processed_dir, "compiled_data")
anno_dir <- file.path(base_dir, "supporting_files/annotation")
dir.create(compiled_data_dir, recursive = TRUE, showWarnings = FALSE)

### Load sample metadata ###
sample_metadata = read.csv(file.path(raw_dir, "sample_metadata/sample_metadata_rna.csv"))

# Update file paths conditionally
sample_metadata <- sample_metadata %>%
  mutate(
    file_path_1 = ifelse(!is.na(file_name_1) & nzchar(file_name_1), file.path(raw_dir, "fastqs", file_name_1), NA),
    file_path_2 = ifelse(!is.na(file_name_2) & nzchar(file_name_2), file.path(raw_dir, "fastqs", file_name_2), NA)
  )

# Define STAR alignment BAM suffixes
alignment_suffix = "Aligned.sortedByCoord.out.bam"
sample_metadata$alignment_bam = paste0(sample_metadata$sample_id, "_", alignment_suffix)

alignment_bam_files_star = sample_metadata$alignment_bam
alignment_sjout_files_star = gsub(paste0(alignment_suffix, "$"), "_SJ.out.tab", sample_metadata$alignment_bam)
alignment_count_files_star = gsub(paste0(alignment_suffix, "$"), "_ReadsPerGene.out.tab", sample_metadata$alignment_bam)

# Define paths for Subread and Salmon
gene_count_files_featurecounts = file.path(preprocessed_dir, "subread", paste0(sample_metadata$sample_id, ".featureCounts.txt"))
sample_metadata$transcript_quantification_files_salmon = file.path(preprocessed_dir, "salmon", sample_metadata$sample_id, "quant.sf")

# Save updated metadata
fst::write_fst(x = sample_metadata, path = file.path(processed_dir, "compiled_data/dtl_sample_metadata.fst"))

# Reload
dtl_sample_metadata = fst::read_fst(file.path(processed_dir, "compiled_data/dtl_sample_metadata.fst"))

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

################################################
### Load annotation data: txdb and ensembl mapping
################################################

if (organism == "homo_sapiens") {
    subfolder = "grch38_gencode_36"
    anno_subpath = file.path(anno_dir, organism, subfolder)
    txdb_url = "https://ulukag01.u.hpc.mssm.edu/supporting_files/annotation/homo_sapiens/grch38_gencode_36/txdb/homo_sapiens_grch38_gencode_36.sqlite"
    annotation_url = "https://ulukag01.u.hpc.mssm.edu/supporting_files/annotation/homo_sapiens/grch38_gencode_36/data/ensembl_transcript_annotations.csv"
} else if (organism == "mus_musculus") {
    subfolder = "grcm38_gencode_M25"
    anno_subpath = file.path(anno_dir, organism, subfolder)
    txdb_url = "https://ulukag01.u.hpc.mssm.edu/supporting_files/annotation/mus_musculus/grcm38_gencode_M25/txdb/mus_musculus_grcm38_gencode_M25.sqlite"
    annotation_url = "https://ulukag01.u.hpc.mssm.edu/supporting_files/annotation/mus_musculus/grcm38_gencode_M25/data/ensembl_transcript_annotations.csv"
} else {
    cat("\n⚠️ You picked a custom organism. Please provide:\n",
        "1. Path to a TxDb SQLite file (*.sqlite)\n",
        "2. Path to Ensembl transcript annotation file (ensembl_transcript_annotations.csv)\n")
    
    cat("\nEnter full path to the TxDb SQLite: ")
    txdb_path <- trimws(readLines(con = stdin(), n = 1))
    cat("\nEnter full path to Ensembl transcript annotation CSV: ")
    ensembl_transcript_annotation_path <- trimws(readLines(con = stdin(), n = 1))
}

if (organism %in% c("homo_sapiens", "mus_musculus")) {
    txdb_path <- file.path(anno_subpath, "txdb", basename(txdb_url))
    ensembl_transcript_annotation_path <- file.path(anno_subpath, "data", basename(annotation_url))
    
    dir.create(dirname(txdb_path), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(ensembl_transcript_annotation_path), recursive = TRUE, showWarnings = FALSE)
    
    if (!file.exists(txdb_path)) {
        message("⬇️  Downloading TxDb file...")
        download.file(url = txdb_url, destfile = txdb_path, mode = "wb")
    }
    if (!file.exists(ensembl_transcript_annotation_path)) {
        message("⬇️  Downloading Ensembl annotation file...")
        download.file(url = annotation_url, destfile = ensembl_transcript_annotation_path, mode = "wb")
    }
}

################################################
### Load transcript to gene mappings
################################################

ensembl_transcript_annotation = read.csv(ensembl_transcript_annotation_path, stringsAsFactors = FALSE)

txdb = AnnotationDbi::loadDb(txdb_path)
k = keys(txdb, keytype = "TXNAME")
tx2gene = AnnotationDbi::select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")

################################################
### Load Salmon transcript quantifications
################################################

txi_tx_salmon = tximport(sample_metadata$transcript_quantification_files_salmon, type = "salmon", txOut = TRUE)
txi_gene_salmon = summarizeToGene(txi_tx_salmon, tx2gene)

################################################
### Save Salmon outputs (.fst and .csv at once)
################################################

save_txi <- function(txi, prefix, id_col, sample_ids) {
  types <- c("abundance" = "tpm", "counts" = "counts", "length" = "length")
  
  for (key in names(types)) {
    value = types[[key]]
    dat = as.data.frame(txi[[key]], stringsAsFactors = FALSE)
    
    colnames(dat) <- sample_ids  # ✨ Always force correct sample IDs
    
    dat[[id_col]] = rownames(dat)
    dat = dat[, c(id_col, setdiff(names(dat), id_col))]  # Move ID first
    
    fst::write_fst(dat, paste0(prefix, value, ".fst"))
    write.csv(dat, paste0(prefix, value, ".csv"), row.names = FALSE)
  }
}

# Save both transcript-level and gene-level
save_txi(txi_tx_salmon, file.path(compiled_data_dir, "dtl_salmon_transcript_"), "transcript_id", sample_metadata$sample_id)
save_txi(txi_gene_salmon, file.path(compiled_data_dir, "dtl_salmon_gene_"), "gene_id", sample_metadata$sample_id)

################################################
### Generate normalized counts using DESeq2
################################################

dtl_salmon_gene_counts = fst::read_fst(file.path(compiled_data_dir, "dtl_salmon_gene_counts.fst"))
element_id = dtl_salmon_gene_counts[, 1, drop = FALSE]
dtl_counts = as.matrix(round(dtl_salmon_gene_counts[, -1]))
colnames(dtl_counts) = sample_metadata$sample_id  # ✨ Set sample IDs again here
mode(dtl_counts) = "integer"

coldata = sample_metadata[, "condition", drop = FALSE]
rownames(coldata) = sample_metadata$sample_id

dds = DESeqDataSetFromMatrix(countData = dtl_counts, colData = coldata, design = ~ condition)

# VST normalization
vsd <- vst(dds, blind = TRUE)
normalized_counts = as.data.frame(assay(vsd))
normalized_counts = cbind(element_id, normalized_counts)
fst::write_fst(normalized_counts, file.path(compiled_data_dir, "dtl_salmon_gene_counts_vsd.fst"))
write.csv(normalized_counts, file.path(compiled_data_dir, "dtl_salmon_gene_counts_vsd.csv"), row.names = FALSE)

# Median ratio normalization
dds = estimateSizeFactors(dds)
normalized_counts_median = as.data.frame(counts(dds, normalized = TRUE))
normalized_counts_median = cbind(element_id, normalized_counts_median)
fst::write_fst(normalized_counts_median, file.path(compiled_data_dir, "dtl_salmon_gene_counts_medianratios.fst"))
write.csv(normalized_counts_median, file.path(compiled_data_dir, "dtl_salmon_gene_counts_medianratios.csv"), row.names = FALSE)

message("✅ Finished saving all compiled data under: ", compiled_data_dir)
