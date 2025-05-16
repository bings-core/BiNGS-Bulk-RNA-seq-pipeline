###################################################
### BiNGS Bulk RNA - Pipeline - Transcriptomics ###
###################################################

# Get code and data directory paths based on the script location
script_path="$(realpath "${BASH_SOURCE[0]}")"
code_dir="$(dirname "$script_path")"
data_dir="${code_dir/code_rna/data_rna}"
supporting_files_dir="${code_dir/code_rna/supporting_files}"

# Define directory where the input files (e.g. raw fastq, sample metadata) will be stored
raw_dir="${data_dir}/raw"
preprocessed_dir="${data_dir}/preprocessed"

# Define directories for each package
fastqc_dir="${preprocessed_dir}/fastqc"  # ml fastqc/0.11.9
trimgalore_dir="${preprocessed_dir}/trimgalore"  # ml trim_galore/0.6.6
multiqc_dir="${preprocessed_dir}/multiqc"

# Create multiqc directory if it doesn't exist
mkdir -p "$multiqc_dir"

# Define path to sample metadata
sample_metadata="${raw_dir}/sample_metadata/sample_metadata_rna.csv"

######################################
### Download MultiQC config file   ###
######################################

multiqc_config="$multiqc_dir/multiqc_config_rnaseq_bulk_rna.yaml"

if [[ ! -f "$multiqc_config" ]]; then
    echo "⬇️  Downloading MultiQC configuration YAML..."
    wget -q --show-progress -P "$multiqc_dir" "https://ulukag01.u.hpc.mssm.edu/supporting_files/multiqc/multiqc_config_rnaseq_bulk_rna.yaml"
    echo "✅ Downloaded MultiQC config: $multiqc_config"
else
    echo "✅ Found existing MultiQC config: $multiqc_config"
fi

###################
### RUN MULTIQC ###
###################

echo "🚀 Running MultiQC to summarize all QC reports..."

module load python/3.7.3
module load openssl/1.0.2

multiqc -f \
    --config "$multiqc_config" \
    --outdir "$multiqc_dir" \
    --filename multiqc_report.html \
    "$preprocessed_dir"

module unload python/3.7.3
module unload openssl/1.0.2

echo "✅ MultiQC run complete!"
echo "📂 Report saved here: $multiqc_dir/multiqc_report.html"
