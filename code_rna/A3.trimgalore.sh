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

# Define directories for each package in the data directory before running script
fastqc_dir="${preprocessed_dir}/fastqc" # ml fastqc/0.11.9
trimgalore_dir="${preprocessed_dir}/trimgalore" # ml trim_galore/0.6.6

mkdir -p "${preprocessed_dir}/lsf"

sample_metadata="${raw_dir}/sample_metadata/sample_metadata_rna.csv"

# Automatically detect column numbers
column_sampleid=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{for (i=1; i<=NF; i++) { gsub(/"/, "", $i); if ($i == "sample_id") print i; }}')
column_file1=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{for (i=1; i<=NF; i++) { gsub(/"/, "", $i); if ($i == "file_name_1") print i; }}')
column_file2=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{for (i=1; i<=NF; i++) { gsub(/"/, "", $i); if ($i == "file_name_2") print i; }}')

# Ensure columns were detected
if [[ -z "$column_sampleid" || -z "$column_file1" || -z "$column_file2" ]]; then
    echo "❌ Error: Could not detect necessary columns in sample metadata CSV. Exiting."
    exit 1
fi

# Define where to save project account
project_account_file="${supporting_files_dir}/project_account.txt"

# Check if project account is already saved
if [[ -f "$project_account_file" ]]; then
    project_account=$(cat "$project_account_file")
    echo "Using saved project account: $project_account"
else
    read -p "Enter your project account (e.g., acc_BiNGS_bulk): " project_account
    echo "$project_account" > "$project_account_file"
    echo "Project account saved to $project_account_file"
fi

#######################
### RUN TRIMGALORE! ###
#######################

# Loop through each row of sample metadata and submit TrimGalore jobs
tail -n +2 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' -v col1="$column_file1" -v col2="$column_file2" -v col3="$column_sampleid" '{
    gsub(/"/, "", $col1);
    gsub(/"/, "", $col2);
    gsub(/"/, "", $col3);
    print $col3, $col1, $col2;
}' | while read -r sample_id file1 file2; do

    # Ensure at least one file exists before submitting job
    if [[ -f "$raw_dir/fastqs/$file1" || -f "$raw_dir/fastqs/$file2" ]]; then
        FILE="$preprocessed_dir/lsf/trimgalore_$sample_id.lsf"

        # Generate the LSF submission script
        cat <<EOM > "$FILE"

#!/bin/bash
#BSUB -J trimgalore_$sample_id
#BSUB -P $project_account
#BSUB -q premium
#BSUB -n 6
#BSUB -W 4:00
#BSUB -R "rusage[mem=1000]"
#BSUB -R "span[hosts=1]"
#BSUB -oo "$preprocessed_dir/lsf/trimgalore_$sample_id.out"
#BSUB -eo "$preprocessed_dir/lsf/trimgalore_$sample_id.err"
#BSUB -L /bin/bash

# Create output directory
mkdir -p "$trimgalore_dir/$sample_id"

# Load module
module load trim_galore/0.6.6

# Run TrimGalore!
if [[ -n "$file1" && -n "$file2" ]]; then # paired-end
    trim_galore --fastqc --paired --cores 6 --gzip \
        --basename "${sample_id}" \
        --output_dir "$trimgalore_dir/$sample_id" \
        "$raw_dir/fastqs/$file1" "$raw_dir/fastqs/$file2"
elif [[ -n "$file1" && -z "$file2" ]]; then # single-end
    trim_galore --fastqc --cores 6 --gzip \
        --basename "${sample_id}" \
        --output_dir "$trimgalore_dir/$sample_id" \
        "$raw_dir/fastqs/$file1"
fi

# Unload module
module unload trim_galore/0.6.6
EOM

        # Submit the job
        bsub < "$FILE"
    else
        echo "Warning: One or both FASTQ files missing for sample: $sample_id"
    fi
done

echo "✅ All TrimGalore! jobs submitted!"
echo "🛠️  Use 'bjobs' to check the status of submitted jobs."
echo "🚀 Wait for all FastQC and TrimGalore! jobs to be finished before continuing with A4."
