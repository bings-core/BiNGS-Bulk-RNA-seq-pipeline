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
star_dir="${preprocessed_dir}/star"
samtools_dir="${preprocessed_dir}/samtools"

sample_metadata="${raw_dir}/sample_metadata/sample_metadata_rna.csv"

# Automatically detect column numbers
column_sampleid=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{for (i=1; i<=NF; i++) {gsub(/"/, "", $i); if ($i == "sample_id") print i;}}')
column_file1=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{for (i=1; i<=NF; i++) {gsub(/"/, "", $i); if ($i == "file_name_1") print i;}}')
column_file2=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{for (i=1; i<=NF; i++) {gsub(/"/, "", $i); if ($i == "file_name_2") print i;}}')

# Ensure columns were detected
if [[ -z "$column_sampleid" || -z "$column_file1" || -z "$column_file2" ]]; then
    echo "❌ Error: Could not detect necessary columns in sample metadata CSV. Exiting."
    exit 1
fi

# Handle project account
project_account_file="${supporting_files_dir}/project_account.txt"
if [[ -f "$project_account_file" ]]; then
    project_account=$(cat "$project_account_file")
    echo "✅ Using saved project account: $project_account"
else
    read -p "Enter your project account (e.g., acc_BiNGS_bulk): " project_account
    echo "$project_account" > "$project_account_file"
    echo "✅ Project account saved to $project_account_file"
fi

#####################
### RUN SAMTOOLS  ###
#####################

mkdir -p "$preprocessed_dir/lsf"

# Loop through each sample and submit Samtools jobs
tail -n +2 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' -v col1="$column_file1" -v col2="$column_file2" -v col3="$column_sampleid" '{
    gsub(/"/, "", $col1);
    gsub(/"/, "", $col2);
    gsub(/"/, "", $col3);
    print $col3, $col1, $col2;
}' | while read -r sample_id file1 file2; do

    FILE="$preprocessed_dir/lsf/samtools_$sample_id.lsf"

    # Generate the LSF submission script
    cat <<EOM > "$FILE"
#!/bin/bash
#BSUB -J samtools_$sample_id
#BSUB -P $project_account
#BSUB -q premium
#BSUB -n 4
#BSUB -W 4:00
#BSUB -R "rusage[mem=10000]"
#BSUB -R "span[hosts=1]"
#BSUB -oo "$preprocessed_dir/lsf/samtools_$sample_id.out"
#BSUB -eo "$preprocessed_dir/lsf/samtools_$sample_id.err"
#BSUB -L /bin/bash

# Create output directory
mkdir -p "$samtools_dir/$sample_id"

# Load Samtools module
module load samtools/1.11

echo "🚀 Processing Sample: $sample_id"
echo "📂 Output Directory: $samtools_dir/$sample_id"

# Find the BAM file
bam_file=\$(ls "$star_dir/$sample_id/"*Aligned.sortedByCoord.out.bam 2>/dev/null)

# Check if BAM exists
if [[ ! -f "\$bam_file" ]]; then
    echo "❌ Error: Aligned BAM file missing for sample $sample_id. Skipping..."
    exit 1
fi

# Run Samtools operations
samtools index "\$bam_file"
samtools stats "\$bam_file" > "$samtools_dir/$sample_id/${sample_id}.stats"
samtools idxstats "\$bam_file" > "$samtools_dir/$sample_id/${sample_id}.idxstats"
samtools flagstat "\$bam_file" > "$samtools_dir/$sample_id/${sample_id}.flagstat"

# Unload Samtools module
module unload samtools/1.11
EOM

    # Submit the job
    bsub < "$FILE"

done

echo "✅ All Samtools jobs submitted!"
echo "🛠️  Use 'bjobs' to check the status of submitted jobs."
echo "🚀 Wait for all Samtools jobs to be finished before continuing with A7."
