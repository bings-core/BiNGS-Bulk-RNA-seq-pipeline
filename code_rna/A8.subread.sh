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
star_dir="${preprocessed_dir}/star"
salmon_dir="${preprocessed_dir}/salmon"
subread_dir="${preprocessed_dir}/subread"
multiqc_dir="${preprocessed_dir}/multiqc"

# Create multiqc directory
mkdir -p "$multiqc_dir"

# Download MultiQC header template
wget -q --show-progress -P "$multiqc_dir" "https://ulukag01.u.hpc.mssm.edu/supporting_files/multiqc/multiqc_featurecounts_biotype_header.txt"

sample_metadata="${raw_dir}/sample_metadata/sample_metadata_rna.csv"

#########################################
### Detect or ask for organism name   ###
#########################################

annotation_dir="${supporting_files_dir}/annotation"
organism_dirs=($(ls -d "$annotation_dir"/*/ 2>/dev/null | xargs -n1 basename))

if [[ ${#organism_dirs[@]} -eq 1 ]]; then
    organism="${organism_dirs[0]}"
    echo "✅ Detected organism: $organism"
elif [[ ${#organism_dirs[@]} -eq 0 ]]; then
    echo "⚠️  No organism folder detected."
    read -p "Please specify your organism name exactly (e.g., mus_musculus, homo_sapiens): " organism
    mkdir -p "$annotation_dir/$organism"
    echo "📁 Created organism annotation folder: $annotation_dir/$organism"
else
    echo "⚠️  Multiple organism folders detected:"
    printf '%s\n' "${organism_dirs[@]}" | nl -w2 -s'. '
    echo ""

    while true; do
        read -p "Please type the correct organism name exactly as shown above: " organism
        if [[ " ${organism_dirs[@]} " =~ " $organism " ]]; then
            echo "✅ Confirmed organism: $organism"
            break
        else
            echo "❌ Invalid input. Please choose one of the detected organism names exactly."
        fi
    done

    echo "🧹 Deleting incorrect organism folders..."
    for dir in "${organism_dirs[@]}"; do
        if [[ "$dir" != "$organism" ]]; then
            rm -rf "$annotation_dir/$dir"
            echo "🗑️  Deleted: $dir"
        fi
    done
    echo "✅ Keeping only: $organism"
fi

organism_path="$annotation_dir/$organism"

#########################################
### Find filtered annotation GTF file ###
#########################################

filtered_gtf_file=$(find "$organism_path" -type f -name "*annotation.filtered.gtf" | head -n 1)

if [[ -z "$filtered_gtf_file" ]]; then
    echo "⚠️  No filtered annotation GTF file found under $organism_path."
    read -p "Please enter the path to your filtered annotation GTF: " filtered_gtf_file
    if [[ ! -f "$filtered_gtf_file" ]]; then
        echo "❗ Provided GTF file does not exist. Please re-run after preparing one."
        exit 1
    fi
else
    echo "✅ Found filtered GTF file: $filtered_gtf_file"
fi

#########################################
### Detect strandedness from Salmon   ###
#########################################

first_log_file=$(find "$salmon_dir" -type f -name "salmon_quant.log" | head -n 1)

if [[ -z "$first_log_file" ]]; then
    echo "❗ No salmon_quant.log file found in $salmon_dir."
    exit 1
fi

salmon_library_type=$(grep "Automatically detected most likely library type as" "$first_log_file" | awk '{print $NF}')

if [[ -z "$salmon_library_type" ]]; then
    echo "❗ Could not automatically detect library type from Salmon log."
    exit 1
fi

salmon_library_type=$(echo "$salmon_library_type" | tr '[:lower:]' '[:upper:]')
echo "🔎 Detected Salmon library type: $salmon_library_type"

case "$salmon_library_type" in
    SF|ISF)
        strandedness="1"
        ;;
    SR|ISR)
        strandedness="2"
        ;;
    U|IU)
        strandedness="0"
        ;;
    *)
        strandedness="0"
        ;;
esac

echo "🧬 Based on detection, suggested strandedness code: $strandedness"

#########################################
### Detect sample metadata columns    ###
#########################################

column_sampleid=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{for (i=1; i<=NF; i++) {gsub(/"/, "", $i); if ($i == "sample_id") print i;}}')
column_file1=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{for (i=1; i<=NF; i++) {gsub(/"/, "", $i); if ($i == "file_name_1") print i;}}')
column_file2=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{for (i=1; i<=NF; i++) {gsub(/"/, "", $i); if ($i == "file_name_2") print i;}}')

if [[ -z "$column_sampleid" || -z "$column_file1" || -z "$column_file2" ]]; then
    echo "❌ Error: Could not detect necessary columns in sample metadata CSV. Exiting."
    exit 1
fi

#########################################
### Handle project account            ###
#########################################

project_account_file="${supporting_files_dir}/project_account.txt"
if [[ -f "$project_account_file" ]]; then
    project_account=$(cat "$project_account_file")
    echo "✅ Using saved project account: $project_account"
else
    read -p "Enter your project account (e.g., acc_BiNGS_bulk): " project_account
    echo "$project_account" > "$project_account_file"
    echo "✅ Project account saved to $project_account_file"
fi

###################
### RUN SUBREAD ###
###################

mkdir -p "$preprocessed_dir/lsf"

# Loop through each sample and submit Subread jobs
tail -n +2 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' -v col1="$column_file1" -v col2="$column_file2" -v col3="$column_sampleid" '{
    gsub(/"/, "", $col1);
    gsub(/"/, "", $col2);
    gsub(/"/, "", $col3);
    print $col3, $col1, $col2;
}' | while read -r sample_id file1 file2; do

    FILE="$preprocessed_dir/lsf/subread_$sample_id.lsf"

    # Generate LSF submission script
    cat <<EOM > "$FILE"
#!/bin/bash
#BSUB -J subread_$sample_id
#BSUB -P $project_account
#BSUB -q premium
#BSUB -n 4
#BSUB -W 4:00
#BSUB -R "rusage[mem=10000]"
#BSUB -R "span[hosts=1]"
#BSUB -oo "$preprocessed_dir/lsf/subread_$sample_id.out"
#BSUB -eo "$preprocessed_dir/lsf/subread_$sample_id.err"
#BSUB -L /bin/bash

# Create output directory
mkdir -p "$subread_dir/featurecounts/$sample_id"

# Load subread module
module load subread/2.0.1

echo "🚀 Processing Sample: $sample_id"
echo "📂 Output Directory: $subread_dir/featurecounts/$sample_id"

bam_file=\$(ls "$star_dir/$sample_id/"*Aligned.sortedByCoord.out.bam 2>/dev/null)

# Check if BAM exists
if [[ ! -f "\$bam_file" ]]; then
    echo "❌ Error: Aligned BAM file missing for sample $sample_id. Skipping..."
    exit 1
fi

# Run featureCounts (single-end or paired-end)
if [[ -z "$file2" ]]; then
    echo "🔹 Running featureCounts for single-end reads..."
    featureCounts -B -C -g gene_type -t exon -T 6 -a "$filtered_gtf_file" -s $strandedness -o "$subread_dir/featurecounts/$sample_id/$sample_id.featureCounts.txt" "\$bam_file"
else
    echo "🔹 Running featureCounts for paired-end reads..."
    featureCounts -B -C -p -g gene_type -t exon -T 6 -a "$filtered_gtf_file" -s $strandedness -o "$subread_dir/featurecounts/$sample_id/$sample_id.featureCounts.txt" "\$bam_file"
fi

# Prepare MultiQC compatible biotype counts
cut -f 1,7 "$subread_dir/featurecounts/$sample_id/$sample_id.featureCounts.txt" | tail -n +3 | cat "$multiqc_dir/multiqc_featurecounts_biotype_header.txt" - > "$subread_dir/featurecounts/$sample_id/$sample_id.biotype_counts_mqc.tsv"

# Unload subread module
module unload subread/2.0.1
EOM

    # Submit the job
    bsub < "$FILE"
done

echo "✅ All Subread jobs submitted!"
echo "🛠️ Use 'bjobs' to check the status of submitted jobs."
echo "🚀 Wait for all Qualimap and Subread jobs to be finished before continuing with A9."
