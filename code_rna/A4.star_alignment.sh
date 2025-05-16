###################################################
### BiNGS Bulk RNA - Pipeline - Transcriptomics ###
###################################################

# Get code and data directory paths based on the script location
script_path="$(realpath "${BASH_SOURCE[0]}")"
code_dir="$(dirname "$script_path")"
data_dir="${code_dir/code_rna/data_rna}"
supporting_files_dir="${code_dir/code_rna/supporting_files}"

# Define directory where the input files (e.g., raw fastq, sample metadata) will be stored
raw_dir="${data_dir}/raw"
preprocessed_dir="${data_dir}/preprocessed"

# Define directories for each package
trimgalore_dir="${preprocessed_dir}/trimgalore"
star_dir="${preprocessed_dir}/star"

sample_metadata="${raw_dir}/sample_metadata/sample_metadata_rna.csv"

# Detect common fastq extension (for user information, doesn't affect STAR)
fastq_extension=$(ls "${raw_dir}/fastqs" | grep -v '\.html$' | awk -F/ '{print $NF}' | awk -F. 'NF>1 {ext="."$NF; for (i=NF-1;i>1;i--) ext="."$i ext; print ext}' | sort | uniq -c | awk 'NR==1 {count=$1} $1!=count {exit} {print $2}')
echo "🔹 Common raw data file extension detected: $fastq_extension"

# Automatically detect column numbers
column_sampleid=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{for (i=1; i<=NF; i++) {gsub(/"/, "", $i); if ($i == "sample_id") print i;}}')
column_file1=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{for (i=1; i<=NF; i++) {gsub(/"/, "", $i); if ($i == "file_name_1") print i;}}')
column_file2=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{for (i=1; i<=NF; i++) {gsub(/"/, "", $i); if ($i == "file_name_2") print i;}}')

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

# Define organism path based on organism name
organism_path="$annotation_dir/$organism$( [[ "$organism" == "homo_sapiens" ]] && echo "/grch38_gencode_36" || ([[ "$organism" == "mus_musculus" ]] && echo "/grcm38_gencode_M25") )"

#############################################
### Check if STAR index and GTF exist     ###
#############################################

expected_star_index="$organism_path/star/2.7.5b/index"
expected_gtf_file="$organism_path/data/*annotation.gtf"

if [[ -d "$expected_star_index" && -n $(ls $expected_gtf_file 2>/dev/null) ]]; then
    echo "✅ STAR index and annotation GTF found for $organism."
else
    echo "❌ STAR index or annotation GTF files are missing for $organism. Please set them up first."
    exit 1
fi

# Set paths for STAR
star_index="$organism_path/star/2.7.5b/index"
gtf_file=$(ls "$organism_path/data/"*annotation.gtf 2>/dev/null | head -n1)

if [[ -z "$gtf_file" ]]; then
    echo "❌ Error: Could not find an annotation.gtf file under $organism_path/data/. Exiting."
    exit 1
fi

################
### RUN STAR ###
################

mkdir -p "$preprocessed_dir/lsf"

# Loop through each sample and submit STAR jobs
tail -n +2 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' -v col1="$column_file1" -v col2="$column_file2" -v col3="$column_sampleid" '{
    gsub(/"/, "", $col1);
    gsub(/"/, "", $col2);
    gsub(/"/, "", $col3);
    print $col3, $col1, $col2;
}' | while read -r sample_id file1 file2; do

    FILE="$preprocessed_dir/lsf/staralign_$sample_id.lsf"

    # Single-end vs paired-end
    if [[ -z "$file2" ]]; then
        trimmed_single="$trimgalore_dir/$sample_id/${sample_id}_trimmed.fq.gz"
        if [[ ! -f "$trimmed_single" ]]; then
            echo "❗ Warning: Missing $trimmed_single. Skipping $sample_id."
            continue
        fi
    else
        trimmed_read1="$trimgalore_dir/$sample_id/${sample_id}_val_1.fq.gz"
        trimmed_read2="$trimgalore_dir/$sample_id/${sample_id}_val_2.fq.gz"
        if [[ ! -f "$trimmed_read1" || ! -f "$trimmed_read2" ]]; then
            echo "❗ Warning: Missing $trimmed_read1 or $trimmed_read2. Skipping $sample_id."
            continue
        fi
    fi

    # Generate LSF submission script
    cat <<EOM > "$FILE"
#!/bin/bash
#BSUB -J staralign_$sample_id
#BSUB -P $project_account
#BSUB -q premium
#BSUB -n 8
#BSUB -W 4:00
#BSUB -R "rusage[mem=10000]"
#BSUB -R "span[hosts=1]"
#BSUB -oo "$preprocessed_dir/lsf/staralign_$sample_id.out"
#BSUB -eo "$preprocessed_dir/lsf/staralign_$sample_id.err"
#BSUB -L /bin/bash

# Create output directory
mkdir -p "$star_dir/$sample_id"

# Load STAR module
module load star/2.7.5b

echo "🚀 Processing Sample: $sample_id"
echo "📂 Output Directory: $star_dir/$sample_id"
EOM

    if [[ -z "$file2" ]]; then
        # Single-end STAR
        cat <<EOM >> "$FILE"
echo "🔹 Running STAR for single-end reads..."
STAR --genomeDir "$star_index" \
     --readFilesIn "$trimmed_single" \
     --runThreadN 8 \
     --outFileNamePrefix "$star_dir/$sample_id/$sample_id" \
     --sjdbGTFfile "$gtf_file" \
     --quantMode GeneCounts \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesCommand zcat
EOM
    else
        # Paired-end STAR
        cat <<EOM >> "$FILE"
echo "🔹 Running STAR for paired-end reads..."
STAR --genomeDir "$star_index" \
     --readFilesIn "$trimmed_read1" "$trimmed_read2" \
     --runThreadN 8 \
     --outFileNamePrefix "$star_dir/$sample_id/$sample_id" \
     --sjdbGTFfile "$gtf_file" \
     --quantMode GeneCounts \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesCommand zcat
EOM
    fi

    # Unload STAR
    echo "module unload star/2.7.5b" >> "$FILE"

    # Submit job
    bsub < "$FILE"

done

echo "✅ All STAR alignment jobs submitted!"
echo "🛠️  Use 'bjobs' to check the status of submitted jobs."
echo "🚀 You can run A5 right away."
