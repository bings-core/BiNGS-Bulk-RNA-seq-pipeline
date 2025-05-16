###################################################
### BiNGS Bulk RNA - Pipeline - Transcriptomics ###
###################################################

# Get code and data directory paths based on the script location
script_path="$(realpath "${BASH_SOURCE[0]}")"
code_dir="$(dirname "$script_path")"
data_dir="${code_dir/code_rna/data_rna}"

# Define directory where the input files (e.g. raw fastq, sample metadata) will be stored
raw_dir="${data_dir}/raw"
sample_metadata="${raw_dir}/sample_metadata/sample_metadata_rna.csv"
# echo "🔎 Checking sample metadata file: $sample_metadata"

####################################
### Arrange your sample metadata ###
####################################

if [[ ! -f "${sample_metadata}" ]]; then
    while true; do
        read -p "Would you like to run this pipeline with the example dataset (GSE228989)? (yes/no): " use_example
        if [[ "$use_example" == "yes" || "$use_example" == "no" ]]; then
            break
        else
            echo "❌ Please answer yes or no."
        fi
    done

    if [[ "$use_example" == "yes" ]]; then
        echo "📁 Downloading the sample metadata and fastqs for GSE228989 bulk RNA-seq data..."
        sleep 3
        # Download sample metadata and fastqs
        wget -r -np -nH --cut-dirs=2 --reject "index.html*" https://ulukag01.u.hpc.mssm.edu/data_rna/raw/sample_metadata/ -P "$raw_dir/"
        wget -r -np -nH --cut-dirs=2 --reject "index.html*" https://ulukag01.u.hpc.mssm.edu/data_rna/raw/fastqs/ -P "$raw_dir/"

        echo "✅ Sample metadata and FASTQs copied to: $raw_dir."

        #######################################
        ### Download annotation ###
        #######################################
        supporting_files_dir="${code_dir/code_rna/supporting_files}"
        mkdir -p "$supporting_files_dir/annotation/mus_musculus"

        echo "📁 Downloading annotation files for mus_musculus..."
        sleep 3
        wget -r -np -nH --cut-dirs=3 --reject "index.html*" \
            https://ulukag01.u.hpc.mssm.edu/supporting_files/annotation/mus_musculus/ \
            -P "$supporting_files_dir/annotation/mus_musculus"

        echo "✅ Annotation files copied to: $supporting_files_dir/annotation/mus_musculus"

        echo "🎯 Example dataset setup complete! You can continue to A2."
        exit 1
    else
        echo "Great, you have another dataset in mind!"
        mkdir -p "${raw_dir}/sample_metadata"
        wget -r -np -nH --cut-dirs=3 --reject "index.html*" \
            https://ulukag01.u.hpc.mssm.edu/supporting_files/sample_metadata/ \
            -P "${raw_dir}/sample_metadata"
        sleep 3
        echo "📄 Here is an example sample metadata file for you to edit: $sample_metadata"
        sleep 3
        echo "✋ Please re-run this script after finalizing your metadata."
        sleep 3
        echo "🧠 While you're working on your metadata, let's start downloading the supporting files."
        sleep 3

        #######################################
        ### Ask for organism ###
        #######################################
        while true; do
            read -p "Which organism is your dataset from? (e.g., mus_musculus, homo_sapiens, drosophila_melanogaster): " organism
            if [[ -n "$organism" ]]; then
                break
            else
                echo "❌ Please enter a valid organism name."
            fi
        done

        supporting_files_dir="${code_dir/code_rna/supporting_files}"

        if [[ "$organism" == "mus_musculus" ]]; then
            mkdir -p "$supporting_files_dir/annotation/mus_musculus"
            echo "📁 Downloading annotation files for mus_musculus..."
            sleep 3
            wget -r -np -nH --cut-dirs=3 --reject "index.html*" \
                https://ulukag01.u.hpc.mssm.edu/supporting_files/annotation/mus_musculus/ \
                -P "$supporting_files_dir/annotation/mus_musculus"

            echo "✅ Annotation files copied to: $supporting_files_dir/annotation/mus_musculus"

        elif [[ "$organism" == "homo_sapiens" ]]; then
            mkdir -p "$supporting_files_dir/annotation/homo_sapiens"
            echo "📁 Downloading annotation files for homo_sapiens..."
            sleep 3
            wget -r -np -nH --cut-dirs=3 --reject "index.html*" \
                https://ulukag01.u.hpc.mssm.edu/supporting_files/annotation/homo_sapiens/ \
                -P "$supporting_files_dir/annotation/homo_sapiens"

            echo "✅ Annotation files copied to: $supporting_files_dir/annotation/homo_sapiens"

        else
            while true; do
                read -p "⚠️  Your data is from an organism other than human or mouse, correct? (yes/no): " confirm_other
                if [[ "$confirm_other" == "yes" || "$confirm_other" == "no" ]]; then
                    break
                else
                    echo "❌ Please answer yes or no."
                fi
            done

            if [[ "$confirm_other" == "yes" ]]; then
                mkdir -p "$supporting_files_dir/annotation/$organism"
                echo "📂 Created directory: $supporting_files_dir/annotation/$organism"
                echo "✋ You will have to provide paths for the STAR index and GTF files manually during the analysis."
            else
                echo "❌ Please re-run and type the organism correctly."
                exit 1
            fi
        fi
        exit 1
    fi

else
    while true; do
        read -p "Are you running this pipeline with the example dataset (GSE228989)? (yes/no): " use_example
        if [[ "$use_example" == "yes" || "$use_example" == "no" ]]; then
            break
        else
            echo "❌ Please answer yes or no."
        fi
    done

    if [[ "$use_example" == "yes" ]]; then
        echo "✅ Files should already be copied to: ${raw_dir}/fastqs."
        echo "✋ If you don't see all fastqs there, delete $raw_dir and start over."
        exit 1
    else
        while true; do
            read -p "Is your sample metadata finalized and are you ready to copy over fastqs? (yes/no): " finalized
            if [[ "$finalized" == "yes" || "$finalized" == "no" ]]; then
                break
            else
                echo "❌ Please answer yes or no."
            fi
        done

        if [[ "$finalized" != "yes" ]]; then
            echo "❌ Please finalize your metadata and re-run the script."
            exit 1
        fi
    fi
fi

#############################################
### Ask for user-supplied FASTQ directory ###
#############################################

read -p "📁 Please enter the full path to your directory containing all FASTQ files: " fastq_dir

# Check that path exists
if [[ ! -d "$fastq_dir" ]]; then
    echo "❌ Directory does not exist: $fastq_dir"
    exit 1
fi

# Run R script to complete metadata
echo "🧬 Processing sample metadata..."
ml R/4.2.0
Rscript "${code_dir}/do_not_run_complete_metadata.R" "$sample_metadata" "$fastq_dir"

if [[ $? -ne 0 ]]; then
    echo "❌ Error: Metadata completion failed. Check the error above."
    exit 1
fi

echo "✅ Metadata updated with full FASTQ paths and replicate info!"

########################################################
### Copy over raw fastq files to your data directory ###
########################################################

if [[ ! -d "${raw_dir}/fastqs" || -z "$(ls -A "${raw_dir}/fastqs" 2>/dev/null)" ]]; then
    mkdir -p "${raw_dir}/fastqs"

    # Auto-detect columns
    column_file1=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{
        for (i=1; i<=NF; i++) {
            gsub(/"/, "", $i);
            if ($i == "file_path_1") print i;
        }
    }')
    column_file2=$(head -1 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' '{
        for (i=1; i<=NF; i++) {
            gsub(/"/, "", $i);
            if ($i == "file_path_2") print i;
        }
    }')

    if [[ -z "$column_file1" || -z "$column_file2" ]]; then
        echo "❌ Error: File path columns not found in sample metadata!"
        exit 1
    fi

    total_files=$(tail -n +2 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' -v col1="$column_file1" -v col2="$column_file2" '{
        print $col1 "\n" $col2
    }' | grep -c '.')

    count=0
    tail -n +2 "$sample_metadata" | awk -v FPAT='([^,]*)|("[^"]+")' -v col1="$column_file1" -v col2="$column_file2" '
    {
        gsub(/"/, "", $col1);
        gsub(/"/, "", $col2);
        print $col1;
        print $col2;
    }' | while IFS= read -r fastq; do
        if [[ -n "$fastq" && -f "$fastq" ]]; then
            cp "$fastq" "${raw_dir}/fastqs/"
        fi
        ((count++))
        echo -ne "📂 Copying FASTQs: $count/$total_files files copied...\r"
    done

    echo -e "\n✅ All FASTQs copied to: ${raw_dir}/fastqs/"
    echo "🎯 You're ready for the next step (A2)!"
else
    echo "⚠️ It looks like you already have a fastqs directory: ${raw_dir}/fastqs"
    echo "🧹 If needed, remove it and start over."
    echo "🎯 Otherwise, you're ready for the next step (A2)!"
fi
