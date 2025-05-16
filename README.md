
# BiNGS Bulk RNA-seq Pipeline  
*An End-to-End Bulk RNA-seq Analysis Pipeline*  
Developed by the [BiNGS Core](https://bings.mssm.edu/) at the Tisch Cancer Institute, Mount Sinai

---

## Table of Contents
- Setup
  - [Option 1: CodeServer (OnDemand) Setup](#option-1-codeserver-ondemand-setup)
  - [Option 2: VSCode Desktop App Setup](#option-2-vscode-desktop-app-setup)
- Running Scripts
  - [Bash Scripts](#bash-scripts)
  - [R Scripts](#r-scripts)
- [Good Practices](#-good-practices)
- [Troubleshooting](#%EF%B8%8F-troubleshooting)
  - ["403 Forbidden" error while running the A1 script](#403-forbidden-error-while-running-the-a1-script)
- [Questions?](#-questions)
- [Citation](#-citation)
  
## 🖥️ Setup on HPC

### Clone the Pipeline Repository

```bash
cd /your/working/directory/
git clone https://github.com/GulayBenguUlukaya/BiNGS-Bulk-RNA-seq-pipeline.git
```

---

## Option 1: CodeServer (OnDemand) Setup

1. **Start a session with:**
   - Queue: `premium`
   - Cores: `10`
   - Memory: `64 GB`
   - Walltime: `4 hrs`
   - ⚠️ Use your **own project account** instead of `acc_BiNGS_bulk`.

2. **Create and enter your project directory:**

```bash
mkdir ./rna_seq_pipeline_test
cd ./rna_seq_pipeline_test
git clone https://github.com/GulayBenguUlukaya/BiNGS-Bulk-RNA-seq-pipeline.git
```

> For steps after cloning, begin with `A1.prepare_files.sh` and follow the prompts.

---

## Option 2: VSCode Desktop App Setup

1. **Start a `screen` session and request interactive job:**

```bash
screen -S rnaseq
```

> ⚠️ Use your **own project account** instead of `acc_BiNGS_bulk`.

```bash
bsub -P <your_project_account> -q premium -n 10 -W 4:00 -R span[hosts=1] -R rusage[mem=64000] -Is /bin/bash
```

2. **Create and enter your project directory:**

```bash
mkdir ./rna_seq_pipeline_test
cd ./rna_seq_pipeline_test
git clone https://github.com/GulayBenguUlukaya/BiNGS-Bulk-RNA-seq-pipeline.git
```

3. **Start running from Step A1.**

> 📘 For more on `screen`, see: [linuxize.com/post/how-to-use-linux-screen](https://linuxize.com/post/how-to-use-linux-screen/)

---

## 🔧 Running Scripts

For running all scripts after A1, you don't need to ask as much resources. You can most likely get away with:
   - Cores: `4`
   - Memory: `6 GB`

### Bash Scripts (A1-9.sh)

```bash
bash /path/to/script.sh
```

### R Scripts (B1/B2/C/D.R)

```bash
ml R/4.2.0
R
source("/path/to/script.R")
```

---

## 🧹 Good Practices

- Do **not** modify pipeline scripts.
- Answer prompts directly in the terminal and press Enter.
- Ensure your sample metadata is properly filled in when prompted by `A1.prepare_files.sh`.

> For single-end datasets, leave `file_name_2` blank.

---

## 🛠️ Troubleshooting

### "403 Forbidden" error while running the A1 script

```text
--2025-04-28 10:48:32--  https://ulukag01.u.hpc.mssm.edu/data_rna/raw/sample_metadata/
Connecting to 172.28.7.1:3128... connected.
Proxy request sent, awaiting response... 403 Forbidden
ERROR 403: Forbidden.
```

**Solution:**  
Check that you can access the following in your browser:

- https://ulukag01.u.hpc.mssm.edu/data_rna/
- https://ulukag01.u.hpc.mssm.edu/supporting_files/

If these links do not load, there may be a temporary server issue at Mount Sinai. Please delete your pipeline folder and try again later.

### Issues with LSF job submission

Please check out the Minerva job submission documentation for troubleshooting lsf job issues: 

👉 https://labs.icahn.mssm.edu/minervalab/documentation/lsf-job-scheduler 

---

## 📫 Questions?

If you encounter an issue or have a feature request, please open an **Issue** on the GitHub repository:  
👉 [https://github.com/GulayBenguUlukaya/BiNGS-Bulk-RNA-seq-pipeline/issues](https://github.com/GulayBenguUlukaya/BiNGS-Bulk-RNA-seq-pipeline/issues)

---

## 📚 Citation

Ulukaya, G. B. (2025). *BiNGS-Bulk-RNA-seq-pipeline: An end-to-end bulk RNA-seq analysis pipeline* (Version 1.0) [Software]. GitHub.  
https://github.com/GulayBenguUlukaya/BiNGS-Bulk-RNA-seq-pipeline

> A lot of time and effort went into developing this pipeline. Please cite this material if it contributed to your research.

---

## © License

These materials were developed by members of the [BiNGS Core](https://bings.mssm.edu/) at Mount Sinai.  
They are open access and licensed for unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.
