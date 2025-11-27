# -------- 1. Conda Environment Setup --------
conda create -n ngs_analysis python=3.10 r-base=4.2.0
conda activate ngs_analysis
conda install -c bioconda fastqc multiqc trim-galore cutadapt samtools subread

# Check software versions
fastqc --version
multiqc --version
trim_galore --version
cutadapt --version

# -------- Optional FastQC Setup --------
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
cd FastQC/
chmod 755 fastqc

# -------- MultiQC Setup --------
conda install -c bioconda multiqc
multiqc --version

conda create -n multiqc_env python=3.10
conda activate multiqc_env
conda install -c bioconda multiqc

# -------- Trim Galore Installation --------
wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.10.zip
unzip 0.6.10.zip
cd TrimGalore-0.6.10
chmod +x trim_galore
sudo ln -s $(pwd)/trim_galore /usr/local/bin/
pip install cutadapt
cutadapt --version

conda install -c bioconda trim-galore

# -------- 2. Data Processing --------
pwd
echo $'\n'; echo "======== listing files ========"; ls; echo $'\n'

# -------- Download Reference Genome --------
mkdir -p ~/Suriya/Align/Ref_Genome
cd ~/Suriya/Align/Ref_Genome

wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna \
     /Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

ls -lh *.gz

# -------- Decompress Reference Files --------
gunzip -c Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz > Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
gunzip -c Homo_sapiens.GRCh38.115.gtf.gz > Homo_sapiens.GRCh38.115.gtf
ls -lh Homo_sapiens.GRCh38*

# -------- 3. Quality Control Commands --------
cd /home/vamouda/Suriya/Data/

fastqc *.fastq.gz -o ../FastQC_Results/

multiqc ../FastQC_Results/ -o ../MultiQC_Report/

# -------- 4. Trimming Reads --------
cd /home/vamouda/Suriya/Data/

find /home/vamouda/Suriya/Data/ -name "*_R1.fastq.gz" | while read r1; do
    r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
    if [[ -f "$r2" ]]; then
        echo "Processing: $(basename $r1) and $(basename $r2)"
        trim_galore --paired --quality 30 --fastqc "$r1" "$r2"
    else
        echo "Warning: Missing paired file for $r1"
    fi
done
