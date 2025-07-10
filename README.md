🧬 Cancer Variant Detection from Metagenomic Datasets
This Streamlit-based web application detects cancer-associated exon variants from uploaded metagenomic sequencing data. The pipeline uses BLAST, BWA, SAMtools, and BCFtools to compute relative abundance and perform variant calling, generating both tabular and heatmap visualizations.

📁 Features
Upload .fasta, .fastq, .fq, .gz, or .zip files.

Detect cancer exon hits using BLASTn.

Perform alignments with BWA and call variants with SAMtools + BCFtools.

Automatically decompress and convert formats via Biopython.

Generates:

🧪 Relative abundance heatmap per gene.

🔬 Mutation count heatmap per gene.

📄 Downloadable CSV summaries.

🖼 PNG downloads of heatmaps.

📦 One-click download of all outputs in a .zip archive.

🖥️ Run Locally (Linux/WSL Recommended)
🔧 Prerequisites
Python 3.8+

makeblastdb, blastn (from BLAST+)

bwa

samtools

bcftools

gzip, zip, unzip

git

📦 System Package Installation (Ubuntu/WSL)
bash
Copy
Edit
sudo apt update && sudo apt install -y \
  bwa samtools bcftools ncbi-blast+ \
  gzip unzip build-essential python3-venv
🚀 Steps to Run
bash
Copy
Edit
# Clone the repository
git clone https://github.com/tomar-neeri/cancer-variant-streamlit.git
cd cancer-variant-streamlit

# (Optional) Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install Python packages
pip install --upgrade pip
pip install -r requirements.txt

# Run the Streamlit app
streamlit run app.py
📍 Visit: http://localhost:8501

🧪 Input File Formats
You can upload multiple files in any of the following formats:

.fasta, .fa, .fna

.fastq, .fq

.gz (compressed FASTQ/FASTA)

.zip (must contain valid sequence files)

For each sample:
Performs BLASTn search against cancer exons.

Aligns reads with BWA and calls variants.

Extracts gene-wise abundance and mutation details.

📊 Outputs
📈 Relative Abundance Heatmap
Gene-wise percent abundance from BLAST hits.

🧬 Mutation Count Heatmap
Number of unique mutations detected per cancer gene.

📁 Downloadable Outputs
File	Description
abundance_summary.csv	Relative abundance per cancer gene
mutation_summary.csv	Mutation counts per cancer gene
abundance_heatmap.png	Heatmap image of relative abundance
mutation_heatmap.png	Heatmap image of mutation counts
outputs_bundle.zip	Contains all of the above + intermediate files

🔬 Internal Workflow
Upload & Decompression
Files are decompressed and converted into .fasta if necessary using Biopython.

BLASTn Search
Queries each sample against cancer_exons.fasta using blastn.

BWA Alignment & Variant Calling
Aligns reads and detects variants with SAMtools and BCFtools.

Visualization
Uses Seaborn to generate annotated heatmaps (PNG & Streamlit embedded).

⚠️ Limitations
❌ Not compatible with Streamlit Cloud
Streamlit Cloud does not support running local binaries like bwa, samtools, or blastn.
✅ Run locally or deploy via GCP, AWS, or Docker for full functionality.

📦 The pipeline assumes pre-formatted exon reference (cancer_exons.fasta) is valid and indexed correctly.

📂 Upload limit set to 2 GB per file.

🔐 Do not upload patient-identifiable data. Intended for research/demonstration purposes only.

📁 File Tree (Important Files)
bash
Copy
Edit
.
├── app.py                      # Main Streamlit application
├── cancer_exons.fasta          # Reference exon database
├── requirements.txt            # Python dependencies
├── Dockerfile                  # (Optional) Docker setup
├── /uploads/                   # Uploaded user files
├── /converted_fastas/          # Converted FASTA files
├── /blast_outputs/             # BLAST result files
├── /alignments/                # SAM/BAM/VCF alignments
├── /variants/                  # Parsed variant files
├── abundance_heatmap.png       # Visualization output
├── mutation_heatmap.png        # Visualization output
└── outputs_bundle.zip          # All output files bundled
🐳 Docker Support (Optional)
bash
Copy
Edit
docker build -t cancer-pipeline .
docker run -p 8501:8501 cancer-pipeline
Note: You must ensure all bioinformatics tools are installed inside the Docker container.

🌐 Deployment Recommendations
For cloud deployment, use:

Google Cloud Platform (GCP) — Compute Engine or App Engine

Amazon AWS — EC2 or ECS (not Lambda)

Microsoft Azure — Container Instances

👨‍🔬 Author
Siddharth Singh Tomar
📧 siddharthsinghtomar166@gmail.com

