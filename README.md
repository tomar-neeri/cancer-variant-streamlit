# 🧬 Cancer Variant Detection from Metagenomic Datasets

This Streamlit web application detects cancer-associated exon variants from uploaded metagenomic sequencing data using BLAST, BWA, Samtools, and BCFtools. It performs alignment, variant calling, and visualization for relative gene abundance and mutation counts.

---

## 📁 Features

- Upload `.fasta`, `.fastq`, `.fq`, `.gz`, or `.zip` files.
- Detect hits using **BLASTn** against a cancer exon reference.
- Perform alignments with **BWA** and variant calling using **Samtools + BCFtools**.
- Generate:
  - Relative abundance heatmaps.
  - Mutation count heatmaps per cancer gene.
  - Downloadable CSV summaries.

---

## 🖥️ Run Locally (WSL/Linux)

### 🔧 Prerequisites

- Python 3.8+
- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- BWA
- Samtools
- BCFtools
- gzip, zip
- Git

### 📦 System Package Installation (Ubuntu/WSL)

```bash
sudo apt update && sudo apt install -y \
    bwa samtools bcftools ncbi-blast+ \
    gzip unzip build-essential python3-venv
```

---

### 🚀 Steps to Run

```bash
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
```

Visit: [http://localhost:8501](http://localhost:8501)

---

## 🧪 Input File Formats

You can upload multiple files in any of the following formats:

- `.fasta`, `.fa`, `.fna`
- `.fastq`, `.fq`
- `.gz` (compressed)
- `.zip` (with valid sequences inside)

Each sample will be processed to compute:

- BLAST hits to known cancer exons
- BWA alignment and variant calling
- Relative gene hit counts and mutations

---

## 📊 Outputs

- **Relative Abundance Heatmap**  
  Cancer gene-wise relative % abundance from BLAST hits.

- **Mutation Count Heatmap**  
  Number of mutations per cancer gene from alignment-based variant calling.

- **Downloadable Files**
  - `abundance_summary.csv` — Gene-wise relative abundance
  - `mutation_summary.csv` — Gene-wise mutation counts

---

## 🧠 Internal Workflow

1. **Upload and format conversion**  
   All uploaded files are decompressed and converted to `.fasta` if needed using Biopython.

2. **BLASTn Search**  
   Against the included `cancer_exons.fasta` database, reporting gene hits.

3. **BWA Alignment + Variant Calling**  
   Variant detection from aligned reads using Samtools & BCFtools.

4. **Visualization**  
   Using Seaborn to generate annotated heatmaps.

---

## ⚠️ Notes

- **Do not upload real patient data.** This is for research/demo purposes.
- App uses local command-line tools (`blastn`, `bwa`, etc.). Ensure paths are correct if modifying.

---

## 🧬 File Tree (Key Files)

```
.
├── app.py                      # Main Streamlit application
├── cancer_exons.fasta          # Reference database of cancer exons
├── requirements.txt            # Python dependencies
├── Dockerfile                  # For Docker deployment (optional)
├── .dockerignore               # Excludes keys and data from Docker context
├── /uploads/                   # Uploaded sample files
├── /converted_fastas/          # Converted .fasta files
├── /blast_outputs/             # BLAST result files
├── /alignments/                # SAM/BAM/VCF intermediate files
├── /variants/                  # Parsed CSV variant results
└── .streamlit/config.toml      # Optional Streamlit UI settings
```

---

## 🐳 Docker Support

```bash
docker build -t cancer-pipeline .
docker run -p 8501:8501 cancer-pipeline
```

> ⚠️ Tools like `bwa` and `samtools` must be available in the Docker image.

---

## 🌐 Deployment Notes

- **Streamlit Cloud** doesn’t support `bwa`, `samtools`, or `blastn`. Use local or GCP/AWS for full functionality.
- For production deployment, consider:
  - GCP with Compute Engine or App Engine
  - AWS EC2 or Lambda + EFS
  - Azure Container Instances

---

## 👨‍🔬 Author

**Siddharth Singh Tomar**  
CSIR-NEERI  
📧 siddharthsinghtomar166@gmail.com

---

## 📜 License

This project is licensed under the [MIT License](LICENSE).

---
