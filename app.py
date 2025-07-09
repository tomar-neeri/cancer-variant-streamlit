import os
import shutil
import gzip
import zipfile
import subprocess
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from pathlib import Path
import glob

# === Tool Paths ===
makeblastdb_path = "makeblastdb"
blastn_path = "blastn"
bwa_path = "bwa"
samtools_path = "samtools"
bcftools_path = "bcftools"

# === File Constants ===
cancer_fasta = "cancer_exons.fasta"
blast_db_prefix = "cancer_exon_db"

# === Directories ===
os.makedirs("uploads", exist_ok=True)
os.makedirs("converted_fastas", exist_ok=True)
os.makedirs("blast_outputs", exist_ok=True)
os.makedirs("alignments", exist_ok=True)
os.makedirs("variants", exist_ok=True)

# === Streamlit Config ===
st.set_page_config(page_title="🧬 Cancer Variant Pipeline", layout="wide")
st.title("🧬 Cancer Variant Detection from Metagenomic datasets")

# === Clean folders function ===
def clean_folders(folders):
    """
    Deletes all files inside the specified folders.
    Args:
        folders (list): List of folder paths as strings.
    """
    for folder in folders:
        files = glob.glob(os.path.join(folder, "*"))
        for f in files:
            try:
                if os.path.isfile(f) or os.path.islink(f):
                    os.unlink(f)  # remove file or link
                elif os.path.isdir(f):
                    shutil.rmtree(f)  # remove directory recursively
            except Exception as e:
                st.warning(f"Could not delete {f}. Reason: {e}")

# === Create BLAST DB ===
def create_blast_db():
    try:
        subprocess.run([
            makeblastdb_path,
            "-in", cancer_fasta,
            "-dbtype", "nucl",
            "-out", blast_db_prefix
        ], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        st.error(f"❌ makeblastdb failed:\n{e.stderr}")
        raise

# === Ensure BWA Index ===
def ensure_bwa_index(fasta_file):
    index_files = [fasta_file + ext for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]]
    if not all(os.path.exists(f) for f in index_files):
        st.info("🔧 Indexing reference with BWA...")
        try:
            subprocess.run([bwa_path, "index", fasta_file], check=True, capture_output=True, text=True)
            st.success("✅ BWA index created.")
        except subprocess.CalledProcessError as e:
            st.error(f"❌ BWA indexing failed:\n{e.stderr}")
            raise

# === Prepare Sample ===
def prepare_sample(uploaded_file, output_dir="uploads", convert_dir="converted_fastas"):
    raw_path = os.path.join(output_dir, uploaded_file.name)
    with open(raw_path, "wb") as f:
        f.write(uploaded_file.getbuffer())

    decompressed_path = None
    fasta_path = None

    # Handle .fq.gz or .fastq.gz
    if uploaded_file.name.endswith((".fq.gz", ".fastq.gz")):
        decompressed_path = raw_path.replace(".gz", "")
        with gzip.open(raw_path, "rt") as f_in, open(decompressed_path, "w") as f_out:
            shutil.copyfileobj(f_in, f_out)
        fasta_path = os.path.join(convert_dir, Path(decompressed_path).stem + ".fasta")
        SeqIO.write(SeqIO.parse(decompressed_path, "fastq"), fasta_path, "fasta")

    # Handle .fastq or .fq
    elif uploaded_file.name.endswith((".fastq", ".fq")):
        fasta_path = os.path.join(convert_dir, Path(raw_path).stem + ".fasta")
        SeqIO.write(SeqIO.parse(raw_path, "fastq"), fasta_path, "fasta")

    # Handle .zip
    elif uploaded_file.name.endswith(".zip"):
        with zipfile.ZipFile(raw_path, 'r') as zip_ref:
            for file in zip_ref.namelist():
                if file.endswith((".fasta", ".fa", ".fna", ".fastq", ".fq")):
                    zip_ref.extract(file, output_dir)
                    decompressed_path = os.path.join(output_dir, file)
                    break
        if decompressed_path is None:
            raise Exception("❌ No valid FASTA or FASTQ inside ZIP.")
        if decompressed_path.endswith((".fastq", ".fq")):
            fasta_path = os.path.join(convert_dir, Path(decompressed_path).stem + ".fasta")
            SeqIO.write(SeqIO.parse(decompressed_path, "fastq"), fasta_path, "fasta")
        else:
            fasta_path = decompressed_path

    # Handle other .gz (e.g., .fasta.gz)
    elif uploaded_file.name.endswith(".gz"):
        decompressed_path = raw_path.replace(".gz", "")
        with gzip.open(raw_path, "rt") as f_in, open(decompressed_path, "w") as f_out:
            shutil.copyfileobj(f_in, f_out)
        fasta_path = decompressed_path

    else:
        fasta_path = raw_path  # fallback for .fasta, .fa, etc.

    # Sample name derived from uploaded filename (without extension)
    sample_name = uploaded_file.name
    for ext in [".gz", ".zip", ".fastq", ".fq", ".fasta", ".fa", ".fna"]:
        sample_name = sample_name.replace(ext, "")
    sample_name = sample_name.replace("/", "_")

    return sample_name, fasta_path


# === BLAST ===
def run_blast(query_fasta, output_file):
    subprocess.run([
        blastn_path,
        "-query", query_fasta,
        "-db", blast_db_prefix,
        "-out", output_file,
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore"
    ], check=True)

# === Relative Abundance ===
def compute_relative_abundance(blast_file, sample_name):
    cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
    df = pd.read_csv(blast_file, sep="\t", names=cols)
    if df.empty:
        return pd.DataFrame(columns=["Cancer_Gene", sample_name])
    df["Gene"] = df["sseqid"].str.extract(r"^([A-Za-z0-9\-]+)_")
    abundance = df["Gene"].value_counts().rename_axis("Cancer_Gene").reset_index(name="Hit_Count")
    abundance["Relative (%)"] = (abundance["Hit_Count"] / abundance["Hit_Count"].sum()) * 100
    return abundance.set_index("Cancer_Gene")[["Relative (%)"]].rename(columns={"Relative (%)": sample_name})

# === Align & Call Variants ===
def align_and_call_variants(fasta_path, sample_name):
    ensure_bwa_index(cancer_fasta)
    sam_path = f"alignments/{sample_name}.sam"
    bam_path = f"alignments/{sample_name}.bam"
    sorted_bam = f"alignments/{sample_name}_sorted.bam"
    raw_vcf = f"variants/{sample_name}_raw.vcf"

    subprocess.run([bwa_path, "mem", cancer_fasta, fasta_path], stdout=open(sam_path, "w"), check=True)
    subprocess.run([samtools_path, "view", "-Sb", sam_path, "-o", bam_path], check=True)
    subprocess.run([samtools_path, "sort", "-o", sorted_bam, bam_path], check=True)
    subprocess.run([samtools_path, "index", sorted_bam], check=True)

    mpileup = subprocess.Popen([
        bcftools_path, "mpileup", "-f", cancer_fasta, sorted_bam, "-Ou"
    ], stdout=subprocess.PIPE)

    subprocess.run([
        bcftools_path, "call", "-mv", "-Ov", "-o", raw_vcf
    ], stdin=mpileup.stdout, check=True)

    return raw_vcf

# === Parse VCF ===
def parse_variants(vcf_file, sample_name):
    variants = []
    with open(vcf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]
            gene = chrom.split("_")[0]
            variants.append((gene, f"{ref}>{alt}", pos))
    df = pd.DataFrame(variants, columns=["Cancer_Gene", f"{sample_name}_Variant", "Position"])
    return df

# === Build Mutation Count Matrix ===
def build_mutation_matrix(variant_dir):
    matrix = {}
    for file in os.listdir(variant_dir):
        if file.endswith("_variants.csv"):
            sample_name = file.replace("_variants.csv", "")
            df = pd.read_csv(os.path.join(variant_dir, file))
            counts = df["Cancer_Gene"].value_counts()
            for gene, count in counts.items():
                if gene not in matrix:
                    matrix[gene] = {}
                matrix[gene][sample_name] = count
    return pd.DataFrame(matrix).T.fillna(0).astype(int)

# === Visualizations ===
def visualize_abundance(abundance_df):
    st.subheader("🧪 Relative Abundance Heatmap")
    plt.figure(figsize=(12, 6))
    sns.heatmap(abundance_df, cmap="magma", annot=True, fmt=".1f")
    plt.title("Relative Abundance of Cancer Genes")
    st.pyplot(plt)

def visualize_mutations(mutation_df):
    st.subheader("🔬 Mutation Count Heatmap")
    plt.figure(figsize=(12, 6))
    sns.heatmap(mutation_df, cmap="viridis", annot=True, fmt="d")
    plt.title("Mutation Count per Gene")
    st.pyplot(plt)

# === Upload Interface ===
with st.expander("📁 Upload FASTA / FASTQ / .gz / .zip (≤ 2 GB each)"):
    uploaded_files = st.file_uploader(
        "Upload one or more sequence files",
        type=["fasta", "fa", "fna", "fastq", "fq", "gz", "zip"],
        accept_multiple_files=True
    )

# === Pipeline Runner ===
if uploaded_files:
    # Clean old data first before processing new files
    clean_folders([
        "uploads",
        "converted_fastas",
        "blast_outputs",
        "alignments",
        "variants"
    ])

    create_blast_db()
    all_abundances = []

    for uploaded_file in uploaded_files:
        try:
            sample_name, fasta_path = prepare_sample(uploaded_file)
            st.info(f"🔬 Processing `{sample_name}`")

            blast_out = f"blast_outputs/{sample_name}_blast.csv"
            variant_csv = f"variants/{sample_name}_variants.csv"

            run_blast(fasta_path, blast_out)
            abundance = compute_relative_abundance(blast_out, sample_name)
            all_abundances.append(abundance)

            vcf = align_and_call_variants(fasta_path, sample_name)
            var_df = parse_variants(vcf, sample_name)
            var_df.to_csv(variant_csv, index=False)

            st.success(f"✅ `{sample_name}` complete.")

        except Exception as e:
            st.error(f"❌ Error processing {uploaded_file.name}: {e}")

    # Visualize
    if all_abundances:
        combined_abundance = pd.concat(all_abundances, axis=1).fillna(0)
        st.download_button("⬇ Download Abundance CSV", combined_abundance.to_csv().encode(), file_name="abundance_summary.csv")
        visualize_abundance(combined_abundance)

        mutation_df = build_mutation_matrix("variants")
        st.download_button("⬇ Download Mutation CSV", mutation_df.to_csv().encode(), file_name="mutation_summary.csv")
        visualize_mutations(mutation_df)
else:
    st.info("📂 Upload files to begin analysis.")
