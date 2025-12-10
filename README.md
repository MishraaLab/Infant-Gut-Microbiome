# Infant-Gut-Microbiome

Preterm infants are prone to gut microbiota dysbiosis, increasing risks of NEC and sepsis. Probiotics may restore microbial balance. This study uses PacBio long-read sequencing to assess probiotic effects in preterm infants, with full-term infants as the reference.

## Overview

This repository contains scripts and codes for analyzing long-read microbiome data from infant gut samples. The pipeline supports PacBio HiFi reads and includes comprehensive tools for quality control, taxonomic classification, diversity analysis, statistical comparisons, and visualization.

## Repository Structure

```
Infant-Gut-Microbiome/
├── scripts/                          # Analysis scripts
│   ├── 01_quality_control.py         # QC and filtering for long reads
│   ├── 02_taxonomic_classification.py # Taxonomic classification (Kraken2/Minimap2)
│   ├── 03_diversity_analysis.py      # Alpha/beta diversity calculations
│   ├── 04_statistical_analysis.py    # Differential abundance and statistics
│   ├── 05_visualization.py           # Generate plots and figures
│   └── pipeline.py                   # Main integrated pipeline
├── data/                             # Data directory
│   ├── raw_data/                     # Raw sequencing data
│   ├── processed_data/               # Processed/filtered data
│   └── results/                      # Analysis results
├── docs/                             # Documentation
├── config.example.ini                # Example configuration file
├── requirements.txt                  # Python dependencies
└── README.md                         # This file
```

## Features

### 1. Quality Control (`01_quality_control.py`)
- Quality assessment using NanoPlot
- Length-based filtering
- Quality score filtering
- Summary statistics (N50, read length distribution, quality scores)

### 2. Taxonomic Classification (`02_taxonomic_classification.py`)
- Kraken2-based classification for metagenomic data
- Minimap2-based classification for 16S rRNA data
- Bracken abundance re-estimation
- BIOM format export

### 3. Diversity Analysis (`03_diversity_analysis.py`)
- Alpha diversity metrics:
  - Species richness
  - Shannon diversity
  - Simpson diversity
  - Pielou's evenness
- Beta diversity metrics:
  - Bray-Curtis dissimilarity
  - Jaccard distance
- Rarefaction curves

### 4. Statistical Analysis (`04_statistical_analysis.py`)
- Differential abundance analysis
- Alpha diversity comparisons (Mann-Whitney U, Kruskal-Wallis)
- Beta diversity comparisons (PERMANOVA-like analysis)
- Multiple testing correction (Benjamini-Hochberg)
- Data normalization (TSS, CSS, CLR)

### 5. Visualization (`05_visualization.py`)
- Taxonomic composition barplots
- Alpha diversity boxplots
- Beta diversity PCoA plots
- Abundance heatmaps with hierarchical clustering
- Differential abundance volcano plots
- Rarefaction curves

### 6. Integrated Pipeline (`pipeline.py`)
- Runs all analysis steps in sequence
- Configurable via INI file
- Automatic logging
- Error handling

## Installation

### Requirements

- Python 3.7+
- External tools (optional, depending on analysis):
  - [Kraken2](https://github.com/DerrickWood/kraken2) - taxonomic classification
  - [Bracken](https://github.com/jenniferlu717/Bracken) - abundance estimation
  - [Minimap2](https://github.com/lh3/minimap2) - alignment-based classification
  - [NanoPlot](https://github.com/wdecoster/NanoPlot) - QC visualization

### Install Python Dependencies

```bash
# Create a virtual environment (recommended)
python3 -m venv microbiome_env
source microbiome_env/bin/activate  # On Windows: microbiome_env\Scripts\activate

# Install Python packages
pip install -r requirements.txt
```

### Install External Tools (via Conda - recommended)

```bash
# Create conda environment with all tools
conda create -n microbiome python=3.9
conda activate microbiome

# Install tools
conda install -c bioconda kraken2 bracken minimap2 nanoplot
conda install -c conda-forge numpy pandas scipy matplotlib seaborn biopython statsmodels scikit-learn
```

## Quick Start

### 1. Prepare Your Data

Place your long-read FASTQ files in the `data/raw_data/` directory:

```bash
cp your_sample.fastq data/raw_data/
```

### 2. Create Configuration File

Copy and edit the example configuration:

```bash
cp config.example.ini config.ini
# Edit config.ini with your paths and parameters
```

### 3. Run the Pipeline

```bash
python scripts/pipeline.py -c config.ini
```

Or run individual steps:

```bash
# Quality control
python scripts/01_quality_control.py -i data/raw_data/sample.fastq -o data/processed_data/filtered.fastq

# Taxonomic classification
python scripts/02_taxonomic_classification.py -i data/processed_data/filtered.fastq -o data/results/classification -d /path/to/kraken2/db

# Diversity analysis
python scripts/03_diversity_analysis.py -i data/results/classification/kraken2_report.txt -o data/results/diversity

# Visualization
python scripts/05_visualization.py -o data/results/visualization --alpha-diversity data/results/diversity/alpha_diversity.tsv
```

## Usage Examples

### Example 1: Basic QC and Classification

```bash
# Run quality control
python scripts/01_quality_control.py \
    -i data/raw_data/sample.fastq \
    -o data/processed_data/filtered.fastq \
    --min-length 1000 \
    --min-quality 10

# Run taxonomic classification
python scripts/02_taxonomic_classification.py \
    -i data/processed_data/filtered.fastq \
    -o data/results/classification \
    -d /path/to/kraken2/standard_db \
    --threads 8
```

### Example 2: Diversity Analysis

```bash
python scripts/03_diversity_analysis.py \
    -i data/results/classification/kraken2_report.txt \
    -o data/results/diversity \
    --format kraken \
    --alpha \
    --rarefaction
```

### Example 3: Statistical Comparison with Metadata

```bash
# Create metadata file (see data/metadata.example.tsv)
# Then run statistical analysis

python scripts/04_statistical_analysis.py \
    -a data/results/abundance_table.tsv \
    -m data/metadata.tsv \
    -o data/results/statistics \
    -g treatment \
    --alpha-diversity data/results/diversity/alpha_diversity.tsv
```

### Example 4: Generate Visualizations

```bash
python scripts/05_visualization.py \
    -a data/results/abundance_table.tsv \
    -m data/metadata.tsv \
    -o data/results/plots \
    --alpha-diversity data/results/diversity/alpha_diversity.tsv \
    --differential-abundance data/results/statistics/differential_abundance.tsv \
    --group treatment
```

## Configuration

The pipeline uses an INI configuration file. Key parameters:

### Paths
- `input_fastq`: Input FASTQ file with long reads
- `output_dir`: Output directory for results
- `abundance_table`: Optional abundance table for statistics
- `metadata`: Optional metadata file for group comparisons

### Pipeline Steps
Enable/disable individual steps:
- `run_qc`: Quality control
- `run_classification`: Taxonomic classification
- `run_diversity`: Diversity analysis
- `run_statistics`: Statistical analysis
- `run_visualization`: Generate plots

### QC Parameters
- `min_length`: Minimum read length (default: 1000 bp)
- `min_quality`: Minimum average quality score (default: 10)

### Classification Parameters
- `database`: Path to Kraken2 database
- `method`: Classification method (kraken2 or minimap2)
- `threads`: Number of CPU threads

See `config.example.ini` for full options.

## Output Structure

```
data/results/
├── qc_results/              # Quality control results
│   ├── filtered_reads.fastq # Filtered reads
│   ├── qc_summary.txt       # QC statistics
│   ├── raw/                 # NanoPlot results (raw data)
│   └── filtered/            # NanoPlot results (filtered data)
├── classification/          # Taxonomic classification
│   ├── kraken2_output.txt   # Per-read classifications
│   ├── kraken2_report.txt   # Aggregated report
│   └── taxonomy.biom        # BIOM format export
├── diversity/               # Diversity analysis
│   ├── alpha_diversity.tsv  # Alpha diversity metrics
│   ├── beta_diversity_*.tsv # Distance matrices
│   └── rarefaction_curves.tsv # Rarefaction data
├── statistics/              # Statistical analysis
│   ├── differential_abundance.tsv
│   ├── alpha_diversity_comparison.tsv
│   └── beta_diversity_comparison.tsv
└── visualization/           # Plots and figures
    ├── taxonomic_barplot.png
    ├── alpha_diversity.png
    ├── beta_diversity_pcoa.png
    ├── abundance_heatmap.png
    ├── differential_abundance.png
    └── rarefaction_curves.png
```

## Database Setup

### Kraken2 Standard Database

```bash
# Download and build standard database (~50 GB)
kraken2-build --standard --threads 8 --db kraken2_standard_db

# Or use pre-built databases
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20210517.tar.gz
tar -xzf k2_standard_20210517.tar.gz
```

### Custom 16S rRNA Database for Minimap2

```bash
# Download SILVA or Greengenes database
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_NR99_tax_silva.fasta.gz
gunzip SILVA_138_SSURef_NR99_tax_silva.fasta.gz

# Index with minimap2
minimap2 -d silva_138.mmi SILVA_138_SSURef_NR99_tax_silva.fasta
```

## Citation

If you use these scripts in your research, please cite:

[Publication details to be added]

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

[License to be added]

## Contact

For questions or issues, please open an issue on GitHub or contact the Mishraa Lab.

## Acknowledgments

This work focuses on understanding the gut microbiome of preterm and full-term infants, with particular emphasis on probiotic effects on microbial community composition and diversity.
