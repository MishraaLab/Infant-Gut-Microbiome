# Data Directory

This directory contains data files for the microbiome analysis pipeline.

## Structure

- `raw_data/`: Place your raw sequencing data (FASTQ files) here
- `processed_data/`: Filtered and processed sequences will be saved here
- `results/`: Analysis results, including classification, diversity, statistics, and visualizations

## Data Not Included

Due to file size limitations, sequencing data files are not included in this repository. 

### Expected Data Format

**Raw Data** (place in `raw_data/`):
- PacBio HiFi reads in FASTQ format
- File naming: `sample_name.fastq` or `sample_name.fastq.gz`
- Quality scores in Phred+33 format

**Metadata** (optional, place in this directory):
- Tab-separated values (TSV) format
- Required columns: `sample` (sample ID)
- Additional columns: group information (e.g., `treatment`, `timepoint`, `age_days`)
- See `metadata.example.tsv` for format

### Example Data Download

For testing the pipeline, you can use publicly available datasets:

1. **PacBio Sequel II Gut Microbiome Data**
   - Available at: NCBI SRA, EBI ENA
   - Search for: "PacBio HiFi microbiome"

2. **Mock Community Data**
   - ZymoBIOMICS Microbial Community Standards
   - Available at: https://www.zymoresearch.com/

### Adding Your Data

```bash
# Copy your FASTQ files
cp /path/to/your/samples/*.fastq raw_data/

# Create metadata file
nano metadata.tsv
```

Then run the pipeline as described in the main README.
