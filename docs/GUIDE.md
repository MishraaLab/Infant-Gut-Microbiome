# Long-Read Microbiome Analysis Guide

## Introduction

This guide provides detailed information about analyzing long-read microbiome data from PacBio sequencing platforms. The scripts are specifically designed for infant gut microbiome studies but can be adapted for other microbiome research.

## Long-Read Advantages

PacBio long-read sequencing offers several advantages for microbiome studies:

1. **Full-length 16S rRNA genes**: Capture complete 16S sequences (~1.5 kb) for better taxonomic resolution
2. **Improved species-level identification**: Longer reads provide more sequence context
3. **Reduced amplification bias**: HiFi circular consensus sequencing provides high accuracy
4. **Better strain discrimination**: Longer reads can distinguish closely related strains

## Workflow Overview

### 1. Quality Control

**Input**: Raw PacBio reads in FASTQ format

**Process**:
- Assess read quality using NanoPlot
- Filter reads by length (recommended: ≥1000 bp for 16S, ≥500 bp for shotgun)
- Filter reads by quality score (recommended: Q≥10)

**Output**: Filtered, high-quality reads

**Key Parameters**:
```bash
--min-length 1000    # Minimum read length in base pairs
--min-quality 10     # Minimum average Phred quality score
```

### 2. Taxonomic Classification

**Input**: Quality-filtered reads

**Methods**:

#### Kraken2 (Recommended for shotgun metagenomics)
- Fast k-mer based classification
- Requires pre-built database (standard, custom, or specialized)
- Can classify to species level

#### Minimap2 (Recommended for 16S amplicon)
- Alignment-based classification
- Use with SILVA or Greengenes database
- More sensitive for 16S sequences

**Output**: Taxonomic assignments and abundance reports

**Key Parameters**:
```bash
--database /path/to/db    # Path to classification database
--method kraken2          # Classification method
--threads 8               # Number of CPU threads
```

### 3. Diversity Analysis

**Input**: Taxonomic abundance tables

**Alpha Diversity Metrics**:
- **Richness**: Number of observed taxa
- **Shannon Index**: Accounts for abundance and evenness
- **Simpson Index**: Probability two random individuals are different species
- **Pielou's Evenness**: How evenly distributed the abundances are

**Beta Diversity Metrics**:
- **Bray-Curtis**: Quantitative measure based on abundance
- **Jaccard**: Qualitative measure based on presence/absence

**Output**: Diversity metrics and distance matrices

### 4. Statistical Analysis

**Input**: Abundance tables + metadata

**Analyses**:
- **Differential abundance**: Identify taxa with significant differences between groups
- **Alpha diversity comparison**: Compare diversity metrics between groups
- **Beta diversity comparison**: Test for differences in community composition

**Statistical Tests**:
- Mann-Whitney U test (2 groups)
- Kruskal-Wallis test (>2 groups)
- PERMANOVA-like analysis for beta diversity
- FDR correction for multiple testing

**Key Parameters**:
```bash
--group treatment         # Metadata column for grouping
--normalize total_sum     # Normalization method
--min-count 10           # Filter taxa with <10 total reads
--min-prevalence 0.1     # Filter taxa present in <10% of samples
```

### 5. Visualization

**Input**: Analysis results

**Plots**:
- **Taxonomic barplots**: Stacked bar charts showing composition
- **Alpha diversity boxplots**: Compare diversity across groups
- **PCoA plots**: Visualize beta diversity patterns
- **Heatmaps**: Hierarchical clustering of abundance data
- **Volcano plots**: Differential abundance results
- **Rarefaction curves**: Assess sampling depth adequacy

## Best Practices

### Sample Preparation
- Maintain consistent DNA extraction protocols
- Include negative controls
- Record metadata thoroughly (treatment, timepoint, etc.)

### Quality Control
- Inspect raw data quality before filtering
- Compare before/after filtering statistics
- Aim for >80% read retention after QC

### Taxonomic Classification
- Use appropriate database for your study:
  - Standard database for general microbiome
  - RefSeq database for higher specificity
  - Custom database for targeted studies
- Consider running multiple classifiers and comparing

### Diversity Analysis
- Check rarefaction curves reach plateau
- Normalize data appropriately for comparisons
- Use appropriate statistical tests (parametric vs non-parametric)

### Statistical Analysis
- Pre-filter low abundance taxa
- Apply appropriate normalization
- Use FDR correction for multiple testing
- Validate results with biological knowledge

## Database Recommendations

### For Shotgun Metagenomics (Kraken2)

**Standard Database** (~50 GB):
- Contains bacteria, archaea, viruses
- Good for general microbiome studies
```bash
kraken2-build --standard --threads 8 --db kraken2_standard_db
```

**RefSeq Database** (~100 GB):
- More comprehensive reference sequences
- Better for species-level identification

**Custom Database**:
- Add specific genomes of interest
- Useful for specialized studies

### For 16S Amplicon (Minimap2)

**SILVA Database**:
- Most comprehensive 16S/18S database
- Regular updates
- Download: https://www.arb-silva.de/

**Greengenes Database**:
- Well-curated bacterial/archaeal database
- Compatible with many tools
- Download: https://greengenes.secondgenome.com/

## Troubleshooting

### Low Classification Rate
- Check database completeness
- Verify input data quality
- Consider using more comprehensive database

### Memory Issues
- Reduce number of threads
- Use smaller database
- Process samples in batches

### Inconsistent Results
- Check for batch effects
- Verify metadata accuracy
- Ensure consistent preprocessing

### Missing Dependencies
- Install via conda: `conda install -c bioconda <tool>`
- Check PATH and environment variables
- Verify tool versions

## Performance Optimization

### Speed
- Use more threads for parallel processing
- Use SSD for database storage
- Pre-index databases

### Memory
- Kraken2 standard database requires ~40 GB RAM
- Reduce database size if needed
- Use `--memory-mapping` flag for Kraken2

### Storage
- Compress intermediate files
- Delete temporary files after analysis
- Use symbolic links for large databases

## Advanced Usage

### Batch Processing
```bash
# Process multiple samples
for sample in data/raw_data/*.fastq; do
    python scripts/01_quality_control.py -i $sample -o data/processed_data/$(basename $sample)
done
```

### Custom Analysis
- Modify scripts for specific needs
- Add custom visualization
- Integrate with other tools

### Integration with Other Pipelines
- Export to BIOM format for QIIME2
- Export to Phyloseq for R analysis
- Export to MicrobiomeAnalyst

## Resources

### Tools Documentation
- Kraken2: https://github.com/DerrickWood/kraken2/wiki
- Bracken: https://ccb.jhu.edu/software/bracken/
- NanoPlot: https://github.com/wdecoster/NanoPlot
- Minimap2: https://github.com/lh3/minimap2

### Tutorials
- PacBio microbiome analysis: https://github.com/PacificBiosciences/pb-metagenomics-tools
- Kraken2 tutorial: https://ccb.jhu.edu/software/kraken2/
- Microbiome data analysis: https://microbiomejournal.biomedcentral.com/

### Publications
- Callahan et al. (2019) "High-resolution sample inference from Illumina amplicon data"
- Sczyrba et al. (2017) "Critical Assessment of Metagenome Interpretation"
- Wagner et al. (2016) "Survey of the human gut microbiome"

## Support

For questions, issues, or suggestions:
1. Check this documentation
2. Search existing GitHub issues
3. Open a new issue with detailed description
4. Contact the Mishraa Lab

## Future Enhancements

Planned features:
- [ ] Support for Nanopore long reads
- [ ] Integration with functional profiling tools
- [ ] Machine learning classification
- [ ] Interactive visualization dashboard
- [ ] Automated report generation
- [ ] Docker containerization
- [ ] Cloud computing support
