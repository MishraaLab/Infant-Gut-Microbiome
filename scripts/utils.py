#!/usr/bin/env python3
"""
Utility functions for microbiome data analysis
Helper functions that can be imported by other scripts
"""

import os
import gzip
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def open_file(file_path, mode='r'):
    """
    Open file, automatically handling gzip compression
    
    Args:
        file_path: Path to file
        mode: File mode ('r', 'w', etc.)
    
    Returns:
        File handle
    """
    if file_path.endswith('.gz'):
        return gzip.open(file_path, mode + 't')
    else:
        return open(file_path, mode)


def count_sequences(fastq_file):
    """
    Count number of sequences in FASTQ file
    
    Args:
        fastq_file: Path to FASTQ file
    
    Returns:
        Number of sequences
    """
    count = 0
    with open_file(fastq_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                count += 1
    return count


def convert_kraken_to_table(kraken_report, output_table):
    """
    Convert Kraken report to simple abundance table
    
    Args:
        kraken_report: Kraken2 report file
        output_table: Output TSV file
    """
    import pandas as pd
    
    data = []
    with open(kraken_report, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                percentage = float(parts[0])
                count = int(parts[1])
                rank = parts[3]
                taxid = parts[4]
                taxon = parts[5].strip()
                
                # Only include species level or higher
                if rank in ['S', 'G', 'F', 'O', 'C', 'P', 'D']:
                    data.append({
                        'taxon': taxon,
                        'rank': rank,
                        'taxid': taxid,
                        'count': count,
                        'percentage': percentage
                    })
    
    df = pd.DataFrame(data)
    df.to_csv(output_table, sep='\t', index=False)
    logger.info(f"Converted Kraken report to table: {output_table}")


def merge_abundance_tables(table_files, output_file):
    """
    Merge multiple abundance tables from different samples
    
    Args:
        table_files: List of abundance table files
        output_file: Output merged table
    """
    import pandas as pd
    
    merged = None
    
    for table_file in table_files:
        sample_name = Path(table_file).stem
        df = pd.read_csv(table_file, sep='\t', index_col=0)
        df.columns = [sample_name]
        
        if merged is None:
            merged = df
        else:
            merged = merged.join(df, how='outer')
    
    # Fill missing values with 0
    merged = merged.fillna(0)
    
    merged.to_csv(output_file, sep='\t')
    logger.info(f"Merged {len(table_files)} tables into {output_file}")


def calculate_read_stats(fastq_file):
    """
    Calculate basic statistics for FASTQ file
    
    Args:
        fastq_file: Path to FASTQ file
    
    Returns:
        Dictionary with statistics
    """
    from Bio import SeqIO
    import numpy as np
    
    lengths = []
    qualities = []
    
    for record in SeqIO.parse(fastq_file, 'fastq'):
        lengths.append(len(record.seq))
        if hasattr(record, 'letter_annotations') and 'phred_quality' in record.letter_annotations:
            qualities.append(np.mean(record.letter_annotations["phred_quality"]))
    
    stats = {
        'num_sequences': len(lengths),
        'total_bases': sum(lengths),
        'mean_length': np.mean(lengths),
        'median_length': np.median(lengths),
        'min_length': min(lengths),
        'max_length': max(lengths),
    }
    
    if qualities:
        stats.update({
            'mean_quality': np.mean(qualities),
            'median_quality': np.median(qualities),
            'min_quality': min(qualities),
            'max_quality': max(qualities)
        })
    
    return stats


def validate_inputs(input_file, file_type='fastq'):
    """
    Validate input file exists and has expected format
    
    Args:
        input_file: Path to input file
        file_type: Expected file type
    
    Returns:
        True if valid, raises exception otherwise
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    if os.path.getsize(input_file) == 0:
        raise ValueError(f"Input file is empty: {input_file}")
    
    # Check file extension
    valid_extensions = {
        'fastq': ['.fastq', '.fq', '.fastq.gz', '.fq.gz'],
        'fasta': ['.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz'],
        'tsv': ['.tsv', '.txt'],
        'csv': ['.csv']
    }
    
    if file_type in valid_extensions:
        if not any(input_file.endswith(ext) for ext in valid_extensions[file_type]):
            logger.warning(f"Input file may not be {file_type} format: {input_file}")
    
    return True


def create_sample_sheet(data_dir, output_file, pattern='*.fastq'):
    """
    Create sample sheet from files in directory
    
    Args:
        data_dir: Directory containing data files
        output_file: Output sample sheet file
        pattern: File pattern to match
    """
    import pandas as pd
    from pathlib import Path
    
    files = list(Path(data_dir).glob(pattern))
    
    samples = []
    for file in files:
        sample_name = file.stem.replace('.fastq', '').replace('.fq', '')
        samples.append({
            'sample': sample_name,
            'file_path': str(file),
            'file_size': file.stat().st_size
        })
    
    df = pd.DataFrame(samples)
    df.to_csv(output_file, sep='\t', index=False)
    
    logger.info(f"Created sample sheet with {len(samples)} samples: {output_file}")


def summarize_classification(kraken_output, kraken_report):
    """
    Generate summary statistics from Kraken classification
    
    Args:
        kraken_output: Kraken2 output file
        kraken_report: Kraken2 report file
    
    Returns:
        Dictionary with summary statistics
    """
    # Count classified/unclassified
    classified = 0
    unclassified = 0
    
    with open(kraken_output, 'r') as f:
        for line in f:
            if line.startswith('C'):
                classified += 1
            elif line.startswith('U'):
                unclassified += 1
    
    total = classified + unclassified
    
    # Count unique taxa from report
    unique_taxa = 0
    with open(kraken_report, 'r') as f:
        for line in f:
            unique_taxa += 1
    
    summary = {
        'total_reads': total,
        'classified_reads': classified,
        'unclassified_reads': unclassified,
        'classification_rate': classified / total if total > 0 else 0,
        'unique_taxa': unique_taxa
    }
    
    return summary


if __name__ == '__main__':
    # If run directly, provide utility functions
    import sys
    import argparse
    
    parser = argparse.ArgumentParser(description='Utility functions for microbiome analysis')
    parser.add_argument('--count-sequences', help='Count sequences in FASTQ file')
    parser.add_argument('--read-stats', help='Calculate read statistics for FASTQ file')
    parser.add_argument('--convert-kraken', nargs=2, metavar=('INPUT', 'OUTPUT'),
                       help='Convert Kraken report to table')
    parser.add_argument('--create-sample-sheet', nargs=2, metavar=('DATA_DIR', 'OUTPUT'),
                       help='Create sample sheet from directory')
    
    args = parser.parse_args()
    
    if args.count_sequences:
        count = count_sequences(args.count_sequences)
        print(f"Number of sequences: {count}")
    
    elif args.read_stats:
        stats = calculate_read_stats(args.read_stats)
        print("=== Read Statistics ===")
        for key, value in stats.items():
            print(f"{key}: {value}")
    
    elif args.convert_kraken:
        convert_kraken_to_table(args.convert_kraken[0], args.convert_kraken[1])
    
    elif args.create_sample_sheet:
        create_sample_sheet(args.create_sample_sheet[0], args.create_sample_sheet[1])
    
    else:
        parser.print_help()
