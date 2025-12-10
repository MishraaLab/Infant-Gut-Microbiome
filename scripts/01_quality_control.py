#!/usr/bin/env python3
"""
Quality Control Script for Long-Read Microbiome Data
Performs QC on PacBio long-read sequences using NanoPlot and custom filters
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def run_nanoplot(input_file, output_dir, file_format='fastq'):
    """
    Run NanoPlot for quality visualization
    
    Args:
        input_file: Path to input sequence file
        output_dir: Output directory for QC results
        file_format: Format of input file (fastq, fasta, bam)
    """
    logger.info(f"Running NanoPlot on {input_file}")
    
    cmd = [
        'NanoPlot',
        '--fastq' if file_format == 'fastq' else '--fasta', input_file,
        '-o', output_dir,
        '--plots', 'dot',
        '--legacy', 'hex',
        '-t', '4'  # threads
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"NanoPlot completed successfully. Results in {output_dir}")
    except subprocess.CalledProcessError as e:
        logger.error(f"NanoPlot failed: {e.stderr}")
        raise
    except FileNotFoundError:
        logger.warning("NanoPlot not found. Skipping visualization step.")


def filter_by_length(input_file, output_file, min_length=1000, max_length=None):
    """
    Filter sequences by length
    
    Args:
        input_file: Input FASTQ file
        output_file: Output filtered FASTQ file
        min_length: Minimum sequence length (default: 1000)
        max_length: Maximum sequence length (optional)
    """
    logger.info(f"Filtering sequences by length (min={min_length}, max={max_length})")
    
    from Bio import SeqIO
    
    count_input = 0
    count_output = 0
    
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, 'fastq'):
            count_input += 1
            seq_len = len(record.seq)
            
            if seq_len >= min_length:
                if max_length is None or seq_len <= max_length:
                    SeqIO.write(record, out_handle, 'fastq')
                    count_output += 1
    
    logger.info(f"Filtered {count_input} sequences -> {count_output} sequences")
    if count_input > 0:
        logger.info(f"Retention rate: {count_output/count_input*100:.2f}%")
    else:
        logger.info("Retention rate: N/A (no input sequences)")


def filter_by_quality(input_file, output_file, min_quality=10):
    """
    Filter sequences by average quality score
    
    Args:
        input_file: Input FASTQ file
        output_file: Output filtered FASTQ file
        min_quality: Minimum average quality score (default: 10)
    """
    logger.info(f"Filtering sequences by quality (min_quality={min_quality})")
    
    from Bio import SeqIO
    import numpy as np
    
    count_input = 0
    count_output = 0
    
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, 'fastq'):
            count_input += 1
            # Calculate average quality score
            avg_quality = np.mean(record.letter_annotations["phred_quality"])
            
            if avg_quality >= min_quality:
                SeqIO.write(record, out_handle, 'fastq')
                count_output += 1
    
    logger.info(f"Filtered {count_input} sequences -> {count_output} sequences")
    if count_input > 0:
        logger.info(f"Retention rate: {count_output/count_input*100:.2f}%")
    else:
        logger.info("Retention rate: N/A (no input sequences)")


def generate_qc_report(input_file, output_file):
    """
    Generate a summary QC report
    
    Args:
        input_file: Input FASTQ file
        output_file: Output report file
    """
    logger.info("Generating QC summary report")
    
    from Bio import SeqIO
    import numpy as np
    
    lengths = []
    qualities = []
    
    for record in SeqIO.parse(input_file, 'fastq'):
        lengths.append(len(record.seq))
        qualities.append(np.mean(record.letter_annotations["phred_quality"]))
    
    with open(output_file, 'w') as f:
        f.write("=== Quality Control Summary ===\n\n")
        f.write(f"Total sequences: {len(lengths)}\n")
        f.write(f"Total bases: {sum(lengths)}\n\n")
        
        f.write("Length Statistics:\n")
        f.write(f"  Mean: {np.mean(lengths):.2f} bp\n")
        f.write(f"  Median: {np.median(lengths):.2f} bp\n")
        f.write(f"  Min: {min(lengths)} bp\n")
        f.write(f"  Max: {max(lengths)} bp\n")
        f.write(f"  N50: {calculate_n50(lengths):.2f} bp\n\n")
        
        f.write("Quality Statistics:\n")
        f.write(f"  Mean Q-score: {np.mean(qualities):.2f}\n")
        f.write(f"  Median Q-score: {np.median(qualities):.2f}\n")
        f.write(f"  Min Q-score: {min(qualities):.2f}\n")
        f.write(f"  Max Q-score: {max(qualities):.2f}\n")
    
    logger.info(f"QC report saved to {output_file}")


def calculate_n50(lengths):
    """Calculate N50 statistic"""
    lengths_sorted = sorted(lengths, reverse=True)
    total_length = sum(lengths_sorted)
    target = total_length / 2
    
    cumsum = 0
    for length in lengths_sorted:
        cumsum += length
        if cumsum >= target:
            return length
    return 0


def main():
    parser = argparse.ArgumentParser(
        description='Quality control for long-read microbiome sequences'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input FASTQ file'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output filtered FASTQ file'
    )
    parser.add_argument(
        '--qc-dir',
        default='qc_results',
        help='Directory for QC visualizations (default: qc_results)'
    )
    parser.add_argument(
        '--min-length',
        type=int,
        default=1000,
        help='Minimum sequence length (default: 1000)'
    )
    parser.add_argument(
        '--max-length',
        type=int,
        help='Maximum sequence length (optional)'
    )
    parser.add_argument(
        '--min-quality',
        type=float,
        default=10,
        help='Minimum average quality score (default: 10)'
    )
    parser.add_argument(
        '--skip-nanoplot',
        action='store_true',
        help='Skip NanoPlot visualization'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.qc_dir, exist_ok=True)
    
    # Run NanoPlot for initial QC
    if not args.skip_nanoplot:
        run_nanoplot(args.input, os.path.join(args.qc_dir, 'raw'), 'fastq')
    
    # Filter by length
    temp_file = args.output + '.temp'
    filter_by_length(args.input, temp_file, args.min_length, args.max_length)
    
    # Filter by quality
    filter_by_quality(temp_file, args.output, args.min_quality)
    
    # Clean up temp file
    if os.path.exists(temp_file):
        os.remove(temp_file)
    
    # Run NanoPlot on filtered data
    if not args.skip_nanoplot:
        run_nanoplot(args.output, os.path.join(args.qc_dir, 'filtered'), 'fastq')
    
    # Generate summary report
    report_file = os.path.join(args.qc_dir, 'qc_summary.txt')
    generate_qc_report(args.output, report_file)
    
    logger.info("Quality control completed successfully!")


if __name__ == '__main__':
    main()
