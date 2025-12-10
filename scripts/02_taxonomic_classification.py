#!/usr/bin/env python3
"""
Taxonomic Classification Script for Long-Read Microbiome Data
Supports multiple classifiers for 16S rRNA and shotgun metagenomic data
"""

import os
import sys
import argparse
import subprocess
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def run_kraken2(input_file, output_file, report_file, database_path, threads=4):
    """
    Run Kraken2 for taxonomic classification
    
    Args:
        input_file: Input FASTQ/FASTA file
        output_file: Output classification file
        report_file: Output report file
        database_path: Path to Kraken2 database
        threads: Number of threads to use
    """
    logger.info("Running Kraken2 taxonomic classification")
    
    cmd = [
        'kraken2',
        '--db', database_path,
        '--threads', str(threads),
        '--output', output_file,
        '--report', report_file,
        '--use-names',
        input_file
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Kraken2 completed. Results: {output_file}, Report: {report_file}")
    except FileNotFoundError:
        logger.error("Kraken2 not found. Please install Kraken2.")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(f"Kraken2 failed: {e.stderr}")
        raise


def run_bracken(kraken_report, output_file, database_path, read_length=1500, level='S', threshold=10):
    """
    Run Bracken for abundance estimation
    
    Args:
        kraken_report: Kraken2 report file
        output_file: Output Bracken file
        database_path: Path to Kraken2 database with Bracken files
        read_length: Read length for Bracken (default: 1500 for long reads)
        level: Taxonomic level (S=species, G=genus, etc.)
        threshold: Minimum number of reads required
    """
    logger.info(f"Running Bracken abundance estimation at level {level}")
    
    cmd = [
        'bracken',
        '-d', database_path,
        '-i', kraken_report,
        '-o', output_file,
        '-r', str(read_length),
        '-l', level,
        '-t', str(threshold)
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"Bracken completed. Results: {output_file}")
    except FileNotFoundError:
        logger.warning("Bracken not found. Skipping abundance re-estimation.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Bracken failed: {e.stderr}")
        raise


def run_minimap2_16s(input_file, output_file, reference_db, threads=4):
    """
    Run Minimap2 for 16S rRNA classification
    
    Args:
        input_file: Input FASTQ file
        output_file: Output SAM/BAM file
        reference_db: Path to 16S reference database
        threads: Number of threads
    """
    logger.info("Running Minimap2 for 16S rRNA classification")
    
    cmd = [
        'minimap2',
        '-ax', 'map-pb',  # PacBio mode
        '-t', str(threads),
        reference_db,
        input_file
    ]
    
    try:
        with open(output_file, 'w') as out:
            subprocess.run(cmd, check=True, stdout=out, stderr=subprocess.PIPE, text=True)
        logger.info(f"Minimap2 completed. Results: {output_file}")
    except FileNotFoundError:
        logger.error("Minimap2 not found. Please install Minimap2.")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(f"Minimap2 failed: {e.stderr}")
        raise


def parse_classification_results(kraken_output, output_summary):
    """
    Parse and summarize classification results
    
    Args:
        kraken_output: Kraken2 output file
        output_summary: Output summary file
    """
    logger.info("Parsing classification results")
    
    classified = 0
    unclassified = 0
    taxa_counts = {}
    
    with open(kraken_output, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                status = parts[0]
                taxon = parts[2] if len(parts) > 2 else 'Unknown'
                
                if status == 'C':
                    classified += 1
                    taxa_counts[taxon] = taxa_counts.get(taxon, 0) + 1
                else:
                    unclassified += 1
    
    total = classified + unclassified
    
    with open(output_summary, 'w') as f:
        f.write("=== Taxonomic Classification Summary ===\n\n")
        f.write(f"Total reads: {total}\n")
        f.write(f"Classified: {classified} ({classified/total*100:.2f}%)\n")
        f.write(f"Unclassified: {unclassified} ({unclassified/total*100:.2f}%)\n\n")
        f.write(f"Unique taxa identified: {len(taxa_counts)}\n\n")
        
        f.write("Top 20 Taxa:\n")
        sorted_taxa = sorted(taxa_counts.items(), key=lambda x: x[1], reverse=True)
        for i, (taxon, count) in enumerate(sorted_taxa[:20], 1):
            f.write(f"{i}. {taxon}: {count} reads\n")
    
    logger.info(f"Summary saved to {output_summary}")


def convert_kraken_to_biom(kraken_report, output_biom):
    """
    Convert Kraken report to BIOM format for downstream analysis
    
    Args:
        kraken_report: Kraken2 report file
        output_biom: Output BIOM file
    """
    logger.info("Converting Kraken report to BIOM format")
    
    cmd = [
        'kraken-biom',
        kraken_report,
        '-o', output_biom,
        '--fmt', 'json'
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"BIOM file created: {output_biom}")
    except FileNotFoundError:
        logger.warning("kraken-biom not found. Skipping BIOM conversion.")
    except subprocess.CalledProcessError as e:
        logger.warning(f"BIOM conversion failed: {e.stderr}")


def main():
    parser = argparse.ArgumentParser(
        description='Taxonomic classification for long-read microbiome sequences'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input FASTQ file (quality-filtered)'
    )
    parser.add_argument(
        '-o', '--output-dir',
        required=True,
        help='Output directory for classification results'
    )
    parser.add_argument(
        '-d', '--database',
        required=True,
        help='Path to Kraken2 database'
    )
    parser.add_argument(
        '--method',
        choices=['kraken2', 'minimap2'],
        default='kraken2',
        help='Classification method (default: kraken2)'
    )
    parser.add_argument(
        '--16s-db',
        help='Path to 16S reference database (for minimap2 method)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Number of threads (default: 4)'
    )
    parser.add_argument(
        '--run-bracken',
        action='store_true',
        help='Run Bracken for abundance re-estimation'
    )
    parser.add_argument(
        '--bracken-level',
        default='S',
        choices=['S', 'G', 'F', 'O', 'C', 'P'],
        help='Bracken taxonomic level (default: S for species)'
    )
    parser.add_argument(
        '--read-length',
        type=int,
        default=1500,
        help='Average read length for Bracken (default: 1500)'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.method == 'kraken2':
        # Run Kraken2
        kraken_output = os.path.join(args.output_dir, 'kraken2_output.txt')
        kraken_report = os.path.join(args.output_dir, 'kraken2_report.txt')
        
        run_kraken2(
            args.input,
            kraken_output,
            kraken_report,
            args.database,
            args.threads
        )
        
        # Parse results
        summary_file = os.path.join(args.output_dir, 'classification_summary.txt')
        parse_classification_results(kraken_output, summary_file)
        
        # Run Bracken if requested
        if args.run_bracken:
            bracken_output = os.path.join(
                args.output_dir,
                f'bracken_output_{args.bracken_level}.txt'
            )
            run_bracken(
                kraken_report,
                bracken_output,
                args.database,
                args.read_length,
                args.bracken_level
            )
        
        # Convert to BIOM
        biom_file = os.path.join(args.output_dir, 'taxonomy.biom')
        convert_kraken_to_biom(kraken_report, biom_file)
        
    elif args.method == 'minimap2':
        if not args.s16_db:
            logger.error("--16s-db required for minimap2 method")
            sys.exit(1)
        
        sam_output = os.path.join(args.output_dir, 'minimap2_alignment.sam')
        run_minimap2_16s(
            args.input,
            sam_output,
            args.s16_db,
            args.threads
        )
    
    logger.info("Taxonomic classification completed successfully!")


if __name__ == '__main__':
    main()
