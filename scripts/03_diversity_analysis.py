#!/usr/bin/env python3
"""
Diversity Analysis Script for Long-Read Microbiome Data
Calculates alpha and beta diversity metrics
"""

import os
import sys
import argparse
import logging
import numpy as np
import pandas as pd
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_taxonomy_table(file_path, file_format='kraken'):
    """
    Load taxonomy abundance table
    
    Args:
        file_path: Path to taxonomy file
        file_format: Format of file (kraken, biom, tsv)
    
    Returns:
        pandas DataFrame with samples as rows and taxa as columns
    """
    logger.info(f"Loading taxonomy table from {file_path}")
    
    if file_format == 'kraken':
        # Parse Kraken report
        data = []
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    percentage = float(parts[0])
                    count = int(parts[1])
                    taxon = parts[5].strip()
                    data.append({'taxon': taxon, 'count': count, 'percentage': percentage})
        
        df = pd.DataFrame(data)
        return df
    
    elif file_format == 'tsv':
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        return df
    
    else:
        logger.error(f"Unsupported file format: {file_format}")
        raise ValueError(f"Unsupported file format: {file_format}")


def calculate_shannon_diversity(abundance_vector):
    """
    Calculate Shannon diversity index
    
    Args:
        abundance_vector: Array of abundance values
    
    Returns:
        Shannon diversity index
    """
    abundance_vector = np.array(abundance_vector)
    abundance_vector = abundance_vector[abundance_vector > 0]
    
    if len(abundance_vector) == 0:
        return 0
    
    proportions = abundance_vector / abundance_vector.sum()
    shannon = -np.sum(proportions * np.log(proportions))
    
    return shannon


def calculate_simpson_diversity(abundance_vector):
    """
    Calculate Simpson diversity index
    
    Args:
        abundance_vector: Array of abundance values
    
    Returns:
        Simpson diversity index (1 - D)
    """
    abundance_vector = np.array(abundance_vector)
    abundance_vector = abundance_vector[abundance_vector > 0]
    
    if len(abundance_vector) == 0:
        return 0
    
    total = abundance_vector.sum()
    if total == 0:
        return 0
    
    proportions = abundance_vector / total
    simpson = 1 - np.sum(proportions ** 2)
    
    return simpson


def calculate_richness(abundance_vector):
    """
    Calculate species richness (number of observed taxa)
    
    Args:
        abundance_vector: Array of abundance values
    
    Returns:
        Number of observed taxa
    """
    abundance_vector = np.array(abundance_vector)
    richness = np.sum(abundance_vector > 0)
    
    return richness


def calculate_evenness(abundance_vector):
    """
    Calculate Pielou's evenness index
    
    Args:
        abundance_vector: Array of abundance values
    
    Returns:
        Pielou's evenness (J')
    """
    shannon = calculate_shannon_diversity(abundance_vector)
    richness = calculate_richness(abundance_vector)
    
    if richness <= 1:
        return 0
    
    evenness = shannon / np.log(richness)
    return evenness


def calculate_alpha_diversity(abundance_table, output_file):
    """
    Calculate multiple alpha diversity metrics
    
    Args:
        abundance_table: DataFrame with samples as rows and taxa as columns
        output_file: Output file for results
    """
    logger.info("Calculating alpha diversity metrics")
    
    results = []
    
    # If single sample (from Kraken report)
    if 'count' in abundance_table.columns:
        counts = abundance_table['count'].values
        
        metrics = {
            'sample': 'sample_1',
            'richness': calculate_richness(counts),
            'shannon': calculate_shannon_diversity(counts),
            'simpson': calculate_simpson_diversity(counts),
            'evenness': calculate_evenness(counts)
        }
        results.append(metrics)
    else:
        # Multiple samples
        for sample in abundance_table.index:
            counts = abundance_table.loc[sample].values
            
            metrics = {
                'sample': sample,
                'richness': calculate_richness(counts),
                'shannon': calculate_shannon_diversity(counts),
                'simpson': calculate_simpson_diversity(counts),
                'evenness': calculate_evenness(counts)
            }
            results.append(metrics)
    
    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False, sep='\t')
    
    logger.info(f"Alpha diversity results saved to {output_file}")
    
    # Print summary
    logger.info("\n=== Alpha Diversity Summary ===")
    for _, row in results_df.iterrows():
        logger.info(f"\nSample: {row['sample']}")
        logger.info(f"  Richness: {row['richness']:.2f}")
        logger.info(f"  Shannon: {row['shannon']:.4f}")
        logger.info(f"  Simpson: {row['simpson']:.4f}")
        logger.info(f"  Evenness: {row['evenness']:.4f}")
    
    return results_df


def calculate_bray_curtis(sample1, sample2):
    """
    Calculate Bray-Curtis dissimilarity
    
    Args:
        sample1: Abundance vector for sample 1
        sample2: Abundance vector for sample 2
    
    Returns:
        Bray-Curtis dissimilarity (0-1)
    """
    sample1 = np.array(sample1)
    sample2 = np.array(sample2)
    
    numerator = np.sum(np.abs(sample1 - sample2))
    denominator = np.sum(sample1 + sample2)
    
    if denominator == 0:
        return 0
    
    return numerator / denominator


def calculate_jaccard(sample1, sample2):
    """
    Calculate Jaccard distance
    
    Args:
        sample1: Abundance vector for sample 1
        sample2: Abundance vector for sample 2
    
    Returns:
        Jaccard distance (0-1)
    """
    sample1 = np.array(sample1) > 0
    sample2 = np.array(sample2) > 0
    
    intersection = np.sum(sample1 & sample2)
    union = np.sum(sample1 | sample2)
    
    if union == 0:
        return 0
    
    jaccard_similarity = intersection / union
    return 1 - jaccard_similarity


def calculate_beta_diversity(abundance_table, output_file, metric='bray_curtis'):
    """
    Calculate beta diversity distance matrix
    
    Args:
        abundance_table: DataFrame with samples as rows and taxa as columns
        output_file: Output file for distance matrix
        metric: Distance metric (bray_curtis, jaccard)
    """
    logger.info(f"Calculating beta diversity using {metric}")
    
    # For single sample from Kraken report, skip beta diversity
    if 'count' in abundance_table.columns:
        logger.warning("Beta diversity requires multiple samples. Skipping.")
        return None
    
    samples = abundance_table.index.tolist()
    n_samples = len(samples)
    
    # Initialize distance matrix
    distance_matrix = np.zeros((n_samples, n_samples))
    
    # Calculate pairwise distances
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            sample1 = abundance_table.iloc[i].values
            sample2 = abundance_table.iloc[j].values
            
            if metric == 'bray_curtis':
                dist = calculate_bray_curtis(sample1, sample2)
            elif metric == 'jaccard':
                dist = calculate_jaccard(sample1, sample2)
            else:
                raise ValueError(f"Unknown metric: {metric}")
            
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist
    
    # Save distance matrix
    distance_df = pd.DataFrame(
        distance_matrix,
        index=samples,
        columns=samples
    )
    distance_df.to_csv(output_file, sep='\t')
    
    logger.info(f"Beta diversity distance matrix saved to {output_file}")
    
    return distance_df


def perform_rarefaction(abundance_table, output_file, max_depth=None, step=100):
    """
    Perform rarefaction analysis
    
    Args:
        abundance_table: DataFrame with samples as rows and taxa as columns
        output_file: Output file for rarefaction curves
        max_depth: Maximum sequencing depth (default: minimum sample size)
        step: Step size for rarefaction
    """
    logger.info("Performing rarefaction analysis")
    
    # Handle single sample from Kraken report
    if 'count' in abundance_table.columns:
        counts = abundance_table['count'].values
        total_reads = int(counts.sum())
        
        if max_depth is None:
            max_depth = total_reads
        
        depths = range(step, min(max_depth, total_reads) + 1, step)
        rarefaction_data = []
        
        for depth in depths:
            # Simple rarefaction simulation
            rarefied = np.random.choice(
                len(counts),
                size=depth,
                replace=True,
                p=counts / counts.sum()
            )
            richness = len(np.unique(rarefied))
            
            rarefaction_data.append({
                'sample': 'sample_1',
                'depth': depth,
                'richness': richness
            })
        
        rarefaction_df = pd.DataFrame(rarefaction_data)
    else:
        # Multiple samples
        rarefaction_data = []
        
        for sample in abundance_table.index:
            counts = abundance_table.loc[sample].values
            total_reads = int(counts.sum())
            
            sample_max_depth = min(max_depth, total_reads) if max_depth else total_reads
            depths = range(step, sample_max_depth + 1, step)
            
            for depth in depths:
                rarefied = np.random.choice(
                    len(counts),
                    size=depth,
                    replace=True,
                    p=counts / counts.sum()
                )
                richness = len(np.unique(rarefied))
                
                rarefaction_data.append({
                    'sample': sample,
                    'depth': depth,
                    'richness': richness
                })
        
        rarefaction_df = pd.DataFrame(rarefaction_data)
    
    rarefaction_df.to_csv(output_file, index=False, sep='\t')
    logger.info(f"Rarefaction results saved to {output_file}")
    
    return rarefaction_df


def main():
    parser = argparse.ArgumentParser(
        description='Diversity analysis for microbiome data'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input taxonomy/abundance file'
    )
    parser.add_argument(
        '-o', '--output-dir',
        required=True,
        help='Output directory for diversity results'
    )
    parser.add_argument(
        '--format',
        choices=['kraken', 'tsv', 'biom'],
        default='kraken',
        help='Input file format (default: kraken)'
    )
    parser.add_argument(
        '--alpha',
        action='store_true',
        default=True,
        help='Calculate alpha diversity (default: True)'
    )
    parser.add_argument(
        '--beta',
        action='store_true',
        help='Calculate beta diversity (requires multiple samples)'
    )
    parser.add_argument(
        '--beta-metric',
        choices=['bray_curtis', 'jaccard'],
        default='bray_curtis',
        help='Beta diversity metric (default: bray_curtis)'
    )
    parser.add_argument(
        '--rarefaction',
        action='store_true',
        help='Perform rarefaction analysis'
    )
    parser.add_argument(
        '--rare-depth',
        type=int,
        help='Maximum rarefaction depth (default: min sample size)'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    abundance_table = load_taxonomy_table(args.input, args.format)
    
    # Calculate alpha diversity
    if args.alpha:
        alpha_output = os.path.join(args.output_dir, 'alpha_diversity.tsv')
        calculate_alpha_diversity(abundance_table, alpha_output)
    
    # Calculate beta diversity
    if args.beta:
        beta_output = os.path.join(args.output_dir, f'beta_diversity_{args.beta_metric}.tsv')
        calculate_beta_diversity(abundance_table, beta_output, args.beta_metric)
    
    # Rarefaction analysis
    if args.rarefaction:
        rare_output = os.path.join(args.output_dir, 'rarefaction_curves.tsv')
        perform_rarefaction(abundance_table, rare_output, args.rare_depth)
    
    logger.info("Diversity analysis completed successfully!")


if __name__ == '__main__':
    main()
