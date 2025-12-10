#!/usr/bin/env python3
"""
Statistical Analysis Script for Microbiome Data
Performs differential abundance analysis and statistical comparisons
"""

import os
import sys
import argparse
import logging
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_abundance_data(file_path):
    """
    Load abundance data from file
    
    Args:
        file_path: Path to abundance file (TSV format)
    
    Returns:
        pandas DataFrame
    """
    logger.info(f"Loading abundance data from {file_path}")
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    return df


def load_metadata(file_path):
    """
    Load sample metadata
    
    Args:
        file_path: Path to metadata file (TSV format)
        Expected columns: sample_id, group, [other metadata]
    
    Returns:
        pandas DataFrame
    """
    logger.info(f"Loading metadata from {file_path}")
    df = pd.read_csv(file_path, sep='\t')
    return df


def normalize_abundance(abundance_df, method='total_sum'):
    """
    Normalize abundance data
    
    Args:
        abundance_df: DataFrame with samples as rows, taxa as columns
        method: Normalization method (total_sum, css, clr)
    
    Returns:
        Normalized DataFrame
    """
    logger.info(f"Normalizing abundance data using {method}")
    
    if method == 'total_sum':
        # Total Sum Scaling (TSS) - convert to relative abundance
        normalized = abundance_df.div(abundance_df.sum(axis=1), axis=0)
    
    elif method == 'css':
        # Cumulative Sum Scaling
        # Simple implementation: normalize by quantile sum
        quantile_sums = abundance_df.sum(axis=1).quantile(0.5)
        normalized = abundance_df.div(abundance_df.sum(axis=1), axis=0) * quantile_sums
    
    elif method == 'clr':
        # Centered Log-Ratio transformation
        # Add pseudocount to avoid log(0)
        pseudocount = 1
        abundance_pseudo = abundance_df + pseudocount
        
        # Calculate geometric mean for each sample
        geo_means = np.exp(np.log(abundance_pseudo).mean(axis=1))
        
        # CLR transformation
        normalized = np.log(abundance_pseudo.div(geo_means, axis=0))
    
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    
    return normalized


def filter_low_abundance(abundance_df, min_count=10, min_prevalence=0.1):
    """
    Filter taxa with low abundance or prevalence
    
    Args:
        abundance_df: DataFrame with samples as rows, taxa as columns
        min_count: Minimum total count across all samples
        min_prevalence: Minimum proportion of samples with non-zero counts
    
    Returns:
        Filtered DataFrame
    """
    logger.info("Filtering low abundance taxa")
    
    # Filter by total count
    total_counts = abundance_df.sum(axis=0)
    passing_count = total_counts >= min_count
    
    # Filter by prevalence
    prevalence = (abundance_df > 0).sum(axis=0) / len(abundance_df)
    passing_prevalence = prevalence >= min_prevalence
    
    # Combined filter
    passing = passing_count & passing_prevalence
    
    filtered_df = abundance_df.loc[:, passing]
    
    logger.info(f"Kept {passing.sum()} / {len(abundance_df.columns)} taxa")
    
    return filtered_df


def compare_alpha_diversity(alpha_df, metadata_df, group_column, output_file):
    """
    Compare alpha diversity between groups
    
    Args:
        alpha_df: DataFrame with alpha diversity metrics
        metadata_df: DataFrame with sample metadata
        group_column: Column name for grouping variable
        output_file: Output file for results
    """
    logger.info(f"Comparing alpha diversity between groups ({group_column})")
    
    # Merge alpha diversity with metadata
    merged = alpha_df.merge(
        metadata_df[['sample', group_column]],
        left_on='sample',
        right_on='sample'
    )
    
    # Get unique groups
    groups = merged[group_column].unique()
    
    if len(groups) < 2:
        logger.warning("Need at least 2 groups for comparison")
        return
    
    results = []
    
    # Compare each metric
    metrics = ['richness', 'shannon', 'simpson', 'evenness']
    
    for metric in metrics:
        if metric not in merged.columns:
            continue
        
        group_data = [merged[merged[group_column] == g][metric].values for g in groups]
        
        # Perform statistical test
        if len(groups) == 2:
            # Mann-Whitney U test (non-parametric)
            stat, pvalue = stats.mannwhitneyu(group_data[0], group_data[1])
            test_name = 'Mann-Whitney U'
        else:
            # Kruskal-Wallis test (non-parametric)
            stat, pvalue = stats.kruskal(*group_data)
            test_name = 'Kruskal-Wallis'
        
        # Calculate means for each group
        group_means = {g: merged[merged[group_column] == g][metric].mean() for g in groups}
        
        results.append({
            'metric': metric,
            'test': test_name,
            'statistic': stat,
            'pvalue': pvalue,
            'significant': pvalue < 0.05,
            **{f'mean_{g}': group_means[g] for g in groups}
        })
    
    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False, sep='\t')
    
    logger.info(f"Alpha diversity comparison saved to {output_file}")
    
    # Print summary
    logger.info("\n=== Alpha Diversity Comparison ===")
    for _, row in results_df.iterrows():
        logger.info(f"\n{row['metric']} ({row['test']}):")
        logger.info(f"  p-value: {row['pvalue']:.4f} {'*' if row['significant'] else ''}")
        for g in groups:
            logger.info(f"  Mean {g}: {row[f'mean_{g}']:.4f}")


def compare_beta_diversity(distance_matrix, metadata_df, group_column, output_file):
    """
    Compare beta diversity between groups using PERMANOVA
    
    Args:
        distance_matrix: Distance matrix DataFrame
        metadata_df: DataFrame with sample metadata
        group_column: Column name for grouping variable
        output_file: Output file for results
    """
    logger.info(f"Comparing beta diversity between groups ({group_column})")
    
    # Simple PERMANOVA-like analysis
    # Note: Full PERMANOVA requires specialized libraries
    
    samples = distance_matrix.index.tolist()
    
    # Get group assignments
    sample_groups = {}
    for sample in samples:
        group = metadata_df[metadata_df['sample'] == sample][group_column].values
        if len(group) > 0:
            sample_groups[sample] = group[0]
    
    # Calculate within-group and between-group distances
    within_distances = []
    between_distances = []
    
    for i, sample1 in enumerate(samples):
        for j, sample2 in enumerate(samples):
            if i >= j:
                continue
            
            if sample1 not in sample_groups or sample2 not in sample_groups:
                continue
            
            distance = distance_matrix.loc[sample1, sample2]
            
            if sample_groups[sample1] == sample_groups[sample2]:
                within_distances.append(distance)
            else:
                between_distances.append(distance)
    
    # Statistical test
    if len(within_distances) > 0 and len(between_distances) > 0:
        stat, pvalue = stats.mannwhitneyu(within_distances, between_distances)
        
        results = {
            'test': 'Mann-Whitney U (within vs between)',
            'statistic': stat,
            'pvalue': pvalue,
            'mean_within': np.mean(within_distances),
            'mean_between': np.mean(between_distances),
            'significant': pvalue < 0.05
        }
        
        # Save results
        results_df = pd.DataFrame([results])
        results_df.to_csv(output_file, index=False, sep='\t')
        
        logger.info(f"Beta diversity comparison saved to {output_file}")
        logger.info(f"\nMean within-group distance: {results['mean_within']:.4f}")
        logger.info(f"Mean between-group distance: {results['mean_between']:.4f}")
        logger.info(f"p-value: {pvalue:.4f} {'*' if pvalue < 0.05 else ''}")
    else:
        logger.warning("Insufficient data for beta diversity comparison")


def differential_abundance_analysis(abundance_df, metadata_df, group_column, output_file):
    """
    Perform differential abundance analysis between groups
    
    Args:
        abundance_df: DataFrame with samples as rows, taxa as columns
        metadata_df: DataFrame with sample metadata
        group_column: Column name for grouping variable
        output_file: Output file for results
    """
    logger.info(f"Performing differential abundance analysis ({group_column})")
    
    # Merge abundance with metadata
    abundance_with_meta = abundance_df.copy()
    abundance_with_meta['sample'] = abundance_with_meta.index
    
    merged = abundance_with_meta.merge(
        metadata_df[['sample', group_column]],
        on='sample'
    )
    
    groups = merged[group_column].unique()
    
    if len(groups) != 2:
        logger.warning("Differential abundance currently only supports 2-group comparison")
        return
    
    results = []
    
    # Test each taxon
    taxa = [col for col in abundance_df.columns]
    
    for taxon in taxa:
        group1_data = merged[merged[group_column] == groups[0]][taxon].values
        group2_data = merged[merged[group_column] == groups[1]][taxon].values
        
        # Mann-Whitney U test
        stat, pvalue = stats.mannwhitneyu(group1_data, group2_data, alternative='two-sided')
        
        # Calculate fold change
        mean1 = np.mean(group1_data) + 1e-10  # pseudocount
        mean2 = np.mean(group2_data) + 1e-10
        fold_change = mean2 / mean1
        log2_fc = np.log2(fold_change)
        
        results.append({
            'taxon': taxon,
            f'mean_{groups[0]}': np.mean(group1_data),
            f'mean_{groups[1]}': np.mean(group2_data),
            'fold_change': fold_change,
            'log2_fold_change': log2_fc,
            'statistic': stat,
            'pvalue': pvalue
        })
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Multiple testing correction (Benjamini-Hochberg)
    from statsmodels.stats.multitest import multipletests
    
    _, results_df['padj'], _, _ = multipletests(
        results_df['pvalue'],
        method='fdr_bh'
    )
    
    results_df['significant'] = results_df['padj'] < 0.05
    
    # Sort by p-value
    results_df = results_df.sort_values('pvalue')
    
    # Save results
    results_df.to_csv(output_file, index=False, sep='\t')
    
    logger.info(f"Differential abundance results saved to {output_file}")
    logger.info(f"Found {results_df['significant'].sum()} significantly different taxa (padj < 0.05)")


def main():
    parser = argparse.ArgumentParser(
        description='Statistical analysis for microbiome data'
    )
    parser.add_argument(
        '-a', '--abundance',
        required=True,
        help='Abundance table file (TSV format)'
    )
    parser.add_argument(
        '-m', '--metadata',
        required=True,
        help='Metadata file (TSV format with sample and group columns)'
    )
    parser.add_argument(
        '-o', '--output-dir',
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '-g', '--group',
        required=True,
        help='Metadata column for group comparison'
    )
    parser.add_argument(
        '--alpha-diversity',
        help='Alpha diversity file for comparison (optional)'
    )
    parser.add_argument(
        '--beta-diversity',
        help='Beta diversity distance matrix for comparison (optional)'
    )
    parser.add_argument(
        '--normalize',
        choices=['total_sum', 'css', 'clr'],
        default='total_sum',
        help='Normalization method (default: total_sum)'
    )
    parser.add_argument(
        '--min-count',
        type=int,
        default=10,
        help='Minimum total count for taxa filtering (default: 10)'
    )
    parser.add_argument(
        '--min-prevalence',
        type=float,
        default=0.1,
        help='Minimum prevalence for taxa filtering (default: 0.1)'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    abundance_df = load_abundance_data(args.abundance)
    metadata_df = load_metadata(args.metadata)
    
    # Filter low abundance taxa
    abundance_filtered = filter_low_abundance(
        abundance_df,
        args.min_count,
        args.min_prevalence
    )
    
    # Normalize abundance
    abundance_normalized = normalize_abundance(abundance_filtered, args.normalize)
    
    # Differential abundance analysis
    da_output = os.path.join(args.output_dir, 'differential_abundance.tsv')
    differential_abundance_analysis(
        abundance_normalized,
        metadata_df,
        args.group,
        da_output
    )
    
    # Alpha diversity comparison
    if args.alpha_diversity:
        alpha_df = pd.read_csv(args.alpha_diversity, sep='\t')
        alpha_output = os.path.join(args.output_dir, 'alpha_diversity_comparison.tsv')
        compare_alpha_diversity(alpha_df, metadata_df, args.group, alpha_output)
    
    # Beta diversity comparison
    if args.beta_diversity:
        beta_df = pd.read_csv(args.beta_diversity, sep='\t', index_col=0)
        beta_output = os.path.join(args.output_dir, 'beta_diversity_comparison.tsv')
        compare_beta_diversity(beta_df, metadata_df, args.group, beta_output)
    
    logger.info("Statistical analysis completed successfully!")


if __name__ == '__main__':
    main()
