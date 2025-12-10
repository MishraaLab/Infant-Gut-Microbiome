#!/usr/bin/env python3
"""
Visualization Script for Microbiome Data
Creates taxonomic composition plots, diversity plots, and heatmaps
"""

import os
import sys
import argparse
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300


def plot_taxonomic_barplot(abundance_df, metadata_df, output_file, top_n=20, group_column=None):
    """
    Create stacked barplot of taxonomic composition
    
    Args:
        abundance_df: DataFrame with samples as rows, taxa as columns
        metadata_df: DataFrame with sample metadata (optional)
        output_file: Output file for plot
        top_n: Number of top taxa to display
        group_column: Column for grouping/coloring samples
    """
    logger.info(f"Creating taxonomic barplot (top {top_n} taxa)")
    
    # Get top N most abundant taxa
    mean_abundance = abundance_df.mean(axis=0)
    top_taxa = mean_abundance.nlargest(top_n).index.tolist()
    
    # Create abundance table with top taxa + "Other"
    plot_data = abundance_df[top_taxa].copy()
    plot_data['Other'] = abundance_df.drop(columns=top_taxa).sum(axis=1)
    
    # Normalize to relative abundance
    plot_data = plot_data.div(plot_data.sum(axis=1), axis=0) * 100
    
    # Sort samples by group if provided
    if group_column and metadata_df is not None:
        merged = plot_data.merge(
            metadata_df[['sample', group_column]],
            left_index=True,
            right_on='sample',
            how='left'
        )
        merged = merged.sort_values(group_column)
        plot_data = merged.drop(columns=['sample', group_column])
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Create color palette
    colors = sns.color_palette("tab20", len(plot_data.columns))
    
    # Plot stacked bar chart
    plot_data.T.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        color=colors,
        width=0.8,
        legend=False
    )
    
    ax.set_xlabel('Sample')
    ax.set_ylabel('Relative Abundance (%)')
    ax.set_title('Taxonomic Composition')
    ax.legend(
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        fontsize=8,
        frameon=True
    )
    
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Taxonomic barplot saved to {output_file}")


def plot_alpha_diversity(alpha_df, metadata_df, output_file, group_column=None):
    """
    Create boxplots for alpha diversity metrics
    
    Args:
        alpha_df: DataFrame with alpha diversity metrics
        metadata_df: DataFrame with sample metadata (optional)
        output_file: Output file for plot
        group_column: Column for grouping samples
    """
    logger.info("Creating alpha diversity plots")
    
    metrics = ['richness', 'shannon', 'simpson', 'evenness']
    available_metrics = [m for m in metrics if m in alpha_df.columns]
    
    if not available_metrics:
        logger.warning("No alpha diversity metrics found")
        return
    
    # Merge with metadata if provided
    if group_column and metadata_df is not None:
        plot_data = alpha_df.merge(
            metadata_df[['sample', group_column]],
            on='sample',
            how='left'
        )
    else:
        plot_data = alpha_df.copy()
        group_column = None
    
    # Create subplots
    n_metrics = len(available_metrics)
    fig, axes = plt.subplots(1, n_metrics, figsize=(4 * n_metrics, 5))
    
    if n_metrics == 1:
        axes = [axes]
    
    for i, metric in enumerate(available_metrics):
        ax = axes[i]
        
        if group_column and group_column in plot_data.columns:
            sns.boxplot(
                data=plot_data,
                x=group_column,
                y=metric,
                ax=ax,
                palette='Set2'
            )
            sns.swarmplot(
                data=plot_data,
                x=group_column,
                y=metric,
                ax=ax,
                color='black',
                alpha=0.5,
                size=3
            )
        else:
            sns.boxplot(
                data=plot_data,
                y=metric,
                ax=ax,
                color='lightblue'
            )
        
        ax.set_title(metric.capitalize())
        ax.set_ylabel(metric.capitalize())
        
        if group_column:
            ax.set_xlabel('Group')
        else:
            ax.set_xlabel('')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Alpha diversity plots saved to {output_file}")


def plot_rarefaction_curves(rarefaction_df, output_file):
    """
    Create rarefaction curves
    
    Args:
        rarefaction_df: DataFrame with rarefaction data
        output_file: Output file for plot
    """
    logger.info("Creating rarefaction curves")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    samples = rarefaction_df['sample'].unique()
    
    for sample in samples:
        sample_data = rarefaction_df[rarefaction_df['sample'] == sample]
        ax.plot(
            sample_data['depth'],
            sample_data['richness'],
            marker='o',
            markersize=3,
            label=sample
        )
    
    ax.set_xlabel('Sequencing Depth')
    ax.set_ylabel('Observed Richness')
    ax.set_title('Rarefaction Curves')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Rarefaction curves saved to {output_file}")


def plot_pcoa(distance_matrix, metadata_df, output_file, group_column=None):
    """
    Create PCoA (Principal Coordinates Analysis) plot
    
    Args:
        distance_matrix: Distance matrix DataFrame
        metadata_df: DataFrame with sample metadata (optional)
        output_file: Output file for plot
        group_column: Column for coloring points
    """
    logger.info("Creating PCoA plot")
    
    from sklearn.manifold import MDS
    
    # Perform MDS (equivalent to PCoA for distance matrices)
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    coords = mds.fit_transform(distance_matrix.values)
    
    # Create DataFrame with coordinates
    pcoa_df = pd.DataFrame(
        coords,
        columns=['PC1', 'PC2'],
        index=distance_matrix.index
    )
    pcoa_df['sample'] = pcoa_df.index
    
    # Merge with metadata if provided
    if group_column and metadata_df is not None:
        pcoa_df = pcoa_df.merge(
            metadata_df[['sample', group_column]],
            on='sample',
            how='left'
        )
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if group_column and group_column in pcoa_df.columns:
        groups = pcoa_df[group_column].unique()
        colors = sns.color_palette("Set2", len(groups))
        
        for i, group in enumerate(groups):
            group_data = pcoa_df[pcoa_df[group_column] == group]
            ax.scatter(
                group_data['PC1'],
                group_data['PC2'],
                c=[colors[i]],
                label=group,
                s=100,
                alpha=0.7,
                edgecolors='black'
            )
    else:
        ax.scatter(
            pcoa_df['PC1'],
            pcoa_df['PC2'],
            s=100,
            alpha=0.7,
            edgecolors='black'
        )
    
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title('Principal Coordinates Analysis (PCoA)')
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.3)
    
    if group_column:
        ax.legend(title='Group')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"PCoA plot saved to {output_file}")


def plot_heatmap(abundance_df, metadata_df, output_file, top_n=50, group_column=None):
    """
    Create heatmap of abundance data
    
    Args:
        abundance_df: DataFrame with samples as rows, taxa as columns
        metadata_df: DataFrame with sample metadata (optional)
        output_file: Output file for plot
        top_n: Number of top taxa to display
        group_column: Column for sample annotation
    """
    logger.info(f"Creating abundance heatmap (top {top_n} taxa)")
    
    # Get top N most abundant taxa
    mean_abundance = abundance_df.mean(axis=0)
    top_taxa = mean_abundance.nlargest(top_n).index.tolist()
    
    # Subset data
    plot_data = abundance_df[top_taxa].copy()
    
    # Log transform (with pseudocount)
    plot_data = np.log10(plot_data + 1)
    
    # Sort samples by group if provided
    if group_column and metadata_df is not None:
        merged = plot_data.merge(
            metadata_df[['sample', group_column]],
            left_index=True,
            right_on='sample',
            how='left'
        )
        merged = merged.sort_values(group_column)
        row_colors = merged[group_column]
        plot_data = merged.drop(columns=['sample', group_column])
    else:
        row_colors = None
    
    # Create clustermap
    fig = sns.clustermap(
        plot_data.T,
        cmap='YlOrRd',
        figsize=(12, 10),
        cbar_kws={'label': 'log10(Abundance + 1)'},
        row_cluster=True,
        col_cluster=True,
        yticklabels=True,
        xticklabels=True
    )
    
    plt.suptitle('Taxonomic Abundance Heatmap', y=0.98)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Heatmap saved to {output_file}")


def plot_differential_abundance(diff_abundance_df, output_file, top_n=20):
    """
    Create volcano plot and barplot for differential abundance
    
    Args:
        diff_abundance_df: DataFrame with differential abundance results
        output_file: Output file for plot
        top_n: Number of top taxa to show in barplot
    """
    logger.info("Creating differential abundance plots")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Volcano plot
    diff_abundance_df['-log10(padj)'] = -np.log10(diff_abundance_df['padj'])
    
    # Define colors based on significance and fold change
    colors = []
    for _, row in diff_abundance_df.iterrows():
        if row['padj'] < 0.05:
            if row['log2_fold_change'] > 0:
                colors.append('red')
            else:
                colors.append('blue')
        else:
            colors.append('gray')
    
    ax1.scatter(
        diff_abundance_df['log2_fold_change'],
        diff_abundance_df['-log10(padj)'],
        c=colors,
        alpha=0.6,
        s=50
    )
    
    ax1.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    ax1.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    ax1.set_xlabel('log2(Fold Change)')
    ax1.set_ylabel('-log10(adjusted p-value)')
    ax1.set_title('Volcano Plot')
    
    # Barplot of top differentially abundant taxa
    sig_taxa = diff_abundance_df[diff_abundance_df['significant']].copy()
    
    if len(sig_taxa) > 0:
        top_sig = sig_taxa.nlargest(top_n, 'log2_fold_change', keep='all')
        
        ax2.barh(
            range(len(top_sig)),
            top_sig['log2_fold_change'],
            color=['red' if x > 0 else 'blue' for x in top_sig['log2_fold_change']]
        )
        ax2.set_yticks(range(len(top_sig)))
        ax2.set_yticklabels(top_sig['taxon'], fontsize=8)
        ax2.set_xlabel('log2(Fold Change)')
        ax2.set_title(f'Top {min(top_n, len(top_sig))} Differentially Abundant Taxa')
        ax2.axvline(x=0, color='black', linestyle='-', alpha=0.3)
    else:
        ax2.text(0.5, 0.5, 'No significant taxa', ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('No Differentially Abundant Taxa')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Differential abundance plots saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Visualization for microbiome data'
    )
    parser.add_argument(
        '-a', '--abundance',
        help='Abundance table file (TSV format)'
    )
    parser.add_argument(
        '-m', '--metadata',
        help='Metadata file (TSV format)'
    )
    parser.add_argument(
        '-o', '--output-dir',
        required=True,
        help='Output directory for plots'
    )
    parser.add_argument(
        '--alpha-diversity',
        help='Alpha diversity file'
    )
    parser.add_argument(
        '--beta-diversity',
        help='Beta diversity distance matrix'
    )
    parser.add_argument(
        '--rarefaction',
        help='Rarefaction curves data'
    )
    parser.add_argument(
        '--differential-abundance',
        help='Differential abundance results'
    )
    parser.add_argument(
        '--group',
        help='Metadata column for grouping'
    )
    parser.add_argument(
        '--top-n',
        type=int,
        default=20,
        help='Number of top taxa to display (default: 20)'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load metadata if provided
    metadata_df = None
    if args.metadata:
        metadata_df = pd.read_csv(args.metadata, sep='\t')
    
    # Taxonomic barplot
    if args.abundance:
        abundance_df = pd.read_csv(args.abundance, sep='\t', index_col=0)
        barplot_file = os.path.join(args.output_dir, 'taxonomic_barplot.png')
        plot_taxonomic_barplot(abundance_df, metadata_df, barplot_file, args.top_n, args.group)
        
        # Heatmap
        heatmap_file = os.path.join(args.output_dir, 'abundance_heatmap.png')
        plot_heatmap(abundance_df, metadata_df, heatmap_file, args.top_n, args.group)
    
    # Alpha diversity plots
    if args.alpha_diversity:
        alpha_df = pd.read_csv(args.alpha_diversity, sep='\t')
        alpha_plot_file = os.path.join(args.output_dir, 'alpha_diversity.png')
        plot_alpha_diversity(alpha_df, metadata_df, alpha_plot_file, args.group)
    
    # Beta diversity plots
    if args.beta_diversity:
        beta_df = pd.read_csv(args.beta_diversity, sep='\t', index_col=0)
        pcoa_file = os.path.join(args.output_dir, 'beta_diversity_pcoa.png')
        plot_pcoa(beta_df, metadata_df, pcoa_file, args.group)
    
    # Rarefaction curves
    if args.rarefaction:
        rare_df = pd.read_csv(args.rarefaction, sep='\t')
        rare_plot_file = os.path.join(args.output_dir, 'rarefaction_curves.png')
        plot_rarefaction_curves(rare_df, rare_plot_file)
    
    # Differential abundance plots
    if args.differential_abundance:
        diff_df = pd.read_csv(args.differential_abundance, sep='\t')
        diff_plot_file = os.path.join(args.output_dir, 'differential_abundance.png')
        plot_differential_abundance(diff_df, diff_plot_file, args.top_n)
    
    logger.info("Visualization completed successfully!")


if __name__ == '__main__':
    main()
