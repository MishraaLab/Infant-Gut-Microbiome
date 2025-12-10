#!/usr/bin/env python3
"""
Main Pipeline for Long-Read Microbiome Data Analysis
Integrates all analysis steps: QC, classification, diversity, statistics, and visualization
"""

import os
import sys
import argparse
import logging
import subprocess
import configparser
from pathlib import Path
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class MicrobiomePipeline:
    """Main pipeline class for microbiome analysis"""
    
    def __init__(self, config_file):
        """
        Initialize pipeline with configuration
        
        Args:
            config_file: Path to configuration file
        """
        self.config = self.load_config(config_file)
        self.scripts_dir = Path(__file__).parent
        self.start_time = datetime.now()
        
        # Create output directories
        self.create_output_dirs()
    
    def load_config(self, config_file):
        """Load configuration from file"""
        logger.info(f"Loading configuration from {config_file}")
        
        config = configparser.ConfigParser()
        config.read(config_file)
        
        return config
    
    def create_output_dirs(self):
        """Create output directory structure"""
        output_dir = self.config['Paths']['output_dir']
        
        subdirs = [
            'qc_results',
            'classification',
            'diversity',
            'statistics',
            'visualization',
            'logs'
        ]
        
        for subdir in subdirs:
            os.makedirs(os.path.join(output_dir, subdir), exist_ok=True)
    
    def run_command(self, cmd, step_name):
        """
        Run a command and log output
        
        Args:
            cmd: Command to run (list)
            step_name: Name of the pipeline step
        """
        logger.info(f"Running {step_name}...")
        logger.info(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            
            logger.info(f"{step_name} completed successfully")
            
            # Log output
            log_file = os.path.join(
                self.config['Paths']['output_dir'],
                'logs',
                f'{step_name}.log'
            )
            with open(log_file, 'w') as f:
                f.write(f"=== {step_name} Log ===\n\n")
                f.write("STDOUT:\n")
                f.write(result.stdout)
                f.write("\n\nSTDERR:\n")
                f.write(result.stderr)
            
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"{step_name} failed!")
            logger.error(f"Error: {e.stderr}")
            
            # Log error
            log_file = os.path.join(
                self.config['Paths']['output_dir'],
                'logs',
                f'{step_name}_error.log'
            )
            with open(log_file, 'w') as f:
                f.write(f"=== {step_name} Error Log ===\n\n")
                f.write("STDOUT:\n")
                f.write(e.stdout)
                f.write("\n\nSTDERR:\n")
                f.write(e.stderr)
            
            return False
    
    def step_quality_control(self):
        """Run quality control step"""
        if not self.config['Pipeline'].getboolean('run_qc', True):
            logger.info("Skipping quality control (disabled in config)")
            return True
        
        input_file = self.config['Paths']['input_fastq']
        output_dir = os.path.join(self.config['Paths']['output_dir'], 'qc_results')
        output_file = os.path.join(output_dir, 'filtered_reads.fastq')
        
        cmd = [
            'python3',
            str(self.scripts_dir / '01_quality_control.py'),
            '-i', input_file,
            '-o', output_file,
            '--qc-dir', output_dir,
            '--min-length', self.config['QC']['min_length'],
            '--min-quality', self.config['QC']['min_quality']
        ]
        
        if self.config['QC'].getboolean('skip_nanoplot', False):
            cmd.append('--skip-nanoplot')
        
        success = self.run_command(cmd, 'quality_control')
        
        if success:
            self.filtered_reads = output_file
        
        return success
    
    def step_taxonomic_classification(self):
        """Run taxonomic classification step"""
        if not self.config['Pipeline'].getboolean('run_classification', True):
            logger.info("Skipping taxonomic classification (disabled in config)")
            return True
        
        input_file = self.filtered_reads
        output_dir = os.path.join(self.config['Paths']['output_dir'], 'classification')
        
        cmd = [
            'python3',
            str(self.scripts_dir / '02_taxonomic_classification.py'),
            '-i', input_file,
            '-o', output_dir,
            '-d', self.config['Classification']['database'],
            '--method', self.config['Classification']['method'],
            '--threads', self.config['Classification']['threads']
        ]
        
        if self.config['Classification'].getboolean('run_bracken', False):
            cmd.append('--run-bracken')
            cmd.extend(['--bracken-level', self.config['Classification']['bracken_level']])
        
        success = self.run_command(cmd, 'taxonomic_classification')
        
        if success:
            self.classification_report = os.path.join(output_dir, 'kraken2_report.txt')
        
        return success
    
    def step_diversity_analysis(self):
        """Run diversity analysis step"""
        if not self.config['Pipeline'].getboolean('run_diversity', True):
            logger.info("Skipping diversity analysis (disabled in config)")
            return True
        
        input_file = self.classification_report
        output_dir = os.path.join(self.config['Paths']['output_dir'], 'diversity')
        
        cmd = [
            'python3',
            str(self.scripts_dir / '03_diversity_analysis.py'),
            '-i', input_file,
            '-o', output_dir,
            '--format', 'kraken',
            '--alpha'
        ]
        
        if self.config['Diversity'].getboolean('calculate_beta', False):
            cmd.append('--beta')
            cmd.extend(['--beta-metric', self.config['Diversity']['beta_metric']])
        
        if self.config['Diversity'].getboolean('rarefaction', False):
            cmd.append('--rarefaction')
        
        success = self.run_command(cmd, 'diversity_analysis')
        
        if success:
            self.alpha_diversity = os.path.join(output_dir, 'alpha_diversity.tsv')
            self.beta_diversity = os.path.join(output_dir, 'beta_diversity_bray_curtis.tsv')
        
        return success
    
    def step_statistical_analysis(self):
        """Run statistical analysis step"""
        if not self.config['Pipeline'].getboolean('run_statistics', False):
            logger.info("Skipping statistical analysis (disabled in config)")
            return True
        
        # Check if metadata is provided
        if 'metadata' not in self.config['Paths']:
            logger.warning("No metadata file provided. Skipping statistical analysis.")
            return True
        
        abundance_file = self.config['Paths'].get('abundance_table')
        metadata_file = self.config['Paths']['metadata']
        output_dir = os.path.join(self.config['Paths']['output_dir'], 'statistics')
        
        if not abundance_file:
            logger.warning("No abundance table provided. Skipping statistical analysis.")
            return True
        
        cmd = [
            'python3',
            str(self.scripts_dir / '04_statistical_analysis.py'),
            '-a', abundance_file,
            '-m', metadata_file,
            '-o', output_dir,
            '-g', self.config['Statistics']['group_column'],
            '--normalize', self.config['Statistics']['normalization']
        ]
        
        if hasattr(self, 'alpha_diversity'):
            cmd.extend(['--alpha-diversity', self.alpha_diversity])
        
        if hasattr(self, 'beta_diversity'):
            cmd.extend(['--beta-diversity', self.beta_diversity])
        
        success = self.run_command(cmd, 'statistical_analysis')
        
        if success:
            self.diff_abundance = os.path.join(output_dir, 'differential_abundance.tsv')
        
        return success
    
    def step_visualization(self):
        """Run visualization step"""
        if not self.config['Pipeline'].getboolean('run_visualization', True):
            logger.info("Skipping visualization (disabled in config)")
            return True
        
        output_dir = os.path.join(self.config['Paths']['output_dir'], 'visualization')
        
        cmd = [
            'python3',
            str(self.scripts_dir / '05_visualization.py'),
            '-o', output_dir
        ]
        
        # Add optional inputs
        if 'abundance_table' in self.config['Paths']:
            cmd.extend(['-a', self.config['Paths']['abundance_table']])
        
        if 'metadata' in self.config['Paths']:
            cmd.extend(['-m', self.config['Paths']['metadata']])
            cmd.extend(['--group', self.config['Statistics']['group_column']])
        
        if hasattr(self, 'alpha_diversity'):
            cmd.extend(['--alpha-diversity', self.alpha_diversity])
        
        if hasattr(self, 'beta_diversity'):
            cmd.extend(['--beta-diversity', self.beta_diversity])
        
        if hasattr(self, 'diff_abundance'):
            cmd.extend(['--differential-abundance', self.diff_abundance])
        
        success = self.run_command(cmd, 'visualization')
        
        return success
    
    def run(self):
        """Run the complete pipeline"""
        logger.info("=" * 60)
        logger.info("Starting Microbiome Analysis Pipeline")
        logger.info("=" * 60)
        
        steps = [
            ('Quality Control', self.step_quality_control),
            ('Taxonomic Classification', self.step_taxonomic_classification),
            ('Diversity Analysis', self.step_diversity_analysis),
            ('Statistical Analysis', self.step_statistical_analysis),
            ('Visualization', self.step_visualization)
        ]
        
        failed_steps = []
        
        for step_name, step_func in steps:
            logger.info(f"\n{'=' * 60}")
            logger.info(f"Step: {step_name}")
            logger.info(f"{'=' * 60}")
            
            try:
                success = step_func()
                if not success:
                    failed_steps.append(step_name)
                    if self.config['Pipeline'].getboolean('stop_on_error', True):
                        logger.error(f"Stopping pipeline due to error in {step_name}")
                        break
            except Exception as e:
                logger.error(f"Exception in {step_name}: {str(e)}")
                failed_steps.append(step_name)
                if self.config['Pipeline'].getboolean('stop_on_error', True):
                    break
        
        # Summary
        elapsed_time = datetime.now() - self.start_time
        
        logger.info("\n" + "=" * 60)
        logger.info("Pipeline Summary")
        logger.info("=" * 60)
        logger.info(f"Total runtime: {elapsed_time}")
        
        if failed_steps:
            logger.warning(f"Failed steps: {', '.join(failed_steps)}")
        else:
            logger.info("All steps completed successfully!")
        
        logger.info(f"Results saved in: {self.config['Paths']['output_dir']}")


def create_default_config(output_file):
    """Create a default configuration file"""
    config = configparser.ConfigParser()
    
    config['Paths'] = {
        'input_fastq': 'data/raw_data/sample.fastq',
        'output_dir': 'data/results',
        'abundance_table': '',  # Optional
        'metadata': ''  # Optional
    }
    
    config['Pipeline'] = {
        'run_qc': 'True',
        'run_classification': 'True',
        'run_diversity': 'True',
        'run_statistics': 'False',
        'run_visualization': 'True',
        'stop_on_error': 'True'
    }
    
    config['QC'] = {
        'min_length': '1000',
        'min_quality': '10',
        'skip_nanoplot': 'False'
    }
    
    config['Classification'] = {
        'database': '/path/to/kraken2/database',
        'method': 'kraken2',
        'threads': '4',
        'run_bracken': 'False',
        'bracken_level': 'S'
    }
    
    config['Diversity'] = {
        'calculate_beta': 'False',
        'beta_metric': 'bray_curtis',
        'rarefaction': 'False'
    }
    
    config['Statistics'] = {
        'group_column': 'treatment',
        'normalization': 'total_sum'
    }
    
    with open(output_file, 'w') as f:
        config.write(f)
    
    logger.info(f"Default configuration file created: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Long-read microbiome analysis pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Create default config
  python pipeline.py --create-config config.ini
  
  # Run pipeline
  python pipeline.py -c config.ini
        """
    )
    parser.add_argument(
        '-c', '--config',
        help='Configuration file'
    )
    parser.add_argument(
        '--create-config',
        help='Create a default configuration file'
    )
    
    args = parser.parse_args()
    
    if args.create_config:
        create_default_config(args.create_config)
        return
    
    if not args.config:
        parser.print_help()
        sys.exit(1)
    
    # Run pipeline
    pipeline = MicrobiomePipeline(args.config)
    pipeline.run()


if __name__ == '__main__':
    main()
