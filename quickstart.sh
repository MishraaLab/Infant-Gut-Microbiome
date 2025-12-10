#!/bin/bash
# Quick start script for long-read microbiome analysis

echo "====================================================="
echo "Infant Gut Microbiome Analysis - Quick Start"
echo "====================================================="
echo ""

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "ERROR: Python 3 is not installed. Please install Python 3.7 or higher."
    exit 1
fi

echo "✓ Python 3 found: $(python3 --version)"

# Check if virtual environment exists
if [ ! -d "microbiome_env" ]; then
    echo ""
    echo "Creating virtual environment..."
    python3 -m venv microbiome_env
    echo "✓ Virtual environment created"
fi

# Activate virtual environment
echo ""
echo "Activating virtual environment..."
source microbiome_env/bin/activate

# Install requirements
echo ""
echo "Installing Python dependencies..."
pip install -q --upgrade pip
pip install -q -r requirements.txt
echo "✓ Dependencies installed"

# Check for external tools
echo ""
echo "Checking for external tools..."

check_tool() {
    if command -v $1 &> /dev/null; then
        echo "  ✓ $1 found"
        return 0
    else
        echo "  ✗ $1 not found (optional)"
        return 1
    fi
}

check_tool "kraken2"
check_tool "bracken"
check_tool "minimap2"
check_tool "NanoPlot"

# Create config file if it doesn't exist
if [ ! -f "config.ini" ]; then
    echo ""
    echo "Creating default configuration file..."
    python3 scripts/pipeline.py --create-config config.ini
    echo "✓ Configuration file created: config.ini"
    echo ""
    echo "⚠ Please edit config.ini with your data paths and database locations."
fi

# Show next steps
echo ""
echo "====================================================="
echo "Setup Complete!"
echo "====================================================="
echo ""
echo "Next steps:"
echo "1. Place your FASTQ files in data/raw_data/"
echo "2. Edit config.ini with your settings"
echo "3. Run the pipeline:"
echo "   python3 scripts/pipeline.py -c config.ini"
echo ""
echo "For help with individual scripts:"
echo "   python3 scripts/01_quality_control.py --help"
echo "   python3 scripts/02_taxonomic_classification.py --help"
echo "   python3 scripts/03_diversity_analysis.py --help"
echo "   python3 scripts/04_statistical_analysis.py --help"
echo "   python3 scripts/05_visualization.py --help"
echo ""
echo "For detailed documentation, see:"
echo "   - README.md"
echo "   - docs/GUIDE.md"
echo ""
echo "====================================================="
