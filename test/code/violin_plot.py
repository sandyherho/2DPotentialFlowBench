#!/usr/bin/env python
"""
Generate violin plots for runtime and memory performance comparison.

Author: Sandy H. S. Herho <sandy.herho@email.ucr.edu>
License: WTFPL
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

def load_performance_data(file_path: str) -> pd.DataFrame:
    """
    Load performance data from CSV file.
    
    Args:
        file_path (str): Path to the CSV file
        
    Returns:
        pd.DataFrame: Loaded performance data
    """
    return pd.read_csv(file_path)

def create_violin_plot(data: pd.DataFrame, 
                      metric: str,
                      output_path: str,
                      y_label: str) -> None:
    """
    Create and save violin plot for the specified metric.
    
    Args:
        data (pd.DataFrame): Performance data
        metric (str): Column name for the metric to plot
        output_path (str): Path to save the plot
        y_label (str): Label for y-axis
    """
    plt.figure(figsize=(10, 8))
    
    # Set style
    sns.set_style("whitegrid")
    sns.set_palette("deep")
    
    # Create violin plot
    ax = sns.violinplot(data=data,
                       x='language',
                       y=metric,
                       cut=0,
                       inner='box')
    
    # Customize plot
    plt.ylabel(y_label, fontsize=20)
    plt.xlabel("Language", fontsize=20)
    #plt.title(f'Distribution of {y_label} by Language')
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Save plot
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    """Main execution function."""
    # Define paths
    input_path = Path("../outputs/data/performance/benchmark_results.csv")
    output_dir = Path("../outputs/figs/performance")
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    data = load_performance_data(input_path)
    
    # Create plots
    create_violin_plot(
        data=data,
        metric='runtime_seconds',
        output_path=str(output_dir / 'runtime_violin.png'),
        y_label='Execution Time (seconds)'
    )
    
    create_violin_plot(
        data=data,
        metric='memory_mb',
        output_path=str(output_dir / 'memory_violin.png'),
        y_label='Memory Usage (MB)'
    )

if __name__ == "__main__":
    main()
