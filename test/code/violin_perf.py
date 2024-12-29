#!/usr/bin/env python

"""
Symmetric Significance Matrix Visualization
=======================================

Creates clean, black and white significance matrices showing symmetric relationships.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os

def create_output_dir():
    """Create output directory if it doesn't exist."""
    os.makedirs("../outputs/figs/performance", exist_ok=True)

def setup_plot_style():
    """Set up clean, minimalist plot style."""
    plt.style.use('default')
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'figure.facecolor': 'white'
    })

def compute_significance_matrix(df, metric, languages):
    """Compute significance matrix with Bonferroni correction."""
    n = len(languages)
    sig_matrix = np.zeros((n, n))
    n_comparisons = (n * (n-1)) // 2  # number of unique comparisons
    alpha = 0.05
    
    for i, lang1 in enumerate(languages):
        for j, lang2 in enumerate(languages[i+1:], i+1):
            data1 = df[df['language'] == lang1][metric]
            data2 = df[df['language'] == lang2][metric]
            
            # Mann-Whitney U test with Bonferroni correction
            _, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
            adjusted_p = min(p_value * n_comparisons, 1.0)
            
            # If significant, mark both symmetric positions
            if adjusted_p < alpha:
                sig_matrix[i, j] = 1
                sig_matrix[j, i] = 1  # symmetric position
    
    return sig_matrix

def create_significance_plot(sig_matrix, languages, title, filename):
    """Create a clean black and white significance matrix plot."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot the matrix
    im = ax.imshow(sig_matrix, cmap='binary', vmin=0, vmax=1)
    
    # Add annotations
    for i in range(len(languages)):
        for j in range(len(languages)):
            if i != j:  # Don't show anything on diagonal
                if sig_matrix[i, j] == 1:
                    ax.text(j, i, "SIG", ha='center', va='center', 
                           color='white', fontweight='bold')
    
    # Customize plot
    ax.set_xticks(np.arange(len(languages)))
    ax.set_yticks(np.arange(len(languages)))
    ax.set_xticklabels(languages)
    ax.set_yticklabels(languages)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    
    # Add title
    #ax.set_title(f"{title}\n(Bonferroni-corrected, α=0.05)")
    
    # Add grid
    ax.set_xticks(np.arange(-.5, len(languages), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(languages), 1), minor=True)
    ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def main():
    # Read the data
    df = pd.read_csv("../outputs/data/performance/benchmark_results.csv")
    languages = sorted(df['language'].unique())
    
    # Set up plot style
    setup_plot_style()
    
    # Create runtime significance matrix
    runtime_sig = compute_significance_matrix(df, 'runtime_seconds', languages)
    create_significance_plot(
        runtime_sig,
        languages,
        'Runtime Significance Matrix',
        '../outputs/figs/performance/runtime_significance.png'
    )
    
    # Create memory significance matrix
    memory_sig = compute_significance_matrix(df, 'memory_mb', languages)
    create_significance_plot(
        memory_sig,
        languages,
        'Memory Usage Significance Matrix',
        '../outputs/figs/performance/memory_significance.png'
    )
    
    print("Significance matrices created! Files saved in ../outputs/figs/performance/")

if __name__ == "__main__":
    create_output_dir()
    main()
