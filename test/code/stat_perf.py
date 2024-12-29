#!/usr/bin/env python

"""
Statistical Analysis of Performance Data
=====================================

Performs robust statistical analysis on runtime and memory usage data
across different language implementations.

Features:
- Non-parametric statistical tests
- Effect size calculations
- Bootstrapped confidence intervals
- Exports results to CSV
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import mannwhitneyu
import os
from itertools import combinations

def create_output_dir():
    """Create output directory if it doesn't exist."""
    os.makedirs("../outputs/data/performance", exist_ok=True)

def bootstrap_mean_ci(data, n_bootstrap=10000, ci=0.95):
    """Calculate bootstrapped confidence intervals for mean."""
    bootstrap_means = []
    for _ in range(n_bootstrap):
        sample = np.random.choice(data, size=len(data), replace=True)
        bootstrap_means.append(np.mean(sample))
    
    ci_lower = np.percentile(bootstrap_means, (1 - ci) * 100 / 2)
    ci_upper = np.percentile(bootstrap_means, (1 + ci) * 100 / 2)
    return ci_lower, ci_upper

def cliffs_delta(x, y):
    """Calculate Cliff's Delta effect size."""
    nx, ny = len(x), len(y)
    dominance = 0
    
    for i in x:
        for j in y:
            if i > j:
                dominance += 1
            elif i < j:
                dominance -= 1
                
    return dominance / (nx * ny)

def analyze_metric(df, metric_name):
    """Perform statistical analysis on a given metric with Bonferroni correction."""
    languages = df['language'].unique()
    results = []
    
    # Descriptive statistics with bootstrapped CIs
    for lang in languages:
        data = df[df['language'] == lang][metric_name]
        ci_lower, ci_upper = bootstrap_mean_ci(data)
        
        results.append({
            'Comparison': f'{lang} (Descriptive)',
            'Metric': metric_name,
            'Mean': float(np.round(np.mean(data), 3)),
            'Median': float(np.round(np.median(data), 3)),
            'Std': float(np.round(np.std(data), 3)),
            'CI_Lower': float(np.round(ci_lower, 3)),
            'CI_Upper': float(np.round(ci_upper, 3)),
            'IQR': float(np.round(stats.iqr(data), 3))
        })
    
    # Pairwise comparisons with Bonferroni correction
    num_comparisons = len(list(combinations(languages, 2)))
    alpha = 0.05  # Original alpha level
    
    for lang1, lang2 in combinations(languages, 2):
        data1 = df[df['language'] == lang1][metric_name]
        data2 = df[df['language'] == lang2][metric_name]
        
        # Mann-Whitney U test
        statistic, pvalue = mannwhitneyu(data1, data2, alternative='two-sided')
        
        # Apply Bonferroni correction
        adjusted_pvalue = min(pvalue * num_comparisons, 1.0)
        
        # Cliff's Delta effect size
        effect_size = cliffs_delta(data1, data2)
        
        results.append({
            'Comparison': f'{lang1} vs {lang2}',
            'Metric': metric_name,
            'Test': 'Mann-Whitney U',
            'Statistic': float(np.round(statistic, 3)),
            'P_Value': float(np.round(pvalue, 3)),
            'Adjusted_P_Value': float(np.round(adjusted_pvalue, 3)),
            'Effect_Size': float(np.round(effect_size, 3)),
            'Effect_Size_Type': "Cliff's Delta",
            'Significant': adjusted_pvalue < alpha
        })
    
    return pd.DataFrame(results)

def main():
    # Read the data
    df = pd.read_csv("../outputs/data/performance/benchmark_results.csv")
    
    # Analyze runtime
    runtime_results = analyze_metric(df, 'runtime_seconds')
    
    # Analyze memory
    memory_results = analyze_metric(df, 'memory_mb')
    
    # Combine results
    all_results = pd.concat([runtime_results, memory_results], ignore_index=True)
    
    # Save results
    output_file = "../outputs/data/performance/statistical_analysis.csv"
    all_results.to_csv(output_file, index=False, float_format='%.3f')
    print(f"Statistical analysis results saved to: {output_file}")

if __name__ == "__main__":
    create_output_dir()
    main()
