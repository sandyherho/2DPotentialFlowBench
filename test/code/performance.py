#!/usr/bin/env python

"""
Performance Benchmark Script
==========================

Measures runtime and memory usage for Python, Julia, and R implementations
of the 2D Laplace equation solver.

Features:
- Performs warmup runs
- Collects runtime and memory usage data
- Saves results to CSV
- Runs each implementation 200 times
"""

import subprocess
import time
import os
import pandas as pd
import psutil
from datetime import datetime

def create_output_dir():
    """Create output directory if it doesn't exist."""
    os.makedirs("../outputs/data/performance", exist_ok=True)

def get_memory_mb(pid):
    """Get memory usage in MB for a given process ID."""
    try:
        process = psutil.Process(pid)
        mem_info = process.memory_info()
        # Using RSS (Resident Set Size) as memory metric
        return mem_info.rss / 1024 / 1024  # Convert bytes to MB
    except:
        return None

def run_command(cmd, language):
    """Run a single command and measure its performance."""
    start_time = time.time()
    
    # Start the process and capture its output
    process = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             shell=True)
    
    # Get the peak memory usage while process is running
    peak_memory = 0
    while process.poll() is None:
        current_memory = get_memory_mb(process.pid)
        if current_memory is not None and current_memory > peak_memory:
            peak_memory = current_memory
        time.sleep(0.1)  # Small delay to prevent excessive CPU usage
    
    # Get the completion time
    execution_time = time.time() - start_time
    
    # Get the process output
    stdout, stderr = process.communicate()
    
    # Check if the process completed successfully
    if process.returncode != 0:
        print(f"Error running {language} implementation:")
        print(stderr.decode())
        return None, None
    
    return execution_time, peak_memory

def benchmark_implementation(cmd, language, num_runs=10, num_warmup=5):
    """Benchmark a single implementation with warmup runs."""
    print(f"\nBenchmarking {language} implementation...")
    results = []
    
    # Perform warmup runs
    print(f"Performing {num_warmup} warmup runs...")
    for i in range(num_warmup):
        run_command(cmd, language)
    
    # Perform actual benchmark runs
    print(f"Performing {num_runs} benchmark runs...")
    for i in range(num_runs):
        execution_time, peak_memory = run_command(cmd, language)
        if execution_time is not None and peak_memory is not None:
            results.append({
                'language': language,
                'run_number': i + 1,
                'runtime_seconds': execution_time,
                'memory_mb': peak_memory
            })
        
        if (i + 1) % 10 == 0:
            print(f"Completed {i + 1}/{num_runs} runs")
    
    return results

def run_benchmarks():
    """Run benchmarks for all implementations."""
    # Define commands for each implementation
    commands = {
        'Python': 'python sim.py',
        'Julia': 'julia sim.jl',
        'R': 'Rscript sim.R'
    }
    
    # Store all results
    all_results = []
    
    # Run benchmarks for each implementation
    for language, cmd in commands.items():
        results = benchmark_implementation(cmd, language)
        all_results.extend(results)
    
    # Convert results to DataFrame
    df = pd.DataFrame(all_results)
    
    # Add timestamp column
    df['timestamp'] = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    
    # Save results
    output_file = f'../outputs/data/performance/benchmark_results.csv'
    df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")
    
    # Print summary statistics
    print("\nSummary Statistics:")
    summary = df.groupby('language').agg({
        'runtime_seconds': ['mean', 'std', 'min', 'max'],
        'memory_mb': ['mean', 'std', 'min', 'max']
    }).round(3)
    print(summary)

if __name__ == "__main__":
    create_output_dir()
    run_benchmarks()
