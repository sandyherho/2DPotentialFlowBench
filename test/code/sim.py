#!/usr/bin/env python

"""
2D Laplace Equation Numerical Simulation
=======================================

This script computes the numerical solution to the 2D Laplace equation with
rectangular obstacles. It implements a finite difference method with iterative
solving and handles boundary conditions including farfield and obstacle boundaries.

Features:
- Finite difference solution of 2D Laplace equation
- Multiple rectangular obstacle handling
- Automatic data saving to ../outputs/data/
- Error tracking for convergence monitoring

Dependencies:
- NumPy: For numerical computations
- OS: For directory handling

Author: Sandy H. S. Herho <sandy.herho@email.ucr.edu>
Date: 12/28/2024
License: WTFPL
Version: 1.0.0
"""

import numpy as np
import os
import pandas as pd

def setup_geometry():
    """Set up the geometric parameters and grid"""
    # Domain dimensions
    Lx, Ly = 20, 11
    nx, ny = 201, 111
    dx = Lx/(nx-1)
    dy = Ly/(ny-1)
    
    # Create coordinate grid
    X, Y = np.meshgrid(np.linspace(0,Lx,nx), np.linspace(0,Ly,ny))
    
    return Lx, Ly, nx, ny, dx, dy, X, Y

def setup_obstacles():
    """Define the positions of rectangular obstacles"""
    return {
        'box_imin': [30, 60],
        'box_imax': [50, 80],
        'box_jmin': [60, 90],
        'box_jmax': [80, 120],
        'box_colors': ['gray', 'gray']
    }

def apply_boundary_conditions(phi, phi_old, obstacles, dx, Uinf):
    """Apply boundary conditions including obstacles and farfield"""
    # Farfield boundaries
    phi[0:-1,0] = phi[0:-1,1] - Uinf*dx
    phi[0:-1,-1] = phi[0:-1,-2] + Uinf*dx
    phi[0,0:-1] = phi[1,0:-1]
    phi[-1,0:-1] = phi[-2,0:-1]
    
    # Obstacles are impenetrable
    for bc in range(len(obstacles['box_imin'])):
        phi[obstacles['box_imin'][bc]:obstacles['box_imax'][bc],
            obstacles['box_jmin'][bc]] = phi[obstacles['box_imin'][bc]:obstacles['box_imax'][bc],
                                           obstacles['box_jmin'][bc]-1]
        phi[obstacles['box_imin'][bc]:obstacles['box_imax'][bc],
            obstacles['box_jmax'][bc]] = phi[obstacles['box_imin'][bc]:obstacles['box_imax'][bc],
                                           obstacles['box_jmax'][bc]+1]
        phi[obstacles['box_imin'][bc],
            obstacles['box_jmin'][bc]+1:obstacles['box_jmax'][bc]-1] = phi[obstacles['box_imin'][bc]-1,
                                                                         obstacles['box_jmin'][bc]+1:obstacles['box_jmax'][bc]-1]
        phi[obstacles['box_imax'][bc],
            obstacles['box_jmin'][bc]+1:obstacles['box_jmax'][bc]-1] = phi[obstacles['box_imax'][bc]+1,
                                                                         obstacles['box_jmin'][bc]+1:obstacles['box_jmax'][bc]-1]
    return phi

def compute_velocities_and_pressure(phi, dx, dy, Uinf):
    """Compute velocities and pressure coefficient from potential"""
    ny, nx = phi.shape
    U = np.zeros([ny,nx])
    V = np.zeros([ny,nx])
    
    # Calculate velocities
    U[1:-1,1:-1] = (phi[1:-1,1:-1] - phi[1:-1,0:-2])/dx
    V[1:-1,1:-1] = (phi[1:-1,1:-1] - phi[0:-2,1:-1])/dy
    
    # Calculate pressure coefficient
    vel_sq = (U[1:-1,1:-1]*U[1:-1,1:-1] + V[1:-1,1:-1]*V[1:-1,1:-1])
    Cp = 1 - vel_sq/(Uinf*Uinf)
    
    return U, V, Cp

def save_2d_array_to_csv(array, filename):
    """Save 2D numpy array to CSV with proper indexing"""
    df = pd.DataFrame(array)
    df.to_csv(filename, index=False)

def save_obstacles_to_csv(obstacles):
    """Save obstacles information to CSV"""
    # Convert each obstacle parameter to a DataFrame row
    data = {
        'parameter': ['box_imin', 'box_imax', 'box_jmin', 'box_jmax', 'box_colors'],
        'value_1': [obstacles['box_imin'][0], obstacles['box_imax'][0], 
                   obstacles['box_jmin'][0], obstacles['box_jmax'][0], 
                   obstacles['box_colors'][0]],
        'value_2': [obstacles['box_imin'][1], obstacles['box_imax'][1], 
                   obstacles['box_jmin'][1], obstacles['box_jmax'][1], 
                   obstacles['box_colors'][1]]
    }
    df = pd.DataFrame(data)
    df.to_csv('../outputs/data/obstacles.csv', index=False)

def run_simulation(numSteps=40000, Uinf=1):
    """Run the main simulation and save results"""
    # Setup
    Lx, Ly, nx, ny, dx, dy, X, Y = setup_geometry()
    obstacles = setup_obstacles()
    
    # Initialize solution arrays
    phi = np.zeros([ny,nx])
    err = np.zeros([numSteps,1])
    
    # Create output directory if it doesn't exist
    os.makedirs('../outputs/data', exist_ok=True)
    
    # Time stepping loop
    for t in range(numSteps):
        phi_old = np.copy(phi)
        
        # Solve Laplace equation
        phi[1:-1,1:-1] = 0.25*(
            phi_old[0:-2,1:-1] + phi_old[2:,1:-1] + 
            phi_old[1:-1,0:-2] + phi_old[1:-1,2:])
        
        # Apply boundary conditions
        phi = apply_boundary_conditions(phi, phi_old, obstacles, dx, Uinf)
        
        # Calculate error
        err[t] = np.linalg.norm(phi-phi_old)
    
    # Compute final velocities and pressure
    U, V, Cp = compute_velocities_and_pressure(phi, dx, dy, Uinf)
    
    # Save results to CSV files
    save_2d_array_to_csv(X, '../outputs/data/X.csv')
    save_2d_array_to_csv(Y, '../outputs/data/Y.csv')
    save_2d_array_to_csv(phi, '../outputs/data/phi.csv')
    save_2d_array_to_csv(U, '../outputs/data/U.csv')
    save_2d_array_to_csv(V, '../outputs/data/V.csv')
    save_2d_array_to_csv(Cp, '../outputs/data/Cp.csv')
    save_2d_array_to_csv(err, '../outputs/data/err.csv')
    save_obstacles_to_csv(obstacles)
    
    print("Simulation completed and results saved to CSV files!")

if __name__ == "__main__":
    run_simulation()
