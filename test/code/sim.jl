#!/usr/bin/env julia

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

Author: Sandy H. S. Herho <sandy.herho@email.ucr.edu>
Date: 12/28/2024
License: WTFPL
Version: 1.0.0
"""

using LinearAlgebra
using DataFrames
using CSV

"""
    setup_geometry()

Set up the geometric parameters and grid.
Returns named tuple containing geometry parameters and grids.
"""
function setup_geometry()
    # Domain dimensions
    Lx, Ly = 20, 11
    nx, ny = 201, 111
    dx = Lx/(nx-1)
    dy = Ly/(ny-1)
    
    # Create coordinate grid (Julia is column-major, transpose for consistency)
    x = range(0, Lx, length=nx)
    y = range(0, Ly, length=ny)
    X = [x[i] for j in 1:ny, i in 1:nx]
    Y = [y[j] for j in 1:ny, i in 1:nx]
    
    return (Lx=Lx, Ly=Ly, nx=nx, ny=ny, dx=dx, dy=dy, X=X, Y=Y)
end

"""
    setup_obstacles()

Define the positions of rectangular obstacles.
Returns named tuple containing obstacle parameters.
"""
function setup_obstacles()
    return (
        box_imin = [30, 60],
        box_imax = [50, 80],
        box_jmin = [60, 90],
        box_jmax = [80, 120],
        box_colors = ["gray", "gray"]
    )
end

"""
    apply_boundary_conditions!(phi, phi_old, obstacles, dx, Uinf)

Apply boundary conditions including obstacles and farfield.
Modifies phi in-place.
"""
function apply_boundary_conditions!(phi, phi_old, obstacles, dx, Uinf)
    # Farfield boundaries
    phi[1:end-1,1] .= @view(phi[1:end-1,2]) .- Uinf*dx
    phi[1:end-1,end] .= @view(phi[1:end-1,end-1]) .+ Uinf*dx
    phi[1,1:end-1] .= @view(phi[2,1:end-1])
    phi[end,1:end-1] .= @view(phi[end-1,1:end-1])
    
    # Obstacles are impenetrable
    for bc in 1:length(obstacles.box_imin)
        i_range = obstacles.box_imin[bc]:obstacles.box_imax[bc]
        j_range = obstacles.box_jmin[bc]:obstacles.box_jmax[bc]
        
        # Apply boundary conditions for each side of obstacle
        phi[i_range, obstacles.box_jmin[bc]] .= 
            @view(phi[i_range, obstacles.box_jmin[bc]-1])
        phi[i_range, obstacles.box_jmax[bc]] .= 
            @view(phi[i_range, obstacles.box_jmax[bc]+1])
        phi[obstacles.box_imin[bc], j_range[2:end-1]] .= 
            @view(phi[obstacles.box_imin[bc]-1, j_range[2:end-1]])
        phi[obstacles.box_imax[bc], j_range[2:end-1]] .= 
            @view(phi[obstacles.box_imax[bc]+1, j_range[2:end-1]])
    end
    return phi
end

"""
    compute_velocities_and_pressure(phi, dx, dy, Uinf)

Compute velocities and pressure coefficient from potential.
Returns named tuple containing U, V velocities and pressure coefficient.
"""
function compute_velocities_and_pressure(phi, dx, dy, Uinf)
    ny, nx = size(phi)
    U = zeros(ny, nx)
    V = zeros(ny, nx)
    
    # Calculate velocities
    U[2:end-1,2:end-1] .= (@view(phi[2:end-1,2:end-1]) .- @view(phi[2:end-1,1:end-2])) ./ dx
    V[2:end-1,2:end-1] .= (@view(phi[2:end-1,2:end-1]) .- @view(phi[1:end-2,2:end-1])) ./ dy
    
    # Calculate pressure coefficient
    vel_sq = @view(U[2:end-1,2:end-1]).^2 .+ @view(V[2:end-1,2:end-1]).^2
    Cp = 1 .- vel_sq ./ (Uinf^2)
    
    return (U=U, V=V, Cp=Cp)
end

"""
    save_2d_array_to_csv(array, filename)

Save array to CSV with Julia suffix.
Handles both 1D and 2D arrays appropriately.
"""
function save_2d_array_to_csv(array, filename)
    # Extract basename and add _jl suffix
    base_name = split(basename(filename), ".")[1]
    new_filename = joinpath("../outputs/data", base_name * "_jl.csv")
    
    # Handle different array dimensions
    if ndims(array) == 1
        # For 1D array, create single-column DataFrame
        df = DataFrame(value = array)
    else
        # For 2D array, create DataFrame with numbered columns
        df = DataFrame(array, [Symbol("X$i") for i in 1:size(array, 2)])
    end
    
    CSV.write(new_filename, df)
end

"""
    save_obstacles_to_csv(obstacles)

Save obstacles information to CSV.
"""
function save_obstacles_to_csv(obstacles)
    df = DataFrame(
        parameter = ["box_imin", "box_imax", "box_jmin", "box_jmax", "box_colors"],
        value_1 = [obstacles.box_imin[1], obstacles.box_imax[1], 
                  obstacles.box_jmin[1], obstacles.box_jmax[1], 
                  obstacles.box_colors[1]],
        value_2 = [obstacles.box_imin[2], obstacles.box_imax[2], 
                  obstacles.box_jmin[2], obstacles.box_jmax[2], 
                  obstacles.box_colors[2]]
    )
    CSV.write("../outputs/data/obstacles_jl.csv", df)
end

"""
    run_simulation(numSteps=40000, Uinf=1)

Run the main simulation and save results.
"""
function run_simulation(numSteps=40000, Uinf=1)
    # Setup
    geom = setup_geometry()
    obstacles = setup_obstacles()
    
    # Initialize solution arrays
    phi = zeros(geom.ny, geom.nx)
    err = zeros(numSteps)
    
    # Create output directory if it doesn't exist
    mkpath("../outputs/data")
    
    # Time stepping loop
    for t in 1:numSteps
        phi_old = copy(phi)
        
        # Solve Laplace equation
        phi[2:end-1,2:end-1] .= 0.25 .* (
            @view(phi_old[1:end-2,2:end-1]) .+ 
            @view(phi_old[3:end,2:end-1]) .+
            @view(phi_old[2:end-1,1:end-2]) .+ 
            @view(phi_old[2:end-1,3:end])
        )
        
        # Apply boundary conditions
        apply_boundary_conditions!(phi, phi_old, obstacles, geom.dx, Uinf)
        
        # Calculate error
        err[t] = norm(phi .- phi_old)
    end
    
    # Compute final velocities and pressure
    results = compute_velocities_and_pressure(phi, geom.dx, geom.dy, Uinf)
    
    # Save results to CSV files
    save_2d_array_to_csv(geom.X, "X.csv")
    save_2d_array_to_csv(geom.Y, "Y.csv")
    save_2d_array_to_csv(phi, "phi.csv")
    save_2d_array_to_csv(results.U, "U.csv")
    save_2d_array_to_csv(results.V, "V.csv")
    save_2d_array_to_csv(results.Cp, "Cp.csv")
    save_2d_array_to_csv(err, "err.csv")
    save_obstacles_to_csv(obstacles)
    
    println("Simulation completed and results saved to CSV files!")
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    run_simulation()
end
