#!/usr/bin/env Rscript

#' 2D Laplace Equation Numerical Simulation
#' ======================================
#'
#' This script computes the numerical solution to the 2D Laplace equation with
#' rectangular obstacles. It implements a finite difference method with iterative
#' solving and handles boundary conditions including farfield and obstacle boundaries.
#'
#' Features:
#' - Finite difference solution of 2D Laplace equation
#' - Multiple rectangular obstacle handling
#' - Automatic data saving to ../outputs/data/
#' - Error tracking for convergence monitoring
#'
#' Author: Sandy H. S. Herho <sandy.herho@email.ucr.edu>
#' Date: 12/28/2024
#' License: WTFPL
#' Version: 1.0.0

library(tidyverse)

#' Set up the geometric parameters and grid
#' @return List containing geometry parameters and grids
setup_geometry <- function() {
  # Domain dimensions
  Lx <- 20
  Ly <- 11
  nx <- 201
  ny <- 111
  dx <- Lx/(nx-1)
  dy <- Ly/(ny-1)
  
  # Create coordinate grid using expand.grid for R equivalent of meshgrid
  x_seq <- seq(0, Lx, length.out = nx)
  y_seq <- seq(0, Ly, length.out = ny)
  grid <- expand.grid(x = x_seq, y = y_seq)
  X <- matrix(grid$x, nrow = ny, ncol = nx)
  Y <- matrix(grid$y, nrow = ny, ncol = nx)
  
  list(
    Lx = Lx, Ly = Ly,
    nx = nx, ny = ny,
    dx = dx, dy = dy,
    X = X, Y = Y
  )
}

#' Define the positions of rectangular obstacles
#' @return List containing obstacle parameters
setup_obstacles <- function() {
  list(
    box_imin = c(30, 60),
    box_imax = c(50, 80),
    box_jmin = c(60, 90),
    box_jmax = c(80, 120),
    box_colors = c("gray", "gray")
  )
}

#' Apply boundary conditions including obstacles and farfield
#' @param phi Matrix of potential values
#' @param phi_old Previous iteration's potential values
#' @param obstacles List of obstacle parameters
#' @param dx Grid spacing in x direction
#' @param Uinf Free stream velocity
#' @return Updated phi matrix
apply_boundary_conditions <- function(phi, phi_old, obstacles, dx, Uinf) {
  # Farfield boundaries
  phi[1:(nrow(phi)-1), 1] <- phi[1:(nrow(phi)-1), 2] - Uinf*dx
  phi[1:(nrow(phi)-1), ncol(phi)] <- phi[1:(nrow(phi)-1), ncol(phi)-1] + Uinf*dx
  phi[1, 1:(ncol(phi)-1)] <- phi[2, 1:(ncol(phi)-1)]
  phi[nrow(phi), 1:(ncol(phi)-1)] <- phi[nrow(phi)-1, 1:(ncol(phi)-1)]
  
  # Obstacles are impenetrable
  for (bc in 1:length(obstacles$box_imin)) {
    i_range <- obstacles$box_imin[bc]:obstacles$box_imax[bc]
    j_range <- (obstacles$box_jmin[bc]+1):(obstacles$box_jmax[bc]-1)
    
    # Apply boundary conditions for each side of obstacle
    phi[i_range, obstacles$box_jmin[bc]] <- 
      phi[i_range, obstacles$box_jmin[bc]-1]
    phi[i_range, obstacles$box_jmax[bc]] <- 
      phi[i_range, obstacles$box_jmax[bc]+1]
    phi[obstacles$box_imin[bc], j_range] <- 
      phi[obstacles$box_imin[bc]-1, j_range]
    phi[obstacles$box_imax[bc], j_range] <- 
      phi[obstacles$box_imax[bc]+1, j_range]
  }
  
  phi
}

#' Compute velocities and pressure coefficient from potential
#' @param phi Matrix of potential values
#' @param dx Grid spacing in x direction
#' @param dy Grid spacing in y direction
#' @param Uinf Free stream velocity
#' @return List containing U, V velocities and pressure coefficient
compute_velocities_and_pressure <- function(phi, dx, dy, Uinf) {
  ny <- nrow(phi)
  nx <- ncol(phi)
  
  U <- matrix(0, nrow = ny, ncol = nx)
  V <- matrix(0, nrow = ny, ncol = nx)
  
  # Calculate velocities
  U[2:(ny-1), 2:(nx-1)] <- (phi[2:(ny-1), 2:(nx-1)] - phi[2:(ny-1), 1:(nx-2)])/dx
  V[2:(ny-1), 2:(nx-1)] <- (phi[2:(ny-1), 2:(nx-1)] - phi[1:(ny-2), 2:(nx-1)])/dy
  
  # Calculate pressure coefficient
  vel_sq <- U[2:(ny-1), 2:(nx-1)]^2 + V[2:(ny-1), 2:(nx-1)]^2
  Cp <- 1 - vel_sq/(Uinf^2)
  
  list(U = U, V = V, Cp = Cp)
}

#' Save matrix to CSV with R suffix
#' @param mat Matrix to save
#' @param filename Base filename
save_matrix_to_csv <- function(mat, filename) {
  # Extract basename without extension and add _R suffix
  base_name <- tools::file_path_sans_ext(basename(filename))
  new_filename <- file.path("../outputs/data", paste0(base_name, "_R.csv"))
  write.csv(mat, new_filename, row.names = FALSE)
}

#' Save obstacles information to CSV
#' @param obstacles List of obstacle parameters
save_obstacles_to_csv <- function(obstacles) {
  # Create data frame from obstacles list
  df <- data.frame(
    parameter = c("box_imin", "box_imax", "box_jmin", "box_jmax", "box_colors"),
    value_1 = c(obstacles$box_imin[1], obstacles$box_imax[1], 
                obstacles$box_jmin[1], obstacles$box_jmax[1], 
                obstacles$box_colors[1]),
    value_2 = c(obstacles$box_imin[2], obstacles$box_imax[2], 
                obstacles$box_jmin[2], obstacles$box_jmax[2], 
                obstacles$box_colors[2])
  )
  write.csv(df, "../outputs/data/obstacles_R.csv", row.names = FALSE)
}

#' Run the main simulation and save results
#' @param numSteps Number of time steps
#' @param Uinf Free stream velocity
run_simulation <- function(numSteps = 40000, Uinf = 1) {
  # Setup
  geom <- setup_geometry()
  obstacles <- setup_obstacles()
  
  # Initialize solution arrays
  phi <- matrix(0, nrow = geom$ny, ncol = geom$nx)
  err <- matrix(0, nrow = numSteps, ncol = 1)
  
  # Create output directory if it doesn't exist
  dir.create("../outputs/data", recursive = TRUE, showWarnings = FALSE)
  
  # Time stepping loop
  for (t in 1:numSteps) {
    phi_old <- phi
    
    # Solve Laplace equation
    phi[2:(geom$ny-1), 2:(geom$nx-1)] <- 0.25 * (
      phi_old[1:(geom$ny-2), 2:(geom$nx-1)] + 
      phi_old[3:geom$ny, 2:(geom$nx-1)] +
      phi_old[2:(geom$ny-1), 1:(geom$nx-2)] + 
      phi_old[2:(geom$ny-1), 3:geom$nx]
    )
    
    # Apply boundary conditions
    phi <- apply_boundary_conditions(phi, phi_old, obstacles, geom$dx, Uinf)
    
    # Calculate error using Frobenius norm
    err[t] <- norm(phi - phi_old, type = "F")
  }
  
  # Compute final velocities and pressure
  results <- compute_velocities_and_pressure(phi, geom$dx, geom$dy, Uinf)
  
  # Save results to CSV files
  save_matrix_to_csv(geom$X, "X.csv")
  save_matrix_to_csv(geom$Y, "Y.csv")
  save_matrix_to_csv(phi, "phi.csv")
  save_matrix_to_csv(results$U, "U.csv")
  save_matrix_to_csv(results$V, "V.csv")
  save_matrix_to_csv(results$Cp, "Cp.csv")
  save_matrix_to_csv(err, "err.csv")
  save_obstacles_to_csv(obstacles)
  
  cat("Simulation completed and results saved to CSV files!\n")
}

# Main execution
if (!interactive()) {
  run_simulation()
}
