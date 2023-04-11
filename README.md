# HOMEDF
Higher-Order Maximum Entropy Distribution Functions Solver

## Overview

HOMEDF is an optimization algorithm that solves the entropy maximization problem. Through 
a Newton search and some numerical integration, it finds the free coefficients of the 
distribution function associated with a given moment state. 

## Main files examples

Two main files are provided as examples. The first one finds the distribution functions for 
a given moment state for the three, five, and seven moment distributions. It outputs the 
closing flux, the eigenvalues of the corresponding flux Jacobian and writes the distribution 
function to a .dat file so that it can be plotted. The second main file solves the five moment 
distributions for a range of moment states and then calculates the closing flux. The closing 
flux is then output to VTK so that it can be visualized.

Note that the search algorithm can only handle non-dimensional moment states. Member functions 
are provided in order to easily non-dimensionalize the state.

## Build

HOMEDF can be built using cmake in a build directory. The required libraries are Eigen, VTK, 
gtest, MPI and OpenMP.
