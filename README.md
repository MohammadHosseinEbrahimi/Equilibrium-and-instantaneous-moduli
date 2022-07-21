# Equilibrium-and-instantaneous-moduli

This script is to analyze EQUILIBRIUM and instantaneous moduli based on multistep stress-relaxation protocol

Script loads the 4 step biomechanical testing data, finds the equilibrium (full relaxation) points, fits a line into the data and calculates the
Youngs modulus for the samples. it also calucaltes instantaneous moduli from each loading step, fits a line and
extracts initial and strain-dependent instantaneous moduli.

some variables like relaxation time, strain, strain rate, the unit of force and displacement, 
number of relaxation points etc, must be modified according to the measurement protocol

This code also includes Correction of modulus based on Hayes correction factor.

 By Mohammadhossein Ebrahimi 14.2.2018
