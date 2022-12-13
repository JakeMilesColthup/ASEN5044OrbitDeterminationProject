# ASEN5044OrbitDeterminationProject

Main.m is the main script to run the entire project for everything past part 1. MainPart1.m is for running part 1 and generating those results.

checkVisibleStations performs the angle checks

computeLinearizedDyn calculates F, G, and Omega

computeLinearizedH calculates H

computeYNL generates nonlinear measurement data using lowercase h

nonLinearOrbit is the function used by ode45

All of the "reformat" functions are ways to rearrange the different storing method of measurement data

The LKF, EKF, and UKF each have a runCombined function to run them with simulated data and withData function to run with the data log

simulateNoisyMeas generates simulated noisy measurements

simulateTruthModel generates truth model data with process noise
