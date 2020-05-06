# Life History and the Evolution of Social Learning

This repository contains the scripts to reproduce all results and plots in 

***Deffner D, McElreath R. 2020. The importance of life history and population regulation for the evolution of social learning.
Phil. Trans. R. Soc. B 20190492. http://dx.doi.org/10.1098/rstb.2019.0492***

Description of files:

- "FigS1Recursions.r" contains numerical simulation of analytical models and plotting script to reproduce population growth curves in Fig. S1 in the ESM

- "FigS2Isoclines.r" calculates and plots isoclines for equilibrium population sizes (Fig. S2)

- "Simulation_Code_single_FigS3.r". Code for running a single simulation with specified parameter combination. Also contains code to create basic adaptation and demographic dynamics plot from simulation output (Fig. S3)

- "Simulation_Code_Cluster.r" provides the code for the main simulation models for all 810 parameter combinations used in the paper

- "Fig2MainResults.r" lets you compute relevant summary statistics from output of "Simulation_Code_Cluster.r" and plot main results for different costs of individual learning (Figs. 2, S4 and S5)

- "Fig3TimingeffRates.r" takes individual-level data from "Simulation_Code_Cluster.r" and computes adaptation levels, number of births as well as effective vital rates for different times after change in the environment (Fig. 3)
