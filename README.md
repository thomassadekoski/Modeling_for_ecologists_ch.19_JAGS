Applied Statistical Modeling for Ecologists (Chapter 19)
Overview
This repository contains R code implementing occupancy models described in Chapter 19 of Applied Statistical Modeling for Ecologists. The scripts demonstrate how imperfect detection influences species distribution inference and how hierarchical occupancy models correct this bias using both maximum likelihood and Bayesian approaches.
This workflow includes:
Simulation of true occupancy and detection processes
Naive logistic regression ignoring detection error
Likelihood-based occupancy modeling using unmarked
Bayesian occupancy modeling using JAGS
Posterior prediction and visualization using ggplot2
The example uses a simulated environmental covariate (humidity index) that influences both occupancy probability and detection probability.
