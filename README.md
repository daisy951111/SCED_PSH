# SCED-PSH: Variable Selection for Competing Risks

This repository provides MATLAB implementations for variable selection and parameter estimation under the proportional subdistribution hazards (PSH) model using the proposed SCED (Smoothly Clipped Elastic Deviation) penalty.

The proposed method aims to achieve simultaneous estimation and variable selection while approximating the L0 norm via a smooth nonconvex penalty, thereby bridging the parsimonious modeling property of best subset selection and the computational tractability of regularization-based approaches.

---

## 📌 Overview

The repository contains code for:

- Variable selection using the SCED penalty
- Parameter estimation under the PSH model
- Tuning parameter selection
- Simulation studies
- Real data analysis
- Extension to high-dimensional settings via coordinate descent

---

## 📂 File Structure

The main files are organized as follows:

- `Sced_estimation.m`  
  Performs variable selection and compares the proposed SCED method with other methods as mentioned in paper.

- `Sced_lqa.m`  
  Main function for parameter estimation based on the local quadratic approximation (LQA) algorithm.

- `Tunning_sced.m`  
  Implements tuning parameter selection (e.g., EBIC-based selection).

- `Real_data.m`  
  Real data analysis for empirical validation of the proposed method.

- `Sced_cd.m`  
  Extension to high-dimensional settings using coordinate descent for simultaneous variable selection and estimation.

---


