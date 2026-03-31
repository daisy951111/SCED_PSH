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

## 📂 File

The main files are organized as follows:

- `sced_lqa.m`  
  Main function for simultaneous variable selection and estimation based on the local quadratic approximation (LQA) algorithm.

- `tunning_sced.m`  
  Implements tuning parameter selection (e.g., EBIC-based selection).

- `real_data.m`  
  Real data analysis for empirical validation of the proposed method.

- `sced_cd.m`  
  Extension to high-dimensional settings using coordinate descent for simultaneous variable selection and estimation.

- `computing_time.m`  
  Conducts Monte Carlo simulations to evaluate the computational efficiency of the proposed SCED method and competing approaches.

---

## 📊 Data
Simulation data are generated within the code
Real data example:
GSE5479 dataset
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5479

---

## 🔖 Citation

If you use this code, please cite:  Variable Selection for Competing Risks Survival Data via the SCED Penalty. (Under review)

---

## 📬 Contact

For questions or comments, please contact:  daisy1111@csu.edu.cn
