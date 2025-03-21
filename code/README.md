## Supplementary R Scripts

This folder contains R scripts used in the "Clustered Hidden Markov Models" (J. Bok & S. Shin). Each script provides detailed information on its main functions, including input/output parameters.

### 1. Auxiliary Functions
- **File:** `Auxiliary_Functions.R`  
- **Description:**  
  Contains auxiliary functions that support the main routines (e.g., Forward-Backward Algorithm, BIC score calculation, etc.).

---

### 2. Main Functions

#### 2-1. EM Algorithm for MLE of CHMM
- **File:** `EM_MLE.R`  
- **Description:**  
  Implements an EM algorithm for the Maximum Likelihood Estimation (MLE) of CHMM via the `EM_MLE` function.

#### 2-2. EM Algorithm for MPLE of CHMM
- **File:** `EM_CHMM.R`  
- **Description:**  
  Implements an EM algorithm for the Maximum Penalized Likelihood Estimation (MPLE) of CHMM.  
  - `EM_MPLE_singleLambda`: EM algorithm using a single penalty parameter.  
  - `EM_MPLE`: EM algorithm using multiple penalty parameters.

#### 2-3. EM Algorithm for OE of Stationary HMM
- **File:** `EM_Oracle.R`  
- **Description:**  
  Implements an EM algorithm for Oracle Estimation (OE) of a stationary Hidden Markov Model (HMM).  
  - `EM_Oracle_Sim{i}_A{j}`: Oracle estimation function for Simulation `{i}`, using the `{j}`-th transition matrix.

#### 2-4. EM Algorithm for MLE and OE of FHMM
- **File:** `EM_FHMM.R`  
- **Description:**  
  Implements an EM algorithm for both the MLE and OE of an FHMM.  
  - `fhmm`: Function for OE of FHMM.

---

### 3. Scripts for Reproducing Numerical Results

All datasets are stored [here](https://github.com/JS-Bok/Clustered-hidden-Markov-models/tree/main/data) and are automatically loaded by each script.

#### 3-1. Reproducing Simulation 1
- **File:** `Reproduce_Simulation1.R`  
- **Description:**  
  Reproduces the numerical results of Simulation 1 using parallel computing.

#### 3-2. Reproducing Simulation 2
- **File:** `Reproduce_Simulation2.R`  
- **Description:**  
  Reproduces the numerical results of Simulation 2 using parallel computing.

#### 3-3. Reproducing Protein Structure Analysis
- **File:** `Reproduce_ProteinStructureAnalysis.R`  
- **Description:**  
  Reproduces the numerical results of the Protein Structure Analysis.

---

### 4. Scripts for Reproducing Figures in Protein Structure Analysis
- **File:** `Code_for_figures.R`  
- **Description:**  
  Generates the figures for the Protein Structure Analysis.  
- **Prerequisite:** Output from `Reproduce_ProteinStructureAnalysis.R` is required.
