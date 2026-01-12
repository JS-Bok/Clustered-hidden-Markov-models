#---------------------------------#
# README: Supplementary R Scripts #
#---------------------------------#

This folder contains R scripts used in 'Clustered Hidden Markov Models*' (J. Bok & S. Shin). The detailed information (inputs and outputs) of the main functions is described in each script.

---

## 1. Auxiliary Functions
- File: 'Auxiliary_Functions.R'
- Description:  
  This file contains auxiliary functions for the main functions (e.g., Forward-Backward Algorithm, BIC score, etc.).

---

## 2. Main Functions

### 2-1. EM Algorithm for MLE of CHMM
- File: 'EM_MLE.R'
- Description:  
  This file contains a function that performs an EM algorithm for the MLE of a CHMM ('EM_MLE').

### 2-2. EM Algorithm for MPLE of CHMM
- File: 'EM_CHMM.R'
- Description:  
  This file contains a function that performs an EM algorithm for the MPLE of a CHMM.  
  - 'EM_MPLE_singleLambda': EM algorithm with a single penalty parameter.  
  - 'EM_MPLE': EM algorithm with multiple penalty parameters.  

### 2-3. EM Algorithm for OE of Stationary HMM
- File: 'EM_Oracle.R'
- Description:  
  This file contains a function that performs an EM algorithm for the Oracle Estimation (OE) of a stationary HMM.  
  - 'EM_Oracle_Sim{i}_A{j}': Oracle estimation function for Simulation '{i}', with the '{j}'-th transition matrix.

### 2-4. EM Algorithm for MLE and OE of FHMM
- File: 'EM_FHMM.R'
- Description:  
  This file contains a function that performs an EM algorithm for the OE of an FHMM.  
  - 'fhmm': Function for OE of FHMM.    

---

## 3. Scripts for Reproducing Numerical Results

### 3-1. Reproducing Simulation 1
- File: 'Reproduce_Simulation1.R'
- Description:  
  This script reproduces the numerical results of Simulation 1 using parallel computing.

### 3-2. Reproducing Simulation 2
- File: 'Reproduce_Simulation2.R'
- Description:  
  This script reproduces the numerical results of Simulation 2 using parallel computing.

### 3-3. Reproducing Protein Structure Analysis
- File: 'Reproduce_ProteinStructureAnalysis.R'
- Description:  
  This script reproduces the numerical results of the Protein Structure Analysis.

---

## 4. Scripts for Reproducing Figures in Protein structure analysis
- File: 'Code_for_figures.R'
- Description:  
  This script generates Figures in Protein structure analysis.  
- Prerequisite: Output of 'Reproduce_ProteinStructureAnalysis.R'.
