# Clustered Hidden Markov Models (CHMMs)
Datasets, codes, and numerical results for reproduce the results in "Clustered hidden Markov models", J. Bok & S. Shin.

## Environment Details
- **Operating System**: Ubuntu 22.04.2 LTS
- **Kernel Version**: 5.19.0-46-generic (x86_64)
- **Programming Language**: R 4.4.1

## Repository Structure

| Directory | Description |
|-----------|-------------|
| **`data/`** | Contains all datasets used in **Simulation 1**, **Simulation 2**, and the **protein-structure analysis**. A file `data/README.md` explains how each dataset was obtained. |
| **`code/`** | R scripts implementing every algorithm evaluated in the paper, including our proposed **EM-ADMM** procedure. Three turnkey scripts—`Reproduce_Simulation1.R`, `Reproduce_Simulation2.R`, and `Reproduce_ProteinStructureAnalysis.R`—recreate all numerical results and figures reported in the manuscript. |
| **`results/`** | Outputs produced by running the reproduction scripts above. |

### Reproduction Workflow

1. **Restore the R environment**
   ```r
   install.packages("renv")  # if not already installed
   renv::restore()

2. **Run the reproduction scripts for numerical results**

   Run `code/Reproduce_Simulation1.R`, `code/Reproduce_Simulation2.R`, and `code/ProteinStructureAnalysis.R`.

3. **Run the reproduction scripts for figures in manuscripts**

   Run `Code_for_figures.R`.

   **Important**: Reproducing **Figure 4** requires a brief manual action. Please follow the instructions provided in the *####Code generating Figure 4####* section of 'Code_for_figures.R' before running that part of the script.

   
