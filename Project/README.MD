# ASEN6412 Final Project

This folder contains the code used to generate all plots and results for my ASEN6412 final project. At a high level, it contains a simplified longitudinal dynamics model of a hypersonic flight vehicle, upon which uncertainty quantification (UQ) and parameter estimation are performed.

## How to Use

The project is split into two phases: forward UQ and inverse UQ.

### Forward UQ

| Script | Description |
|---|---|
| `RunForwardUQ.m` | Runs a small batch MC to screen for bimodality, then runs stochastic collocation to estimate the mean and variance of the flight vehicle. |
| `RunGridDependenceStudy.m` | Runs a small study to ensure the relative error tolerance is acceptable for the RK45 solver. |
| `RunMC_SensitivityAnalysis.m` | Runs a large batch MC simulation that is then used to generate a PCE. The PCE is then used to calculate total sensitivity Sobol indices. |

### Inverse UQ

| Script | Description |
|---|---|
| `RunInverseUq_DiscreteBayes.m` | Estimates both constant and linear `Cm_alpha` using a discrete Bayes filter approach. |
| `RunInverseUQ_NLS.m` | Estimates a polynomial model of `Cm_alpha` using non-linear least squares. |
| `RunNLS_IC_SensitivityStudy.m` | Checks the sensitivity of the NLS solution to the initial guess for the `Cm_alpha` model. |

## Other Files of Relevance

- `ModelParams.m`: Contains the baseline model parameters and is the master source of truth.
- `ForwardUqPerturbationVars.m`: Defines the distributions for uncertain parameters.
- `Startup.m`: Initializes dependencies. It is run by default in all scripts.

## Dependencies

### MATLAB Toolboxes

- Parallel Computing
- Symbolic Math
- Statistics and Machine Learning

### `../Library`

- KLE solver
- parfor waitbar (MATLAB File Exchange)
- truncated random normal generator (MATLAB File Exchange) — **no longer used**
- Sparse Smolyak grid generator

## Folders Overview

| Folder | Contents |
|---|---|
| `GaussianProcess` | KLE function |
| `Helpers` | Struct updater |
| `Images` | Figures saved for report generation |
| `InverseUQ` | Discrete Bayes filter and NLS implementation |
| `Models` | Dynamics model of the flight vehicle |
| `MonteCarlo` | MC-related functionality |
| `Plotting` | Plotting scripts |
| `PolyChaos` | PCE basis function generation and support |
| `Statistics` | Truncation function for a normal PDF |
| `StochasticCollocation` | Same as MonteCarlo, but for SC |
| `Testing` | Testing scripts, many deprecated |
