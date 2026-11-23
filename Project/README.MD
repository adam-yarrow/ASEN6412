This folder contains the code that was used to generate all plots and results for my
ASEN6412 final project. At a high level this folder contains a simplified longitudinal
dynamics model of a hypersonic flight vehicle upon which uncertainty quantification and 
parameter estimation is performed.

How to Use:
The project is split into two different phases, forward UQ and inverse UQ.
For forward UQ the following can be executed:
	- RunForwardUQ.m
		This runs a small batch MC to screen for bimodality. Then runs stochastic collocation
		to estimate the mean and variance of the flight vehicle.
	- RunGridDependenceStudy.m
		Runs a small study to ensure the relative error tolerance is acceptable for the 
		RK45 solver.
	- RunMC_SensitivityAnalysis.m
		Runs a large batch MC simulation that is then used to generate a PCE. The PCE is 
		then used to calculate total sensitivity SOBOL indices.
	
For inverse UQ the following can be executed:
	- RunInverseUq_DiscreteBayes.m
		Estimates both constant and linear Cm_alpha using a discrete bayes filter approach.
	- RunInverseUQ_NLS.m
		Estimates a polynomial model of Cm_alpha using non-linear least squares.
	- RunNLS_IC_SensitivityStudy.m
		Used to check the sensitivity of the NLS solution to the initial guess for the Cm_alpha
		model.	

Other files of relevance:
	- ModelParams.m:
		This file contains the baseline model parameters and is the master source of truth.
	- ForwardUqPerturbationVars.m:
		This file defines the distributions for uncertain parameters.
	- Startup.m
		This file initializes dependencies. It is ran by default in all scripts.

Dependencies:
	- MATLAB Toolboxes
		- Parallel Computing
		- Symbolic Math
		- Statistics and Machine Learning
	- ../Library 
		- KLE solver
		- parfor waitbar (MATLAB file exchange)
		- truncated random normal generator (MATLAB file exchange) NO LONGER USED
		- Sparse Smolyak grid generator
	
Folders Overview:
	-GaussianProcess
		KLE function
	- Helpers
		struct updater
	- Images
		Figures that were saved for report generation
	- InverseUQ
		Contains discrete bayes filter and NLS implementation
	- Models
		Contains the dynamics model of the flight vehicle
	- MonteCarlo
		Contains MC related functionality
	- Plotting
		Plotting scripts
	- PolyChaos
		PCE basis function generation and suppport
	- Statistics	
		truncation function for a normal pdf
	- StochasticCollocation
		Same as MC but for SC 
	- Testing
		Testing scripts, many deprecated.