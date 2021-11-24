Functional data analysis applied to vertical ground reaction forces obtained from countermovement jumps.

This code supports my draft paper entitled 'The Effects of Curve Registration on Vertical Ground Reaction Force curves of the Countermovement Jump'.  

The code generates 384 linear models that predict jump height, peak power or classify the jumps into those performed with or without arm swing. Each model is based on different data preprocessing operations for:
- curve length standardisation (padding or time normalisation)
- landmark registration (16 combinations of four landmarks, including no registration)
- continuous registration (either applied or not)
- feature extraction (Functional PCA or Analysis of Characteristiing Phases)

The MATLAB library for  functional data analysis can be found here: https://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/

The analysis is run by jumpsFDA.m, loading raw data (file available on request).

The code main loop iterates through the following operations:
- data partitioning
- length standardisation
- functional smoothing
- landmark registration
- continuous registration
- decomposition analysis
- linear modelling 

Other files are dedicated to generating plots of the results.
