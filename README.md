Functional data analysis applied to vertical ground reaction forces obtained from countermovement jumps.

This code supports my draft paper entitled 'The Effects of Curve Registration on Vertical Ground Reaction Force curves of the Countermovement Jump'.  

The code generates 1024 linear models that predict jump height, peak power or classify the jumps into those performed with or without arm swing. Each model is based on different data preprocessing operations for:
- curve length standardisation (padding or time normalisation)
- landmark registration (16 combinations of four landmarks, including no registration)
- continuous registration (either applied or not)
- feature extraction (Functional PCA or Analysis of Characteristiing Phases)
- component rotation (unrotated/varimax)

The MATLAB library for  functional data analysis can be found here: https://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/

The analysis is run by jumpsFDA.m. The main loop iterates through the following operations:
- data partitioning
- length standardisation
- functional smoothing
- landmark registration
- continuous registration
- decomposition analysis
- linear modelling 

jumpsFDA.m runs the full analysis, using a raw datafile that is not available yet - see below. 

jumpsFDA.m outputs a Matlab data file containing:
- decomp (2x2x16x2 cell array) amplitude/phase decomposition  
- fdPar (2x2 cell array) FDA smoothing parameters, dataset by length standardisation
- isValid (2x2x16 cell array) logical array indicating valid registrations
- models (2x2x16x2 cell array) structures holding the fitted models
- name (2x2x16x2 string array) names of the datasets
- results (2x2x16x2 cell array) tables holding model component scores and model outputs (originally used in SAS models)
- vgrfACP (2x2x16x2 cell array) structures holding the ACP unrotated and varimax scores
- vgrfFd (2x2x16x2 cell array) functional data objects for the smoothed VGRF curves
- vgrfPCA (2x2x16x2 cell array) structures holding the functional PCA, unrotated and varimax
- warpFd (2x2x16x2 cell array) functional data objects for the time-warp curves

The 4D arrays has this structure
- 2 datasets:CMJ_NA & CMJ_A 
- 2 length standardisations: PAD & LTN
- 16 landmark registration: 0000-1111
- 2 continuous registrations: - & C


jumpsFDA.m reads the raw datafile, 'compactjumpsdata.mat', which contains anonymised data. It is currently unavailable until ethics approval is granted. It has the following structure:
- bwall: 64x16 double (bodyweights)
- grf.raw: 64x16 cell (raw 1D VGRF data of variable lengths)
- grf.initiation: 64x16 double (jump initiaion index)
- grf.takeoff: 64x16 doubel (take-off index)
- jumpOrder: 104x16 (master list specifying the the jump type for each trial)
- nJumpsPerSubject: 64x1 double (how many jumps for each subject)
- sDataID: 1x64 double (lookup table linking array index to subject ID)
- sJumpID: 1x104 double (lookup table linking trial index to jump order ID)


Other files are dedicated to generating plots of the results.
