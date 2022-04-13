% ************************************************************************
% Script: jumpsFDA
% Purpose:  Perform batches of functional data analysis on the jumps time series
%
%
% ************************************************************************

clear;

% ************************************************************************
%   Read data file
% ************************************************************************

% get the repository's code path
path = fileparts( which('jumpsFDA.m') );
path = [path '/../'];
cd(path);

% read the processed data file
load( 'data/compactjumpdata.mat' );


% ************************************************************************
%   Constants
% ************************************************************************

g = 9.812; % acceleration due to gravity

tFreq = 1; % sampling frequency

nSets = 2; % number of data sets (1 = CMJ(No arms); 2 = CMJ(All))
nStd = 2; % number of methods to standardise length
nModels = 4; % number of models
nLMReg = 16; % number of landmark registration combinations including none
nCTReg = 2; % number of continuous registrations (applied or none)

% presets based on not exceeding jump performance error threshold
preset.nBasis = [ 130 105; 130 105 ];
preset.lambda = [ 1E0 1E0; 1E0 1E0 ];

% ************************************************************************
%   Command switches
% ************************************************************************

options.doInitialFit = false; % whether to fit with max flexibility first

options.doCheckFit = false; % whether to check the goodness of fit
options.doResidualCheck = false; % for this check, whether to check the residuals
options.doSortCheck = false; % for this check, whether to sort from largest error to smallest
options.doUseWtFunc = false; % whether to use the weighted functional data object
options.doPlotDerivatives = false; % whether to plot the function and its first two derivatives

options.reg.doCurvePlots = false; % flag whether to show registration curve plots
options.reg.doRegHistogram = false; % flag whether to show the spread of landmark points
options.reg.doRegPlots = false; % flag whether to show registration curve plots

options.reg.doRemoveFaulty = true; % whether to remove faulty registrations

% ************************************************************************
%   Baseline settings
% ************************************************************************

setup.data.tFreq = 1; % time intervals per second
setup.data.initial = 1; % initial padding value
setup.data.threshold1 = 0.08; % fraction of BW for first detection
setup.data.threshold2 = 0.025; % fraction of BW for sustained low threshold
setup.data.prctileLimit = 90; % no outliers are beyond this limit

setup.Fd.basisOrder = 4; % 5th order for a basis expansion of quartic splines
setup.Fd.penaltyOrder = 2; % roughness penalty
setup.Fd.names = [{'Time (ms)'},{'Jumps'},{'GRF (BW)'}]; % axes names
setup.Fd.tolerance = 0.001; % performance measure error tolerance

setup.reg.nBasis = 13; % numbers of bases for registration
setup.reg.basisOrder = 3; % time warping basis order for registration
setup.reg.wLambda = 1E0; % roughness penalty for time warp 1E-2
setup.reg.XLambda = 1E3; % roughness penalty to prevent wiggles in y
setup.reg.convCriterion = 0.001; % smallest change in C 
setup.reg.maxIterations = 4; % maximum iterations to prevent infinite loop
setup.reg.nBasisTotalWarp = 203; % for re-smoothing total warp
setup.reg.allMustBeValid = false; % all curves must be valid for registration acceptance

setup.reg.lm.grfmin = false; % use VGRF minimum as a landmark?
setup.reg.lm.pwrmin = false; % use Power minimum as a landmark?
setup.reg.lm.pwrcross = false; % use Power crossing point as a landmark?
setup.reg.lm.pwrmax = false; % use Power maximum as a landmark?

setup.reg.faultCriterion = 'RelativeArea'; % after vs before area ratio
setup.reg.faultZScore = 3.5; % fault threshold

setup.pca.nComp = 15; % number of PCA components to be retained
setup.pca.nCompWarp = 5; % number of PCA components to be retained

setup.models.nRepeats = 4; % number of repetitions of CV
setup.models.nFolds = 5; % number of CV folds for each repetition
setup.models.seed = 12345; % random seed for reproducibility
setup.models.spec = 'linear'; % type of GLM
setup.models.upper = 'linear'; % linear model without interactions
setup.models.interactions = false; % interactions between ampl and warp
setup.models.criterion = 'bic'; % predictor selection criterion
setup.models.RSqMeritThreshold = 0.7; % merit threshold for stepwise selection

setup.filename = 'results/jumpsAnalysis4.mat'; % where to save the analysis


% ************************************************************************
%   Extract data
% ************************************************************************

% exclude jumps from subjects in the second data collection
subjectExclusions = find( ismember( sDataID, ...
            [ 14, 39, 68, 86, 87, 11, 22, 28, 40, 43, 82, 88, 95, 97, ...
              100, 121, 156, 163, 196 ] ) );

% exclude specific jumps with excessive double movements
jumpExclusions = [  6104, 6114, 6116, ...
                    9404, 9411, 9413, 9416, ...
                    0101, 0103, 0106, 0111, 0114 ];

[ rawData, refSet, typeSet ] =  extractVGRFData( ... 
                                    grf, bwall, nJumpsPerSubject, ...
                                    sDataID, sJumpID, jumpOrder, ...
                                    subjectExclusions, jumpExclusions, ...
                                    setup.data );


% ************************************************************************
%   Determine smoothing levels
% ************************************************************************

perf = cell( nSets );
tSpan = cell( nSets, nStd );
vgrfData = cell( nSets, nStd );

for i = 1:nSets

    seriesLengths = cellfun( @length, rawData{i} );
    padLen = max( seriesLengths );
    normLen = fix(median( seriesLengths ));
   
    for j = 1:nStd

       switch j
           case 1              
               % truncate pre-padded series to this length
               vgrfData{i,j} = padData( rawData{i}, ...
                                        padLen, ...
                                        setup.data.initial );              
               
           case 2              
               % truncate times series and time normalise
               vgrfData{i,j} = timeNormData( rawData{i}, ...
                                         normLen );
                                                                         
       end
       
       % compute jump performances from the truncated raw data
       if j == 1
           perf{i} = jumpperf( vgrfData{i,j} );
       end
       
       % store time span
       fixLen = size( vgrfData{i,j}, 1 );
       tSpan{i,j} = -fixLen+1:0;
       
       if options.doInitialFit
            % one basis function per data point
            setup.Fd.nBasis = fixLen + setup.Fd.basisOrder + 2;
            disp(['Smoothing Validation: ' num2str([i j])]);
            disp(['# Points = ' num2str(fixLen)]);

            while setup.Fd.nBasis > 0
                setup.Fd.nBasis = input('# Basis Functions = ');
    
                if setup.Fd.nBasis > 0
                    validateSmoothing(  vgrfData{i,j}, ...
                                        tSpan{i,j}, ...
                                        setup.Fd, ...
                                        perf{i} );
                end
            end
       end
       
    end
    
end



% ************************************************************************
%   Begin functional data analysis
% ************************************************************************

vgrfFd = cell( nSets, nStd, nLMReg, nCTReg );
warpFd = cell( nSets, nStd, nLMReg, nCTReg );
fdPar = cell( nSets, nStd );

regIter = zeros( nSets, nStd, nLMReg, nCTReg );
decomp = cell( nSets, nStd, nLMReg, nCTReg );
isValid = cell( nSets, nStd, nLMReg, nCTReg );

name = strings( nSets, nStd, nLMReg, nCTReg );

vgrfPCA = cell( nSets, nStd, nLMReg, nCTReg );
vgrfACP = cell( nSets, nStd, nLMReg, nCTReg );

results = cell( nSets, nStd, nLMReg, nCTReg );

models = cell( nSets, nStd, nLMReg, nCTReg );

%load( setup.filename );
%load( 'results/jumpsAnalysis3.mat' );
%vgrfFd{1,1,16,1} = [];
%warpFd{1,1,16,1} = [];

% set random seed for reproducibility
rng( setup.models.seed );

for i = 1:nSets
    
    % setup partitioning for all models using this data set   
    partitions = kFoldSubjectCV(  refSet{i}(:,1), ...
                                  setup.models.nRepeats, ...
                                  setup.models.nFolds );
    
    for j = 1:nStd

       % use presets specific to padding or time normalisation
       fdSetup = setup.Fd;
       fdSetup.nBasis = preset.nBasis(i,j);
       fdSetup.lambda = preset.lambda(i,j);

       % generate the smooth curves
       if isempty( vgrfFd{i,j,1,1} ) || isempty( fdPar{i,j} )
           [ vgrfFd{i,j,1,1}, fdPar{i,j} ] = smoothVGRF(  ...
                                                vgrfData{i,j}, ...
                                                tSpan{i,j}, ...
                                                fdSetup, ...
                                                options );
       end      
       
       for k = 1:nLMReg
           
           for l = 1:nCTReg
               
               % encode the processing procedure
               [ name{i,j,k,l}, setup.reg.lm ] = encodeProc( i, j, k, l );
               disp(['*************** ' name{i,j,k,l} ' ***************']);
               
               % setup reference data in case rows have to be removed
               if l==1 
                   % create a baseline for the first registration
                   baselineFd = vgrfFd{i,j,1,1};
                   ref = refSet{i};
                   type = typeSet{i};
                   jperf = perf{i};
                   part = partitions;
               end
               
               if l == 1 % first 'l' loop
                   if k > 1 && isempty( vgrfFd{i,j,k,l} )
                       % landmark registration required
                       % applied to unregistered curves
                       [ vgrfFd{i,j,k,l}, warpFd{i,j,k,l}, ...
                           regIter(i,j,k,l), isValidReg ] = registerVGRF( ...
                                         tSpan{i,j}, ...
                                         vgrfFd{i,j,1,1}, ...
                                         'Landmark', ...
                                         setup.reg );                                    
                   else
                       isValidReg = true( size(type) );
                   end
                   
               else % second 'l' loop
                   if isempty( vgrfFd{i,j,k,l} )
                       % continuous registration required
                       % applied to prior-registered curves
                       [ vgrfFd{i,j,k,l}, warpFd{i,j,k,l}, ...
                           regIter(i,j,k,l), isValidReg ] = registerVGRF( ...
                                        tSpan{i,j}, ...
                                        vgrfFd{i,j,k,1}, ...
                                        'Continuous', ...
                                        setup.reg, ...
                                        warpFd{i,j,k,1} );
                   end
               end
               
               isValid{i,j,k,l} = isValidReg;
               if ~all( isValidReg )
                   % remove faulty registrations
                   vgrfFd{i,j,k,l} = selectFd( vgrfFd{i,j,k,l}, isValidReg );
                   warpFd{i,j,k,l} = selectFd( warpFd{i,j,k,l}, isValidReg );
                   baselineFd = selectFd( baselineFd, isValidReg );
                   ref = ref( isValidReg, : );
                   type = type( isValidReg );
                   jperf.JHtov = jperf.JHtov( isValidReg );
                   jperf.JHwd = jperf.JHwd( isValidReg );
                   jperf.PP = jperf.PP( isValidReg );
                   part = part( isValidReg, : );
               end

               if k > 1
                   % perform a decomposition analysis
                   decomp{i,j,k,l} = regDecomp( ...
                           baselineFd, ...
                           vgrfFd{i,j,k,l}, ...
                           warpFd{i,j,k,l} );
               end
               
               if isempty( vgrfPCA{i,j,k,l} )

                   % run principal component analsyis
                   vgrfPCA{i,j,k,l} = pcaVGRF( vgrfFd{i,j,k,l}, ...
                                              fdPar{i,j}, ...
                                              warpFd{i,j,k,l}, ...
                                              setup.pca.nComp, ...
                                              setup.pca.nCompWarp );

               end

               if isempty( vgrfACP{i,j,k,l} )

                   % run analysis of characterising phases
                   vgrfACP{i,j,k,l} = acpVGRF( tSpan{i,j}, ...
                                             vgrfFd{i,j,k,l}, ...
                                             warpFd{i,j,k,l}, ...
                                             vgrfPCA{i,j,k,l} );

               end                      

               % generate output tables
               if isempty( results{i,j,k,l} ) 
                    results{i,j,k,l} = outputTable( name{i,j,k,l}, ...
                                            ref, ...
                                            type, ...
                                            jperf, ...
                                            vgrfPCA{i,j,k,l}, ...
                                            vgrfACP{i,j,k,l}, ...
                                            setup );
               end

               % fit models to the data
               if isempty( models{i,j,k,l} )
                   models{i,j,k,l} = fitVGRFModels( ...
                                        results{i,j,k,l}, ...
                                        part, setup, ...
                                        models{i,j,k,l} );
               end
               
           end

           % store data
           save( setup.filename, ...
                 'decomp', 'fdPar', 'name', 'vgrfFd', 'warpFd', ...
                 'isValid', 'regIter', 'vgrfPCA', 'vgrfACP', ...
                 'models' ); 
         

       end
       
  
       
    end
    
end


% ************************************************************************
%   Compile and save the results table
% ************************************************************************


longResults = compileResults( results );
longPerformance = compileResults( models, 'perf' );
longInclude = compileResults( models, 'incl' );
longTStat = compileResults( models, 'tStat' );
longCoeffRSq = compileResults( models, 'coeffRSq' );
longDecomp = compileResults( decomp );
faultyCountWOA = squeeze(cellfun( @sum, isValid(1,:,:) ))';
faultyCountALL = squeeze(cellfun( @sum, isValid(2,:,:) ))';


   
