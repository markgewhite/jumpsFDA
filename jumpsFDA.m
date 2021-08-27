% ************************************************************************
% Script: jumpsFDA
% Purpose:  Perform batches of functional data analysis on the jumps time series
%
%
% ************************************************************************

clear;

% ************************************************************************
%     Setup file paths
% ************************************************************************

if ismac
    rootpath = '/Users/markgewhite/Google Drive/PhD/Studies/Jumps';
else
    rootpath = 'C:\Users\markg\Google Drive\PhD\Studies\Jumps';
end

datapath = [ rootpath '\Data\Processed\Training' ];

if ismac
    datapath = strrep( datapath, '\', '/') ;
end


% ************************************************************************
%   Read data file
% ************************************************************************

cd(datapath);

% read the processed data file
load(fullfile(datapath,'compactjumpdata.mat'));


% ************************************************************************
%   Constants
% ************************************************************************

g = 9.812; % acceleration due to gravity

% tLength = 10000; % (5000) duration in milliseconds %%%% FILTER %%%% 
tFreq = 1; % sampling frequency

nSets = 3; % number of data sets (1 = All; 2 = CMJ(NA); 3 = CMJ(A))
nStd = 2; % number of methods to standardise length
nModels = 5; % number of models
nProc = (2^4-1)+2; % number of curves processes

% presets based on approx 1E-3 mean GCV error
%preset.nBasis = [ 670 140; 530 135; 660 150 ];
%preset.lambda = [ 1E2 1E2; 1E2 1E2; 1E2 1E2 ];

% presets based on minimising GCV error
%preset.nBasis = [ 1180 350; 1180 335; 950 375 ];
%preset.lambda = [ 1E0 1E0; 1E0 1E0; 1E0 1E0 ];

% presets based on not exceed jump performance error threshold
preset.nBasis = [ 190 80; 190 80; 190 80 ];
preset.lambda = [ 1E3 1E3; 1E3 1E3; 1E3 1E3 ];

% ************************************************************************
%   Command switches
% ************************************************************************

options.test = 'None'; % type of test to perform


options.doFiltering = false; % whether to do low-pass filtering on the GRF data
options.doInitialFit = false; % whether to fit with max flexibility first

options.doCheckFit = false; % whether to check the goodness of fit
options.doResidualCheck = false; % for this check, whether to check the residuals
options.doSortCheck = false; % for this check, whether to sort from largest error to smallest
options.doUseWtFunc = false; % whether to use the weighted functional data object
options.doPlotDerivatives = false; % whether to plot the function and its first two derivatives

options.reg.doCurvePlots = false; % flag whether to show registration curve plots
options.reg.doRegHistogram = false; % flag whether to show the spread of landmark points
options.reg.doRegPlots = false; % flag whether to show registration curve plots


% ************************************************************************
%   Baseline settings
% ************************************************************************

setup.data.tFreq = 1; % time intervals per second
setup.data.sampleFreq = 1000; % sampling frequency
setup.data.cutoffFreq = 10; % 15 Hz cut-off frequency for filtering
setup.data.padding = 500; % milliseconds of padding for filtering
setup.data.form = 'Vertical'; % data representation
setup.data.initial = 1; % initial padding value
setup.data.threshold1 = 0.08; % primary detection threshold 
setup.data.threshold2 = 0.01; % secondary detection threshold
setup.data.sustained = 100; % milliseconds for secondary threshold 

setup.Fd.basisOrder = 4; % 5th order for a basis expansion of quartic splines
setup.Fd.penaltyOrder = 2; % roughness penalty
setup.Fd.lambda = 1E2 ; % roughness penalty
setup.Fd.names = [{'Time (ms)'},{'Jumps'},{'GRF (BW)'}]; % axes names
setup.Fd.tolerance = 0.001; % performance measure error tolerance

setup.reg.nIterations = 2; % Procrustes iterations
setup.reg.nBasis = 10; % numbers of bases for registration
setup.reg.basisOrder = 3; % time warping basis order for registration
setup.reg.wLambda = 1E-2; % roughness penalty for time warp
% *** CHANGE? ***
setup.reg.XLambda = 1E0; % roughness penalty to prevent wiggles in y

setup.reg.lm.grfmin = false; % use VGRF minimum as a landmark?
setup.reg.lm.grfcross = false; % use VGRF crossing point as a landmark?
setup.reg.lm.pwrmin = false; % use Power minimum as a landmark?
setup.reg.lm.pwrcross = false; % use Power crossing point as a landmark?
setup.reg.lm.pwrmax = false; % use Power maximum as a landmark?

setup.pca.nComp = 15; % number of PCA components to be retained
setup.pca.nCompWarp = 8; % number of PCA components to be retained

setup.models.nRepeats = 2; % number of repetitions of CV
setup.models.nFolds = 5; % number of CV folds for each repetition
setup.models.seed = 12345; % random seed for reproducibility
setup.models.spec = 'linear'; % type of GLM
setup.models.upper = 'linear'; % linear model without interactions
setup.models.criterion = 'bic'; % predictor selection criterion
setup.models.RSqMeritThreshold = 0.8; % merit threshold for stepwise selection

setup.filename = fullfile(datapath,'jumpsAnalysis2.mat'); % where to save the analysis



% ************************************************************************
%   Extract data
% ************************************************************************

% exclude jumps from subjects in the second data collection
subjectExclusions = find( ismember( sDataID, ...
            [ 14, 39, 68, 86, 87, 11, 22, 28, 40, 43, 82, 88, 95, 97, ...
              100, 121, 156, 163, 196 ] ) );

% specific jumps that should be excluded
jumpExclusions = [3703 3113 2107 2116 0503 0507 6010 1109];

[ rawData, refSet, typeSet ] =  extractVGRFData( ... 
                                    grf, bwall, nJumpsPerSubject, ...
                                    sDataID, sJumpID, jumpOrder, ...
                                    subjectExclusions, jumpExclusions );

                                
% ************************************************************************
%   Smooth data - first cut
% ************************************************************************


tSpan0 = cell( nSets, 1 );
for i = 1:nSets
    
    % standardised length with padding
    maxLen = max( cellfun( @length, rawData{i} ) );
    tSpan0{i} = -maxLen+1:0; % time domain in milliseconds 
    
    rawData{i} = padData( rawData{i}, ...
                           maxLen, ...
                           setup.data.initial );

    % filter data
    if options.doFiltering
        rawData{i} = filterVGRF( rawData{i}, ...
                                    setup.data.sampleFreq, ...
                                    setup.data.cutoffFreq, ...
                                    setup.data.padding );
    end
                       
end 


% ************************************************************************
%   Determine smoothing levels
% ************************************************************************

perf = cell( nSets );
tSpan = cell( nSets, nStd );
vgrfData{i} = cell( nSets, nStd );

for i = 1:nSets
   
    % determine jump initation
    tStart = jumpInit( rawData{i}, ...
                       setup.data.threshold1, ...
                       setup.data.threshold2, ...
                       setup.data.sustained );
                   
   
    for j = 1:nStd

       switch j
           case 1
               % pad out to longest series
               fixStart = min( tStart );
               
               % truncate pre-padded series to this length
               vgrfData{i,j} = rawData{i}( fixStart:end, : );              
               
           case 2
               % time normalising to median series length
               fixStart = fix( median( tStart ) );
               
               % truncate times series and time normalise
               vgrfData{i,j} = timeNormData( rawData{i}, ...
                                         tStart, ...
                                         fixStart );
                                                                         
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
            validateSmoothing(  vgrfData{i,j}, ...
                                tSpan{i,j}, ...
                                setup.Fd, ...
                                perf{i} );
            pause;
       end
       
    end
    
end



% ************************************************************************
%   Begin functional data analysis
% ************************************************************************

vgrfFd = cell( nSets, nStd, nProc );
warpFd = cell( nSets, nStd, nProc );
fdPar = cell( nSets, nStd );
decomp = cell( nSets, nStd, nProc );
name = strings( nSets, nStd, nProc );

vgrfPCA = cell( nSets, nStd, nProc );
vgrfACP = cell( nSets, nStd, nProc );



results = cell( nSets, nStd, nProc );
models = cell( nSets, nStd, nProc );
performance = cell( nSets, nStd, nProc );

load( setup.filename );



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
       if isempty( vgrfFd{i,j,1} ) || isempty( fdPar{i,j} )
           [ vgrfFd{i,j,1}, fdPar{i,j} ] = smoothVGRF(  ...
                                                vgrfData{i,j}, ...
                                                tSpan{i,j}, ...
                                                fdSetup, ...
                                                options );
           name{i,j,1} = [ num2str(i) '-' num2str(j) 'VGRF----' ];
       end
             
       % compute jump performances from the padded curves
       %if j == 1
       %     perf{i} = jumpperf_fd( vgrfFd{i,j,1} );
       %end
       
       
       for k = 1:nProc
           
           if k > 1
               % processing includes regression
               
               % select landmarks to find
               kbin = dec2bin( k-2, 4 );
               setup.reg.lm.grfmin = (kbin(1)=='1');
               setup.reg.lm.pwrmin = (kbin(2)=='1');
               setup.reg.lm.pwrcross = (kbin(3)=='1');
               setup.reg.lm.pwrmax = (kbin(4)=='1');
               name{i,j,k} = [ num2str(i) '-' num2str(j) 'VGRF' kbin ];

               % register the curves
               if isempty( vgrfFd{i,j,k} ) ...
                       || isempty( warpFd{i,j,k} ) ...
                       || isempty( decomp{i,j,k} )
                   
                    [ vgrfFd{i,j,k}, warpFd{i,j,k}, decomp{i,j,k} ] = ...
                                  registerVGRF( tSpan{i,j}, ...
                                                vgrfFd{i,j,1}, ...
                                                setup.reg, ...
                                                options.reg );
               
               end
                                                                      
           end
           
           if isempty( vgrfPCA{i,j,k} )
                                
               % run principal component analsyis
               vgrfPCA{i,j,k} = pcaVGRF( vgrfFd{i,j,k}, ...
                                          fdPar{i,j}, ...
                                          warpFd{i,j,k}, ...
                                          setup.pca.nComp, ...
                                          setup.pca.nCompWarp );
                 
           end
           
           if isempty( vgrfACP{i,j,k} )
               
               % run analysis of characterising phases
               vgrfACP{i,j,k} = acpVGRF( tSpan{i,j}, ...
                                         vgrfFd{i,j,k}, ...
                                         warpFd{i,j,k}, ...
                                         vgrfPCA{i,j,k} );
                             
           end
           
           % store data
           %save( setup.filename, ...
           %             'vgrfFd', 'warpFd', 'decomp', 'name', 'fdPar', ...
           %             'vgrfPCA', 'vgrfACP' );
                    
                    
                                        
           % generate output tables
           if isempty( results{i,j,k} ) 
                    results{i,j,k} = outputTable( name{i,j,k}, ...
                                        refSet{i}, ...
                                        typeSet{i}, ...
                                        perf{i}, ...
                                        vgrfPCA{i,j,k}, ...
                                        vgrfACP{i,j,k}, ...
                                        setup );
           end
                                    
           % fit models to the data
           models{i,j,k} = fitVGRFModels( ...
                                    results{i,j,k}, ...
                                    partitions, setup, ...
                                    models{i,j,k} );
                                                  
                                       
                                        
       end
       
       
    end
    
end


% ************************************************************************
%   Compile and save the results table
% ************************************************************************

save( setup.filename, ...
      'decomp', 'fdPar', 'name', 'vgrfFd', 'warpFd', ...
      'vgrfPCA', 'vgrfACP', ...
      'models' );

longResults = compileResults( results );
longPerformance = compileResults( performance );
longInclude = compileResults( models, 'incl' );
longCoeff = compileResults( models, 'coeff' );
longTStat = compileResults( models, 'tStat' );
longIncludeW = compileResults( models, 'inclW' );
longCoeffW = compileResults( models, 'coeffW' );
longTStatW = compileResults( models, 'tStatW' );


   
