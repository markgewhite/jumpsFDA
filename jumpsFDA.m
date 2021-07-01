% ************************************************************************
% Script: jumpsFDA
% Purpose:  Perform batches of functional data analysis on the jumps time series
%
%
% ************************************************************************

clear;

% ************************************************************************
%     MASTER OPTION
% ************************************************************************

dataset = 'Training';


% ************************************************************************
%     Setup file paths
% ************************************************************************

if ismac
    rootpath = '/Users/markgewhite/Google Drive/PhD/Studies/Jumps';
else
    rootpath = 'C:\Users\markg\Google Drive\PhD\Studies\Jumps';
end

switch dataset 
    case 'Training'
        datapath = [ rootpath '\Data\Processed\Training' ];
    case 'Testing'
        datapath = [ rootpath '\Data\Processed\Testing' ];
    case 'All'
        datapath = [ rootpath '\Data\Processed\All' ];
end

if ismac
    datapath = strrep( datapath, '\', '/') ;
end


% ************************************************************************
%   Constants
% ************************************************************************

g = 9.812; % acceleration due to gravity

tLength = 5000; % (5000) duration in milliseconds %%%% FILTER %%%% 
tFreq = 1; % sampling frequency

nMeasures = 1; % number of measures under analysis
vgrfID = 1; % vertical ground reaction force identifier
vgrfWoaID = 2; % VGRF subset for jumps without arms
vgrfWaID = 3; % VGRF subset for jumps with arms

subjectRefID = 1; % identifier for subject ID index in ref
jumpRefID = 2; % identifier for jump ID in ref


% ************************************************************************
%   Command switches
% ************************************************************************

options.test = 'Registration'; % type of test to perform

options.doTimeNormalisation = false; % whether to do time normalisation
options.doTruncation = false; % whether to truncate time series at take-off
options.doFiltering = false; % whether to do low-pass filtering on the GRF data
options.doCrossValidation = false; % whether to do cross validation to find best value for lambda
options.doCrossValidation2 = false; % whether to do cross validation to find best value for number of bases
options.doCrossValidation3 = false; % whether to do validation based on jump height
options.doDifferentiation = false; % whether to take the first derivative of the smoothed function

options.doCheckFit = false; % whether to check the goodness of fit
options.doResidualCheck = false; % for this check, whether to check the residuals
options.doSortCheck = true; % for this check, whether to sort from largest error to smallest
options.doUseWtFunc = false; % whether to use the weighted functional data object
options.doPlotDerivatives = false; % whether to plot the function and its first two derivatives

options.lm.doReg0 = true; % flag whether to do landmark registration
options.lm.doReg0Calc = true; % flag whether to do the calculation (or read from a file)
options.lm.doReg1 = false; % flag whether to do landmark registration
options.lm.doReg1Calc = false; % flag whether to do the calculation (or read from a file)
options.ct.doReg = false; % flag whether to do continuous registration
options.ct.doRegCalc = false; % flag whether to do the calculation (or read from a file)

options.lm.doCurvePlots = false; % flag whether to show registration curve plots
options.lm.doRegHistogram = false; % flag whether to show the spread of landmark points
options.lm.doRegPlots = false; % flag whether to show registration curve plots

options.pca.doShowComponents = false; % flag whether to show plots of PCA components
options.pca.doVarimaxRotation = true; % flag whether to do a vaximax rotation on PCA components
options.pca.doComponentsValidation = false; % flag whether to assess performance fit
options.pca.doCentreFunctions = true; % flag whether to centre about the mean curve


% ************************************************************************
%   Baseline settings
% ************************************************************************

if options.doTruncation
    preLength = tLength; % time duration before take-off
    postLength = 0; % time duration after take-off
else
    preLength = tLength; % time duration before take-off
    postLength = tLength; % time duration after take-off
end   
tSpan = -preLength:1/tFreq:postLength; % time domian in seconds 

setup.data.tFreq = tFreq; % time intervals per second
setup.data.tLength = tLength; % duration of the extracted dataset 
setup.data.maxLength = tFreq*tLength+1; % maximum possible duration
setup.data.maxLength1 = tFreq*preLength+1;
setup.data.maxLength2 = tFreq*postLength+1;
setup.data.form = 'Vertical'; % data representation
setup.data.cutoffFreq = 15; % 15 Hz cut-off frequency for filtering
setup.data.initial = 1; % initial padding value

setup.Fd.basis = 400; % number of bases
setup.Fd.basisOrder = 5; % 5th order for a basis expansion of quartic splines
setup.Fd.penaltyOrder = 3; % roughness penalty
setup.Fd.lambda = 1E0 ; % roughness penalty
setup.Fd.names = [{'Time (ms)'},{'Jumps'},{'GRF (BW)'}]; % axes names

setup.lm.basis = 10; % numbers of bases for landmark registration
setup.lm.basisOrder = 1; % time warping basis order for landmark registration
setup.lm.lambda = 1E0; % roughness penalty for landmark registration
setup.lm.yLambda = 1E0; % roughness penalty to prevent wiggles in y
setup.lm.searchLambda = 1E5; % heavily smoothed for LM search
setup.lm.filename = fullfile(datapath,'registrationLM.mat'); % where to save the analysis

setup.lm.set.grfmin = false; % use VGRF minimum as a landmark?
setup.lm.set.grfcross = false; % use VGRF crossing point as a landmark?
setup.lm.set.pwrmin = false; % use Power minimum as a landmark?
setup.lm.set.pwrcross = false; % use Power crossing point as a landmark?
setup.lm.set.pwrmax = false; % use Power maximum as a landmark?

setup.ct.basis = 50; % numbers of bases for continuous registration
setup.ct.basisOrder = 4; % time warping basis order for continuous registration
setup.ct.lambda = 1E-1; % roughness penalty for continuous registration
setup.ct.filename = fullfile(datapath,'registrationCT.mat'); % where to save the analysis

setup.pca.nComponents = 15; % number of PCA components to be retained


% ************************************************************************
%   Read data file
% ************************************************************************


% read the results file
load(fullfile(datapath,'processedjumpdata.mat'));


% ************************************************************************
%   Read data file
% ************************************************************************

cd(datapath);

% read the processed data file
load(fullfile(datapath,'processedjumpdata.mat'));


% nSubjects = 1; % **** OVERRIDE ****


% ************************************************************************
%   Normalise force
% ************************************************************************

% VGRF data (vertical jumps only)
nTotal = sum( nJumpsPerSubject );
subjectExclusions = find( sDataID==100 ); % [2,3];
jumpExclusions = [3703 3113 2107 2116 0503 0507 6010 1109];
vgrfData = cell( nTotal, 1 );
vgrfRef = zeros( nTotal, 2 );
withArms = false( nTotal, 1 );
initiation = zeros( nTotal, 1 );
takeoff = zeros( nTotal, 1 );
flightTime = zeros( nTotal, 1 );
k = 0;
for i = 1:nSubjects
    for j = 1:nJumpsPerSubject(i)
        jump = jumpOrder( sJumpID==sDataID(i), j );
        jumpID = sDataID(i)*100+j;
        duration = grf.takeoff(i,j)-grf.initiation(i,j)+1;
        if jump{1}(1) == 'V' ...
                && duration < setup.data.maxLength ...
                && ~ismember( i, subjectExclusions ) ...
                && ~ismember( jumpID, jumpExclusions )
            k = k+1;
            if options.doTruncation
                vgrfData{ k } = ...
                    grf.raw{i,j,totalID}( ...
                        grf.initiation(i,j):grf.takeoff(i,j), vAxisID ) ...
                               / bwall(i,j);
            else
                vgrfData{ k } = ...
                    grf.raw{i,j,totalID}(:,vAxisID) ...
                               / bwall(i,j);
            end
            vgrfRef( k, subjectRefID ) = i;
            vgrfRef( k, jumpRefID ) = j;
            withArms( k ) = (length(jump{1}) == 2);
            initiation( k ) = grf.initiation(i,j);
            takeoff( k ) = grf.takeoff(i,j);
            flightTime( k ) = jumpflighttime( ...
                                    grf.raw{i,j,totalID}(:,vAxisID), 1 );
        end
    end
end
vgrfData = vgrfData(1:k);
vgrfRef = vgrfRef(1:k,:);
withArms = withArms(1:k,:);
initiation = initiation(1:k);
takeoff = takeoff(1:k);
flightTime = flightTime(1:k);

% combine all datasets together
curveSet = { vgrfData, vgrfData(~withArms), vgrfData(withArms) };
curveIDSet = { vgrfRef, vgrfRef(~withArms,:), vgrfRef(withArms,:) };
curveTOSet = { takeoff, takeoff(~withArms), takeoff(withArms) };
curveFTSet = { flightTime, flightTime(~withArms), flightTime(withArms) };


% ************************************************************************
%   Performance the functional data analysis
% ************************************************************************

jumpperf.bwall = bwall;

grfFd = cell(nMeasures,1);
grfFdParams = cell(nMeasures,1);
grfPCA = cell(nMeasures,1);
warpFd = cell(nMeasures,1);
LM = cell(nMeasures,1);
name = cell( nMeasures, 1 );
results = cell(2*nMeasures,2);
m = 0;
for i = 1:nMeasures
    
    tic;
    
    % specify the registration setup
    switch options.test
        case 'Registration'
            ib = dec2bin( i-1, 5 );
            options.doTimeNormalisation = (ib(1)=='1');
            if options.doTimeNormalisation
                setup.Fd.basis = 256;
            else
                setup.Fd.basis = 450;
            end
            options.lm.doReg0 = not(strcmp( ib(2:5), '0000' ));
            setup.lm.set.grfmin = (ib(5)=='1');
            setup.lm.set.pwrmin = (ib(4)=='1');
            setup.lm.set.pwrcross = (ib(3)=='1');
            setup.lm.set.pwrmax = (ib(2)=='1');
            name{ i } = [ 'VGRF' ib ];
        case 'LMorder'
            setup.lm.basisOrder = i-1;
            name{ i } = [ 'VGRF Order' num2str(i-1) ];
        case 'LMlambda'
            exponent = i-2;
            setup.lm.lambda = 10^exponent;
            name{ i } = [ 'VGRF Lambda 10E' num2str(exponent) ];
        case 'LMnbasis'
            setup.lm.basis = i*5;
            name{ i } = [ 'VGRF nBases ' num2str(setup.lm.basis) ];
    end
    
    disp(['MEASURE: ' char(name(i))]);
       
    % perform the functional data analysis
    [grfFd{i}, grfFdParams{i}, grfPCA{i}, warpFd{i}, LM{i}] = ...
                                    performFDA( ...
                                                vgrfData(~withArms), ...
                                                tSpan, ...
                                                takeoff, ...
                                                setup, ...
                                                options, ...
                                                @findGRFlandmarks, ...
                                                []);  
    
    % recalculate jump performances
    if not(options.doTimeNormalisation)
        [newperf.height, newperf.peakPower] = ...
                            jumpperf_pca(tSpan,grfPCA{i}.unrotated);
    else
        newperf.height = zeros( nSubjects, nJumps );
        newperf.peakPower = zeros( nSubjects, nJumps );
    end
    
    % generate the output tables

    % 1: unrotated PCA scores
    fullname = [num2str(i,'%02d') '-' name{i} 'UA'];
    m = m+1;
    [ results{m,1}, results{m,2} ] = outputPCAtables( ...
                    fullname, vgrfRef, sDataID, jumpOrder, sJumpID, ...
                    grfFd{i}, ...
                    grfPCA{i}.unrotated, jumpperf, newperf );
                
    % 2: time domain and time warp PCA scores
    if options.lm.doReg0 || options.ct.doReg
        fullname = [num2str(i,'%02d') '-' name{i} 'UW'];
        m = m+1;
        [ results{m,1}, results{m,2} ] = outputPCAtables( ...
                        fullname, vgrfRef, sDataID, jumpOrder, sJumpID, ...
                        grfFd{i}, ...
                        grfPCA{i}.warp, jumpperf, newperf );
    end        
        
    % 3: varimax PCA scores
    if options.pca.doVarimaxRotation
        if not(options.doTimeNormalisation)
            [newperf.height, newperf.peakPower] = ...
                                jumpperf_pca(tSpan,grfPCA{i}.varimax);
        else
            newperf.height = zeros( nSubjects, nJumps );
            newperf.peakPower = zeros( nSubjects, nJumps );
        end
        fullname = [num2str(i,'%02d') '-' name{i} 'VA'];
        m = m+1;
        [ results{m,1}, results{m,2} ] = outputPCAtables( ...
                        fullname, vgrfRef, sDataID, jumpOrder, sJumpID, ...
                        grfFd{i}, ...
                        grfPCA{i}.varimax, jumpperf, newperf );
                    
        if options.lm.doReg0 || options.ct.doReg
            fullname = [num2str(i,'%02d') '-' name{i} 'VW'];
            m = m+1;
            [ results{m,1}, results{m,2} ] = outputPCAtables( ...
                            fullname, vgrfRef, sDataID, jumpOrder, sJumpID, ...
                            grfFd{i}, ...
                            grfPCA{i}.warp, jumpperf, newperf );
        end  
    end
        
    duration = toc;
    disp(['Processing time = ' num2str(duration)]);
    
end
results = results(1:m,:);


% ************************************************************************
%   Save the results files
% ************************************************************************

% find the broadest width to the results tables
maxN = 0;
for i = 1:m
    maxN = max(maxN,size(results{i,1},2));
end

% standardise the results table widths
for i = 1:m
    if size(results{i,1},2)<maxN
        rows = size(results{i,1},1);
        cols = size(results{i,1},2);
        results{i,1} = [results{i,1} array2table(zeros(rows,maxN-cols))];
        for j = cols-9:maxN-10
            results{i,1}.Properties.VariableNames{10+j} = char("PCA"+j);
        end
    end
end

% combine the results tables together
stdResults = results{1,1};
longResults = results{1,2};
for i = 2:m
    stdResults = outerjoin(stdResults,results{i,1},'MergeKeys',true);
    longResults = outerjoin(longResults,results{i,2},'MergeKeys',true);
end
    
stdName = fullfile(datapath,'All PCA-ACP Results.csv');
writetable(stdResults,stdName,'WriteRowNames', true);
disp(['Output file saved: ' stdName]);

longName = strcat(stdName(1:end-4),'-long.csv');
writetable(longResults,longName,'WriteRowNames', true);
disp(['Output file saved: ' longName]);


% ************************************************************************
%   Save the Matlab variables
% ************************************************************************

save(fullfile(datapath,'GRFFeatures-Test.mat'), ...
         'grfFd', 'grfFdParams', ...
         'LM', 'warpFd', ...
         'grfPCA', ...
         'curveSet', 'curveIDSet', 'curveTOSet', 'curveFTSet', ...
         'stdResults');
disp(['Matlab file saved: ' fullfile(datapath,'GRFFeatures-Test.mat')]);
   
