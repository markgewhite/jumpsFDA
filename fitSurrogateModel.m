% ************************************************************************
% Function: fitSurrogateModel
% Purpose:  Fit surrogate model to the GLMs' performance.
%           This is specific to the accuracy figures.
%
% Parameters:
%       perf: cell array of model performance tables
%
% Outputs:
%       model: surrogate model
%
% ************************************************************************

function [sModel, sModelInt] = fitSurrogateModel( glmModels )

% assemble tables for both jump types together 
perf1 = compileResults( glmModels(1,:,:), 'perf' ); 
% retain only classification results
retain1 = strcmp( perf1.Model, 'jumpType' );

% assemble tables for jumps without arm swing
perf2 = compileResults( glmModels(2,:,:), 'perf' );
% retain only regression results
retain2 = strcmp( perf2.Model, 'JHtov' ) ...
           | strcmp( perf2.Model, 'JHwd' ) ...
           | strcmp( perf2.Model, 'PP' );

% combine
perf = [ perf1( retain1, : ); perf2( retain2, : ) ];

% exclude ACPU and ACPV
retain3 = strcmp( perf.PredictorSet, 'PCAU' ) ...
           | strcmp( perf.PredictorSet, 'ACPV' );
perf = perf( retain3, : );

% create long table so there is one line for train and test
nRows = size( perf, 1 );
perf = [ perf; perf ];

isTrain = 1:2*nRows <= nRows;

% assign classes to each of the four blocks
perf.Set( isTrain  ) = "Train";
perf.Set( ~isTrain  ) = "Test";
perf.Outcome( isTrain ) = perf.TrainRSq( isTrain );
perf.Outcome( ~isTrain ) = perf.TestRSq( ~isTrain );

% add signal length standardisation method 
usedLLN = extract( perf.Dataset, 3 )=='1';
perf.Method( usedLLN ) = "LLN";
perf.Method( ~usedLLN ) = "PAD";

% extract the registration from the dataset label
perf.Registration = extractBetween( perf.Dataset, 8, 11 );


% retain only landmark registration (remove continuous reg)
isLMReg = ~strcmp( perf.Registration, "0000" );
perf = perf( isLMReg, : );

% re-badge no registration
isNoReg = strcmp( perf.Registration, "----" );
perf.Registration( isNoReg ) = "0000";

% remove unwanted fields
perf.Dataset = [];
perf.TrainRSq = [];
perf.TrainRSq_NoWarp = [];
perf.TestRSq = [];
perf.TestRSq_NoWarp = [];


% now fit the model
sModel = fitglm( perf, 'ResponseVar', 'Outcome' );





end



