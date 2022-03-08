% ************************************************************************
% Function: fitMetaModels
% Purpose:  Fit meta models to the GLMs' performance.
%
% Parameters:
%       perf: cell array of model performance tables
%       type: type of model, either 'Linear' or 'Interactions'
%       predSplit: whether to split out the predictor types
%
% Outputs:
%       model: surrogate model
%
% ************************************************************************

function [metaJHtov, metaJHwd, metaPP, metaCL, metaAll ] = ...
                            fitMetaModels( glmModels, type, predSplit )

if nargin < 3
    predSplit = false;
end
if nargin < 2
    type = 'linear';
end

% assemble tables for regression models
% -------------------------------------
mdlOut1 = compileResults( glmModels(1,:,:,:), 'perf' );

% retain only regression results
mdlOut1.Model = categorical( mdlOut1.Model );
retain1 = mdlOut1.Model=='JHtov' ...
            | mdlOut1.Model=='JHwd' | mdlOut1.Model=='PP';
mdlOut1 = mdlOut1( retain1, : );

% create common outcome variable for train and test
mdlRegression = setCommonOutcome( mdlOut1, 'RMSE' );

% add other factors and remove others
mdlRegression = addFactors( mdlRegression, predSplit );
mdlRegression = removeUnwanted( mdlRegression );
mdlRegression.LogLikelihood = [];


% assemble tables for classification model
% ----------------------------------------
mdlOut2 = compileResults( glmModels(2,:,:,:), 'perf' ); 
% retain only classification results
mdlOut2.Model = categorical( mdlOut2.Model );
retain2 = mdlOut2.Model=='jumpType';
mdlOut2 = mdlOut2( retain2, : );

% create common outcome variable for train and test
mdlClassification = setCommonOutcome( mdlOut2, 'Accuracy' );

% add other factors and remove others
mdlClassification = addFactors( mdlClassification, predSplit );
mdlClassification = removeUnwanted( mdlClassification );
mdlClassification.LogLikelihood = [];


% create the outcome specific models
% ----------------------------------
mdlJHtov = mdlRegression( mdlRegression.Model=='JHtov', : );
mdlJHtov.Model = [];
disp('Fitting JHtov meta model:');
metaJHtov = fitModel( mdlJHtov, type, 'Outcome' );

mdlJHwd = mdlRegression( mdlRegression.Model=='JHwd', : );
mdlJHwd.Model = [];
disp('Fitting JHwd meta model:');
metaJHwd = fitModel( mdlJHwd, type, 'Outcome' );

mdlPP = mdlRegression( mdlRegression.Model=='PP', : );
mdlPP.Model = [];
disp('Fitting PP meta model:');
metaPP = fitModel( mdlPP, type, 'Outcome' );

mdlCL = mdlClassification;
mdlCL.Model = [];
disp('Fitting classification meta model:');
metaCL = fitModel( mdlCL, type, 'Outcome' );


% assemble table covering all models
% ----------------------------------
mdlAll = [ mdlOut1; mdlOut2 ];
mdlAll = addFactors( mdlAll, predSplit );
mdlAll = removeUnwanted( mdlAll );

% ensure all LL's are positive for a gamma distribution
mdlAll.LogLikelihood = mdlAll.LogLikelihood - min(mdlAll.LogLikelihood) ...
    + 5000;

% fit log likelihood model
disp('Fitting LL meta model:');
metaAll = fitModel( mdlAll, type, 'LogLikelihood' );


end


function mdl = fitModel( data, type, response )

switch type
    case 'linear'
        mdl = fitglm( data, ...
                         'ResponseVar', response, ...
                         'Distribution', 'Gamma', ...
                         'Link', 'identity' );
    case 'interactions'
        mdl = stepwiseglm( data, 'interactions', ...
                         'ResponseVar', response, ...
                         'Distribution', 'Gamma', ...
                         'Link', 'identity', ...
                         'Criterion', 'BIC' );

end

end



function T = setCommonOutcome( T, measure )

% create long table so there is one line for train and test
nRows = size( T, 1 );
T = [ T; T ];

isTrain = 1:2*nRows <= nRows;
measureTrain = ['Train' measure];
measureTest = ['Test' measure];

T.Set( isTrain  ) = "Train";
T.Set( ~isTrain  ) = "Test";
T.Outcome( isTrain ) = T.(measureTrain)( isTrain );
T.Outcome( ~isTrain ) = T.(measureTest)( ~isTrain );

end


function T = addFactors( T, predSplit )

% add signal length standardisation method 
usedPAD = extract( T.Dataset, 3 )=='1';
T.Method( usedPAD ) = "PAD";
T.Method( ~usedPAD ) = "LTN";
T.Method = categorical( T.Method, {'PAD','LTN'}, 'Ordinal', true );

if predSplit
    % extract component type
    pca = extractBefore( T.PredictorSet, 4 )=='PCA';
    T.Feature( pca ) = "PCA";
    T.Feature( ~pca ) = "ACP";
    T.Feature = categorical( T.Feature, {'PCA','ACP'}, 'Ordinal', true );
    
    % extract rotation employed
    varimax = extract( T.PredictorSet, 4 )=='V';
    T.Rotation( ~varimax ) = "Unrotated";
    T.Rotation( varimax ) = "Varimax";
    T.Rotation = categorical( T.Rotation, ...
                            {'Unrotated','Varimax'}, 'Ordinal', true );

    T.PredictorSet = [];
else
    % provide alternative representation
    T.PredictorSet = categorical( T.PredictorSet, ...
                   {'PCAU', 'PCAV', 'ACPU', 'ACPV'}, 'Ordinal', true );
end

% extract landmark registration from the dataset label
T.RegLM = extractBetween( T.Dataset, 8, 11 );
T.RegLM = categorical( T.RegLM );

% extract continuous registration from the dataset label
T.RegC = extractAfter( T.Dataset, 11 );
T.RegC = categorical( T.RegC );

end


function T = removeUnwanted( T )

% remove unwanted fields
T.Dataset = [];
T.TrainRSq = [];
T.TestRSq = [];
T.TrainRMSE = [];
T.TestRMSE = [];
T.TrainAccuracy = [];
T.TestAccuracy = [];
T.AIC = [];
T.ExplByWarp = [];

end


