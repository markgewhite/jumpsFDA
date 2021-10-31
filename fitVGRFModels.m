% ************************************************************************
% Function: fitVGRFModels
% Purpose:  Fit models to the VGRF data and evaluate their performance
%
% Parameters:
%       results: results table from the functional data analysis
%       trnSelect: training partitions to use (test partitions are inverse)
%       setup: parameters to use
%
% Outputs:
%       models: cell array of structure listings of model performance 
%
% ************************************************************************

function models = fitVGRFModels( results, trnSelect, setup, models )


% setup partitions
nPartitions = size( trnSelect, 2 );
tstSelect = ~trnSelect;

% setup datasets and predictors
nPred = setup.pca.nComp + setup.pca.nCompWarp;
dataset = results.Dataset(1,:);

% find out if the data set has two arm conditions 
classifier = length(unique( results.Arms ))==2;

% define models
modelName = { 'JHtov', 'JHwd', 'PP', 'jumpType' };
modelDist = { 'normal', 'normal', 'normal', 'binomial' };
modelLink = { 'identity', 'identity', 'identity', 'logit' };
outcomeIdx = [ 5, 6, 7, 4 ];
if classifier
    chosenModels = 4;
else
    chosenModels = [ 1 2 3 ];
end
nModels = length( chosenModels );
    

% define predictor subsets
pred1st = 8;
predSetName = { 'PCAU', 'PCAV', 'ACPU', 'ACPV' };
predIdxSet{1} = pred1st : pred1st+nPred-1;
predIdxSet{2} = pred1st+nPred : pred1st+2*nPred-1;
predIdxSet{3} = pred1st+2*nPred : pred1st+3*nPred-1;
predIdxSet{4} = pred1st+3*nPred : pred1st+4*nPred-1;
chosenPredSets = [ 1 4 ]; % only unrotated PCA and varimax ACP
nPredSets = length( chosenPredSets );

if isempty( models )       

    % define model performance table
    varNames = { 'Dataset', 'Model', 'PredictorSet', ...
                    'TrainRSq', 'TrainRSq_NoWarp', ...
                    'TestRSq', 'TestRSq_NoWarp', ...
                    'TestRMSE', 'TestAccuracy', ...
                    'TestSensitivity' , 'TestSpecificity' };
    varTypes = [repmat( {'string'}, 1, 3) ...
                 repmat( {'double'}, 1, 8) ];
    models.perf = table( ...
                'Size', [ nModels*nPredSets length(varNames) ], ...
                'VariableNames', varNames, ...
                'VariableTypes', varTypes );
     
    % define model definition tables
    varNames = [ { 'Dataset', 'Model', 'PredictorSet' } ...
                    results.Properties.VariableNames(8:7+nPred) ];
    varTypes = [ { 'string', 'string', 'string' } ...
                    repmat( {'double'}, 1, nPred) ];
    models.incl = table( ...
                'Size', [ nModels*nPredSets length(varNames) ], ...
                'VariableNames', varNames, ...
                'VariableTypes', varTypes );
    models.coeff = models.incl;
    models.tStat = models.incl; 
    models.coeffRSq = models.incl;

    newFit = true;
    
else
    
    newFit = false;

end



m = 0;
for i = chosenModels

    for j = chosenPredSets
        
        m = m + 1;
        
        if newFit 
            % carry out model fitting without stepwise selection
            doStepwise = false;
                                
        else
            % check the fixed model above a performance threshold
            if models.perf.TestRSq(m) > setup.models.RSqMeritThreshold
                % fit a new model with stepwise selection
                doStepwise = true;
            else
                % skip this iteration and move to the next model
                continue
            end
            
        end
        
        if newFit
            % setup record    
            models.perf = setupRecord( models.perf, m, ...
                                dataset, modelName{i}, predSetName{j} );
            models.incl = setupRecord( models.incl, m, ...
                                dataset, modelName{i}, predSetName{j} );
            models.tStat = setupRecord( models.tStat, m, ...
                                dataset, modelName{i}, predSetName{j} );
            models.coeffRSq = setupRecord( models.coeffRSq, m, ...
                                dataset, modelName{i}, predSetName{j} );
        end
            
        predIdxAll = predIdxSet{j};
               
        nPW = length(predIdxAll) - setup.pca.nCompWarp;
        predIdxAmpl = predIdxAll( 1:nPW );
        
        noWarpPred = sum( table2array( ...
                results(:, predIdxAll(nPW+1:end)) ), 'all' )==0;
        if noWarpPred
            predIdxAll = predIdxAmpl;
        end
            
        % setup records of fold-specific model fits
        nP = length( predIdxAll );
        trnRSq = zeros( nPartitions, 2 );
        tstRSq = zeros( nPartitions, 2 );
        inModel = false( nPartitions, nP );
        tStat = zeros( nPartitions, nP );
        coeffRSq = zeros( nPartitions, nP );
        tstRMSE = zeros( nPartitions, 1 );
        tstAcc = zeros( nPartitions, 1 );
        tstSen = zeros( nPartitions, 1 );
        tstSpc = zeros( nPartitions, 1 );

        
        for k = 1:nPartitions          
      
            % create training and testing sets
            trnData = results( trnSelect(:,k), ...
                                [ predIdxAll outcomeIdx(i) ] );
            tstData = results( tstSelect(:,k), ...
                                [ predIdxAll outcomeIdx(i) ] );
    
            % fit the full model with all predictors
            [ trnRSq(k,1), tstRSq(k,1), ...
                inModel(k,:), tStat(k,:), coeffRSq(k,:), ...
                tstRMSE(k), tstAcc(k), tstSen(k), tstSpc(k)] = ...
                          fitModel( ...
                                    trnData, tstData, ...
                                    modelDist{i}, modelLink{i}, ...
                                    doStepwise, setup.models );

            
            if ~noWarpPred
                % create the revised data sets without warp predictors
                trnData = results( trnSelect(:,k), ...
                                    [ predIdxAmpl outcomeIdx(i) ] );
                tstData = results( tstSelect(:,k), ...
                                    [ predIdxAmpl outcomeIdx(i) ] );

                % fit the model with only amplitude predictors
                [ trnRSq(k,2), tstRSq(k,2) ] = ...
                                fitModel( ...
                                        trnData, tstData, ...
                                        modelDist{i}, modelLink{i}, ...
                                        doStepwise, setup.models );
                                    
            end
            
        end
                                       
        
        % average performance across all partitions
        models.perf.TrainRSq(m) = mean( trnRSq(:,1) );            
        models.perf.TestRSq(m) = mean( tstRSq(:,1) );
        models.perf.TrainRSq_NoWarp(m) = mean( trnRSq(:,2) );            
        models.perf.TestRSq_NoWarp(m) = mean( tstRSq(:,2) );
        
        models.perf.TestRMSE(m) = mean( tstRMSE );
        models.perf.TestAccuracy(m) = mean( tstAcc );
        models.perf.TestSensitivity(m) = mean( tstSen );
        models.perf.TestSpecificity(m) = mean( tstSpc );
        
        % count number of variables included
        models.incl( m, 4:nP+3 ) = array2table(sum( inModel ));
        
        % average coefficients that were included
        models.tStat( m, 4:nP+3 ) = array2table(mean(tStat,'omitnan'));
        models.coeffRSq( m, 4:nP+3 ) = array2table(mean(coeffRSq,'omitnan'));

             
        disp(['Model: type = ' modelName{i} ...
          '; vars = ' predSetName{j} ...
          '; train rsq = ' num2str(models.perf.TrainRSq(m)) ...
          '; test rsq = ' num2str(models.perf.TestRSq(m)) ...
          '; -- No Warp Model: train rsq = ' ...
                            num2str(models.perf.TrainRSq_NoWarp(m)) ...
          '; test rsq = ' num2str(models.perf.TestRSq_NoWarp(m)) ]);
        
    end
    
end


end



function [ trnRSq, tstRSq, inModel, tStat, explRSq, ...
            tstRMSE, tstAcc, tstSen, tstSpc ] = ...
                                fitModel( trnData, tstData, ...
                                          distFn, linkFn, ...
                                          doStepwise, setup )
    warning( 'off', 'all' );
    % fit model
    if doStepwise
        mdl =stepwiseglm( trnData, ...
                          setup.spec, ...
                          'Distribution', distFn, ...
                          'Link', linkFn, ...
                          'Upper', setup.upper, ...
                          'Criterion', setup.criterion, ...
                          'Verbose', 0);

    else
        mdl = fitglm( trnData, ...
                      setup.spec, ...
                      'Distribution', distFn, ...
                      'Link', linkFn );

    end
    
    % get ground truth
    trnY = table2array( trnData(:,end) );
    tstY = table2array( tstData(:,end) );

    % make predictions
    trnYHat = predict( mdl, trnData );
    tstYHat = predict( mdl, tstData );
    
    % get model fit
    trnRSq = corr( trnY, trnYHat )^2;
    tstRSq = corr( tstY, tstYHat )^2;

    % compute errors
    if strcmp(linkFn,'logit')       
        % is classifier model
        tstYHat = round( tstYHat, 0 );
        tstClassPerf = classperf( tstY, tstYHat );
        tstAcc = tstClassPerf.CorrectRate;
        tstSen = tstClassPerf.Sensitivity;
        tstSpc = tstClassPerf.Specificity;
        tstRMSE = 0;

    else
        % is regression model
        tstRMSE = sqrt( mean( (tstYHat-tstY).^2) );
        tstAcc = 0;
        tstSen = 0;
        tstSpc = 0;

    end
    
    % record model coefficients
    inModel = mdl.VariableInfo.InModel(1:end-1);
    nCoeff = length( inModel );
    tStat = zeros( nCoeff, 1 );
    j = 1;
    for i = 1:nCoeff
        if inModel(i)
            j = j+1;
            tStat(i) = mdl.Coefficients.tStat(j);
        else
            tStat(i) = NaN;
        end
    end
    
    % obtain explained variances
    explRSq = zeros( nCoeff, 1 );
    for i = 1:nCoeff
        if inModel(i)
            mdl1 = removeTerms( mdl, ...
                            mdl.VariableInfo.Properties.RowNames{i} );
            explRSq(i) = mdl.Rsquared.Ordinary - mdl1.Rsquared.Ordinary;
        end
    end
    
    warning( 'on', 'all' );
    
end



function T = setupRecord( T, m, dataset, modelName, predSetName )

T.Dataset(m) = dataset;
T.Model(m) = modelName;
T.PredictorSet(m) = predSetName;
            
end