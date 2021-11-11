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

% define model tables
varNamesModels = { 'Dataset', 'Model', 'PredictorSet' };
varNamesPerf = { 'LogLikelihood', 'AIC', ...
                'TrainRSq', 'TestRSq', ...
                'TrainRMSE', 'TestRMSE', ...
                'TrainAccuracy', 'TestAccuracy', ...
                'ExplByWarp' };
varNamesDef = results.Properties.VariableNames(8:7+nPred);
nVarModels = length( varNamesModels );
nVarPerf = length( varNamesPerf );
nVarDef = length( varNamesDef );
            
if isempty( models )       

    % define model performance table
    varTypes = [ repmat({'string'}, 1, nVarModels) ...
                 repmat({'double'}, 1, nVarPerf) ];
    models.perf = table( ...
                'Size', [ nModels*nPredSets nVarModels+nVarPerf ], ...
                'VariableNames', [varNamesModels varNamesPerf], ...
                'VariableTypes', varTypes );
    models.stderr = models.perf;
     
    % define model definition tables
    varTypes = [ repmat({'string'}, 1, nVarModels) ...
                 repmat({'double'}, 1, nVarDef) ];
    models.incl = table( ...
                'Size', [ nModels*nPredSets nVarModels+nVarDef ], ...
                'VariableNames', [varNamesModels varNamesDef], ...
                'VariableTypes', varTypes );
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
            models.stderr = setupRecord( models.stderr, m, ...
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
            
        % setup ararys for fold-specific model fits
        nPred = length( predIdxAll );
        perf = table( ...
                'Size', [ nPartitions nVarPerf ], ...
                'VariableNames', varNamesPerf, ...
                'VariableTypes', repmat({'double'}, 1, nVarPerf) );
        inModel = zeros( nPartitions, nPred );
        tStat = zeros( nPartitions, nPred );
        coeffRSq = zeros( nPartitions, nPred );

        for k = 1:nPartitions          
      
            % create training and testing sets
            trnData = results( trnSelect(:,k), ...
                                [ predIdxAll outcomeIdx(i) ] );
            tstData = results( tstSelect(:,k), ...
                                [ predIdxAll outcomeIdx(i) ] );
    
            % fit the model for this partition
            [ perf(k,:), inModel(k,:), tStat(k,:), coeffRSq(k,:) ] = ...
                          fitModel( perf(k,:), ...
                                    trnData, tstData, ...
                                    modelDist{i}, modelLink{i}, ...
                                    doStepwise, setup.models );
                                            
        end
                                              
        % average performance across all partitions
        models.perf( m, 4:nVarPerf+3 ) = array2table(mean( table2array(perf) ));
        
        % calculate the standard error in this estimate
        models.stderr( m, 4:nVarPerf+3 ) = ...
            array2table( std( table2array(perf) )/sqrt( nPartitions ) );
        
        % count number of variables included
        models.incl( m, 4:nPred+3 ) = array2table(sum( inModel ));
        
        % average coefficients that were included
        models.tStat( m, 4:nPred+3 ) = array2table(mean( tStat,'omitnan'));
        models.coeffRSq( m, 4:nPred+3 ) = array2table(mean( coeffRSq, 'omitnan' ));

        if classifier
            disp(['Model: type = ' modelName{i} ...
              '; vars = ' predSetName{j} ...
              '; loglikelihood = ' num2str(models.perf.LogLikelihood(m),'%.1f') ...
              '; train fit = ' num2str(models.perf.TrainAccuracy(m),'%.3f') ...
                ' ' char(177) ' ' num2str(models.stderr.TrainAccuracy(m),'%.3f') ...
              '; test fit = ' num2str(models.perf.TestAccuracy(m),'%.3f') ...
                ' ' char(177) ' ' num2str(models.stderr.TestAccuracy(m),'%.3f')]);
        else
            disp(['Model: type = ' modelName{i} ...
              '; vars = ' predSetName{j} ...
              '; loglikelihood = ' num2str(models.perf.LogLikelihood(m),'%.1f') ...
              '; train fit = ' num2str(models.perf.TrainRSq(m),'%.3f') ...
                ' ' char(177) ' ' num2str(models.stderr.TrainRSq(m),'%.3f') ...
              '; test fit = ' num2str(models.perf.TestRSq(m),'%.3f') ...
                ' ' char(177) ' ' num2str(models.stderr.TestRSq(m),'%.3f')]);
        end
        
    end
    
end


end



function [ perf, inModel, tStat, explRSq ] = ...
                                fitModel( perf, trnData, tstData, ...
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
    perf.LogLikelihood = mdl.LogLikelihood;
    perf.AIC = mdl.ModelCriterion.AIC;
    
    perf.TrainRSq = corr( trnY, trnYHat )^2;
    perf.TestRSq = corr( tstY, tstYHat )^2;

    % compute errors
    if strcmp(linkFn,'logit')       
        % is classifier model
        trnYHat = round( trnYHat, 0 );
        trnClassPerf = classperf( trnY, trnYHat );
        perf.TrainAccuracy = trnClassPerf.CorrectRate;
        perf.TrainRMSE = 0;      
        
        tstYHat = round( tstYHat, 0 );
        tstClassPerf = classperf( tstY, tstYHat );
        perf.TestAccuracy = tstClassPerf.CorrectRate;
        perf.TestRMSE = 0;

    else
        % is regression model
        perf.TrainRMSE = sqrt( mean( (trnYHat-trnY).^2) );
        perf.TrainAccuracy = 0;
                
        perf.TestRMSE = sqrt( mean( (tstYHat-tstY).^2) );
        perf.TestAccuracy = 0;


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
    perf.ExplByWarp = 0;
    for i = 1:nCoeff
        if inModel(i)
            mdl1 = removeTerms( mdl, ...
                            mdl.VariableInfo.Properties.RowNames{i} );
            explRSq(i) = mdl.Rsquared.Ordinary - mdl1.Rsquared.Ordinary;
            if mdl.VariableInfo.Properties.RowNames{i}(6) == 'W'
                perf.ExplByWarp = perf.ExplByWarp + explRSq(i);
            end
        end
    end
    
    warning( 'on', 'all' );
    
end



function T = setupRecord( T, m, dataset, modelName, predSetName )

T.Dataset(m) = dataset;
T.Model(m) = modelName;
T.PredictorSet(m) = predSetName;
            
end