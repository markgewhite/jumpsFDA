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
doClassifiers = length(unique( results.Arms ))==2;

% define models
nModels = 3 + doClassifiers;
modelName = { 'JHtov', 'JHwd', 'PP', 'jumpType' };
modelDist = { 'normal', 'normal', 'normal', 'binomial' };
modelLink = { 'identity', 'identity', 'identity', 'logit' };
outcomeIdx = [ 5, 6, 7, 4 ];

% define predictor subsets
Pred1st = 8;
nPredSets = 4;
predSetName = { 'PCAU', 'PCAV', 'ACPU', 'ACPV' };
predIdxSet{1} = Pred1st : Pred1st+nPred-1;
predIdxSet{2} = Pred1st+nPred : Pred1st+2*nPred-1;
predIdxSet{3} = Pred1st+2*nPred : Pred1st+3*nPred-1;
predIdxSet{4} = Pred1st+3*nPred : Pred1st+4*nPred-1;

if isempty( models )       

    % define model performance table
    varNames = { 'Dataset', 'Model', 'PredictorSet', ...
                    'TrainRSq', 'TrainRSq_NoWarp', ...
                    'TestRSq', 'TestRSq_NoWarp' };
    varTypes = { 'string', 'string', 'string', ...
                    'double', 'double', 'double', 'double' };
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
    models.inclW = models.incl;
    models.coeff = models.incl;
    models.coeffW = models.incl;     
    models.tStat = models.incl; 
    models.tStatW = models.incl;
    models.coeffRSq = models.incl;
    models.coeffRSqW = models.incl;

    newFit = true;
    
else
    
    newFit = false;

end



m = 12; %%%%% TEMP
for i = 4 % 1:nModels

    for j = [1 4 ] % 1:nPredSets
        
        %m = m + 1;
        % TEMPORARY
        if j==1
            m = m + 1;
        else
            m = m + 3;
        end
        
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
            models.inclW = setupRecord( models.inclW, m, ...
                                dataset, modelName{i}, predSetName{j} );
            models.coeff = setupRecord( models.coeff, m, ...
                                dataset, modelName{i}, predSetName{j} );
            models.coeffW = setupRecord( models.coeffW, m, ...
                                dataset, modelName{i}, predSetName{j} );
            models.tStat = setupRecord( models.tStat, m, ...
                                dataset, modelName{i}, predSetName{j} );
            models.tStatW = setupRecord( models.tStatW, m, ...
                                dataset, modelName{i}, predSetName{j} );
            models.coeffRSq = setupRecord( models.coeffRSq, m, ...
                                dataset, modelName{i}, predSetName{j} );
            models.coeffRSqW = setupRecord( models.coeffRSqW, m, ...
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
        inModel = false( nPartitions, nP, 2 );
        coeff = zeros( nPartitions, nP, 2 );
        tStat = zeros( nPartitions, nP, 2 );
        coeffRSq = zeros( nPartitions, nP, 2 );
        
        for k = 1:nPartitions          
      
            % create training and testing sets
            trnData = results( trnSelect(:,k), ...
                                [ predIdxAll outcomeIdx(i) ] );
            tstData = results( tstSelect(:,k), ...
                                [ predIdxAll outcomeIdx(i) ] );
    
            % fit the full model with all predictors
            [ trnRSq(k,1), tstRSq(k,1), inModel(k,:,1), ...
                coeff(k,:,1), tStat(k,:,1), coeffRSq(k,:,1) ] = ...
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
                [ trnRSq(k,2), tstRSq(k,2), ...
                    inModel(k,1:nPW,2), coeff(k,1:nPW,2), ...
                    tStat(k,1:nPW,2), coeffRSq(k,1:nPW,2) ] = ...
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
        
        % count number of variables included
        models.incl( m, 4:nP+3 ) = array2table(sum( inModel(:,:,1) ));
        models.inclW( m, 4:nP+3 ) = array2table(sum( inModel(:,:,2) ));
        
        % average coefficients that were included
        models.coeff( m, 4:nP+3 ) = array2table(mean(coeff(:,:,1),'omitnan'));
        models.coeffW( m, 4:nP+3 ) = array2table(mean(coeff(:,:,2),'omitnan'));
        models.tStat( m, 4:nP+3 ) = array2table(mean(tStat(:,:,1),'omitnan'));
        models.tStatW( m, 4:nP+3 ) = array2table(mean(tStat(:,:,2),'omitnan'));
        models.coeffRSq( m, 4:nP+3 ) = array2table(mean(coeffRSq(:,:,1),'omitnan'));
        models.coeffRSqW( m, 4:nP+3 ) = array2table(mean(coeffRSq(:,:,2),'omitnan'));

             
        disp(['Model: type = ' modelName{i} ...
          '; vars = ' predSetName{j} ...
          '; train error = ' num2str(models.perf.TrainRSq(m)) ...
          '; test error = ' num2str(models.perf.TestRSq(m)) ...
          '; -- No Warp Model: train error = ' ...
                            num2str(models.perf.TrainRSq_NoWarp(m)) ...
          '; test error = ' num2str(models.perf.TestRSq_NoWarp(m)) ]);
        
    end
    
end


end



function [ trnRSq, tstRSq, inModel, coeff, tStat, RSq ] = ...
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
    warning( 'on', 'all' );
    
    % get ground truth
    trnY = table2array( trnData(:,end) );
    tstY = table2array( tstData(:,end) );

    % make predictions
    trnYHat = predict( mdl, trnData );
    tstYHat = predict( mdl, tstData );

    % compute errors
    if strcmp(linkFn,'logit')
        
        % is classifier model
        trnRSq = accuracy( trnY, trnYHat );            
        tstRSq = accuracy( tstY, tstYHat ); 

    else
        % is regression model
        trnRSq = rsquared( trnY, trnYHat );
        tstRSq = rsquared( tstY, tstYHat );

    end
    
    % record model coefficients
    inModel = mdl.VariableInfo.InModel(1:end-1);
    nCoeff = length( inModel );
    coeff = zeros( nCoeff, 1 );
    tStat = zeros( nCoeff, 1 );
    j = 1;
    for i = 1:nCoeff
        if inModel(i)
            j = j+1;
            coeff(i) = mdl.Coefficients.Estimate(j);
            tStat(i) = mdl.Coefficients.tStat(j);
        else
            coeff(i) = NaN;
            tStat(i) = NaN;
        end
    end
    
    % obtain explained variances
    RSq = zeros( nCoeff, 1 );
    for i = 1:nCoeff
        if inModel(i)
            mdl1 = removeTerms( mdl, ...
                            mdl.VariableInfo.Properties.RowNames{i} );
            RSq(i) = mdl.Rsquared.Ordinary - mdl1.Rsquared.Ordinary;
        end
    end

end



function rsq = rsquared( y, yHat )

rsq = corr( y, yHat )^2;

end



function correct = accuracy( y, yHat )

c = y > 0.5;
cHat = yHat > 0.5;

isMatch = c==cHat;

correct = sum( isMatch )/length( y );

end



function T = setupRecord( T, m, dataset, modelName, predSetName )

T.Dataset(m) = dataset;
T.Model(m) = modelName;
T.PredictorSet(m) = predSetName;
            
end