% ************************************************************************
% Function: plotCoeffRSq
% Purpose:  Create a panel plot of explained variance by coefficient
%           with PAD models on the left and LTN models on the right
%           for the outcome variables specified
%
% Parameters:
%       models : cell array of model output structures
%       outcomeVars : string array of outcome variables to present
%                     ["JHtov", "JHwd", "PP", "CL"]
%
% Outputs:
%       one figure
%
% ************************************************************************

function fig = plotCoeffRSq( models, outcomeVars )

% compile table of coefficients to work with
C = compileResults( models, 'coeffRSq' );

% compile table of model performances 
P = compileResults( models, 'perf' );
% reverse accuracy to make lowest the best
P.TrainAccuracy = -P.TrainAccuracy;
P.TestAccuracy = -P.TestAccuracy;

% create categories
C.Model = categorical( C.Model );
C.PredictorSet = categorical( C.PredictorSet );

usedPAD = extract( C.Dataset, 3 )=='1';
C.Norm( usedPAD ) = "PAD";
C.Norm( ~usedPAD ) = "LTN";
C.Norm = categorical( C.Norm );

C.Dataset = categorical( C.Dataset );



% retain the unrotated PCA subset
retain = C.PredictorSet=='PCAU';
C = C( retain, : );
P = P( retain, : );

% determine how many plots
nVars = length( outcomeVars );

% create figure
fig = figure;
fig.Units = 'centimeters';
fig.Position(3) = 18.0; % width
fig.Position(4) = 2.0 + 6*nVars; % height

% iterate through outcome variables
for i = 1:nVars

    % identify the appropriate error field and y-limit
    switch outcomeVars(i)
        case {"JHtov", "JHwd"}
            errFld = 9;
            yMax = 50;
            name = "Jump Height";
        case "PP"
            errFld = 9;
            yMax = 70;
            name = "Peak Power";
        case "jumpType"
            errFld = 11;
            yMax = 30;
            name = "Classification";
    end

    % work with the PAD dataset for this outcome variable
    subset = (C.Model==outcomeVars(i) & C.Norm=='PAD');
    % find the best model
    [~, bestIdx ] = min( table2array(P( subset, errFld )) );
    datasets = P.Dataset( subset );
    bestDataset = datasets( bestIdx );
    
    % display the plot
    ax = subplot( nVars, 2, i*2-1 );
    plotRSqBoxPlot( ax, C( subset, : ), bestDataset, yMax, [i==nVars true] );
    title( ax, strcat( name, ": PAD") );
    text( ax, -0.2, 1.1, ['(' char(64+i*2-1) ')'], 'Units', 'normalized' );

    % work with the LTN dataset for this outcome variable
    subset = (C.Model==outcomeVars(i) & C.Norm=='LTN');
    % find the best model
    [~, bestIdx ] = min( table2array(P( subset, errFld )) );
    datasets = P.Dataset( subset );
    bestDataset = datasets( bestIdx );
    
    % display the plot
    ax = subplot( nVars, 2, i*2 );
    plotRSqBoxPlot( ax, C( subset, : ), bestDataset, yMax, [i==nVars false] );
    title( ax, strcat( name, ": LTN") );
    text( ax, -0.2, 1.1, ['(' char(64+i*2) ')'], 'Units', 'normalized' );

end

end


function plotRSqBoxPlot( ax, data, best, yMax, showLabels )

nAmp = 5;
nPha = 1;

% present spread of values
fields = [ 4:4+nAmp-1 19:19+nPha-1 ];
R = table2array(data(:,fields))*100; % first N components

boxObj = boxchart( ax, R );
boxObj.MarkerStyle = 'x';
boxObj.MarkerSize = 4;

% overlay the best model
hold( ax, 'on' );
X = 1:(nAmp+nPha);
Y = table2array( data( data.Dataset==best, fields) )*100;
scObj = scatter( ax, X, Y, 20, [0.6350 0.0780 0.1840], 'filled', 'o' );
for i = 1:length(X)
    text( ax, X(i)+0.25, Y(i), [num2str( Y(i), '%.1f' ) '%'], ...
        'VerticalAlignment', 'baseline', 'FontSize', 8 );
    if i<length(X)
        ax.XTickLabel{i} = [ 'PCA' num2str(i) ];
    else
        ax.XTickLabel{i} = [ 'PCAW' num2str(i-nAmp) ];
    end
end
hold( ax, 'off' );

name = strcat( "Best: ", extractAfter(best, 7) );
legend( ax, scObj, name, 'Location', 'best' );

ax.FontName = 'Arial';
ax.TickDir = 'out';
ylim( ax, [0 yMax] );
if showLabels(1)
    xlabel( ax, 'Model Predictors' );
end
if showLabels(2)
    ylabel( ax, 'Expl Variance, \DeltaR^{2} (%)' );
end

end

