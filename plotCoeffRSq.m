% ************************************************************************
% Function: plotCoeffRSq
% Purpose:  Create a panel plot of explained variance by coefficient
%
% Parameters:
%       models : cell array of model output structures
%
% Outputs:
%       one figure
%
% ************************************************************************

function fig = plotCoeffRSq( models )

% compile table to work with
T = compileResults( models, 'coeffRSq' );

% create categories
T.Model = categorical( T.Model );
T.PredictorSet = categorical( T.PredictorSet );

usedPAD = extract( T.Dataset, 3 )=='1';
T.Norm( usedPAD ) = "PAD";
T.Norm( ~usedPAD ) = "LTN";
T.Norm = categorical( T.Norm );

T.Dataset = categorical( T.Dataset );

% retain a subset
retain1 = T.Model=='JHtov' | T.Model=='PP' | T.Model=='jumpType';
retain2 = T.PredictorSet=='PCAU';
T = T( retain1 & retain2, : );

% create figure
fig = figure;
fig.Position(2) = fig.Position(2) - fig.Position(4)*0.5;
fig.Position(3) = fig.Position(3)*1.75;
fig.Position(4) = fig.Position(4)*1.5;

% JH-PAD panel 
ax = subplot( 3, 2, 1 );
subset = (T.Model=='JHtov' & T.Norm=='PAD');
plotRSqBoxPlot( ax, T( subset, : ), "1-1VGRF0001-", 50, [false true] );
title( ax, 'JH: PAD' );

% JH-LTN panel 
ax = subplot( 3, 2, 2 );
subset = (T.Model=='JHtov' & T.Norm=='LTN');
plotRSqBoxPlot( ax, T( subset, : ), "1-2VGRF0001-", 50, [false false] );
title( ax, 'JH: LTN' );

% PP-PAD panel 
ax = subplot( 3, 2, 3 );
subset = (T.Model=='PP' & T.Norm=='PAD');
plotRSqBoxPlot( ax, T( subset, : ), "1-1VGRF0000-", 70, [false true] );
title( ax, 'PP: PAD' );

% PP-LTN panel 
ax = subplot( 3, 2, 4 );
subset = (T.Model=='PP' & T.Norm=='LTN');
plotRSqBoxPlot( ax, T( subset, : ), "1-2VGRF0000-", 70, [false false] );
title( ax, 'PP: LTN' );

% CL-PAD panel 
ax = subplot( 3, 2, 5 );
subset = (T.Model=='jumpType' & T.Norm=='PAD');
plotRSqBoxPlot( ax, T( subset, : ), "2-1VGRF0010-", 40, [true true] );
title( ax, 'CL: PAD' );

% CL-LTN panel 
ax = subplot( 3, 2, 6 );
subset = (T.Model=='jumpType' & T.Norm=='LTN');
plotRSqBoxPlot( ax, T( subset, : ), "2-2VGRF0010-", 40, [true false] );
title( ax, 'CL: LTN' );

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
    ylabel( ax, 'Expl Var, \DeltaR^{2} (%)' );
end

end

