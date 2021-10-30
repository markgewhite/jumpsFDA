% ************************************************************************
% Function: plotTopModels
% Purpose:  Create a panel plot of top models
%
% Parameters:
%       data: table used in ANOVA
%       nTop: number of top models to show
%
% Outputs:
%       one figure
%
% ************************************************************************

function figObj = plotTopModels( data, nTop )

figObj = figure;

% ignore training data
data = data( data.Set=="Test", : );

% calculate the outcome gain from using warp
noReg = (data.Registration=="0000");
data.Amplitude = data.Outcome_NoWarp;
data.Amplitude( noReg ) = data.Outcome( noReg );
data.Warp = data.Outcome - data.Outcome_NoWarp;
data.Warp( noReg ) = 0;

% combine method, predictor set and registration in one field
data.Processing = data.Method.*data.PredictorSet.*data.Registration;

% sort rows in descending order of outcome
data = sortrows( data, { 'Model', 'Outcome' }, ...
                       { 'ascend', 'descend' } );

% plot jump height (take-off velocity)
ax(1) = subplot( 2,2,1 );
plotChart( ax(1), data, 'JHtov', nTop, [false true] );

% plot jump height (work done)
ax(2) = subplot( 2,2,2 );
plotChart( ax(2), data, 'JHwd', nTop, [false false] );

% plot peak power
ax(3) = subplot( 2,2,3 );
plotChart( ax(3), data, 'PP', nTop, [true true] );

% plot classification
ax(4) = subplot( 2,2,4 );
plotChart( ax(4), data, 'Classify', nTop, [true false] );

% sort out size and positioning
for i = 1:length(ax)
    subPos = ax(i).Position;
    ax(i).Position = [ subPos(1) subPos(2)+0.08 subPos(3) subPos(4)-0.07 ];
    ax(i).YAxis.Label.Position(1) = -0.8;
    ax(i).XAxis.Label.Position(2) = -50;
end


finalisePlot( figObj, [0 100] );



end


function plotChart( ax, data, model, nTop, showLabels )

% extract subset
dataM = data( data.Model == model, : );

% select the top models
y = table2array(dataM( 1:nTop, [8 9] ));

% only retain the top categories
x = cellstr( dataM.Processing( 1:nTop ) );
x = categorical( strrep( x, ' ', '-' ) );
x = reordercats( x, string(x) );

% plot the chart

barObj = bar( ax, x, y, 'stacked', 'LineWidth', 0.75 );
barObj.FaceAlpha = 0.4;
barObj.EdgeColor = barObj.FaceColor;

% annotate values
valStr = string(compose( '%.1f%%', barObj(1).YData + barObj(2).YData ));

yTips = repelem( 30, nTop );
xTips = barObj(1).XEndPoints;
text( xTips, yTips, valStr, 'HorizontalAlignment', 'right', ...
                            'VerticalAlignment', 'middle', ...
                            'Rotation', 90, ...
                            'FontSize', 7 );

if showLabels(1) 
    xlabel( ax, 'Data Processing' );
end
if showLabels(2)
    ylabel( ax, 'R^{2} / Acc (%)' ); 
end
title( ax, model );
legend( ax, 'Ampl', 'Warp' );

end

