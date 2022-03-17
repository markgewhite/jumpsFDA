% ************************************************************************
% Function: plotComponentOverlay
% Purpose:  Plot selected FPC as an overlay
%
% Parameters:
%       ax: figure axes
%       i: component number
%       pca: functional FPCA data structure
%       coeff: model coefficients (to provide sign)
%       name: title for the plot
%       showLabels: logical stating whether to show x and y axis labels
%       showLegend: whether to show the legend
%       
%
% Output:
%       Plots and text
%
% ************************************************************************

function plotRef = plotComponentOverlay( ax, i, pca, coeff, name, lineSpec )

% use the unrotated components
pca = pca.unrotated; 

% do calculations
tSpan = -1000:5:0;

% obtain the mean and the time-varying adjustment for the component
y0 = eval_fd( tSpan, pca.meanfd );
yd = sqrt( pca.values(i) )*eval_fd( tSpan, pca.harmfd(i) );

yPlus = y0+sign(coeff(i))*yd;
yMinus = y0-sign(coeff(i))*yd;
    
% draw mean and border lines
hold( ax, 'on' );
plot( ax, tSpan, y0, lineSpec, 'LineWidth', 0.5 );
plot( ax, tSpan, yPlus, lineSpec, 'LineWidth', 1.0 );
plotRef = plot( ax, tSpan, yMinus, lineSpec, 'LineWidth', 1.0, ...
            'DisplayName', name );
hold( ax, 'off' );

end