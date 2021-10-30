% ************************************************************************
% Function: plotComponent
% Purpose:  Plot selected FPC
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

function plotComponent( ax, i, pca, coeff, name, showLabels, showLegend  )


% use the unrotated components
pca = pca.unrotated; 

% do calculations
tSpan = -1100:10:0;

% obtain the mean and the time-varying adjustment for the component
y0 = eval_fd( tSpan, pca.meanfd );
yd = sqrt( pca.values(i) )*eval_fd( tSpan, pca.harmfd(i) );

yPlus = y0+sign(coeff(i))*yd;
yMinus = y0-sign(coeff(i))*yd;

% work out boundary of shaded area
tRev = [ tSpan, fliplr(tSpan) ];
yPlusRev = [ yPlus; flipud(y0) ];
yMinusRev = [ yMinus; flipud(y0) ];
    
% draw shaded regions for plus and minus
hold on;
plotRef(1) = fill( ax, tRev, yPlusRev, 'r', 'FaceAlpha', 0.4, ...
                'DisplayName', '+' );
plotRef(2) = fill( ax, tRev, yMinusRev, 'b', 'FaceAlpha', 0.4, ...
                'DisplayName', '-' );

% draw a border line
plot( ax, tSpan, yPlus, 'r' );
plot( ax, tSpan, yMinus, 'b' );

% draw the mean line
plot( ax, tSpan, y0, 'k' );  

hold off;

if showLegend
    legend( plotRef, 'Location', 'northwest' );
end
    
xlim( ax, [-1100 0] );
ytickformat( ax, '%.1f' );

title( ax, [name num2str(i)], ...
        'FontSize', 9, 'FontWeight', 'normal' );

if showLabels(1)
    xlabel( ax, 'Time (ms)');
end
if showLabels(2)
    ylabel( ax, 'VGRF (BW)');
end

end