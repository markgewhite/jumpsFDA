% ************************************************************************
% Function: plotWarpComponent
% Purpose:  Plot selected warp (temporal) FPC
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
% Output:
%       Plots and text
%
% ************************************************************************

function plotWarpComponent2( ax, i, pca, coeff, name, showLabels, showLegend  )

% use the unrotated components
pcaA = pca.unrotated;
pcaW = pca.warp;

% do calculations
tSpan = -1100:10:0;

% get the mean curve's points when even time spacing
t0 = eval_fd( tSpan, pcaW.meanfd );

% get the warp component's SD (variance squared rooted) 
td = sqrt( pcaW.values(i) )*eval_fd( tSpan, pcaW.harmfd(i) );

tPlus = t0+sign(coeff(i))*td;
tMinus= t0-sign(coeff(i))*td;

tRev = [ tSpan, fliplr(tSpan) ];
tPlusRev = [ tPlus; flipud(t0) ];
tMinusRev = [ tMinus; flipud(t0) ];

y0 = eval_fd( t0, pcaA.meanfd );
yPlus = eval_fd( tPlus, pcaA.meanfd );
yMinus = eval_fd( tMinus, pcaA.meanfd );
yPlusRev = eval_fd( tPlusRev, pcaA.meanfd );
yMinusRev = eval_fd( tMinusRev, pcaA.meanfd );

% draw shaded regions for plus and minus
hold on;
plotRef(1) = fill( ax, tRev, yPlusRev, 'b', 'FaceAlpha', 0.4, ...
                'DisplayName', '+' );
plotRef(2) = fill( ax, tRev, yMinusRev, 'r', 'FaceAlpha', 0.4, ...
                'DisplayName', '-' );

% draw a border line
plot( ax, tSpan, yPlus, 'b' );
plot( ax, tSpan, yMinus, 'r' );
% draw the mean line
plot( ax, tSpan, y0, 'k' );  

hold off;

if showLegend
    legend( plotRef, 'Location', 'northwest' );
end

xlim( ax, [-1100 0] );
ytickformat( ax, '%.1f' );

title( ax, [name num2str(i) 'W'], ...
        'FontSize', 9, 'FontWeight', 'normal' );

if showLabels(1)
    xlabel( ax, 'Time (ms)');
end
if showLabels(2)
    ylabel( ax, 'VGRF (BW)');
end

end