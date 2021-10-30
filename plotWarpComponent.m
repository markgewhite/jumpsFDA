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

function plotWarpComponent( ax, i, pca, coeff, name, showLabels, showLegend  )

% use the unrotated components
pcaA = pca.unrotated;
pcaW = pca.warp;

% do calculations
tSpan = -1100:10:0;

% get the mean curve's points when even time spacing
y0 = eval_fd( tSpan, pcaA.meanfd );
t0 = eval_fd( tSpan, pcaW.meanfd );

% get the warp component's SD (variance squared rooted) 
td = sqrt( pcaW.values(i) )*eval_fd( tSpan, pcaW.harmfd(i) );

% calculate the adjustment to the mean 
% note: polarity reversed because time points are adjusted
%       if y points were adjusted not such reverse required
%        - see plotWarpComponent2.m
%       but that requires evaluating fd outside of the defined range
tPlus = t0-sign(coeff(i))*td;
tMinus= t0+sign(coeff(i))*td;

% work out boundary of shaded area
yRev = [ y0; flipud(y0) ];
tPlusRev = [ tPlus; flipud(t0) ];
tMinusRev = [ tMinus; flipud(t0) ];
    
% draw shaded regions for plus and minus
hold on;
plotRef(1) = fill( ax, tPlusRev, yRev, 'b', 'FaceAlpha', 0.4, ...
                'DisplayName', '+' );
plotRef(2) = fill( ax, tMinusRev, yRev, 'r', 'FaceAlpha', 0.4, ...
                'DisplayName', '-' );

% draw a border line
plot( ax, tPlus, y0, 'b' );
plot( ax, tMinus, y0, 'r' );

% draw the mean line
plot( ax, t0, y0, 'k' );  

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