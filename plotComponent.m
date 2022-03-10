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
%       showLabels: logical stating whether to show x and y axis labels%       
%
% Output:
%       Plots and text
%
% ************************************************************************

function plotRef = plotComponent( ax, i, pca, coeff, name, showLabels  )

% use the unrotated components
pca = pca.unrotated; 

% do calculations
tSpan = -950:10:0;

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
hold( ax, 'on' );
plotRef(1) = fill( ax, tRev, yPlusRev, 'b', 'FaceAlpha', 0.25, ...
                'EdgeColor', 'none', ...
                'DisplayName', [name ': +ve'] );
plotRef(2) = fill( ax, tRev, yMinusRev, 'b', 'FaceAlpha', 0.50, ...
                'EdgeColor', 'none', ...
                'DisplayName', [name ': -ve'] );

% draw a border line
plot( ax, tSpan, yPlus, 'b' );
plot( ax, tSpan, yMinus, 'b' );  

hold( ax, 'off' );
    
xlim( ax, [-1100 0] );
ytickformat( ax, '%.1f' );

title( ax, [name(1:3) ' FPC' num2str(i)], ...
        'FontSize', 9, 'FontWeight', 'normal' );

if showLabels(1)
    xlabel( ax, 'Time (ms)');
end
if showLabels(2)
    ylabel( ax, 'VGRF (BW)');
end

end