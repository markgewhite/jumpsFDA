% ************************************************************************
% Function: plotCurves
% Purpose:  Create a panel plot of VGRF curves with different registrations
%
% Parameters:
%       vgrfFd : cell array of all registrations
%       warpFd : cell array of matching warping functions
%
% Outputs:
%       one figure
%
% ************************************************************************

function [ fig1, fig2 ] = plotCurves( vgrfFd, warpFd )

% plot the VGRF curves
fig1 = figure;
fig1.Units = 'centimeters';
fig1.Position(3) = 18.0;
fig1.Position(4) = 10.0;

ax = subplot( 2, 3, 1 );
fd = vgrfFd{ 1, 1, 1, 1 }; % WOA-PAD0000-
plotCurveSet( ax, fd, [ false true ], 1 );
title( ax, 'PAD0000-' );

ax = subplot( 2, 3, 2 );
fd = vgrfFd{ 1, 1, 3, 1 }; % WOA-PAD0010-
plotCurveSet( ax, fd, [ false false ], 2 );
title( ax, 'PAD0010-' );

ax = subplot( 2, 3, 3 );
fd = vgrfFd{ 1, 1, 3, 2 }; % WOA-PAD0010C
plotCurveSet( ax, fd, [ false false ], 3 );
title( ax, 'PAD0010C' );

ax = subplot( 2, 3, 4 );
fd = vgrfFd{ 1, 2, 1, 1 }; % WOA-LTN0000-
plotCurveSet( ax, fd, [ true true ], 4 );
title( ax, 'LTN0000-' );

ax = subplot( 2, 3, 5 );
fd = vgrfFd{ 1, 2, 3, 1 }; % WOA-LTN0010-
plotCurveSet( ax, fd, [ true false ], 5 );
title( ax, 'LTN0010-' );

ax = subplot( 2, 3, 6 );
fd = vgrfFd{ 1, 2, 3, 2 }; % WOA-LTN0010C
plotCurveSet( ax, fd, [ true false], 6 );
title( ax, 'LTN0010C' );


% plot the warp curves
fig2 = figure;
fig2.Units = 'centimeters';
fig2.Position(3) = 18.0;
fig2.Position(4) = 5.0;

ax = subplot( 1, 4, 1 );
fd = warpFd{ 1, 1, 3, 1 }; % WOA-PAD0010- warp
plotCurveSet( ax, fd, [ true true ], 1 );
title( ax, 'PAD0010- (Warp)' );

ax = subplot( 1, 4, 2 );
fd = warpFd{ 1, 1, 3, 2 }; % WOA-PAD0010C warp
plotCurveSet( ax, fd, [ true false ], 2 );
title( ax, 'PAD0010C (Warp)' );

ax = subplot( 1, 4, 3 );
fd = warpFd{ 1, 2, 3, 1 }; % WOA-LTN0010- warp
plotCurveSet( ax, fd, [ true false ], 3 );
title( ax, 'LTN0010- (Warp)' );

ax = subplot( 1, 4, 4 );
fd = warpFd{ 1, 2, 3, 2 }; % WOA-LTN0010C warp
plotCurveSet( ax, fd, [ true false ], 3 );
title( ax, 'LTN0010C (Warp)' );

end


function plotCurveSet( ax, fd, showLabel, p )

%  extract basis information
names = getnames( fd );
isVGRF = (strcmp(names{3}(end-7:end), 'GRF (BW)' ));

basis = getbasis( fd );
tRng = getbasisrange( basis );
if isVGRF
    tRng(1) = max( -1500, tRng(1) );
end

nBasis = getnbasis( basis );
nPts = max( [101, 10*nBasis+1] );
tSpan = linspace( tRng(1), tRng(2), nPts )';

% compute points for all curves
y = eval_fd( tSpan, fd );

% plot the curves
nCurves = size( y, 2 );
hold( ax, 'on' );
for i = 1:nCurves
    lineObj = plot( ax, tSpan, y(:,i), 'LineWidth', 0.25 );
    lineObj.Color(4) = 0.6; % transparency
end

if isVGRF
    % plot the mean
    plot( ax, tSpan, mean(y,2), 'k', 'LineWidth', 2 );
else
    % plot the diagonal/fixed time line
    plot( ax, tSpan, tSpan, '--k', 'LineWidth', 2 );
end

hold( ax, 'off' );
xlim( ax, tRng );
if isVGRF
    ax.YLim = [0 3];
    ax.YTick = [0 1 2 3];
else
    ax.YLim = tRng;
    axis( ax, 'square' );
    ax.YTick = ax.XTick;
end
if showLabel(1)
    xlabel( ax, 'Time (ms)' );
end
if showLabel(2)
    if isVGRF
        ylabel( ax, 'VGRF (BW)' );
        ax.YAxis.TickLabelFormat = '%.1f';
        ax.YAxis.Label.Units = 'normalized';
        ax.YAxis.Label.Position(1) = -0.20;
    else
        ylabel( ax, 'Warp Time (ms)' );
        ax.YAxis.Label.Units = 'normalized';
        ax.YAxis.Label.Position(1) = -0.30;
    end
end

ax.TickDir = 'out';
ax.FontName = 'Arial';
ax.FontSize = 8;

text( ax, -0.25, 1.15, ['(' char(64+p) ')'], ...
        'Units', 'normalized', 'FontSize', 9 );

end