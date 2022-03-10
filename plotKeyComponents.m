% ************************************************************************
% Function: plotKeyComponents
% Purpose:  Create a panel plot of selected FPCs (hard coded)
%
% Parameters:
%       vgrfFd: array of FPC analyses
%       models: array of models
%
% Outputs:
%       one figure
%
% ************************************************************************

function figObj = plotKeyComponents( vgrfPCA )

figObj = figure;
figObj.Position(2) = figObj.Position(2)-figObj.Position(4);
figObj.Position(3) = figObj.Position(3)*1.5;
figObj.Position(4) = figObj.Position(4)*2;

% constants
woaID = 1;

padID = 1;
ltnID = 2;

nRows = 4;
nCols = 2;
ax = gobjects( nRows*nCols, 1 );
plotRef = gobjects( 4, 1 );

% plot first four unregistered PAD components down left hand column
axes = 1:2:7;
coeff = [ 1 1 1 1 ]; % positive model coefficients for all four components
[ ax(axes), plotRef(1:2) ] = plotFPCSet( [], nRows, nCols, axes, ...
                    'PAD0000-', coeff, vgrfPCA{ woaID, padID, 1, 1 } );
% plot 0001- overlay
[ ~, plotRef(3) ] = plotFPCSet( ax(axes), nRows, nCols, axes, ...
                'PAD0001-', coeff, vgrfPCA{ woaID, padID, 2, 1 }, 'k-' );
% plot 0010- overlay
[ ~, plotRef(4) ] = plotFPCSet( ax(axes), nRows, nCols, axes, ...
                'PAD0010-', coeff, vgrfPCA{ woaID, padID, 3, 1 }, 'r-' );
% display legend
lgd = legend( ax(1), plotRef, 'Location', 'northwest' );
lgl.Box = 'off';


% plot first four unregistered PAD components down right hand column
axes = 2:2:8;
[ ax(axes), plotRef(1:2) ] = plotFPCSet( [], nRows, nCols, axes, ...
                    'LTN0000-', coeff, vgrfPCA{ woaID, ltnID, 1, 1 } );
% plot 0001- overlay
[ ~, plotRef(3) ] = plotFPCSet( ax(axes), nRows, nCols, axes, ...
                'LTN0001-', coeff, vgrfPCA{ woaID, ltnID, 2, 1 }, 'k-' );
% plot 0010- overlay
[ ~, plotRef(4) ] = plotFPCSet( ax(axes), nRows, nCols, axes, ...
                'LTN0010-', coeff, vgrfPCA{ woaID, ltnID, 3, 1 }, 'r-' );  
% display legend
lgd = legend( ax(2), plotRef, 'Location', 'northwest' );
lgl.Box = 'off';     
                    
% sort out size and positioning
%for i = 1:length(ax)
%    subPos = ax(i).Position;
%    vOffset = 0.01*(nCols-ceil(i/nCols));
%    ax(i).Position = [ subPos(1) subPos(2)+vOffset ...
%                                    subPos(3) subPos(4)-0.05 ];
%    ax(i).XAxis.Label.Units = 'normalized';
%    ax(i).YAxis.Label.Units = 'normalized';
%    ax(i).XAxis.Label.Position(2) = -0.25;
%    ax(i).YAxis.Label.Position(1) = -0.3;
%end
                    
% finalise plot
finalisePlot( figObj, [0 2.5], 'northwest' );


end


function [ ax, plotRef ] = plotFPCSet( ax, axRows, axCols, axRng, name, ...
                            coeff, vgrfPCA, lineSpec, useWarp )

if nargin<9
    useWarp = false;
end
if nargin<8
    lineSpec = 'k-';
end

nComp = length(axRng);
newPlot = isempty( ax );
if newPlot
   ax = gobjects( nComp, 1 );
end

for i = 1:nComp
    
    if newPlot
        ax(i) = subplot( axRows, axCols, axRng(i) );
        showLabels = [ axRng(i) > (axRows-1)*axCols ...
                            mod( axRng(i), axCols )==1 ];
    else
        showLabels = [ false false ];
    end
    
    if ~useWarp
        if newPlot
            plotRef = plotComponent( ax(i), i, vgrfPCA, coeff, ...
                            name, showLabels );
        else
            plotRef = plotComponentOverlay(  ax(i), i, vgrfPCA, ...
                            coeff, name, lineSpec );
        end
    else
        plotWarpComponent( ax(i), i, vgrfPCA, coeff, ...
                            name, showLabels, false );
    end
            
end


end

