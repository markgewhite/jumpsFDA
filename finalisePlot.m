% ************************************************************************
% Script:   finalisePlot
% Purpose:  Finalise a figure to standard setup
%
%
% ************************************************************************

function finalisePlot( fig, yLimits, legPos )

if nargin<3
    legPos = 'Best';
end

multiLimits = size( yLimits, 1 )>1;

ax = fig.Children;

a = 0;
p = 0;
for i = length(ax):-1:1 % reverse order
    
    
    switch class(ax(i))
        
        case 'matlab.graphics.axis.Axes'
            a = a+1;
            axYLim = yLimits( 1+(a-1)*multiLimits, : );
            ylim( ax(i), axYLim );
            if ax(i).InnerPosition(3) < 2*ax(i).InnerPosition(4)
                labelX = -0.25;
                labelY = 1.15;
            else
                labelX = -0.20;
                labelY = 1.15;
            end
    
            p = p+1;
            text( ax(i), labelX, labelY, ['(' char(64+p) ')'], ...
                        'Units', 'normalized', ...
                        'FontName', 'Arial', ...
                        'FontSize', 9 );

            ax(i).XLabel.Units = 'normalized';
            ax(i).XLabel.Position(2) = -0.22;

            ax(i).YLabel.Units = 'normalized';
            ax(i).YLabel.Position(1) = -0.12;

            ax(i).FontName = 'Arial';
            ax(i).FontSize = 8;
            ax(i).Box = false;
            ax(i).TickDir = 'out';
            
            
        case 'matlab.graphics.illustration.Legend'
            
            ax(i).Location = legPos;
            ax(i).Box = 'off';
            
    end

end