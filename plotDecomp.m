% ************************************************************************
% Function: plotDecomp
% Purpose:  Create a panel plot of decomposition RSq results
%
% Parameters:
%       data: table of decomposition results
%
% Outputs:
%       one figure
%
% ************************************************************************

function figObj = plotDecomp( data )

% extract data
T = table( data.Registration( data.Method=="PAD" ), ...
                        'VariableNames', {'Registration'} );
T.RSqPAD = data.rSq( data.Method=="PAD" )*100;
T.RSqLTN = data.rSq( data.Method=="LTN" )*100;

figObj = figure;

blueColour = [0 0.45 0.74 ]; % blue for PAD decompositions
redColour = [0.85 0.33 0.10 ]; % red for LTN decompositions

% create box chart
ax = subplot(1,1,1);
barObj = bar( ax, T.Registration, ...
                [ T.RSqPAD T.RSqLTN ], ...
                'BarWidth', 1, 'LineWidth', 0.75 );
barObj(1).EdgeColor = blueColour;
barObj(2).EdgeColor = redColour;
barObj(1).FaceColor = blueColour;
barObj(2).FaceColor = redColour;
barObj(1).FaceAlpha = 0.2;
barObj(2).FaceAlpha = 0.2;

ylabel( ax, 'Proportion Phase Variance (%)');
xlabel( ax, 'Registration' );
    
legend( ax, 'PAD', 'LTN' );

% finalise
finalisePlot( figObj, [0 100], 'northwest' );

% re-align the registration ticks
ax.XTickLabelRotation = 45;


end