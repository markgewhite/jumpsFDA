% ************************************************************************
% Function: plotDecompScatter
% Purpose:  Create a scatter plot of the registration decompositions
%
% Parameters:
%       decomp : cell array of decomposition analysis structures
%
% Outputs:
%       one figure
%
% ************************************************************************

function fig = plotDecompScatter( decomp )

% initialise
darkBlue = [0.0000 0.4470 0.7410];
lightBlue = [0.3010 0.7450 0.9330];
darkRed = [0.6350 0.0780 0.1840];
lightRed = [0.9290 0.6940 0.1250];
xCol1 = 55;
xCol2 = 90;
xOffset = 4;
fontSize = 8;

labelsShown = {'LTN1110C', 'LTN1100C', 'PAD0001C', 'PAD0000C', ...
                'PAD1110-', 'LTN1110-', 'PAD0001-', 'LTN1000-', ...
                'LTN1111C', 'LTN0001-', 'LTN1111-' };

% create the data set
n = numel( decomp )/2-2;
normStr = {'PAD', 'LTN'};
CTStr = {'-', 'C'};
D.norm = strings( n, 1);
D.LM = strings( n, 1 );
D.CT = strings( n, 1 );
D.pt = zeros( n, 2 );
i = 0;
for j = 1:2 % PAD & LTN
    for l = 1:2 % CT
        for k = 1:16 % LM
            d = decomp{1,j,k,l};
            if ~isempty(d)
                i = i+1;
                D.norm(i) = normStr(j);
                D.LM(i) = dec2bin( k-1, 4 );
                D.CT(i) = CTStr(l);
                D.pt(i,:) = [ d.ampVar d.phaVar ];
            end
        end
    end
end

fig = figure;
ax = gca;
hold( ax, 'on' );
for i = 1:n
    switch D.norm(i)
        case 'PAD'
            marker = 'o'; % circle
            blue = darkBlue;
            red = darkRed;
        case 'LTN'
            marker = 's'; % square
            blue = lightBlue;
            red = lightRed;
    end
    switch D.CT(i)
        case '-'
            colour = blue;
        case 'C'
            colour = red;
    end
    
    sz = 20 + length( strfind(D.LM(i),'1') )*15;
    if marker=='s'
        sz = sz*sqrt(2);
    end

    scatter( ax, D.pt(i,1), D.pt(i,2), sz, colour, 'filled', marker, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.25 );
    
    label = strcat( D.norm(i), D.LM(i), D.CT(i) );
    if any(strcmp( label, labelsShown ))
        text( ax, D.pt(i,1)+xOffset, D.pt(i,2), label, ...
            'FontName', 'Arial', 'FontSize', fontSize );
    end

end

plot( ax, [0 100], [0 0], 'k:', 'LineWidth', 0.5 );

% manual legend
scatter( ax, xCol1, 160, 50*sqrt(2), lightRed, ...
    'filled', 's', 'MarkerEdgeColor', 'k', 'LineWidth', 0.25 );
text( ax, xCol1+xOffset, 160, 'LTN - With CT', ...
    'FontName', 'Arial', 'FontSize', fontSize );
scatter( ax, xCol1, 150, 50, darkRed, ...
    'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 0.25 );
text( ax, xCol1+xOffset, 150, 'PAD - With CT', ...
    'FontName', 'Arial', 'FontSize', fontSize );
scatter( ax, xCol1, 140, 50*sqrt(2), lightBlue, ...
    'filled', 's', 'MarkerEdgeColor', 'k', 'LineWidth', 0.25 );
text( ax, xCol1+xOffset, 140, 'LTN - No CT', ...
    'FontName', 'Arial', 'FontSize', fontSize );
scatter( ax, xCol1, 130, 50, darkBlue, ...
    'filled', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 0.25 );
text( ax, xCol1+xOffset, 130, 'PAD - No CT', ...
    'FontName', 'Arial', 'FontSize', fontSize );

text( ax, xCol2+5, 180, 'No. Landmarks', ...
    'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', fontSize );
scatter( ax, xCol2+5, 170, 20, [0 0 0 ], 'o', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.25 );
text( ax, xCol2, 170, '0', 'FontName', 'Arial', 'FontSize', fontSize );
scatter( ax, xCol2+5, 160, 20+15, [0 0 0 ], 'o', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.25 );
text( ax, xCol2, 160, '1', 'FontName', 'Arial', 'FontSize', fontSize );
scatter( ax, xCol2+5, 150, 20+30, [0 0 0 ], 'o', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.25 );
text( ax, xCol2, 150, '2', 'FontName', 'Arial', 'FontSize', fontSize );
scatter( ax, xCol2+5, 140, 20+45, [0 0 0 ], 'o', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.25 );
text( ax, xCol2, 140, '3', 'FontName', 'Arial', 'FontSize', fontSize );
scatter( ax, xCol2+5, 130, 20+60, [0 0 0 ], 'o', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.25 );
text( ax, xCol2, 130, '4', 'FontName', 'Arial', 'FontSize', fontSize );

hold( ax, 'off');

% formatting
ax.FontSize = fontSize;
ax.FontName = 'Arial';
ax.TickDir = 'out';
xlim( ax, [0 100] );
ylim( ax, [-20 180] );
xlabel( ax, 'Amplitude Variance' );
ylabel( ax, 'Phase Variance' );



end