% ************************************************************************
% Function: plotLLEffects
% Purpose:  Create a panel plot of the effects on log likelihood
%           of different landmark registrations. 
%           The data must come from SAS's meta model (GLIMMIX)
%           a non-linear mixed model, too advanced for MATLAB to do
%
% Parameters:
%       data: table used in ANOVA
%
% Outputs:
%       one figure
%
% ************************************************************************

function figObj = plotLLEffects( data )

sigLevel = 80.2; % Chi-squared for df=61, p=0.05

data.Effect = categorical( data.Effect );

colours = [ 0.0000 0.4470 0.7410; ...
            0.6350 0.0780 0.1840; ...
            0.4660 0.6740 0.1880 ];

figObj = figure;

% LTN bar chart
filter = (data.Effect=='Norm');
x = 1:sum(filter);
ax(1) = subplot(5,6,1);
barObj{1} = bar(  ax(1), x, table2array(data(filter, 6:8)) );
ax(1).XTickLabel = 'LTN';

% CT bar chart
filter = (data.Effect=='CTReg') & (data.JH_Estimate~=0);
x = 1:sum(filter);
ax(2) = subplot(5,6,2);
barObj{2} = bar(  ax(2), x, table2array(data(filter, 6:8)) );
ax(2).XTickLabel = 'CT';
                     
% ACP bar chart
filter = (data.Effect=='Predictor') & (data.JH_Estimate~=0);
ax(3) = subplot(5,6,3);
barObj{3} = bar(  ax(3), x, table2array(data(filter, 6:8)) );
ax(3).XTickLabel = 'ACP';                    

% LTN*CT bar chart
filter = (data.Effect=='Norm*CTReg') & (data.JH_Estimate~=0);
ax(4) = subplot(5,6,4);
barObj{4} = bar(  ax(4), x, table2array(data(filter, 6:8)) );
ax(4).XTickLabel = 'LTN\timesCT';

% LTN*ACP bar chart
filter = (data.Effect=='Norm*Predictor') & (data.JH_Estimate~=0);
ax(5) = subplot(5,6,5);
barObj{5} = bar(  ax(5), x, table2array(data(filter, 6:8)) );
ax(5).XTickLabel = 'LTN\timesACP';

% CT*ACP bar chart
filter = (data.Effect=='CTReg*Predictor') & (data.JH_Estimate~=0);
ax(6) = subplot(5,6,6);
barObj{6} = bar(  ax(6), x, table2array(data(filter, 6:8)) );
ax(6).XTickLabel = 'CT\timesACP';

% LM box plot
filter = (data.Effect=='LMReg') & (data.JH_Estimate~=0);
x = 1:sum(filter);
reg = categorical(data.LMReg( filter ));
ax(7) = subplot(5,1,2);
barObj{7} = bar(  ax(7), x, table2array(data(filter, 6:8)) );
ax(7).XTickLabel = reg;
xlabel( ax(7), 'LM' );

% LTN*LM box plot
filter = (data.Effect=='Norm*LMReg') & (data.JH_Estimate~=0);
ax(8) = subplot(5,1,3);
barObj{8} = bar(  ax(8), x, table2array(data(filter, 6:8)) );
ax(8).XTickLabel = reg;
xlabel( ax(8), 'LM \times LTN' );

% CT*LM box plot
filter = (data.Effect=='LMReg*CTReg') & (data.JH_Estimate~=0);
ax(9) = subplot(5,1,4);
barObj{9} = bar(  ax(9), x, table2array(data(filter, 6:8)) );
ax(9).XTickLabel = reg;
xlabel( ax(9), 'LM \times CT' );

% ACP*LM box plot
filter = (data.Effect=='LMReg*Predictor') & (data.JH_Estimate~=0);
ax(10) = subplot(5,1,5);
barObj{10} = bar(  ax(10), x, table2array(data(filter, 6:8)) );
ax(10).XTickLabel = reg;
xlabel( ax(10), 'LM \times ACP' );

% format and all charts
for i = 1:length(ax)
    
    % significance threshold
    hold( ax(i), 'on' );
    if i<=6 
        w = 2;
    else
        w = length(reg)+1;
    end
    plot( ax(i), [0 w], [ sigLevel sigLevel ], 'k--', 'LineWidth', 0.5);
    plot( ax(i), [0 w], [ -sigLevel -sigLevel ], 'k--', 'LineWidth', 0.5);
    hold( ax(i), 'off' );
    xlim( ax(i), [0.5 w-0.5] );
    
    % lines and colours
    for j = 1:3
        barObj{i}(j).CData = colours(j,:);
    end
       
    % axis labels
    if i==1 || i>=7 
        ylabel( ax(i), '-2LL');
        yticks( ax(i), [-400 0 400] );
    end
    if i>=2 && i<=6
        yticks( ax(i), [] );
    end
    
    % legend
    if i==10
        legend( ax(i), 'JH', 'PP', 'CL', '\itp < 0.05', ...
                       'Location', 'southeast', ...
                       'Orientation', 'horizontal', ...
                       'Box', 'off' );
    end
    
end


% finalise
finalisePlot( figObj, [-400 400] );

end


function x = shiftOrder( x )

x = [ x(end,:); x(1:end-1,:) ];

end