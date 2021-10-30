% ************************************************************************
% Function: plotRSq
% Purpose:  Create a panel plot of RSq results
%
% Parameters:
%       data: table used in ANOVA
%
% Outputs:
%       one figure
%
% ************************************************************************

function figObj = plotRSq( data, method )

if nargin > 1
    % retain only specific method
    data = data( data.Method==method, : );
end

% split the data into regression and classification
data1 = data( ismember(data.Model, {'JHtov', 'JHwd', 'PP' }), : );
data2 = data( data.Model=='Classify', : );



figObj = figure;

c = 3;
colour = [0 0 1]; % blue for amplitude components
ax = gobjects( 8, 1 );
barObj = gobjects( 8, 2 );
yLimits = zeros( 8, 2 );
for i = 1:8
    
    if i~=6
        c = c+1; % next component
    else
        c = c+16; % switch to first warp component
        colour = [1 0 0]; % red for warp components
    end
    
    % summarise the regression data
    data1Comp = groupsummary( data1, 'Registration', {'Mean', 'std'}, ...
                             data1.Properties.VariableNames{c} );
    % add in the classification data
    data1Comp = [ data1Comp data2(:,c) ];
    data1Comp.Properties.VariableNames{3} = 'RegMean';
    data1Comp.Properties.VariableNames{4} = 'RegSD';
    data1Comp.Properties.VariableNames{5} = 'Classify';
    
    % create box chat for registration
    ax(i) = subplot(4,2,i);
    barObj(i,:) = bar(  ax(i), data1Comp.Registration, ...
                               table2array(data1Comp(:,[3 5])), ...
                               'FaceColor', colour, ...
                               'BarWidth', 1 );
    barObj(i,1).FaceAlpha = 0.6;
    barObj(i,2).FaceAlpha = 0.2;

    
    % add error bars for regression
    hold on;
    errObj(i) = errorbar( barObj(i,1).XEndPoints, ...
                          barObj(i,1).YEndPoints, data1Comp.RegSD, ...
                          'LineStyle', 'none', ...
                          'Color', [0 0 0], ...
                          'CapSize', 2 );
    hold off;
     
    if mod(i,2)==1 %odd - show label on left hand side only
        ylabel( ax(i), '\DeltaR^{2} (%)');
    end
    if i>=7 % add x axis labels for the bottom row
        xlabel( ax(i), 'Registration' );
    end
    if i<=5 % amplitude component
        title( ax(i), ['Ampl' num2str(i)] );
        yLimits(i,:) = [0 50];
    else
        title( ax(i), ['Warp' num2str(i-5)] );
        yLimits(i,:) = [0 10];
    end
    if i==7
        % special annotation
        text( ax(i), barObj(i,2).XEndPoints(end)-0.3, 9, ...
                     num2str( data1Comp.Classify(end), '%.1f%%' ), ...
                     'HorizontalAlignment', 'right', ...
                     'FontName', 'Arial', ...
                     'FontSize', 8 );
    end
    
    legend( 'Regression', 'Classification' );
    
end


% sort out size and positioning
for i = 1:8
    subPos = ax(i).Position;
    vOffset = (ceil((9-i)/2)-1)*0.015;
    ax(i).Position = [ subPos(1) subPos(2)+vOffset subPos(3)+0.04 subPos(4)-0.02 ];
    ax(i).XAxis.Label.Units = 'normalized';
    ax(i).YAxis.Label.Units = 'normalized';
    ax(i).XAxis.Label.Position(2) = -0.5;
    ax(i).YAxis.Label.Position(1) = -0.1;
end


% finalise
finalisePlot( figObj, yLimits, 'northwest' );

% re-align the registration ticks
for i = 1:8
    ax(i).XTickLabelRotation = 90;
end


end