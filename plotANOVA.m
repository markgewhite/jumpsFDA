% ************************************************************************
% Function: plotANOVA
% Purpose:  Create a panel plot of ANONA interactions
%
% Parameters:
%       data: table used in ANOVA
%
% Outputs:
%       one figure
%
% ************************************************************************

function figObj = plotANOVA( data )

figObj = figure;

% model*method interaction box plot
ax(1) = subplot(3,3,1);
boxObj{1} = boxchart(  ax(1), data.Model, data.Outcome, ...
                            'GroupByColor', data.Method, ...
                            'MarkerStyle', 'none', ...
                            'LineWidth', 0.75 );
xlabel( ax(1), 'Model');
ylabel( ax(1), 'R^{2} / Acc (%)');
legend( ax(1), 'Location', 'best' );

% model*predictors interaction box plot
ax(2) = subplot(3,3,2);
boxObj{2} = boxchart(  ax(2), data.Model, data.Outcome, ...
                            'GroupByColor', data.PredictorSet, ...
                            'MarkerStyle', 'none', ...
                            'LineWidth', 0.75 );
xlabel( ax(2), 'Model');
legend( ax(2), 'Location', 'best' );

% model*partition plot
ax(3) = subplot(3,3,3);
boxObj{3} = boxchart(  ax(3), data.Model, data.Outcome, ...
                            'GroupByColor', data.Set, ...
                            'MarkerStyle', 'none', ...
                            'LineWidth', 0.75 );
xlabel( ax(3), 'Model');
legend( ax(3), 'Location', 'best' );

% method*registration plot
ax(4) = subplot(3,1,2);
boxObj{4} = boxchart(  ax(4), data.Registration, data.Outcome, ...
                            'GroupByColor', data.Method, ...
                            'MarkerStyle', 'none', ...
                            'LineWidth', 0.75 );
xlabel( ax(4), 'Registration');
ylabel( ax(4), 'R^{2} / Acc (%)'); 
legend( ax(4), 'Location', 'best' );

% predictors*registration plot
ax(5) = subplot(3,1,3);
boxObj{5} = boxchart(  ax(5), data.Registration, data.Outcome, ...
                            'GroupByColor', data.PredictorSet, ...
                            'MarkerStyle', 'none', ...
                            'LineWidth', 0.75 );
xlabel( ax(5), 'Registration');
ylabel( ax(5), 'R^{2} / Acc (%)');
legend( ax(5), 'Location', 'best' );

% sort out size and positioning
for i = 1:length(ax)
    subPos = ax(i).Position;
    switch i
        case {1, 2, 3}
            ax(i).Position = [ subPos(1) subPos(2)+0.05 subPos(3) subPos(4)-0.03 ];
            ax(i).YAxis.Label.Position(1) = -0.4;
        case 4
            ax(i).Position = [ subPos(1) subPos(2)+0.025 subPos(3) subPos(4)-0.03 ];
            ax(i).YAxis.Label.Position(1) = -0.5;
        case 5
            ax(i).Position = [ subPos(1) subPos(2)+0.00 subPos(3) subPos(4)-0.03 ];
            ax(i).YAxis.Label.Position(1) = -0.5;
    end
end


% finalise
finalisePlot( figObj, [0 100] );

end