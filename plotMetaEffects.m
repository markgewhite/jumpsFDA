% ************************************************************************
% Function: plotMetaEffects
% Purpose:  Create a panel plot of the effects from meta models
%           either 'Features' or 'Regs'
%
% Parameters:
%       metaModels  : fitted meta models
%       type        : 'Features' or 'Regs'
%       scaling     : y-scaling factor (optional)
%       unit        : units to display (char array) 
%
% Outputs:
%       one figure
%
% ************************************************************************

function [figObj, preds] = ...
                plotMetaEffects( metaModels, type, scaling, units )

% parsing
isMulti = iscell(metaModels);
if isMulti
    % multiple models to plot
    nModels = length( metaModels );
    if ~iscell(units) 
        error('No multiple units');
    elseif length(units)~=nModels
        error('Number of units does not match number of models.');
    end

else
    % just one, but make it a cell
    nModels = 1;
    temp = metaModels;
    clear metaModels;
    metaModels = {temp};
    temp = units;
    clear units;
    units = {temp};
end

switch type
    case 'Features'
        nRows = nModels;
        nextModelFcn = @(m) m;
        legendFcn = @(s) [s(1:3) '(' s(4) ')'];
        ordinate = 'Length Normalisation';
    case 'Regs'
        nRows = 2*nModels;
        nextModelFcn = @(m) fix((m+1)/2);
        legendFcn = @(s) [s(1:7) '(' s(8) ')'];
        ordinate = 'Landmark Registration';
end

if nargin<3
    scaling = 1;
end
if nargin<4
    units = repmat("?", nModels);
end


% setup figure
figObj = figure;
ax = gobjects( nRows, 1 );
barObj = gobjects( nRows, 4 );

for i = 1:nRows

    % get meta model predictions
    m = nextModelFcn( i );
    preds = metaModelPredictions( metaModels{m}, type );

    % prepare data for bar chart
    switch type
        case 'Features'
            v1 = 1;
            v2 = 4;
        case 'Regs'
            h = mod(i-1,2);
            v1 = 8*h+1;
            v2 = 8*(h+1);
    end

    vNames = preds.Properties.VariableNames(v1:v2);
    x = categorical( vNames, vNames, 'Ordinal', true );
    y = scaling*table2array(preds)';
    y = y(v1:v2, :);
    z = preds.Properties.RowNames;
    z = cellfun( legendFcn, z, 'UniformOutput', false );
    
    % draw bar chart
    ax(i) = subplot( nRows, 1, i);
    barObj(i,:) = bar(  ax(i), x, y );
    ylabel( ax(i), units(m) );
    if i==1
        xlabel( ax(i), ordinate );
        legend( ax(i), z, 'Location', 'northwest', ...
                          'Orientation', 'horizontal' );
    end

    % finalise
    set( ax(i), 'FontName', 'Arial' );
    set( ax(i), 'FontSize', 9 );
    set( ax(i), 'Box', false );
    set( ax(i), 'TickDir', 'out' );
    

end


end