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
        width = 9.0;
        height = 10.0;
        legendFcn = @(s) [s(1:3) '(' s(4) ')'];
        legendLoc = 'northwest';
        legendCols = 1;
        ordinate = 'Length Normalisation';
    case 'Regs'
        width = 19.0;
        height = 12.0;
        legendFcn = @(s) [s(1:7) '(' s(8) ')'];
        legendLoc = 'north';
        legendCols = 4;
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
figObj.Units = 'centimeters';
figObj.Position(3) = width;
figObj.Position(4) = height;

ax = gobjects( nModels, 1 );
barObj = gobjects( nModels, 4 );

for i = 1:nModels

    % get meta model predictions
    preds = metaModelPredictions( metaModels{i}, type );

    vNames = preds.Properties.VariableNames;
    x = categorical( vNames, vNames, 'Ordinal', true );
    y = scaling*table2array(preds)';
    z = preds.Properties.RowNames;
    z = cellfun( legendFcn, z, 'UniformOutput', false );
    
    % draw bar chart
    ax(i) = subplot( nModels, 1, i);
    barObj(i,:) = bar(  ax(i), x, y );
    ylabel( ax(i), units(i) );
    if i==1
        legend( ax(i), z, 'Location', legendLoc );
        ax(i).Legend.NumColumns = legendCols;
        ax(i).Legend.Position(2) = 0.95;
    end
    if i==nModels
        xlabel( ax(i), ordinate );
    end

    % finalise
    set( ax(i), 'YGrid', 'on' );
    set( ax(i), 'FontName', 'Arial' );
    set( ax(i), 'FontSize', 9 );
    set( ax(i), 'Box', false );
    set( ax(i), 'TickDir', 'out' );
    

end


end