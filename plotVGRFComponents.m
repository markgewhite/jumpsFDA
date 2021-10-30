% ************************************************************************
% Function: plotVGRFComponents
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

function figObj = plotVGRFComponents( vgrfPCA, models, figName )

figObj = figure;

% constants
allID = 1;
woaID = 2;

padID = 1;
llnID = 2;

unregID = 1;
regCTID = 2;

switch figName

    case 'PAD-LLN' % PAD-LLN comparison
        nRows = 2;
        nCols = 4;
        ax = gobjects( nRows*nCols, 1 );
        vSizeAdj = -0.05;

        % plot unregistered FPCs for WOA and PAD
        coeff = getCoeff( models{ woaID, padID, unregID }, 'PP', false );
        ax(1:4) = plotFPCSet( nRows, nCols, 1:4, 'PAD', ...
                                coeff, vgrfPCA{ woaID, padID, unregID } );

        % plot unregistered FPCs for WOA and LLN
        coeff = getCoeff( models{ woaID, llnID, unregID }, 'PP', false );
        ax(5:8) = plotFPCSet( nRows, nCols, 5:8, 'LLN', ...
                                coeff, vgrfPCA{ woaID, llnID, unregID } );

                            
    case 'CTReg' % continuous registered comparison
        nRows = 4;
        nCols = 4;
        ax = gobjects( nRows*nCols, 1 );
        vSizeAdj = -0.02;
        
        % plot registered amplitude FPCs for WOA and PAD
        coeff = getCoeff( models{ woaID, padID, regCTID }, 'PP', false );
        ax(1:4) = plotFPCSet( nRows, nCols, 1:4, 'PAD', ...
                        coeff, vgrfPCA{ woaID, padID, regCTID } );
                            
        % plot registered warp FPCs for WOA and PAD
        coeff = getCoeff( models{ woaID, padID, regCTID }, 'PP', true );
        ax(5:8) = plotFPCSet( nRows, nCols, 5:8, 'PAD', ...
                        coeff, vgrfPCA{ woaID, padID, regCTID }, true );

        % plot registered amplitude FPCs for WOA and LLN
        coeff = getCoeff( models{ woaID, llnID, regCTID }, 'PP', false );
        ax(9:12) = plotFPCSet( nRows, nCols, 9:12, 'LLN', ...
                        coeff, vgrfPCA{ woaID, llnID, regCTID } );
                            
        % plot registered warp FPCs for WOA and LLN
        coeff = getCoeff( models{ woaID, llnID, regCTID }, 'PP', true );
        ax(13:16) = plotFPCSet( nRows, nCols, 13:16, 'LLN', ...
                        coeff, vgrfPCA{ woaID, llnID, regCTID }, true );
                            
                            
    case 'BestPAD' % best PAD registrations
        nRows = 4;
        nCols = 4;
        ax = gobjects( nRows*nCols, 1 );
        vSizeAdj = -0.02;
        
        % plot registered amplitude FPCs for WOA and PAD-0001
        coeff = getCoeff( models{ woaID, padID, 3 }, 'PP', false );
        ax(1:2) = plotFPCSet( nRows, nCols, 1:2, 'PAD-0001-', ...
                        coeff, vgrfPCA{ woaID, padID, 3 } );
                            
        % plot registered warp FPCs for WOA and PAD-0001
        coeff = getCoeff( models{ woaID, padID, 3 }, 'PP', true );
        ax(3:4) = plotFPCSet( nRows, nCols, 3:4, 'PAD-0001-', ...
                        coeff, vgrfPCA{ woaID, padID, 3 }, true );

                    
        % plot registered amplitude FPCs for WOA and PAD-0010
        coeff = getCoeff( models{ woaID, padID, 4 }, 'PP', false );
        ax(5:6) = plotFPCSet( nRows, nCols, 5:6, 'PAD-0010-', ...
                        coeff, vgrfPCA{ woaID, padID, 4 } );
                            
        % plot registered warp FPCs for WOA and PAD-0010
        coeff = getCoeff( models{ woaID, padID, 4 }, 'PP', true );
        ax(7:8) = plotFPCSet( nRows, nCols, 7:8, 'PAD-0010-', ...
                        coeff, vgrfPCA{ woaID, padID, 4 }, true );
                    
        
        % plot registered amplitude FPCs for WOA and PAD-0100
        coeff = getCoeff( models{ woaID, padID, 6 }, 'PP', false );
        ax(9:10) = plotFPCSet( nRows, nCols, 9:10, 'PAD-0100-', ...
                        coeff, vgrfPCA{ woaID, padID, 6 } );
                            
        % plot registered warp FPCs for WOA and and PAD-0100
        coeff = getCoeff( models{ woaID, padID, 6 }, 'PP', true );
        ax(11:12) = plotFPCSet( nRows, nCols, 11:12, 'PAD-0100-', ...
                        coeff, vgrfPCA{ woaID, padID, 6 }, true );
                         
        
        % plot registered amplitude FPCs for WOA and PAD-1001
        coeff = getCoeff( models{ woaID, padID, 11 }, 'PP', false );
        ax(13:14) = plotFPCSet( nRows, nCols, 13:14, 'PAD-0100-', ...
                        coeff, vgrfPCA{ woaID, padID, 11 } );
                            
        % plot registered warp FPCs for WOA and and PAD-1001
        coeff = getCoeff( models{ woaID, padID, 11 }, 'PP', true );
        ax(15:16) = plotFPCSet( nRows, nCols, 15:16, 'PAD-0100-', ...
                        coeff, vgrfPCA{ woaID, padID, 11 }, true );
                    
                    
                    
    case 'UserSelect' % user selected
        nRows = 2;
        nCols = 4;
        ax = gobjects( nRows*nCols, 1 );

        while true

            i = input('Registration ID = ');
            regStr = dec2bin( i-2, 4 );
            
            figure;
            
            % plot registered amplitude FPCs for WOA
            coeff = getCoeff( models{ woaID, padID, i }, 'PP', false );
            ax(1:4) = plotFPCSet( nRows, nCols, 1:4, ['PAD' regStr '-'], ...
                            coeff, vgrfPCA{ woaID, padID, i } );
                            
            % plot registered warp FPCs for WOA and PAD-0011
            coeff = getCoeff( models{ woaID, padID, i }, 'PP', true );
            ax(5:8) = plotFPCSet( nRows, nCols, 5:8, ['PAD' regStr '-'], ...
                            coeff, vgrfPCA{ woaID, padID, i }, true );
                        
        end
                   
end
        
                    
% sort out size and positioning
for i = 1:length(ax)
    subPos = ax(i).Position;
    vOffset = 0.01*(nCols-ceil(i/nCols));
    ax(i).Position = [ subPos(1) subPos(2)+vOffset ...
                                    subPos(3) subPos(4)+vSizeAdj ];
    ax(i).XAxis.Label.Units = 'normalized';
    ax(i).YAxis.Label.Units = 'normalized';
    ax(i).XAxis.Label.Position(2) = -0.25;
    ax(i).YAxis.Label.Position(1) = -0.3;
end
                    
% finalise plot
finalisePlot( figObj, [0 2.5], 'northwest' );


end


function ax = plotFPCSet( axRows, axCols, axRng, name, ...
                            coeff, vgrfPCA, useWarp )

if nargin<7
    useWarp = false;
end

nComp = length(axRng);
ax = gobjects( nComp, 1 );
for i = 1:nComp
    
    showLegend = (axRng(i)==1);
    showLabels = [ axRng(i) > (axRows-1)*axCols ...
                    mod( axRng(i), axCols )==1 ];
                
    ax(i) = subplot( axRows, axCols, axRng(i) );
    
    if ~useWarp
        plotComponent( ax(i), i, vgrfPCA, coeff, ...
                            name, showLabels, showLegend );
    else
        plotWarpComponent( ax(i), i, vgrfPCA, coeff, ...
                            name, showLabels, showLegend );
    end
            
end


end


function c = getCoeff( models, modelType, isWarp )

% extract model coefficients
cTbl = models.coeff;
m = ((cTbl.Model==modelType) & (cTbl.PredictorSet=="PCAU"));

% find the block of warp coefficients
for i = 4:size( cTbl, 2 )
    if cTbl.Properties.VariableNames{i}(6) == 'W'
        break
    end
end

% extract the relevant block
if isWarp
    c = table2array( cTbl(m,i:end) );
else
    c = table2array( cTbl(m,4:i-1) );
end

end