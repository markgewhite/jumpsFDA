% ************************************************************************
% Function: outputTable
% Purpose:  Generate the results table for a given data set
%
% Parameters:
%       dataset: title of this data set
%       ref: identifiers
%       withArms: jump type
%       jumpperf: jump performances structure
%       pca: pca structure
%       acp: acp structure
%
% Outputs:
%       results: table with results in wide format
%
% ************************************************************************


function results = outputTable( dataset, ref, ...
                                withArms, perf, ...
                                pca, acp, ...
                                setup )
               
N = size( ref, 1 );
isWarp = isfield( pca, 'warp' );

nComp = setup.pca.nComp;
nWarp = setup.pca.nCompWarp;

nVar = 6 + (nComp+nWarp)*4;
data = zeros( N, nVar );

% build output table

% classifiers
data( :, 1:2 ) = ref; % subject ID + trial number
data( :, 3 ) = withArms; % jump type

% performance metrics
data( :, 4 ) = perf.JHtov; % jump height (take-off velocity)
data( :, 5 ) = perf.JHwd; % jump height (work done)
data( :, 6 ) = perf.PP; % peak power


% PCA components (allowing possibility that they may be reduced)
c(1) = 7;
data( :, c(1):c(1)+nComp-1 ) = pca.unrotated.harmscr(:,1:nComp);
c(2) = c(1)+nComp;
if isWarp
    data( :, c(2):c(2)+nWarp-1 ) = pca.warp.harmscr(:,1:nWarp);
end
c(3) = c(2)+nWarp;

data( :, c(3):c(3)+nComp-1 ) = pca.varimax.harmscr(:,1:nComp);
c(4) = c(3)+nComp;
if isWarp
    data( :, c(4):c(4)+nWarp-1 ) = pca.warpVarimax.harmscr(:,1:nWarp);
end
c(5) = c(4)+nWarp;


% ACP components
data( :, c(5):c(5)+nComp-1 ) = acp.unrotated(:,1:nComp);
c(6) = c(5)+nComp;
if isWarp
    data( :, c(6):c(6)+nWarp-1 ) = acp.warp(:,1:nWarp);
end
c(7) = c(6)+nWarp;

data( :, c(7):c(7)+nComp-1 ) = acp.varimax(:,1:nComp);
c(8) = c(7)+nComp;
if isWarp
    data( :, c(8):c(8)+nWarp-1 ) = acp.warpVarimax(:,1:nWarp);
end


% convert numeric array into results table

results = array2table(data);
results = [ table(repmat( dataset, N, 1 )) results ];

results.Properties.VariableNames{1} = 'Dataset';
results.Properties.VariableNames{2} = 'SubjectID';
results.Properties.VariableNames{3} = 'Trial';
results.Properties.VariableNames{4} = 'Arms';
results.Properties.VariableNames{5} = 'JH_TOV';
results.Properties.VariableNames{6} = 'JH_WD';
results.Properties.VariableNames{7} = 'PP';

for i = 1:nComp
    results.Properties.VariableNames{ c(1)+i } = char("PCA_U"+i);
    results.Properties.VariableNames{ c(3)+i } = char("PCA_V"+i);
    results.Properties.VariableNames{ c(5)+i } = char("ACP_U"+i);
    results.Properties.VariableNames{ c(7)+i } = char("ACP_V"+i);
end

for i = 1:nWarp
    results.Properties.VariableNames{ c(2)+i } = char("PCA_UW"+i);
    results.Properties.VariableNames{ c(4)+i } = char("PCA_VW"+i);
    results.Properties.VariableNames{ c(6)+i } = char("ACP_UW"+i);
    results.Properties.VariableNames{ c(8)+i } = char("ACP_VW"+i);
end


end
