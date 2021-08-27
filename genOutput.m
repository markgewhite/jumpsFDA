% ************************************************************************
% Function: genOutput
% Purpose:  Generate the output table for SAS
%
% Parameters:
%       name: title of this data set
%       vgrfPCA: pca structure
%       vgrfACP: acp structure
%       refSet: identifiers
%       withArms: jump type (logical)
%       jumpperf: jump performances structure
%
% Output:
%       output: formatted table
%
% ************************************************************************

function output = genOutput( name, ...
                             refSet, withArms, perf, ...
                             vgrfPCA, warpPCA, ...
                             vgrfACP, warpACP )
 

output = cell( 8, 1 );                   
                         
% 1: unrotated PCA scores
fullname = [name 'UA_PCA'];
output{1} = outputTable( fullname, refSet, withArms, perf, ...
                         vgrfPCA.unrotated );

% 2: varimax PCA scores
fullname = [name 'VA_PCA'];
output{2} = outputTable( fullname, refSet, vgrfPCA.varimax, perf );

% 3: unrotated ACP scores
fullname = [name 'UA_ACP'];
output{3} = outputTable( fullname, refSet, vgrfACP.unrotated, perf );

% 4: varimax ACP scores
fullname = [name 'VA_ACP'];
output{4} = outputTable( fullname, refSet, vgrfACP.varimax, perf );

if ~isempty( warpPCA )
    
    % 5: time warp PCA scores
    fullname = [name 'UW_PCA'];
    output{5} = outputTable( fullname, refSet, warpPCA.unrotated, perf );
    
    % 6: time warp PCA scores
    fullname = [name 'VW_PCA'];
    output{6} = outputTable( fullname, refSet, warpPCA.varimax, perf );

end

if ~isempty( warpACP )
    
    % 7: time warp ACP scores
    fullname = [name 'UW_ACP'];
    output{7} = outputTable( fullname, refSet, warpPCA.unrotated, perf );
    
    % 8: time warp PCA scores
    fullname = [name 'VW_ACP'];
    output{8} = outputTable( fullname, refSet, warpPCA.varimax, perf );
    
end        


end
