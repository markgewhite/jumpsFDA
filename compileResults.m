% ************************************************************************
% Function: compileResults
% Purpose:  Compile a single table of results from a cell array
%
% Parameters:
%       results: results table from the functional data analysis
%       field: field identifying sub-table
%
% Outputs:
%       longResults: full assembled table
%
% ************************************************************************

function longT = compileResults( results, field )

isSubTable = nargin==2;

[ nSets, nStd, nProc ] = size( results );

longT = table([]);
for i = 1:nSets
    for j = 1:nStd
        for k = 1:nProc

            T = results{i,j,k};
            if isSubTable
                if isfield( T, field )
                    T = results{i,j,k}.(field);
                end
            end
            
            if istable( T )
                if isempty( longT )
                    longT = T;
                else
                    longT = [ longT; T ]; %#ok<AGROW>
                end
            end
                
        end
    end
end

end
  