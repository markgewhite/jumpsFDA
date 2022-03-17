% ************************************************************************
% Function: jumpPerfDescr
% Purpose:  Jump performance descriptive statistics
%
% Parameters:
%       results: cell array of results tables
%
% Outputs:
%       Group summary statistics
%
% ************************************************************************

function jumpPerfDescr( results )

% assemble tables for both jump types together 
T = compileResults( results ); 

% identify the subset of jumps without arms
woaIdx = T.Dataset(:,1)=='1';

% and jumps with arms
waIdx = ~woaIdx & T.Arms==1;

disp('Jumps without arms:');
disp([ 'JH = ' num2str( mean(T.JH_TOV(woaIdx)), '%.3f') ...
       ' +/- ' num2str( std(T.JH_TOV(woaIdx)), '%.3f') ]);
disp([ 'PP = ' num2str( mean(T.PP(woaIdx)), '%.3f') ...
       ' +/- ' num2str( std(T.PP(woaIdx)), '%.3f') ]);
   
disp('Jumps with arms:');
disp([ 'JH = ' num2str( mean(T.JH_TOV(waIdx)), '%.3f') ...
       ' +/- ' num2str( std(T.JH_TOV(waIdx)), '%.3f') ]);
disp([ 'PP = ' num2str( mean(T.PP(waIdx)), '%.3f') ...
       ' +/- ' num2str( std(T.PP(waIdx)), '%.3f') ]);

end



