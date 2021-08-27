% ************************************************************************
% Function: validateSmoothing
% Purpose:  Validate the functional smoothing
%
% Parameters:
%       X: curve matrix
%       timespan: time series of points
%       setup: smoothing settings
%       perfCrit: jump performance criterion measures
%
% Output:
%       Graphs
%
% ************************************************************************

function validateSmoothing( X, timespan, setup, perfCrit )

N = size( X, 2 );

tStart = timespan(1);
tEnd = timespan(end);
basis = create_bspline_basis( [tStart,tEnd], setup.nBasis, setup.basisOrder);

% Find minimum GCV value of lambda
% search for the best value for lambda, the roughness penalty
logLambda   = -8:1:8;
gcvSave = zeros( length(logLambda), size(X,3) );
dfSave  = zeros( length(logLambda), 1 );
jh1Save  = zeros( length(logLambda), 1 );
jh2Save  = zeros( length(logLambda), 1 );
ppSave  = zeros( length(logLambda), 1 );

for i = 1:length(logLambda)
    
    % set smoothing parameters
    XFdPari = fdPar( basis, setup.penaltyOrder, 10^logLambda(i) );
    
    % perform smoothing
    [XFd, dfi, gcvi] = smooth_basis( timespan, X, XFdPari );
    
    % determine mean GCV and degrees of freedom
    gcvSave(i,:) = sqrt( sum( gcvi )/N ); 
    dfSave(i)  = dfi;
    
    % compute jump performances
    perfi = jumpperf_fd( XFd );

    % determine error with reference to lightest smoothing
    jh1Save(i) = sqrt( sum(( perfi.JHtov-perfCrit.JHtov ).^2)/N );
    jh2Save(i) = sqrt( sum(( perfi.JHwd-perfCrit.JHwd ).^2)/N );
    ppSave(i) = sqrt( sum(( perfi.PP-perfCrit.PP ).^2)/N );
    
end

%  plot the results for GCV and DF
figure;

plot( logLambda, log10(gcvSave), 'k-o' );
ylabel('\fontsize{13} log_{10}( GCV ))');
hold on;

yyaxis right;
plot( logLambda, log10(dfSave), 'r-o' );
ylabel('\fontsize{13} log_{10}( DF ))');

xlabel('\fontsize{13} log_{10}(\lambda)');

% plot the results for jump performance errors
figure;

plot( logLambda, log10(jh1Save), 'k-o' );
hold on;
jhTol1 = log10( setup.tolerance*mean( perfCrit.JHtov ) );
plot( [logLambda(1) logLambda(end)], [jhTol1 jhTol1], 'k--' );

plot( logLambda, log10(jh2Save), 'b-o' );
jhTol2 = log10( setup.tolerance*mean( perfCrit.JHwd ) );
plot( [logLambda(1) logLambda(end)], [jhTol2 jhTol2], 'b--' );


ylabel('\fontsize{13} log_{10}( RMSE( JH ) )');

yyaxis right;
plot( logLambda, log10(ppSave), 'r-o' );

ppTol = log10( setup.tolerance*mean( perfCrit.PP ) );
plot( [logLambda(1) logLambda(end)], [ppTol ppTol], 'r--' );

ylabel('\fontsize{13} log_{10}( RMSE( PP ) )');

xlabel('\fontsize{13} log_{10}(\lambda)');

end