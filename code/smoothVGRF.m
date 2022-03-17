% ************************************************************************
% Function: smoothVGRF
% Purpose:  Smooth the VGRF data
%
% Parameters:
%       X: raw data (of a standardised length)
%       t: time series of points
%       setup: smoothing settings
%       options: relevant options
%
% Output:
%       XFd: smoothed curves
%       XFdPar: smoothing parameters
%
% ************************************************************************

function [ XFd, XFdPar ] = smoothVGRF( X, t, setup, option )

if nargin < 4
    option.doCheckFit = false;
end

% set up a basis and fdPar object for smoothing the data
basis = create_bspline_basis( [ t(1), t(end) ], ...
                              setup.nBasis, ...
                              setup.basisOrder);
                          
XFdPar = fdPar( basis, setup.penaltyOrder, setup.lambda ); 

% create the FD object

[ XFd, df, gcv ] = smooth_basis( t, X, XFdPar);
XFd = putnames( XFd, setup.names);

disp(['Degrees of freedom = ',num2str(df)])
disp(['Generalized cross-validation = ',num2str(sum(gcv,2)')])

% check the goodness of fit
if option.doCheckFit
    %  plot each curve along with the data
    clf;
    plotfit_fd( X, t, XFd, option.doResidualCheck, option.doSortCheck );
end


end