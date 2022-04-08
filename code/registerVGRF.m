% ************************************************************************
% Function: registerVGRF
% Purpose:  Perform curve registration on VGRF data
%
% Parameters:
%       t: time span
%       XFd: smoothed curves
%       setup: smoothing settings
%       warpFd0: prior warp function (optional)
%
% Output:
%       XFdReg: registered smoothed curves
%
% ************************************************************************


function [ XFdReg, warpFd ] = registerVGRF( t, XFd, type, setup, warpFd0 )

% initialise
monotonic = true;
N = size( getcoef( XFd ), 2 );

XFdReg = XFd;
if nargin < 5 || isempty( warpFd0 )
    warpT = repmat( t, N, 1 )';
else
    warpT = eval_fd( t, warpFd0 );
    warpT = max( min( warpT, t(end) ), t(1) );
end

% use a Procustes style loop
hasConverged = false;
prevRSq = 0;
while ~hasConverged   
    
    switch type
        
        case 'Landmark'
        
            % *** landmark registration ***
            
            % locate landmarks
            lm = findGRFlandmarks( t, XFdReg, setup.lm );

            % sort into order (in necessary)
            [ lm.mean, lm.order ] = sort( lm.mean, 'ascend' );
            lm.case = lm.case( :, lm.order );

            disp(['Landmark means = ' num2str( lm.mean )]);
            disp(['Landmark SDs   = ' num2str( std( lm.case ) )]);

            wBasis = create_bspline_basis( [t(1),t(end)], ...
                                           setup.nBasis, ...
                                           setup.basisOrder, ...
                                           [ t(1) lm.mean t(end) ] );

            wFdReg = fd( zeros( setup.nBasis, 1), wBasis );
            wFdRegPar = fdPar( wFdReg, 1, setup.wLambda );

            [ XFdReg, warpFd ] = landmarkreg( ...
                                        XFdReg, lm.case, lm.mean, ...
                                        wFdRegPar, monotonic, setup.XLambda );
            
                                
        case 'Continuous'
                                    
            % *** continuous registration  ***

            wBasis = create_bspline_basis( [t(1),t(end)], ...
                                           setup.nBasis, ...
                                           setup.basisOrder );

            wFdReg = fd( zeros( setup.nBasis, 1), wBasis );
            wFdRegPar = fdPar( wFdReg, 1, setup.wLambda );

            XMeanFd = mean( XFdReg );

            [ XFdReg, warpFd ] = register_fd( XMeanFd, XFdReg, wFdRegPar );                       
                                   
        otherwise
            error([ 'Unrecognised registration type: ' type ]);
            
    end
  
    
    % update time warp to maintain link with original
    for j = 1:N
        % warp the previous warp
        warpT(:,j) = eval_fd( warpT(:,j), warpFd(j) );
        % separate the points evenly
        warpT(:,j) = interp1( t, warpT(:,j), t, 'spline', 'extrap' );
        % impose limits in case of over/underflow
        warpT(:,j) = max( min( warpT(:,j), t(end) ), t(1) );
    end

    % re-smooth the warping curve using a more extensive basis
    % without landmarks
    wBasis = create_bspline_basis( [t(1),t(end)], ...
                                               4*setup.nBasis, ...
                                               setup.basisOrder ); 
    wFdPar = fdPar( wBasis, 1, setup.wLambda );
    warpFd = smooth_basis( t, warpT, wFdPar );
    
    % check on progress
    decomp = regDecomp( XFd, XFdReg, warpFd );

    hasConverged = (abs(decomp.rSq-prevRSq) < setup.convCriterion);

    prevRSq = decomp.rSq;

end

if strcmp( type, 'Landmark' )
    lm = findGRFlandmarks( t, XFdReg, setup.lm );
    disp(['Landmark means  = ' num2str( lm.mean )]);
    disp(['Landmark SDs    = ' num2str( std( lm.case ) )]);
end

end