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
nProcustes = setup.nIterations;
N = size( getcoef( XFd ), 2 );

XFdReg = XFd;
if nargin < 5
    warpT = repmat( t, N, 1 );
else
    warpT = eval_fd( repmat( t, N, 1 ), warpFd0 );
end

% use a Procustes style loop
for i = 1:nProcustes   
    
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

            [ XFdReg, warpFd ] = register_fd( ...
                                                XMeanFd, XFdReg, wFdRegPar) ;                       
                                   
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
    warpFd = smooth_basis( t, warpT, wFdRegPar );
        
end


end