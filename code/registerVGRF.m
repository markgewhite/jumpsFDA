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
%       finalXFdReg: registered smoothed curves
%       finalWarpFd: time warping curves
%       i: number of iterations performed
%       finalValidity: array indicating which curves are valid
%
% ************************************************************************


function [ finalXFdReg, finalWarpFd, i, finalValidity ] = ...
                        registerVGRF( tSpan, XFd, type, setup, warpFd0 )

% initialise
N = size( getcoef( XFd ), 2 );

XFdReg = XFd;

% initialize the total warp time series to high resolution
if nargin < 5 || isempty( warpFd0 )
    totWarpT = repmat( tSpan, N, 1 )';
else
    totWarpT = eval_fd( tSpan, warpFd0 );
    totWarpT = max( min( totWarpT, tSpan(end) ), tSpan(1) );
end

% use a Procustes style loop
hasConverged = false;
prevC = 1;
isValid = true( 1, N );
i = 0;
while ~hasConverged && all(isValid) && i < setup.maxIterations 
    
    switch type
        
        case 'Landmark'
        
            % *** landmark registration ***
            lm = findGRFlandmarks( tSpan, XFdReg, setup.lm );
            disp(['Landmark means  = ' num2str( lm.mean )]);
            disp(['Landmark SDs    = ' num2str( std( lm.case ) )]);

            % sort landmarks into order (if necessary)
            [ lm.mean, lm.order ] = sort( lm.mean, 'ascend' );
            lm.case = lm.case( :, lm.order );

            wBasis = create_bspline_basis( [tSpan(1),tSpan(end)], ...
                                           setup.nBasis, ...
                                           setup.basisOrder, ...
                                           [ tSpan(1) lm.mean tSpan(end) ] );

            wFdReg = fd( zeros( setup.nBasis, 1), wBasis );
            wFdRegPar = fdPar( wFdReg, setup.basisOrder-2, setup.wLambda );

            [ XFdReg, warpFd ] = landmarkreg( ...
                                        XFdReg, lm.case, lm.mean, ...
                                        wFdRegPar, true, setup.XLambda );
            
                                
        case 'Continuous'
                                    
            % *** continuous registration  ***

            wBasis = create_bspline_basis( [tSpan(1),tSpan(end)], ...
                                           setup.nBasis, ...
                                           setup.basisOrder );

            wFdReg = fd( zeros( setup.nBasis, 1), wBasis );
            wFdRegPar = fdPar( wFdReg, setup.basisOrder-2, setup.wLambda );

            XMeanFd = mean( XFdReg );

            [ XFdReg, warpFd ] = register_fd( XMeanFd, XFdReg, wFdRegPar );                       
                                   
        otherwise
            error([ 'Unrecognised registration type: ' type ]);
            
    end
  
    
    % update time warp to maintain link with original
    totWarpT = updateTotalWarp( tSpan, totWarpT, warpFd );

    % convert the total warp time series into smooth functions
    % checking for monotonicity
    [ warpFd, isMonotonic ] = resmoothWarp( tSpan, totWarpT, setup.rewarp );

    % check landmarks as indication of validity even if not landmark reg
    hasValidLandmarks = validateLandmarks( tSpan, XFdReg );
    
    % check for convergence
    decomp = regDecomp( XFd, XFdReg, warpFd );
    hasConverged = (abs(prevC-decomp.c) < setup.convCriterion);
    prevC = decomp.c;

    % determine overall validity
    isValid = isMonotonic & hasValidLandmarks;
    
    i = i + 1;
    if ~(setup.allMustBeValid && ~all( isValid ))
        % accept registration
        finalWarpFd = warpFd;
        finalXFdReg = XFdReg;
        finalValidity = isValid;
    end
    if ~all(isValid)
        disp('Warp functions not all valid.');
        disp( ['Faulty registrations (non-monotonic) = ' ...
                    num2str(sum(~isMonotonic)) ] );
        disp( ['Faulty registrations (landmarks) = ' ...
                    num2str(sum(~hasValidLandmarks)) ] );
    end

end

if strcmp( type, 'Landmark' )
    % locate landmarks
    lm = findGRFlandmarks( tSpan, ...
                           selectFd( finalXFdReg, isMonotonic ), ...
                           setup.lm );
    disp(['Final landmark means = ' num2str( lm.mean )]);
    disp(['Final landmark SDs   = ' num2str( std( lm.case ) )]);
end



end


function totWarpT = updateTotalWarp( t, totWarpT, newWarpFd )

    N = size( totWarpT, 2 );
    for i = 1:N
        % warp the previous warp
        totWarpT(:,i) = eval_fd( totWarpT(:,i), newWarpFd(i) );
        % separate the points evenly
        totWarpT(:,i) = interp1( t, totWarpT(:,i), t, 'spline', 'extrap' );
        % impose limits in case of over/underflow
        totWarpT(:,i) = max( min( totWarpT(:,i), t(end) ), t(1) );
    end

end


function [ warpFd, isMonotonic ] = resmoothWarp( t, warpT, setup )

    % re-smooth the warping curve using a more extensive basis
    % to ensure the closest adherence to the interpolated points
    % using the minimal number of basis functions required 

    b = 4 + setup.basisOrder;
    allMonotonic = false;
    while ~allMonotonic && b < setup.nBasisMax

        % double number of basis functions
        b = (b-setup.basisOrder)*2 + setup.basisOrder;

        % smooth time series with the trial basis
        wBasis = create_bspline_basis( [t(1),t(end)], b, ...
                                       setup.basisOrder ); 
        wFdPar = fdPar( wBasis, ...
                        setup.basisOrder-2, ...
                        setup.wLambda );
        warpFd = smooth_basis( t, warpT, wFdPar );
    
        % compute the first derivative
        warpDT = eval_fd( t, warpFd, 1 );
        % check monotonicty by curve
        isMonotonic = all( warpDT > 0 );
        allMonotonic = all( isMonotonic );
        
    end
    disp(['Resmooth warp: number of basis functions = ' num2str(b)]);

end


function isValid = validateLandmarks( t, XFdReg )

    % find all possible landmarks
    setup.grfmin = true;
    setup.pwrmin = true;
    setup.pwrcross = true;
    setup.pwrmax = true;
    lm = findGRFlandmarks( t, XFdReg, setup );

    [ nCurves, nLMs ]= size( lm.case );
    isValid = true( nCurves, 1 );

    % check for any extreme outliers
    for i = 1:nLMs
        % compute mean and SD with outliers excluded
        limits = prctile( lm.case(:,i), [2.5 97.5] );
        outliers = (lm.case(:,i)<limits(1) | lm.case(:,i)>limits(2));
        sd = std( lm.case(~outliers, i) );
        avg = mean( lm.case(~outliers, i) );
        % check if within max Z score 
        % (very large because SD computed without outliers)
        isValid = isValid & (abs(lm.case(:,i)-avg)/sd < 50);
    end

    % check if any landmarks are out of sequence
    for i = 1:nCurves
        [~, order] = sort( lm.case(i,:) );
        isValid(i) = isValid(i) && all(diff(order)==1);
    end

    isValid = isValid';

end