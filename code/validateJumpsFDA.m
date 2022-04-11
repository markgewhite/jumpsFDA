% ************************************************************************
% Script: validateJumpsFDA
% Purpose:  Check the jumpsFDA analysis
%
% ************************************************************************

% check that the time warps are monotonic

[nSets, nStd, nLMReg, nCTReg] = size( warpFd );

for i = 1:nSets
    for j = 1:nStd
        basisFd = getbasis( warpFd{i,j,2,1} );
        tRange = getbasisrange( basisFd );
        tSpanWarp = [ tRange(1) getbasispar( basisFd ) tRange(2) ];

        for k = 2:nLMReg
            for l = 1:nCTReg
                warpDerivPts = eval_fd( tSpanWarp, warpFd{i,j,k,l}, 1 );
                if any( warpDerivPts<0, 'all' )
                    disp( ['Dataset: Set = ' num2str(i) ...
                           '; Norm = ' num2str(j) ...
                           '; LMReg = ' num2str(k) ...
                           '; CTReg = ' num2str(l) ] );
                    minima = min(warpDerivPts);
                    disp( ['Number of curves = ' ...
                        num2str(sum(minima<0)) ]);
                    disp( ['Most negative = ' ...
                        num2str(min(minima))] );
                    disp( ' ' );
                end
            end
        end
    end
end


