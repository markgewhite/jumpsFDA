% ************************************************************************
% Script: validateJumpsFDA
% Purpose:  Check the jumpsFDA analysis
%
% ************************************************************************

function validateJumpsFDA( warpFd, showAll )
% check that the time warps are monotonic

if nargin<2
    showAll = false;
end

[nSets, nStd, nLMReg, nCTReg] = size( warpFd );

for i = 1:nSets
    for j = 1:nStd
        if isempty( warpFd{i,j,2,1} )
            continue
        end

        basisFd = getbasis( warpFd{i,j,2,1} );
        tRange = getbasisrange( basisFd );
        tSpanWarp = [ tRange(1) getbasispar( basisFd ) tRange(2) ];

        maxPts = 0;
        for k = 2:nLMReg
            for l = 1:nCTReg
                if isempty( warpFd{i,j,k,l} )
                    continue
                end
                warpDerivPts = eval_fd( tSpanWarp, warpFd{i,j,k,l}, 1 );
                nPts = size(warpDerivPts,2);
                if nPts > maxPts
                    maxPts = nPts;
                end
                if any( warpDerivPts<0, 'all' ) || showAll
                    disp( ['Dataset: Set = ' num2str(i) ...
                           '; Norm = ' num2str(j) ...
                           '; LMReg = ' num2str(k) ...
                           '; CTReg = ' num2str(l) ] );
                    disp( ['Number of curves = ' ...
                        num2str() ]);
                    minima = min(warpDerivPts);
                    disp( ['Number of non-monotonic curves = ' ...
                        num2str(sum(minima<0)) ]);
                    disp( ['Most negative = ' ...
                        num2str(min(minima))] );
                    disp( ' ' );
                end
                if ~showAll && nPts < maxPts
                    disp( ['Dataset: Set = ' num2str(i) ...
                           '; Norm = ' num2str(j) ...
                           '; LMReg = ' num2str(k) ...
                           '; CTReg = ' num2str(l) ] );
                    disp( ['Missing curves = ' ...
                        num2str(maxPts-nPts) ]);
                end


            end
        end
    end
end

end


