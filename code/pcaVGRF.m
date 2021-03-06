% ************************************************************************
% Function: pcaVGRF
% Purpose:  Perform functional principal component analysis on VGRF data
%
% Parameters:
%       XFd: smoothed curves
%       XFdPar: smoothing parameters
%       warpFd: smoothed time warp (optional)
%       nRetained: number of components to retain
%       nRetainedWarp: number of components to retain from the warp
%
% Output:
%       pca: PCA structure
%
% ************************************************************************

function pca = pcaVGRF( XFd, XFdPar, warpFd, nRetained, nRetainedWarp )


pca.unrotated = pca_fd( XFd, ...
                        nRetained, ...
                        XFdPar );
                    
pca.varimax = varmx_pca( pca.unrotated );

pca.unrotated = compactPCA( pca.unrotated );
pca.varimax = compactPCA( pca.varimax );

if ~isempty( warpFd )
    
    pca.warp = pca_fd(  warpFd, ...
                        nRetainedWarp, ...
                        XFdPar );
                    
    pca.warpVarimax= varmx_pca( pca.warp );
    
    pca.warp = compactPCA( pca.warp );
    pca.warpVarimax = compactPCA( pca.warpVarimax );
    
    % re-scale warp harmonic scores from milliseconds to seconds
    pca.warp.harmscr = pca.warp.harmscr/1000;
    pca.warpVarimax.harmscr = pca.warpVarimax.harmscr/1000;
    
end

disp(['Total explained variance = ' ...
                    num2str( sum(pca.unrotated.varprop) ) ]);
                

end


function pca = compactPCA( pca )

% strip out fields that are unnecessary in order to save memory

pca = rmfield( pca, 'fdhatfd' );


end