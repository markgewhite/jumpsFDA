% ************************************************************************
% Function: acpVGRF
% Purpose:  Perform analysis of characteristing phases on VGRF data
%
% Parameters:
%       tSpan: time span
%       XFd: smoothed curves
%       warpFd: smoothed time warp (optional)
%       XPCA: pca structure
%
% Output:
%       ss: ACP similarity scores
%
% ************************************************************************

function ss = acpVGRF( tSpan, XFd, warpFd, XPCA )


ss.unrotated = acp( tSpan, XFd, XPCA.unrotated );

ss.varimax = acp( tSpan, XFd, XPCA.varimax );

if ~isempty( warpFd )
    
   ss.warp = acp( tSpan, XFd, XPCA.warp );

   ss.warpVarimax = acp( tSpan, XFd, XPCA.warp );
    
end
                

end
