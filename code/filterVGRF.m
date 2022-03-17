% ************************************************************************
% Function: filterVGRF
% Purpose:  Filter the VGRF raw data using a low-pass filter
%
% Parameters:
%       X: raw data (of a standardised length)
%       fs: sampling frequency
%       fs: cut-off frequency
%
% Output:
%       Xf: filtered data
%
% ************************************************************************


function Xf = filterVGRF( X, fs, fc, padding )

N = size( X, 2 );

Xpad = [ ones( padding, N ); X; zeros( padding, N ) ];

[ butterB, butterA ] = butter( 4, fc/(fs/2), 'low');

Xf = filter( butterB, butterA, Xpad );

Xf = Xf( padding+1:end-padding, : );

trialErrors = sqrt( sum( (Xf-X).^2, 1 )/size(X,1) );
meanError = sum( trialErrors ) / N;
disp(['Mean VGRF RMSE = ' num2str(meanError) ]);

end


