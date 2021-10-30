% ************************************************************************
% Function: statsDataDecomp
% Purpose:  Prepare the data set for statistical analysis in SAS
%           for the registration decompositional data.
%           This applies only to jumps without arm swing for simplicity.
%
% Parameters:
%       decomp: cell array of decomposition structrures
%
% Outputs:
%       data: table for statistical analysis
%
% ************************************************************************

function data = statsDataDecomp( decomp )

% assemble tables, one method at a time
data1 = compileResults( decomp(2,1,:) ); 
data1.Method(:) = "PAD";
nRows = size( data1, 1 );
data1.Registration(:) = string( dec2bin( 0:nRows-1, 4 ) );
data1.Registration(1) = "CCCC";

data2 = compileResults( decomp(2,2,:) ); 
data2.Method(:) = "LTN";
nRows = size( data2, 1 );
data2.Registration(:) = string( dec2bin( 0:nRows-1, 4 ) );
data2.Registration(1) = "CCCC";

% combine
data = [ data1; data2 ];

% define category orders
data.Method = categorical( data.Method, ...
            {'PAD', 'LTN'}, 'Ordinal', true );
data.Registration = categorical( data.Registration, ...
        {'0001', '0010', '0011', '0100', '0101', '0110', '0111', ...
         '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111', ...
         'CCCC' }, 'Ordinal', true);

end



