% ************************************************************************
% Function: statsDataCoeff
% Purpose:  Prepare the data set for statistical analysis in SAS
%           for the analysis by component
%
% Parameters:
%       perf: cell array of model performance tables
%
% Outputs:
%       perf: table for statistical analysis
%
% ************************************************************************

function data = statsDataCoeff( glmModels )

% assemble tables for both jump types together 
data1 = compileResults( glmModels(1,:,:), 'coeffRSq' ); 
% retain only classification results
retain1 = strcmp( data1.Model, 'jumpType' );
data1.Model( retain1 ) = "Classify";

% assemble tables for jumps without arm swing
perf2 = compileResults( glmModels(2,:,:), 'coeffRSq' );
% retain only regression results
retain2 = strcmp( perf2.Model, 'JHtov' ) ...
           | strcmp( perf2.Model, 'JHwd' ) ...
           | strcmp( perf2.Model, 'PP' );

% combine
data = [ data1( retain1, : ); perf2( retain2, : ) ];

% retain only PCAU (uncorrelated)
retain3 = strcmp( data.PredictorSet, 'PCAU' );
data = data( retain3, : );
data.PredictorSet = extractBetween( data.PredictorSet, 1, 3 );

% add signal length standardisation method 
usedPAD = extract( data.Dataset, 3 )=='1';
data.Method( usedPAD ) = "PAD";
data.Method( ~usedPAD ) = "LTN";

% extract the registration from the dataset label
data.Registration = extractBetween( data.Dataset, 8, 11 );


% re-badge continuous registration
isCTReg = strcmp( data.Registration, "0000" );
data.Registration( isCTReg ) = "CCCC";

% re-badge no registration
isNoReg = strcmp( data.Registration, "----" );
data.Registration( isNoReg ) = "0000";

% make Outcome a percentage
data(:,4:26) = array2table( table2array(data(:,4:26))*100 );

% define category orders
data.Model = categorical( data.Model, ...
            {'JHtov', 'JHwd', 'PP', 'Classify'}, 'Ordinal', true );
data.PredictorSet = categorical( data.PredictorSet, ...
            {'PCA'}, 'Ordinal', true );
data.Method = categorical( data.Method, ...
            {'PAD', 'LTN'}, 'Ordinal', true );
data.Registration = categorical( data.Registration, ...
        {'0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111', ...
         '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111', ...
         'CCCC' }, 'Ordinal', true);

end



