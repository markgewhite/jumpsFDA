% ************************************************************************
% Function: statsData
% Purpose:  Prepare the data set for statistical analysis in SAS
%
% Parameters:
%       perf: cell array of model performance tables
%
% Outputs:
%       perf: table for statistical analysis
%
% ************************************************************************

function perf = statsData( glmModels )

% assemble tables for both jump types together 
perf1 = compileResults( glmModels(1,:,:), 'perf' ); 

% assemble tables for jumps without arm swing
perf2 = compileResults( glmModels(2,:,:), 'perf' );

% combine
perf = [ perf1; perf2 ];
perf.PredictorSet = extractBetween( perf.PredictorSet, 1, 3 );

% create long table so there is one line for train and test
nRows = size( perf, 1 );
perf = [ perf; perf ];

isTrain = 1:2*nRows <= nRows;

% assign classes to each of the four blocks
perf.Set( isTrain  ) = "Train";
perf.Set( ~isTrain  ) = "Test";
perf.Outcome( isTrain ) = perf.TrainRSq( isTrain );
perf.Outcome( ~isTrain ) = perf.TestRSq( ~isTrain );
perf.Outcome_NoWarp( isTrain ) = perf.TrainRSq_NoWarp( isTrain );
perf.Outcome_NoWarp( ~isTrain ) = perf.TestRSq_NoWarp( ~isTrain );


% add signal length standardisation method 
usedPAD = extract( perf.Dataset, 3 )=='1';
perf.Method( usedPAD ) = "PAD";
perf.Method( ~usedPAD ) = "LTN";

% extract the registration from the dataset label
perf.LMReg = extractBetween( perf.Dataset, 8, 11 );
perf.CTReg = extractAfter( perf.Dataset, 11 );

% make Outcome a percentage
perf.Outcome = perf.Outcome*100;
perf.Outcome_NoWarp = perf.Outcome_NoWarp*100;

% remove unwanted fields
perf.Dataset = [];
perf.TrainRSq = [];
perf.TrainRSq_NoWarp = [];
perf.TestRSq = [];
perf.TestRSq_NoWarp = [];

% define category orders
perf.Model = categorical( perf.Model, ...
            {'JHtov', 'JHwd', 'PP', 'Classify'}, 'Ordinal', true );
perf.PredictorSet = categorical( perf.PredictorSet, ...
            {'PCA', 'ACP'}, 'Ordinal', true );
perf.Set = categorical( perf.Set, ...
            {'Train', 'Test'}, 'Ordinal', true );
perf.Method = categorical( perf.Method, ...
            {'PAD', 'LTN'}, 'Ordinal', true );
perf.LMReg = categorical( perf.LMReg, ...
        {'0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111', ...
         '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111', ...
         }, 'Ordinal', true);
perf.CTReg = categorical( perf.CTReg, ...
            {'-', 'C'}, 'Ordinal', true );

end



