% ************************************************************************
% Function: extractVGRFData
% Purpose:  Extract vertical jumps from the first data collection
%           Data from both force plates have already been added together
%           and the only the vertical component retained.
%
% Parameters:
%       grf: ground reaction force data structure
%       bwall: all bodyweights (by jump)
%       sDataID: assignment of subject ID
%       sJumpID: assignment of jumps (due to discards)
%       jumpOrder: array identifying the jump type
%       subjectExclusions: which subjects to exclude
%       jumpExclusions: which specific jumps to exclude
%
% Output:
%       curveSet: extracted VGRF data structure
%       IDSet: associated identifiers
%       typeSet: jump type
%
% ************************************************************************


function [ curveSet, IDSet, typeSet ] =  extractVGRFData( ...
                                    grf, bwall, nJumpsPerSubject, ...
                                    sDataID, sJumpID, jumpOrder, ...
                                    subjectExclusions, jumpExclusions )

detectionThreshold = 0.05; % proportion of bodyweight
tMinStable = 250; % ms - minimum period below detection threshold
zScoreLimit = 3; % outliers are beyond this limit

newDetectionMethod = true;
visualCheck = false;

nSubjects = length( sDataID );

nTotal = sum( nJumpsPerSubject );

% VGRF data (vertical jumps only)
vgrfData = cell( nTotal, 1 );
vgrfRef = zeros( nTotal, 2 );
withArms = false( nTotal, 1 );
typErr = zeros( nTotal, 1 );
fltTime = zeros( nTotal, 1 );
k = 0;
c = 0;
for i = 1:nSubjects
    
    for j = 1:nJumpsPerSubject(i)
        
        jump = jumpOrder( sJumpID==sDataID(i), j );
        jumpID = sDataID(i)*100+j;

        if jump{1}(1) == 'V' ...
                && ~ismember( i, subjectExclusions ) ...
                && ~ismember( jumpID, jumpExclusions )

            c = c+1;
            if newDetectionMethod
                % find jump initiation and jump take-off
                [tStart, tEnd, valid] = demarcateJump( grf.raw{i,j}, ...
                                                bwall(i,j), ...
                                                detectionThreshold, ...
                                                tMinStable );
    
                % check validity
                vgrfBW = grf.raw{i,j}/bwall(i,j);
                valid = valid && ...
                    all(abs( vgrfBW(1:tStart)-1 ) < detectionThreshold);

            else
                % original method
                tStart = grf.initiation(i,j);
                tEnd = grf.takeoff(i,j);
                vgrfBW = grf.raw{i,j}/bwall(i,j);
                valid = true;
            end


            if valid

                if visualCheck && i==42
                   plot( grf.raw{i,j}( tStart:tEnd )/bwall(i,j) );
                   disp( num2str(jumpID) );
                   pause;
                end

                k = k+1;
            
                % store the VGRF data in bodyweight units
                vgrfData{ k } = grf.raw{i,j}( tStart:tEnd )/bwall(i,j);
    
                vgrfRef( k, : ) = [ i j ];
                withArms( k ) = (length(jump{1}) == 2);
                typErr( k ) = std( vgrfData{k}(1:100) );

                fltTime( k ) = jumpFlightTime( grf.raw{i,j}, 1 );

            end
                
        end


        
    end
end

% trim the arrays
vgrfData = vgrfData(1:k);
vgrfRef = vgrfRef(1:k,:);
withArms = withArms(1:k,:);
typErr = typErr(1:k);
fltTime = fltTime(1:k);

% check lengths
len = cellfun( @length, vgrfData(1:k) );

disp(['Detection Threshold = ' num2str(detectionThreshold) ' BW']);
disp(['Min Stable Period   = ' num2str(tMinStable) ' ms']);

outliers = (abs(len-mean(len)))/std(len) > 3.5;
disp(['Outliers (Z-score > ' num2str(zScoreLimit) ') = ' ...
             num2str(sum(outliers)) ' {' num2str(len(outliers)') '}']);
disp('Excluded');

% exclude outliers 
vgrfData = vgrfData(~outliers);
vgrfRef = vgrfRef(~outliers,:);
withArms = withArms(~outliers,:);
typErr = typErr(~outliers);
fltTime = fltTime(~outliers);
len = len(~outliers);
k = k - sum(outliers);

% report
disp(['Valid Jumps = ' num2str(k)]);
disp(['Proportion with Arms = ' num2str(100*sum(withArms)/k, '%.0f') '%']);
disp(['Subjects = ' num2str( length(unique(vgrfRef(:,1))) ) ]);
disp(['Median Duration = ' num2str( median( len ) ) ]);
disp(['Mode Duration   = ' num2str( mode( len ) ) ]);
disp(['5th Percentile  = ' num2str( fix(prctile( len, 5 ) )) ]);
disp(['95th Percentile = ' num2str( ceil(prctile( len, 95 ) )) ]);
disp(['Min Duration    = ' num2str( min( len ) ) ]);
disp(['Max Duration    = ' num2str( max( len ) ) ]);



disp(['Typical VGRF Error = ' num2str(mean(typErr)) ' BW']);
disp(['Flight Time = ' num2str(mean(fltTime), '%.0f') ' +/- ' ...
                        num2str(std(fltTime), '%.0f') ' ms' ]);

% combine all datasets together
curveSet = { vgrfData(~withArms), vgrfData };
IDSet = { vgrfRef(~withArms,:), vgrfRef };
typeSet = { false(length(IDSet{1}),1), withArms };

end