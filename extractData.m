% ************************************************************************
% Function: extractData
% Purpose:  Extract vertical jumps from the first data collection
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


function [ curveSet, IDSet, typeSet ] =  extractData( ...
                                    grf, bwall, nJumpsPerSubject, ...
                                    sDataID, sJumpID, jumpOrder, ...
                                    subjectExclusions, jumpExclusions )

% constants
totalID = 3; % combined data from both force platforms 
vAxisID = 3; % vertical axis

subjectRefID = 1; % identifier for subject ID index in ref
jumpRefID = 2; % identifier for jump ID in ref

nSubjects = length( sDataID );

% nSubjects = 2; % OVERRIDE %
nTotal = sum( nJumpsPerSubject );

% VGRF data (vertical jumps only)
vgrfData = cell( nTotal, 1 );
vgrfRef = zeros( nTotal, 2 );
withArms = false( nTotal, 1 );
initiation = zeros( nTotal, 1 );
takeoff = zeros( nTotal, 1 );
flightTime = zeros( nTotal, 1 );
typErr = zeros( nTotal, 1 );
k = 0;

for i = 1:nSubjects
    for j = 1:nJumpsPerSubject(i)
        
        jump = jumpOrder( sJumpID==sDataID(i), j );
        jumpID = sDataID(i)*100+j;
        
        if jump{1}(1) == 'V' ...
                && ~ismember( i, subjectExclusions ) ...
                && ~ismember( jumpID, jumpExclusions )
            
            k = k+1;
            vgrfData{ k } = ...
                    grf.raw{i,j,totalID}( 1:grf.takeoff(i,j), vAxisID ) ...
                               / bwall(i,j);

            vgrfRef( k, subjectRefID ) = i;
            vgrfRef( k, jumpRefID ) = j;
            withArms( k ) = (length(jump{1}) == 2);
            initiation( k ) = grf.initiation(i,j);
            typErr( k ) = std( vgrfData{k}(1:100) );
            takeoff( k ) = grf.takeoff(i,j);
            flightTime( k ) = jumpflighttime( ...
                                    grf.raw{i,j,totalID}(:,vAxisID), 1 );
        end
        
    end
end

vgrfData = vgrfData(1:k);
vgrfRef = vgrfRef(1:k,:);
withArms = withArms(1:k,:);
initiation = initiation(1:k);
takeoff = takeoff(1:k);
flightTime = flightTime(1:k);
typErr = typErr(1:k);

disp(['Typical VGRF Error = ' num2str(mean(typErr)) ' BW']);

% combine all datasets together
curveSet = { vgrfData, vgrfData(~withArms), vgrfData(withArms) };
IDSet = { vgrfRef, vgrfRef(~withArms,:), vgrfRef(withArms,:) };
typeSet = { withArms, false(length(IDSet{2}),1), true(length(IDSet{3}),1) };

end