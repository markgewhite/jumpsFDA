% ************************************************************************
% Function: demarcatejump
% Purpose:  Identify the start and end times of a jump
%           by examining the vertical GRF
%           Start point is when the vGRF departs from the mean value
%           in the first second by 5 times the standard deviation
%           End point is when the vGRF reaches (< 10 N)
%
% Parameters:
%   vgrf: array with vertical GRF
%   bw: bodyweight
%
% Output:
%   t1: array index when jump starts
%   t2: array index at takeoff when jump ends
%
% ************************************************************************


function [ t1, t2 ] = demarcateJump( vgrf, bw, threshold, tMinStable )

% constants
mingrf = 10; % Newtons
plotVGRF = false;

if plotVGRF
    colours = lines(7);
end

% find take-off - the end demarcation pointv
t2 = find( abs(vgrf) < mingrf, 1 );

% normalise to bodyweights
vgrf = (vgrf-bw)/bw;

if plotVGRF
    % draw the VGRF curve
    plot( vgrf(1:t2),'LineWidth',1,'Color','b');
    hold on;
    % draw threshold lines
    plot( [1 t2], [threshold threshold], ...
                '--', 'LineWidth', 1.5, 'Color', 'k' );
    plot( [1 t2], [-threshold -threshold], ...
                '--', 'LineWidth', 1.5, 'Color', 'k' );
    % setup the axes
    xlim([t2-2500,t2]);
    ylim([-0.1 0.1]);
    xlabel('Time (ms)');
    ylabel('VGRF (BW) Centred');
    grid on;
end

% find the VGRF minimum prior to takeoff than is below 0 (BW)
[ ~, tMin ] = findpeaks( -vgrf(1:t2), ...
                        'MinPeakHeight', -2*threshold, ...
                        'MinPeakWidth', 20 );
tMin = tMin( end );

% work backwards to find where vGRF falls to < lower threshold
% ensure this is a stable low period 
t1 = tMin;
t3 = tMin;
while (t1-t3) < tMinStable
    % find the next point backwards where VGRF is within 2.5% of BW
    t1 = find( abs(vgrf(1:t3)) < threshold, 1, 'Last' );

    if plotVGRF
        % plot stability box
        box = patch( [t1 t1-tMinStable t1-tMinStable t1], ...
              [threshold threshold -threshold -threshold], ...
              colours(3,:), 'FaceAlpha', 0.2, 'EdgeAlpha', 0 );
    end
    
    % find there, find where the VGRF rises back above this threshold
    t3 = find( abs(vgrf(1:t1)) > threshold, 1, 'Last' );
    if isempty(t3)
        t3 = 1;
    end

end

if isempty(t1)
    t1 = 1;
else
    [~, t0] = min( abs(vgrf(t1-tMinStable+1:t1)) );
    t1 = max( t0 + t1 - tMinStable, 1 );
    %t1 = max( t1 - fix(0.5*tMinStable), 1 );
end



if plotVGRF

    % plot the initiation point on the VGRF curve
    plot( t1, vgrf(t1), 'o', 'MarkerFaceColor', 'k', ...
                             'MarkerEdgeColor', 'r', ...
                             'MarkerSize', 8, ...
                             'LineWidth', 1 );

    hold off;
    pause;

end
   
end
