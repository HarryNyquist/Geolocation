
%% Figure 2, TDOA Geometry (recreation of figure 11.1b)
% Initialize Detector/Source Locations
x_sensor1 = [0; 0];
x_sensor2 = [120e3; 0];
x_sensor3 = [0; 100e3];
x_sensor4 = [120e3; 110e3];

x_source = [0; 500e3];


% Compute Ranges

r1 = rng(x_sensor1,x_source);
r2 = rng(x_sensor2,x_source);
r3 = rng(x_sensor3,x_source);
r4  = rng(x_sensor4,x_source);

% Find Isochrones
xiso1 = drawIsochrone(x_sensor1,x_sensor2,r2-r1,15000,650e3);
xiso2 = drawIsochrone(x_sensor2,x_sensor3,r3-r2,15000,650e3);
xiso3 = drawIsochrone(x_sensor1,x_sensor4,r4-r1,15000,620e3);

% Draw Figure
fig2 = figure();hold on;

% Isochrones
plot(xiso1(1,:),xiso1(2,:),':','DisplayName','Isochrone');
hiso2=plot(xiso2(1,:),xiso2(2,:),':');
hiso3=plot(xiso3(1,:),xiso3(2,:),':');
excludeFromLegend(hiso2);

% Isochrone Labels
text(mean([x_sensor1(1),x_sensor2(1)]),mean([x_sensor1(2),x_sensor2(2)])-.2,'$TDOA_{1,2}$');
text(mean([x_sensor2(1),x_sensor3(1)])+.3,mean([x_sensor2(2),x_sensor3(2)]),'$TDOA_{2,3}$');
text(mean([x_sensor1(1),x_sensor4(1)])+.3,mean([x_sensor1(2),x_sensor4(2)]),'$TDOA_{1,4}$')
% Position Markers
%hiso2=scatter([x_sensor1(1),x_sensor2(1),x_sensor3(1)],[x_sensor1(2),x_sensor2(2),x_sensor3(2)],'k-','LineWidth',1);
%utils.excludeFromLegend(hiso2);
scatter([x_sensor1(1),x_sensor2(1),x_sensor3(1)],[x_sensor1(2),x_sensor2(2),x_sensor3(2)],'^','filled','DisplayName','Sensors');
% hdl=scatter(x_sensor2(1),x_sensor2(2),'^','filled');
% utils.excludeFromLegend(hdl);
% hdl=scatter(x_sensor3(1),x_sensor3(2),'^','filled');
% utils.excludeFromLegend(hdl);
scatter(x_source(1),x_source(2),'s','filled','DisplayName','Transmitter');

% Position Labels
text(x_sensor1(1)+.05,x_sensor1(2)-.1,'$S_1$');
text(x_sensor2(1)+.05,x_sensor2(2)-.1,'$S_2$');
text(x_sensor3(1)+.05,x_sensor3(2)-.1,'$S_3$');
text(x_sensor4(1)+.05,x_sensor4(2)-.1,'$S_4$');
% text(x_source(1)+.05,x_source(2)+.05,'$T_1$');

% Adjust Axes
xlim([-200e3 250e3]);
ylim([-500e3 550e3]);
legend('Location','SouthEast');

%setPlotStyle(gca,{'clean','equal','widescreen','tight'});
%exportPlot(fig2,[prefix '2']);

%-----------------------------------------------------------------------

function r = rng(x1,x2)
% r = rng(x1,x2)
%
% Computes the range between two N-dimensional position vectors, using
% the Euclidean (L2) norm.
%
% Inputs:
%
%   x1          NxM1 matrix of N-dimensional positions
%   x2          NxM2 matrix of N-dimensional positions
%
% Outputs:
%
%   r           M1xM2 matrix of pair-wise ranges
%
% Nicholas O'Donoughue
% 1 July 2019

% Find the input dimensions and test for compatibility
[N1,M1] = size(x1);
[N2,M2] = size(x2);
if N1~=N2
    fprintf('Error; first dimension of x1 and x2 must match.\n');
end

% Reshape the inputs
x1 = reshape(x1,[N1,M1]);
x2 = reshape(x2,[N1,1,M2]);

% Compute the Euclidean Norm
r = reshape(sqrt(sum(abs(bsxfun(@minus,x1,x2)).^2,1)),[M1,M2]);

end

%------------------------------------------------------------------------
function dr = rngDiff(x0,x1,x2)
% dr = rngDiff(x0,x1,x2)
%
% Computes the difference in range from the reference input (x0) to each
% of the input vectors x1 and x2.  Difference is taken as the range from 
% x0 to each column in x2 minus the range from x0 to each column in x1, 
% pair-wise.
%
% The first dimension of all three inputs must match.  The dimensions of
% x1 and x2 determine the dimensions of the output dr.
%
% Inputs:
%
%   x0      Nx1 reference position
%   x1      NxM1 vector of test positions
%   x2      NxM2 vector of test positions
%
% Outputs:
%
%   dr      M1xM2 matrix of range differences
%
% Nicholas O'Donoughue
% 1 July 2019

% Compute the range from the reference position to each of the set of
% test positions
r1 = utils.rng(x0,x1); % 1xM1
r2 = utils.rng(x0,x2); % 1xM2

% Take the difference, with appropriate dimension reshaping
dr = bsxfun(@minus,r2(:)',r1(:));
end
%-------------------------------------------------------------------------
function iso = drawIsochrone(x1,x2,rdiff,numPts,maxOrtho)
% iso = drawIsochrone(x1,x2,rdiff,numPts,maxOrtho)
%
% Finds the isochrone with the stated range difference from points x1
% and x2.  Generates an arc with 2*numPts-1 points, that spans up to
% maxOrtho distance from the intersection line of x1 and x2
%
% Inputs:
%   x1          Sensor position 1
%   x2          Sensor position 2
%   rdiff       Desired isochrone range difference
%   numPts      Number of points to draw
%   maxOrtho    Maximum offset (orthogonal to the line connecting x1 and
%               x2) at which points should be drawn
%
% Outputs:
%   iso         2 x numPts array of isochrone coordinates
%
% Nicholas O'Donoughue
% 1 July 2019

% Generate pointing vectors u and v in rotated coordinate space
%  u = unit vector from x1 to x2
%  v = unit vector orthogonal to u
rotMat = [0,1; -1,0];
R = utils.rng(x1,x2);
u = (x2-x1)/R;
v = rotMat*u;

xProj = [u v];

% Position of reference points in uv-space
x1uv = [0;0];
x2uv = [R;0];

% Initialize isochrone positions in uv-space
vv = linspace(0,maxOrtho,numPts);
uu = zeros(size(vv));
uu(1) = (R-rdiff)/2;
xuv = [uu; vv];

% Integrate over points, search for isochrone position
maxIter = 10000;

for i=1:numPts
    if i>1
        xuv(1,i) = xuv(1,i-1); % Initialize u position with previous value
    end
    
    offset = R;
    numIter = 1;
    
    while offset > R/10000 && numIter <= maxIter
        numIter = numIter + 1;
        
        rdiff0 = utils.rngDiff(xuv(:,i),x1uv,x2uv);
        offset = rdiff0-rdiff;
        
        xuv(1,i) = xuv(1,i)+offset;
    end
end

% Isochrone is symmetric about u axis
    % Flip the u axis, flip and negate v
xuv = [fliplr([xuv(1,2:end);-xuv(2,2:end)]),xuv];

% Convert to x/y space and re-center at origin
iso = xProj*xuv + x1;
    
end
%------------------------------------------------------------------------
function excludeFromLegend(graphObj)
% excludeFromLegend(graphObj)
%
% Modifies the supplied graphics object such that it will not appead in a
% legend; useful for plotting annotations that do not belong in the legend,
% such as dashed lines between a label and the appropriate line or
% marker.
%
% The handles are returned from the plotting command that was used to
% generate them, such as
%   h = plot(...)
% or
%   h = text(...)
%
% Alternatively, an axis object can be searched with the findall or findobj
% commands, such as:
%   h1 = findall(gcf);
%   h2 = findall(gcf,'Type','text');
%   h3 = findobj(gcf,'Type','text');
%
% If multiple graph objects are input the command is applied to each.
%
% Inputs:
%   graphObj        Graphics object handles (scalar or array)
%
% Nicholas O'Donoughue
% 1 July 2019

if numel(graphObj)>1
    arrayfun(@(x) utils.excludeFromLegend(x),graphObj);
    return;
end

% Graphics Objects have an Annotation property.
% Within the Annotation object, there is a LegendInformation Property.
% The LegendInformation object has an IconDisplayStyle property.
%
% Set that IconDisplayProperty to 'off' and a graphics object will be
% excluded from legends
graphObj.Annotation.LegendInformation.IconDisplayStyle = 'off';

end
%------------------------------------------------------------------------
function setPlotStyle(ax,styleName)
% setPlotStyle(ax,styleName)
%
% Preset script for various plot styles, to handle quick setting of
% common parameters, such as "clean" for drawings with no need to mark the
% axes.
%
% If styleName is a cell array of strings, each string is interpreted, in
% order.  This gives preference to the last element in the array, which may
% override settings called by prior elements.
%
% Where possible, the visibility of an axis element is adjusted, rather
% than its presence; to allow for toggling by later commands.
%
% Inputs:
%
%   ax          Axis handle to modify
%   styleName   String or cell array of strings.  Valid options include:
%                   'clean'     strips axes
%                   'notick'    removes numbers and ticks, but leaves axes
%                   'equal'     sets scale equal for x and y axes
%                   'widescreen'sets the plot dimensions to a 16:9 format
%                   'tight'     minimizes white space outside axes
%
% Nicholas O'Donoughue
% 1 July 2019


% Iterate over comma-separated lists; in order.
if iscell(styleName)
    cellfun(@(x) utils.setPlotStyle(ax,x),styleName);
    return;
end

switch styleName
    case 'clean'
        % Remove the tick marks and grid
%        grid(ax,'off');
%        set(ax,'xtick',[]);
%        set(ax,'ytick',[]);
%        set(ax,'xticklabel',[]);
%        set(ax,'yticklabel',[]);
        ax.Visible ='off';
        
        if ~isempty(ax.Legend)
            ax.Legend.Box = 'off';
        end
    
    case 'notick'
        set(ax,'xticklabel',[]);
        set(ax,'yticklabel',[]);
        
    case 'box only'
        set(ax,'xticklabel',[]);
        set(ax,'yticklabel',[]);
        grid(ax,'off');
        
    case 'equal'
        axis(ax,'equal');
        
        % Set figure size to square
        %posVec = get(ax.Parent,'Position');
        %posVec(4)=posVec(3);
        %set(ax.Parent,'Position',posVec);
        
        % Set axes position to square
        %posVec = get(ax,'OuterPosition');
        %posVec(4) = posVec(3);
        %set(ax,'OuterPosition',posVec);
    case 'widescreen'
        % 16:9 format
        
        % Set figure size
        posVec = get(ax.Parent,'Position');
        posVec(4) = 9*posVec(3)/16;
        set(ax.Parent,'Position',posVec);
        
        % Set axis position
        posVec = get(ax,'OuterPosition');
        posVec(4) = 9*posVec(3)/16;
        set(ax,'OuterPosition',posVec);
    case 'tight'
        outerpos = ax.OuterPosition;
        ti = ax.TightInset; 
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
end
end


%-------------------------------------------------------------------------
function exportPlot(fig,fnm)
% exportPlot(fig,fnm)
%
% Exports a plot to file using the desired resolution and format for
% inclusion in a textbook.
%
% Currently exports both .Fig and .EPS files with 1200 dpi resolution.
%
% Inputs:
%   fig         Figure handle to export
%   fnm         Filename, including directory
%
% Nicholas O'Donoughue
% 1 July 2019

fig.InvertHardcopy = 'off';

saveas(fig,[fnm '.fig']);
print(fig,fnm,'-depsc','-r1200','-painters');
print(fig,fnm,'-dpng','-r1200');

end

%------------------------------------------------------------------------

