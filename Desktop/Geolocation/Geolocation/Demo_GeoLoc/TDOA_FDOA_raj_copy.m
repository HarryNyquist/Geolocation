
%% Figure 4, Hybrid Geometry (recreation of figure 13.1)
%x_source = [2;4]; % Transmitter/source
x_source = [50e3; 500e3];
v_sensor = [10, 10; 10, 10; 10,10; 10,10]';
x_sensor = [0, 0; 80e3, 0; 80e3, 120e3; 0, 120e3]';



fig4=figure;
plot(x_source(1),x_source(2),'k^','DisplayName','Source');
hold on;

% Position Markers
plot(x_sensor(1,:),x_sensor(2,:),'ko','DisplayName','Sensors');

% Draw velocity arrows
% drawArrow(x_sensor(1,1)+[0 v_sensor(1,1)]/4,x_sensor(2,1)+[0 v_sensor(2,1)]/4);
% drawArrow(x_sensor(1,2)+[0 v_sensor(1,2)]/4,x_sensor(2,2)+[0 v_sensor(2,2)]/4);
% drawArrow(x_sensor(1,3)+[0 v_sensor(1,3)]/4,x_sensor(2,3)+[0 v_sensor(2,3)]/4);
% drawArrow(x_sensor(1,4)+[0 v_sensor(1,4)]/4,x_sensor(2,4)+[0 v_sensor(2,4)]/4);
% % -- Direction of Arrival
% 
% % Compute Ranges
% r = utils.rng(x_sensor,x_source);
% 
% % Error Values
% epsang = 5*pi/180;
% 
% % Find AOA 
% lob = x_source - x_sensor;
% psi = atan2(lob(2,:),lob(1,:));
% 
% xaoa1 = x_sensor(:,1) + [0 cos(psi(1));0 sin(psi(1))]*5*r(1);
% xaoap1 = x_sensor(:,1) + [0 cos(psi(1)+epsang);0 sin(psi(1)+epsang)]*5*r(1);
% xaoam1 = x_sensor(:,1) + [0 cos(psi(1)-epsang);0 sin(psi(1)-epsang)]*5*r(1);
% lobFill1 = cat(2,xaoap1,fliplr(xaoam1),xaoap1(:,1));
% 
% xaoa2 = x_sensor(:,2) + [0 cos(psi(2));0 sin(psi(2))]*5*r(2);
% xaoap2 = x_sensor(:,2) + [0 cos(psi(2)+epsang);0 sin(psi(2)+epsang)]*5*r(2);
% xaoam2 = x_sensor(:,2) + [0 cos(psi(2)-epsang);0 sin(psi(2)-epsang)]*5*r(2);
% lobFill2 = cat(2,xaoap2,fliplr(xaoam2),xaoap2(:,1));
% 
% % LOBs
% plot(xaoa1(1,:),xaoa1(2,:),'k-','DisplayName','Line of Bearing');
% hdl2=plot(xaoa2(1,:),xaoa2(2,:),'k-');
% utils.excludeFromLegend(hdl2);

% -- Time Difference of Arrival
% Initialize Detector/Source Locations

% Isochrones
dr1 = rngDiff(x_source,x_sensor(:,1),x_sensor(:,2)); %between 1&2
dr2 = rngDiff(x_source,x_sensor(:,1),x_sensor(:,3));% between 1&3
dr3 = rngDiff(x_source,x_sensor(:,1),x_sensor(:,4));% between 1&4
%dr3 = rngDiff(x_source,x_sensor(:,1),x_sensor(:,3));% between 1&3
xiso1 = drawIsochrone(x_sensor(:,1),x_sensor(:,2),dr1,15000,800e3);%between 1&2
xiso2 = drawIsochrone(x_sensor(:,1),x_sensor(:,3),dr2,15000,800e3);% between 1&3
xiso3 = drawIsochrone(x_sensor(:,1),x_sensor(:,4),dr3,15000,800e3);% between 1&3
plot(xiso1(1,:),xiso1(2,:),'m--','LineWidth',1.3,'DisplayName','Line of Constant TDOA-12'); hold on;
hiso2=plot(xiso2(1,:),xiso2(2,:),'r:','LineWidth',1.3,'DisplayName','Line of Constant TDOA-13'); hold on;
hiso3=plot(xiso3(1,:),xiso3(2,:),'b:','LineWidth',1.3,'DisplayName','Line of Constant TDOA-14'); hold on;
%excludeFromLegend(hiso2);
% text(.5,-.5,'TDOA Solution');

% Isochrone Labels
% text(mean([x_sensor1(1),x_sensor2(1)]),mean([x_sensor1(2),x_sensor2(2)])-.2,'$TDOA_{1,2}$');
% text(mean([x_sensor2(1),x_sensor3(1)])+.3,mean([x_sensor2(2),x_sensor3(2)]),'$TDOA_{2,3}$');

% -- Frequency Difference of Arrival

% Draw isodoppler line
 text(x_sensor(1,1)-.2,x_sensor(2,1)-.2,'$S_1$','FontSize',10);
 text(x_sensor(1,2)-.2,x_sensor(2,2)-.2,'$S_2$','FontSize',10);
 text(x_sensor(1,3)-.2,x_sensor(2,3)-.2,'$S_3$','FontSize',10);
 text(x_sensor(1,4)-.2,x_sensor(2,4)-.2,'$S_3$','FontSize',10);
% Draw velocity arrows
 drawArrow(x_sensor(1,1)+[0 v_sensor(1,1)]/4,x_sensor(2,1)+[0 v_sensor(2,1)]/4);
 drawArrow(x_sensor(1,2)+[0 v_sensor(1,2)]/4,x_sensor(2,2)+[0 v_sensor(2,2)]/4);
 drawArrow(x_sensor(1,3)+[0 v_sensor(1,3)]/4,x_sensor(2,3)+[0 v_sensor(2,3)]/4);

% Draw isodoppler line S12
vdiff12 = utils.dopDiff(x_source,[0 0]',x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),3e8);
x_isodop12 = fdoa.drawIsodop(x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),vdiff12,1000,800e3);
plot(x_isodop12(1,:),x_isodop12(2,:),'-.','LineWidth',1.3,'DisplayName','Line of Constant FDOA-12');
%text(-1,2.7,'$S_{12}$ Solution','FontSize',10);

% Draw isodoppler line S23
vdiff23 = utils.dopDiff(x_source,[0 0]',x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),3e8);
x_isodop23 = fdoa.drawIsodop(x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),vdiff23,1000,800e3);
hh=plot(x_isodop23(1,:),x_isodop23(2,:),'g-.','LineWidth',1.3,'DisplayName','Line of Constant FDOA-23');

%plot(x,y,'--gs', 'LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b',...'MarkerFaceColor',[0.5,0.5,0.5])
%excludeFromLegend(hh);
%text(1.5,.85,'$S_{23}$ Solution','FontSize',10);
% vdiff1 = dopDiff(x_source,[0 0]',x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),1.62025e8);% between 1&2
% vdiff2 = dopDiff(x_source,[0 0]',x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),1.62025e8);% between 2&3
% vdiff3 = dopDiff(x_source,[0 0]',x_sensor(:,3),v_sensor(:,3),x_sensor(:,4),v_sensor(:,4),1.62025e8);% between 3&4
% x_isodop1 = drawIsodop(x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),vdiff1,1500,600e3);
% x_isodop2 = drawIsodop(x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),vdiff2,1500,600e3);
% x_isodop3 = drawIsodop(x_sensor(:,3),v_sensor(:,3),x_sensor(:,4),v_sensor(:,4),vdiff3,1500,600e3);
% plot(x_isodop1(1,:),x_isodop1(2,:),'-.','DisplayName','Line of Constant FDOA-12');
% %plot(x_isodop2(1,:),x_isodop2(2,:),'-.','DisplayName','Line of Constant FDOA-23');
% %plot(x_isodop3(1,:),x_isodop3(2,:),'-.','DisplayName','Line of Constant FDOA-34');
% % text(.5,4.1,'FDOA Solution');

xlim([-300e3 350e3]);
ylim([-600e3 750e3]);
% ylim([-5 5]);
% xlim([-2 4]);
legend('Location','SouthWest');

%setPlotStyle(gca,{'clean','widescreen','tight'});
%exportPlot(fig4,[prefix '4']);
%--------------------------------------------------------------------------

% %% Figure 3, FDOA Geometry (recreation of figure 12.1)
% 
% % Draw Geometry
% 
% 
% text(x_sensor(1,1)-.2,x_sensor(2,1)-.2,'$S_1$','FontSize',10);
% text(x_sensor(1,2)-.2,x_sensor(2,2)-.2,'$S_2$','FontSize',10);
% text(x_sensor(1,3)-.2,x_sensor(2,3)-.2,'$S_3$','FontSize',10);
% 
% % Draw velocity arrows
% drawArrow(x_sensor(1,1)+[0 v_sensor(1,1)]/4,x_sensor(2,1)+[0 v_sensor(2,1)]/4);
% drawArrow(x_sensor(1,2)+[0 v_sensor(1,2)]/4,x_sensor(2,2)+[0 v_sensor(2,2)]/4);
% drawArrow(x_sensor(1,3)+[0 v_sensor(1,3)]/4,x_sensor(2,3)+[0 v_sensor(2,3)]/4);
% 
% % Draw isodoppler line S12
% vdiff12 = utils.dopDiff(x_source,[0 0]',x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),3e8);
% x_isodop12 = fdoa.drawIsodop(x_sensor(:,1),v_sensor(:,1),x_sensor(:,2),v_sensor(:,2),vdiff12,1000,800e3);
% plot(x_isodop12(1,:),x_isodop12(2,:),'-.','DisplayName','Line of Constant FDOA');
% text(-1,2.7,'$S_{12}$ Solution','FontSize',10);
% 
% % Draw isodoppler line S23
% vdiff23 = utils.dopDiff(x_source,[0 0]',x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),3e8);
% x_isodop23 = fdoa.drawIsodop(x_sensor(:,2),v_sensor(:,2),x_sensor(:,3),v_sensor(:,3),vdiff23,1000,800e3);
% hh=plot(x_isodop23(1,:),x_isodop23(2,:),'-.');
% excludeFromLegend(hh);
% text(1.5,.85,'$S_{23}$ Solution','FontSize',10);



%--------------------------------------------------------------------------

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

function hdl = drawArrow(x,y,maxHeadSize)
% hdl = drawArrow(x,y,maxHeadSize)
%
% Draw an arrow on the current axes, with the provided coordinates.  Used
% as an alternative to the annotation function.  Uses the quiver function
% to draw the arrows.
%
% INPUTS:
%   x               2-element vector of x coordinate start/end positions
%   y               2-element vector of y coordinate start/end positions
%   maxHeadSize     Maximum head size, in points
%
% OUTPUTS:
%   hdl             Handle object for arrow drawn
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 3 || isempty(maxHeadSize)
    maxHeadSize=5;
end

hdl = quiver(x(1),     y(1),...
            x(2)-x(1),y(2)-y(1),...
            0,'MaxHeadSize',maxHeadSize );
hdl.LineStyle='-';
hdl.Color=[0 0 0];

%% Remove the arrow from legend entries
excludeFromLegend(hdl);
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
function ddop = dopDiff(x0,v0,x1,v1,x2,v2,f)
% ddop = dopDiff(x0,v0,x1,v1,x2,v2,f)
%
% Computes the difference in Doppler shift between reference and test 
% sensor pairs and a source.  The two sets of sensor positions (x1 and x2)
% must have the same size.  Corresponding pairs will be compared.
%
% INPUTS:
%   x0      Position vector for N sources (nDim x N), in m
%   v0      Velocity vector for N sources (nDim x N), in m/s
%   x1      Position vector for M reference sensors (nDim x M), in m
%   v1      Velocity vector for M reference sensors (nDim x M), in m/s
%   x2      Position vector for M test sensors (nDim x M), in m
%   v2      Velocity vector for M test sensors (nDim x M), in m/s
%   f       Carrier frequency, in Hertz
%
% OUTPUTS:
%   ddop    Different Doppler shift (N x M), in Hertz
%
% Nicholas O'Donoughue
% 1 July 2019

% Compute Doppler velocity from reference to each set of test positions
dop1 = dop(x0,v0,x1,v1,f); % N x M
dop2 = dop(x0,v0,x2,v2,f); % N x M

% Doppler difference
ddop = dop2 - dop1;
end
%-------------------------------------------------------------------------

function fd = dop(x1,v1,x2,v2,f)
% fd = dop(x1,v1,x2,v2,f)
%
% Given source and sensor at position x1 and x2 with velocity v1 and v2,
% compute the Doppler velocity shift
%
% INPUTS
%   x1      Position vector of N sources (nDim x N), in m
%   v1      Velocity vector of N sources (nDim x N), in m/s
%   x2      Position vector of M sensors (nDim x M), in m
%   v2      Velocity vector of M sensors (nDim x M), in m/s
%   f       Carrier frequency, in Hertz
%
% OUTPUTS
%   fd      Doppler shift for each source, sensor pair (N x M), in Hertz
%
% Nicholas O'Donoughue
% 1 July 2019

% Reshape inputs
[nDim,N] = size(x1);
[nDim2,M] = size(x2);

x1 = reshape(x1',N,1,nDim);
v1 = reshape(v1',[],1,nDim);
x2 = reshape(x2',1,M,nDim2);
v2 = reshape(v2',1,[],nDim2);

if nDim~=nDim2
    fprintf('Error: input dimensions do not match.');
    fd = [];
    return
end

% Unit vector from x1 to x2
u12 = (x2-x1)./sqrt(sum(abs(x2-x1).^2,3));
u21 = -u12;

% x1 velocity towards x2
vv1 = sum(v1.*u12,3);
vv2 = sum(v2.*u21,3);

% Sum of combined velocity
v = vv1 + vv2;

% Convert to Doppler
c = 3e8;
fd = f .* (1+ v./c);
end
%------------------------------------------------------------------------

function xy_iso = drawIsodop(x1,v1,x2,v2,vdiff,numPts,maxOrtho)
% xy_iso = drawIsodop(x1,v1,x2,v2,vdiff,numPts,maxOrtho)
%
% Finds the isochrone with the stated range rate difference from points x1
% and x2.  Generates an arc with 2*numPts-1 points, that spans up to
% maxOrtho distance from the intersection line of x1 and x2
%
% Inputs:
%   x1          Position of first sensor (Ndim x 1) [m]
%   v1          Velocity vector of first sensor (Ndim x 1) [m/s]
%   x2          Position of second sensor (Ndim x 1) [m]
%   v2          Velocity vector of second sensor (Ndim x 1) [m/s]
%   vdiff       Desired velocity difference [m/s]
%   numPts      Number of points to compute
%   maxOrtho    Maximum offset from line of sight between x1 and x2 [m]
%
% Outputs:
%   x_iso       2 x numPts iso doppler contour [m]
%
% Nicholas O'Donoughue
% 1 July 2019

% Set frequency to 3e8, so that c/f_0 is unity, and output of utils.dopDiff
% is velocity difference [m/s]
f_0 = 3e8;

%% Set up test points
xx_vec = maxOrtho(:).*repmat(linspace(-1,1,numPts),2,1);
x_vec = xx_vec(1,:);
y_vec = xx_vec(2,:);
[XX,YY] = meshgrid(xx_vec(1,:),xx_vec(2,:));
x_plot = [XX(:),YY(:)]'; % numPts^2 x 2

df_plot = dopDiff(x_plot,[0 0]',x1,v1,x2,v2,f_0);
    % numPts^2 x (N-1)
    
    
%% Compute contour
fig00=figure;
[cc,hh] = contour(x_vec,y_vec,reshape(df_plot,numPts,numPts),[vdiff vdiff]);
% Close the figure generated
close(fig00);

x_iso=cc(1,:);
y_iso=cc(2,:);

%% Filter points out of bounds
out_of_bounds = abs(x_iso) > maxOrtho | abs(y_iso) > maxOrtho;
x_iso = x_iso(~out_of_bounds);
y_iso = y_iso(~out_of_bounds);

xy_iso = cat(1, x_iso(:)', y_iso(:)');
end
%-----------------------------------------------------------------------



