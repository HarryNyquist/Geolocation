function [fig_geo,fig_err] = ex10_1()
% [fig_geo,fig_err] = ex10_1()
%
% Executes Example 10.1 and generates two figures
%
% INPUTS
%   none
%
% OUTPUTS
%   fig_geo     figure handle for geographic layout
%   fig_err     figure handle for error as a function of iteration
%
% Nicholas O'Donoughue
% 1 July 2019ans

% Define sensor positions
x_sensor = 30e3*[-1 0;0 0;1 0]';
N=size(x_sensor,2);

% Define source position
x_source = [15,100]'*1e3;

% Grab a noisy measurement
psi_act = triang.measurement(x_sensor,x_source);
C_psi = (1*pi/180)^2*eye(N);
psi = psi_act + sqrt(C_psi)*randn(N,1);
% r = 0:(1.2*max(utils.rng(x_source,x_sensor))); -- unused

% Compute Ranges
% r0 = utils.rng(x_sensor(:,1),x_source); -- unused
r1 = utils.rng(x_sensor(:,2),x_source);
r2 = utils.rng(x_sensor(:,3),x_source);
%r3 = utils.rng(x_sensor(:,4),x_source);
% Error Values
epsang = 2*sqrt(diag(C_psi));

% Find AOA 
% lob0 = x_source - x_sensor(:,1); -- unused
% xaoa0 = x_sensor(:,1) + [0 cos(psi(1));0 sin(psi(1))]*5*r1; -- unused
xaoap0 = x_sensor(:,1) + [0 cos(psi_act(1)+epsang(1));0 sin(psi_act(1)+epsang(1))]*5*r1;
xaoam0 = x_sensor(:,1) + [0 cos(psi_act(1)-epsang(1));0 sin(psi_act(1)-epsang(1))]*5*r1;
lobFill0 = cat(2,xaoap0,fliplr(xaoam0),xaoap0(:,1));

% lob1 = x_source - x_sensor(:,2); -- unused
% xaoa1 = x_sensor(:,2) + [0 cos(psi(2));0 sin(psi(2))]*5*r1; -- unused
xaoap1 = x_sensor(:,2) + [0 cos(psi_act(2)+epsang(2));0 sin(psi_act(2)+epsang(2))]*5*r1;
xaoam1 = x_sensor(:,2) + [0 cos(psi_act(2)-epsang(2));0 sin(psi_act(2)-epsang(2))]*5*r1;
lobFill1 = cat(2,xaoap1,fliplr(xaoam1),xaoap1(:,1));

% lob2 = x_source - x_sensor(:,3); -- unused
% xaoa2 = x_sensor(:,3) + [0 cos(psi(3));0 sin(psi(3))]*5*r2; -- unused
xaoap2 = x_sensor(:,3) + [0 cos(psi_act(3)+epsang(3));0 sin(psi_act(3)+epsang(3))]*5*r2;
xaoam2 = x_sensor(:,3) + [0 cos(psi_act(3)-epsang(3));0 sin(psi_act(3)-epsang(3))]*5*r2;
lobFill2 = cat(2,xaoap2,fliplr(xaoam2),xaoap2(:,1));

% % lob2 = x_source - x_sensor(:,4); -- unused
% % xaoa2 = x_sensor(:,3) + [0 cos(psi(3));0 sin(psi(3))]*5*r2; -- unused
% xaoap3 = x_sensor(:,4) + [0 cos(psi_act(4)+epsang(4));0 sin(psi_act(4)+epsang(4))]*5*r3;
% xaoam3 = x_sensor(:,4) + [0 cos(psi_act(4)-epsang(4));0 sin(psi_act(4)-epsang(4))]*5*r3;
% lobFill3 = cat(2,xaoap3,fliplr(xaoam3),xaoap3(:,1));

% Geometric Solutions
% x_centroid = triang.centroid(x_sensor,psi);
% x_incenter = triang.angle_bisector(x_sensor,psi);

% Iterative Methods
% epsilon = .01; % Stopping condition -- unused
nMC = 50;
numIters=200;
% alpha=.3;  % Gradient Descent Line Search parameter -- unused
% beta=.6;   % Gradient Descent Line Search parameter -- unused
% C_psi = eye(N);
psi = cat(2,psi,psi_act + sqrt(C_psi)*randn(N,nMC-1)); % preserve prior estimate
    % so that geometric and iterative solutions use same inputs for plot

x_ml = zeros(2,nMC);          
x_bf = zeros(2,nMC);
x_centroid = zeros(2,nMC);
x_incenter = zeros(2,nMC);
x_ls_full = zeros(2,numIters,nMC);
x_grad_full = zeros(2,numIters,nMC);
x_ls = zeros(2,nMC);
x_gd = zeros(2,nMC);
fprintf('Conducting MC trial for triangulation error...\n');
x_initial = [50;30e3];
for idx = 1:nMC
    if mod(idx,floor(nMC/80))==0
        fprintf('.');
    end
    %x_ml(:,idx) = triang.mlSoln(x_sensor,psi(:,idx),C_psi,x_initial,[50;50],.1);
    %x_bf(:,idx) = triang.bfSoln(x_sensor,psi(:,idx),C_psi,x_initial,[50;50],.1);
    %x_centroid(:,idx) = triang.centroid(x_sensor,psi(:,idx));
    %x_incenter(:,idx) = triang.angle_bisector(x_sensor,psi(:,idx));
    [x_ls(:,idx),x_ls_full(:,:,idx)] = triang.lsSoln(x_sensor,psi(:,idx),C_psi,x_initial,[],numIters,true,[]);
    [x_gd(:,idx),x_grad_full(:,:,idx)] = triang.gdSoln(x_sensor,psi(:,idx),C_psi,x_initial,[],[],[],numIters,true,[]);
end
%err_ml = x_source - x_ml;
%err_bf = x_source - x_bf;
%err_cnt = x_source - x_centroid;
%err_inc = x_source - x_incenter;
%err_ls = x_source - x_ls_full;
%err_grad = x_source - x_grad_full;
err_ls = x_source - x_ls;
err_grad = x_source - x_gd;
fprintf('done.\n');

x_gd_avg = mean(x_gd, 2);
x_ls_avg = mean(x_ls, 2);
    
theta = rad2deg(atan2(x_source(2)-x_sensor(2,2),x_source(1)-x_sensor(1,2)));
%theta_ml = rad2deg(atan2(x_ml_avg(2),x_ml_avg(1)));
theta_gd = rad2deg(atan2(x_gd_avg(2)-x_sensor(2,2),x_gd_avg(1)-x_sensor(1,2)));
theta_ls = rad2deg(atan2(x_ls_avg(2)-x_sensor(2,2),x_ls_avg(1)-x_sensor(1,2)));
% 
r_actual = norm(x_source - x_sensor(:,2));
%r_ml = norm(x_ml_avg - x_tdoa(:,1));
r_gd = norm(x_gd_avg - x_sensor(:,2));
r_ls = norm(x_ls_avg - x_sensor(:,2));

% Display Results
fprintf('\nActual Source Location is: (%.2f, %.2f)\n', x_source);
%fprintf('\nAverage ML solution: (%.2f, %.2f)', x_ml_avg);
fprintf('\nAverage GD solution: (%.2f, %.2f)', x_gd_avg);
fprintf('\nAverage LS solution: (%.2f, %.2f)\n', x_ls_avg);

fprintf('\n Actal r and theta is %.2f,%.2f ', r_actual, theta)
%fprintf('\n r and theta MLE pos estimate is: %.2f,%.2f', r_ml, theta_ml);
fprintf('\n r and theta from GD pos estimate is: %.2f,%.2f', r_gd,theta_gd);
fprintf('\n r and theta from LS pos estimate is: %.2f,%.2f', r_ls, theta_ls);

% Compute Deviations
x_ls_dev = [];
x_ls_dev_x = [];
x_ls_dev_y = [] ;

x_ls_dev = x_source - x_ls;
% Compute absolute deviations (magnitude only, no sign)
abs_dev_x = abs(x_ls_dev(1,:));
abs_dev_y = abs(x_ls_dev(2,:));

% Plot bar charts
figure;
% X deviations
subplot(1,2,1);
%bar(abs_dev_x, 'b');
histogram(x_ls_dev(1,:), 'BinWidth', 500,'FaceColor', 'b' )
title('Deviations in X from source');
xlabel('deviation in meters');
ylabel('Num of deviations');
legend('Location','NorthWest');
grid on;

% Y deviations
subplot(1,2,2);
%bar(abs_dev_y, 'r');
histogram(x_ls_dev(2,:), 'BinWidth', 1000,'FaceColor', 'g')
title('Deviations in Y from source');
xlabel('deviation in meters');
ylabel('Num of deviations');
legend('Location','NorthWest');
grid on;

% XY deviations
figure;
%subplot(2,2,3)
plot(x_ls_dev(1,:), x_ls_dev(2,:), '*r','DisplayName','XY errors' );
%title('XY-Deviations plot from True Source position');
xlabel('X deviation from source(m)');
ylabel('Y-deviation from source(m)');
legend('Location','NorthWest');
grid on;

% --- Compute Euclidean Distances ---
% Deviations from True Source
distances_source = sqrt(x_ls_dev(1,:).^2 + x_ls_dev(2,:).^2);

% --- Compute CEP (50%) and 80% Radius ---
% For True Source Deviations
sorted_dist_source = sort(distances_source);
CEP50_source = prctile(sorted_dist_source, 50);  % 50% CEP
CEP80_source = prctile(sorted_dist_source, 80);  % 80% CEP

% --- Plot Scatter with CEP Circles ---
theta_dev = linspace(0, 2*pi, 100);

% (2) Deviations from True Source
%subplot(1,2,2);
figure;
%subplot(2,2,4);
scatter(x_ls_dev(1,:), x_ls_dev(2,:), '.b','DisplayName','XY errors');
hold on;

% Plot CEP50 (50% circle)
x_CEP50_source = CEP50_source * cos(theta_dev);
y_CEP50_source = CEP50_source * sin(theta_dev);
plot(x_CEP50_source, y_CEP50_source, 'r-', 'LineWidth', 2,'DisplayName','CEP50(50%)' );

% Plot CEP80 (80% circle)
x_CEP80_source = CEP80_source * cos(theta_dev);
y_CEP80_source = CEP80_source * sin(theta_dev);
plot(x_CEP80_source, y_CEP80_source, 'g--', 'LineWidth', 2, 'DisplayName','CEP80(80%)');

%title('Deviations from True Source with CEP');
xlabel('X Error(m)');
ylabel('Y Error(m)');
legend('Location','NorthWest');
grid on;
axis equal;

% --- Display CEP Values ---

fprintf('\n--- CEP (True Source Deviations) ---\n');
fprintf('CEP50 (50%% radius): %.4f\n', CEP50_source);
fprintf('CEP80 (80%% radius): %.4f\n', CEP80_source);
%% First subfigure
% Draw Figure
fig_geo = figure();hold on;

% Uncertainty Intervals
h = fill(lobFill0(1,:),lobFill0(2,:),'b','FaceAlpha',.15,'EdgeColor','none');
utils.excludeFromLegend(h);

h = fill(lobFill1(1,:),lobFill1(2,:),'b','FaceAlpha',.15,'EdgeColor','none');
utils.excludeFromLegend(h);

h = fill(lobFill2(1,:),lobFill2(2,:),'b','FaceAlpha',.15,'EdgeColor','none');
utils.excludeFromLegend(h);

% h = fill(lobFill3(1,:),lobFill3(2,:),'b','FaceAlpha',.15,'EdgeColor','none');
% utils.excludeFromLegend(h);


% Position Markers
hh=plot(x_sensor(1,:),x_sensor(2,:),'go','MarkerSize', 8,'DisplayName','Sensors');
utils.excludeFromLegend(hh);

% Position Labels
text(x_sensor(1,1)+155,x_sensor(2,1)-3,'S_0');
text(x_sensor(1,2)+155,x_sensor(2,2)-3,'S_1');
text(x_sensor(1,3)+155,x_sensor(2,3)-3,'S_2');
%text(x_sensor(1,4)+155,x_sensor(2,4)-3,'S_3');
% Plot the points and lobs
plot(x_source(1),x_source(2),'r*','MarkerSize', 8,'DisplayName','Transmitter');

% Geometric Solutions
%plot(x_centroid(1),x_centroid(2),'k*','DisplayName','Centroid');
%plot(x_incenter(1),x_incenter(2),'k+','DisplayName','Incenter');
%plot(x_ml(1,1),x_ml(2,1),'kv','DisplayName','Maximum Likelihood');
%plot(x_bf(1,1),x_bf(2,1),'ko','DisplayName','BestFix');

% Iterative Solutions
plot(x_initial(1),x_initial(2),'m*','MarkerSize', 7,'DisplayName','Initial Estimate');
plot(x_ls_full(1,:,10),x_ls_full(2,:,10),'b--','LineWidth',1.9 ,'DisplayName','Least Squares');
%plot(x_ls_full(1,numIters,10),x_ls_full(2,numIters,10),'ms','MarkerSize', 4 ,'DisplayName','LS-solution');
%utils.excludeFromLegend(hh);
plot(x_grad_full(1,:,10),x_grad_full(2,:,10),'g--','LineWidth',2.2,'DisplayName','Grad Descent');
%plot(x_grad_full(1,numIters,10),x_grad_full(2,numIters,10),'g^','MarkerSize', 4,'DisplayName','GD-solution' );
plot(x_gd_avg(1), x_gd_avg(2), 'g+','MarkerSize', 7,'DisplayName','Avg DG Sol' );
plot(x_ls_avg(1), x_ls_avg(2), 'bo','MarkerSize', 8,'DisplayName','Avg LS Sol' );
%utils.excludeFromLegend(hh);
%legend('Location','NorthWest');
xlabel('[m]');
ylabel('[m]');

% Compute CRLB and Error Ellipse
err_crlb = triang.computeCRLB(x_sensor,x_source,C_psi);
crlb_cep50 = utils.computeCEP50(err_crlb); % [km]
crlb_ellipse = utils.drawErrorEllipse(x_source,err_crlb,100,90);
plot(crlb_ellipse(1,:),crlb_ellipse(2,:),'r','LineWidth',1.0,'DisplayName','90% Error Ellipse');
%utils.excludeFromLegend(h);
%text(-20,45,'90\% Error Ellipse','FontSize',10);
%h=plot([1 11],[45 45],'k-','LineWidth',.5);
legend('Location','NorthEast');
utils.excludeFromLegend(h);
xlim([-40 50]*1e3);
ylim([-10 160]*1e3);

% Annotation Arrows
% annotation(fig_geo,'arrow',[0.635451388888888 0.709236111111111],...
%     [0.313924768518518 0.425601851851852]);
% annotation(fig_geo,'arrow',[0.559027777777778 0.506944444444444],...
%     [0.278645833333333 0.420572916666667]);
% 
% grid off;
% 
% utils.setPlotStyle(gca,{'tight'});


% %% Second subfigure
% fig_err=figure;
% loglog(1:numIters,cep50_ls,'m:','LineWidth',0.8,'DisplayName','Least-Squares');hold on;
% %text(1.2,4,'Least-Squares','FontSize',10);
% plot(1:numIters,cep50_grad,'b--','LineWidth',0.8,'DisplayName','Gradient Descent');
% %text(2.5,15,'Gradient Descent','FontSize',10);
% %plot(1:numIters,cep50_ml*ones(1,numIters),'DisplayName','Max Likelihood');
% %plot(1:numIters,cep50_bf*ones(1,numIters),'DisplayName','BestFix');
% %plot(1:numIters,cep50_cnt*ones(1,numIters),'DisplayName','Centroid');
% %text(15,2.6,'Centroid','FontSize',10);
% %plot(1:numIters,cep50_inc*ones(1,numIters),'DisplayName','Incenter');
% %text(3,35,'Incenter','FontSize',10);
% plot(1:numIters,crlb_cep50*ones(1,numIters),'g-.','LineWidth',0.8,'DisplayName','CRLB');
% %text(1.2,1.8,'CRLB','FontSize',10);
% 
% xlabel('Iteration Number');
% ylabel('$CEP_{50}$ [m]');
% legend('location','NorthEast');
% xlim([1 150]);
% ylim([1 100]*1e3);
% 
% utils.setPlotStyle(gca,{'widescreen','tight'});
