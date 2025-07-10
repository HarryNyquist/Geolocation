function fig = TDOA_LS_GD()
% fig=book2_ex2_2()
%
% Executes Example 2.2 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 22 April 2021

% Set up sensor and target coordinates
x_source = [3; 200]*1e3;

x_tdoa = [0, 25, 65, 100; 0, 0, 0 , 0]*1e3;
%x_fdoa = [0, 40, 80; 0, 0, 0]*1e3;

%v_fdoa = [3, 1, 1; 1, 1, 0]*100 ;%*sqrt(.5)*100; % 100 m/s, at -45 deg heading
%v_fdoa = [1, 1; 1, 1]*100 ;%*sqrt(.5)*100; % 100 m/s, at -45 deg heading


err_time = (20e-9)/sqrt(2); % 20 ns timing error
err_r = err_time * utils.constants.c;
cov_r = (err_r)^2*eye(size(x_tdoa,2)) ;% m^2, double for the combination of test/ref msmts
cov_r_out = utils.resampleCovMtx(cov_r,1);

freq_err = 10/sqrt(2); % Hz
f0 = 1e9; % Hz
%rr_err = freq_err * utils.constants.c/f0; % (m/s)
%cov_rr = rr_err^2*eye(size(x_fdoa,2)) ;% (m/s)^2
%cov_rr_out = utils.resampleCovMtx(cov_rr,1);

% Hybrid measurement and combined covariance matrix
%z = hybrid.measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source);
z = tdoa.measurement(x_tdoa, x_source);
cov_z = blkdiag(cov_r); % raw sensor measurements

cov_z_out = blkdiag(cov_r_out) % sensor pairs

%cov_z_out = [35.9502, 0, 0,  0; 0, 35.9502,  0, 0; 0,  0, 2.2469, 0; 0, 0, 0,  2.2469];
%cov_z_out = [35.9502, 0;  0, 2.2469];
% Generate Random Noise
L = chol(cov_z_out,'lower'); % Cholesky decomposition of the covariance matrix
%size(L,1);
% noise = L*randn(size(L,2),1);
% 
% % Noisy Measurements
% zeta = z + noise;

%% ML Search Parameters
x_ctr = [5; 150]*1e3;
grid_size = [90e3; 250e3];
grid_res = 400;  % meters, grid resolution

%% GD and LS Search Parameters
x_init = [5; 110]*1e3;
epsilon = grid_res;
max_num_iterations = 150;
force_full_calc = true;
plot_progress = false;

% %% ML Soln
% 
% x_ml = hybrid.mlSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, cov_z, ...
%        x_ctr, grid_size, epsilon);
% 
% 
% 
% theta = rad2deg(atan2(x_source(2),x_source(1)));
% theta_ml = rad2deg(atan2(x_ml(2),x_ml(1)));
% 
% %% GD Soln
% [x_gd, x_gd_full] = hybrid.gdSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, ...
%                     cov_z, x_init, [], [], epsilon, ... 
%                     max_num_iterations, force_full_calc, plot_progress);
%                 
% 
% 
% %% LS Soln
% [x_ls, x_ls_full] = hybrid.lsSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, ...
%                     cov_z, x_init, epsilon, max_num_iterations, ...
%                     force_full_calc, plot_progress);

%-------------------------------------------------------------------------------------------------------------------------------------------
% Monte Carlo Simulations Setup
numSimulations = 100; % Number of Monte Carlo runs

% Initialize result storage
x_ml_results = zeros(2, numSimulations);
x_gd_results = zeros(2, numSimulations);
x_ls_results = zeros(2, numSimulations);

for i = 1:numSimulations
    % Generate Random Noise
    noise = L * randn(size(L,2),1);
    zeta_noisy = z + noise;
    
    % ML Solution
    %x_ml_results(:, i) = tdoa.mlSoln_TF(x_tdoa, zeta_noisy, cov_z, x_ctr, grid_size, epsilon);
    
    % GD Solution
    x_gd_results(:, i) = tdoa.gdSoln(x_tdoa, zeta_noisy, cov_z, x_init, [], [], epsilon, max_num_iterations, force_full_calc, plot_progress,1);
                             % gdSoln(x_tdoa,rho,C,x_init,alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress,ref_idx)
    % LS Solution
    x_ls_results(:, i) = tdoa.lsSoln(x_tdoa, zeta_noisy, cov_z, x_init, epsilon, max_num_iterations, force_full_calc, plot_progress,1);
                            % lsSoln(x_tdoa,rho,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress,ref_idx)
end
% Compute Averages
%x_ml_avg = mean(x_ml_results, 2);
x_gd_avg = mean(x_gd_results, 2);
x_ls_avg = mean(x_ls_results, 2);
    
theta = rad2deg(atan2(x_source(2),x_source(1)));
%theta_ml = rad2deg(atan2(x_ml_avg(2),x_ml_avg(1)));
theta_gd = rad2deg(atan2(x_gd_avg(2),x_gd_avg(1)));
theta_ls = rad2deg(atan2(x_ls_avg(2),x_ls_avg(1)));
% 
r_actual = norm(x_source - x_tdoa(:,1));
%r_ml = norm(x_ml_avg - x_tdoa(:,1));
r_gd = norm(x_gd_avg - x_tdoa(:,1));
r_ls = norm(x_ls_avg - x_tdoa(:,1));

% Display Results
fprintf('\nActual Source Location is: (%.2f, %.2f)\n', x_source);
%fprintf('\nAverage ML solution: (%.2f, %.2f)', x_ml_avg);
fprintf('\nAverage GD solution: (%.2f, %.2f)', x_gd_avg);
fprintf('\nAverage LS solution: (%.2f, %.2f)\n', x_ls_avg);

fprintf('\n Actal r and theta is %.2f,%.2f ', r_actual, theta)
%fprintf('\n r and theta MLE pos estimate is: %.2f,%.2f', r_ml, theta_ml);
fprintf('\n r and theta from GD pos estimate is: %.2f,%.2f', r_gd,theta_gd);
fprintf('\n r and theta from LS pos estimate is: %.2f,%.2f', r_ls, theta_ls);

% Plot Monte Carlo Results
figure;
hold on;
%plot(x_source(1), x_source(2), 'kx', 'DisplayName','Target');

%plot(x_aoa(1), x_aoa(2), 'ko', 'DisplayName','AOA Sensor');
plot(x_tdoa(1, :), x_tdoa(2, :), 'ks', 'DisplayName','TDOA Sensor');


hold on;
%plot(x_ml_results(1,:), x_ml_results(2,:), 'bv', 'MarkerSize', 6, 'DisplayName', 'ML Estimates');
plot(x_gd_results(1,:), x_gd_results(2,:), 'rx', 'MarkerSize', 6, 'DisplayName', 'GD Estimates');
plot(x_ls_results(1,:), x_ls_results(2,:), 'g*', 'MarkerSize', 6, 'DisplayName', 'LS Estimates');
plot(x_source(1), x_source(2), 'ko', 'MarkerSize', 8, 'DisplayName', 'True Source');

grid on;
ylim([-20 250]*1e3);
xlim([-20 120]*1e3);
%zlim([-20 100]*1e3);
caxis([-20 0]);
set(gca,'ydir','normal');
legend('Location','NorthEast');

xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Monte Carlo Simulation Results');

% Compute Deviations
x_ls_dev = [];
x_ls_dev_x_dir = [];
x_ls_dev_y_dir = [] ;
%Compute Deviations-- Vectorized version
x_ls_dev = x_ls_avg - x_ls_results;
x_ls_dev_x_dir = x_ls_avg(1) - x_ls_results(1,:); % Deviations from mean solution
x_ls_dev_y_dir = x_ls_avg(2) - x_ls_results(2,:); % Deviations from mean solution
x_ls_dev_x_dir_source = x_source(1) - x_ls_results(1,:);  % Deviations from True solution
x_ls_dev_y_dir_source = x_source(2) - x_ls_results(2,:); % Deviations from True solution
% Compute absolute deviations (magnitude only, no sign)
abs_dev_x = abs(x_ls_dev_x_dir);
abs_dev_y = abs(x_ls_dev_y_dir);

abs_dev_x_source = abs(x_ls_dev_x_dir_source);
abs_dev_y_source = abs(x_ls_dev_y_dir_source);

% Plot bar charts
figure;

% X deviations
subplot(2,2,1);
%bar(abs_dev_x, 'b');
histogram(x_ls_dev_x_dir_source, 'BinWidth', 20,'FaceColor', 'b')
title('Deviations in X from source');
xlabel('deviation in meters');
ylabel('Num of deviations');
grid on;

% Y deviations
subplot(2,2,2);
%bar(abs_dev_y, 'r');
histogram(x_ls_dev_y_dir_source, 'BinWidth', 20,'FaceColor', 'g')
title('Deviations in Y from source');
xlabel('deviation in meters');
ylabel('Num of deviations');
grid on;

% XY deviations
subplot(2,2,3);
plot(x_ls_dev_x_dir_source, x_ls_dev_y_dir_source, '*r');
title('XY-Deviations plot from True Source position');
xlabel('X deviation from source pos');
ylabel('Y-deviation from source pos');
grid on;

% --- Compute Euclidean Distances ---
% Deviations from Mean
distances_mean = sqrt(x_ls_dev_x_dir.^2 + x_ls_dev_y_dir.^2);
% Deviations from True Source
distances_source = sqrt(x_ls_dev_x_dir_source.^2 + x_ls_dev_y_dir_source.^2);

% --- Compute CEP (50%) and 80% Radius ---
% For Mean Deviations
sorted_dist_mean = sort(distances_mean);
CEP50_mean = prctile(sorted_dist_mean, 50);  % 50% CEP
CEP80_mean = prctile(sorted_dist_mean, 80);  % 80% CEP

% For True Source Deviations
sorted_dist_source = sort(distances_source);
CEP50_source = prctile(sorted_dist_source, 50);  % 50% CEP
CEP80_source = prctile(sorted_dist_source, 80);  % 80% CEP

% --- Plot Scatter with CEP Circles ---
theta_dev = linspace(0, 2*pi, 100);


% (2) Deviations from True Source
subplot(2,2,4);
scatter(x_ls_dev_x_dir_source, x_ls_dev_y_dir_source, 'b', 'filled');
hold on;

% Plot CEP50 (50% circle)
x_CEP50_source = CEP50_source * cos(theta_dev);
y_CEP50_source = CEP50_source * sin(theta_dev);
plot(x_CEP50_source, y_CEP50_source, 'r-', 'LineWidth', 2);

% Plot CEP80 (80% circle)
x_CEP80_source = CEP80_source * cos(theta_dev);
y_CEP80_source = CEP80_source * sin(theta_dev);
plot(x_CEP80_source, y_CEP80_source, 'g--', 'LineWidth', 2);

title('Deviations from True Source with CEP');
xlabel('X Deviation');
ylabel('Y Deviation');
legend('Deviations', 'CEP50 (50%)', 'CEP80 (80%)', 'Location', 'best');
grid on;
axis equal;

% --- Display CEP Values ---

fprintf('\n--- CEP (True Source Deviations) ---\n');
fprintf('CEP50 (50%% radius): %.4f\n', CEP50_source);
fprintf('CEP80 (80%% radius): %.4f\n', CEP80_source);

%-------------------------------------------------------------------------------------