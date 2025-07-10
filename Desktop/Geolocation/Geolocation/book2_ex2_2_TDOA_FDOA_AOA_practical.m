function fig = book2_ex2_2()
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
x_source = [3; 50]*1e3;

x_aoa = [4; 0]*1e3;
x_tdoa = [0, 12, 12, 0; 0, 0, 12, 12]*1e3;
x_fdoa = [0, 12, 12, 0; 0, 0, 12, 12]*1e3;
v_fdoa = [1, 1, 1, 1; -1, -1, -1, -1]*sqrt(.5)*100; % 100 m/s, at -45 deg heading

% Error Covariance Matrix
err_aoa = 5; % deg
cov_psi = (err_aoa*pi/180)^2; % rad^2

err_time = 20e-9; % 20 ns timing error
err_r = err_time * utils.constants.c;
cov_r = (err_r)^2*eye(size(x_tdoa,2)); % m^2, double for the combination of test/ref msmts
cov_r_out = utils.resampleCovMtx(cov_r,1);

freq_err = 5; % Hz
f0 = 1e9; % Hz
rr_err = freq_err * utils.constants.c/f0; % (m/s)
cov_rr = rr_err^2*eye(size(x_fdoa,2)); % (m/s)^2
cov_rr_out = utils.resampleCovMtx(cov_rr,1);

% Hybrid measurement and combined covariance matrix
z = hybrid.measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source);
cov_z = blkdiag(cov_psi, cov_r, cov_rr); % raw sensor measurements
cov_z_out = blkdiag(cov_psi, cov_r_out, cov_rr_out); % sensor pairs

% Generate Random Noise
L = chol(cov_z_out,'lower'); % Cholesky decomposition of the covariance matrix
size(L,1);
noise = L*randn(size(L,2),1);

% Noisy Measurements
zeta = z + noise;

%% ML Search Parameters
x_ctr = [5; 25]*1e3;
grid_size = [20e3; 60e3];
grid_res = 400;  % meters, grid resolution

%% GD and LS Search Parameters
x_init = [1; 1]*1e3;
epsilon = grid_res;
max_num_iterations = 1500;
force_full_calc = true;
plot_progress = false;

%% ML Soln

x_ml = hybrid.mlSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, cov_z, ...
       x_ctr, grid_size, epsilon);

% numSimulations = 10;
% 
% % Initialize the result storage
% x_ml_results = zeros(2, numSimulations); % Assuming x_ml is a 3x1 vector
% 
% % Run Monte Carlo simulations
% for i = 1:numSimulations
%     % Generate random noise (adjust according to your problem)
%     %noise = mvnrnd(zeros(size(zeta)), cov_z)'; % Assuming Gaussian noise
%     
%     % Run the hybrid.mlSoln function with noisy input
%     x_ml_results(:, i) = hybrid.mlSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, cov_z, ...
%        x_ctr, grid_size, epsilon);
% end
% 
% % Compute the average result
% x_ml_avg = mean(x_ml_results, 2);
% 
% % Display results
% disp('A single x_ml is:');
% disp(x_ml)
% disp('\n Average x_ml over 100 Monte Carlo simulations:');
% disp(x_ml_avg);

theta = rad2deg(atan2(x_source(2),x_source(1)));
theta_ml = rad2deg(atan2(x_ml(2),x_ml(1)));

%% GD Soln
[x_gd, x_gd_full] = hybrid.gdSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, ...
                    cov_z, x_init, [], [], epsilon, ... 
                    max_num_iterations, force_full_calc, plot_progress);
                
% numSimulations = 10;
% 
% % Initialize the result storage
% x_gd_results = zeros(2, numSimulations); % Assuming x_gd is a 3x1 vector
% 
% % Run Monte Carlo simulations
% for i = 1:numSimulations
%     % Generate random noise (adjust according to your problem)
%     %noise = mvnrnd(zeros(size(zeta)), cov_z)'; % Assuming Gaussian noise
%     
%     % Run the hybrid.mlSoln function with noisy input
%     x_gd_results(:, i) = hybrid.gdSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, ...
%                     cov_z, x_init, [], [], epsilon, ... 
%                     max_num_iterations, force_full_calc, plot_progress);
% end
% 
% % Compute the average result
% x_gd_avg = mean(x_gd_results, 2);
% disp('\n Average x_gd over 10 Monte Carlo simulations:');
% disp(x_gd_avg);
% disp('\n the 10 x_gd are:');
% disp(x_gd_results);

%% LS Soln
[x_ls, x_ls_full] = hybrid.lsSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, ...
                    cov_z, x_init, epsilon, max_num_iterations, ...
                    force_full_calc, plot_progress);

%% Plot Result
theta = rad2deg(atan2(x_source(2),x_source(1)));
theta_ml = rad2deg(atan2(x_ml(2),x_ml(1)));
theta_gd = rad2deg(atan2(x_gd(2),x_gd(1)));
theta_ls = rad2deg(atan2(x_ls(2),x_ls(1)));

r_actual = norm(x_source - x_tdoa(:,1));
r_ml = norm(x_ml - x_tdoa(:,1));
r_gd = norm(x_gd - x_tdoa(:,1));
r_ls = norm(x_ls - x_tdoa(:,1));

fprintf('\n Actual Source location is: (%.2f,%.2f)', x_source);
fprintf('\n ML solution is: (%.2f, %.2f)', x_ml);
fprintf('\n GD solution is: (%.2f,%.2f)', x_gd);
fprintf('\n LS solution is: (%.2f,%.2f)', x_ls);

fprintf('\n Actal r and theta is %.2f,%.2f ', r_actual, theta)
fprintf('\n r and theta MLE pos estimate is: %.2f,%.2f', r_ml, theta_ml);
fprintf('\n r and theta from GD pos estimate is: %.2f,%.2f', r_gd,theta_gd);
fprintf('\n r and theta from LS pos estimate is: %.2f,%.2f', r_ls, theta_ls);

fig=figure;
plot(x_source(1), x_source(2), 'kx', 'DisplayName','Target');
hold on;
plot(x_aoa(1), x_aoa(2), 'ko', 'DisplayName','AOA Sensor');
plot(x_tdoa(1, :), x_tdoa(2, :), 'ks', 'DisplayName','TDOA Sensor');
plot(x_fdoa(1, :), x_fdoa(2, :), 'k^', 'DisplayName','FDOA Sensor');
utils.drawArrow(x_fdoa(1,1)+[0 v_fdoa(1,1)],x_fdoa(2,1)+[0 v_fdoa(2,1)]);
utils.drawArrow(x_fdoa(1,2)+[0 v_fdoa(1,2)],x_fdoa(2,2)+[0 v_fdoa(2,2)]);

plot(x_ml(1), x_ml(2), 'v', 'DisplayName', 'ML Solution', 'MarkerSize',12);
hdl=plot(x_gd_full(1,:), x_gd_full(2,:), '-.', 'LineWidth',2);
utils.excludeFromLegend(hdl);
plot(x_gd(1),x_gd(2),'-.x','DisplayName','GD Solution','Color',hdl.Color, 'MarkerSize',13);
hdl=plot(x_ls_full(1,:), x_ls_full(2,:), '-', 'LineWidth', 2);
utils.excludeFromLegend(hdl);
plot(x_ls(1), x_ls(2), '-*','DisplayName','LS Solution','Color',hdl.Color,  'MarkerSize', 12);

grid on;
ylim([0 70]*1e3);
xlim([-5 30]*1e3);
caxis([-20 0]);
set(gca,'ydir','normal');
legend('Location','NorthEast');
%utils.setPlotStyle(gca,{'widescreen','tight'});
