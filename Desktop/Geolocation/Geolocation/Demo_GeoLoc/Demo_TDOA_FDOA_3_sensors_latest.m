function fig = Demo_TDOA_FDOA_2_sensors()
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

%x_aoa = [4; 0]*1e3;
x_tdoa = [0, 40, 80; 0, 0, 0]*1e3;
x_fdoa = [0, 40, 80; 0, 0, 0]*1e3;
v_fdoa = [3, 1, 1; 1, 1, 0]*100 ;%*sqrt(.5)*100; % 100 m/s, at -45 deg heading

% Error Covariance Matrix
% err_aoa = 35; % deg
% cov_psi = (err_aoa*pi/180)^2; % rad^2

err_time = 20e-9; % 20 ns timing error
err_r = err_time * utils.constants.c;
cov_r = (err_r)^2*eye(size(x_tdoa,2)) ;% m^2, double for the combination of test/ref msmts
cov_r_out = utils.resampleCovMtx(cov_r,1);

freq_err = 5; % Hz
f0 = 1e9; % Hz
rr_err = freq_err * utils.constants.c/f0; % (m/s)
cov_rr = rr_err^2*eye(size(x_fdoa,2)) ;% (m/s)^2
cov_rr_out = utils.resampleCovMtx(cov_rr,1);

% Hybrid measurement and combined covariance matrix
%z = hybrid.measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source);
z = hybrid.measurement_TF(x_tdoa, x_fdoa, v_fdoa, x_source);
cov_z = blkdiag(cov_r, cov_rr); % raw sensor measurements

cov_z_out = blkdiag(cov_r_out, cov_rr_out); % sensor pairs

% cov_z_out = [35.9502, 0, 0,  0; 0, 35.9502,  0, 0; 0,  0, 2.2469, 0; 0, 0, 0,  2.2469];
% Generate Random Noise
L = chol(cov_z_out,'lower'); % Cholesky decomposition of the covariance matrix
%size(L,1);
% noise = L*randn(size(L,2),1);
% 
% % Noisy Measurements
% zeta = z + noise;

%% ML Search Parameters
x_ctr = [5; 80]*1e3;
grid_size = [90e3; 250e3];
grid_res = 600;  % meters, grid resolution

%% GD and LS Search Parameters
x_init = [15; 50]*1e3;
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
numSimulations = 1; % Number of Monte Carlo runs

% Initialize result storage
x_ml_results = zeros(2, numSimulations);
x_gd_results = zeros(2, numSimulations);
x_ls_results = zeros(2, numSimulations);

for i = 1:numSimulations
    % Generate Random Noise
    noise = L * randn(size(L,2),1);
    zeta_noisy = z + noise;
    
    % ML Solution
    x_ml_results(:, i) = hybrid.mlSoln_TF(x_tdoa, x_fdoa, v_fdoa, zeta_noisy, cov_z, x_ctr, grid_size, epsilon);
    
    % GD Solution
    x_gd_results(:, i) = hybrid.gdSoln_TF(x_tdoa, x_fdoa, v_fdoa, zeta_noisy, cov_z, x_init, [], [], epsilon, max_num_iterations, force_full_calc, plot_progress);
    
    % LS Solution
    x_ls_results(:, i) = hybrid.lsSoln_TF(x_tdoa, x_fdoa, v_fdoa, zeta_noisy, cov_z, x_init, epsilon, max_num_iterations, force_full_calc, plot_progress);
end

% Compute Averages
x_ml_avg = mean(x_ml_results, 2);
x_gd_avg = mean(x_gd_results, 2);
x_ls_avg = mean(x_ls_results, 2);

theta = rad2deg(atan2(x_source(2),x_source(1)));
theta_ml = rad2deg(atan2(x_ml_avg(2),x_ml_avg(1)));
theta_gd = rad2deg(atan2(x_gd_avg(2),x_gd_avg(1)));
theta_ls = rad2deg(atan2(x_ls_avg(2),x_ls_avg(1)));
% 
r_actual = norm(x_source - x_tdoa(:,1));
r_ml = norm(x_ml_avg - x_tdoa(:,1));
r_gd = norm(x_gd_avg - x_tdoa(:,1));
r_ls = norm(x_ls_avg - x_tdoa(:,1));

% Display Results
fprintf('\nActual Source Location is: (%.2f, %.2f)\n', x_source);
fprintf('\nAverage ML solution: (%.2f, %.2f)', x_ml_avg);
fprintf('\nAverage GD solution: (%.2f, %.2f)', x_gd_avg);
fprintf('\nAverage LS solution: (%.2f, %.2f)\n', x_ls_avg);

fprintf('\n Actal r and theta is %.2f,%.2f ', r_actual, theta)
fprintf('\n r and theta MLE pos estimate is: %.2f,%.2f', r_ml, theta_ml);
fprintf('\n r and theta from GD pos estimate is: %.2f,%.2f', r_gd,theta_gd);
fprintf('\n r and theta from LS pos estimate is: %.2f,%.2f', r_ls, theta_ls);

% Plot Monte Carlo Results
figure;
hold on;
%plot(x_source(1), x_source(2), 'kx', 'DisplayName','Target');

%plot(x_aoa(1), x_aoa(2), 'ko', 'DisplayName','AOA Sensor');
plot(x_tdoa(1, :), x_tdoa(2, :), 'ks', 'DisplayName','TDOA Sensor');
plot(x_fdoa(1, :), x_fdoa(2, :), 'k^', 'DisplayName','FDOA Sensor');
utils.drawArrow(x_fdoa(1,1)+[0 v_fdoa(1,1)],x_fdoa(2,1)+[0 v_fdoa(2,1)]);
utils.drawArrow(x_fdoa(1,2)+[0 v_fdoa(1,2)],x_fdoa(2,2)+[0 v_fdoa(2,2)]);

hold on;
plot(x_ml_results(1,:), x_ml_results(2,:), 'bv', 'MarkerSize', 6, 'DisplayName', 'ML Estimates');
plot(x_gd_results(1,:), x_gd_results(2,:), 'rx', 'MarkerSize', 6, 'DisplayName', 'GD Estimates');
plot(x_ls_results(1,:), x_ls_results(2,:), 'g*', 'MarkerSize', 6, 'DisplayName', 'LS Estimates');
plot(x_source(1), x_source(2), 'ko', 'MarkerSize', 8, 'DisplayName', 'True Source');

grid on;
ylim([-20 250]*1e3);
xlim([-20 100]*1e3);
caxis([-20 0]);
set(gca,'ydir','normal');
legend('Location','NorthEast');

xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Monte Carlo Simulation Results');



% fprintf('\n Actual Source location is: (%.2f,%.2f)', x_source);

%-------------------------------------------------------------------------------------
%clc; clear; close all;

%% Parameters
f = 162.025e6; % Signal frequency in Hz
lambda = 3e8 / f; % Wavelength
k = 2 * pi / lambda; % Wavenumber

source_pos = [3; 200]*1e3; % Source position (x;y) in meters
sensor_pos = [0 1000 2000; 0 0 0]; % Sensor positions (columns)

sigma_phase = 1; % Phase noise standard deviation in degrees
snapshots = 20; % Number of snapshots
coarse_angle = theta_ls ; % I have taken theta_ls for Coarse angle
%atand(source_pos(2) / source_pos(1)); % True angle
angle_error = 0.1; % Coarse estimate error in degrees

%% Compute Phase Differences
num_sensors = size(sensor_pos,2);
distances = sqrt(sum((sensor_pos - source_pos).^2, 1)); % Distance from source
phases = mod(k * distances, 2*pi); % True phase at each sensor
%phases_noisy = repmat(phases.', 1, snapshots) + deg2rad(sigma_phase) * randn(num_sensors, snapshots); % Add noise
% Step-by-step phase noise addition
phases = phases(:); % Ensure phases is a column vector
noise = deg2rad(sigma_phase) * randn(num_sensors, snapshots); % Generate noise matrix
phases_noisy = repmat(phases, 1, snapshots) + noise; % Add noise

%% Construct Covariance Matrix
X = exp(1j * phases_noisy); % Simulated received signal
R = (X * X') / snapshots ; % Sample covariance matrix

%% Apply MUSIC Algorithm
[U, ~, ~] = svd(R); % Eigenvalue decomposition
Un = U(:,2:end); % Noise subspace
angles = coarse_angle + linspace(-0.5, 0.5, 100); % Search space around coarse estimate
Pmusic = zeros(size(angles));

for i = 1:length(angles)
    steering_vector = exp(1j * k * sensor_pos(1,:) * cosd(angles(i))); % Use column-wise positions
    Pmusic(i) = 1 / (steering_vector * (Un * Un') * steering_vector');
end

Pmusic = abs(Pmusic) / max(abs(Pmusic)); % Normalize

%% Plot MUSIC Spectrum & Identify Peaks
figure;
plot(angles, 10*log10(Pmusic), 'LineWidth', 2);
xlabel('Angle (degrees)'); ylabel('Power (dB)');
grid on; title('MUSIC Spectrum');

%% Identify and Resolve Ambiguities
[peaks, locs] = findpeaks(Pmusic, 'SortStr', 'descend');  % Find peaks
angles_peaks = angles(locs);  % Map indices to angles
% angles_peaks = angles(locs);
% resolved_angle = angles_peaks(abs(angles_peaks - coarse_angle) < angle_error)
% Select angles within the error range of the coarse estimate
valid_idx = abs(angles_peaks - coarse_angle) < angle_error;
resolved_angle_MUSIC = angles_peaks(valid_idx)
resolved_peaks = peaks(valid_idx);

% Find the final resolved angle with the highest peak value
if ~isempty(resolved_peaks) 
    [~, max_idx] = max(resolved_peaks);  
    final_resolved_angle = resolved_angle_MUSIC(max_idx);
else
    final_resolved_angle = NaN; % No valid peak found
    warning('No peak found within the specified angle error range.');
end

disp(['Final Resolved Angle: ', num2str(final_resolved_angle)]);

hold on;
plot(angles_peaks, 10*log10(peaks), 'ro', 'MarkerSize', 8);
plot(resolved_angle_MUSIC, 10*log10(Pmusic(locs(abs(angles_peaks - coarse_angle) < angle_error))), 'gs', 'MarkerSize', 10);
legend('MUSIC Spectrum', 'Ambiguities', 'Resolved Angle');
% 
% 



%------------------------------------------------------------------------------------------------------------------------------------------------
% % Code for AOA
% x_source = [3; 200]*1e3;
% x_aoa = [0, 1; 0, 0]*1e3;
% err_delta_phi = 1*pi/180; % error in phase estimation 1deg
% f = 162.025e6;   c = 3*1e8;
% lam = c/f ;
% d = 1e3  ;
% d_lam = d/lam ;
% true_AOA = atan2((x_source(2,1)-x_aoa(2,1)),(x_source(1,1)-x_aoa(1,1))) ; % in radians
% phase_diff = 2*pi*d_lam*sin(pi/2 - true_AOA) ;
% wrapped_phase_diff = mod(phase_diff + pi, 2*pi) - pi ;
% psi_coarse = theta_ls; %(pi/180)*89.18 ;
% psi_est = [];
% k = floor(d_lam);
% n= k ;
% phi_diff_est = [];
% noise_delta_phi = 1*pi/180*randn(1);
% if mod(k,2) == 0
%     for i = 1:k
%         %fprintf("st")
%         phi_diff_est(i) = wrapped_phase_diff + 2*pi*((i-1) -n/2) ;
%         phi_diff_est(i) = phi_diff_est(i) + noise_delta_phi ;
%         
%     end
% else
%     for i = 1:k
%         %fprintf("st")
%         phi_diff_est(i) = wrapped_phase_diff + 2*pi*((i-1) -(n-1)/2) ;
%         phi_diff_est(i) = wrapped_phase_diff + noise_delta_phi;
%         
%     end
% end
% 
% 
% for i = 1:k
%     psi_est(i) = pi/2 - asin(phi_diff_est(i)/(2*pi*d_lam));
% end
% 
% resolved_psi_est = [];
% angle_error = 0.07; % Coarse estimate error in degrees
% for i = 1:k
%     if abs(psi_coarse -psi_est(i))< angle_error*pi/180 ;
%       resolved_psi_est = [resolved_psi_est, psi_est(i)*180/pi]; % Append value correctly
%     end
% end
% 
% resolved_psi_est
% 
% fprintf('\n Actal R and Angle is TRUE:(%.2f,\t%.2f) ', r_actual, theta);
% fprintf('\n R and Angle  from MLE: (%.2f,\t%.2f)', r_ml, theta_ml);
% fprintf('\n R and Angle from GD: (%.2f,\t%.2f)', r_gd,theta_gd);
% %fprintf('\n r and theta from LS pos estimate is: (%.2f,\t%.2f'), r_ls, theta_ls);
% fprintf('\n R from LS and Angle from AOA estimateis is: (%.2f,\t%.2f)', r_gd,theta_AOA);




