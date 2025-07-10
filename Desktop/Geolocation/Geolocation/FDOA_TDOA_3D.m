function fig = FDOA_TDOA_3D()
% Set up sensor and target coordinates
x_source = [10; 500; 500]*1e3;

x_tdoa = [0, 40, 80; 0, 0, 0; 0, 0, 0]*1e3; %three TDOA sensors, each at (0,0), (40k, 0), (80k, 0)
x_fdoa = [0, 40, 80; 0, 0, 0; 0, 0, 0]*1e3; 
v_fdoa = [3, 1, 1; 1, 1, 0; 2, 1, 1]*100 

err_time = (20e-9)/sqrt(2); % 20 ns timing error
err_r = err_time * utils.constants.c;
cov_r = (err_r)^2*eye(size(x_tdoa,2)) ;% m^2, double for the combination of test/ref msmts
cov_r_out = utils.resampleCovMtx(cov_r,1);

freq_err = 5/sqrt(2); % Hz
f0 = 1e9; % Hz
rr_err = freq_err * utils.constants.c/f0; % (m/s)
cov_rr = rr_err^2*eye(size(x_fdoa,2)) ;% (m/s)^2
cov_rr_out = utils.resampleCovMtx(cov_rr,1);

z = hybrid.measurement_TF(x_tdoa, x_fdoa, v_fdoa, x_source);
cov_z = blkdiag(cov_r, cov_rr); % raw sensor measurements

cov_z_out = blkdiag(cov_r_out, cov_rr_out) % sensor pairs

% Generate Random Noise
L = chol(cov_z_out,'lower'); % Cholesky decomposition of the covariance matrix

%% GD and LS Search Parameters
x_init = [5; 50; 50]*1e3;
epsilon = 200; %meters, grid resolution
max_num_iterations = 200;
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
x_ml_results = zeros(3, numSimulations);
x_gd_results = zeros(3, numSimulations);
x_ls_results = zeros(3, numSimulations);
for i = 1:numSimulations
    % Generate Random Noise
    noise = L * randn(size(L,2),1);
    zeta_noisy = z + noise;
    
    % ML Solution
    %x_ml_results(:, i) = hybrid.mlSoln_TF(x_tdoa, x_fdoa, v_fdoa, zeta_noisy, cov_z, x_ctr, grid_size, epsilon);
    
    % GD Solution
    x_gd_results(:, i) = hybrid.gdSoln_TF(x_tdoa, x_fdoa, v_fdoa, zeta_noisy, cov_z, x_init, [], [], epsilon, max_num_iterations, force_full_calc, plot_progress);
    
    % LS Solution
    x_ls_results(:, i) = hybrid.lsSoln_TF(x_tdoa, x_fdoa, v_fdoa, zeta_noisy, cov_z, x_init, epsilon, max_num_iterations, force_full_calc, plot_progress);
end
% Compute Averages
%x_ml_avg = mean(x_ml_results, 2);
x_gd_avg = mean(x_gd_results, 2);
x_ls_avg = mean(x_ls_results, 2);

 
% Display Results
fprintf('\nActual Source Location is: (%.2f, %.2f, %.2f)\n', x_source);
%fprintf('\nAverage ML solution: (%.2f, %.2f)', x_ml_avg);
fprintf('\nAverage GD solution: (%.2f, %.2f, %.2f)', x_gd_avg);
fprintf('\nAverage LS solution: (%.2f, %.2f, %.2f)\n', x_ls_avg);

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
plot(x_ls_dev_x_dir, x_ls_dev_y_dir, '*b');
title('XY-deviations plot from Mean');
xlabel('X deviation(m)');
ylabel('Y-deviation(m)');
grid on;

% XY deviations
subplot(2,2,4);
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

% (1) Deviations from Mean
figure;
subplot(1,2,1);
scatter(x_ls_dev_x_dir, x_ls_dev_y_dir, 'b', 'filled');
hold on;

% Plot CEP50 (50% circle)
x_CEP50_mean = CEP50_mean * cos(theta_dev);
y_CEP50_mean = CEP50_mean * sin(theta_dev);
plot(x_CEP50_mean, y_CEP50_mean, 'r-', 'LineWidth', 2);

% Plot CEP80 (80% circle)
x_CEP80_mean = CEP80_mean * cos(theta_dev);
y_CEP80_mean = CEP80_mean * sin(theta_dev);
plot(x_CEP80_mean, y_CEP80_mean, 'g--', 'LineWidth', 2);

title('Deviations from Mean with CEP');
xlabel('X Deviation');
ylabel('Y Deviation');
legend('Deviations', 'CEP50 (50%)', 'CEP80 (80%)', 'Location', 'best');
grid on;
axis equal;

% (2) Deviations from True Source
subplot(1,2,2);
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
fprintf('--- CEP (Mean Deviations) ---\n');
fprintf('CEP50 (50%% radius): %.4f\n', CEP50_mean);
fprintf('CEP80 (80%% radius): %.4f\n', CEP80_mean);

fprintf('\n--- CEP (True Source Deviations) ---\n');
fprintf('CEP50 (50%% radius): %.4f\n', CEP50_source);
fprintf('CEP80 (80%% radius): %.4f\n', CEP80_source);

%-------------------------------------------------------------------------------------
%clc; clear; close all;

%% Parameters
% f = 162.025e6; % Signal frequency in Hz
% lambda = 3e8 / f; % Wavelength
% k = 2 * pi / lambda; % Wavenumber
% 
% source_pos = [10; 500]*1e3;% Source position (x;y) in meters 
% %dont know why this was declared as a seperate variable.
% 
% sensor_pos = [0 1000 2000; 0 0 0]; % Sensor positions (columns)
% 
% sigma_phase = 1; % Phase noise standard deviation in degrees
% snapshots = 20; % Number of snapshots
% coarse_angle = theta_ls ; % I have taken theta_ls for Coarse angle
% %atand(source_pos(2) / source_pos(1)); % True angle
% angle_error = 0.1; % Coarse estimate error in degrees
% 
% %% Compute Phase Differences
% num_sensors = size(sensor_pos,2);
% distances = sqrt(sum((x_ls_avg - source_pos).^2, 1)); % Distance from source
% phases = mod(k * distances, 2*pi); % True phase at each sensor
% %phases_noisy = repmat(phases.', 1, snapshots) + deg2rad(sigma_phase) * randn(num_sensors, snapshots); % Add noise
% % Step-by-step phase noise addition
% phases = phases(:); % Ensure phases is a column vector
% noise = deg2rad(sigma_phase) * randn(num_sensors, snapshots); % Generate noise matrix
% phases_noisy = repmat(phases, 1, snapshots) + noise; % Add noise
% 
% %% Construct Covariance Matrix
% X = exp(1j * phases_noisy); % Simulated received signal
% R = (X * X') / snapshots ; % Sample covariance matrix
% 
% %% Apply MUSIC Algorithm
% [U, ~, ~] = svd(R); % Eigenvalue decomposition
% Un = U(:,2:end); % Noise subspace
% angles = coarse_angle + linspace(-0.5, 0.5, 100); % Search space around coarse estimate
% Pmusic = zeros(size(angles));
% 
% for i = 1:length(angles)
%     steering_vector = exp(1j * k * sensor_pos(1,:) * cosd(angles(i))); % Use column-wise positions
%     Pmusic(i) = 1 / (steering_vector * (Un * Un') * steering_vector');
% end
% 
% Pmusic = abs(Pmusic) / max(abs(Pmusic)); % Normalize
% 
% %% Plot MUSIC Spectrum & Identify Peaks
% figure;
% plot(angles, 10*log10(Pmusic), 'LineWidth', 2);
% xlabel('Angle (degrees)'); ylabel('Power (dB)');
% grid on; title('MUSIC Spectrum');
% 
% %% Identify and Resolve Ambiguities
% [peaks, locs] = findpeaks(Pmusic, 'SortStr', 'descend');  % Find peaks
% angles_peaks = angles(locs);  % Map indices to angles
% % angles_peaks = angles(locs);
% % resolved_angle = angles_peaks(abs(angles_peaks - coarse_angle) < angle_error)
% % Select angles within the error range of the coarse estimate
% valid_idx = abs(angles_peaks - coarse_angle) < angle_error;
% resolved_angle_MUSIC = angles_peaks(valid_idx)
% resolved_peaks = peaks(valid_idx);
% 
% % Find the final resolved angle with the highest peak value
% if ~isempty(resolved_peaks) 
%     differences = abs(resolved_angle_MUSIC - theta_ls);
%     % Find index of the minimum difference
%     [~, idx] = min(differences);
%     % Select the fine measurement closest to coarse measurement
%     final_measurement = resolved_angle_MUSIC(idx);
% 
% else
%     final_resolved_angle = NaN; % No valid peak found
%     warning('No peak found within the specified angle error range.');
% end
% 
% disp(['Final Resolved Angle: ', num2str(final_measurement)]);
% 
% hold on;
% plot(angles_peaks, 10*log10(peaks), 'ro', 'MarkerSize', 8);
% plot(resolved_angle_MUSIC, 10*log10(Pmusic(locs(abs(angles_peaks - coarse_angle) < angle_error))), 'gs', 'MarkerSize', 10);
% legend('MUSIC Spectrum', 'Ambiguities', 'Resolved Angle');
% % 
% % 
% 
% 
% end

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



