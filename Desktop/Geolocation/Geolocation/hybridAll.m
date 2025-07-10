x_source = [10; 500]*1e3; %source at (100, 200) km

tdoa_sensors = [0, 40, 80; 0, 0, 0]*1e3; % sensors at (30,12), (50,27) km 
%just randomly chosen numbers btw

%we'll put the final sensor as the reference. 
err_time = (20e-9)/sqrt(2);
err_r = err_time * utils.constants.c;
%cov_r = (err_r)^2*eye(size(tdoa_sensors,2)) ;%NxN covariance matrix of individual sensor measurements
cov_t = (err_time)^2*eye(size(tdoa_sensors,2));

%experiments go brrrrr

cov_r = [(err_r)^2, 0, 0; 0, 20*(err_r)^2, 0; 0, 0, 50*(err_r)^2];
%for the above, CEP50 ~ 5100m and CEP90 ~ 12600m


%cov_r = [50*(err_r)^2, 0, 0; 0, 20*(err_r)^2, 0; 0, 0, (err_r)^2];
%on 22-06-25, CEP50 ~ 5200m and CEP90 ~ 12600m 



%cov_r = [5*(err_r)^2, 0, 0; 0, 5*(err_r)^2, 0; 0, 0, 5*(err_r)^2];
%for the above, CEP50 ~ 2500m and CEP90 ~ 6000m

%cov_r = [2*(err_r)^2, 0, 0; 0, 2*(err_r)^2, 0; 0, 0, 2*(err_r)^2];
%for the above, CEP50 ~ 1550m and CEP90 ~ 4000m

%cov_r = [1e-12, 0, 0; 0, (err_r)^2, 0; 0, 0, 1e-12];
%for the above, MATLAB gives up due to badly conditioned matrix.

%cov_r = [1e-7, 0, 0; 0, 2*(err_r)^2, 0; 0, 0, 1e-7]; %err_r ~ 4
%for the above, CEP50 ~ 635m and CEP90 ~ 1520m

%cov_r = [2*(err_r)^2, 0, 0; 0, 1e-7, 0; 0, 0, 1e-7];
%for the above, CEP50 ~ 635m and CEP90 ~ 1580m

%cov_r = [1e-9, 0, 0; 0, 1e-9, 0; 0, 0, 1e-9];
%for the above, CEP50 ~ 0.01m and CEP90 ~ 0.02m

%this concludes the sanity check, looks to be good.

cov_r_out = utils.resampleCovMtx(cov_r,size(tdoa_sensors,2));
%this is to make final sensor as reference.

L = chol(cov_r_out, 'lower'); %find L s.t cov_r_out = LL*, where L* is herimitian 

noiseless_measurement = tdoa.measurement(tdoa_sensors, x_source);

max_num_iterations = 100;
force_full_calc = true;
plot_progress = false;
%we only care about LS solution either ways, soo:dd
x_init = [5; 50]*1e3;
epsilon = 100; %m;
numSimulations = 1000;
x_ls_results = zeros(2, numSimulations);
for i = 1:numSimulations
    noise = L * randn(size(L,2),1);
    zeta_noisy = noiseless_measurement + noise;
    % LS Solution
    x_ls_results(:, i) = tdoa.lsSoln(tdoa_sensors, zeta_noisy, cov_r, x_init);
    %for a sanity check
    %Update: Sanity check workssss!!!
    %x_ls_results(:, i) = tdoa.lsSoln(tdoa_sensors, noiseless_measurement, cov_r, x_init);
end
% Compute Averages
x_ls_avg = mean(x_ls_results, 2);
% Display Results
fprintf('\nActual Source Location is: (%.2f, %.2f)\n', x_source);
fprintf('\nAverage LS solution: (%.2f, %.2f)\n', x_ls_avg);

% --- Compute Deviations (you already have these) ---
x_ls_dev_x_dir     = x_ls_avg(1)          - x_ls_results(1,:);
x_ls_dev_y_dir     = x_ls_avg(2)          - x_ls_results(2,:);
x_ls_dev_x_src_dir = x_source(1)          - x_ls_results(1,:);
x_ls_dev_y_src_dir = x_source(2)          - x_ls_results(2,:);

% --- Euclidean Distances ---
dist_mean   = sqrt(x_ls_dev_x_dir.^2     + x_ls_dev_y_dir.^2);
dist_source = sqrt(x_ls_dev_x_src_dir.^2 + x_ls_dev_y_src_dir.^2);

% --- Compute CEP50 & CEP90 ---
CEP50_mean   = prctile(dist_mean,   50);
CEP90_mean   = prctile(dist_mean,   90);
CEP50_source = prctile(dist_source, 50);
CEP90_source = prctile(dist_source, 90);

% --- Scatter + CEP Circles ---
theta = linspace(0,2*pi,200);

figure;
% (1) From Mean
subplot(1,2,1);
scatter(x_ls_dev_x_dir, x_ls_dev_y_dir, 10, 'filled'); hold on;
plot( CEP50_mean*cos(theta),   CEP50_mean*sin(theta),   'r-',  'LineWidth',1.5 );
plot( CEP90_mean*cos(theta),   CEP90_mean*sin(theta),   'm--', 'LineWidth',1.5 );
axis equal; grid on;
title('Deviations from Mean with CEP');
xlabel('X deviation (m)'); ylabel('Y deviation (m)');
legend('Errors','CEP50','CEP90','Location','best');

% (2) From True Source
subplot(1,2,2);
scatter(x_ls_dev_x_src_dir, x_ls_dev_y_src_dir, 10, 'filled'); hold on;
plot( CEP50_source*cos(theta), CEP50_source*sin(theta), 'r-',  'LineWidth',1.5 );
plot( CEP90_source*cos(theta), CEP90_source*sin(theta), 'm--', 'LineWidth',1.5 );
axis equal; grid on;
title('Deviations from True Source with CEP');
xlabel('X deviation (m)'); ylabel('Y deviation (m)');
legend('Errors','CEP50','CEP90','Location','best');

% --- Print CEPs ---
fprintf('\n--- CEP (Mean Deviations) ---\n');
fprintf('CEP50: %.2f m\n', CEP50_mean);
fprintf('CEP90: %.2f m\n', CEP90_mean);

fprintf('\n--- CEP (True Source Deviations) ---\n');
fprintf('CEP50: %.2f m\n', CEP50_source);
fprintf('CEP90: %.2f m\n', CEP90_source);

% crlb = tdoa.computeCRLB(tdoa_sensors, x_source, cov_t);
% cep50 = utils.computeCEP50(crlb);
% fprintf('CEP50 calculated from book function: %.4f m\n', cep50);


theta_ls = rad2deg(atan2(x_ls_avg(2),x_ls_avg(1))); 
theta_true = rad2deg(atan2(x_source(2), x_source(1)));
coarse_angle = theta_ls; %this is in degrees, just so we r clear.
f = 162.056e6;
lambda = utils.constants.c/f;
k = 2*pi/lambda;
aoa_sensors = [0, 1000, 2000; 0, 0, 0];
d = sqrt((aoa_sensors(1,1) - aoa_sensors(2,1))^2 + (aoa_sensors(1,2) - aoa_sensors(2,2))^2);
%it should be a UNIFORM LINEAR ARRAY!! 
distances = sqrt(sum((x_source - aoa_sensors).^2, 1)); % Distance from source
phases = mod(k * distances, 2*pi); % True phase at each sensor
phases = phases(:);
sigma_phase = 0.01; %1 degree std dev.
snapshots = 20;
num_sensors = size(aoa_sensors, 2);
noise = deg2rad(sigma_phase) * randn(num_sensors, snapshots);
%here, both noise and X are 3*20 vectors
phases_noisy = repmat(phases, 1, snapshots) + noise;
x    = exp(1j * phases_noisy);
angle_error = 0.1;
%For some reason, couldnt make the built in functions work.
% d_lam = d/lambda;
% v = array.make_steering_vector(d_lam, num_sensors);
% 
% min_psi = coarse_angle - 0.5;
% max_psi = coarse_angle + 0.5;
% 
% P_music = array.music_with_bounds(x, v, 1, min_psi, max_psi, 100);
% P_music = abs(P_music)/max(abs(P_music));
% angle_error = 0.1;
% % Reconstruct the angle axis:
% angles = linspace(min_psi, max_psi, numel(P_music));


%Fine, I'll do it myself
R = (x * x') / snapshots ; % Sample covariance matrix

%% Apply MUSIC Algorithm
[U, ~, ~] = svd(R); % Eigenvalue decomposition
Un = U(:,2:end); % Noise subspace
angles = coarse_angle + linspace(-0.5, 0.5, 100); % Search space around coarse estimate
P_music = zeros(size(angles));

for i = 1:length(angles)
    steering_vector = exp(1j * k * aoa_sensors(1,:) * cosd(angles(i))); % Use column-wise positions
    P_music(i) = 1 / (steering_vector * (Un * Un') * steering_vector');
end

P_music = abs(P_music) / max(abs(P_music)); % Normalize

%% Plot MUSIC Spectrum & Identify Peaks
figure;
plot(angles, 10*log10(P_music), 'LineWidth', 2);
xlabel('Angle (degrees)'); ylabel('Power (dB)');
grid on; title('MUSIC Spectrum');

%% Identify and Resolve Ambiguities
[peaks, locs] = findpeaks(P_music, 'SortStr', 'descend');  % Find peaks
angles_peaks = angles(locs);  % Map indices to angles
% angles_peaks = angles(locs);
% resolved_angle = angles_peaks(abs(angles_peaks - coarse_angle) < angle_error)
% Select angles within the error range of the coarse estimate
valid_idx = abs(angles_peaks - coarse_angle) < angle_error;
resolved_angle_MUSIC = angles_peaks(valid_idx);
resolved_peaks = peaks(valid_idx);

% Find the final resolved angle with the highest peak value
if ~isempty(resolved_peaks) 
    differences = abs(resolved_angle_MUSIC - theta_ls);
    % Find index of the minimum difference
    [~, idx] = min(differences);
    % Select the fine measurement closest to coarse measurement
    final_measurement = resolved_angle_MUSIC(idx);

else
    final_resolved_angle = NaN; % No valid peak found
    warning('No peak found within the specified angle error range.');
end

disp(['Final Resolved Angle: ', num2str(final_measurement)]);

hold on;
plot(angles_peaks, 10*log10(peaks), 'ro', 'MarkerSize', 8);
plot(resolved_angle_MUSIC, 10*log10(P_music(locs(abs(angles_peaks - coarse_angle) < angle_error))), 'gs', 'MarkerSize', 10);
legend('MUSIC Spectrum', 'Ambiguities', 'Resolved Angle');
% 
% 





%----------------------------------------------------------------
% f = 162.025e6; % Signal frequency in Hz
% lambda = 3e8 / f; % Wavelength
% k = 2 * pi / lambda; % Wavenumber
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
% distances = sqrt(sum((x_source - sensor_pos).^2, 1)); % Distance from source
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
% resolved_angle_MUSIC = angles_peaks(valid_idx);
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








%------------------------------------------------------------------



