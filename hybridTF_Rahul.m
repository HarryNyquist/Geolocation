x_source = [10; 500]*1e3;

x_tdoa = [0, 40, 80; 0, 0, 0]*1e3; % sensors at (30,12), (50,27) km 
x_fdoa = [0, 40, 80; 0, 0, 0]*1e3; 
v_fdoa = [3, 1, 1; 1, 1, 0]*100;

err_time = (20e-9)/sqrt(2);
err_r = err_time * utils.constants.c;
%err_r = 1e-6;
cov_r = (err_r)^2*eye(size(x_tdoa,2)) ;% m^2, double for the combination of test/ref msmts
cov_r_out = utils.resampleCovMtx(cov_r,1);

freq_err = 5/sqrt(2);
%freq_err = 1e-5;
f0 = 1e9; % 1 GHz
rr_err = freq_err * utils.constants.c/f0; % (m/s)
%rr_err = 1e-6
cov_rr = rr_err^2*eye(size(x_fdoa,2)) ;% (m/s)^2
cov_rr_out = utils.resampleCovMtx(cov_rr,1);

noiseless_measurement = hybrid.measurement_TF(x_tdoa, x_fdoa, v_fdoa, x_source);
cov_z = blkdiag(cov_r, cov_rr); 
cov_z_out = blkdiag(cov_r_out, cov_rr_out); 

% Generate Random Noise
L = chol(cov_z_out,'lower'); 

x_init = [5; 50]*1e3;
epsilon = 100;
max_num_iterations = 100;
force_full_calc = true;
plot_progress = false;

%-------------------------------------------------------------------------------------------------------------------------------------------
% Monte Carlo Simulations Setup
numSimulations = 100; % Number of Monte Carlo runs

x_ls_results = zeros(2, numSimulations);
for i = 1:numSimulations
    noise = L * randn(size(L,2),1);
    zeta_noisy = noiseless_measurement + noise;
    x_ls_results(:, i) = hybrid.lsSoln_TF(x_tdoa, x_fdoa, v_fdoa, zeta_noisy, cov_z, x_init, epsilon, max_num_iterations, force_full_calc, plot_progress);
    %x_ls_results(:, i) = hybrid.lsSoln_TF(x_tdoa, x_fdoa, v_fdoa, noiseless_measurement, cov_z, x_init, epsilon, max_num_iterations, force_full_calc, plot_progress);
end

x_ls_avg = mean(x_ls_results, 2);   
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
sigma_phase = 1; %1 degree std dev.
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
[N,M] = size(x);
C = zeros(N,N);
for idx_m = 1:M
    C = C + x(:,idx_m)*x(:,idx_m)';
end
C = C/M;

% Perform Eigendecomposition of C
[U,Lam] = eig(C);
lam = diag(Lam);

% Sort the eigenvalues
[~,idx_sort] = sort(abs(lam),'descend');
U_sort = U(:,idx_sort);
lam_sort = lam(idx_sort);
Un = U_sort(:,2:end); % Noise subspace
angles = coarse_angle + linspace(-0.5, 0.5, 100); % Search space around coarse estimate
%how to choose this search space is another question atp.

P_music = zeros(size(angles));

for i = 1:length(angles)
    steering_vector = exp(1j * k * aoa_sensors(1,:) * cosd(angles(i)))/sqrt(N); % Use column-wise positions
    P_music(i) = 1 ./ abs(steering_vector * (Un * Un') * steering_vector');
    %changes in MUSIC made to reflect book codes.
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
    %changed so that it marks the max value as the peak. 
    [~, idx] = max(resolved_peaks);  
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




