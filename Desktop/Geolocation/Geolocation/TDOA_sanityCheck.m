x_source = [10; 500]*1e3; %source at (100, 200) km

%tdoa_sensors = [30, 50, 85; 12, 27, 105]*1e3; % sensors at (30,12), (50,27) km 
tdoa_sensors = [0, 40, 80; 0, 0, 0];

err_time = 20e-9;
err_time = 20e-9/sqrt(2); 
err_r = err_time * utils.constants.c;
cov_t = (err_time)^2 * eye(size(tdoa_sensors,2));
cov_r = (err_r)^2 * eye(size(tdoa_sensors, 2)); 
cov_t_out = utils.resampleCovMtx(cov_t, 1);
cov_r_out = utils.resampleCovMtx(cov_r, 1);

noiseless_measurement = tdoa.measurement(tdoa_sensors, x_source);
L = chol(cov_r_out,'lower');

x_init = [5; 50]*1e3;
%x_init = [50; 50]*1e3;
max_num_iterations = 200;
force_full_calc = true;
plot_progress = false;
epsilon = 600; %m
numSimulations = 100;

x_ls_results = zeros(2, numSimulations);
for i = 1:numSimulations
    noise = L * randn(size(L,2),1);
    %basically, cov_r_out = L*L';
    %noise = L*z, where z is a standard normal vector (2D)
    %cov(noise) = E(noise*noise') = L*E(z*z')*L' = L*L' = cov_r_out
    %So what this basically does is, it creates a vector with cov_r_out as
    %covariance and 0 mean.
    zeta_noisy = noiseless_measurement + noise;
    % LS Solution
    %x_ls_results(:, i) = tdoa.lsSoln(tdoa_sensors, zeta_noisy, cov_r, x_init,1);
    %for a sanity check
    %Update: Sanity check workssss!!!
    x_ls_results(:, i) = tdoa.lsSoln(tdoa_sensors, noiseless_measurement, cov_r, x_init);
end
x_ls_avg = mean(x_ls_results, 2);
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
fprintf('The LS solution is: %.4f\n', x_ls_avg);

% crlb = tdoa.computeCRLB(tdoa_sensors, x_source, cov_t);
% cep50 = utils.computeCEP50(crlb);
% fprintf('CEP50 calculated from book function: %.4f m\n', cep50);