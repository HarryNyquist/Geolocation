function fig = Demo_FDOA_3_sensors_only()

% Source location
x_source = [10; 500]*1e3;

% FDOA sensor locations and velocities
x_fdoa = [0, 40, 80; 0, 0, 0]*1e3;
v_fdoa = [3, 1, 1; 1, 1, 0]*100;

% Frequency error
freq_err = 5/sqrt(2); % Hz
f0 = 1e9; % Hz
rr_err = freq_err * utils.constants.c/f0; % m/s
cov_rr = rr_err^2 * eye(size(x_fdoa, 2)); % covariance of raw sensor msmt
cov_rr_out = utils.resampleCovMtx(cov_rr, 1); % pairwise covariance

% Measurement and noise
z = fdoa.measurement(x_fdoa, v_fdoa, x_source);
cov_z = cov_rr;
cov_z_out = cov_rr_out;

L = chol(cov_z_out, 'lower'); % For generating Gaussian noise

% Search parameters
x_ctr = [5; 80]*1e3;
grid_size = [90e3; 250e3];
grid_res = 600;

x_init = [15; 50]*1e3;
epsilon = grid_res;
max_num_iterations = 200;
force_full_calc = true;
plot_progress = false;

% Monte Carlo setup
numSimulations = 100;
x_gd_results = zeros(2, numSimulations);
x_ls_results = zeros(2, numSimulations);

for i = 1:numSimulations
    noise = L * randn(size(L,2), 1);
    zeta_noisy = z + noise;

    x_gd_results(:, i) = fdoa.gdSoln(x_fdoa, v_fdoa, zeta_noisy, cov_z, ...
        x_init, [], [], epsilon, max_num_iterations, force_full_calc, plot_progress);

    x_ls_results(:, i) = fdoa.lsSoln(x_fdoa, v_fdoa, zeta_noisy, cov_z, ...
        x_init, epsilon, max_num_iterations, force_full_calc, plot_progress);
end

% Averages
x_gd_avg = mean(x_gd_results, 2);
x_ls_avg = mean(x_ls_results, 2);

theta = rad2deg(atan2(x_source(2), x_source(1)));
theta_gd = rad2deg(atan2(x_gd_avg(2), x_gd_avg(1)));
theta_ls = rad2deg(atan2(x_ls_avg(2), x_ls_avg(1)));

r_actual = norm(x_source - x_fdoa(:,1));
r_gd = norm(x_gd_avg - x_fdoa(:,1));
r_ls = norm(x_ls_avg - x_fdoa(:,1));

% Display
fprintf('\nActual Source Location: (%.2f, %.2f)\n', x_source);
fprintf('\nAverage GD solution: (%.2f, %.2f)', x_gd_avg);
fprintf('\nAverage LS solution: (%.2f, %.2f)\n', x_ls_avg);
fprintf('\nActual r, theta: %.2f, %.2f\n', r_actual, theta);
fprintf('\nGD estimate r, theta: %.2f, %.2f\n', r_gd, theta_gd);
fprintf('\nLS estimate r, theta: %.2f, %.2f\n', r_ls, theta_ls);

% Plot Monte Carlo Results
figure;
hold on;
plot(x_fdoa(1, :), x_fdoa(2, :), 'k^', 'DisplayName', 'FDOA Sensor');
utils.drawArrow(x_fdoa(1,1)+[0 v_fdoa(1,1)], x_fdoa(2,1)+[0 v_fdoa(2,1)]);
utils.drawArrow(x_fdoa(1,2)+[0 v_fdoa(1,2)], x_fdoa(2,2)+[0 v_fdoa(2,2)]);
plot(x_gd_results(1,:), x_gd_results(2,:), 'rx', 'DisplayName', 'GD Estimates');
plot(x_ls_results(1,:), x_ls_results(2,:), 'g*', 'DisplayName', 'LS Estimates');
plot(x_source(1), x_source(2), 'ko', 'MarkerSize', 8, 'DisplayName', 'True Source');

grid on;
ylim([-20 250]*1e3);
xlim([-20 120]*1e3);
set(gca,'ydir','normal');
legend('Location','NorthEast');
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('FDOA-Only Monte Carlo Simulation Results');

% Deviations
x_ls_dev_x_dir_source = x_source(1) - x_ls_results(1,:);
x_ls_dev_y_dir_source = x_source(2) - x_ls_results(2,:);

% Histograms
figure;
subplot(2,2,1);
histogram(x_ls_dev_x_dir_source, 'BinWidth', 20, 'FaceColor', 'b');
title('Deviations in X from Source');
xlabel('Deviation (m)'); ylabel('Count'); grid on;

subplot(2,2,2);
histogram(x_ls_dev_y_dir_source, 'BinWidth', 20, 'FaceColor', 'g');
title('Deviations in Y from Source');
xlabel('Deviation (m)'); ylabel('Count'); grid on;

subplot(2,2,3);
plot(x_ls_dev_x_dir_source, x_ls_dev_y_dir_source, '*r');
title('XY-Deviations from True Source');
xlabel('X Deviation'); ylabel('Y Deviation'); grid on;

% CEP
distances_source = sqrt(x_ls_dev_x_dir_source.^2 + x_ls_dev_y_dir_source.^2);
CEP50_source = prctile(distances_source, 50);
CEP80_source = prctile(distances_source, 80);

theta_dev = linspace(0, 2*pi, 100);
x_CEP50 = CEP50_source * cos(theta_dev);
y_CEP50 = CEP50_source * sin(theta_dev);
x_CEP80 = CEP80_source * cos(theta_dev);
y_CEP80 = CEP80_source * sin(theta_dev);

subplot(2,2,4);
scatter(x_ls_dev_x_dir_source, x_ls_dev_y_dir_source, 'b', 'filled');
hold on;
plot(x_CEP50, y_CEP50, 'r-', 'LineWidth', 2);
plot(x_CEP80, y_CEP80, 'g--', 'LineWidth', 2);
title('CEP Circles (from True Source)');
xlabel('X Deviation'); ylabel('Y Deviation'); grid on; axis equal;

fprintf('\n--- CEP (True Source Deviations) ---\n');
fprintf('CEP50: %.2f m\n', CEP50_source);
fprintf('CEP80: %.2f m\n', CEP80_source);
end