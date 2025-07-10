function [fig_error, fig_scatter] = sensor_count_plots()
% Generates two plot sets compatible with pre-R2019b MATLAB versions

% Parameters
nMC = 80;                % Monte Carlo trials
numIters = 300;          % Iterations
sensor_counts = 3;     % Sensor counts to test
sensor_spacing = 30e3;   % Distance between sensors
lambda = 1e-3;           % Regularization parameter

% Source position (kept constant)
x_source = [15, 120]' * 1e3;

% Initialize storage
avg_gd_errors = zeros(size(sensor_counts));
avg_ls_errors = zeros(size(sensor_counts));
all_gd_solutions = cell(length(sensor_counts),1);
all_ls_solutions = cell(length(sensor_counts),1);

%% Main simulation loop
for n_idx = 1:length(sensor_counts)
    N = sensor_counts(n_idx);
    
    % Sensor positions with small y-variation
    x_sensor = zeros(2, N);
    x_sensor(1,:) = linspace(0, (N-1)*sensor_spacing, N);
    x_sensor(2,:) = randn(1,N)*sensor_spacing*0.1;
    
    % Measurement setup
    C_psi = (2*pi/180)^2 * eye(N);
    psi_act = triang.measurement(x_sensor, x_source);
    psi = psi_act + sqrt(C_psi) * randn(N, nMC);
    
    % Initial guess near centroid
    x_initial = [5e3; 60e3];
    
    % Storage for this sensor count
    x_gd = zeros(2, nMC);
    x_ls = zeros(2, nMC);
    
    fprintf('Processing %d sensors...\n', N);
    for idx = 1:nMC
        [x_gd(:,idx), ~] = triang.gdSoln(x_sensor, psi(:,idx), C_psi, ...
                          x_initial, [], lambda, [], numIters, false, []);
        [x_ls(:,idx), ~] = triang.lsSoln(x_sensor, psi(:,idx), C_psi, ...
                          x_initial, lambda, numIters, false, []);
    end
    
    % Store solutions and errors
    all_gd_solutions{n_idx} = x_gd;
    all_ls_solutions{n_idx} = x_ls;
    avg_gd_errors(n_idx) = mean(sqrt(sum((x_gd - x_source).^2, 1)));
    avg_ls_errors(n_idx) = mean(sqrt(sum((x_ls - x_source).^2, 1)));
end

%% Plot 1: Error vs Number of Sensors
fig_error = figure('Position', [100 100 800 400]);

%plot(sensor_counts, avg_gd_errors/1e3, 'g-o', 'LineWidth', 2, ...
     %'MarkerFaceColor', 'g', 'DisplayName', 'GD');
%hold on;
plot(sensor_counts, avg_ls_errors/1e3, 'm-s', 'LineWidth', 2, ...
     'MarkerFaceColor', 'm', 'DisplayName', 'LS');
grid on;
xlabel('Number of Sensors');
ylabel('Mean Error (km)');
title('(a) Positioning Error vs Sensor Count');
legend('Location', 'NorthEast');
xticks(sensor_counts);

% subplot(1,2,2);
% %semilogy(sensor_counts, avg_gd_errors/1e3, 'g-o', 'LineWidth', 2, ...
%          'MarkerFaceColor', 'g', 'DisplayName', 'GD');
% %hold on;
% semilogy(sensor_counts, avg_ls_errors/1e3, 'm-s', 'LineWidth', 2, ...
%          'MarkerFaceColor', 'm', 'DisplayName', 'LS');
% grid on;
% xlabel('Number of Sensors');
% ylabel('Mean Error (km) - Log Scale');
% title('(b) Error Trend (Log Scale)');
% xticks(sensor_counts);

%% Plot 2: Scatter Plots for Each Sensor Count (using subplot instead of tiledlayout)
fig_scatter = figure('Position', [100 100 1200 800]);
rows = 2;
cols = ceil(length(sensor_counts)/rows);

for n_idx = 1:length(sensor_counts)
    subplot(rows, cols, n_idx);
    N = sensor_counts(n_idx);
    
    % Plot GD solutions
    %scatter(all_gd_solutions{n_idx}(1,:)/1e3, all_gd_solutions{n_idx}(2,:)/1e3, ...
            %40, 'g', 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'GD');
    %hold on;
    
    % Plot LS solutions
    scatter(all_ls_solutions{n_idx}(1,:)/1e3, all_ls_solutions{n_idx}(2,:)/1e3, ...
            40, 'm', 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'LS');
    
    % Plot true source
    plot(x_source(1)/1e3, x_source(2)/1e3, 'rp', 'MarkerSize', 12, ...
         'MarkerFaceColor', 'r', 'DisplayName', 'True Source');
    
    % Formatting
    grid on;
    title(sprintf('%d Sensors', N));
    xlabel('X (km)');
    ylabel('Y (km)');
    axis equal;
    
    % Add error ellipses if function exists
    if exist('error_ellipse', 'file')
        cov_gd = cov(all_gd_solutions{n_idx}');
        error_ellipse(cov_gd, x_source/1e3, 'conf', 0.95, 'style', 'g--');
        cov_ls = cov(all_ls_solutions{n_idx}');
        error_ellipse(cov_ls, x_source/1e3, 'conf', 0.95, 'style', 'm--');
    end
    
    if n_idx == 1
        legend('Location', 'best');
    end
end

%% Display results
fprintf('\n=== Mean Positioning Errors ===\n');
fprintf('Sensors\t GD Error (km)\t LS Error (km)\n');
for n_idx = 1:length(sensor_counts)
    fprintf('%d\t %.2f\t\t %.2f\n', sensor_counts(n_idx), ...
            avg_gd_errors(n_idx)/1e3, avg_ls_errors(n_idx)/1e3);
end
end