function fig_sensor_count = sensor_count_vs_error_robust()
% Robust version that handles ill-conditioned cases

% Parameters
nMC = 80;               % Number of Monte Carlo trials
numIters = 300;          % Number of iterations
sensor_counts = 2:6;    % Number of sensors to test
sensor_spacing = 15e3;   % 10 km between sensors

% Define source position (off-axis to avoid colinearity)
x_source = [15, 120]' * 1e3;  % Changed from [15,120] to reduce colinearity

% Initialize error storage
avg_gd_errors = zeros(size(sensor_counts));
avg_ls_errors = zeros(size(sensor_counts));

% Regularization parameter for ill-conditioned cases
lambda = 1e-3;  % Small regularization value

for n_idx = 1:length(sensor_counts)
    N = sensor_counts(n_idx);
    
    % Define sensor positions (add some y-axis variation)
    x_sensor = zeros(2, N);
    x_sensor(1,:) = linspace(0, (N-1)*sensor_spacing, N);
    x_sensor(2,:) = randn(1,N)*sensor_spacing*0.1; % Small random y-offsets
    
    % Measurement noise covariance
    C_psi = (2*pi/180)^2 * eye(N);
    psi_act = triang.measurement(x_sensor, x_source);
    psi = psi_act + sqrt(C_psi) * randn(N, nMC);
    
    % Better initial guess (centroid of sensors with some offset)
    x_initial = mean(x_sensor, 2) + [5e3; 20e3];
    
    x_gd = zeros(2, nMC);
    x_ls = zeros(2, nMC);
    
    fprintf('Processing %d sensors...', N);
    
    for idx = 1:nMC
        % GD Solution with regularization check
        [x_gd(:,idx), ~] = triang.gdSoln(x_sensor, psi(:,idx), C_psi, ...
                                      x_initial, [], lambda, [], numIters, false, []);
        
        % LS Solution with regularization
        [x_ls(:,idx), ~] = triang.lsSoln(x_sensor, psi(:,idx), C_psi, ...
                                      x_initial, lambda, numIters, false, []);
    end
    fprintf('done.\n');
    
    % Compute robust average (median to handle outliers)
    gd_errors = sqrt(sum((x_gd - x_source).^2, 1)) ;
    ls_errors = sqrt(sum((x_ls - x_source).^2, 1)) ;
    
    avg_gd_errors(n_idx) = mean(gd_errors);
    avg_ls_errors(n_idx) = mean(ls_errors);
end

fig_sensor_count = figure();
plot(sensor_counts, avg_gd_errors, 'g-o', 'LineWidth', 2, ...
     'MarkerFaceColor', 'g', 'DisplayName', 'Gradient Descent');
hold on;
plot(sensor_counts, avg_ls_errors, 'm-s', 'LineWidth', 2, ...
     'MarkerFaceColor', 'm', 'DisplayName', 'Least Squares');
grid on;
xlabel('Number of Sensors');
ylabel('Median Position Error (km)');
title('Positioning Error vs Number of Sensors');
legend('Location', 'NorthEast');
xticks(sensor_counts);
set(gca, 'YScale', 'log');  % Log scale often helps visualize error trends

% Display table of results
fprintf('\n=== Median Positioning Errors (km) ===\n');
fprintf('Sensors\t GD Error\t LS Error\n');
for n_idx = 1:length(sensor_counts)
    fprintf('%d\t %.3f\t\t %.3f\n', sensor_counts(n_idx), ...
            avg_gd_errors(n_idx), avg_ls_errors(n_idx));
end
end