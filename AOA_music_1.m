clc; clear; close all;

%% Parameters
f = 162.025e6; % Signal frequency in Hz
lambda = 3e8 / f; % Wavelength
k = 2 * pi / lambda; % Wavenumber

source_pos = [3000; 200000]; % Source position (x;y) in meters
sensor_pos = [0 500 1000; 0 0 0]; % Sensor positions (columns)

sigma_phase = 1; % Phase noise standard deviation in degrees
snapshots = 20; % Number of snapshots
coarse_angle = 89.12; %atand(source_pos(2) / source_pos(1)); % True angle
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
angles = coarse_angle + linspace(-1, 1, 1000); % Search space around coarse estimate
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
[peaks, locs] = findpeaks(Pmusic, 'SortStr', 'descend');
angles_peaks = angles(locs);
%resolved_angle = angles_peaks(abs(angles_peaks - coarse_angle) < angle_error);
valid_indices = abs(angles_peaks - coarse_angle) < angle_error;
disp(valid_indices);
resolved_angle = angles_peaks(valid_indices)
hold on;
plot(angles_peaks, 10*log10(peaks), 'ro', 'MarkerSize', 8);
plot(resolved_angle, 10*log10(Pmusic(locs(abs(angles_peaks - coarse_angle) < angle_error))), 'gs', 'MarkerSize', 10);
legend('MUSIC Spectrum', 'Ambiguities', 'Resolved Angle');
