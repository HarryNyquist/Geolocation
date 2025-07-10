% MATLAB Implementation of MUSIC Algorithm for Linear Array of 4 Elements

% Parameters
c = 3e8; % Speed of light (m/s)
f = 300e6; % Frequency (300 MHz)
lambda = c / f; % Wavelength (m)
d = lambda / 2; % Inter-element spacing (half-wavelength)
theta_range = -90:0.1:90; % Range of angles to search (degrees)
num_elements = 4; % Number of antenna elements in the linear array
num_signals = 2; % Number of incoming signals
phase_error = 5; % Phase measurement error in degrees

% Array Geometry (Linear Array)
array_positions = (0:num_elements-1) * d; % Positions of the array elements along the x-axis

% Simulated Signals
theta_true = [-30, 36]; % True DOAs of the signals (degrees)
A = exp(1i * 2 * pi * array_positions' * sind(theta_true) / lambda); % Steering matrix (4x2)
noise_power = 0.1; % Noise power
signal_power = [1; 1]; % Signal power (2x1)
X = A * diag(sqrt(signal_power)) + sqrt(noise_power) * (randn(num_elements, 1) + 1i * randn(num_elements, 1)); % Received signal (4x1)

% Add phase measurement error
phase_error_matrix = exp(1i * deg2rad(phase_error) * randn(size(X)));
X = X .* phase_error_matrix;

% MUSIC Algorithm
R = X * X' / size(X, 2); % Covariance matrix (4x4)
[E, D] = eig(R); % Eigenvalue decomposition
[~, idx] = sort(diag(D), 'descend'); % Sort eigenvalues in descending order
E = E(:, idx); % Sort eigenvectors accordingly
En = E(:, num_signals+1:end); % Noise subspace eigenvectors (4x2)

% MUSIC Spectrum
music_spectrum = zeros(size(theta_range));
for i = 1:length(theta_range)
    a = exp(1i * 2 * pi * array_positions' * sind(theta_range(i)) / lambda); % Steering vector (4x1)
    denominator = abs(a' * (En * En') * a); % Ensure this is a scalar
    music_spectrum(i) = 1 / denominator; % MUSIC spectrum
end

% Plot MUSIC Spectrum
figure;
plot(theta_range, 10 * log10(abs(music_spectrum)));
xlabel('Angle (degrees)');
ylabel('MUSIC Spectrum (dB)');
title('MUSIC Spectrum for Linear Array of 4 Elements');
grid on;

% Find the peaks in the MUSIC spectrum
[~, peak_indices] = findpeaks(abs(music_spectrum), 'SortStr', 'descend', 'NPeaks', num_signals);
estimated_thetas = theta_range(peak_indices);
fprintf('Estimated DOAs: %.2f degrees and %.2f degrees\n', estimated_thetas(1), estimated_thetas(2));