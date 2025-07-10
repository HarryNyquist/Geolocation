clc; clear; close all;

%% Parameters
fc = 300e6;               % Carrier frequency (Hz)
c = 3e8;                  % Speed of light (m/s)
lambda = c / fc;          % Wavelength (m)
d = 2000 * lambda;        % Inter-element spacing (m)
M = 4;                    % Number of elements (2x2 array flattened)
angles = [-30, 20];       % DOA angles in degrees
SNR = 20;                 % Signal-to-Noise Ratio (dB)
Nsources = length(angles);
Nsnapshots = 200;         % Number of snapshots

%% Generate Steering Matrix
array_positions = [0, d, d, 0; 0, 0, d, d]; % 2x2 Grid
array_positions = array_positions(:); % Flatten to column vector
A = zeros(M, Nsources) + 1j * zeros(M, Nsources); % Ensure complex allocation
for k = 1:Nsources
    theta = angles(k) * pi / 180;
    A(:, k) = exp(-1j * (2 * pi / lambda) * array_positions * sin(theta));
end

%% Debugging Output
disp('Size of A:'), disp(size(A));
disp('Size of array_positions:'), disp(size(array_positions));
disp('Size of lambda:'), disp(size(lambda));
disp('Size of theta:'), disp(size(theta));

%% Generate Signals
signals = randn(Nsources, Nsnapshots) + 1j * randn(Nsources, Nsnapshots);
received_signal = A * signals;

%% Add Noise
noise = (randn(M, Nsnapshots) + 1j * randn(M, Nsnapshots)) * 10^(-SNR/20);
X = received_signal + noise;

%% Apply Phase Measurement Error
phase_error = (1 * pi / 180) * (randn(M, 1)); % 1-degree error per channel
X = X .* exp(1j * phase_error);

%% Covariance Matrix and Eigen Decomposition
Rxx = (X * X') / Nsnapshots;
[E, D] = eig(Rxx);
[eigenvalues, idx] = sort(diag(D), 'descend');
E = E(:, idx);
En = E(:, Nsources+1:end); % Noise subspace

%% MUSIC Spectrum Computation
theta_scan = -90:0.1:90;
Pmusic = zeros(size(theta_scan));
for i = 1:length(theta_scan)
    sv = exp(-1j * (2 * pi / lambda) * array_positions * sin(theta_scan(i) * pi / 180));
    Pmusic(i) = 1 / (sv' * (En * En') * sv);
end
Pmusic = abs(Pmusic);
Pmusic = 10 * log10(Pmusic / max(Pmusic));

%% Plot MUSIC Spectrum
figure;
plot(theta_scan, Pmusic, 'LineWidth', 2);
xlabel('Angle (degrees)');
ylabel('MUSIC Spectrum (dB)');
title('MUSIC Algorithm Spectrum');
grid on;
axis([-90 90 -40 0]);
