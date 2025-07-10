f      = 1e9;             % Carrier frequency (Hz)
lambda = utils.constants.c / f;           % Wavelength (m)
M      = 200;             % Number of snapshots
SNR_dB = 20;              % Signal-to-noise ratio per sensor (dB)

%% 1) True source
psi_true = 30 * pi/180;               % True azimuth (radians)
x_source = [cos(psi_true); sin(psi_true)];  % Unit‐direction vector

%% 2) Sensor array (3 sensors in 2D)
sensor_pos = [0,   40,  80;            % x-coordinates (m)
              0,    0,   0];           % y-coordinates (m)
N = size(sensor_pos,2);

%% 3) Build steering‐vector function handle
%    v(psi) returns an N×1 vector of exp(j*2π/λ * (sensor·u(psi)))
v = @(psi) exp( ...
          1j * 2*pi/lambda * ...
          ( sensor_pos(1,:).*cos(psi) + sensor_pos(2,:).*sin(psi) )'  ...
        );

%% 4) Simulate snapshots x = v(ψ_true)*s + noise
s = exp(1j*2*pi*rand(1,M));            % random QPSK‐like phases
A = v(psi_true) / sqrt(N);             % normalized steering
signal_power = 1;
noise_power  = signal_power / (10^(SNR_dB/10));
noise_sigma  = sqrt(noise_power/2);

x = A * s + noise_sigma*(randn(N,M) + 1j*randn(N,M));

%% 5) Run MUSIC
D       = 1;       % number of sources
max_psi = 90*pi/180;
N_pts   = 361;     % one‐degree resolution
[P,psi_vec] = array.music(x, v, D, max_psi, N_pts);

%% 6) Plot spatial spectrum
figure;
plot(psi_vec*180/pi, 10*log10(P), 'LineWidth',1.5);
hold on;
xlabel('Angle \psi (deg)');
ylabel('MUSIC Spectrum (dB)');
title('MUSIC Spatial Spectrum');
grid on;

% mark true angle
[~, idx_true] = min(abs(psi_vec - psi_true));
stem(psi_true*180/pi, 10*log10(P(idx_true)), 'r','filled','LineWidth',1.5);
legend('P(\psi)','True \psi');

fprintf('True DOA = %.2f°\n', psi_true*180/pi);
