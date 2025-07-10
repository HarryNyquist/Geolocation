%% Figure 9/11 - Beamscan Images
% Set up array
N = 3;
f = 162.025e6; % Signal frequency in Hz
lam = 3e8 / f; % Wavelength
d = 1000 ; 
k = 2 * pi /lam; % Wavenumber
d_lam = d/lam ;
% d_psi = .89/((N-1)*d_lam);
% du = sin(d_psi);
%v = array.make_steering_vector(d_lam,N);
v = @(psi) exp(1i*2*pi*d_lam*(0:N-1).*sin(psi(:))).';

v_dot = @(psi) (1i*2*pi*d_lam*(0:N-1).*cos(psi(:))).'.*v(psi);

% % Set up source signals
% spacing = 1;
% u1 = du*spacing/2;
%u2 = -u1;
source_pos = [3; 200]*1e3; % Source position (x;y) in meters
sensor_pos = [0 1000 2000; 0 0 0]; % Sensor positions (columns)


psi_actual = atan2(source_pos(2,:),source_pos(1,:));
psi1  = pi/2 - psi_actual ;%asin(u1);
%psi2 = asin(u2);

v1 = v(psi1); % N x numel(u1)
%v2 = v(psi2); % N x numel(u2)

% Generate snapshots
M = 30;
s1 = sqrt(.5)*(randn(1,M)+1i*randn(1,M));
%s2 = sqrt(.5)*(randn(1,M)+1i*randn(1,M));
x1 = v1.*s1;
%x2 = v2.*s2;

% Add noise
snr_db = 15;
snr = 10.^(snr_db/10);
n = sqrt(1./(snr*2))*(randn(size(x1))+1i*randn(size(x1)));
x = x1 + n ; %+ x2;

psi_LS = (90 -89.20)*pi/180 ; 
fig9=figure;hold on;
[P,psi_vec] = array.beamscan_modified(x,v,psi_LS,501);


P_music = array.music_modified(x,v,2,psi_LS,501);
P_music = P_music./max(abs(P_music(:)));
set(gca,'LineStyleOrderIndex',3);
%plot(sin(psi_vec),10*log10(abs(P_music)),'DisplayName','MUSIC');
plot((psi_vec)*180/pi,10*log10(abs(P_music)),'DisplayName','MUSIC');
utils.setPlotStyle(gca,{'widescreen'});
%utils.exportPlot(fig9,[prefix '11']);




