%% Figure 9/11 - Beamscan Images
% Set up array
N = 25;
d_lam = .5;
d_psi = .89/((N-1)*d_lam);
du = sin(d_psi);
v = array.make_steering_vector(d_lam,N);

% Set up source signals
spacing = 1;
u1 = du*spacing/2;
u2 = -u1;

psi1 = asin(u1);
psi2 = asin(u2);

v1 = v(psi1); % N x numel(u1)
v2 = v(psi2); % N x numel(u2)

% Generate snapshots
M = 30;
s1 = sqrt(.5)*(randn(1,M)+1i*randn(1,M));
s2 = sqrt(.5)*(randn(1,M)+1i*randn(1,M));
x1 = v1.*s1;
x2 = v2.*s2;

% Add noise
snr_db = 10;
snr = 10.^(snr_db/10);
n = sqrt(1./(snr*2))*(randn(size(x1))+1i*randn(size(x1)));
x = x1 + x2 + n;


% Generate Beamscan images
[P,psi_vec] = array.beamscan(x,v,pi/2,1001);
P_mvdr = array.beamscan_mvdr(x,v,pi/2,1001);

% Scale outputs
P = P/max(P(:));
P_mvdr = P_mvdr/max(P_mvdr(:));

fig9=figure;hold on;
plot(sin(psi_vec),10*log10(P),'LineWidth',1.5,'DisplayName','Beamscan');
set(gca,'LineStyleOrderIndex',2);
plot(sin(psi_vec),10*log10(P_mvdr),'LineWidth',1.25,'DisplayName','MVDR');
%h1=plot(u1*[1 1],[-60 0],'k--');
%h2=plot(u2*[1 1],[-60 0],'k--');
%utils.excludeFromLegend([h1 h2]);
ylim([-30 0]);
xlabel('u');
ylabel('P [dB]');
legend('Location','NorthEast');
utils.setPlotStyle(gca,{'widescreen'});
utils.exportPlot(fig9,[prefix '9']);

P_music = array.music(x,v,2,pi/2,1001);
P_music = P_music./max(abs(P_music(:)));
set(gca,'LineStyleOrderIndex',3);
plot(sin(psi_vec),10*log10(abs(P_music)),'DisplayName','MUSIC');
utils.setPlotStyle(gca,{'widescreen'});
utils.exportPlot(fig9,[prefix '11']);

%% Figure 10 - Beamscan Example Images

fig10=ex8_1;
utils.exportPlot(fig10,[prefix '10']);