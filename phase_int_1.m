% % Code for AOA
% x_source = [3; 200]*1e3;
% x_aoa = [0, 1; 0, 0]*1e3;
% err_delta_phi = 1*pi/180; % error in phase estimation 1deg
% f = 162.025e6;   c = 3*1e8;
% lam = c/f ;
% d = 1e3  ;
% d_lam = d/lam 
% true_AOA = atan2((x_source(2,1)-x_aoa(2,1)),(x_source(1,1)-x_aoa(1,1)))  % in radians
% phase_diff = 2*pi*d_lam*sin(pi/2 - true_AOA) 
% wrapped_phase_diff = mod(phase_diff + pi, 2*pi) - pi
% psi_coarse = (pi/180)*89.18 ;
% psi_est = [];
% k = floor(d_lam);
% n= k
% phi_diff_est = [];
% noise_delta_phi = 1*pi/180*randn(1);
% if mod(k,2) == 0
%     for i = 1:k
%         %fprintf("st")
%         phi_diff_est(i) = wrapped_phase_diff + 2*pi*((i-1) -n/2) ;
%         phi_diff_est(i) = phi_diff_est(i) + noise_delta_phi ;
%         
%     end
% else
%     for i = 1:k
%         %fprintf("st")
%         phi_diff_est(i) = wrapped_phase_diff + 2*pi*((i-1) -(n-1)/2) ;
%         phi_diff_est(i) = wrapped_phase_diff + noise_delta_phi;
%         
%     end
% end
% 
% 
% for i = 1:k
%     psi_est(i) = pi/2 - asin(phi_diff_est(i)/(2*pi*d_lam));
% end
% 
% resolved_psi_est = [];
% for i = 1:k
%     if abs(psi_coarse -psi_est(i))< 0.1*pi/180 ;
%       resolved_psi_est = [resolved_psi_est, psi_est(i)*180/pi]; % Append value correctly
%     end
% end
% 
% resolved_psi_est
    
% numSimulations = 20;
% for i = 1:numSimulations
%     noise_delta_phi = 1*pi/180*randn(1);
%     
%     phi_est(i) = wrapped_phase_diff + noise_delta_phi;
%     psi_est(i) = 90 - (180/pi)*asin(phi_est(i)/(2*pi*d_lam));
% end
% 
% psi_est = mean(psi_est);

%Constants
num_trials = 1;  % Number of Monte Carlo iterations
x_source = [3; 200] * 1e3;
x_aoa = [0, 1; 0, 0] * 1e3;
err_delta_phi = 1 * pi / 180; % Error in phase estimation (1 degree)
f = 162.025e6;   
c = 3e8;  
lam = c / f;
d = 1e3;  
d_lam = d / lam;
true_AOA = atan2((x_source(2,1) - x_aoa(2,1)), (x_source(1,1) - x_aoa(1,1)));  % Radians
phase_diff = 2 * pi * d_lam * sin(pi/2 - true_AOA);
wrapped_phase_diff = mod(phase_diff + pi, 2*pi) - pi;
psi_coarse = (pi/180) * 89.18;
k = floor(d_lam);
n = k;

% Store resolved psi estimates
resolved_psi_est_all = zeros(num_trials, k);  

for trial = 1:num_trials
    phi_diff_est = zeros(1, k);  % Preallocate array
    noise_delta_phi = (1 * pi / 180) * randn(1);  % Gaussian noise (1° std)

    % Compute phase difference estimates
    if mod(k, 2) == 0
        for i = 1:k
            phi_diff_est(i) = wrapped_phase_diff + 2 * pi * ((i - 1) - n / 2);
            phi_diff_est(i) = phi_diff_est(i) + noise_delta_phi;
        end
    else
        for i = 1:k
            phi_diff_est(i) = wrapped_phase_diff + 2 * pi * ((i - 1) - (n - 1) / 2);
            phi_diff_est(i) = phi_diff_est(i) + noise_delta_phi;
        end
    end

    % Estimate psi values
    psi_est = pi/2 - asin(phi_diff_est / (2 * pi * d_lam));

    % Resolve psi estimates within threshold
    resolved_psi_est = [];
    for i = 1:k
        if abs(psi_coarse - psi_est(i)) < (0.06 * pi / 180)
            resolved_psi_est = [resolved_psi_est, psi_est(i) * 180 / pi];  % Convert to degrees
        end
    end
    resolved_psi_est
    % Store results for this trial
    if ~isempty(resolved_psi_est)
        resolved_psi_est_all(trial, 1:length(resolved_psi_est)) = resolved_psi_est;
    end
end

% % Display results
% disp('Resolved Psi Estimates (in degrees) from Monte Carlo simulations:');
% disp(resolved_psi_est_all);
