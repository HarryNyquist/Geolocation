% fig=book2_ex2_2()
%
% Executes Example 2.2 from Practical Geolocation for Electronic Warfare
% with MATLAB.
%
% INPUTS
%   none
%
% OUTPUTS
%   fig         figure handle
%
% Nicholas O'Donoughue
% 22 April 2021

% Set up sensor and target coordinates

%------------------------------------------------------------------------

% -----------------------------------------------------------------------

x_source = [3; 3]*1e3;

x_aoa = [4; 0]*1e3;
x_tdoa = [1, 3; 0, .5]*1e3;
x_fdoa = [0, 0; 1, 2]*1e3;
v_fdoa = [1, 1; -1, -1]*sqrt(.5)*300; % 300 m/s, at -45 deg heading
c = 299792458;% Speed of Light

% Error Covariance Matrix
err_aoa = 0.5; % 0.5 deg
cov_psi = (err_aoa*pi/180)^2; % rad^2

err_time = 20e-9; % 20 ns timing error
err_r = err_time * c; %utils.constants.c;
cov_r = (err_r)^2*eye(size(x_tdoa,2)); % m^2, double for the combination of test/ref msmts
cov_r_out = resampleCovMtx(cov_r,1);

freq_err = 5; % 5 Hz
f0 = 1e9; % Hz
rr_err = freq_err * c/f0; % (m/s)
cov_rr = rr_err^2*eye(size(x_fdoa,2)); % (m/s)^2
cov_rr_out = resampleCovMtx(cov_rr,1);

% Hybrid measurement and combined covariance matrix
z = hybrid_measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source);
cov_z = blkdiag(cov_psi, cov_r, cov_rr) % raw sensor measurements
cov_z_out = blkdiag(cov_psi, cov_r_out, cov_rr_out) % sensor pairs

% Generate Random Noise
L = chol(cov_z_out,'lower') % Cholesky decomposition of the covariance matrix
noise = L*randn(size(L,2),1)

% Noisy Measurements
zeta = z + noise;

%% ML Search Parameters
x_ctr = [2.5; 2.5]*1e3;
grid_size = [5e3; 5e3];
grid_res = 25;  % meters, grid resolution

%% GD and LS Search Parameters
x_init = [1; 1]*1e3;
epsilon = grid_res;
max_num_iterations = 100;
force_full_calc = true;
plot_progress = false;

%% ML Soln
x_ml = hybrid_mlSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, cov_z, ...
       x_ctr, grid_size, epsilon);

%% GD Soln
[x_gd, x_gd_full] = hybrid_gdSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, ...
                    cov_z, x_init, [], [], epsilon, ... 
                    max_num_iterations, force_full_calc, plot_progress);

%% LS Soln
[x_ls, x_ls_full] = hybrid_lsSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, zeta, ...
                    cov_z, x_init, epsilon, max_num_iterations, ...
                    force_full_calc, plot_progress);

%% Plot Result

fig=figure;
plot(x_source(1), x_source(2), 'kx', 'DisplayName','Target');
hold on;
plot(x_aoa(1), x_aoa(2), 'ko', 'DisplayName','AOA Sensor');
plot(x_tdoa(1, :), x_tdoa(2, :), 'ks', 'DisplayName','TDOA Sensor');
plot(x_fdoa(1, :), x_fdoa(2, :), 'k^', 'DisplayName','FDOA Sensor');
utils.drawArrow(x_fdoa(1,1)+[0 v_fdoa(1,1)],x_fdoa(2,1)+[0 v_fdoa(2,1)]);
utils.drawArrow(x_fdoa(1,2)+[0 v_fdoa(1,2)],x_fdoa(2,2)+[0 v_fdoa(2,2)]);

plot(x_ml(1), x_ml(2), 'v', 'DisplayName', 'ML Solution');
hdl=plot(x_gd_full(1,:), x_gd_full(2,:), '-.');
utils.excludeFromLegend(hdl);
plot(x_gd(1),x_gd(2),'-.+','DisplayName','GD Solution','Color',hdl.Color);
hdl=plot(x_ls_full(1,:), x_ls_full(2,:), '-');
utils.excludeFromLegend(hdl);
plot(x_ls(1), x_ls(2), '-*','DisplayName','LS Solution','Color',hdl.Color);

grid on;
ylim([0 4]*1e3);
xlim([-0.5 5.5]*1e3);
caxis([-20 0]);
set(gca,'ydir','normal');
legend('Location','NorthEast');
utils.setPlotStyle(gca,{'widescreen','tight'});

%------------------------------------------------------------------------
function cov_out = resampleCovMtx(cov, test_idx, ref_idx, test_wts, ref_wts)
% cov_out = resampleCovMtx(cov, test_idx, ref_idx, test_wts, ref_wts)
%
% Resample a covariance matrix based on a set of reference and test 
% indices.  This assumes a linear combination of the test and reference 
% vectors.  The output is an n_pair x n_pair covariance matrix for the 
% n_pair linear combinations.
%
% The measurements can be optionally weighted, to over-emphasize some
% measurements and de-emphasize others.
%
% In the resampled covariance matrix, the i,j-th entry is given
%    [Cout]_ij = [C]_bibj + [C]_aiaj - [C]_aibj - [C]_biaj
%       where:  ai, aj are the i-th and j-th reference indices
%               bi, bj are the i-th and j-th test indices
%               C is the input covariance matrix
%
% Any indices that are NaN will be ignored, to represent a single-sensor
% measurement, such as AoA (which does not need a reference sensor
% measurement against which to compare), or a noise-free measurement with
% no error.
%
% If the third input, ref_idx, is missing or empty, then the second input,
% test_idx, will be passed to utils.parseReferenceSensor to generate
% matching test and reference vectors.
%
% INPUTS:
%   cov         NxN covariance matrix of individual sensor measurements
%   test_idx    n_pair x 1 vector of test sensor indices
%   ref_idx     n_pair x 1 vector of reference sensor indices [Optional]
%   test_wts    Optional n_pair x 1 vector of test measurement weights
%   ref_wts     Optional n_pair x 1 vector of reference measurement weights
%
% OUTPUTS:
%   cov_out     n_pair x n_pair output covariance matrix of sensor
%               measurement pairs
%
% Nicholas O'Donoughue
% 25 May 2021

%% Input handling
% Parse array sizes and indices
n_sensor = size(cov, 1);

% Handle test/reference inputs
if nargin < 3 || isempty(ref_idx)
    [test_idx, ref_idx] = parseReferenceSensor(test_idx, n_sensor);
end

% Parse output size
n_test = numel(test_idx);
n_ref = numel(ref_idx);
n_out = max(n_test, n_ref);

if n_test > 1 && n_ref > 1 && n_test ~= n_ref
    error(strcat("Error calling covariance matrix resample. Reference and", ...
                 " test vectors must have the same shape."))
end

if any(test_idx > n_sensor) || any(ref_idx > n_sensor)
	error(strcat("Error calling covariance matrix resample. Indices exceed", ...
                 " the dimensions of the covariance matrix."))
end

% Parse sensor weights
do_test_wt = ~(nargin < 4 || isempty(test_wts));
do_ref_wt = ~(nargin < 5 || isempty(ref_wts));

if do_test_wt
    n_test_wt = numel(test_wts);
end

if do_ref_wt
    n_ref_wt = numel(ref_wts);
end

% Initialize output
cov_out = zeros(n_out, n_out);

a_i_wt = 1;
a_j_wt = 1;
b_i_wt = 1;
b_j_wt = 1;

% Step through reference sensors
for idx_row = 1:n_out
    % Parse sensor indices.  The mod commands seamlessly handle scalar
    % inputs
    a_i = test_idx(1+mod(idx_row-1, n_test));
    b_i = ref_idx(1+mod(idx_row-1, n_ref));

    % Parse sensor weights
    if do_test_wt
        a_i_wt = test_wts(1+mod(idx_row-1, n_test_wt));
    end

    if do_ref_wt
        b_i_wt = ref_wts(1+mod(idx_row-1, n_ref_wt));
    end

    for idx_col =1:n_out
        % Parse sensor indices.  The mod commands seamlessly handle scalar
        % inputs.
        a_j = test_idx(1+mod(idx_col-1, n_test));
        b_j = ref_idx(1+mod(idx_col-1, n_ref));

        if do_test_wt
            a_j_wt = test_wts(1+mod(idx_col-1, n_test_wt));
        end
    
        if do_ref_wt
            b_j_wt = ref_wts(1+mod(idx_col-1, n_ref_wt));
        end
        
        % Parse Input covariances
        if isnan(b_i) || isnan(b_j)
            cov_bibj = 0;
        else
            cov_bibj = cov(b_i, b_j);
        end
        if isnan(a_i) || isnan(a_j)
            cov_aiaj = 0;
        else
            cov_aiaj = cov(a_i, a_j);
        end
        if isnan(a_i) || isnan(b_j)
            cov_aibj = 0;
        else
            cov_aibj = cov(a_i, b_j);
        end
        if isnan(b_i) || isnan(a_j)
            cov_biaj = 0;
        else
            cov_biaj = cov(b_i, a_j);
        end
        
        %  [Cout]_ij = [C]_bibj + [C]_aiaj - [C]_aibj - [C]_biaj
        cov_out(idx_row, idx_col) = b_i_wt * b_j_wt * cov_bibj + ...
                                    a_i_wt * a_j_wt * cov_aiaj - ...
                                    a_i_wt * b_j_wt * cov_aibj - ...
                                    b_i_wt * a_j_wt * cov_biaj;

    end
end
end

%--------------------------------------------------------------------------
function [test_idx_vec, ref_idx_vec] = parseReferenceSensor(ref_idx, num_sensors)

%% Default Behavior
if isempty(ref_idx)
    % Default behavior is to generate all possible sensor pairs
    test_idx_vec = 1:num_sensors-1;
    ref_idx_vec = num_sensors * ones(size(test_idx_vec));
elseif strcmpi(ref_idx,'full')
    % Do the full set of N(N-1)/2 pairs
    full_set = nchoosek(1:num_sensors, 2);
    test_idx_vec = full_set(:,2)';
    ref_idx_vec = full_set(:,1)';
elseif isscalar(ref_idx)
    % Scalar reference provided, use all others as test sensors
    test_idx_vec = setdiff(1:num_sensors, ref_idx);
    ref_idx_vec = ref_idx * ones(size(test_idx_vec));
else
    % Explicit sensor pairs provided, parse them
    test_idx_vec = ref_idx(1, :);
    ref_idx_vec = ref_idx(2, :);
end
        
if size(test_idx_vec) ~= size(ref_idx_vec)
    warning('utils/parseReferenceSensor.m generated unequal test and reference vectors.  Check for bugs.');
end

return
end
%-------------------------------------------------------------------------
function z = hybrid_measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source, tdoa_ref_idx, fdoa_ref_idx, do2Daoa, alpha_aoa, alpha_tdoa, alpha_fdoa)
% z = measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source, tdoa_ref_idx, 
%                                                fdoa_ref_idx, do2Daoa)
%
% Computes hybrid measurements, for AOA, TDOA, and FDOA sensors.
%
% INPUTS:
%   x_aoa       nDim x nAOA array of sensor positions
%   x_tdoa      nDim x nTDOA array of TDOA sensor positions
%   x_fdoa      nDim x nFDOA array of FDOA sensor positions
%   v_fdoa      nDim x nFDOA array of FDOA sensor velocities
%   x_source    nDim x nSource array of source positions
%   tdoa_ref_idx    Index for reference TDOA sensor or 2 x nPair set of
%                   TDOA sensor pairing indices [optional]
%   fdoa_ref_idx    Index for reference FDOA sensor or 2 x nPair set of
%                   FDOA sensor pairing indices [optional]
%   do2Daoa     Boolean flag, if true then 2D AOA measurements will be
%               generated (azimuth and elevation)
%   alpha_aoa   (Optional) nAOA x 1 (or nAOA x 2 if 2DAOR) vector of AOA 
%                   bias terms [default=0]
%   alpha_tdoa  (Optional) nTDOA x 1 vector of range bias terms 
%                   (time bias * c) [default=0]
%   alpha_fdoa  (Optional) nFDOA x 1 vector of range-rate bias terms 
%                   (freq bias * c/f0) [default=0]
%
% OUTPUTS:
%   z           nAoa + nTDOA + nFDOA - 2 x nSource array of measurements
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 6 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end
if nargin < 7 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

if nargin < 8 || ~exist('do2Daoa','var')
    do2Daoa = true;
end

if nargin < 9 || ~exist('alpha_aoa','var')
    alpha_aoa = [];
end

if nargin < 10 || ~exist('alpha_tdoa','var')
    alpha_tdoa = [];
end

if nargin < 11 || ~exist('alpha_fdoa','var')
    alpha_fdoa = [];
end

% Construct component measurements
if ~isempty(x_aoa)
    z_a =measurement(x_aoa, x_source, do2Daoa, alpha_aoa);
else
    z_a = [];
end
if ~isempty(x_tdoa)
    z_t =measurement(x_tdoa, x_source,tdoa_ref_idx, alpha_tdoa);
else
z_t = [];
end
if ~isempty(x_fdoa)
    z_f =fdoa.measurement(x_fdoa, v_fdoa, x_source,fdoa_ref_idx, alpha_fdoa);
else
    z_f = [];
end

% Combine into a single data vector
z = [z_a; z_t; z_f];
end

%-------------------------------------------------------------------------
function psi = measurement(x_sensor, x_source, do2DAoA, alpha)
% psi = measurement(x_sensor, x_source, do2DAoA)
%
% Computes angle of arrival measurements.  If the input is 2D, then only
% the azimuth angle is report (atan(y/x)).  If the input is 3D, and the
% flag do2DAoA is missing or set to true, then the output contains two 
% sets of measurements, with azimuth reported first in psi(1:nSensor, :)
% and elevation reported second in psi(nSensor+1:end, :).
%
% INPUTS:
%   x_sensor    nDim x nSensor array of sensor positions
%   x_source    nDim x nSource array of source positions
%   do2DAoA     Boolean flag to activate 2D (az/el) AOA measurements
%               [default=True]
%   alpha       (Optional) nSensor x 1 vector of AOA biases (nSensor x 2 if
%               do2DAoA is true).  [default = 0]
%
% OUTPUTS:
%   psi         nSensor x nSource array of AOA measurements
%               (2*nSensor x nSource if 2D AOA measurements are generated)
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 3 || isempty(do2DAoA)
    do2DAoA = true;
end

if nargin < 4 || isempty(alpha)
    alpha_az = 0;
    alpha_el = 0;
else
    alpha_az = alpha(:,1);
    if do2DAoA && size(alpha,2) > 1
        alpha_el = alpha(:,2);
    else
        alpha_el = 0;
    end
end

[nDim1,nSensor] = size(x_sensor);
[nDim2,nSource] = size(x_source);

if nDim1~=nDim2
    error('First dimension of all inputs must match');
end
if nDim1 < 2
    error('Must have at least two dimensions (x and y)');
end
nDim = nDim1;

dx = reshape(x_source,nDim,1,nSource) - reshape(x_sensor,nDim,nSensor);
        % nDim x nSensor x nSource
        
az = reshape(atan2(dx(2,:,:),dx(1,:,:)),nSensor,nSource) + alpha_az;

if nDim >2 && do2DAoA
    ground_rng = sqrt(sum(abs(dx(1:2,:,:)).^2,1));
    el = reshape(atan2(dx(3,:,:),ground_rng),nSensor,nSource) + alpha_el;
    
    psi = cat(1, az, el);
else
    psi = az;
end
end

%--------------------------------------------------------------------------
function rdoa = tdoa_measurement(x_sensor, x_source, ref_idx, alpha)
% rdoa = measurement(x_sensor, x_source, ref_idx)
%
% Computes range difference measurements, using the
% final sensor as a common reference for all TDOA measurements.
%
% INPUTS:
%   x_sensor    nDim x nSensor array of sensor positions
%   x_source    nDim x nSource array of source positions
%   ref_idx     Either a scalar index for which sensor is the reference,
%               or a 2 x nPairing matrix of sensor pairing indices
%   alpha       nSensor x 1 vector of range bias terms
%
% OUTPUTS:
%   rdoa        nSensor -1 x nSource array of RDOA measurements
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
[nDim1,nSensor] = size(x_sensor);
[nDim2,nSource] = size(x_source);
if nDim1~=nDim2
    error('First dimension of all inputs must match');
end

if nargin < 3 || ~exist('ref_idx','var') || isempty(ref_idx)
    ref_idx = nSensor;
end

[test_idx_vec, ref_idx_vec] = parseReferenceSensor(ref_idx, nSensor);

if nargin < 4 || ~exist('alpha','var') || isempty(alpha)
    rdoa_bias = 0;
else
    % Parse the TDOA bias
    rdoa_bias = (alpha(test_idx_vec) - alpha(ref_idx_vec));
end

% Compute range from each source to each sensor
dx = reshape(x_source,nDim1,1,nSource) - reshape(x_sensor,nDim1,nSensor);
R = reshape(sqrt(sum(abs(dx).^2,1)),nSensor,nSource); % nSensor x nSource

% Compute range difference for each pair of sensors
rdoa = R(test_idx_vec,:) - R(ref_idx_vec,:) + rdoa_bias;
end

%-------------------------------------------------------------------------
function rrdoa = fdoa_measurement(x_sensor, v_sensor, x_source, ref_idx, alpha)
% rrdoa = measurement(x_sensor, v_sensor, x_source, ref_idx)
%
% Computed range rate difference measurements, using the
% final sensor as a common reference for all FDOA measurements.
%
% INPUTS:
%   x_sensor    nDim x nSensor array of sensor positions
%   v_sensor    nDim x nSensor array of sensor velocities
%   x_source    nDim x nSource array of source positions
%   ref_idx         Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings
%   alpha       nSensor x 1 array of range-rate bias terms
%
% OUTPUTS:
%   rrdoa       nSensor -1 x nSource array of RRDOA measurements [m/s]
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
[nDim1,nSensor1] = size(x_sensor);
[nDim2,nSensor2] = size(v_sensor);
[nDim3,nSource] = size(x_source);
if nDim1~=nDim2 || nSensor1 ~=nSensor2
    error('First two inputs must have macthing size');
end

if nDim1~=nDim3
    error('First dimension of all inputs must match');
end

if nargin < 3 || ~exist('ref_idx','var') || isempty(ref_idx)
    ref_idx = nSensor1;
end

% Parse Reference Sensor
n_sensor = size(x_sensor, 2);
[test_idx_vec, ref_idx_vec] = parseReferenceSensor(ref_idx, n_sensor);

% Parse FDOA Bias
if nargin < 4 || ~exist('alpha','var') || isempty(alpha)
    rrdoa_bias = 0;
else
    rrdoa_bias = alpha(test_idx_vec) - alpha(ref_idx_vec);
end

% Compute distance from each source position to each sensor
dx = reshape(x_source,nDim1,1,nSource) - reshape(x_sensor,nDim1,nSensor1);
R = sqrt(sum(abs(dx).^2,1)); % 1 x nSensor x nSource

% Compute range rate from range and velocity
rr = reshape(sum(v_sensor.*dx./R,1),nSensor1,nSource); % nSensor x nSource

% Apply reference sensors to compute range rate difference for each sensor
% pair
rrdoa = rr(test_idx_vec,:) - rr(ref_idx_vec,:) + rrdoa_bias;
end

%-------------------------------------------------------------------------
function [x_est,A,x_grid] = hybrid_mlSoln(x_aoa,x_tdoa,x_fdoa,v_fdoa,zeta,C,x_ctr,search_size,epsilon,tdoa_ref_idx,fdoa_ref_idx)
% [x_est,A,x_grid] = mlSoln(x_aoa,x_tdoa,x_fdoa,v_fdoa,zeta,C,x_ctr,...
%                           search_size,epsilon,tdoa_ref_idx,fdoa_ref_idx)
%
% Construct the ML Estimate by systematically evaluating the log
% likelihood function at a series of coordinates, and returning the index
% of the maximum.  Optionally returns the full set of evaluated
% coordinates, as well.
%
% INPUTS:
%   x_aoa           AOA sensor positions [m]
%   x_tdoa          TDOA sensor positions [m]
%   x_fdoa          FDOA sensor positions [m]
%   v_fdoa          FDOA sensor velocities [m/s]
%   zeta            Combined measurement vector
%   C               Combined measurement error covariance matrix
%   x_ctr           Center of search grid [m]
%   search_size     2-D vector of search grid sizes [m]
%   epsilon         Desired resolution of search grid [m]
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA
%
% OUTPUTS:
%   x_est           Estimated source position [m]
%   A               Likelihood computed across the entire set of candidate
%                   source positions
%   x_grid          Candidate source positions
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 10 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 11 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

% Set up function handle
ell = @(x) loglikelihood(x_aoa,x_tdoa,x_fdoa,v_fdoa,zeta,C,x,tdoa_ref_idx,fdoa_ref_idx);

% Call the util function
[x_est,A,x_grid] = utils_mlSoln(ell,x_ctr,search_size,epsilon);
end

%-------------------------------------------------------------------------
function ell = loglikelihood(x_aoa,x_tdoa,x_fdoa,v_fdoa,zeta,C,x_source,tdoa_ref_idx,fdoa_ref_idx)
% function ell = loglikelihood(x_aoa,x_tdoa,x_fdoa,zeta,C,x_source,...
%                                               tdoa_ref_idx,fdoa_ref_idx)
%
% Computes the Log Likelihood for Hybrid sensor measurement (AOA, TDOA, and
% FDOA), given the received measurement vector zeta, covariance matrix C, 
% and set of candidate source positions x_source.
%
% INPUTS:
%   x_aoa       nDim x nAOA vector of AOA sensor positions
%   x_tdoa      nDim x nTDOA vector of TDOA sensor positions
%   x_fdoa      nDim x nFDOA vector of FDOA sensor positions
%   v_fdoa      nDim x nFDOA vector of FDOA sensor velocities
%   zeta        Combined AOA/TDOA/FDOA measurement vector
%   C           Combined AOA/TDOA/FDOA measurement covariance matrix
%   x_source    Candidate source positions
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA measurements
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA measurements
%
% OUTPUTS:
%   ell         Log-likelihood evaluated at each position x_source.
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 8 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 9 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

n_dim = size(x_source,1);
n_source_pos = size(x_source,2);
ell = zeros(n_source_pos,1);

% Resample the covariance matrix
n_aoa = size(x_aoa,2);
n_tdoa = size(x_tdoa,2);
n_fdoa = size(x_fdoa,2);

% Determine if AOA measurements are 1D (azimuth) or 2D (az/el)
if n_aoa > 0
    assert(size(C,1) == n_aoa + n_tdoa + n_fdoa || size(C,1) == 2*n_aoa + n_tdoa + n_fdoa,'Unable to determine if AOA measurements are 1D or 2D');
    do2DAoA = size(C,1) == 2*n_aoa + n_tdoa + n_fdoa;
    if do2DAoA, n_aoa = 2*n_aoa; end
end

% Parse the TDOA and FDOA reference indices together
[tdoa_test_idx_vec, tdoa_ref_idx_vec] = parseReferenceSensor(tdoa_ref_idx,n_tdoa);
[fdoa_test_idx_vec, fdoa_ref_idx_vec] = parseReferenceSensor(fdoa_ref_idx,n_fdoa);
test_idx_vec = cat(2,1:n_aoa, n_aoa + tdoa_test_idx_vec, n_aoa + n_tdoa + fdoa_test_idx_vec);
ref_idx_vec = cat(2,nan(1,n_aoa), n_aoa + tdoa_ref_idx_vec, n_aoa + n_tdoa + fdoa_ref_idx_vec);

% For now, we assume the AOA is independent of TDOA/FDOA
C_tilde = resampleCovMtx(C, test_idx_vec, ref_idx_vec);

% Ensure the covariance matrix is invertible
C_tilde = ensureInvertible(C_tilde);

% Pre-compute covariance matrix inverses
do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C_tilde,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C_tilde);
end

for idx_source = 1:n_source_pos
    x_i = x_source(:,idx_source);
    
    % Generate the ideal measurement matrix for this position
    z = hybrid_measurement(x_aoa,x_tdoa,x_fdoa,v_fdoa, x_i, tdoa_ref_idx, fdoa_ref_idx);
    
    % Evaluate the measurement error
    err = (z - zeta);

    % Compute the scaled log likelihood
    if do_decomp
        ell(idx_source) = -err'/C_d*err;
    else
        ell(idx_source) = -err'*C_inv*err;
    end
end
end

%-----------------------------------------------------------------------
function [x_est,A, x_grid] = utils_mlSoln(ell,x_ctr,search_size,epsilon)
% function [x_est,A,x_grid] = mlSoln(ell,x_ctr,search_size,epsilon)
%
% Execute ML estimation through brute force computational methods.
%
% INPUTS:
%   ell          Function handle for the likelihood of a given position
%                must accept x_ctr (and similar sized vectors) as the sole
%                input.
%   x_ctr        Center position for search space (x, x/y, or z/y/z).
%   search_size  Search space size (same units as x_ctr)
%   epsilon      Search space resolution (same units as x_ctr)
%
% OUTPUTS:
%   x_est        Estimated minimum
%   A            Likelihood computed at each x position in the search space
%   x_grid       Set of x positions for the entire search space (M x N) for
%                N=1, 2, or 3.
%
% Nicholas O'Donoughue
% 1 July 2019

[x_set, x_grid] = make_nd_grid(x_ctr, search_size, epsilon);

% rearrange to a matrix, where each column is
% Evaluate the likelihood function at each coordinate in the search space
A = ell(x_set);

% Find the peak
[~,idx_pk] = max(A(:));
x_est = x_set(:,idx_pk);
            

%--------------------------------------------------------------------------
function [x_set, dims] = make_nd_grid(ctr,max_offset,spacing)
% x_set = make_nd_grid(ctr,max_offset,spacing)
%
% Accepts a center value, maximum offset, and spacing vectors, of arbitrary
% length n_dim, and generates a search grid over n_dim dimensions.
%
% The returned matrix x_set has n_dim columns and N rows, where N is the
% product of the search vector in each dimension.
%    n_elements = 1 + 2*max_offset./spacing.
%    N = prod(n_elements)
%
% Inputs ctr, max_offset, and spacing must either be scalar, or have a
% common number of elements.  They are assumed to be vectors (shape is not
% preserved).
%
% Includes some error checking on the total array size allowable; which is
% currently set at no more than 10% of MATLAB's maximum array size, to be
% conservative.
%
% INPUTS:
%   ctr         Scalar or n_dim vector of center point for each dimension
%   max_offset  Scalar or n_dim vector of maximum search space offset from
%               ctr for each dimension
%   epsilon     Scalar or n_dim vector of search space resolution for each
%               dimension
%
% OUTPUTS:
%   x_set       n_dim x n_pt matrix of coordinates for n_pt test points
%   x_grid      n_dim x 1 cell array, each bearing the principal vector for
%               one dimension (to be used in plotting commands)
%
% Nicholas O'Donoughue
% 7 Nov 2021

%% Parse Inputs and Check for Max Array Size
n_dims = numel(ctr);

if numel(max_offset)==1
    max_offset = max_offset*ones(n_dims,1);
end

if numel(spacing)==1
    spacing = spacing * ones(n_dims,1);
end

assert(n_dims == numel(max_offset) && ...
       n_dims == numel(spacing),...
       'Search space dimensions do not match across specification of the center, search_size, and epsilon.');

n_elements = 1+ 2*max_offset(:)./spacing(:);

% Check Search Size
[user, ~] = memory;
maxElements = user.MaxPossibleArrayBytes/10; % Set a conservative limit,
                                             % since we'll need at least 2
assert(prod(n_elements) < maxElements, 'Search size is too large; MATLAB is likely to crash or become unresponsive.');

%% Initialize search space

% dims is a cell array of dimensions, each of which contains a vector of
% grid points along that dimension
dims = arrayfun(@(x,x_mx,n) x + x_mx * linspace(-1,1,n),  ctr(:), max_offset(:), n_elements(:), 'UniformOutput',false);

% Use meshgrid expansion; each element of dims_full is now a full n_dim 
% dimensioned grid for one of the elements of x
[x_grid{1:numel(dims)}] = ndgrid(dims{:}); %use comma-separated list expansion on both sides

% Rearrange to an n_dim x N matrix
x_grid_vec = cellfun(@(A) reshape(A,1,[]), x_grid, 'UniformOutput',false);
x_set = cell2mat(x_grid_vec(:)); % n_dim x N
end
end
%------------------------------------------------------------------------
function cov_out = ensureInvertible(cov, epsilon)
% function cov_out = ensureInvertible(cov, epsilon)
%
% Check the input matrix for invertibility by finding the eigenvalues and
% checking that they are all >= a small value (epsilon).
%
% If any of the eigenvalues are too small, then a diagonal loading term
% is applied to ensure that the matrix is positive definite (all
% eigenvalues are >= epsilon).
%
% INPUTS:
%   cov     Input square covariance matrix.  If the input has >2
%           dimensions, then the process is repeated across the extra
%           dimensions.
%
%   epsilon (Optional) Specifies the minimum eigenvalue to use when
%           checking for invertibility. Default = 1e-10
%
% OUTPUTS:
%   cov_out Output covariance matrix, guaranteed to be invertible.
%
% Nicholas O'Donoughue
% 1 June 2021

% Check for epsilon input
if nargin < 2 || isempty(epsilon)
    epsilon = 1e-20;
end

% Check input dimensions
sz = size(cov);
assert(numel(sz) > 1, 'Input must have at least two dimensions.');
assert(sz(1) == sz(2), 'First two dimensions of input matrix must be equal.');
dim = sz(1);
if numel(sz) > 2
    n_matrices = prod(sz(3:end));
else
    n_matrices = 1;
end

cov_out = zeros(size(cov));
for idx_matrix = 1:n_matrices
   % Check min eigenvalue
   if min(eig(cov(:,:,idx_matrix))) < epsilon
      % We need to add a diagonal loading term, determine appropriate size
      d = epsilon;
      while min(eig(cov(:,:,idx_matrix) + d * eye(dim))) < epsilon
          d = d * 10;
      end
      
      % Apply diagonal loading
      cov_out(:,:,idx_matrix) = cov(:,:,idx_matrix) + d * eye(dim);
   else
       % No diagonal loading
       cov_out(:,:,idx_matrix) = cov(:,:,idx_matrix);
   end
end
end

%------------------------------------------------------------------------
function [x,x_full] = hybrid_gdSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, z,C,x_init,alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress,tdoa_ref_idx,fdoa_ref_idx)
% [x,x_full] = gdSoln(x_aoa, x_tdoa, x_fdoa, v_fdoa, z,C,x_init,alpha,...
%            beta,epsilon,max_num_iterations,force_full_calc,plot_progress)
%
% Computes the gradient descent solution for hybrid AOA, TDOA, and
% FDOA processing.
%
% Inputs:   
%   x_aoa               AOA sensor positions
%   x_tdoa              TDOA sensor positions
%   x_fdoa              FDOA sensor positions
%   v_fdoa              FDOA sensor velocities
%   z                   Measurement vector
%   C                   Combined error covariance matrix
%   x_init              Initial estimate of source position [m]
%   alpha               Backtracking line search parameter
%   beta                Backtracking line search parameter
%   epsilon             Desired position error tolerance (stopping 
%                       condition)
%   max_num_iterations  Maximum number of iterations to perform
%   force_full_calc     Boolean flag to force all iterations (up to
%                       max_num_iterations) to be computed, regardless
%                       of convergence (DEFAULT = False)
%   plot_progress       Boolean flag dictacting whether to plot
%                       intermediate solutions as they are derived 
%                       (DEFAULT = False).
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA
% Outputs:
%   x               Estimated source position
%   x_full          Iteration-by-iteration estimated source positions
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 15 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

if nargin < 14 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 13 || ~exist('plot_progress','var')
    plot_progress = false;
end

if nargin < 12 || ~exist('force_full_calc','var')
    force_full_calc = false;
end

if nargin < 11 || ~exist('max_num_iterations','var')
    max_num_iterations = [];
end

if nargin < 10 || ~exist('epsilon','var')
    epsilon = [];
end

if nargin < 9 || ~exist('beta','var')
    beta = [];
end

if nargin < 8 || ~exist('alpha','var')
    alpha = [];
end

% Initialize measurement error and Jacobian function handles
y = @(x) z - hybrid_measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x, tdoa_ref_idx, fdoa_ref_idx);
J = @(x) hybrid.jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x, tdoa_ref_idx, fdoa_ref_idx);

% Resample the covariance matrix
n_aoa = size(x_aoa,2);
n_tdoa = size(x_tdoa,2);
n_fdoa = size(x_fdoa,2);

% Parse the TDOA and FDOA reference indices together
[tdoa_test_idx_vec, tdoa_ref_idx_vec] = parseReferenceSensor(tdoa_ref_idx,n_tdoa);
[fdoa_test_idx_vec, fdoa_ref_idx_vec] = parseReferenceSensor(fdoa_ref_idx,n_fdoa);
test_idx_vec = cat(2,tdoa_test_idx_vec, n_tdoa + fdoa_test_idx_vec);
ref_idx_vec = cat(2,tdoa_ref_idx_vec, n_tdoa + fdoa_ref_idx_vec);

% For now, we assume the AOA is independent of TDOA/FDOA
C_aoa = C(1:n_aoa, 1:n_aoa);
C_tfdoa = C(n_aoa+1:end, n_aoa+1:end);
C_tilde = blkdiag(C_aoa, resampleCovMtx(C_tfdoa, test_idx_vec, ref_idx_vec));

[x,x_full] = utils_gdSoln(y,J,C_tilde,x_init,alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress);
end
%------------------------------------------------------------------------
function [x,x_full] = utils_gdSoln(y,J,C,x_init,alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress)
% [x,x_full] = gdSoln(y,J,C,x_init,alpha,beta,epsilon,max_num_iterations,force_full_calc,plot_progress)
%
% Computes the gradient descent solution for localization given the 
% provided measurement and Jacobian function handles, and measurement 
% error covariance.
%
% Inputs:
%   
%   y               Measurement vector function handle (accepts n_dim 
%                   vector of source position estimate, responds with error 
%                   between received and modeled data vector)
%   J               Jacobian matrix function handle (accepts n_dim vector
%                   of source position estimate, and responds with n_dim x
%                   n_sensor Jacobian matrix
%   C               Measurement error covariance matrix
%   x_init          Initial estimate of source position
%   alpha           Backtracking line search parameter
%   beta            Backtracking line search parameter
%   epsilon         Desired position error tolerance (stopping condition)
%   max_num_iterations  Maximum number of LS iterations to perform
%   force_full_calc Forces all max_num_iterations to be executed
%   plot_progress   Binary flag indicating whether to plot error/pos est
%                   over time
%
% Outputs:
%   x               Estimated source position
%   x_full          Iteration-by-iteration estimated source positions
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
n_dims = numel(x_init);

if nargin < 9 || isempty(plot_progress)
    % By default, do not plot
    plot_progress = false;
end

if nargin < 8 || isempty(force_full_calc)
    % If not specified, don't force the solver to run until it hits
    % max_num_iterations
    force_full_calc = false;
end

if nargin < 7 || isempty(max_num_iterations)
    % Number of iterations not specified, default is 10,000
    max_num_iterations = 10000;
end

if nargin < 6 || isempty(epsilon)
    % Maximum error not specified; default to .1 (distance units)
    epsilon = 1e-6;
end

if nargin < 5 || isempty(beta)
    beta = .8;
end

if nargin < 4 || isempty(alpha)
    alpha = .3;
end

% Initialize loop
iter = 1;
error = Inf;
x_full = zeros(n_dims,max_num_iterations);
x_prev = x_init;
x_full(:,1) = x_prev;

% Ensure the covariance matrix is invertible
if exist('+utils/ensureInvertible.m','file')
    C = ensureInvertible(C);
end

% Pre-compute covariance matrix inverses
do_decomp = exist('OCTAVE_VERSION','builtin')==0 && ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C);
end

% Cost Function for Gradient Descent
if do_decomp
    f = @(x) y(x)'/C_d*y(x);
else
    f = @(x) y(x)'*C_inv*y(x);
end

% Initialize Plotting
if plot_progress
    fig_plot=figure;
    subplot(211);
    xlabel('Iteration Number');
    ylabel('Change in Position Estimate');
    hold on;
    set(gca,'yscale','log')
    subplot(212);
    xlabel('x');
    ylabel('y');
end

% Divergence Detection
num_expanding_iters = 0;
max_num_expanding_iters = 5;
prev_error = Inf;

% Loop until either the desired tolerance is achieved or the maximum
% number of iterations have occurred
while iter < max_num_iterations && (force_full_calc || error >= epsilon)
    iter = iter+1;

    % Evaluate Residual and Jacobian Matrix
    y_i = y(x_prev);
    J_i = J(x_prev);
    
    % Compute Gradient and Cost function
    if do_decomp
        grad = -2*J_i/C_d*y_i;
    else
        grad = -2*J_i*C_inv*y_i;
    end
    
    % Descent direction is the negative of the gradient
    del_x = -grad/norm(grad);
    
    % Compute the step size
    t = utils_backtrackingLineSearch(f,x_prev,grad,del_x,alpha,beta);
    
    % Update x position
    x_full(:,iter) = x_prev + t*del_x;
    
    % Update variables
    x_prev = x_full(:,iter);
    error = t;
    
    if plot_progress
        figure(fig_plot);
        subplot(211);
        plot(iter,error,'.');
        subplot(212);
        plot(x_full(1,1:iter),x_full(2,1:iter),'-+');
    end
    
    % Check for divergence
    if error <= prev_error
        num_expanding_iters = 0;
    else
        num_expanding_iters = num_expanding_iters + 1;
        if num_expanding_iters >= max_num_expanding_iters
            % Divergence detected
            x_full(:,iter:end) = NaN;%repmat(x_full(:,iter),1,maxIters-iter+1);
            break;
        end
    end
    prev_error = error;
end

% Bookkeeping
if ~force_full_calc
    x_full = x_full(:,1:iter);
end
x = x_full(:,iter);
end

%------------------------------------------------------------------------
function t = utils_backtrackingLineSearch(f,x,grad,del_x,alpha,beta)
% t = backtrackingLineSearch(f,x,grad,del_x,alpha,beta)
%
% Performs backtracking line search according to algorithm 9.2 of
% Stephen Boyd's, Convex Optimization
%
% Inputs:
%
%   f       Function handle to minimize
%   x       Current estimate of x
%   grad    Gradient of f at x
%   del_x   Descent direction
%   alpha   Constant between 0 and 0.5
%   beta    Constant between 0 and 1
%
% Outputs:
%
%   t       Optimal step size for the current estimate x.
%
% Nicholas O'Donoughue
% 1 July 2019

% Initialize the search parameters and direction
t = 1;
startingVal = f(x);
slope = grad'*del_x;

% Make sure the starting value is large enough
while f(x+t*del_x) < startingVal+alpha*t*slope
    t = 2*t;
end

% Conduct the backtracking line search
while f(x+t*del_x) > startingVal+alpha*t*slope
    t = beta*t;
end
end

%-------------------------------------------------------------------------
function [x,x_full] = hybrid_lsSoln(x_aoa,x_tdoa,x_fdoa,v_fdoa,z,C,x_init,epsilon,max_num_iterations, force_full_calc, plot_progress, tdoa_ref_idx, fdoa_ref_idx)
% [x,x_full] = lsSoln(x0,rho,C,xs_init,epsilon,numIterations,
%                                      force_full_calc, plot_progress)
%
% Computes the least square solution for combined AOA, TDOA, and
% FDOA processing.
%
% Inputs:
%   x_aoa               AOA sensor positions
%   x_tdoa              TDOA sensor positions
%   x_fdoa              FDOA sensor positions
%   v_fdoa              FDOA sensor velocities
%   z                   Measurement vector
%   C                   Combined Error Covariance Matrix
%   x_init             Initial estimate of source position [m]
%   epsilon             Desired estimate resolution [m]
%   max_num_iterations  Maximum number of iterations to perform
%   force_full_calc     Boolean flag to force all iterations (up to
%                       max_num_iterations) to be computed, regardless
%                       of convergence (DEFAULT = False)
%   plot_progress       Boolean flag dictacting whether to plot
%                       intermediate solutions as they are derived 
%                       (DEFAULT = False).
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA
%
% Outputs:
%   x               Estimated source position
%   x_full          Iteration-by-iteration estimated source positions
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
if nargin < 13 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

if nargin < 12 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 11 || ~exist('plot_progress','var')
    plot_progress = false;
end

if nargin < 10 || ~exist('force_full_calc','var')
    force_full_calc = false;
end

if nargin < 9 || ~exist('max_num_iterations','var')
    max_num_iterations = [];
end

if nargin < 8 || ~exist('epsilon','var')
    epsilon = [];
end

% Initialize measurement error and Jacobian function handles
y = @(x) z - hybrid_measurement(x_aoa, x_tdoa, x_fdoa, v_fdoa, x, tdoa_ref_idx, fdoa_ref_idx);
J = @(x) hybrid_jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x, tdoa_ref_idx, fdoa_ref_idx);

% Resample the covariance matrix
n_aoa = size(x_aoa,2);
n_tdoa = size(x_tdoa,2);
n_fdoa = size(x_fdoa,2);

% Determine if AOA measurements are 1D (azimuth) or 2D (az/el)
assert(size(C,1) == n_aoa + n_tdoa + n_fdoa || size(C,1) == 2*n_aoa + n_tdoa + n_fdoa,'Unable to determine if AOA measurements are 1D or 2D');
do2DAoA = size(C,1) == 2*n_aoa + n_tdoa + n_fdoa;
if do2DAoA, n_aoa = 2*n_aoa; end

% Parse the TDOA and FDOA reference indices together
[tdoa_test_idx_vec, tdoa_ref_idx_vec] = parseReferenceSensor(tdoa_ref_idx,n_tdoa);
[fdoa_test_idx_vec, fdoa_ref_idx_vec] = parseReferenceSensor(fdoa_ref_idx,n_fdoa);
test_idx_vec = cat(2,tdoa_test_idx_vec, n_tdoa + fdoa_test_idx_vec);
ref_idx_vec = cat(2,tdoa_ref_idx_vec, n_tdoa + fdoa_ref_idx_vec);

% For now, we assume the AOA is independent of TDOA/FDOA
C_aoa = C(1:n_aoa, 1:n_aoa);
C_tfdoa = C(n_aoa+1:end, n_aoa+1:end);
C_tilde = blkdiag(C_aoa, resampleCovMtx(C_tfdoa, test_idx_vec, ref_idx_vec));

% Call the generic Least Square solver
[x,x_full] = utils_lsSoln(y,J,C_tilde,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress);

end
%------------------------------------------------------------------------

function [J,Jv] = jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,fdoa_ref_idx, v_source)
% [J,Jv] = jacobian(x_aoa, x_tdoa, x_fdoa, v_fdoa, x_source,tdoa_ref_idx,...
%                                                  fdoa_ref_idx, v_source)
%
% Returns the Jacobian matrix for hybrid set of AOA, TDOA, and FDOA
% measurements.
%
% If the target is moving, as specified by an optional fifth input
% v_source, then the Jacobian is provided with respect to both the target
% position and velocity.  This is only necessary if the geolocation
% algorithm is also solving for target velocity.  If target velocity is
% assumed known, or is not being estimated, then the source velocity can be
% subtracted from sensor velocity.
%
% INPUTS:
%   x_aoa           AOA sensor positions
%   x_tdoa          TDOA sensor positions
%   x_fdoa          FDOA sensor positions
%   v_fdoa          FDOA sensor velocities
%   x_source        Candidate source positions
%   tdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for TDOA measurements
%   fdoa_ref_idx    Scalar index of reference sensor, or nDim x nPair
%                   matrix of sensor pairings for FDOA measurements
%   v_source        [Optional] nDim x nSource vector of source velocities
%                   Target assumed stationary if not provided.
%
% OUTPUTS:
%   J               nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position
%   Jv              nDim x nMeasurement x nSource matrix of Jacobians,
%                   one for each candidate source position (if v_source is
%                   provided).
%
% Nicholas O'Donoughue
% 1 July 2019

if nargin < 6 || ~exist('tdoa_ref_idx','var')
    tdoa_ref_idx = [];
end

if nargin < 7 || ~exist('fdoa_ref_idx','var')
    fdoa_ref_idx = [];
end

% Compute Jacobian for AOA measurements
if ~isempty(x_aoa)
    J_aoa = triang.jacobian(x_aoa, x_source);
else
    J_aoa = [];
end

% Compute Jacobian for TDOA measurements
if ~isempty(x_tdoa)
    J_tdoa= tdoa.jacobian(x_tdoa, x_source, tdoa_ref_idx);
else
    J_tdoa = [];
end

% Compute Jacobian for FDOA measurements
if ~isempty(x_fdoa) && ~isempty(v_fdoa)
    do_vel_jacobian = nargin > 4 && exist('v_source','var') && ~isempty(v_source);

    if do_vel_jacobian
        [J_fdoa, J_fdoa_v] = fdoa.jacobian(x_fdoa, v_fdoa, x_source, fdoa_ref_idx, v_source);
    else
        J_fdoa = fdoa.jacobian(x_fdoa, v_fdoa, x_source, fdoa_ref_idx);
        J_fdoa_v = [];
    end
else
    do_vel_jacobian = false;

    J_fdoa = [];
    J_fdoa_v = [];
end

%% Combine component Jacobians
J = [J_aoa, J_tdoa, J_fdoa];

if do_vel_jacobian
    num_rows = size(J_fdoa_v,1);
    num_zero_cols = size(J_aoa,2) + size(J_tdoa,2);
    Jv = cat(2,zeros(num_rows, num_zero_cols), J_fdoa_v);
else
    Jv = [];
end
end

%-------------------------------------------------------------------------
function [x,x_full] = utils_lsSoln(y,J,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress)
% [x,x_full] = lsSoln(y,J,C,x_init,epsilon,max_num_iterations,force_full_calc,plot_progress)
%
% Computes the least square solution for TDOA processing.
%
% Inputs:
%   
%   y               Measurement vector function handle (accepts n_dim 
%                   vector of source position estimate, responds with error 
%                   between received and modeled data vector)
%   J               Jacobian matrix function handle (accepts n_dim vector
%                   of source position estimate, and responds with n_dim x
%                   n_sensor Jacobian matrix
%   C               Measurement error covariance matrix
%   x_init          Initial estimate of source position
%   epsilon         Desired position error tolerance (stopping condition)
%   max_num_iterations  Maximum number of LS iterations to perform
%   force_full_calc Forces all max_num_iterations to be calculated
%   plot_progress   Binary flag indicating whether to plot error/pos est
%                   over time
%
% Outputs:
%   x               Estimated source position
%   x_full          Iteration-by-iteration estimated source positions
%
% Nicholas O'Donoughue
% 1 July 2019

% Parse inputs
n_dims = numel(x_init);

if nargin < 8 || isempty(plot_progress)
    % By default, do not plot
    plot_progress = false;
end

if nargin < 7 || isempty(force_full_calc)
    % If not specified, don't force the solver to run until it hits
    % max_num_iterations
    force_full_calc = false;
end

if nargin < 6 || isempty(max_num_iterations)
    % Number of iterations not specified, default is 10,000
    max_num_iterations = 10000;
end

if nargin < 5 || isempty(epsilon)
    % Maximum error not specified; default to .1 (distance units)
    epsilon = 1e-6;
end

% Initialize loop
iter = 1;
error = Inf;
x_full = zeros(n_dims,max_num_iterations);
x_prev = x_init;
x_full(:,1) = x_prev;
    
% For now, we assume the AOA is independent of TDOA/FDOA
C = ensureInvertible(C);

% Pre-compute covariance matrix inverses
do_decomp = ~verLessThan('MATLAB','9.3');
if do_decomp
    % Starging in R2017b, MATLAB released the DECOMPOSITION function,
    % which can decompose matrices for faster computation of left- and
    % right-division in for loops.
    C_d = decomposition(C,'chol');
else
    % If DECOMPOSITION is unavailable, let's precompute the pseudo-inverse.
    C_inv = pinv(C);
end

% Initialize Plotting
if plot_progress
    figure;
    xlabel('Iteration Number');
    ylabel('Change in Position Estimate');
    hold on;
    set(gca,'yscale','log');
end

% Divergence Detection
num_expanding_iters = 0;
max_num_expanding_iters = 10;
prev_error = Inf;

% Loop until either the desired tolerance is achieved or the maximum
% number of iterations have occurred
while iter < max_num_iterations && (force_full_calc || error >= epsilon)
    iter = iter+1;

    % Evaluate Residual and Jacobian Matrix
    y_i = y(x_prev);
    J_i = J(x_prev);
            
    % Compute delta_x^(i), according to 13.18
    if do_decomp
        % Check for invertibility
        jc = J_i/C_d;
        jcj = jc*J_i';
        if cond(jcj) > 10000
            % Ill-conditioned, apply diagonal loading
            diag_ldng = 1e-10*eye(size(J_i,1));
            jcj = jcj + diag_ldng;
        end
        delta_x = jcj\jc*y_i;
    else
        jc = J_i*C_inv;
        jcj = jc*J_i';
        delta_x = jcj\jc*y_i;
    end
    
    % Update predicted location
    x_full(:,iter) = x_prev + delta_x;
    
    % Update variables
    x_prev = x_full(:,iter);
    error = norm(delta_x);

    if plot_progress
        plot(iter,error,'.');
    end
    
    % Check for divergence
    if error <= prev_error
        num_expanding_iters = 0;
    else
        num_expanding_iters = num_expanding_iters + 1;
        if num_expanding_iters >= max_num_expanding_iters
            % Divergence detected
            break;
        end
    end
    prev_error = error;
end

% Bookkeeping
if ~force_full_calc
    x_full = x_full(:,1:iter);
end
x = x_full(:,iter);
end

%-------------------------------------------------------------------------


%---------------------------------------------------------------------------