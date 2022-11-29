%% N ellipsoid optimization 
% varying xyz axis lengthes, locations, and rotations
% freq = 1, 1.5, 2 GHz
% optimize in dB scale
clc;clear;close all
global N
N = 16;

%% Import aircraft RCS data (Cut out extra data)
global rcsAbs
rcs1 = importdata('f1.txt');
rcs2 = importdata('f1.5.txt');
rcs3 = importdata('f2.txt');
theta_data = rcs1.data(1:648,1);
phi_data = rcs1.data(1:648,2);
rcsAbs = [rcs1.data(1:648,6) rcs2.data(1:648,6) rcs3.data(1:648,6)];   % phi direction
rcsPhase = [rcs1.data(1:648,7) rcs2.data(1:648,7) rcs3.data(1:648,7)]; % phi direction

%% Construct a list of direction unit vectors with theta_data and phi_data
global direction
rho = 1;
r_xy = rho .* sind(theta_data);
x  = r_xy  .* cosd(phi_data);
y  = r_xy  .* sind(phi_data);
z  = rho .* cosd(theta_data);
direction = [x'; y'; z'];  % 3, 18*36

%% Construct the 2D matrix of theta/phi angles (theta first, then phi)
global THETA PHI
[THETA, PHI] = ndgrid((-180:10:170), (-90:10:80));

%% Optimization
p0 = repmat([0; 0; 0; 1; 1; 1; 0; 0], N, 1); % location_N, rxyz_ell_N, theta/phi_offset(rotation)_N
lb = repmat([-10; -10; -10; 1e-5; 1e-5; 1e-5; -45; -45], N, 1);
ub = repmat([10; 10; 10; 15; 15; 15; 45; 45], N, 1);
errorFun = @(p) sqrt(mean((rcsAbs - computeRCS(p)).^2, 'all'));
options = optimoptions('fmincon','Display','iter','Algorithm','sqp')%,'MaxFunctionEvaluations', 10000000,'MaxIterations',1000000);
tic
[p, fval] = fmincon(errorFun, p0, [], [], [], [], lb, ub, [], options);
load('opt_3freq.mat')  %comment out when computing optimization
computationTime = toc;
disp(['Initial RMSE = ', num2str(errorFun(p0))])
disp(['Optimal RMSE = ', num2str(errorFun(p))])

%% Save Result
% save('opt_3freq','p')  %comment out when not computing optimization

%% Plot RCS comparison
% 2D plot
plot2dRCScomparison_3freq(rcsAbs, computeRCS(p), N)

% 3D polar plot
plot3dRCScomparison_3freq(rcsAbs, computeRCS(p), N)

% 2D polar plot (fix phi= 0&-90)
plot2dRCSpolarComparison_3freq(rcsAbs, computeRCS(p), N)

% Plot ellipsoid model
plotEllipsoidModel(p, N)

% Plot ellipsoid model (closer)
plotEllipsoidModel_closer(p, N)

%%
function rcsModel = computeRCS(p) % Input:N*8 variables, Output: RCS table   
    global N THETA PHI direction
    locations = zeros(N, 3);
    rcs_models = zeros(36, 18, N);
    
    % RCS table for each ellipsoid 
    for k = 1:N
        locations(k, :) = [p(1 + (k-1)*8), p(2 + (k-1)*8), p(3 + (k-1)*8)];
        angle_offsets = [p(7 + (k-1)*8), p(8+ (k-1)*8)];       
        rx_ell = p(4 + (k-1)*8);
        ry_ell = p(5 + (k-1)*8);
        rz_ell = p(6 + (k-1)*8);
        phi_ell = PHI - angle_offsets(2);
        theta_ell = THETA - angle_offsets(1);
        sin_phi2 = sind(phi_ell).^2;
        cos_phi2 = cosd(phi_ell).^2;
        sin_theta2 = sind(theta_ell).^2;
        cos_theta2 = cosd(theta_ell).^2;
        rcs_models(:, :, k) = (pi * rx_ell^2 * ry_ell^2 * rz_ell^2) ./ (rx_ell^2*cos_phi2.*sin_theta2 + ry_ell^2*sin_phi2.*sin_theta2 + rz_ell^2.*cos_theta2).^2 ; 
    end
    
    c = physconst('LightSpeed');
    fc = [1e9; 1.5e9; 2e9];
    reflection = sqrt(rcs_models);
    reflection_reshaped = reshape(reflection, [36*18, 16]);
    delay = 2*locations*direction;  %location(N,3) * dir(3,18*36)
    rcsModel1 = pow2db(abs(sum(reflection_reshaped' .* exp(1i*2*pi*fc(1)/c*delay), 1)).^2);
    rcsModel2 = pow2db(abs(sum(reflection_reshaped' .* exp(1i*2*pi*fc(2)/c*delay), 1)).^2);
    rcsModel3 = pow2db(abs(sum(reflection_reshaped' .* exp(1i*2*pi*fc(3)/c*delay), 1)).^2);
    rcsModel =[rcsModel1' rcsModel2' rcsModel3'];
end