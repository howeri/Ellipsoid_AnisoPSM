%% N sphere optimization 
% varying sphere size and locations
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
 
%% Optimization
p0 = repmat([0; 0; 0; 1], N, 1); % location_N, r_N
lb = repmat([-10; -10; -10; 1e-5], N, 1);
ub = repmat([10; 10; 10; 15], N, 1);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp')%,'MaxFunctionEvaluations', 10000000,'MaxIterations',1000000);
tic
errorFun = @(p) sqrt(mean((computeRCS(p)-rcsAbs).^2, 'all'));
[p, fval] = fmincon(errorFun, p0, [], [], [], [], lb, ub, [], options);
load('opt_3freq_16.mat')  %comment out when computing optimization
computationTime = toc;
disp(['Initial RMSE = ', num2str(errorFun(p0))])
disp(['Optimal RMSE = ', num2str(errorFun(p))])

%% Save Result
% save('opt_3freq','p') %comment out when not computing optimization

%% Plot RCS comparison
% 2D plot
plot2dRCScomparison_3freq(rcsAbs, computeRCS(p), N) % remember to change the legend

% 3D polar plot
plot3dRCScomparison_3freq(rcsAbs, computeRCS(p), N) % remember to change the legend

% 2D polar plot (fix phi= 0&-90)
plot2dRCSpolarComparison_3freq(rcsAbs, computeRCS(p), N) % remember to change the legend

% Plot sphere model
plotSphereModel(p, N)
 
% Plot sphere model (closer look)
plotSphereModel_closer(p, N)

%%
function rcsModel = computeRCS(p) % Input:N*4 variables, Output: RCS Table 
    global N direction
    locations = [p(1:4:((N-1)*4)+1) p(2:4:((N-1)*4)+2) p(3:4:((N-1)*4)+3)];  % N,3
    radius = p(4:4:((N-1)*4)+4);  % N,1   
    c = physconst('LightSpeed');
    fc = [1e9; 1.5e9; 2e9];      
    sphereReflection = sqrt(pi*radius.^2);
    delay = 2*locations*direction;  %location(N,3) * dir(3,18*36)
    rcsModel1 = pow2db(abs(sphereReflection' * exp(1i*2*pi*fc(1)/c*delay)).^2); %alpha'(1,N) * delay(N, 18*36)
    rcsModel2 = pow2db(abs(sphereReflection' * exp(1i*2*pi*fc(2)/c*delay)).^2);
    rcsModel3 = pow2db(abs(sphereReflection' * exp(1i*2*pi*fc(3)/c*delay)).^2); 
    rcsModel = [rcsModel1' rcsModel2' rcsModel3'];
end