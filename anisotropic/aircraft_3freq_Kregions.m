%% N ellipsoids optimization 
% varying xyz axis lengthes, locations, and rotations
% angle increment is not paramatized yet
% freq = 1, 1.5, 2 GHz
% optimize in dB scale
% Predefine K regions and seperate the optimization problem
% clc;clear;close all
global N
global K
N = 16;
K = 4;

%% Import rcs data (Cut out extra data)
rcs1 = importdata('f1.txt');
rcs2 = importdata('f1.5.txt');
rcs3 = importdata('f2.txt');
theta_data = rcs1.data(:,1);
phi_data = rcs1.data(:,2);
rcsAbs = [rcs1.data(:,6) rcs2.data(:,6) rcs3.data(:,6)];   % phi direction
rcsPhase = [rcs1.data(:,7) rcs2.data(:,7) rcs3.data(:,7)]; % phi direction
theta_data = theta_data(1:648);
phi_data = phi_data(1:648);
rcsAbs = rcsAbs(1:648, :);
rcsPhase = rcsPhase(1:648, :);

%% Construct theta/phi & azimuth/elevatoin table
global theta phi az el
[theta, phi] = meshgrid(deg2rad(-180:10:170), deg2rad(-90:10:80));
% From theta&phi to azmuth&elevation
el = NaN(18,36);
az = NaN(18,36);
[ind1, ind2] = find(0<=theta & theta<deg2rad(180));
el(ind1, ind2) = deg2rad(90) - theta(ind1, ind2);
az(ind1, ind2) = phi(ind1, ind2);
[ind1, ind2] = find(deg2rad(-180)<=theta & theta<0 & 0>=phi);
el(ind1, ind2) = deg2rad(90) - (- theta(ind1, ind2));
az(ind1, ind2) = phi(ind1, ind2)+deg2rad(180);
[ind1, ind2] = find(deg2rad(-180)<=theta & theta<0 & 0<phi);
el(ind1, ind2) = deg2rad(90) - (- theta(ind1, ind2));
az(ind1, ind2) = phi(ind1, ind2)-deg2rad(180);

%% Construct region logic 
global logic_regionK  % 18*36*K
logic_regionK = false(18, 36, K);
range_theta = (180-(-180))/K;
for k=1:K
     logic_regionK(:,:,k) = theta>=deg2rad(-180+(k-1)*range_theta) & theta<deg2rad(-180+k*range_theta); % 2nd inequality is due to -180~"170"
end
% Plot visualization of region setting
% plotNregionVisualization(logic_regionK, K)

%% Extract rcsAbs for regionK
rcsAbs_r = cell(K,1);
for k=1:K
    ind = find(logic_regionK(:,:,k)');
    rcsAbs_r{k} = rcsAbs(ind, :);
end

%% Optimization K times
p0 = repmat([0; 0; 0; 1; 1; 1; 0; 0], N, 1);
lb = repmat([-10; -10; -10; 1e-5; 1e-5; 1e-5; -45; -45], N, 1);
ub = repmat([10; 10; 10; 15; 15; 15; 45; 45], N, 1);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp')%,'MaxFunctionEvaluations', 10000000,'MaxIterations',1000000);
cost = cell(K, 1);
p_r = zeros(N*8, K);
fval_r = zeros(K, 1);
cost_r = zeros(K, 1);
time_r = zeros(K, 1);
load('opt_PredefineKregions_16each.mat')  %comment out when computing optimization
for k=1:K
    cost{k} =  @(p) sqrt(mean((rcsAbs_r{k} - computeRCS(p, k)).^2, 'all'));
    tic
%     [p_r(:,k), fval_r(k)] = fmincon(cost{k}, p0, [], [], [], [], lb, ub, [], options);
    cost_r(k) = cost{k}(p_r(:,k));
    time_r(k) = toc;
    disp(['Initial RMSE at region', num2str(k),' = ', num2str(cost{k}(p0))])
    disp(['Optimal RMSE at region', num2str(k), ' = ', num2str(cost{k}(p_r(:,k)))])
end
costTotal = sqrt(sum(cost_r.^2)/K);
computationTime = sum(time_r);

%% Combine results from K regions
RCS_model = zeros(648, 3);
for k=1:K
    ind = find(logic_regionK(:,:,k)');
    RCS_model(ind,:) = computeRCS(p_r(:,k), k);
end

%% Save Result
% save('opt_PredefineKregions','p_r')  %comment out when not computing optimization

%% Plot RCS comparison
% 2D plot
plot2dRCScomparison_3freq(rcsAbs, RCS_model, N)

% 3D polar plot
plot3dRCScomparison_3freq(rcsAbs, RCS_model, N)

% 2D polar plot (fix phi= 0&-90)
plot2dRCSpolarComparison_3freq(rcsAbs, RCS_model, N)

% Plot ellipsoid model
for k=1:K
    plotEllipsoidModel(p_r(:,k), N)
    title(strcat('Ellipsoid Model at Region',num2str(k)))
end

% Plot ellipsoid model (closer)
for k=1:K
    plotEllipsoidModel_closer(p_r(:,k), N)
    title(strcat('Ellipsoid Model at Region',num2str(k)))
end

%%
function rcs_model = computeRCS(p, region) % Input:N*8*K variables, Output:Ellipoid model analytical solution table   
    global N theta phi az el logic_regionK
    locations = zeros(N, 3);
    r_ells = zeros(N, 3);
    angle_offsets = zeros(N, 2);
    rcs_models = zeros(18, 36, N);
    
    % RCS table for each ellipsoid 
    for k = 1:N
        locations(k, :) = [p(1 + (k-1)*8), p(2 + (k-1)*8), p(3 + (k-1)*8)];
        r_ells(k, :) = [p(4 + (k-1)*8), p(5 + (k-1)*8), p(6 + (k-1)*8)];
        angle_offsets(k, :) = deg2rad([p(7 + (k-1)*8), p(8+ (k-1)*8)]);       
        rx_ell = r_ells(k, 1);
        ry_ell = r_ells(k, 2);
        rz_ell = r_ells(k, 3);
        phi_ell = phi - angle_offsets(k, 2);
        theta_ell = theta - angle_offsets(k, 1);
        sin_phi_s = sin(phi_ell).^2;
        cos_phi_s = cos(phi_ell).^2;
        rcs_models(:, :, k) = (pi * rx_ell^2 * ry_ell^2 * rz_ell^2) ./ (rx_ell^2 * cos_phi_s .* (sin(theta_ell).^2) + ry_ell^2 * sin_phi_s .* (sin(theta_ell).^2) + rz_ell^2 .* (cos(theta_ell).^2)).^2 ; 
    end
    
    % Combined RCS value of f1
    c = physconst('LightSpeed');
    fc = 1e9;    
    lam = c/fc;
    extrcs1 = zeros(18,36);
    for k = 1:18
        sv = steervec(locations'/lam, [az(k, :); el(k, :)]);
        rcs_models_reshaped = squeeze(rcs_models(k,:,:));
        if N==1
            temp1 = sqrt(rcs_models_reshaped).*(sv.^2); % "squeeze" doesn't effect 2D matrix       
        elseif N>=2
            temp1 = sqrt(rcs_models_reshaped').*(sv.^2);
        end
        extrcs1(k,:) = abs(sum(temp1, 1)).^2;
    end        
    
    % Combined RCS value of f2
    fc = 1.5e9;
    lam = c/fc;
    extrcs2 = zeros(18,36);
    for k = 1:18
        sv = steervec(locations'/lam, [az(k, :); el(k, :)]);
        rcs_models_reshaped = squeeze(rcs_models(k,:,:));
        if N==1
            temp1 = sqrt(rcs_models_reshaped).*(sv.^2); % "squeeze" doesn't effect 2D matrix       
        elseif N>=2
            temp1 = sqrt(rcs_models_reshaped').*(sv.^2);
        end
        extrcs2(k,:) = abs(sum(temp1, 1)).^2;
    end 
    
    % Combined RCS value of f3
    fc = 2e9;
    lam = c/fc;
    extrcs3 = zeros(18,36);
    for k = 1:18
        sv = steervec(locations'/lam, [az(k, :); el(k, :)]);
        rcs_models_reshaped = squeeze(rcs_models(k,:,:));
        if N==1
            temp1 = sqrt(rcs_models_reshaped).*(sv.^2); % "squeeze" doesn't effect 2D matrix       
        elseif N>=2
            temp1 = sqrt(rcs_models_reshaped').*(sv.^2);
        end
        extrcs3(k,:) = abs(sum(temp1, 1)).^2;
    end 
    
    % Output K region result
    rcs_model_f1 = pow2db(reshape(extrcs1',[18*36,1]));
    rcs_model_f2 = pow2db(reshape(extrcs2',[18*36,1]));
    rcs_model_f3 = pow2db(reshape(extrcs3',[18*36,1]));    
    ind = find(logic_regionK(:,:,region)');
    rcs_model = [rcs_model_f1(ind) rcs_model_f2(ind) rcs_model_f3(ind)];
end
