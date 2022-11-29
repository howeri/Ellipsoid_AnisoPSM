function plot3dRCScomparison_3freq(rcsAbs, RCS_model, N)
    figure
    subplot(1,2,1)
    [theta, phi] = meshgrid(deg2rad(-180:10:170), deg2rad(-90:10:90));
    el = NaN(19,36);
    az = NaN(19,36);
    [ind1, ind2] = find(0<=theta & theta<deg2rad(180));
    el(ind1, ind2) = deg2rad(90) - theta(ind1, ind2);
    az(ind1, ind2) = phi(ind1, ind2);
    [ind1, ind2] = find(deg2rad(-180)<=theta & theta<0 & 0>=phi);
    el(ind1, ind2) = deg2rad(90) - (- theta(ind1, ind2));
    az(ind1, ind2) = phi(ind1, ind2)+deg2rad(180);
    [ind1, ind2] = find(deg2rad(-180)<=theta & theta<0 & 0<phi);
    el(ind1, ind2) = deg2rad(90) - (- theta(ind1, ind2));
    az(ind1, ind2) = phi(ind1, ind2)-deg2rad(180);
    R_true = zeros(19, 36);
    R_model = zeros(19, 36);
    rcsAbs_add = [rcsAbs(:, 1); rcsAbs(1, 1); rcsAbs(36:-1:2, 1)]; % Add extra data to complete the polar plot
    RCS_model_final_add = [RCS_model(:, 1); RCS_model(1, 1); RCS_model(36:-1:2, 1)]; % Add extra data to complete the polar plot
    for k=1:19
        R_true(k, :) = (rcsAbs_add( 1+(k-1)*36 : 36*k));
        R_model(k, :) = (RCS_model_final_add( 1+(k-1)*36 : 36*k));
    end       
    minthresh = max(R_true(:))-100;           % Define floor (dB)
    R_true(R_true<minthresh) = minthresh;     % Replace -Inf with minthresh
    r_true = R_true-minthresh;                % Radius must be positive
    [x,y,z] = sph2cart(az,el,r_true);
    surf(x,y,z,R_true);
    axis equal
    axis off
    title('True RCS from CST at f=1GHz')
    hold on
    subplot(1,2,2)
    R_model(R_model<minthresh) = minthresh;
    r = R_model-minthresh;              
    [x,y,z] = sph2cart(az,el,r);
    surf(x,y,z,R_model);
    title(strcat('Modeled RCS using ',num2str(N),' ellipsoid(s) at f=1GHz'))
    axis equal
    axis off

    figure
    subplot(1,2,1)
    [theta, phi] = meshgrid(deg2rad(-180:10:170), deg2rad(-90:10:90));
    el = NaN(19,36);
    az = NaN(19,36);
    [ind1, ind2] = find(0<=theta & theta<deg2rad(180));
    el(ind1, ind2) = deg2rad(90) - theta(ind1, ind2);
    az(ind1, ind2) = phi(ind1, ind2);
    [ind1, ind2] = find(deg2rad(-180)<=theta & theta<0 & 0>=phi);
    el(ind1, ind2) = deg2rad(90) - (- theta(ind1, ind2));
    az(ind1, ind2) = phi(ind1, ind2)+deg2rad(180);
    [ind1, ind2] = find(deg2rad(-180)<=theta & theta<0 & 0<phi);
    el(ind1, ind2) = deg2rad(90) - (- theta(ind1, ind2));
    az(ind1, ind2) = phi(ind1, ind2)-deg2rad(180);
    R_true = zeros(19, 36);
    R_model = zeros(19, 36);
    rcsAbs_add = [rcsAbs(:, 2); rcsAbs(1, 2); rcsAbs(36:-1:2, 2)]; % Add extra data to complete the polar plot
    RCS_model_final_add = [RCS_model(:, 2); RCS_model(1, 2); RCS_model(36:-1:2, 2)];  % Add extra data to complete the polar plot
    for k=1:19
        R_true(k, :) = (rcsAbs_add( 1+(k-1)*36 : 36*k));
        R_model(k, :) = (RCS_model_final_add( 1+(k-1)*36 : 36*k));
    end
    R_true(R_true<minthresh) = minthresh;     % Replace -Inf with minthresh
    r_true = R_true-minthresh;                % Radius must be positive
    [x,y,z] = sph2cart(az,el,r_true);
    surf(x,y,z,R_true);
    axis equal
    axis off
    title('True RCS from CST at f=1.5GHz')
    hold on
    subplot(1,2,2)
    R_model(R_model<minthresh) = minthresh;
    r = R_model-minthresh;              
    [x,y,z] = sph2cart(az,el,r);
    surf(x,y,z,R_model);
    title(strcat('Modeled RCS using ',num2str(N),' ellipsoid(s) at f=1.5GHz'))
    axis equal
    axis off

    figure
    subplot(1,2,1)
    [theta, phi] = meshgrid(deg2rad(-180:10:170), deg2rad(-90:10:90));
    el = NaN(19,36);
    az = NaN(19,36);
    [ind1, ind2] = find(0<=theta & theta<deg2rad(180));
    el(ind1, ind2) = deg2rad(90) - theta(ind1, ind2);
    az(ind1, ind2) = phi(ind1, ind2);
    [ind1, ind2] = find(deg2rad(-180)<=theta & theta<0 & 0>=phi);
    el(ind1, ind2) = deg2rad(90) - (- theta(ind1, ind2));
    az(ind1, ind2) = phi(ind1, ind2)+deg2rad(180);
    [ind1, ind2] = find(deg2rad(-180)<=theta & theta<0 & 0<phi);
    el(ind1, ind2) = deg2rad(90) - (- theta(ind1, ind2));
    az(ind1, ind2) = phi(ind1, ind2)-deg2rad(180);
    R_true = zeros(19, 36);
    R_model = zeros(19, 36);
    rcsAbs_add = [rcsAbs(:, 3); rcsAbs(1, 3); rcsAbs(36:-1:2, 3)]; % Add extra data to complete the polar plot
    RCS_model_final_add = [RCS_model(:, 3); RCS_model(1, 3); RCS_model(36:-1:2, 3)];  % Add extra data to complete the polar plot
    for k=1:19
        R_true(k, :) = (rcsAbs_add( 1+(k-1)*36 : 36*k));
        R_model(k, :) = (RCS_model_final_add( 1+(k-1)*36 : 36*k));
    end
    R_true(R_true<minthresh) = minthresh;     % Replace -Inf with minthresh
    r_true = R_true-minthresh;                % Radius must be positive
    [x,y,z] = sph2cart(az,el,r_true);
    surf(x,y,z,R_true);
    axis equal
    axis off
    title('True RCS from CST at f=2GHz')
    hold on
    subplot(1,2,2)
    R_model(R_model<minthresh) = minthresh;
    r = R_model-minthresh;              
    [x,y,z] = sph2cart(az,el,r);
    surf(x,y,z,R_model);
    title(strcat('Modeled RCS using ',num2str(N),' ellipsoid(s) at f=2GHz'))
    axis equal
    axis off
end

