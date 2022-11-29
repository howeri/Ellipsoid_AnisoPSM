function plot3dRCScomparison(rcsAbs, RCS_model, N)
    figure
    [theta2, phi2] = meshgrid(deg2rad(-180:10:170), deg2rad(-90:10:90)); % Add extra data to complete the polar plot
    el2 = NaN(19,36);
    az2 = NaN(19,36);
    [ind1, ind2] = find(0<=theta2 & theta2<deg2rad(180));
    el2(ind1, ind2) = deg2rad(90) - theta2(ind1, ind2);
    az2(ind1, ind2) = phi2(ind1, ind2);
    [ind1, ind2] = find(deg2rad(-180)<=theta2 & theta2<0 & 0>=phi2);
    el2(ind1, ind2) = deg2rad(90) - (- theta2(ind1, ind2));
    az2(ind1, ind2) = phi2(ind1, ind2)+deg2rad(180);
    [ind1, ind2] = find(deg2rad(-180)<=theta2 & theta2<0 & 0<phi2);
    el2(ind1, ind2) = deg2rad(90) - (- theta2(ind1, ind2));
    az2(ind1, ind2) = phi2(ind1, ind2)-deg2rad(180);
    R_true = zeros(19, 36);
    R_model = zeros(19, 36);
    rcsAbs_added = [rcsAbs; rcsAbs(1); rcsAbs(36:-1:2)]; % Add extra data to complete the polar plot
    RCS_model_opt_added = [RCS_model; RCS_model(1); RCS_model(36:-1:2)];  % Add extra data to complete the polar plot
    for k=1:19
        R_true(k, :) = rcsAbs_added( 1+(k-1)*36 : 36*k);
        R_model(k, :) = RCS_model_opt_added( 1+(k-1)*36 : 36*k);
    end
    minthresh = max(R_true(:))-100;           % Define floor (dB)
    R_true(R_true<minthresh) = minthresh;     % Replace -Inf with minthresh
    r_true = R_true-minthresh;                % Radius must be positive
    [x,y,z] = sph2cart(az2,el2,r_true);
    surf(x,y,z,R_true);
    axis equal
    axis off
    title('True RCS from CST')
    set(gca,'FontSize',20)
    hold on
    figure
    minthresh = max(R_model(:))-100;
    R_model(R_model<minthresh) = minthresh;
    r = R_model-minthresh;              
    [x,y,z] = sph2cart(az2,el2,r);
    surf(x,y,z,R_model);
    title(strcat('Modeled RCS using ',num2str(N),' ellipsoid(s)'))
    set(gca,'FontSize',20)
    axis equal
    axis off
end

