function plot2RegionVisualization(logic_region2)
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
    
    figure
    [x,y,z] = sph2cart(az, el, ones(18,36));
    scatter3(x(~logic_region2),y(~logic_region2),z(~logic_region2),500,'MarkerFaceColor', 'r')
    hold on
    scatter3(x(logic_region2),y(logic_region2),z(logic_region2),500,'MarkerFaceColor', 'b')
    title('Region Visualization')
    legend('Region1','Region2')
    axis equal
    set(gca,'FontSize',20) 
end

