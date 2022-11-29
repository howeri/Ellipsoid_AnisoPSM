function plotNregionVisualization(logic_regionK, K)
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
    legend_char = '';
    for k=1:K
        scatter3(x(logic_regionK(:,:,k)), y(logic_regionK(:,:,k)), z(logic_regionK(:,:,k)), 500, 'filled')
        hold on
        legend_char = strcat(legend_char, sprintf('"Region%d', k), '",');
    end
    title('Region Visualization')
    eval(append('legend(',legend_char(1:end-1),')'))
    axis equal
    set(gca,'FontSize',20) 
end

