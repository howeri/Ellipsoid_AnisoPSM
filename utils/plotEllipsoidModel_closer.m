function plotEllipsoidModel_closer(p, N)
    figure
    axis equal
    scatter3(0,0,0,500,'+','MarkerFaceColor', 'b')
    hold on
    for k = 1:N
        [x, y, z] = ellipsoid(p(1+(k-1)*8), p(2+(k-1)*8), p(3+(k-1)*8), p(4+(k-1)*8)/10, p(5+(k-1)*8)/10, p(6+(k-1)*8)/10);  %scale the radius accordingly
        S = surf(x, y, z); 
        rotate(S, [0, 0, 1], p(8+(k-1)*8))
        rotate(S, [-tan(p(8+(k-1)*8)), 1, 0], p(7+(k-1)*8))
    end
    alpha(.3)
    colormap copper
    axis equal
    axis off
    title('Ellipsoid Model (closer look)')
    set(gca,'FontSize',20)
end

