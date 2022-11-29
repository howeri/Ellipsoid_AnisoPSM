function plotEllipsoidModel(p, N)
    figure
    AC = stlread('f35.stl');
    trimesh(AC,'FaceColor','none','EdgeColor','k')
    axis equal
    hold on
    for k = 1:N
        [x, y, z] = ellipsoid(p(1+(k-1)*8)+50, p(2+(k-1)*8)+10, p(3+(k-1)*8)+40, p(4+(k-1)*8)/14.8*120, p(5+(k-1)*8)/14.8*120, p(6+(k-1)*8)/14.8*120);
        S = surf(x, y, z);
        rotate(S, [0, 0, 1], p(8+(k-1)*8))
        rotate(S, [-tan(p(8+(k-1)*8)), 1, 0], p(7+(k-1)*8))
    end
    alpha(.3)
    colormap copper
    axis off
    axis equal
    title('Ellipsoid Model')
    set(gca,'FontSize',20)
end

