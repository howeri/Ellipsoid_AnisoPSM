function plotSphereModel(p, N)
    figure
    AC = stlread('f35.stl');
    trimesh(AC,'FaceColor','none','EdgeColor','k')
    axis equal
    hold on
    for k = 1:N
        [x, y, z] = sphere;
        surf(x*p(4+(k-1)*4)/14.8*120 + 50+p(1+(k-1)*4), y*p(4+(k-1)*4)/14.8*120 + 10+p(2+(k-1)*4), z*p(4+(k-1)*4)/14.8*120 + 40+p(3+(k-1)*4));
    end
    alpha(.3)
    colormap copper
    axis equal
    axis off
    title('Sphere Model')
    set(gca,'FontSize',20)
end

