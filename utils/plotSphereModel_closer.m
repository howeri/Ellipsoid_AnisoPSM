function plotSphereModel_closer(p, N)
    figure
    axis equal
    scatter3(0,0,0,100,'+','MarkerFaceColor', 'b')
    hold on        
    for k = 1:N
        scatter3(p(1+(k-1)*4), p(2+(k-1)*4), p(3+(k-1)*4), 200, p(4+(k-1)*4), 'filled');  %scale the radius accordingly
    end
    title('Sphere Model (closer look)')
    axis equal
    set(gca,'FontSize',20)
    
    figure
    axis equal
    scatter3(0,0,0,100,'+','MarkerFaceColor', 'b')
    hold on        
    for k = 1:N
        [x, y, z] = sphere;
        surf(x*p(4+(k-1)*4)/10 + p(1+(k-1)*4), y*p(4+(k-1)*4)/10 + p(2+(k-1)*4), z*p(4+(k-1)*4)/10 + p(3+(k-1)*4));  %scale the radius accordingly
    end          
    alpha(.3)
    colormap copper
    axis equal
    axis off
    title('Sphere Model (closer look)')
    set(gca,'FontSize',20)
end

