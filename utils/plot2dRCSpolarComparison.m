function plot2dRCSpolarComparison(rcsAbs, RCS_model, N)
    figure
    theta3 = deg2rad(-180:10:180);
    maxRCS = max([rcsAbs RCS_model], [], 'all');
    minRCS = min([rcsAbs RCS_model], [], 'all');
    subplot(1,2,1)
    rcs_true = [rcsAbs(325:360); rcsAbs(325)];
    rcs_model = [RCS_model(325:360); RCS_model(325)];
    polarplot(theta3, rcs_true, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    hold on
    polarplot(theta3, rcs_model, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    set(gca,'FontSize',20)
    legend('True RCS from CST',strcat('Modeled RCS using ', num2str(N), ' ellipsoid(s)'))
    title('RCS Plot at \phi=0(m^2 dB)')
    subplot(1,2,2)
    rcs_true = [rcsAbs(1:36); rcsAbs(1)];
    rcs_model = [RCS_model(1:36); RCS_model(1)];
    polarplot(theta3, rcs_true, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    hold on
    polarplot(theta3,rcs_model, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    set(gca,'FontSize',20)
    legend('True RCS from CST',strcat('Modeled RCS using ', num2str(N), ' ellipsoid(s)'))
    title('RCS Plot at \phi=-90(m^2 dB)')
end

