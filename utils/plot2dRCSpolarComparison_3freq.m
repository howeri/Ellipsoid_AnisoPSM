function plot2dRCSpolarComparison_3freq(rcsAbs, RCS_model, N)
    figure
    theta = deg2rad(-180:10:180);
    maxRCS = max([rcsAbs(:, 1) RCS_model(:,1)],[],'all');
    minRCS = min([rcsAbs(:, 1) RCS_model(:,1)],[],'all');
    subplot(1,2,1)
    rcs_true = [rcsAbs(325:360, 1); rcsAbs(325, 1)];
    rcs_model = [RCS_model(325:360, 1); RCS_model(325, 1)];
    polarplot(theta,rcs_true, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    hold on
    polarplot(theta,rcs_model, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    set(gca,'FontSize',20)
    legend('True RCS from CST', strcat('Modeled RCS using ',num2str(N),' ellipsoid(s)'))
    title('RCS Plot at \phi=0 at f=1GHz (m^2 dB)')
    subplot(1,2,2)
    rcs_true = [rcsAbs(1:36, 1); rcsAbs(1, 1)];
    rcs_model = [RCS_model(1:36, 1); RCS_model(1, 1)];
    polarplot(theta,rcs_true, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    hold on
    polarplot(theta,rcs_model, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    set(gca,'FontSize',20)
    legend('True RCS from CST', strcat('Modeled RCS using ',num2str(N),' ellipsoid(s)'))
    title('RCS Plot at \phi=-90 at f=1GHz (m^2 dB)')

    figure
    theta = deg2rad(-180:10:180);
    maxRCS = max([rcsAbs(:, 2) RCS_model(:,2)],[],'all');
    minRCS = min([rcsAbs(:, 2) RCS_model(:,2)],[],'all');
    subplot(1,2,1)
    rcs_true = [rcsAbs(325:360, 2); rcsAbs(325, 2)];
    rcs_model = [RCS_model(325:360, 2); RCS_model(325, 2)];
    polarplot(theta,rcs_true, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    hold on
    polarplot(theta,rcs_model, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    set(gca,'FontSize',20)
    legend('True RCS from CST', strcat('Modeled RCS using ',num2str(N),' ellipsoid(s)'))
    title('RCS Plot at \phi=0 at f=1.5GHz (m^2 dB)')
    subplot(1,2,2)
    rcs_true = [rcsAbs(1:36, 2); rcsAbs(1, 2)];
    rcs_model = [RCS_model(1:36, 2); RCS_model(1, 2)];
    polarplot(theta,rcs_true, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    hold on
    polarplot(theta,rcs_model, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    set(gca,'FontSize',20)
    legend('True RCS from CST', strcat('Modeled RCS using ',num2str(N),' ellipsoid(s)'))
    title('RCS Plot at \phi=-90 at f=1.5GHz (m^2 dB)')

    figure
    theta = deg2rad(-180:10:180);
    maxRCS = max([rcsAbs(:, 3) RCS_model(:,3)],[],'all');
    minRCS = min([rcsAbs(:, 3) RCS_model(:,3)],[],'all');
    subplot(1,2,1)
    rcs_true = [rcsAbs(325:360, 3); rcsAbs(325, 3)];
    rcs_model = [RCS_model(325:360, 3); RCS_model(325, 3)];
    polarplot(theta,rcs_true, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    hold on
    polarplot(theta,rcs_model, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    set(gca,'FontSize',20)
    legend('True RCS from CST', strcat('Modeled RCS using ',num2str(N),' ellipsoid(s)'))
    title('RCS Plot at \phi=0 at f=2GHz (m^2 dB)')
    subplot(1,2,2)
    rcs_true = [rcsAbs(1:36, 3); rcsAbs(1, 3)];
    rcs_model = [RCS_model(1:36, 3); RCS_model(1, 3)];
    polarplot(theta,rcs_true, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    hold on
    polarplot(theta,rcs_model, 'Linewidth',4)
    rlim([minRCS-10 maxRCS+10])
    set(gca,'FontSize',20)
    legend('True RCS from CST', strcat('Modeled RCS using ',num2str(N),' ellipsoid(s)'))
    title('RCS Plot at \phi=-90 at f=2GHz (m^2 dB)')
end

