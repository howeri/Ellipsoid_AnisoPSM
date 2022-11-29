function plot2dRCScomparison_3freq(rcsAbs, RCS_model, N)
    figure
    subplot(3,1,1)
    hold on
    plot(rcsAbs(:,1), 'LineWidth',4)
    plot(RCS_model(:,1), 'LineWidth',2)
    set(gca,'FontSize',10)
    title('RCS Value Comparison at f=1GHz')
    xlabel('Angle Points Indices (\phi(-90~80) x \theta(-180~170))')  %(outer loop=\phi (-90~80), inner loop= \theta (-180~170))')
    ylabel('RCS(dB m^2)')
    for k=1:17
        line([36*k, 36*k],[-40, 60],'Color','black','LineStyle','--')
    end
    legend('True RCS from CST',strcat('Modeled RCS using ',num2str(N),' ellipsoid(s)'))
    xlim([1 648])
    xticks(0:36:648)
    subplot(3,1,2)
    hold on
    plot(rcsAbs(:,2), 'LineWidth',4)
    plot(RCS_model(:,2), 'LineWidth',2)
    set(gca,'FontSize',10)
    title('RCS Value Comparison at f=1.5GHz')
    xlabel('Angle Points Indices (\phi(-90~80) x \theta(-180~170))')  %(outer loop=\phi (-90~80), inner loop= \theta (-180~170))')
    ylabel('RCS(dB m^2)')
    for k=1:17
        line([36*k, 36*k],[-40, 60],'Color','black','LineStyle','--')
    end
    legend('True RCS from CST',strcat('Modeled RCS using ',num2str(N),' ellipsoid(s)'))
    xlim([1 648])
    xticks(0:36:648)
    subplot(3,1,3)
    hold on
    plot(rcsAbs(:,3), 'LineWidth',4)
    plot(RCS_model(:,3), 'LineWidth',2)
    set(gca,'FontSize',10)
    title('RCS Value Comparison at f=2GHz')
    xlabel('Angle Points Indices (\phi(-90~80) x \theta(-180~170))')  %(outer loop=\phi (-90~80), inner loop= \theta (-180~170))')
    ylabel('RCS(dB m^2)')
    for k=1:17
        line([36*k, 36*k],[-40, 60],'Color','black','LineStyle','--')
    end
    legend('True RCS from CST',strcat('Modeled RCS using ',num2str(N),' ellipsoid(s)'))
    xlim([1 648])
    xticks(0:36:648)
end
