function plot2dRCScomparison(rcsAbs, rcsModel, N)
    figure
    hold on
    plot(rcsAbs, 'LineWidth',4)
    plot(rcsModel, 'LineWidth',2)
    set(gca,'FontSize',20)
    title('RCS Value Comparison')
    xlabel('Angle Points Indices (\phi(-90~80) x \theta(-180~170))')  %(outer loop=\phi (-90~80), inner loop= \theta (-180~170))')
    ylabel('RCS(dB m^2)')
    xlim([1 648])
    for k=1:17
        line([36*k, 36*k],[-20, 50],'Color','black','LineStyle','--')
    end
    legend('True RCS from CST',strcat('Modeled RCS using ', num2str(N), ' ellipsoid(s)'))
    xticks(0:36:648)
end

