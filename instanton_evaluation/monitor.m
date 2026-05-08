function mon = monitor(it,xf,E,mod_dx,delta,s)

    figure(1)
    set(gcf,'position',[10 45 650 600])


    subplot(2,2,1)
    hold on
    plot(it,xf,'.b')
    title('$x_1(s=1)$','Interpreter','Latex','FontSize',12)
    drawnow;

    %plot(s,x(2,:),'-.')
    %plot(it,alpha,'s')
    subplot(2,2,2)
    hold on;
    plot(s,mod_dx,'-.')
    title('Speed','Interpreter','Latex','FontSize',12)
    drawnow;

    subplot(2,2,3)
    hold on;
    plot(it,delta,'.k')
    title('Error','Interpreter','Latex','FontSize',12)

    subplot(2,2,4)
    hold on;
    plot(s,E,'-.r')
    title('$\textrm{Hamiltonian}$','Interpreter','Latex','FontSize',12)
%     hold off
    mon = 1;
end
