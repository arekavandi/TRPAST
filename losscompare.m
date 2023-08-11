clc
clear all
close all

x=linspace(-6,6,100);
a=[0.98 0.85 0.5];
mu=2;
for j=1:3
for i=1:length(x)
    out(i,j)=(1/(1-a(j)))*(1-exp(0.5*(a(j)-1)*(x(i)^2)));
    outl2(i)=(x(i))^2;
    outl1(i)=abs(x(i));
end
end
figure;
plot((x),(outl2),'LineWidth',2)
hold on
plot((x),(out(:,1)),'LineWidth',2)
hold on
plot((x),(out(:,2)),'LineWidth',2)
hold on
plot((x),(out(:,3)),'LineWidth',2)
hold on
plot((x),(outl1),'LineWidth',3)
grid on
legend({'$\ell_2$-norm','$\ell_{\alpha}$ with $\alpha=0.98$','$\ell_{\alpha}$ with $\alpha=0.85$','$\ell_{\alpha}$ with $\alpha=0.5$','$\ell_1$-norm'}, ...
        'Interpreter', 'LaTeX')
    xlabel('$norm(z)$', 'Interpreter', 'LaTeX')
    ylabel('Loss Value', 'Interpreter', 'LaTeX')
    title('', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')
% print('losscompare','-depsc') 