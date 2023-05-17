clc;
close all;
data = load('integral_polar.txt');
psi = data(:,1);
psi = psi/180*pi;
int = data(:,7);
intMod = sqrt(data(:,7).^2 + data(:,10).^2);

data = load('stPhase_polar.txt');
stPhase = data(:,7);
stPhaseMod = sqrt(data(:,7).^2 + data(:,10).^2);

polarplot(psi, intMod, psi, stPhaseMod,'--', 'lineWidth', 3);
thetalim([0 180]);


% legItem1 = '$|\mathbf{u}_{ss}|$ integral';
% legItem2 = '$|\mathbf{u}_{ss}|$ asymptotics';
% leg1 = legend(legItem1, legItem2);
% set(leg1, 'Interpreter', 'Latex');


grid on;
set(gcf,'color','w');

legend('интеграл','асимптотика','Location','northwest');
title('$\mathbf{u}_{ss}$', 'Interpreter', 'Latex', 'FontSize',75);
set(get(gca,'title'),'Position',[0 0])
set(gca, 'GridAlpha', 1);

set(findall(gcf,'type','text'), 'FontSize', 14, 'FontName', 'Times')
set(gca,'FontSize', 14, 'FontName', 'Times')