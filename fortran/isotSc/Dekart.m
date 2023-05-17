clc;
close all;
fh = figure();
fh.WindowState = 'maximized';

data = load('integral.txt');
x = data(:,3);
int = data(:,10);
intMod = sqrt(data(:,7).^2 + data(:,10).^2);

data = load('stPhase.txt');
stPhase = data(:,10);
stPhaseMod = sqrt(data(:,7).^2 + data(:,10).^2);

error = abs(intMod-stPhaseMod)./intMod;

plot(x, int, x, stPhase,'--', 'lineWidth', 1.5);

%plot(x, error, 'lineWidth', 1.5);

% legItem1 = '$|\mathbf{u}_{ss}|$ integral';
% legItem2 = '$|\mathbf{u}_{ss}|$ asymptotics';
% leg1 = legend(legItem1, legItem2);
% set(leg1, 'Interpreter', 'Latex');


%grid on;
set(gcf,'color','w');

legend('интеграл','асимптотика','Location','northwest');
%title('$\mathbf{u}^+_{ss}$', 'Interpreter', 'Latex', 'FontSize',75);
xlabel('$x$', 'FontSize',22, 'Interpreter', 'LaTeX');
%set(get(gca,'title'),'Position',[0 0])
set(gca, 'GridAlpha', 1);

set(findall(gcf,'type','text'), 'FontSize', 65, 'FontName', 'Times New Roman')
set(gca,'FontSize', 65, 'FontName', 'Times New Roman')
pbaspect([1.62 1 1])
set(gcf,'units','centimeters','position',[0,0,6.2,3.8])
%ylim([0.0000005 inf])
%xlim([1 50])