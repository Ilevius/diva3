clc;
close all;
fh = figure();
fh.WindowState = 'maximized';

data = load('COMSOL\WITH reflection dataset.txt');
xwithRef = data(:,1);
withRef = data(:,3);

data = load('COMSOL\No reflection dataset.txt');
xnoRef = data(:,1);
noRef = data(:,3);

data = load("C:\Users\tiama\OneDrive\–абочий стол\IMMI\DIVA3\fortran\isotSc\integral.txt");
xint = data(:,3)./1000;
int = data(:,10)./1000;



 %plot(xwithRef, withRef, xnoRef, noRef,'--', 'lineWidth', 1.5);
% legend('с отражением','без отражени€','Location','northwest');

plot(xwithRef, withRef-noRef, xint, int, 'lineWidth', 1.5);


%plot(x, error, 'lineWidth', 1.5);

% legItem1 = '$|\mathbf{u}_{ss}|$ integral';
% legItem2 = '$|\mathbf{u}_{ss}|$ asymptotics';
% leg1 = legend(legItem1, legItem2);
% set(leg1, 'Interpreter', 'Latex');


%grid on;
set(gcf,'color','w');


%title('$\mathbf{u}^+_{ss}$', 'Interpreter', 'Latex', 'FontSize',75);
%xlabel('$x$', 'FontSize',22, 'Interpreter', 'LaTeX');
%set(get(gca,'title'),'Position',[0 0])
set(gca, 'GridAlpha', 1);

set(findall(gcf,'type','text'), 'FontSize', 65, 'FontName', 'Times New Roman')
set(gca,'FontSize', 65, 'FontName', 'Times New Roman')
%pbaspect([1.62 1 1])
%set(gcf,'units','centimeters','position',[0,0,6.2,3.8])
%ylim([0.0000005 inf])
%xlim([1 50])