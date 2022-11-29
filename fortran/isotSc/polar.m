clc;
close all;
data = load('integral.txt');
psi = data(:,1);
psi = psi/180*pi;
int = data(:,4);
intMod = sqrt(data(:,4).^2 + data(:,7).^2);

data = load('stPhase.txt');
stPhase = data(:,4);
stPhaseMod = sqrt(data(:,4).^2 + data(:,7).^2);

polarplot(psi, intMod, psi, stPhaseMod,'--', 'lineWidth', 3);
thetalim([0 180]);
legend('integral','asymptotics');
title('u_{pp}');