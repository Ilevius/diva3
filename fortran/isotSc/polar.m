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
legend('integral','asymptotics');
title('u_{pp}');