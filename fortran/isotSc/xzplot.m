clc;
close all;
data = load('integral.txt');
x = data(:,3);
int = data(:,5);
intMod = sqrt(data(:,7).^2 + data(:,10).^2);

data = load('stPhase.txt');
stPhase = data(:,5);
stPhaseMod = sqrt(data(:,7).^2 + data(:,10).^2);

plot(x, intMod, x, stPhaseMod,'--', 'lineWidth', 3);
legend('integral','asymptotics');
title('u_{pp}');