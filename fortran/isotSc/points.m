close all;
clc;
data = load('points.txt');
x = data(:,1);
z = data(:,2);
xRe = data(:,3);
zRe = data(:,4);
plot(x, z, xRe, zRe, 'x')