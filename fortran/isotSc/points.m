close all;
data = load('points.txt');
x = data(:,1);
z = data(:,2);
plot(x, z, 'x')