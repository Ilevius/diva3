close all;
data = load('test.txt');
x = data(:,1);
y = data(:,2);
z = data(:,4);
plot(x, y, x, z, 'o');
