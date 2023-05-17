close all;
data = load('test.txt');
x = data(:,1);
y = data(:,2);
z = data(:,4);
plot(x, y, x, z, 'o', 'lineWidth', 2);
legend('left hand', 'right hand');
xlabel('\alpha')
grid on;
set(gcf,'color','w');
set(gca,'FontSize', 65, 'FontName', 'Times New Roman')
