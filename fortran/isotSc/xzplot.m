clc;
close all;

data = load('integral.txt');
x = data(:,3);
intU = data(:,7);
intW = data(:,10);
intMod = sqrt(data(:,7).^2 + data(:,10).^2);

fieldCode = data(1,11);
cp1 = data(1,12); cp2 = data(1,13);
cs1 = data(1,14); cs2 = data(1,15);
rho1 = data(1,16); rho2 = data(1,17);
w = data(1,19);
lambda = min(cp1, cp2)*2*pi/w;
xlabelContent = "$x, \lambda = "+lambda+"$";
titleContent = "$ \omega = " + w + ", cp(1) = " + cp1 + ", cp(2) = " + cp2;
titleContent = titleContent + ", cs(1) = " + cs1+", cs(2) = " + cs2+ ", \rho(1) = " + rho1+", \rho(2) = " + rho2 + "$";






data = load('stPhase.txt');
stPhaseU = data(:,7);
stPhaseW = data(:,10);
stPhaseMod = sqrt(data(:,7).^2 + data(:,10).^2);




if fieldCode == 11 
    subfield = '_{pp}';
    pictureName = 'upp';
elseif fieldCode == 12
    subfield = '_{ps}';  
    pictureName = 'ups';
elseif fieldCode == 21
    subfield = '_{sp}';
    pictureName = 'usp';
elseif fieldCode == 22
    subfield = '_{ss}';
    pictureName = 'uss';
else
    title('check field code!');
end

pictureName = pictureName+ "w"+w+"cp(1)"+cp1+"cp(2)"+cp2+"cs(1)"+cs1+"cs(2"+cs2+"rho(1)"+rho1+"rho(2)"+rho2;

pictureName = strrep(pictureName, '.', ',')

body = "\begin{figure}[h]\centering\includegraphics[scale=0.38]{"+pictureName+".eps";
body = body + "}\caption{";
body = body +  titleContent + "}\end{figure}";

workFile = fopen('plots.tex','a');
fprintf(workFile, '%s', body);
fclose(workFile);


f = figure;
f.WindowState = 'maximized';

figure(f);
%                                       overal module
subplot(2,2,1);
field = "\textbf{u}";
yContent = "$|"+field + subfield+"|$";
plot(x, intMod, x, stPhaseMod,'--', 'lineWidth', 3);
legend('integral','asymptotics');
%title(titleContent,'interpreter','latex');
xlabel(xlabelContent, 'FontSize',22, 'Interpreter', 'LaTeX');
ylabel(yContent, 'FontSize',22, 'Interpreter', 'LaTeX');
set(gca, 'FontSize',22);

%                                   horisontal component module
subplot(2,2,2);
field = "u";
yContent = "$|"+field + subfield+"|$";
plot(x, intU, x, stPhaseU,'--', 'lineWidth', 3);
legend('integral','asymptotics');
%title(titleContent,'interpreter','latex');
xlabel(xlabelContent, 'FontSize',22, 'Interpreter', 'LaTeX');
ylabel(yContent, 'FontSize',22, 'Interpreter', 'LaTeX');
set(gca, 'FontSize',22);


%                                       vertical component module
subplot(2,2,3);
field = "w";
yContent = "$|"+field + subfield+"|$";
plot(x, intW, x, stPhaseW,'--', 'lineWidth', 3);
legend('integral','asymptotics');
%title(titleContent,'interpreter','latex');
xlabel(xlabelContent, 'FontSize',22, 'Interpreter', 'LaTeX');
ylabel(yContent, 'FontSize',22, 'Interpreter', 'LaTeX');
set(gca, 'FontSize',22);set(gca, 'FontSize',22);


grid on;
set(gcf,'color','w');

saveas(gcf, pictureName, 'epsc');



