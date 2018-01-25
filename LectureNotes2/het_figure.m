clear all;
clc;
FEfig1 = figure('visible','off');
fplot(@(x) x, [1, 5])
hold on
fplot(@(x) .5*x+  1, [1, 5],'-.r')
hold on
fplot(@(x) .5*x+1.5, [1, 5],'-.r')
hold on
fplot(@(x) .5*x+  2, [1, 5],'-.r')
xlabel('X')
ylabel('y')
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'XTick',[])
set(gca,'YTick',[])
legend('Pooled OLS line','FE lines with individual intercepts')
title('OLS vs. Fixed Effects: Overstated Case','FontWeight','bold')
exportfig(FEfig1, 'C:/Users/Tyler & Nichole/Documents/My Dropbox/Teaching/DukeMatlab/Lectures2/FE_figure1.eps')

FEfig2 = figure('visible','off');
fplot(@(x) .5*x+1, [1, 5])
hold on
fplot(@(x) x-   1, [1, 5],'-.r')
hold on
fplot(@(x) x  -.5, [1, 5],'-.r')
hold on
fplot(@(x) x     , [1, 5],'-.r')
xlabel('X')
ylabel('y')
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'XTick',[])
set(gca,'YTick',[])
legend('OLS line','FE line')
legend('Pooled OLS line','FE lines with individual intercepts')
title('OLS vs. Fixed Effects: Understated Case','FontWeight','bold')
exportfig(FEfig2, 'C:/Users/Tyler & Nichole/Documents/My Dropbox/Teaching/DukeMatlab/Lectures2/FE_figure2.eps')