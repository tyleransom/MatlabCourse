%%% plot (and save) results
% reservation wage as f(lambdas)
rw1 = figure('visible','off');
surf(res_wage1);
hold on
mesh(bmat);
xlabel('\lambda^e')
ylabel('\lambda^u')
zlabel('reservation (log) wage')
% title('Reservation Wage as a function of \lambda^u and \lambda^e',... 
%   'FontWeight','bold')
set(gca,'XTickLabel',{'0','0.25','0.50','0.75','1.00'})
set(gca,'YTickLabel',{'0','0.25','0.50','0.75','1.00'})
% txstr(1) = {'The flat plane denotes the value of (log) b'};
% text(0,2,b-.5,txstr,'HorizontalAlignment','center')
% print(rw1,'-dps','res_wage.ps')
exportfig(rw1, 'C:/Users/Tyler & Nichole/Documents/My Dropbox/Dissertation/Writing/res_wage.eps')
hold off