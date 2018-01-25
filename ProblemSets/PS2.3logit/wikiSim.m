clear all;
clc;

tic;
int_hat_quad = quad(@(x) exp(-x.^2),-10,1);
quad_time=toc;

n = 2*[1e1 1e2 1e3 1e4 1e5 1e6 1e7]';
int_hat = zeros(size(n));
int_time = zeros(size(n));
for i = 1:length(n)
    tic;
    area = 11;
    grid = [11*rand(n(i),1)-1 rand(n(i),1)];
    below = (grid(:,2)<=exp(-(grid(:,1)).^2))&(grid(:,2)>0);
    above = (grid(:,2)>=exp(-(grid(:,1)).^2))&(grid(:,2)<0);
    int_hat(i) = area*(sum(below)-sum(above))/length(grid);
    int_time(i)= toc;
end

% general formula for any function over any interval is
%
% grid = [U(xmin to xmax) U(ymin to ymax)]
% area  = (xmax-xmin)*(ymax-ymin)
% below = (y<=f(x))&(y>0), or
% below = (grid(:,2)<=f(grid(:,1))&(grid(:,2)>0)
% above = (y>=f(x))&(y<0), or
% above = (grid(:,2)>=f(grid(:,1))&(grid(:,2)<0)
% integral_hat = (area of box)*(sum(below)-sum(above))/length(grid)