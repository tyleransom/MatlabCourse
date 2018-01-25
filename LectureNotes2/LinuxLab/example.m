clear all;
clc;
delete example.diary
diary('example.diary')
x = 1;
y = sqrt(2);
z = x*y;
% if z>=sqrt(2)
%     error('z should be equal to sqrt(2)')
% end
% keyboard
X = rand(15,3);
Y = rand(3,5);
Z = X*Y
% exit
h = @(x) x.^2 + 2.*x +3;

Z = product2(X,Y,h,x,y,z)

h(2)
diary off