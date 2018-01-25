clear all;
clc;
% bsxfun practice
N        = 8984;
step     = .03;
variance = 2*rand(N,1)+.25;
mean     = 4*rand(N,1)-2;
a1       = mean-5*sqrt(variance);
b1       = mean+5*sqrt(variance);
G        = floor((b1-a1)./step);

grid    = zeros(max(G),N);
weights = zeros(max(G),N);
tic
for i=1:N
	[grid(1:G(i),i),weights(1:G(i),i)]=lgwt(G(i),a1(i),b1(i));
end
disp('time to loop through quadrature weights')
toc

f   = normpdf(grid',repmat(mean,1,max(G)),repmat(variance,1,max(G)));
int = sum(normpdf(grid').*weights',2);
Ex  = sum((grid').*f.*weights',2);