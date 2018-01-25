function [like,grad,hess]=normalMLEhessian(b,X,Y)
beta      = b(1:end-1);
wagesigma = b(end);
n         = length(Y);
like = -sum(-.5*log(2*pi)-.5*log(wagesigma^2)-.5*((Y-X*beta)./wagesigma).^2);

grad = zeros(size(b));
grad(1:end-1) = -X'*(Y-X*beta)./wagesigma^2;
grad(end)     = n./wagesigma - ((Y-X*beta)'*(Y-X*beta))./wagesigma^3;

hess = zeros(length(b));
hess(1:length(b)-1,1:length(b)-1) = (X'*X)./wagesigma.^2;
hess(length(b),1:length(b)-1) = (2*(Y-X*beta)'*X)./wagesigma.^3;
hess(1:length(b)-1,length(b)) = (2*X'*(Y-X*beta))./wagesigma.^3;
hess(length(b),length(b)) = (-n)./wagesigma.^2 + (3./wagesigma.^4)*(Y-X*beta)'*(Y-X*beta);
end