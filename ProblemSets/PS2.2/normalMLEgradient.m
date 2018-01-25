function [like,grad]=normalMLEgradient(b,X,Y)
beta      = b(1:end-1);
wagesigma = b(end);
n         = length(Y);

like = -sum(-.5*log(2*pi)-.5*log(wagesigma^2)-.5*((Y-X*beta)./wagesigma).^2);

grad = zeros(size(b));
grad(1:end-1) = -X'*(Y-X*beta)./wagesigma^2;
grad(end)     = n./wagesigma - ((Y-X*beta)'*(Y-X*beta))./wagesigma^3;
end