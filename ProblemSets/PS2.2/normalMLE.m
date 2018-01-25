function [like]=normalMLE(b,X,Y)
beta      = b(1:end-1);
wagesigma = b(end);
like = -sum(-.5*log(2*pi)-.5*log(wagesigma^2)-(.5*(Y-X*beta)./wagesigma).^2);
end