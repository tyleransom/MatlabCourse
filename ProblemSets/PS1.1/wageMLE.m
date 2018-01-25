function [like]=wageMLE(b,Y,X)
beta      = b(1:end-1);
wagesigma = b(end);
% wagesigma = exp(wagesigma);
% ll = log(normpdf(Y-X*beta,0,wagesigma));
like = sum(log(wagesigma^2)+((Y-X*beta)/wagesigma).^2);