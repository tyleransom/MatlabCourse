function [like,grad,hess] = logit_gradient(b,X,y)
dem = 1+exp(X*b);
like = -sum((y==1).*(X*b) + (y==2) - log(dem));
P1 = exp(X*b)./dem;
grad = zeros(size(b));
grad(1:size(X,2),1)     = -X'*((y==1)-P1);
hess = zeros(length(b));
hess(1:size(X,2),1:size(X,2))         = -X'*(P1)*(1-P1)'*X;