function [like,grad] = mlogit_gradient(b,X,y)
b1 = b(1+0*size(X,2):1*size(X,2));
b2 = b(1+1*size(X,2):2*size(X,2));
b3 = b(1+2*size(X,2):3*size(X,2));
b4 = b(1+3*size(X,2):4*size(X,2));
b5 = b(1+4*size(X,2):5*size(X,2));
dem = 1+exp(X*b1)+exp(X*b2)+exp(X*b3)+exp(X*b4)+exp(X*b5);
like = -sum((y==1).*(X*b1) + (y==2).*(X*b2) + (y==3).*(X*b3) + (y==4).*(X*b4) + (y==5).*(X*b5) + (y==6) - log(dem));
P1 = exp(X*b1)./dem;
P2 = exp(X*b2)./dem;
P3 = exp(X*b3)./dem;
P4 = exp(X*b4)./dem;
P5 = exp(X*b5)./dem;
grad = zeros(size(b));
grad(1+0*size(X,2):1*size(X,2),1) = -X'*((y==1)-P1);
grad(1+1*size(X,2):2*size(X,2),1) = -X'*((y==2)-P2);
grad(1+2*size(X,2):3*size(X,2),1) = -X'*((y==3)-P3);
grad(1+3*size(X,2):4*size(X,2),1) = -X'*((y==4)-P4);
grad(1+4*size(X,2):5*size(X,2),1) = -X'*((y==5)-P5);