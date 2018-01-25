function like = mlogit(b,X,y)
b1 = b(1:size(X,2));
b2 = b(1+size(X,2):end);
dem = 1+exp(X*b1)+exp(X*b2);
like = -sum((y==1).*(X*b1) + (y==2).*(X*b2) + (y==3) - log(dem));