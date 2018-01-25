function like = binaryMLE(b,X,y,type)
% BINARYMLE(b,X,y,type) estimates a simple logit or probit model for
% regressors in X and dependent variable y
% 'type' is a string variable (either 'logit' or 'probit') which dictates
% the distributional assumption of the error term
switch type
    case {'logit',{}}
        logP0 =      -log(1+exp(X*b));
        logP1 = (X*b)-log(1+exp(X*b));
    case 'probit'
        logP0 = log(1-normcdf(X*b,0,1));
        logP1 = log(  normcdf(X*b,0,1));
        % try to trap errors where Matlab tries to evaluate log(0):
        logP0(1-normcdf(X*b,0,1)==0)=log(1e-220);
        logP1(  normcdf(X*b,0,1)==0)=log(1e-220);
end
like = -sum((y==0).*logP0 + (y==1).*logP1);
