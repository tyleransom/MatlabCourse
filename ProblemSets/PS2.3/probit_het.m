function like = probit_het(b,X,y)
%PROBIT_HET Computes random effects for a binary probit
% [like] = probit_het(b,X,y)
% Estimates parameters b of a random-effects probit given data X and y
% Note that the std dev of the random effect is the last element of b
% Estimation method is a Matlab built-in quadrature program
T = size(X,3);
like_i = zeros(length(X),1);
for i=1:length(X)
    like_i(i) = quadv(@(z) prod(normcdf(ctranspose(ctranspose(squeeze(X(i,:,:)))*b(1:end-1))+z*ones(1,T)).^(y(i,:)==1) .* (1-normcdf(ctranspose(ctranspose(squeeze(X(i,:,:)))*b(1:end-1))+z*ones(1,T))).^(y(i,:)==0),2).*normpdf(z,0,b(end)),-5,5);
end

like = -sum(log(like_i));
end

