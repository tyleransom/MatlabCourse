function like = logit_het_disc(b,X,y)
%LOGIT_HET_GRID Computes discrete random effects for a binary logit
% [like] = logit_het(b,X,y)
% Estimates parameters b of a random-effects logit given data X and y
% Note that the std dev of the random effect is the last element of b
% Heterogeneity has finite number of types
[N,K,T]=size(X);
beta  = b(1:end-1);
pi    = b(end);
% alpha = b(end-1);
alpha = 1;

% set up the multiplication X*beta 
X2 = permute(X, [1 3 2]);
Xb  = reshape(X2,N*T,K)*beta;
Xb1 = reshape(Xb,N,T);


%%% test code
% test = beta'*squeeze(X(1,:,:));
% First row of Xb1 should be the same as test
%%%

like_1 = prod(sigmoid(Xb1+alpha).^(y==1) .* (1-sigmoid(Xb1+alpha)).^(y==0),2);
like_2 = prod(sigmoid(Xb1+0    ).^(y==1) .* (1-sigmoid(Xb1+0    )).^(y==0),2);

like = -sum(log(pi*like_1+(1-pi)*like_2));
end

