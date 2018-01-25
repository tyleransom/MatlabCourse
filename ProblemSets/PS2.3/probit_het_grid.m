function like = probit_het_grid(b,X,y)
%PROBIT_HET_GRID Computes random effects for a binary probit
% [like] = probit_het(b,X,y)
% Estimates parameters b of a random-effects probit given data X and y
% Note that the std dev of the random effect is the last element of b
% Integration method is composite Simpson's Rule (by hand)
[N,K,T]=size(X);

% set up the grid matrices for integration
% grid should be NxTxG
step = .03;
grid = [-3:step:3]'*ones(1,T);
G = size(grid,1);
grid = grid';
grid = reshape(grid,[1 T G]);
grid = repmat(grid,[N 1 1]);

% set up the multiplication X*beta
X2 = permute(X, [1 3 2]);
Xb  = reshape(X2,N*T,K)*b(1:end-1);
Xb1 = reshape(Xb,[N T 1]);
Xb1 = repmat(Xb1,[1 1 G]);
y   = reshape(y,[N T 1]);
y   = repmat(y,[1 1 G]);

% form matrix of normal pdf weights (i.e. f(x) in an integrand)
% these weights matrices should be of size NxG
normpdf_weights = normpdf(squeeze(grid(:,1,:)),0,b(end));
% form matrix of composite simpson's rule weights
simpson_weights = ones(size(squeeze(grid(:,1,:))));
for g=1:G
    if mod(g,2)==0;
        simpson_weights(:,g)=2;
    else
        simpson_weights(:,g)=4;
    end
end
simpson_weights(:,1)=ones(N,1);
simpson_weights(:,end)=ones(N,1);
simpson_weights = (step/3)*simpson_weights;

if sum(simpson_weights.*normpdf_weights,2)-1>1e-3
    error('normdpdf doesn`t integrate to 1')
end

like = -sum(log(sum(squeeze(prod(normcdf(Xb1+grid).^((y==1)) .* (1-normcdf(Xb1+grid)).^(y==0),2)).*normpdf_weights.*simpson_weights,2)),1);
end

