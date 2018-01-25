%%% OLStest.m (to be used with OLS1.m, a simple OLS function)
clear all;
clc;
% initialize the number of observations and the X matrix
n = 15000;
X = [ones(n,1) 17+2*randn(n,1)];
% Generate Y from Y = Xb + e, setting b=[2;1]:
Y = X*[2; 1] + randn(n,1);
% Do optimization
[beta1,SSE] = fminsearch('OLS1',.5*ones(2,1),[],X,Y);
[beta2,SSE] = fminunc('OLS1',.5*ones(2,1),[],X,Y);
[beta3,SSE] = fmincon('OLS1',.5*ones(size(X,2),1),[],[],[],[],-3*ones(size(X,2),1),.5*ones(size(X,2),1),[],[],X,Y);