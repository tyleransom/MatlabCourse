clear all;
clc;
load nlsy97
options = optimset('Disp','iter-detailed','MaxFunEvals',1e12,'MaxIter',1e6);
options2 = optimset('Disp','iter-detailed','MaxFunEvals',1e12,'MaxIter',1e6,'GradObj','on','DerivativeCheck','on','LargeScale','off');
options3 = optimset('Disp','iter-detailed','MaxFunEvals',1e12,'MaxIter',1e6,'GradObj','on','DerivativeCheck','on','LargeScale','on','Hessian','on');
rand('seed',1234); randn('seed',1234);
%% Problem 1(a) -- Reshape the data into long format
N = length(ID);
T = 5;
activityt = activity';
log_waget = log_wage';
hgct = hgc';
expert = exper';
Diplomat = Diploma';
AAt = AA';
BAt = BA';
malet = repmat(male,1,T)';
AFQTt = repmat(AFQT,1,T)';
Mhgct = repmat(Mhgc,1,T)';
IDt = repmat(ID,1,T)';
IDl = IDt(:);
IDprime  = [1:numel(ID)]';
IDprimet = repmat(IDprime,1,T)';
IDprimel = IDprimet(:);
X = [ones(N*T,1) malet(:) AFQTt(:) Mhgct(:) hgct(:) expert(:) Diplomat(:) AAt(:) BAt(:)];
y = activityt(:);
K = size(X,2);
%% Problem 1(b) -- estimate mlogit using fminunc
bstart = mnrfit(X(:,2:end),y); % I cheated here and used the mnrfit starting values (for computational speed)
% bstart = rand(2*size(X,2),1);
[bunc,lunc,~,~,gunc,hunc] = fminunc('mlogit',bstart(:),options,X,y);
SEunc = sqrt(diag(inv(hunc)));
%% add gradient
bstart = rand(2*size(X,2),1);
[buncg,luncg,~,~,guncg,huncg] = fminunc('mlogit_gradient',bstart(:),options2,X,y);
SEuncg = sqrt(diag(inv(huncg)));
[bunch,lunch,~,~,gunch,hunch] = fminunc('mlogit_gradient',bstart(:),options3,X,y);
SEunch = sqrt(diag(inv(hunch)));
%% binomial hessian test
% bstart = rand(size(X,2),1);
% [buncb,luncb,~,~,guncb,huncb] = fminunc('logit_gradient',bstart(:),options2,X,2-(y==1|y==2));