clear all;
clc;
load nlsy97m
options = optimset('Disp','iter-detailed','MaxFunEvals',1e12,'MaxIter',1e6);
options2 = optimset('Disp','iter-detailed','MaxFunEvals',1e12,'MaxIter',1e6,'GradObj','on','DerivativeCheck','on');
options3 = optimset('Disp','iter-detailed','MaxFunEvals',1e12,'MaxIter',1e6,'GradObj','on','Hessian','on','DerivativeCheck','on');
rand('seed',1234); randn('seed',1234);
%% Reshape the data into long format
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
y = activityt(:);
log_wagel = log_waget(:);
X = [ones(N*T,1) malet(:) AFQTt(:) Mhgct(:) hgct(:) expert(:) Diplomat(:) AAt(:) BAt(:)];
K = size(X,2);
%% Problem 2(a)
tabulate(activityt(:));
%% Problem 2(b) -- estimate mlogit using fminunc
bstart = rand(5*size(X,2),1);
tic;
[bng,lng,~,~,~,hng] = fminunc('mlogit_gradient',bstart(:),options,X,y);
SEbng = sqrt(diag(inv(hng)));
disp('Time spent for mlogit MLE (no gradient)')
toc;
resultng = [bng SEbng]
%% Problem 2(c) -- estimate mlogit using fminunc (and provide gradient)
tic;
[bg,lg,~,~,~,hg] = fminunc('mlogit_gradient',bstart(:),options2,X,y);
hg = full(hg);
SEbg = sqrt(diag(inv(hg)));
disp('Time spent for mlogit MLE (user-supplied gradient)')
toc;
resultg = [bg SEbg]
%% Problem 3(a) -- log wage normal MLE
subset = ~isnan(log_waget(:));
X = [ones(N*T,1) malet(:) AFQTt(:) Mhgct(:) hgct(:) expert(:) Diplomat(:) AAt(:) BAt(:) y==2 y==3];
X1 = X(subset,:);
y1 = log_wagel(subset);
bols = (X1'*X1)\(X1'*y1);
[b,l,~,~,~,h]=fminunc('normalMLEconstraint',-.45*rand(length(bols)+1,1),options,X1,y1);
b(2)   = .15*tanh(b(2))+.05; 
b(10)  = -(b(10)).^2;
b(end) = exp(b(end));

A = ones(length(bols)+1,1);
A(2)   = .15*(sech(b(2)).^2);
A(10)  = -2*(b(10));
A(end) = exp(b(end));
A = diag(A);

SE = sqrt(diag(A'/(h)*A));

resultMLE = [b SE]
%% Problem 4(a) -- log wage normal MLE with gradient
b0 = rand(length(bols)+1,1);
[bmleng,lmleng,~,~,~,hmleng]=fminunc('normalMLEgradient',b0,options,X1,y1);
SEmleng = sqrt(diag(inv(hmleng)));
resultMLEng = [bmleng SEmleng]
[bmleg,lmleg,~,~,~,hmleg]   =fminunc('normalMLEgradient',b0,options2,X1,y1);
hmleg = full(hmleg);
SEmleg = sqrt(diag(inv(hmleg)));
resultMLEg = [bmleg SEmleg]
%% Problem 4(b) -- log wage normal MLE with hessian
[bmlenh,lmlenh,~,~,~,hmlenh]=fminunc('normalMLEhessian',b0,options,X1,y1);
SEmlenh = sqrt(diag(inv(hmlenh)));
resultMLEnh = [bmlenh SEmlenh]
[bmleh,lmleh,~,~,~,hmleh]   =fminunc('normalMLEhessian',b0,options3,X1,y1);
SEmleh = sqrt(diag(inv(hmleh)));
resultMLEh = [bmleh SEmleh]
save PS2results resultg lg resultng lng resultMLE l resultMLEng lmleng resultMLEg lmleg resultMLEnh lmlenh resultMLEh lmleh