clear all;
clc;
load nlsw88
subset = (~isnan(wage).*~isnan(union).*~isnan(grade)==1);
y = wage(subset);
n = sum(subset);
X = [ones(length(wage),1) ttl_exp grade union];
X = X(subset,:);
[beta,~,~,~,stats]=regress(y,X);
true_parms = [-4.064795;.2635245;.608633;.9885717;];
regress_output = [beta;sqrt(stats(4))];
B = diag((y-X*beta).^2);
S = (X'*B*X);
Avar_GMM = inv((X'*X)/S*(X'*X));
se_GMM = sqrt(diag(Avar_GMM));

Avar_OLS_Hayashi = (X'*X)\S/(X'*X);
se_OLS_Hayashi = sqrt(diag(Avar_OLS_Hayashi));

se_Stata        = sqrt(diag(stats(4)*inv(X'*X)));
se_Stata_robust = sqrt(diag((X'*X)\(X'*B*X)/(X'*X)));

[beta true_parms];
[se_GMM se_OLS_Hayashi se_Stata se_Stata_robust];
% my_answers = b1;

results = cell(size(X,2)+1,4);
results{1,2}='Coefficient';
results{1,3}='Standard Error';
results{1,4}='t-stat';
varnames ={'Variable';'Intercept';'experience';'education';'union';};
for i=1:size(X,2)+1
    results{i,1}=varnames{i,1};
end
for i=1:size(X,2)
    results{i+1,2}=beta(i);
    results{i+1,3}=se_Stata_robust(i);
    results{i+1,4}=beta(i)/se_Stata_robust(i);
end

results
