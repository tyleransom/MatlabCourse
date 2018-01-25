clear all;
clc;
load carbig
format short
y  = Acceleration;
X  = [ones(size(Acceleration)) Cylinders Displacement Weight (Model_Year<76)];
b  = X\y;
se = sqrt(diag(((y-X*b)'*(y-X*b))/(size(X,1)-size(X,2))*((X'*X)\eye(size(X,2)))));
t  = b./se;
p  = 2*(1-tcdf(abs(t),size(X,1)-size(X,2)));
%%%%%%%
estimates  = num2cell([b se t p]);
estimates1 = [b se t p];
varnames   = cellstr(char('Intercept','Cylinders','Displacement','Weight','Older Model'));
headers    = cellstr(char('Variable','Coefficient','Standard Error','T-stat','p-value'))';
result     = cat(1,headers,cat(2,varnames,estimates));
%%%%%%%
xlswrite('estimation_results.xlsx',result,1);
matrix2latex(estimates1, 'estimation_results.tex','rowLabels', varnames, 'columnLabels', headers, 'alignment', 'c','format', '%7.3f');
