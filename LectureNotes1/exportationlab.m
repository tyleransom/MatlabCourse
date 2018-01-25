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
% headers    = cellstr(char('Variable','Coefficient','Standard Error','T-stat','p-value'))';
headers    = cellstr(char('Variable','Coef.','Std Err','T-stat','p-value'))';
result     = cat(1,headers,cat(2,varnames,estimates));
%%%%%%%
xlswrite('estimation_results.xlsx',result,1);
matrix2latex(estimates1, 'estimation_results.tex','rowLabels', varnames, 'columnLabels', headers, 'alignment', 'c','format', '%7.3f');

% filename = 'tester.dat';
% fid = fopen(filename, 'w');
% fprintf(fid, '%s %s %s %s %s\n', result{1,:});
% for row=2:size(result,1)
%     fprintf(fid, '%s %5.4f %5.4f %5.4f %5.4f\n', result{row,:});
% end
% fclose(fid);

str1=fprintf('%14s %7s %7s %7s %7s\n', result{1,1},result{1,2},result{1,3},result{1,4},result{1,5});
for row=2:size(result,1)
    strk=fprintf('%14s %+06.4f %+06.4f %+06.4f %+06.4f\n', result{row,1},result{row,2},result{row,3},result{row,4},result{row,5});
end