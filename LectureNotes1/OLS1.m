function SSE = OLS1(beta,X,Y)
SSE = (Y-X*beta)'*(Y-X*beta);
end