function y = qstar(params,cost,cutoff)
% Compute OPTIMAL QUANTITY SUPPLIED given a cutoff and marginal cost draw.

% This function takes in a 5x1 vector of parameters (PARAMS) representing

%      The mean of the log-cost distribution
%      The variance of the log-cost distribution
%      The variance of the log-disturbance (assumed to be mean zero)
%      The fixed entry cost
%      The inverse demand intercept

% as well as a vector of costdraws and a scalar cutoff, which is the
% presumed entry threshold (enter if s<=cutoff) for competitors


mu = params(1);
var_c = params(2);
var_a = params(3);
f = params(4);
P = params(5);
var_s = var_a + var_c;


H = @(x) exp(mu + var_c*(log(x) - mu)/var_s + 0.5*var_c*var_a/var_s) .* exp(-0.5 * (log(x) - mu).^2 / var_s) .* x.^(-1) * (2 * pi * var_s)^(-.5);
% Expected cost given some signal (E[C|S=s]) times the distribution of
% signals, f(s) together form the integrand of E[C|S<=sstar]


ec_ent = @(sstar) quadgk(H,0,sstar)/logncdf(sstar,mu,sqrt(var_s));
% Integrate from s=0 to s=sstar.
% Since there is no closed form, use quadrature.

% Use quadgk, which works best in the limit.  For example, using quadl
% works up until sstar of about 10 million, whereas quadgk works up until
% sstar of about 1 quadrillion.  It also works for Inf, while quadl does
% not.


EQnot = @(sstar) logncdf(sstar,mu,sqrt(var_s))*(2*P/2 - ec_ent(sstar))/(1 + logncdf(sstar,mu,sqrt(var_s)));
% Expected quantity supplied by competitors is a function of sstar.

qstarf = @(c,sstar) max(0,P/2 - 0.5*EQnot(sstar) - c/2);
% My quantity, given that I've entered, as a function of my true cost


y = qstarf(cost,cutoff);


end
