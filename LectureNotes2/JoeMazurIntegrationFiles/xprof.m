function xprof = argx(params,S_star)
% Compute EXPECTED PROFIT given a signal, s, and entry decision, star.

% This function takes in a 5x1 vector of parameters (PARAMS) representing

%      The mean of the log-cost distribution
%      The variance of the log-cost distribution
%      The variance of the log-disturbance (assumed to be mean zero)
%      The fixed entry cost
%      The inverse demand intercept

% as well as a 2x1 vector of signal parameters (S) representing

%      The signal received prior to making the entry decision
%      The presumed entry threshold (enter if s<=star) for competitors


mu = params(1);
var_c = params(2);
var_a = params(3);
f = params(4);
P = params(5);
var_s = var_a + var_c;


s = S_star(1);
star = S_star(2);



ec_s = @(x) exp(mu + var_c*(log(x) - mu)/var_s + 0.5*var_c*var_a/var_s);
% Expected cost given some signal, E[C|S=s]

f_s = @(x) exp(-0.5 * (log(x) - mu).^2 / var_s) .* x.^(-1) * (2 * pi * var_s)^(-.5);
% Distribution of signals, f(s) 


H = @(x) exp(mu + var_c*(log(x) - mu)/var_s + 0.5*var_c*var_a/var_s) .* exp(-0.5 * (log(x) - mu).^2 / var_s) .* x.^(-1) * (2 * pi * var_s)^(-.5);
% Put them together to form the integrand of E[C|S<=sstar]

ec_ent = @(sstar) quadgk(H,0,sstar)/logncdf(sstar,mu,sqrt(var_s));
% Integrate from s=0 to s=sstar.
% Since there is no closed form, use quadrature.

% Use quadgk, which works best in the limit.  For example, using quadl
% works up until sstar of about 10 million, whereas quadgk works up until
% sstar of about 1 quadrillion.  It also works for Inf, while quadl does
% not.



EQnot = @(sstar) logncdf(sstar,mu,sqrt(var_s))*(2*P/2 - ec_ent(sstar))/(1 + logncdf(sstar,mu,sqrt(var_s)));
% Expected quantity supplied by competitors is a function of sstar.

qstar = @(c,sstar) max(0,P/2 - 0.5*EQnot(sstar) - c/2);
% My quantity, given that I've entered, as a function of my true cost

instprofit = @(c,sstar) (P - qstar(c,sstar) - EQnot(sstar))*qstar(c,sstar);
% Instantaneous profit

diffx = @(s,sstar,c) instprofit(c,sstar)*lognpdf(c,mu + var_c/var_s * (log(s) - mu),sqrt(var_c*var_a/var_s));
% Instantaneous profit for a given c times conditional distribution of c
% given s, all varying with sstar.  The function above is the shorthand
% version of what will be PROFX below.


EQstar = EQnot(star);

profx = @(x) (P - [P/2 - 0.5*EQstar - x./2] - EQstar).*[P/2 - 0.5*EQstar - x./2].*lognpdf(x,mu + var_c/var_s * (log(s) - mu),sqrt(var_c*var_a/var_s));
% Long form for instantaneous profit [(P - qstar - EQnot)*qstar] times density

xprof = quadgk(profx,0,Inf) - f;



end
