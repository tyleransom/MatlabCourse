% Joe Mazur
% Simple Selective Entry Game




mu=5
sigmac=1
sigmaa=0.5
sigmas=sigmaa + sigmac


ec_s = @(s) exp(mu + sigmac*(log(s) - mu)/sigmas + 0.5*sigmac*sigmaa/sigmas)
% Expected cost given some signal, E[C|S=s]

f_s = @(s) exp(-0.5 * (log(s) - mu).^2 / sigmas) .* s.^(-1) * (2 * pi * sigmas)^(-.5)
% Distribution of signals, f(s) 


H = @(x) exp(mu + sigmac*(log(x) - mu)/sigmas + 0.5*sigmac*sigmaa/sigmas) .* exp(-0.5 * (log(x) - mu).^2 / sigmas) .* x.^(-1) * (2 * pi * sigmas)^(-.5)
% Put them together to form the integrand of E[C|S<=sstar]

ec_ent = @(sstar) quadgk(H,0,sstar)/logncdf(sstar,mu,sigmas)
% Integrate from s=0 to s=sstar.
% Since there is no closed form, use quadrature.

% Use quadgk, which works best in the limit.  For example, using quadl
% works up until sstar of about 10 million, whereas quadgk works up until
% sstar of about 1 quadrillion.  It also works for Inf, while quadl does
% not.



EQnot = @(sstar) logncdf(sstar,mu,sigmas)*(2 - ec_ent(sstar))/(1 + logncdf(sstar,mu,sigmas))
% Expected quantity supplied by competitors is a function of sstar.

qstar = @(c,sstar) 1 - 0.5*EQnot(sstar) - c/2
% My quantity, given that I've entered, as a function of my true cost

instprofit = @(c,sstar) (2 - qstar(c,sstar) - EQnot(sstar))*qstar(c,sstar)
% Instantaneous profit

diffx = @(s,sstar,c) instprofit(c,sstar)*lognpdf(c,mu + sigmac/sigmas * (log(s) - mu),sigmac*sigmaa/sigmas)
% Instantaneous profit for a given c times conditional distribution of c
% given s, all varying with sstar.


% Below I've basically just given it an sstar and an s, and then defined
% the expected profit function to be the integral of instantaneous profit
% over the conditional distribution of c, which is set for a given s.

sstar = 100
s = 5

EQstar = EQnot(sstar)

profx = @(x) (2 - [1 - 0.5*EQstar - x./2] - EQstar).*[1 - 0.5*EQstar - x./2].*lognpdf(x,mu + sigmac/sigmas * (log(s) - mu),sigmac*sigmaa/sigmas)




