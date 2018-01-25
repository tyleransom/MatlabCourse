% Joe Mazur
% Cournot Competition with Selective Entry
% Simplified Model with 3 Identical Firms

clear
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explanation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Three firms compete in quantities with noisy signals of their marginal
% costs.  Marginal cost, C, is lognormally distributed i.i.d.  Signal, S,
% is equal to C*A, where A is a lognormally distributed i.i.d. disturbance
% term that is independent of C.

% The parameters of the underlying normal distributions for the logs of
% C,A, and S are as follows:

mu_c   = 2;
var_c = 1;

mu_a   = 0;
var_a = 1;

mu_s   = mu_c;
var_s = var_a + var_c;

% The timeline is as follows:

% 1) Firms first observe the set of potential entrants, their
%    characteristics, and the features of the market.  In this case, firms
%    are identical up to their draws from the common cost and signal
%    distributions.  There is a fixed entry cost, f, common to all firms.

% 2) Firms receive their marginal cost signal, S.

% 3) Each firm decides whether or not to enter the market.

% 4) Each firm observes the true marginal cost draw, C, and chooses
%    quantity, q, still uninformed about other firms' entry or quantity
%    decisions, and given a known market demand curve.

% 5) Other firms' entry and quantity decisions are observed, and market
%    transactions take place.

% The parameters of the market are as follows:

f = 5;         %common fixed entry cost
P = 50;        %inverse demand intercept



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Consider firm i's decision at stage 4.  Given that i has entered, and
% that he now knows his true Ci, how will he choose q?  He will maximize
% expected profit E[q*(2-Q-Ci)].  Solving this through and using symmetry
% gives us optimal quantity for firm i as a function of i's cost draw and
% the expected quantity for competiting firms:

% q = P/2 - 0.5*Qn - c/2 ,
%            where
% Qn = (F(S_star) * (2*P/2 - E[Cj|Sj<S_star]))/(1 + F(S_star)) ,
%            and
% F(.) is the CDF of S
%
% Note that E[Cj|Sj<S_star] depends on the conditional distribution of Cj
% given Sj=S.  Consider momentarily the natural logs of C, S, and A (c, s,
% and a, respectively).  Given that A and C are independent, the
% distribution of s is the sum of two independent normals.  S is, however,
% not independent of C or A.  The correlation between c and s is rho:

rho = sqrt(var_c/var_s);

% The bivariate conditional distribution of c given s is
% Normal(mu_c + sqrt(var_c/var_s)*rho*(s - mu_s) , (1-rho^2)*var_c)

% We can use the mean and variance above as parameters in the corresponding
% lognormal distribution, and also to get E[Cj|Sj=S].

% The expression of interest, E[Cj|Sj<S_star], is the integral of
% E[Cj|Sj=S] * PDF(S) from 0 to S_star, divided by the probability
% that Sj < S_star.

% The integral in question has no closed form solution, so we'll use
% quadrature to numerically integrate.  Resulting values for
% E[Cj|Sj<S_star] will depend on the presumed S_star.  Moving forward,
% calculating expected profit conditional on a signal, S, will require the
% flexibility to change both S and S_star.  Moreover, expected profit will
% be the integral of instantaneous profit given some value of C, which
% depends crucially on E[Cj|Sj<S_star], times the conditional distribution
% of C given S, over the whole support of C (0,Inf).

% The XPROF function will calculate expected profit for a given S and a
% presumed S_star cutoff point of competing firms.

param = [mu_c,var_c,var_a,f,P]';
% Designate the vector of parameters to pass through to XPROF

S_star = [1,5]';
% Just an example of the S_star vector, which sets arbitrary values of the
% signal and the presumed sstar cutoff of competitors.
%      S_star(1) = signal, "s"
%      S_star(2) = presumed competitors' cutoff, "star"

xprof(param,S_star);
% Calculate expected profit for arbitrary values of S and S_Star.


fprintf('Model Parameters\n\n  Mean (log)MC:              %2.1f\n  Variance (log)MC:          %2.1f\n',mu_c,var_c)
fprintf('  Variance (log)Disurbance:  %2.1f\n',var_a)

fprintf('  Fixed Entry Cost:         %2.1f\n  Inverse Demand Intercept: %2.1f\n\n',f,P)



cutoff = fsolve(@(t) xprof(param,[t,t]'), 20)
% Solve for the cutoff value that returns zero expected profit such that
% the equilibrium condition that firms' entry decisions are identical
% holds.

% The syntax above parameterizes the XPROF function by creating an
% anonymous function @(t) that evaluates XPROF with PARAMS == param
% (defined above) and with s == star == t.  FSOLVE tries different values
% of t, starting at 0.01, to make XPROF == 0.

check = xprof(param,[cutoff, cutoff]')
% Double-check that expected profit is zero when the signal and presumed
% competitors' threshold value both equal cutoff.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal Consistency Check Using Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


costs = random('logn',mu_c,sqrt(var_c),10000,1);
signals = costs.*random('logn',0,sqrt(var_a),10000,1);

costmean_analytical = exp(mu_c + 0.5*var_c);
costmean_simulated = mean(costs);

signalmean_analytical = exp(mu_s+0.5*var_s);
signalmean_simulated = mean(signals);



entry = signals < cutoff;
nonentry = signals >= cutoff;
margentry1 = signals >= (cutoff - 3);
margentry2 = margentry1 - nonentry;
entryprob=mean(entry);

qstars = qstar(param,costs,cutoff);

% The following lines compute expected cost, given entry, for assumed
% values of the cutoff point.  See the function m-files qstar.m or xprof.m
% for more detail.
H = @(x) exp(mu_c + var_c*(log(x) - mu_s)/var_s + 0.5*var_c*var_a/var_s) .* exp(-0.5 * (log(x) - mu_s).^2 / var_s) .* x.^(-1) * (2 * pi * var_s)^(-.5);

ec_ent = @(sstar) quadgk(H,0,sstar)/logncdf(sstar,mu_s,sqrt(var_s));
ec_Nent = @(sstar) quadgk(H,sstar,Inf)/(1-logncdf(sstar,mu_s,sqrt(var_s)));

ECost_Entry = ec_ent(cutoff);
ECost_NoEntry = ec_Nent(cutoff);

AvgCost_Entry = mean(costs(find(entry)));
AvgCost_NoEntry = mean(costs(find(nonentry)));

EQnot = logncdf(cutoff,mu_s,sqrt(var_s))*(2*P/2 - ec_ent(cutoff))/(1 + logncdf(cutoff,mu_s,sqrt(var_s)));
% Expected quantity supplied by competitors is a function of the cutoff.



Eprofit_entry = xprof(param,[cutoff, cutoff]')





Eqstar_entry = P/2 - EQnot/2 - ECost_Entry/2;
Avgqstar = mean(qstars(find(entry)));



EQuantity = EQnot*3/2;
%Also = 3*logncdf(cutoff,mu_s,sqrt(var_s))*(P/2 - (1/2)*EQnot - (1/2)*ECost_Entry)

AvgQuantity = 3*mean(entry.*qstars);

EPrice = P - EQuantity;
AvgPrice = P - AvgQuantity;


fprintf('Expected marginal cost: %5.3f\n',costmean_analytical)
fprintf('Average simulated marginal cost: %5.3f\n\n',costmean_simulated)
fprintf('Expected marginal cost signal: %5.3f\n',signalmean_analytical)
fprintf('Average simulated marginal cost signal: %5.3f\n\n',signalmean_simulated)
fprintf('Probability of entry (analytical): %5.3f\n',logncdf(cutoff,mu_s,sqrt(var_s)))
fprintf('Probability of entry (sample): %5.3f\n\n',entryprob)

fprintf('Expected cost for entering firms: %5.3f\n',ECost_Entry)
fprintf('Expected cost for non-entering firms: %5.3f\n\n',ECost_NoEntry)

fprintf('Average simulated cost for entering firms: %5.3f\n',AvgCost_Entry)
fprintf('Average simulated cost for non-entering firms: %5.3f\n\n',AvgCost_NoEntry)

fprintf('Expected "qstar" for entering firms: %5.3f\n',Eqstar_entry)
fprintf('Average "qstar" for entering firms: %5.3f\n\n',Avgqstar)

%fprintf('Expected profit for entering firms: %5.3f\n',Eprofit_entry)
%fprintf('Average profit for entering firms: %5.3f\n\n',Avgprofit)

fprintf('Expected price and quantity (p,q): (%5.3f,%5.3f)\n',EPrice,EQuantity)
fprintf('Average simulated price and quantity (p,q): (%5.3f,%5.3f)\n',AvgPrice,AvgQuantity)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost Distributions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The bivariate conditional distribution of c given s is
% Normal(mu_c + sqrt(var_c/var_s)*rho*(log(s) - mu_s) , (1-rho^2)*var_c)

cond_dist = @(c,s) lognpdf(c,mu_c + sqrt(var_c/var_s)*rho*(log(s) - mu_s) , sqrt((1-rho^2)*var_c));
s_dist    = @(s)   lognpdf(s,mu_s, sqrt(var_s));

joint     = @(c,s) cond_dist(c,s).*s_dist(s);

F = @(c) quadgk(@(x) joint(c,x),0,cutoff)/logncdf(cutoff,mu_s, sqrt(var_s));
%Entrants

G = @(c) quadgk(@(x) joint(c,x),cutoff,Inf)/(1-logncdf(cutoff,mu_s, sqrt(var_s)));
%Non-Entrants

J = @(c) cond_dist(c,cutoff);
%Marginal Entrants


% Quick integral check

% test1=0
% for i=[0:0.1:1000];
%    test2 = test1 + 0.1*H(i);
%    test1=test2;
% end
% test2

%subplot(3,1,1)
%fplot(@(x) lognpdf(x,mu_c,sqrt(var_c)),[0 40])
%subplot(3,1,2)
%fplot(F,[0 40])
%subplot(3,1,3)
%fplot(G,[0 40])

Y = @(c) [lognpdf(c,mu_c,sqrt(var_c)),F(c),J(c),G(c)];
fplot(Y,[0 40])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEXT STEPS

% For next Thursday, have PROOF that the data works, and that the simple
% version of my program is working.

% Then, generalize the model to N firms, and generalize the market demand
% equation, replacing 2 with some parameter.

% If time permits, try GMM estimation using simulated values for P, Qi, and
% entry.