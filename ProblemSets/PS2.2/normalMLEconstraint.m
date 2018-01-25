function ll = normalMLEconstraint(theta,x,y)
% b = theta(1:end-1);
theta(2)   = .15*tanh(theta(2))+.05; 
theta(10)  = -(theta(10)).^2;
theta(end) = exp(theta(end));

ll = -sum(-.5*log(2*pi*(theta(end)^2))-.5*((y-x*theta(1:end-1))/theta(end)).^2);
end