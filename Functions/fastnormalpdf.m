function [y] = fastnormalpdf(X,mu,sigma)
%evaluate X as a normal with para(1) mu and para(2) sigma
y=(1 / (sqrt(2*pi)*sigma)) .* exp(-.5*((X-mu)/sigma).^2);

end

