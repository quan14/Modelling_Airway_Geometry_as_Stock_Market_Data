function [y] = faststudenttpdf(X,mu,sigma,df)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
y=(gamma((df+1)/2)/(sigma*sqrt(df*pi)*gamma(df/2)))*(((df+((X-mu)/sigma).^2)/df).^(-(df+1)/2));

end

