function [y] = uniformpdf(X,a,b)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
if X<a || X>b
    y=0;
else
    y=1/(b-a);
end

end

