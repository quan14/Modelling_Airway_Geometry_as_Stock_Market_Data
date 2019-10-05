function [logprior] = logprior(para,hyper,nameofdist)
%Computes the log prior for different distrubtions
if strcmp(nameofdist,'normal')==1
    logprior=log(fastnormalpdf(para(1),hyper(1),sqrt((para(2)^2)/hyper(2))))+log(inversegampdf((para(2)^2),0.5*hyper(3),0.5*hyper(3)*hyper(4)^2));%log(chi2pdf(hyper(3)*(hyper(4)^2)/(para(2)^2),hyper(3))*((hyper(3)*hyper(4)^2)/((para(2)^2)^2)));
elseif strcmp(nameofdist,'studentt')==1
    logprior=log(fastnormalpdf(para(1),hyper(1),sqrt((para(2)^2)/hyper(2))))+log(inversegampdf((para(2)^2),0.5*hyper(3),0.5*hyper(3)*hyper(4)^2))+log(uniformpdf(para(3),hyper(5),hyper(6)));%log(chi2pdf(hyper(3)*(hyper(4)^2)/(para(2)^2),hyper(3))*((hyper(3)*hyper(4)^2)/((para(2)^2)^2)));
    %location scale t
end

end

