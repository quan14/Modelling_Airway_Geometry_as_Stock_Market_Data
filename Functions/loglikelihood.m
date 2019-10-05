function [ loglike ] = loglikelihood(para,Xsegment,nameofdist)
%computes the loglikelihood for different distributions.
if strcmp(nameofdist,'normal')==1
%     makeloglike=makedist('normal',para(1),para(2));
%     loglike=sum(log(pdf(makeloglike,Xsegment)));
    loglike=sum(log(fastnormalpdf(Xsegment,para(1),para(2))));
elseif strcmp(nameofdist,'studentt')==1
    loglike=sum(log(faststudenttpdf(Xsegment,para(1),para(2),para(3))));
    %loglikelihood of location scale t distribution
end


end

