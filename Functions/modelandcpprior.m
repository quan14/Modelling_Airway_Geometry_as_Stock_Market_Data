function [rprior,gprior,mprior] = modelandcpprior(X,rpriordist,rhyper,rscale,gpriordist,ghyper,gscale,scalefactor,autoscale,locationr,locationg,mpriordist,mhyper,kmax)
%Obtain the prior for the model given the distribution mentioned
if strcmp(rpriordist,'uniform')==1
    %     rprior=(1/(size(X,1)-1))*ones(size(X,1)-1,1);
    rprior=rscale*(1/(size(X,1)-1))*ones(size(X,2),size(X,1)-1);
    
    %% My version - I'm assuming the poir is logistic function
    
    %I need to assume that X is 1D
%     rprior = Logistic_function_maker(length(X),...
%         length(X)*0.5,1,0.5);
    
    %% This is mine
    
    rprior = Return_normalise_dis(rprior);
    %       rprior(9)=1;
elseif strcmp(rpriordist,'single')==1
    rprior=zeros(size(X,1)-1,1);
    rprior(locationr)=1;
else
    disp('Error for changepoint prior')
end
if strcmp(gpriordist,'uniform')==1
    if autoscale==1;
        gprior=(1/(size(X,1)-1))^(scalefactor)*ones(size(X,2),size(X,1)-1);     %%%try a prior where gprior is proportional to the number of series
    else
        gprior=gscale*(1/(size(X,1)-1))*ones(size(X,2),size(X,1)-1);
    end
elseif strcmp(gpriordist,'single')==1
    gprior=zeros(size(X,1)-1,1);
    gprior(locationg)=1;
else
    disp('Error in Global changepoint prior');
end
if strcmp(mpriordist,'binomial')==1
    mprior=binopdf(0:kmax-1,size(X,1)-1,mhyper/(size(X,1)-1));
else
    disp('Error for model prior')
end

end

