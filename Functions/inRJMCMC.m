function [newmodel,newcp,newclcp,newicp,newpara] = inRJMCMC(X,likelihood,numberofpara,N,kmax,inglcl,ijump,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,mprior,hyper,model,cp,clcp,icp,para,v,jacobian,numofclusters,currentcluster)
%Independent RJMCMC, Only one series selected, v is the selected series
%resample, move, add, delete
ijumpto=ijump; %copy the jump probs (for the jump prob)
ijumpfrom=ijump; % copy the jump probs ( for the return probs)
%alter ijumpto prob if in beginning or end model
if model==1; %if we are in bottom model, only add or resample
    ijumpto(2)=0; %move =0 
    ijumpto(4)=0; %delete =0
    ijumpto=ijumpto/sum(ijumpto); %normalise
elseif model==kmax; %if we are in top model, only can resample, move, delete
    ijumpto(3)=0; %add =0
    ijumpto=ijumpto/sum(ijumpto); %normalise
end
u=rand(1); %uniform random number
% u=7/12;
if u<ijumpto(1);%resample
    [newmodel,newcp,newclcp,newicp,newpara]=inresample(X,likelihood,numberofpara,N,kmax,inglcl,ijump,ijumpto,ijumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,mprior,hyper,model,cp,clcp,icp,para,v,jacobian,numofclusters,currentcluster);
elseif u>ijumpto(1) && u<sum(ijumpto(1:2)); %move changepoint
    [newmodel,newcp,newclcp,newicp,newpara]=inmove(X,likelihood,numberofpara,N,kmax,inglcl,ijump,ijumpto,ijumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,mprior,hyper,model,cp,clcp,icp,para,v,jacobian,numofclusters,currentcluster);
elseif u>sum(ijumpto(1:2)) && u<sum(ijumpto(1:3)); %add changepoint
    if model+1==kmax %if we add and reach top model
        ijumpfrom(3)=0; %add =0
        ijumpfrom=ijumpfrom/sum(ijumpfrom);
    end
    [newmodel,newcp,newclcp,newicp,newpara]=inadd(X,likelihood,numberofpara,N,kmax,inglcl,ijump,ijumpto,ijumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,mprior,hyper,model,cp,clcp,icp,para,v,jacobian,numofclusters,currentcluster);
else %delete
    if model-1==1; %if we delete and reach bottom model
        ijumpfrom(2)=0; %move =0
        ijumpfrom(4)=0; %delete =0
        ijumpfrom=ijumpfrom/sum(ijumpfrom);
    end
    [newmodel,newcp,newclcp,newicp,newpara]=indelete(X,likelihood,numberofpara,N,kmax,inglcl,ijump,ijumpto,ijumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,mprior,hyper,model,cp,clcp,icp,para,v,jacobian,numofclusters,currentcluster);
end

end

