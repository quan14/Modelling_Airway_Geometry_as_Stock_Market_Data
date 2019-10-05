function [newmodel,newcp,newgcp,newicp,newpara,dc,newi2gcount,newg2icount] = glRJMCMC(X,likelihood,numberofpara,N,kmax,inglcl,ijump,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,gprior,mprior,hyper,model,cp,gcp,icp,para,jacobian,i2gcount,g2icount,numofclusters,clcp,currentcluster)
%Independent RJMCMC, Only one series selected
%move, add, delete, globaltoindep,indeptoglobal
newi2gcount=i2gcount;
newg2icount=g2icount;
dc=0;
gjumpto=gjump;
gjumpfrom=gjump;
%alter ijumpto prob if in beginning or end model
if sum(find(model==1))>0; %if we are in bottom model for any series
    gjumpto(1)=0; %move global=0 
    gjumpto(3)=0; %delete global=0
    gjumpto(4)=0; %global to independent =0
    gjumpto=gjumpto/sum(gjumpto);
end %split these two as if we have kmax =2 and we have the situation model = 1, 2
if sum(find(model==kmax))>0; %if we are in top model for any series
    gjumpto(2)=0; %add =0
    gjumpto=gjumpto/sum(gjumpto);
end
if sum(gcp,2)==0 %if we have no gcp 
    gjumpto(1)=0; %move global=0 
    gjumpto(3)=0; %delete global=0
    gjumpto(4)=0; %globaltoindep =0
    gjumpto=gjumpto/sum(gjumpto);
end
%%
% %select a series if we want to go independent to global
% icptotal=sum(icp,2); %total of independent changepoints in each series
% icpcheck=find(icptotal>0); %find which series has independent changepoints
% if isempty(icpcheck) %if we have no independent changepoints at all
%     gjumpto(5)=0; %independent to global =0
%     gjumpto=gjumpto/sum(gjumpto);
% else
%%
v=randi([1 size(X,2)]); %select a series
if sum(icp(1,:,v))==0 %if we have no independent changepoints for series v
    gjumpto(5)=0; %independent to global =0
    gjumpto=gjumpto/sum(gjumpto);
end
if sum(find(model==kmax))>0 %if any series is in top model
    maxmodel=find(model==kmax); % find all series in max model
    if isempty(find(maxmodel==v, 1)) %if we have selected a series that is not in top model but we know that one series is in top model, then we set independent to global move as zero.
        gjumpto(5)=0; %independent to global =0
        gjumpto=gjumpto/sum(gjumpto);
        %%
%     newmodel=model;
%     newcp=cp;
%     newgcp=gcp;
%     newicp=icp;
%     newpara=para;
%     return
    end
    %%
% %     if we selected series like this, then we would not have 1/size(X,2)
% probability
% %     icpindex=randi([1 size(icpcheck,1)]); %pick a series
% %     v=icpcheck(icpindex); %select a series
% %     if sum(find(model==kmax))>0 %if any series is in top model
% %         maxmodel=find(model==kmax); % find all series in max model
% %         if isempty(find(maxmodel==v, 1)) %if we have selected a series that is not in top model but we know that one series is in top model, then we reject the whole move.
% %             newmodel=model;
% %             newcp=cp;
% %             newgcp=gcp;
% %             newicp=icp;
% %             newpara=para;
% %             return
% %         end
% %     end 
% end
end

%rare case if we set kmax to 2 and we get model = 1,2
if isnan(gjumpto)
    newmodel=model;
    newcp=cp;
    newgcp=gcp;
    newicp=icp;
    newpara=para;
    return
end
    
u=rand(1);
if u<gjumpto(1);%move
    [newmodel,newcp,newgcp,newicp,newpara]=glmove(X,likelihood,numberofpara,N,kmax,inglcl,ijump,gjumpto,gjumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,gprior,mprior,hyper,model,cp,gcp,icp,para,jacobian);
       dc=1;
elseif u>gjumpto(1) && u<sum(gjumpto(1:2)); %add changepoint
    if sum(model+1==kmax)>0 %if we add and reach top model for any series
        gjumpfrom(2)=0; %add =0
        gjumpfrom=gjumpfrom/sum(gjumpfrom);
    end
    [newmodel,newcp,newgcp,newicp,newpara]=gladd(X,likelihood,numberofpara,N,kmax,inglcl,ijump,gjumpto,gjumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,gprior,mprior,hyper,model,cp,gcp,icp,para,jacobian);
    dc=2;
elseif u>sum(gjumpto(1:2)) && u<sum(gjumpto(1:3)); %delete changepoint
    if sum(model-1==1)>0; %if we delete and reach bottom model
        gjumpfrom(1)=0; %move =0
        gjumpfrom(3)=0; %delete =0
        gjumpfrom(4)=0; % globaltoindep=0
        gjumpfrom=gjumpfrom/sum(gjumpfrom);
    end
    [newmodel,newcp,newgcp,newicp,newpara]=gldelete(X,likelihood,numberofpara,N,kmax,inglcl,ijump,gjumpto,gjumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,gprior,mprior,hyper,model,cp,gcp,icp,para,jacobian);
    dc=3;
elseif u>sum(gjumpto(1:3)) && u<sum(gjumpto(1:4)); %global to independent
    [newmodel,newcp,newgcp,newicp,newpara,newg2icount]=gltoindep(X,likelihood,numberofpara,N,kmax,inglcl,ijump,gjumpto,gjumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,gprior,mprior,hyper,model,cp,gcp,icp,para,jacobian,i2gcount,g2icount);
    dc=4;
else %independent to global
    [newmodel,newcp,newgcp,newicp,newpara,newi2gcount]=indeptogl(X,likelihood,numberofpara,N,kmax,inglcl,ijump,gjumpto,gjumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,gprior,mprior,hyper,model,cp,gcp,icp,para,v,jacobian,i2gcount,g2icount);
    dc=5;
end

end

