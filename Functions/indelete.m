function [newmodel,newcp,newclcp,newicp,newpara] = indelete(X,likelihood,numberofpara,N,kmax,inglcl,ijump,ijumpto,ijumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,mprior,hyper,model,cp,clcp,icp,para,series,jacobian,numofclusters,currentcluster);
%delete an independent changepoint from a series, series is the series
%selected
newmodel=model;
newcp=cp;
newclcp=clcp;
newicp=icp;
newpara=para;

%need possiblemoves for the reverse jump of adding in all possible places
%+1 for the deleted cp
possiblelocations=clcp+icp;
possiblemoves=find(possiblelocations==0);
%%%

if sum(clcp,2)==(model-1)
%     disp('Error in picking model to delete')
    return %leave function
end

if sum(icp,2)==0 %if no icp
    return
end

r=randi([1 sum(icp,2)]);
icptodelete=find(icp==1); %find possible icp that can be deleted
curcp=icptodelete(r); %pick one 
left=(find(cp(1,:)==curcp))-1; %find the left segment
% curcp=9;
% possiblemoves=[];
% left=1;

if strcmp(likelihood,'normal')
    q1=normrnd(0,adproposalsd(1));
    q2=normrnd(0,adproposalsd(2));
    
    newmu=mean(X(cp(left)+1:cp(left+2)))+q1;
    newsig=std(X(cp(left)+1:cp(left+2)))^2+q2;
    
    while newsig<=0;
        q2=normrnd(0,adproposalsd(2));
        newsig=std(X(cp(left)+1:cp(left+2)))^2+q2;
%         disp('Error in delete sig')
    end
    
    post=loglikelihood([newmu,sqrt(newsig)],X(cp(left)+1:cp(left+2)),likelihood)+logprior([newmu,sqrt(newsig)],hyper,likelihood);
    cur=loglikelihood(para(1,left,:),X(cp(left)+1:cp(left+1)),likelihood)+loglikelihood(para(1,left+1,:),X(cp(left+1)+1:cp(left+2)),likelihood)+logprior(para(1,left,:),hyper,likelihood)+logprior(para(1,left+1,:),hyper,likelihood)+log(rprior(series,curcp));
    
    newmodelprior=log(mprior(model-1));
    oldmodelprior=log(mprior(model));
    
    newtooldproposal=log(fastnormalpdf(para(1,left,1)-mean(X(cp(left)+1:cp(left+1))),0,adproposalsd(1)))+log(fastnormalpdf((para(1,left,2)^2)-std(X(cp(left)+1:cp(left+1)))^2,0,adproposalsd(2)))+log(fastnormalpdf(para(1,left+1,1)-mean(X(cp(left+1)+1:cp(left+2))),0,adproposalsd(1)))+log(fastnormalpdf((para(1,left+1,2)^2)-std(X(cp(left+1)+1:cp(left+2)))^2,0,adproposalsd(2)));
    oldtonewproposal=log(fastnormalpdf(q1,0,adproposalsd(1)))+log(fastnormalpdf(q2,0,adproposalsd(2)));
    
    jumpnum=log(ijumpfrom(3)*inglcl(1))+log(1/(size(possiblemoves,2)+1)); %prob of selecting add move and to add a cp at the same location 
    jumpden=log(ijumpto(4)*inglcl(1))+log(1/(sum(icp,2))); %prob of selecting delete move and selecting the right cp to delete
    
    v=rand(1);
    a=min(exp(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal+jumpnum-jumpden),1);
%     if model==4;
%         cell(Aw)
%     end
    if v<=a
        newmodel=model-1;
        newcp(left+1)=[];
        newicp(curcp)=0;
        if left==1
            if left+2<=model
                newpara=newpara(1,left+2:model,:);
            else
                newpara=[];
            end
            if isempty(newpara)
                newpara=cat(3,newmu,sqrt(newsig));
            else
                newpara=[cat(3,newmu,sqrt(newsig)),newpara];
            end
        elseif left==(model-1)
            newpara=newpara(1,1:left-1,:);
            newpara=[newpara,cat(3,newmu,sqrt(newsig))];
        else
            newpara=[newpara(1,1:left-1,:),cat(3,newmu,sqrt(newsig)),newpara(1,left+2:model,:)];
        end
    end
elseif strcmp(likelihood,'studentt')==1
    q1=normrnd(0,adproposalsd(1));
    q2=normrnd(0,adproposalsd(2));
    
    newmu=mean(X(cp(left)+1:cp(left+2)))+q1;
    newsig=std(X(cp(left)+1:cp(left+2)))^2+q2;
    newdf=(para(1,left,3)+para(1,left+1,3))/2;%deterministic proposal for df but because of this, we need a jacobian!!1
    %%%jacobian is calculated as 1/2
    while newsig<=0;
        q2=normrnd(0,adproposalsd(2));
        newsig=std(X(cp(left)+1:cp(left+2)))^2+q2;
%         disp('Error in delete sig')
    end
    
    post=loglikelihood([newmu,sqrt(newsig),newdf],X(cp(left)+1:cp(left+2)),likelihood)+logprior([newmu,sqrt(newsig),newdf],hyper,likelihood);
    cur=loglikelihood(para(1,left,:),X(cp(left)+1:cp(left+1)),likelihood)+loglikelihood(para(1,left+1,:),X(cp(left+1)+1:cp(left+2)),likelihood)+logprior(para(1,left,:),hyper,likelihood)+logprior(para(1,left+1,:),hyper,likelihood)+log(rprior(series,curcp));
    
    newmodelprior=log(mprior(model-1));
    oldmodelprior=log(mprior(model));
    
    newtooldproposal=log(fastnormalpdf(para(1,left,1)-mean(X(cp(left)+1:cp(left+1))),0,adproposalsd(1)))+log(fastnormalpdf((para(1,left,2)^2)-std(X(cp(left)+1:cp(left+1)))^2,0,adproposalsd(2)))+log(fastnormalpdf(para(1,left+1,1)-mean(X(cp(left+1)+1:cp(left+2))),0,adproposalsd(1)))+log(fastnormalpdf((para(1,left+1,2)^2)-std(X(cp(left+1)+1:cp(left+2)))^2,0,adproposalsd(2)))+log(fastnormalpdf(para(1,left,3)-newdf,0,adproposalsd(3)));
    oldtonewproposal=log(fastnormalpdf(q1,0,adproposalsd(1)))+log(fastnormalpdf(q2,0,adproposalsd(2)));
    
    jumpnum=log(ijumpfrom(3)*inglcl(1))+log(1/(size(possiblemoves,2)+1)); %prob of selecting add move and to add a cp at the same location 
    jumpden=log(ijumpto(4)*inglcl(1))+log(1/(sum(icp,2))); %prob of selecting delete move and selecting the right cp to delete
    
    logjacobian=log(1/jacobian); %reverse of the jacobian calculated in the forward jump
    
    v=rand(1);
    a=min(exp(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal+jumpnum-jumpden+logjacobian),1);
%     if model==4;
%         cell(Aw)
%     end
    if v<=a
        newmodel=model-1;
        newcp(left+1)=[];
        newicp(curcp)=0;
        if left==1
            if left+2<=model
                newpara=newpara(1,left+2:model,:);
            else
                newpara=[];
            end
            if isempty(newpara)
                newpara=cat(3,newmu,sqrt(newsig),newdf);
            else
                newpara=[cat(3,newmu,sqrt(newsig),newdf),newpara];
            end
        elseif left==(model-1)
            newpara=newpara(1,1:left-1,:);
            newpara=[newpara,cat(3,newmu,sqrt(newsig),newdf)];
        else
            newpara=[newpara(1,1:left-1,:),cat(3,newmu,sqrt(newsig),newdf),newpara(1,left+2:model,:)];
        end
    end
else
    disp('Need to update likelihood')
    disp('Error in delete')
end
    for j=1:newmodel(1)
        if newcp(1,j+1)-newcp(1,j)<minimumdistance
            cell(aw);
        end
    end
end

