function [newmodel,newcp,newclcp,newicp,newpara] = inmove(X,likelihood,numberofpara,N,kmax,intogl,ijump,ijumpto,ijumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,mprior,hyper,model,cp,clcp,icp,para,series,jacobian,numofclusters,currentcluster);
%Move each changepoint in a series, series is the series selected
newmodel=model;
newcp=cp;
newclcp=clcp;
newicp=icp;
newpara=para;

%check if the model consists of just global changepoints
if sum(clcp,2)==(model-1) %if the number of global changepoint is equal to the number of cahngepoints in the current model i.e no independent changepoints
%     disp('Error in picking model to move')
    return %leave function
end

if sum(icp,2)==0 %no icp 
    return
end

%select a changepoint to move
r=randi([1 sum(icp,2)]); %pick an index
icptomove=find(icp==1); %find possible icp that can be moved
curcp=icptomove(r); %pick one 
left=(find(cp(1,:)==curcp))-1; %find the left segment
procp=curcp+((-1)^(randi([0 1])))*poissrnd(proposalmove,1); %propose a new cp

if procp<=1+minimumdistance||procp>=size(X,1)-1-minimumdistance %can add merge move here
    %dont really need to say its a error as this will get a -inf in the
    %prior for the new changepoint
%     disp('movecp error - out of bounds');
    return %leave function
end
if procp<=cp(1,left)+minimumdistance||procp>=cp(1,left+2)-minimumdistance
%     disp('movecp error - inside minimum distance');
    return %leave function
end

if strcmp(likelihood,'normal')==1    
    q1=normrnd(0,rmproposalsd(1));
    q2=normrnd(0,rmproposalsd(1));
    q3=normrnd(0,rmproposalsd(2));
    q4=normrnd(0,rmproposalsd(2));
    
    newleftmu=para(1,left,1)+q1;
    newrightmu=para(1,left+1,1)+q2;
    newleftsig=para(1,left,2)^2+q3;
    newrightsig=para(1,left+1,2)^2+q4;
    
    while newleftsig<=0||newrightsig<=0
        q3=normrnd(0,rmproposalsd(2));
        q4=normrnd(0,rmproposalsd(2));
        newleftsig=para(1,left,2)^2+q3;
        newrightsig=para(1,left+1,2)^2+q4;
%         disp('Error in move sigma in move cp')
    end
    
    post=loglikelihood([newleftmu,sqrt(newleftsig)],X(cp(left)+1:procp),likelihood)+loglikelihood([newrightmu,sqrt(newrightsig)],X(procp+1:cp(left+2)),likelihood)+logprior([newleftmu,sqrt(newleftsig)],hyper,likelihood)+logprior([newrightmu,sqrt(newrightsig)],hyper,likelihood)+log(rprior(series,procp));
    cur=loglikelihood(para(1,left,:),X(cp(left)+1:cp(left+1)),likelihood)+loglikelihood(para(1,left+1,:),X(cp(left+1)+1:cp(left+2)),likelihood)+logprior(para(1,left,:),hyper,likelihood)+logprior(para(1,left+1,:),hyper,likelihood)+log(rprior(series,curcp));

    v=rand(1);
    a=min(exp(post-cur),1); % no need for proposal distribution as symmetrical
    if v<=a
        newcp(1,left+1)=procp;
        newpara(1,left,:)=[newleftmu,sqrt(newleftsig)];
        newpara(1,left+1,:)=[newrightmu,sqrt(newrightsig)];
        newicp(curcp)=0;
        newicp(procp)=1;
    end
elseif strcmp(likelihood,'studentt')==1
    q1=normrnd(0,rmproposalsd(1));
    q2=normrnd(0,rmproposalsd(1));
    q3=normrnd(0,rmproposalsd(2));
    q4=normrnd(0,rmproposalsd(2));
    q5=normrnd(0,rmproposalsd(3));
    q6=normrnd(0,rmproposalsd(3));
    
    newleftmu=para(1,left,1)+q1;
    newrightmu=para(1,left+1,1)+q2;
    newleftsig=para(1,left,2)^2+q3;
    newrightsig=para(1,left+1,2)^2+q4;
    newleftdf=para(1,left,3)+q5;
    newrightdf=para(1,left+1,3)+q6;
    
    %might not need this, prior should take care of this
    while newleftsig<=0||newrightsig<=0
        q3=normrnd(0,rmproposalsd(2));
        q4=normrnd(0,rmproposalsd(2));
        newleftsig=para(1,left,2)^2+q3;
        newrightsig=para(1,left+1,2)^2+q4;
%         disp('Error in move sigma in move cp')
    end
    while newleftdf<2||newrightdf<2
        q5=normrnd(0,rmproposalsd(3));
        q6=normrnd(0,rmproposalsd(3));
        newleftdf=para(1,left,3)+q5;
        newrightdf=para(1,left+1,3)+q6;    
    end
    
    post=loglikelihood([newleftmu,sqrt(newleftsig),newleftdf],X(cp(left)+1:procp),likelihood)+loglikelihood([newrightmu,sqrt(newrightsig),newrightdf],X(procp+1:cp(left+2)),likelihood)+logprior([newleftmu,sqrt(newleftsig),newleftdf],hyper,likelihood)+logprior([newrightmu,sqrt(newrightsig),newrightdf],hyper,likelihood)+log(rprior(series,procp));
    cur=loglikelihood(para(1,left,:),X(cp(left)+1:cp(left+1)),likelihood)+loglikelihood(para(1,left+1,:),X(cp(left+1)+1:cp(left+2)),likelihood)+logprior(para(1,left,:),hyper,likelihood)+logprior(para(1,left+1,:),hyper,likelihood)+log(rprior(series,curcp));

    v=rand(1);
    a=min(exp(post-cur),1); % no need for proposal distribution as symmetrical
    if v<=a
        newcp(1,left+1)=procp;
        newpara(1,left,:)=[newleftmu,sqrt(newleftsig),newleftdf];
        newpara(1,left+1,:)=[newrightmu,sqrt(newrightsig),newrightdf];
        newicp(curcp)=0;
        newicp(procp)=1;
    end
else
    disp('Need to update Likelihood')
    disp('Error move')
end
    %%%check if we have an minimum distance error
    for j=1:newmodel(1)
        if newcp(1,j+1)-newcp(1,j)<minimumdistance
            cell(aw);
        end
    end

end

