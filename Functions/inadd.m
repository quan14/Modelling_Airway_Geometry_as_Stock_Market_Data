function [newmodel,newcp,newclcp,newicp,newpara] = inadd(X,likelihood,numberofpara,N,kmax,inglcl,ijump,ijumpto,ijumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,mprior,hyper,model,cp,clcp,icp,para,series,jacobian,numofclusters,currentcluster);
%Add changepoint to current mode for one series,series is the series
%selected
newmodel=model;
newcp=cp;
newclcp=clcp;
newicp=icp;
newpara=para;

possiblelocations=clcp+icp; %find where all the cps are
if sum(possiblelocations==0)==0
%     disp('Error is picking a free place to add cp')
    return;
end

% possiblemoves=randi([1 size(X,1)-1]);
possiblemoves=find(possiblelocations==0); %find all the empty locations
procp=possiblemoves(randi([1 size(possiblemoves,2)])); %pick a place to put a cp
% procp=160;
% possiblemoves=1;
% procp=9;

left=find(cp(1,:)<procp,1,'last'); %find left of cp
if procp-cp(left)<=minimumdistance||cp(left+1)-procp<=minimumdistance % if inside minimum distance
%     disp('Error in minimum distance of adding cp')
    return
end

if strcmp(likelihood,'normal')==1
    q1=normrnd(0,adproposalsd(1));
    q2=normrnd(0,adproposalsd(1));
    q3=normrnd(0,adproposalsd(2));
    q4=normrnd(0,adproposalsd(2));
    
    newleftmu=mean(X(cp(left)+1:procp))+q1;
    newrightmu=mean(X(procp+1:cp(left+1)))+q2;
    newleftsig=std(X(cp(left)+1:procp))^2+q3;
    newrightsig=std(X(procp+1:cp(left+1)))^2+q4;
    
    while newleftsig<=0||newrightsig<=0
        q3=normrnd(0,adproposalsd(2));
        q4=normrnd(0,adproposalsd(2));
        newleftsig=std(X(cp(left)+1:procp))^2+q3;
        newrightsig=std(X(procp+1:cp(left+1)))^2+q4;
%         disp('Error in sigma in add cp')
    end
    
    post=loglikelihood([newleftmu,sqrt(newleftsig)],X(cp(left)+1:procp),likelihood)+loglikelihood([newrightmu,sqrt(newrightsig)],X(procp+1:cp(left+1)),likelihood)+logprior([newleftmu,sqrt(newleftsig)],hyper,likelihood)+logprior([newrightmu,sqrt(newrightsig)],hyper,likelihood)+log(rprior(series,procp));
    cur=loglikelihood(para(1,left,:),X(cp(left)+1:cp(left+1)),likelihood)+logprior(para(1,left,:),hyper,likelihood);
    
    newmodelprior=log(mprior(model+1));
    oldmodelprior=log(mprior(model));
    
    newtooldproposal=log(fastnormalpdf(para(1,left,1)-mean(X(cp(left)+1:cp(left+1))),0,adproposalsd(1)))+log(fastnormalpdf((para(1,left,2)^2)-std(X(cp(left)+1:cp(left+1)))^2,0,adproposalsd(2)));
    oldtonewproposal=log(fastnormalpdf(q1,0,adproposalsd(1)))+log(fastnormalpdf(q2,0,adproposalsd(1)))+log(fastnormalpdf(q3,0,adproposalsd(2)))+log(fastnormalpdf(q4,0,adproposalsd(2)));
    
    jumpnum=log(ijumpfrom(4)*inglcl(1))+log(1/(sum(icp,2)+1));%prob of selection delete move and selection right cp to delete  
    jumpden=log(ijumpto(3)*inglcl(1))+log(1/(size(possiblemoves,2)));%prob of selection add move and the prob of picking where the cp goes
    
    v=rand(1);
    a=min(exp(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal+jumpnum-jumpden),1);
%     cell(aw)
    %     if model==2
%         cell(aw)
%     end
    if v<=a
        newmodel=model+1;
        newcp=sort([cp,procp]);
        newicp(procp)=1;
        if left==1
            if left+1>model
                newpara=[];
            else
                newpara=newpara(1,left+1:model,:);
            end
            if isempty(newpara)
                newpara=[cat(3,newleftmu,sqrt(newleftsig)),cat(3,newrightmu,sqrt(newrightsig))];
            else
                newpara=[cat(3,newleftmu,sqrt(newleftsig)),cat(3,newrightmu,sqrt(newrightsig)),newpara];
            end
        elseif left==(model)
            newpara=newpara(1,1:left-1,:);
            newpara=[newpara,cat(3,newleftmu,sqrt(newleftsig)),cat(3,newrightmu,sqrt(newrightsig))];
        else
            newpara=[newpara(1,1:left-1,:),cat(3,newleftmu,sqrt(newleftsig)),cat(3,newrightmu,sqrt(newrightsig)),newpara(1,left:model-1,:)];
        end
    end
elseif strcmp(likelihood,'studentt')==1
    q1=normrnd(0,adproposalsd(1));
    q2=normrnd(0,adproposalsd(1));
    q3=normrnd(0,adproposalsd(2));
    q4=normrnd(0,adproposalsd(2));
    q5=normrnd(0,adproposalsd(3));%%%we use a deterministic reverse jump, which causes a jacobian to arise
    %%% the jacobian is calculated to be 2
    
    newleftmu=mean(X(cp(left)+1:procp))+q1;
    newrightmu=mean(X(procp+1:cp(left+1)))+q2;
    newleftsig=std(X(cp(left)+1:procp))^2+q3;
    newrightsig=std(X(procp+1:cp(left+1)))^2+q4;
    newleftdf=para(1,left,3)+q5;
    newrightdf=para(1,left,3)-q5;
    
    %might not need this as prior should take care of this
    while newleftsig<=0||newrightsig<=0
        q3=normrnd(0,adproposalsd(2));
        q4=normrnd(0,adproposalsd(2));
        newleftsig=std(X(cp(left)+1:procp))^2+q3;
        newrightsig=std(X(procp+1:cp(left+1)))^2+q4;
%         disp('Error in sigma in add cp')
    end
    while newleftdf<2||newrightdf<2
        q5=normrnd(0,adproposalsd(3));
        newleftdf=para(1,left,3)+q5;
        newrightdf=para(1,left,3)-q5;    
    end
    
    post=loglikelihood([newleftmu,sqrt(newleftsig),newleftdf],X(cp(left)+1:procp),likelihood)+loglikelihood([newrightmu,sqrt(newrightsig),newrightdf],X(procp+1:cp(left+1)),likelihood)+logprior([newleftmu,sqrt(newleftsig),newleftdf],hyper,likelihood)+logprior([newrightmu,sqrt(newrightsig),newrightdf],hyper,likelihood)+log(rprior(series,procp));
    cur=loglikelihood(para(1,left,:),X(cp(left)+1:cp(left+1)),likelihood)+logprior(para(1,left,:),hyper,likelihood);
    
    newmodelprior=log(mprior(model+1));
    oldmodelprior=log(mprior(model));
    
    newtooldproposal=log(fastnormalpdf(para(1,left,1)-mean(X(cp(left)+1:cp(left+1))),0,adproposalsd(1)))+log(fastnormalpdf((para(1,left,2)^2)-std(X(cp(left)+1:cp(left+1)))^2,0,adproposalsd(2))); %deterministic reverse jump fro degress of freedom
    oldtonewproposal=log(fastnormalpdf(q1,0,adproposalsd(1)))+log(fastnormalpdf(q2,0,adproposalsd(1)))+log(fastnormalpdf(q3,0,adproposalsd(2)))+log(fastnormalpdf(q4,0,adproposalsd(2)))+log(fastnormalpdf(q5,0,adproposalsd(3)));
    
    jumpnum=log(ijumpfrom(4)*inglcl(1))+log(1/(sum(icp,2)+1));%prob of selection delete move and selection right cp to delete  
    jumpden=log(ijumpto(3)*inglcl(1))+log(1/(size(possiblemoves,2)));%prob of selection add move and the prob of picking where the cp goes
    
    logjacobian=log(jacobian);%%% calculated from the deterministic reverse jump
    
    v=rand(1);
    a=min(exp(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal+jumpnum-jumpden+logjacobian),1);
%     cell(aw)
    %     if model==2
%         cell(aw)
%     end
    if v<=a
        newmodel=model+1;
        newcp=sort([cp,procp]);
        newicp(procp)=1;
        if left==1
            if left+1>model
                newpara=[];
            else
                newpara=newpara(1,left+1:model,:);
            end
            if isempty(newpara)
                newpara=[cat(3,newleftmu,sqrt(newleftsig),newleftdf),cat(3,newrightmu,sqrt(newrightsig),newrightdf)];
            else
                newpara=[cat(3,newleftmu,sqrt(newleftsig),newleftdf),cat(3,newrightmu,sqrt(newrightsig),newrightdf),newpara];
            end
        elseif left==(model)
            newpara=newpara(1,1:left-1,:);
            newpara=[newpara,cat(3,newleftmu,sqrt(newleftsig),newleftdf),cat(3,newrightmu,sqrt(newrightsig),newrightdf)];
        else
            newpara=[newpara(1,1:left-1,:),cat(3,newleftmu,sqrt(newleftsig),newleftdf),cat(3,newrightmu,sqrt(newrightsig),newrightdf),newpara(1,left:model-1,:)];
        end
    end
else
    disp('Need to update likelihood')
    disp('Error add')

end
end

