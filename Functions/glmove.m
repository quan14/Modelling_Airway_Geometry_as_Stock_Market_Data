function [newmodel,newcp,newgcp,newicp,newpara] = glmove(X,likelihood,numberofpara,N,kmax,intogl,ijump,gjumpto,gjumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,gprior,mprior,hyper,model,cp,gcp,icp,para,jacobian);
%Move each changepoint in a series
newmodel=model;
newcp=cp;
newgcp=gcp;
newicp=icp;
newpara=para;

%check if the model consists of just global changepoints
% if sum(gcp,2)==(model-1) %if the number of global changepoint is equal to the number of cahngepoints in the current model i.e no independent changepoints
% %     disp('Error in picking model to move')
%     return %leave function
% end
%select a changepoint to move
r=randi([1 sum(gcp,2)]);
gcptomove=find(gcp==1); %find possible gcp that can be moved
curcp=gcptomove(r); %pick one 
left=zeros(1,size(X,2));
for j=1:size(X,2)
    left(j)=(find(cp{1,j}==curcp))-1; %find the left segment
end
procp=curcp+((-1)^(randi([0 1])))*poissrnd(proposalmove,1);

if procp<1+minimumdistance||procp>size(X,1)-1-minimumdistance 
    %dont really need to say its a error as this will get a -inf in the
    %prior for the new changepoint
%     disp('movecp error - out of bounds');
    return %leave function
end

for j=1:size(X,2)
    if procp<cp{1,j}(1,left(j))+1+minimumdistance||procp>cp{1,j}(1,left(j)+2)-1-minimumdistance %can add merge move here or flip
    %     disp('movecp error - inside minimum distance');
        return %leave function
    end
end

if strcmp(likelihood,'normal') 
    q1=zeros(1,size(X,2));
    q2=zeros(1,size(X,2));
    q3=zeros(1,size(X,2));
    q4=zeros(1,size(X,2));
    newleftmu=zeros(1,size(X,2));
    newrightmu=zeros(1,size(X,2));
    newleftsig=zeros(1,size(X,2));
    newrightsig=zeros(1,size(X,2));
    post=zeros(1,size(X,2));
    cur=zeros(1,size(X,2));
    for j=1:size(X,2)
        q1(j)=normrnd(0,rmproposalsd(1));
        q2(j)=normrnd(0,rmproposalsd(1));
        q3(j)=normrnd(0,rmproposalsd(2));
        q4(j)=normrnd(0,rmproposalsd(2));

        newleftmu(j)=para{1,j}(1,left(j),1)+q1(j);
        newrightmu(j)=para{1,j}(1,left(j)+1,1)+q2(j);
        newleftsig(j)=para{1,j}(1,left(j),2)^2+q3(j);
        newrightsig(j)=para{1,j}(1,left(j)+1,2)^2+q4(j);

        while newleftsig(j)<=0||newrightsig(j)<=0
            q3(j)=normrnd(0,rmproposalsd(2));
            q4(j)=normrnd(0,rmproposalsd(2));
            newleftsig(j)=para{1,j}(1,left(j),2)^2+q3(j);
            newrightsig(j)=para{1,j}(1,left(j)+1,2)^2+q4(j);
%             disp('Error in move sigma in move cp')
        end

        post(j)=loglikelihood([newleftmu(j),sqrt(newleftsig(j))],X(cp{1,j}(left(j))+1:procp,j),likelihood)+loglikelihood([newrightmu(j),sqrt(newrightsig(j))],X(procp+1:cp{1,j}(left(j)+2),j),likelihood)+logprior([newleftmu(j),sqrt(newleftsig(j))],hyper,likelihood)+logprior([newrightmu(j),sqrt(newrightsig(j))],hyper,likelihood);
        cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+loglikelihood(para{1,j}(1,left(j)+1,:),X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood)+logprior(para{1,j}(1,left(j)+1,:),hyper,likelihood);
    end
    v=rand(1);
    a=min(exp(sum(post-cur)+log(gprior(procp))-log(gprior(curcp))),1);
    %%%check if proposal in -inf
    if isinf(sum(post-cur)+log(gprior(procp))-log(gprior(curcp)))
        if isinf(sum(post))
%            disp('move')
        else
        cell(aw)
        end
    end
    if v<=a
        for j=1:size(X,2)
            newcp{1,j}(1,left(j)+1)=procp;
            newpara{1,j}(1,left(j),:)=[newleftmu(j),sqrt(newleftsig(j))];
            newpara{1,j}(1,left(j)+1,:)=[newrightmu(j),sqrt(newrightsig(j))];
        end
        newgcp(curcp)=0;
        newgcp(procp)=1;
    end
elseif strcmp(likelihood,'studentt')==1
    q1=zeros(1,size(X,2));
    q2=zeros(1,size(X,2));
    q3=zeros(1,size(X,2));
    q4=zeros(1,size(X,2));
    q5=zeros(1,size(X,2));
    q6=zeros(1,size(X,2));
    
    newleftmu=zeros(1,size(X,2));
    newrightmu=zeros(1,size(X,2));
    newleftsig=zeros(1,size(X,2));
    newrightsig=zeros(1,size(X,2));
    newleftdf=zeros(1,size(X,2));
    newrightdf=zeros(1,size(X,2));
    
    post=zeros(1,size(X,2));
    cur=zeros(1,size(X,2));
    for j=1:size(X,2)
        q1(j)=normrnd(0,rmproposalsd(1));
        q2(j)=normrnd(0,rmproposalsd(1));
        q3(j)=normrnd(0,rmproposalsd(2));
        q4(j)=normrnd(0,rmproposalsd(2));
        q5(j)=normrnd(0,rmproposalsd(3));
        q6(j)=normrnd(0,rmproposalsd(3));

        newleftmu(j)=para{1,j}(1,left(j),1)+q1(j);
        newrightmu(j)=para{1,j}(1,left(j)+1,1)+q2(j);
        newleftsig(j)=para{1,j}(1,left(j),2)^2+q3(j);
        newrightsig(j)=para{1,j}(1,left(j)+1,2)^2+q4(j);
        newleftdf(j)=para{1,j}(1,left(j),3)+q5(j);
        newrightdf(j)=para{1,j}(1,left(j)+1,3)+q6(j);

        while newleftsig(j)<=0||newrightsig(j)<=0
            q3(j)=normrnd(0,rmproposalsd(2));
            q4(j)=normrnd(0,rmproposalsd(2));
            newleftsig(j)=para{1,j}(1,left(j),2)^2+q3(j);
            newrightsig(j)=para{1,j}(1,left(j)+1,2)^2+q4(j);
%             disp('Error in move sigma in move cp')
        end
        while newleftdf(j)<=2||newrightdf(j)<=2
            q5(j)=normrnd(0,rmproposalsd(3));
            q6(j)=normrnd(0,rmproposalsd(3));
            newleftdf(j)=para{1,j}(1,left(j),3)+q5(j);
            newrightdf(j)=para{1,j}(1,left(j)+1,3)+q6(j);
        end

        post(j)=loglikelihood([newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)],X(cp{1,j}(left(j))+1:procp,j),likelihood)+loglikelihood([newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)],X(procp+1:cp{1,j}(left(j)+2),j),likelihood)+logprior([newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)],hyper,likelihood)+logprior([newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)],hyper,likelihood);
        cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+loglikelihood(para{1,j}(1,left(j)+1,:),X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood)+logprior(para{1,j}(1,left(j)+1,:),hyper,likelihood);
    end
    v=rand(1);
    a=min(exp(sum(post-cur)+log(gprior(procp))-log(gprior(curcp))),1);
    %%%check if proposal in -inf
    if isinf(sum(post-cur)+log(gprior(procp))-log(gprior(curcp)))
        if isinf(sum(post))
%            disp('move')
        else
        cell(aw)
        end
    end
    if v<=a
        for j=1:size(X,2)
            newcp{1,j}(1,left(j)+1)=procp;
            newpara{1,j}(1,left(j),:)=[newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)];
            newpara{1,j}(1,left(j)+1,:)=[newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)];
        end
        newgcp(curcp)=0;
        newgcp(procp)=1;
    end
else
    disp('Need to update Likelihood')
    disp('Error move')
end

end

