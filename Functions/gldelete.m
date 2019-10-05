function [newmodel,newcp,newgcp,newicp,newpara] = gldelete(X,likelihood,numberofpara,N,kmax,inglcl,ijump,gjumpto,gjumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,gprior,mprior,hyper,model,cp,gcp,icp,para,jacobian);
%delete an independent changepoint from a series
newmodel=model;
newcp=cp;
newgcp=gcp;
newicp=icp;
newpara=para;

%need possiblemoves for the reverse jump of adding in all possible places
%+1 for the deleted cp
possiblelocations=gcp+sum(icp,3);
possiblemoves=find(possiblelocations==0);
%%%

% if sum(gcp,2)==(model-1) %only needed for indelete
% %     disp('Error in picking model to delete')
%     return %leave function
% end
left=zeros(1,size(X,2));
r=randi([1 sum(gcp,2)]);
gcptodelete=find(gcp==1); %find possible gcp that can be moved
curcp=gcptodelete(r); %pick one 
for j=1:size(X,2)
    left(j)=(find(cp{1,j}(1,:)==curcp))-1; %find the left segment
end
% curcp=9;
% possiblemoves=[];
% left=1;

if strcmp(likelihood,'normal')
    q1=zeros(1,size(X,2));
    q2=zeros(1,size(X,2));
    newmu=zeros(1,size(X,2));
    newsig=zeros(1,size(X,2));
    post=zeros(1,size(X,2));
    cur=zeros(1,size(X,2));
    newmodelprior=zeros(1,size(X,2));
    oldmodelprior=zeros(1,size(X,2));
    newtooldproposal=zeros(1,size(X,2));
    oldtonewproposal=zeros(1,size(X,2));
%     jumpnum=zeros(1,size(X,2));
%     jumpden=zeros(1,size(X,2));
    for j=1:size(X,2)
        q1(j)=normrnd(0,adproposalsd(1));
        q2(j)=normrnd(0,adproposalsd(2));

        newmu(j)=mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j))+q1(j);
        newsig(j)=std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j))^2+q2(j);

        while newsig(j)<=0;
            q2(j)=normrnd(0,adproposalsd(2));
            newsig(j)=std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j))^2+q2(j);
%             disp('Error in delete sig')
        end

        post(j)=loglikelihood([newmu(j),sqrt(newsig(j))],X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j),likelihood)+logprior([newmu(j),sqrt(newsig(j))],hyper,likelihood);
        cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+loglikelihood(para{1,j}(1,left(j)+1,:),X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood)+logprior(para{1,j}(1,left(j)+1,:),hyper,likelihood);

        newmodelprior(j)=log(mprior(model(j)-1));
        oldmodelprior(j)=log(mprior(model(j)));

        newtooldproposal(j)=log(fastnormalpdf(para{1,j}(1,left(j),1)-mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j)),0,adproposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j),2)^2)-std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j))^2,0,adproposalsd(2)))+log(fastnormalpdf(para{1,j}(1,left(j)+1,1)-mean(X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j)),0,adproposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j)+1,2)^2)-std(X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j))^2,0,adproposalsd(2)));
        oldtonewproposal(j)=log(fastnormalpdf(q1(j),0,adproposalsd(1)))+log(fastnormalpdf(q2(j),0,adproposalsd(2)));

    end
    jumpnum=log(gjumpfrom(2)*(inglcl(2)))+log(1/(size(possiblemoves,2)+1)); %prob of selecting add move and to add a cp at the same location 
    jumpden=log(gjumpto(3)*(inglcl(2)))+log(1/(sum(gcp,2))); %prob of selecting delete move and selecting the right cp to delete
    v=rand(1);
    a=min(exp(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal)+jumpnum-jumpden-log(gprior(curcp))),1);

    %%%check if proposal is -inf
    if isinf(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal)+jumpnum-jumpden-log(gprior(curcp)))
        if isinf(sum(post))
%             disp('delete')
        else
        cell(aw)
        end
    end
    %     if model==4;
%         cell(Aw)
%     end
    if v<=a
        for j=1:size(X,2)
            newmodel(j)=model(j)-1;
            newcp{1,j}(left(j)+1)=[];
            if left(j)==1
                if left(j)+2<=model(j)
                    newpara{1,j}=newpara{1,j}(1,left(j)+2:model(j),:);
                else
                    newpara{1,j}=[];
                end
                if isempty(newpara{1,j})
                    newpara{1,j}=cat(3,newmu(j),sqrt(newsig(j)));
                else
                    newpara{1,j}=[cat(3,newmu(j),sqrt(newsig(j))),newpara{1,j}];
                end
            elseif left(j)==(model(j)-1)
                newpara{1,j}=newpara{1,j}(1,1:left(j)-1,:);
                newpara{1,j}=[newpara{1,j},cat(3,newmu(j),sqrt(newsig(j)))];
            else
                newpara{1,j}=[newpara{1,j}(1,1:left(j)-1,:),cat(3,newmu(j),sqrt(newsig(j))),newpara{1,j}(1,left(j)+2:model(j),:)];
            end
        end
        newgcp(curcp)=0;
    end
elseif strcmp(likelihood,'studentt')==1
    q1=zeros(1,size(X,2));
    q2=zeros(1,size(X,2));
    newmu=zeros(1,size(X,2));
    newsig=zeros(1,size(X,2));
    newdf=zeros(1,size(X,2));
    post=zeros(1,size(X,2));
    cur=zeros(1,size(X,2));
    logjacobian=zeros(1,size(X,2));
    newmodelprior=zeros(1,size(X,2));
    oldmodelprior=zeros(1,size(X,2));
    newtooldproposal=zeros(1,size(X,2));
    oldtonewproposal=zeros(1,size(X,2));
%     jumpnum=zeros(1,size(X,2));
%     jumpden=zeros(1,size(X,2));
    for j=1:size(X,2)
        q1(j)=normrnd(0,adproposalsd(1));
        q2(j)=normrnd(0,adproposalsd(2));

        newmu(j)=mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j))+q1(j);
        newsig(j)=std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j))^2+q2(j);
        newdf(j)=(para{1,j}(1,left(j),3)+para{1,j}(1,left(j)+1,3))/2;        

        while newsig(j)<=0;
            q2(j)=normrnd(0,adproposalsd(2));
            newsig(j)=std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j))^2+q2(j);
%             disp('Error in delete sig')
        end

        post(j)=loglikelihood([newmu(j),sqrt(newsig(j)),newdf(j)],X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j),likelihood)+logprior([newmu(j),sqrt(newsig(j)),newdf],hyper,likelihood);
        cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+loglikelihood(para{1,j}(1,left(j)+1,:),X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood)+logprior(para{1,j}(1,left(j)+1,:),hyper,likelihood);

        newmodelprior(j)=log(mprior(model(j)-1));
        oldmodelprior(j)=log(mprior(model(j)));

        newtooldproposal(j)=log(fastnormalpdf(para{1,j}(1,left(j),1)-mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j)),0,adproposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j),2)^2)-std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j))^2,0,adproposalsd(2)))+log(fastnormalpdf(para{1,j}(1,left(j)+1,1)-mean(X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j)),0,adproposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j)+1,2)^2)-std(X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j))^2,0,adproposalsd(2)))+log(fastnormalpdf(para{1,j}(1,left(j),3)-newdf(j),0,adproposalsd(3)));
        oldtonewproposal(j)=log(fastnormalpdf(q1(j),0,adproposalsd(1)))+log(fastnormalpdf(q2(j),0,adproposalsd(2)));

    end
    jumpnum=log(gjumpfrom(2)*(inglcl(2)))+log(1/(size(possiblemoves,2)+1)); %prob of selecting add move and to add a cp at the same location 
    jumpden=log(gjumpto(3)*(inglcl(2)))+log(1/(sum(gcp,2))); %prob of selecting delete move and selecting the right cp to delete
    
    logjacobian=log(1./jacobian);
    
    v=rand(1);
    a=min(exp(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal+logjacobian)+jumpnum-jumpden-log(gprior(curcp))),1);

    %%%check if proposal is -inf
    if isinf(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal)+jumpnum-jumpden-log(gprior(curcp)))
        if isinf(sum(post))
%             disp('delete')
        else
        cell(aw)
        end
    end
    %     if model==4;
%         cell(Aw)
%     end
    if v<=a
        for j=1:size(X,2)
            newmodel(j)=model(j)-1;
            newcp{1,j}(left(j)+1)=[];
            if left(j)==1
                if left(j)+2<=model(j)
                    newpara{1,j}=newpara{1,j}(1,left(j)+2:model(j),:);
                else
                    newpara{1,j}=[];
                end
                if isempty(newpara{1,j})
                    newpara{1,j}=cat(3,newmu(j),sqrt(newsig(j)),newdf(j));
                else
                    newpara{1,j}=[cat(3,newmu(j),sqrt(newsig(j)),newdf(j)),newpara{1,j}];
                end
            elseif left(j)==(model(j)-1)
                newpara{1,j}=newpara{1,j}(1,1:left(j)-1,:);
                newpara{1,j}=[newpara{1,j},cat(3,newmu(j),sqrt(newsig(j)),newdf(j))];
            else
                newpara{1,j}=[newpara{1,j}(1,1:left(j)-1,:),cat(3,newmu(j),sqrt(newsig(j)),newdf(j)),newpara{1,j}(1,left(j)+2:model(j),:)];
            end
        end
        newgcp(curcp)=0;
    end    
else
    disp('Need to update likelihood')
    disp('Error in delete')
end

end

