function [newmodel,newcp,newgcp,newicp,newpara] = gladd(X,likelihood,numberofpara,N,kmax,inglcl,ijump,gjumpto,gjumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,gprior,mprior,hyper,model,cp,gcp,icp,para,jacobian);
%Add changepoint to current mode for one series
newmodel=model;
newcp=cp;
newgcp=gcp;
newicp=icp;
newpara=para;

possiblelocations=gcp+sum(icp,3);
if sum(possiblelocations==0)==0
    disp('Error is picking a free place to add cp')
    return;
end

% possiblemoves=randi([1 size(X,1)-1]);
possiblemoves=find(possiblelocations==0);
procp=possiblemoves(randi([1 size(possiblemoves,2)]));
% possiblemoves=1;
% procp=9;

left=zeros(1,size(X,2));
for j=1:size(X,2)
    left(j)=find(cp{1,j}(1,:)<procp,1,'last');
    if procp-(cp{1,j}(left(j))+1)<minimumdistance||cp{1,j}(left(j)+1)-1-procp<minimumdistance
%     disp('Error in minimum distance of adding cp')
        return
    end
end

if strcmp(likelihood,'normal')==1
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
    newmodelprior=zeros(1,size(X,2));
    oldmodelprior=zeros(1,size(X,2));
    newtooldproposal=zeros(1,size(X,2));
    oldtonewproposal=zeros(1,size(X,2));
%     jumpnum=zeros(1,size(X,2));
%     jumpden=zeros(1,size(X,2));
    
    for j=1:size(X,2)
        q1(j)=normrnd(0,adproposalsd(1));
        q2(j)=normrnd(0,adproposalsd(1));
        q3(j)=normrnd(0,adproposalsd(2));
        q4(j)=normrnd(0,adproposalsd(2));

        newleftmu(j)=mean(X(cp{1,j}(left(j))+1:procp,j))+q1(j);
        newrightmu(j)=mean(X(procp+1:cp{1,j}(left(j)+1),j))+q2(j);
        newleftsig(j)=std(X(cp{1,j}(left(j))+1:procp,j))^2+q3(j);
        newrightsig(j)=std(X(procp+1:cp{1,j}(left(j)+1),j))^2+q4(j);

        while newleftsig(j)<=0||newrightsig(j)<=0
            q3(j)=normrnd(0,adproposalsd(2));
            q4(j)=normrnd(0,adproposalsd(2));
            newleftsig(j)=std(X(cp{1,j}(left(j))+1:procp,j))^2+q3(j);
            newrightsig(j)=std(X(procp+1:cp{1,j}(left(j)+1),j))^2+q4(j);
%             disp('Error in sigma in add cp')
        end

        post(j)=loglikelihood([newleftmu(j),sqrt(newleftsig(j))],X(cp{1,j}(left(j))+1:procp,j),likelihood)+loglikelihood([newrightmu(j),sqrt(newrightsig(j))],X(procp+1:cp{1,j}(left(j)+1),j),likelihood)+logprior([newleftmu(j),sqrt(newleftsig(j))],hyper,likelihood)+logprior([newrightmu(j),sqrt(newrightsig(j))],hyper,likelihood);
        cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood);

        newmodelprior(j)=log(mprior(model(j)+1));
        oldmodelprior(j)=log(mprior(model(j)));

        newtooldproposal(j)=log(fastnormalpdf(para{1,j}(1,left(j),1)-mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j)),0,adproposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j),2)^2)-std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j))^2,0,adproposalsd(2)));
        oldtonewproposal(j)=log(fastnormalpdf(q1(j),0,adproposalsd(1)))+log(fastnormalpdf(q2(j),0,adproposalsd(1)))+log(fastnormalpdf(q3(j),0,adproposalsd(2)))+log(fastnormalpdf(q4(j),0,adproposalsd(2)));

    end
%     jumpnum=log(gjumpfrom(3)*(1-intogl))+log(1/(sum(gcp,2)+1))+log(gjumpfrom(4)*(1-intogl))+size(X,2)*log(1-proportion1or0);%prob of selection delete move and selection right cp to delete and global to independent when we delete all 
     if gjump(4)==0; %if we set global to independent to zero
        jumpnum=log(gjumpfrom(3)*(inglcl(2))*(1/(sum(gcp,2)+1)));%+(gjumpfrom(4)*(1-intogl))*(1-proportion1or0)^(size(X,2)));%prob of selection delete move and selection right cp to delete and global to independent when we delete all
     else
        jumpnum=log(gjumpfrom(3)*(inglcl(2))*(1/(sum(gcp,2)+1))+(gjumpfrom(4)*(inglcl(2)))*(1-proportion1or0)^(size(X,2)));%prob of selection delete move and selection right cp to delete and global to independent when we delete all 
     end
     jumpden=log(gjumpto(2)*(inglcl(2)))+log(1/(size(possiblemoves,2)));%prob of selection add move and the prob of picking where the cp goes
    v=rand(1);
    a=min(exp(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal)+jumpnum-jumpden+log(gprior(procp))),1);

    %%%Check proposal is not -inf
    if isinf(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal)+jumpnum-jumpden+log(gprior(procp)));
        if isinf(sum(post))
%             disp('add');
        else
        cell(aw)
        end
    end
    %     cell(aw)
    %     if model==2
%         cell(aw)
%     end
    if v<=a
        for j=1:size(X,2)
            newmodel(j)=model(j)+1;
            newcp{1,j}=sort([cp{1,j},procp]);
            if left(j)==1 %if in first segment
                if left(j)+1>model %i.e if we are in model 1
                    newpara{1,j}=[];
                else %we copy the segments above the segment we cut into two
                    newpara{1,j}=newpara{1,j}(1,left(j)+1:model(j),:);
                end
                if isempty(newpara{1,j}) %i.e we move from model one to two
                    newpara{1,j}=[cat(3,newleftmu(j),sqrt(newleftsig(j))),cat(3,newrightmu(j),sqrt(newrightsig(j)))];
                else %add the two new segment para
                    newpara{1,j}=[cat(3,newleftmu(j),sqrt(newleftsig(j))),cat(3,newrightmu(j),sqrt(newrightsig(j))),newpara{1,j}];
                end
            elseif left(j)==(model(j)) %if we split the last segment of the model
                newpara{1,j}=newpara{1,j}(1,1:left(j)-1,:);
                newpara{1,j}=[newpara{1,j},cat(3,newleftmu(j),sqrt(newleftsig(j))),cat(3,newrightmu(j),sqrt(newrightsig(j)))];
            else %add the two new segments between the unaltered left and right segments
                newpara{1,j}=[newpara{1,j}(1,1:left(j)-1,:),cat(3,newleftmu(j),sqrt(newleftsig(j))),cat(3,newrightmu(j),sqrt(newrightsig(j))),newpara{1,j}(1,left(j):model(j)-1,:)];
            end
        end
        newgcp(procp)=1;
    end
elseif strcmp(likelihood,'studentt')==1
    q1=zeros(1,size(X,2));
    q2=zeros(1,size(X,2));
    q3=zeros(1,size(X,2));
    q4=zeros(1,size(X,2));
    q5=zeros(1,size(X,2)); %use a deterministic return so only one random variable needed for df
    
    newleftmu=zeros(1,size(X,2));
    newrightmu=zeros(1,size(X,2));
    newleftsig=zeros(1,size(X,2));
    newrightsig=zeros(1,size(X,2));
    newleftdf=zeros(1,size(X,2));
    newrightdf=zeros(1,size(X,2));
    
    logjacobian=zeros(1,size(X,2));
    
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
        q2(j)=normrnd(0,adproposalsd(1));
        q3(j)=normrnd(0,adproposalsd(2));
        q4(j)=normrnd(0,adproposalsd(2));
        q5(j)=normrnd(0,adproposalsd(3));

        newleftmu(j)=mean(X(cp{1,j}(left(j))+1:procp,j))+q1(j);
        newrightmu(j)=mean(X(procp+1:cp{1,j}(left(j)+1),j))+q2(j);
        newleftsig(j)=std(X(cp{1,j}(left(j))+1:procp,j))^2+q3(j);
        newrightsig(j)=std(X(procp+1:cp{1,j}(left(j)+1),j))^2+q4(j);
        newleftdf(j)=para{1,j}(1,left(j),3)+q5(j);
        newrightdf(j)=para{1,j}(1,left(j),3)-q5(j);
        
        while newleftsig(j)<=0||newrightsig(j)<=0
            q3(j)=normrnd(0,adproposalsd(2));
            q4(j)=normrnd(0,adproposalsd(2));
            newleftsig(j)=std(X(cp{1,j}(left(j))+1:procp,j))^2+q3(j);
            newrightsig(j)=std(X(procp+1:cp{1,j}(left(j)+1),j))^2+q4(j);
%             disp('Error in sigma in add cp')
        end
        while newleftdf(j)<=2||newrightdf(j)<=2
            q5(j)=normrnd(0,adproposalsd(3));
            newleftdf(j)=para{1,j}(1,left(j),3)+q5(j);
            newrightdf(j)=para{1,j}(1,left(j),3)-q5(j);
        end
        

        post(j)=loglikelihood([newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)],X(cp{1,j}(left(j))+1:procp,j),likelihood)+loglikelihood([newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)],X(procp+1:cp{1,j}(left(j)+1),j),likelihood)+logprior([newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)],hyper,likelihood)+logprior([newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)],hyper,likelihood);
        cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood);

        newmodelprior(j)=log(mprior(model(j)+1));
        oldmodelprior(j)=log(mprior(model(j)));

        newtooldproposal(j)=log(fastnormalpdf(para{1,j}(1,left(j),1)-mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j)),0,adproposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j),2)^2)-std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j))^2,0,adproposalsd(2))); %deterministic reversal move for df
        oldtonewproposal(j)=log(fastnormalpdf(q1(j),0,adproposalsd(1)))+log(fastnormalpdf(q2(j),0,adproposalsd(1)))+log(fastnormalpdf(q3(j),0,adproposalsd(2)))+log(fastnormalpdf(q4(j),0,adproposalsd(2)))+log(fastnormalpdf(q5(j),0,adproposalsd(3)));
        
        logjacobian(j)=log(jacobian(j)); %determinant of block matrix is product of diagonals
    end
%     jumpnum=log(gjumpfrom(3)*(1-intogl))+log(1/(sum(gcp,2)+1))+log(gjumpfrom(4)*(1-intogl))+size(X,2)*log(1-proportion1or0);%prob of selection delete move and selection right cp to delete and global to independent when we delete all 
     if gjump(4)==0; %if we set global to independent to zero
        jumpnum=log(gjumpfrom(3)*(inglcl(2))*(1/(sum(gcp,2)+1)));%+(gjumpfrom(4)*(1-intogl))*(1-proportion1or0)^(size(X,2)));%prob of selection delete move and selection right cp to delete and global to independent when we delete all
     else
        jumpnum=log(gjumpfrom(3)*(inglcl(2))*(1/(sum(gcp,2)+1))+(gjumpfrom(4)*(inglcl(2)))*(1-proportion1or0)^(size(X,2)));%prob of selection delete move and selection right cp to delete and global to independent when we delete all 
     end
     jumpden=log(gjumpto(2)*(inglcl(2)))+log(1/(size(possiblemoves,2)));%prob of selection add move and the prob of picking where the cp goes
    v=rand(1);
    a=min(exp(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal+logjacobian)+jumpnum-jumpden+log(gprior(procp))),1);

    %%%Check proposal is not -inf
    if isinf(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal)+jumpnum-jumpden+log(gprior(procp)));
        if isinf(sum(post))
%             disp('add');
        else
        cell(aw)
        end
    end
    %     cell(aw)
    %     if model==2
%         cell(aw)
%     end
    if v<=a
        for j=1:size(X,2)
            newmodel(j)=model(j)+1;
            newcp{1,j}=sort([cp{1,j},procp]);
            if left(j)==1 %if in first segment
                if left(j)+1>model %i.e if we are in model 1
                    newpara{1,j}=[];
                else %we copy the segments above the segment we cut into two
                    newpara{1,j}=newpara{1,j}(1,left(j)+1:model(j),:);
                end
                if isempty(newpara{1,j}) %i.e we move from model one to two
                    newpara{1,j}=[cat(3,newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)),cat(3,newrightmu(j),sqrt(newrightsig(j)),newrightdf(j))];
                else %add the two new segment para
                    newpara{1,j}=[cat(3,newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)),cat(3,newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)),newpara{1,j}];
                end
            elseif left(j)==(model(j)) %if we split the last segment of the model
                newpara{1,j}=newpara{1,j}(1,1:left(j)-1,:);
                newpara{1,j}=[newpara{1,j},cat(3,newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)),cat(3,newrightmu(j),sqrt(newrightsig(j)),newrightdf(j))];
            else %add the two new segments between the unaltered left and right segments
                newpara{1,j}=[newpara{1,j}(1,1:left(j)-1,:),cat(3,newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)),cat(3,newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)),newpara{1,j}(1,left(j):model(j)-1,:)];
            end
        end
        newgcp(procp)=1;
    end
else
    disp('Need to update likelihood')
    disp('Error add')

end

