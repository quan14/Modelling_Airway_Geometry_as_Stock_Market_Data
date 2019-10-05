function [newmodel,newcp,newclcp,newicp,newpara] = inresample(X,likelihood,numberofpara,N,kmax,intogl,ijump,ijumpto,ijumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,mprior,hyper,model,cp,clcp,icp,para,series,jacobian,numofclusters,currentcluster)
%resample new parameters for the existing changepoint model, series is the
%series selected
newmodel=model;
newcp=cp;
newclcp=clcp;
newicp=icp;
newpara=para;

if strcmp(likelihood,'normal')==1
    q=zeros(1,model,numberofpara);
    for j=1:model
        q(1,j,1)=normrnd(para(1,j,1),rmproposalsd(1));
        q(1,j,2)=normrnd(para(1,j,2)^2,rmproposalsd(2));
        
        while q(1,j,2)<=0;
            q(1,j,2)=normrnd(para(1,j,2)^2,rmproposalsd(2));
%             disp('Error in resample sig')
        end
        
        post=loglikelihood([q(1,j,1),sqrt(q(1,j,2))],X(cp(1,j)+1:cp(j+1)),likelihood)+logprior([q(1,j,1),sqrt(q(1,j,2))],hyper,likelihood);
        cur=loglikelihood(para(1,j,:),X(cp(1,j)+1:cp(j+1)),likelihood)+logprior(para(1,j,:),hyper,likelihood);
        
        a=min(exp(post-cur),1); %no need for the reverseble jump as we use a symmetric proposal
        u=rand(1);
%         if model==8            
%             if para(1,8,2)^2>358
%                 cell(aw);
%             end
%         end
        if u<=a;
            newpara(1,j,1)=q(1,j,1);
            newpara(1,j,2)=sqrt(q(1,j,2));
        end
    end
elseif strcmp(likelihood,'studentt')==1
    q=zeros(1,model,numberofpara);
    for j=1:model
        q(1,j,1)=normrnd(para(1,j,1),rmproposalsd(1));
        q(1,j,2)=normrnd(para(1,j,2)^2,rmproposalsd(2));
        q(1,j,3)=normrnd(para(1,j,3),rmproposalsd(3));
        
        %might not need these checks, as the prior should take care of this
        while q(1,j,2)<=0;
            q(1,j,2)=normrnd(para(1,j,2)^2,rmproposalsd(2));
%             disp('Error in resample sig')
        end
        while q(1,j,3)<=2 %if the degrees of freedom is less than 2
            q(1,j,3)=normrnd(para(1,j,3),rmproposalsd(3));
        end
        
        post=loglikelihood([q(1,j,1),sqrt(q(1,j,2)),q(1,j,3)],X(cp(1,j)+1:cp(j+1)),likelihood)+logprior([q(1,j,1),sqrt(q(1,j,2)),q(1,j,3)],hyper,likelihood);
        cur=loglikelihood(para(1,j,:),X(cp(1,j)+1:cp(j+1)),likelihood)+logprior(para(1,j,:),hyper,likelihood);
        
        a=min(exp(post-cur),1); %no need for the reverseble jump as we use a symmetric proposal
        u=rand(1);
%         if model==8            
%             if para(1,8,2)^2>358
%                 cell(aw);
%             end
%         end
        if u<=a;
            newpara(1,j,1)=q(1,j,1);
            newpara(1,j,2)=sqrt(q(1,j,2));
            newpara(1,j,3)=q(1,j,3);
        end
    end   
else
    disp('Need to update Likelihood')
    disp('Error Resample')
end
end

