function [newmodel,newcp,newgcp,newicp,newpara,newg2icount] = gltoindep(X,likelihood,numberofpara,N,kmax,inglcl,ijump,gjumpto,gjumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,gprior,mprior,hyper,model,cp,gcp,icp,para,jacobian,i2gcount,g2icount)
%Add changepoint to current mode for one series
newg2icount=g2icount;
newmodel=model;
newcp=cp;
newgcp=gcp;
newicp=icp;
newpara=para;
if strcmp(likelihood,'normal')==1
    
    r=randi([1 sum(gcp,2)]); %pick a index for the number of global changepoints
    gcptoselect=find(gcp==1); %find all the global changepoints
    curcp=gcptoselect(r); %pick a specific global changepoint
    
    procp=zeros(1,size(X,2));
    left=zeros(1,size(X,2));
    q1=zeros(1,size(X,2));
    q2=zeros(1,size(X,2));
    q3=zeros(1,size(X,2));
    q4=zeros(1,size(X,2));
    newleftmu=zeros(1,size(X,2));
    newrightmu=zeros(1,size(X,2));
    newleftsig=zeros(1,size(X,2));
    newrightsig=zeros(1,size(X,2));
    newmu=zeros(1,size(X,2));
    newsig=zeros(1,size(X,2));
    post=zeros(1,size(X,2));
    cur=zeros(1,size(X,2));
    newmodelprior=zeros(1,size(X,2));
    oldmodelprior=zeros(1,size(X,2));
    newtooldproposal=zeros(1,size(X,2));
    oldtonewproposal=zeros(1,size(X,2));
    % jumpnum=zeros(1,size(X,2));
    % jumpden=zeros(1,size(X,2));
    moveordel=zeros(1,size(X,2)); %move is 1, delete is 0
    possiblemoves=zeros(1,size(X,2)); %number of icp available for each series
    
        for j=1:size(X,2) %for eaach series
            left(j)=find(cp{1,j}(1,:)==curcp)-1; %find the left segment of the chosen global cp i.e cp = left+1
            u=rand(1); %pick a uniform random number
            if u<=proportion1or0 %if less than proportion 1 or 0, we keep the changepoint
                moveordel(j)=1;%keep track of moving cp
                possiblemoves(j)=sum(icp(:,:,j),2); %number of current icps but doesnt change when you move cp
                procp(j)=curcp+randi([-gltoinbox gltoinbox]); %choose a new location for the icp
                if procp(j)-cp{1,j}(left(j))<=minimumdistance||cp{1,j}(left(j)+2)-procp(j)<=minimumdistance %if inside minimum distance
%                     disp('Error in minimum distance of global to independent box')
                    return
                end
                
                %effectively move cp
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
%                     disp('Error in move sigma in move cp')
                end
                post(j)=loglikelihood([newleftmu(j),sqrt(newleftsig(j))],X(cp{1,j}(left(j))+1:procp,j),likelihood)+loglikelihood([newrightmu(j),sqrt(newrightsig(j))],X(procp+1:cp{1,j}(left(j)+2),j),likelihood)+logprior([newleftmu(j),sqrt(newleftsig(j))],hyper,likelihood)+logprior([newrightmu(j),sqrt(newrightsig(j))],hyper,likelihood)+log(rprior(j,procp(j)));
                cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+loglikelihood(para{1,j}(1,left(j)+1,:),X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood)+logprior(para{1,j}(1,left(j)+1,:),hyper,likelihood);
                %%%no need for proposal terms as we use symmetric proposals
            else
                moveordel(j)=0; %keep track of deleting a cp
                procp(j)=0; 
                %i.e delete the changepoint
                q1(j)=normrnd(0,adproposalsd(1));
                q2(j)=normrnd(0,adproposalsd(2));

                newmu(j)=mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j))+q1(j);
                newsig(j)=std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j))^2+q2(j);

                while newsig(j)<=0;
                    q2(j)=normrnd(0,adproposalsd(2));
                    newsig(j)=std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j))^2+q2(j);
%                     disp('Error in delete sig')
                end

                post(j)=loglikelihood([newmu(j),sqrt(newsig(j))],X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j),likelihood)+logprior([newmu(j),sqrt(newsig(j))],hyper,likelihood);
                cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+loglikelihood(para{1,j}(1,left(j)+1,:),X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood)+logprior(para{1,j}(1,left(j)+1,:),hyper,likelihood);

                newmodelprior(j)=log(mprior(model(j)-1));
                oldmodelprior(j)=log(mprior(model(j)));

                newtooldproposal(j)=log(fastnormalpdf(para{1,j}(1,left(j),1)-mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j)),0,adproposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j),2)^2)-std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j))^2,0,adproposalsd(2)))+log(fastnormalpdf(para{1,j}(1,left(j)+1,1)-mean(X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j)),0,adproposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j)+1,2)^2)-std(X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j))^2,0,adproposalsd(2)));
                oldtonewproposal(j)=log(fastnormalpdf(q1(j),0,adproposalsd(1)))+log(fastnormalpdf(q2(j),0,adproposalsd(2)));
            end
        end
% need to sum number of series that still have a icp for reverse move
    reverseprob=0;
    indexformove=find(moveordel==1); %find all series that we keep a changepoint
    for z=1:sum(moveordel)
        reverseprob=reverseprob+(1/size(X,2))*(1/(possiblemoves(indexformove(z))+1))*(1/(2*gltoinbox+1));
    end 

    if sum(moveordel)==0 %if we delete the global cp, reverse move is add gcp
        possiblelocations=gcp+sum(icp,3);
        locationstoaddglobal=find(possiblelocations==0);
        jumpnum=log(gjumpfrom(2)*(inglcl(2)))+log(1/(size(locationstoaddglobal,2))+1); %prob of selecting add move and to add a cp at the same location 
        jumpden=log(gjumpto(4)*(inglcl(2)))+log(1/(sum(gcp,2)))+(sum(moveordel)*log(proportion1or0))+((size(X,2)-sum(moveordel))*log(1-proportion1or0)); %prob of selecting delete move and selecting the right cp to delete
    else  
        jumpnum=log(gjumpfrom(5)*(inglcl(2)))+log(reverseprob);  %prob of selecting add move and to add a cp at the same location 
        jumpden=log(gjumpto(4)*(inglcl(2)))+log(1/(sum(gcp,2)))+sum(moveordel)*log(1/(2*gltoinbox+1))+(sum(moveordel)*log(proportion1or0))+((size(X,2)-sum(moveordel))*log(1-proportion1or0)); %prob of selecting delete move and selecting the right cp to delete
    end
    jacobian=0;
    
    a=min(exp(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal)+jumpnum-jumpden+jacobian-log(gprior(curcp))),1);
    %%%check if proposal is -inf
    if isinf(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal)+jumpnum-jumpden+jacobian-log(gprior(curcp)))
        if isinf(sum(post))
%             disp('gltoindep')
        else
        cell(aw)
        end
    end
    v=rand(1);
    if v<=a
        newg2icount=g2icount+1;
        for j=1:size(X,2)
            if moveordel(j)==1  %if we move cp
                newcp{1,j}(1,left(j)+1)=procp(j);
                newpara{1,j}(1,left(j),:)=[newleftmu(j),sqrt(newleftsig(j))];
                newpara{1,j}(1,left(j)+1,:)=[newrightmu(j),sqrt(newrightsig(j))];
                newicp(1,procp(j),j)=1;
            else %if we delete
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
        end
        newgcp(curcp)=0;
    end
    for j=1:newmodel(1)
        if newcp{1,1}(1,j+1)-newcp{1,1}(1,j)<minimumdistance
            cell(aw);
        end
    end
elseif strcmp(likelihood,'studentt')
    
    r=randi([1 sum(gcp,2)]); %pick a index for the number of global changepoints
    gcptoselect=find(gcp==1); %find all the global changepoints
    curcp=gcptoselect(r); %pick a specific global changepoint
    
    procp=zeros(1,size(X,2));
    left=zeros(1,size(X,2));
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
    newmu=zeros(1,size(X,2));
    newsig=zeros(1,size(X,2));
    newdf=zeros(1,size(X,2));
    post=zeros(1,size(X,2));
    cur=zeros(1,size(X,2));
    newmodelprior=zeros(1,size(X,2));
    oldmodelprior=zeros(1,size(X,2));
    newtooldproposal=zeros(1,size(X,2));
    oldtonewproposal=zeros(1,size(X,2));
    % jumpnum=zeros(1,size(X,2));
    % jumpden=zeros(1,size(X,2));
    moveordel=zeros(1,size(X,2)); %move is 1, delete is 0
    possiblemoves=zeros(1,size(X,2)); %number of icp available for each series
    
        for j=1:size(X,2) %for eaach series
            left(j)=find(cp{1,j}(1,:)==curcp)-1; %find the left segment of the chosen global cp i.e cp = left+1
            u=rand(1); %pick a uniform random number
            if u<=proportion1or0 %if less than proportion 1 or 0, we keep the changepoint
                moveordel(j)=1;%keep track of moving cp
                possiblemoves(j)=sum(icp(:,:,j),2); %number of current icps but doesnt change when you move cp
                procp(j)=curcp+randi([-gltoinbox gltoinbox]); %choose a new location for the icp
                if procp(j)-cp{1,j}(left(j))<=minimumdistance||cp{1,j}(left(j)+2)-procp(j)<=minimumdistance %if inside minimum distance
%                     disp('Error in minimum distance of global to independent box')
                    return
                end
                
                %effectively move cp
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
%                     disp('Error in move sigma in move cp')
                end
                
                while newleftdf(j)<=2||newrightdf(j)<=2
                    q5(j)=normrnd(0,rmproposalsd(3));
                    q6(j)=normrnd(0,rmproposalsd(3));
                    newleftdf(j)=para{1,j}(1,left(j),3)+q5(j);
                    newrightdf(j)=para{1,j}(1,left(j)+1,3)+q6(j);
                end
                
                post(j)=loglikelihood([newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)],X(cp{1,j}(left(j))+1:procp,j),likelihood)+loglikelihood([newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)],X(procp+1:cp{1,j}(left(j)+2),j),likelihood)+logprior([newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)],hyper,likelihood)+logprior([newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)],hyper,likelihood)+log(rprior(j,procp(j)));
                cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+loglikelihood(para{1,j}(1,left(j)+1,:),X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood)+logprior(para{1,j}(1,left(j)+1,:),hyper,likelihood);
                %%%no need for proposal terms as we use symmetric proposals
            else
                moveordel(j)=0; %keep track of deleting a cp
                procp(j)=0; 
                %i.e delete the changepoint
                q1(j)=normrnd(0,adproposalsd(1));
                q2(j)=normrnd(0,adproposalsd(2));

                newmu(j)=mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j))+q1(j);
                newsig(j)=std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j))^2+q2(j);
                newdf(j)=(para{1,j}(1,left(j),3)+para{1,j}(1,left(j)+1,3))/2;

                while newsig(j)<=0;
                    q2(j)=normrnd(0,adproposalsd(2));
                    newsig(j)=std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j))^2+q2(j);
%                     disp('Error in delete sig')
                end

                post(j)=loglikelihood([newmu(j),sqrt(newsig(j)),newdf(j)],X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+2),j),likelihood)+logprior([newmu(j),sqrt(newsig(j)),newdf(j)],hyper,likelihood);
                cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+loglikelihood(para{1,j}(1,left(j)+1,:),X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood)+logprior(para{1,j}(1,left(j)+1,:),hyper,likelihood);

                newmodelprior(j)=log(mprior(model(j)-1));
                oldmodelprior(j)=log(mprior(model(j)));

                newtooldproposal(j)=log(fastnormalpdf(para{1,j}(1,left(j),1)-mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j)),0,adproposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j),2)^2)-std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j))^2,0,adproposalsd(2)))+log(fastnormalpdf(para{1,j}(1,left(j)+1,1)-mean(X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j)),0,adproposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j)+1,2)^2)-std(X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j))^2,0,adproposalsd(2)))+log(fastnormalpdf(para{1,j}(1,left(j),3)-newdf(j),0,adproposalsd(3)));
                oldtonewproposal(j)=log(fastnormalpdf(q1(j),0,adproposalsd(1)))+log(fastnormalpdf(q2(j),0,adproposalsd(2)));
            end
        end
% need to sum number of series that still have a icp for reverse move
    reverseprob=0;
    indexformove=find(moveordel==1); %find all series that we keep a changepoint
    for z=1:sum(moveordel)
        reverseprob=reverseprob+(1/size(X,2))*(1/(possiblemoves(indexformove(z))+1))*(1/(2*gltoinbox+1));
    end 

    if sum(moveordel)==0 %if we delete the global cp, reverse move is add gcp
        possiblelocations=gcp+sum(icp,3);
        locationstoaddglobal=find(possiblelocations==0);
        jumpnum=log(gjumpfrom(2)*(inglcl(2)))+log(1/(size(locationstoaddglobal,2))+1); %prob of selecting add move and to add a cp at the same location 
        jumpden=log(gjumpto(4)*(inglcl(2)))+log(1/(sum(gcp,2)))+(sum(moveordel)*log(proportion1or0))+((size(X,2)-sum(moveordel))*log(1-proportion1or0)); %prob of selecting delete move and selecting the right cp to delete
    else  
        jumpnum=log(gjumpfrom(5)*(inglcl(2)))+log(reverseprob);  %prob of selecting add move and to add a cp at the same location 
        jumpden=log(gjumpto(4)*(inglcl(2)))+log(1/(sum(gcp,2)))+sum(moveordel)*log(1/(2*gltoinbox+1))+(sum(moveordel)*log(proportion1or0))+((size(X,2)-sum(moveordel))*log(1-proportion1or0)); %prob of selecting delete move and selecting the right cp to delete
    end
%     jacobian=0;
    if sum(moveordel)==size(X,2)
        logjacobian=0; %i.e dont delete cp so no need for jacobian
    else
        logjacobian=zeros(1,size(X,2));
        indexdelete=find(moveordel==0);
        for c=1:size(X,2)-sum(moveordel)
            logjacobian(indexdelete(c))=log(1./jacobian(indexdelete(c)));
        end
    end
    
    a=min(exp(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal+logjacobian)+jumpnum-jumpden-log(gprior(curcp))),1);
    %%%check if proposal is -inf
    if isinf(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal+logjacobian)+jumpnum-jumpden-log(gprior(curcp)))
        if isinf(sum(post))
%             disp('gltoindep')
        else
        cell(aw)
        end
    end
    v=rand(1);
    if v<=a
        newg2icount=g2icount+1;
        for j=1:size(X,2)
            if moveordel(j)==1  %if we move cp
                newcp{1,j}(1,left(j)+1)=procp(j);
                newpara{1,j}(1,left(j),:)=[newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)];
                newpara{1,j}(1,left(j)+1,:)=[newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)];
                newicp(1,procp(j),j)=1;
            else %if we delete
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
        end
        newgcp(curcp)=0;
    end
    for j=1:newmodel(1)
        if newcp{1,1}(1,j+1)-newcp{1,1}(1,j)<minimumdistance
            cell(aw);
        end
    end   
else
    disp('Need to update likelihood')
    disp('Error add')

end

