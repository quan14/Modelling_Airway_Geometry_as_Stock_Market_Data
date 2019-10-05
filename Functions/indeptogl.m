function [newmodel,newcp,newgcp,newicp,newpara,newi2gcount] = indeptogl(X,likelihood,numberofpara,N,kmax,inglcl,ijump,gjumpto,gjumpfrom,gjump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,gltoinbox,proportion1or0,rprior,gprior,mprior,hyper,model,cp,gcp,icp,para,v,jacobian,i2gcount,g2icount);
%Take a independent cp and make it into a global, v is the chosen series 
newi2gcount=i2gcount;
newmodel=model;
newcp=cp;
newgcp=gcp;
newicp=icp;
newpara=para;
if strcmp(likelihood,'normal')==1
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
    addormove=zeros(1,size(X,2));
    leftorright=zeros(1,size(X,2));
    possiblemoves=zeros(1,size(X,2));
    
    %pick a cp in the selected series
    r=randi([1 sum(icp(1,:,v))]);
    gposition=find(icp(1,:,v)==1);
    curcp=gposition(r);
    %generate a global cp
    procp=curcp+randi([-gltoinbox gltoinbox]);
%     %now for each series, check if this is within a minimum distance 
%     for j=1:size(X,2)
%         left(j)=find(cp{1,j}(1,:)<procp,1,'last');
%         if procp-(cp{1,j}(left(j))+1)<minimumdistance||cp{1,j}(left(j)+1)-1-procp<minimumdistance
%         disp('Error in minimum distance of intogl')
%             return
%         end
%     end    
    %calculate the acceptance
    
    for j=1:size(X,2)
        checkbox=find(cp{1,j}(1,:)<=procp+gltoinbox & cp{1,j}(1,:)>=procp-gltoinbox, 1);
        if checkbox==1 %start cp
            checkbox=[];
        elseif checkbox==model(j)+1 %end cp
            checkbox=[];
        end
        if isempty(checkbox)  %if there is no changepoint inside this region (find(x,1) just takes one value from the find)
            left(j)=find(cp{1,j}(1,:)<=procp+gltoinbox,1,'last'); %find the left segment below the proposed changepoint + the gltoinbox
            %we need to add a changepoint to this series
            addormove(j)=1; %1 means add
            %%%check if we are in minimum distance when we propose a new
            %%%global
            if procp-(cp{1,j}(left(j))+1)<minimumdistance||cp{1,j}(left(j)+1)-1-procp<minimumdistance
%                 disp('Error in minimum distance of intogl')
                return
            end
            %%% also check if we reach max model
            if model(j)==kmax
%                 disp('at max model in indeptogl')
                return
            end
        else %we move a cp 
            %%%check if we have two cps inside gltoinbox
            if numel(find(cp{1,j}(1,:)<=procp+gltoinbox & cp{1,j}(1,:)>=procp-gltoinbox))>1 %if we have two cps inside the gltoinbox
                %%%We reject move as we will breach minimum distance
                return
            end
            left(j)=find(cp{1,j}(1,:)<=procp+gltoinbox,1,'last')-1; %find the left segment below the proposed changepoint + the gltoinbox
            if procp==cp{1,j}(1,left(j)+1) %%No need to do anything
            else 
                %%% check if the cp to be absorb is left or right and if it is
                %%% a global changepoint
                if procp>cp{1,j}(1,left(j)+1) %procp-cp{1,j}(1,left(j)+1)<=gltoinbox %if the left is to be absorb
                    leftorright(j)=1; %keep track of if it is left or right, left is 1
                    %need to check if gcp or icp
                    if left(j)~=1 % we are not in boundary case                   
%                         if gcp(cp{1,j}(1,left(j))+1)==1 %%%%% Not sure about this %i.e the changepoint that we want to absorb is global, then we reject move
                        if gcp(cp{1,j}(1,left(j)+1))==1 %i.e the changepoint that we want to absorb is global, then we reject move
%                             disp('unable to absorb 1');
%                             cell(aw)
                            return
                        end
                        %%%also need to check if we now are in minimum distance; if
                        %%%we absorb left, then we bring it closer to right
                        if cp{1,j}(left(j)+2)-procp<minimumdistance
%                             disp('Error in minimum distance (right) of intogl')
%                             cell(aw)
                            return
                        end
                    else %if we are in boundary case, 0 is not a gcp or icp
                        %%%also need to check if we now are in minimum distance; if
                        %%%we absorb left, then we bring it closer to right
                        if cp{1,j}(left(j)+2)-procp<minimumdistance
%                             disp('Error in minimum distance (right) of intogl')
%                             cell(aw)
                            return
                        end
                    end
                else %the right is to be absorbed
                    %leftorright(j)=0 is right
                    %need to check if gcp or icp
                    if gcp(cp{1,j}(1,left(j)+1))==1 %i.e the changepoint that we want to absorb is global, then we reject move
%                         disp('unable to absorb 2');
%                         cell(aw)
                        return
                    end
                     %%%also need to check if we now are in minimum distance; if
                    %%%we absorb right, then we bring it closer to left
                    if procp-cp{1,j}(left(j))<minimumdistance
%                         disp('Error in minimum distance (left) of intogl')
                        return
                    end
                end
            end
            addormove(j)=0; %0 means move
        end
        if addormove(j)==0 %i.e we move cp
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
                post(j)=loglikelihood([newleftmu(j),sqrt(newleftsig(j))],X(cp{1,j}(left(j))+1:procp,j),likelihood)+loglikelihood([newrightmu(j),sqrt(newrightsig(j))],X(procp+1:cp{1,j}(left(j)+2),j),likelihood)+logprior([newleftmu(j),sqrt(newleftsig(j))],hyper,likelihood)+logprior([newrightmu(j),sqrt(newrightsig(j))],hyper,likelihood);
                cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+loglikelihood(para{1,j}(1,left(j)+1,:),X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood)+logprior(para{1,j}(1,left(j)+1,:),hyper,likelihood)+log(rprior(j,procp)); %%%should be existing changepoints
        else %if we add a cp
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
%                 disp('Error in sigma in add cp')
            end

            post(j)=loglikelihood([newleftmu(j),sqrt(newleftsig(j))],X(cp{1,j}(left(j))+1:procp,j),likelihood)+loglikelihood([newrightmu(j),sqrt(newrightsig(j))],X(procp+1:cp{1,j}(left(j)+1),j),likelihood)+logprior([newleftmu(j),sqrt(newleftsig(j))],hyper,likelihood)+logprior([newrightmu(j),sqrt(newrightsig(j))],hyper,likelihood);
            cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood);

            newmodelprior(j)=log(mprior(model(j)+1));
            oldmodelprior(j)=log(mprior(model(j)));

            newtooldproposal(j)=log(fastnormalpdf(para{1,j}(1,left(j),1)-mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j)),0,adproposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j),2)^2)-std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j))^2,0,adproposalsd(2)));
            oldtonewproposal(j)=log(fastnormalpdf(q1(j),0,adproposalsd(1)))+log(fastnormalpdf(q2(j),0,adproposalsd(1)))+log(fastnormalpdf(q3(j),0,adproposalsd(2)))+log(fastnormalpdf(q4(j),0,adproposalsd(2)));
        end
    end
        %%%reverse move is selecting the gltoin move, selecting the
        %%%proposed global cp, the correct proportions and the gltoinbox
        jacobian=0;
        
        jumpnum=log(gjumpfrom(4)*(inglcl(2)))+log(1/(sum(gcp,2)+1))+(size(X,2)-sum(addormove))*log(1/(2*gltoinbox+1))+sum(addormove)*log(1-proportion1or0)+(size(X,2)-sum(addormove))*log(proportion1or0);
        jumpden=log(gjumpto(5)*(inglcl(2)))+log(1/size(X,2))+log(1/sum(icp(1,:,v)))+log(1/(2*gltoinbox+1));
        a=min(exp(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal)+jumpnum-jumpden+jacobian+log(gprior(procp))),1);
        %%%chec if proposal is -inf
        if isinf(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal)+jumpnum-jumpden+jacobian+log(gprior(procp)))
            if isinf(sum(post))
%                 disp('indeptogl')
            else
            cell(aw)
            end
        end
        
        z=rand(1);
        if z<=a
            newi2gcount=i2gcount+1;
            newicp(1,curcp,v)=0;
             for j=1:size(X,2)
                 if addormove(j)==1  %if we add to get gcp
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
                else %if we move to get gcp
                    if leftorright(j)==1 %i.e we move left icp
                        newicp(1,cp{1,j}(1,left(j)+1),j)=0;
                        newcp{1,j}(1,left(j)+1)=procp;
                        newpara{1,j}(1,left(j),:)=[newleftmu(j),sqrt(newleftsig(j))];
                        newpara{1,j}(1,left(j)+1,:)=[newrightmu(j),sqrt(newrightsig(j))];
                    else %we move right icp
                        newicp(1,cp{1,j}(1,left(j)+1),j)=0;
                        newcp{1,j}(1,left(j)+1)=procp;
                        newpara{1,j}(1,left(j),:)=[newleftmu(j),sqrt(newleftsig(j))];
                        newpara{1,j}(1,left(j)+1,:)=[newrightmu(j),sqrt(newrightsig(j))];
                    end
                 end
             end
             newgcp(procp)=1;
        end
%         if addormove(j)==0; %i.e we only move cps
%             left(j)=find(cp{1,j}(1,:)<=procp+gltoinbox,1,'last')-1; %find the left segment below the proposed changepoint + the gltoinbox
%         else
%             left(j)=find(cp{1,j}(1,:)<=procp+gltoinbox,1,'last'); %find the left segment below the proposed changepoint + the gltoinbox
%         end
        %%%by construction, we can only have 
%         if procp(j)-cp{1,j}(left(j))<=gltoinbox
%             %absorb the cp inside i.e move cp
%             addormove(j)=0; %means moved
%         end
%         if cp{1,j}(left(j)+1)-procp(j)<=gltoinbox
%             %absorb the cp i.e move cp
%             addormove(j)=0; %means moved
%         end
        %%%%%need to check if the absorbed cp is independent or global
%         if procp-cp{1,j}(1,left(j))<=gltoinbox %if the proposed global cp is next to a cp on the left that is inside the global to indep box
%             if find(icp(1,cp{1,j}(left(j)),j))==1 %i.e if the cp to the left is an icp
%                 newicp(1,cp{1,j}(left(j)),j)=0; %set the left icp to zero
%             end
%         elseif procp-cp{1,j}(1,left(j)+1)<=gltoinbox %if the proposed global cp is next to a cp on the right that is inside the global to indep box
%             if find(icp(1,cp{1,j}(left(j)+1),j))==1 %i.e if the cp to the right is an icp
%                 newicp(1,cp{1,j}(left(j)+1),j)=0; %set the right icp to zero
%             end
%         end
        %ignore the absorption of cp, just add a cp
%         q1(j)=normrnd(0,proposalsd(1));
%         q2(j)=normrnd(0,proposalsd(1));
%         q3(j)=normrnd(0,proposalsd(2));
%         q4(j)=normrnd(0,proposalsd(2));
% 
%         newleftmu(j)=mean(X(cp{1,j}(left(j))+1:procp))+q1(j);
%         newrightmu(j)=mean(X(procp+1:cp{1,j}(left(j)+1)))+q2(j);
%         newleftsig(j)=std(X(cp{1,j}(left(j))+1:procp))^2+q3(j);
%         newrightsig(j)=std(X(procp+1:cp{1,j}(left(j)+1)))^2+q4(j);
%         
%         while newleftsig(j)<=0||newrightsig(j)<=0
%             q3(j)=normrnd(0,proposalsd(2));
%             q4(j)=normrnd(0,proposalsd(2));
%             newleftsig(j)=std(X(cp{1,j}(left(j))+1:procp))^2+q3(j);
%             newrightsig(j)=std(X(procp+1:cp{1,j}(left(j)+1)))^2+q4(j);
%             disp('Error in sigma in add cp for intogl')
%         end
%         
%         post(j)=loglikelihood([newleftmu(j),sqrt(newleftsig(j))],X(cp{1,j}(left(j))+1:procp),likelihood)+loglikelihood([newrightmu(j),sqrt(newrightsig(j))],X(procp+1:cp{1,j}(left(j)+1)),likelihood)+logprior([newleftmu(j),sqrt(newleftsig(j))],hyper,likelihood)+logprior([newrightmu(j),sqrt(newrightsig(j))],hyper,likelihood);
%         cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1)),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood);
% 
%         newmodelprior(j)=log(mprior(model(j)+1));
%         oldmodelprior(j)=log(mprior(model(j)));
% 
%         newtooldproposal(j)=log(fastnormalpdf(para{1,j}(1,left(j),1)-mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1))),0,proposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j),2)^2)-std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1)))^2,0,proposalsd(2)));
%         oldtonewproposal(j)=log(fastnormalpdf(q1(j),0,proposalsd(1)))+log(fastnormalpdf(q2(j),0,proposalsd(1)))+log(fastnormalpdf(q3(j),0,proposalsd(2)))+log(fastnormalpdf(q4(j),0,proposalsd(2)));
    for j=1:newmodel(1)
        if newcp{1,1}(1,j+1)-newcp{1,1}(1,j)<minimumdistance
            cell(aw);
        end
    end
    
elseif strcmp(likelihood,'studentt')
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
    addormove=zeros(1,size(X,2));
    leftorright=zeros(1,size(X,2));
    possiblemoves=zeros(1,size(X,2));
    
    %pick a cp in the selected series
    r=randi([1 sum(icp(1,:,v))]);
    gposition=find(icp(1,:,v)==1);
    curcp=gposition(r);
    %generate a global cp
    procp=curcp+randi([-gltoinbox gltoinbox]);
%     %now for each series, check if this is within a minimum distance 
%     for j=1:size(X,2)
%         left(j)=find(cp{1,j}(1,:)<procp,1,'last');
%         if procp-(cp{1,j}(left(j))+1)<minimumdistance||cp{1,j}(left(j)+1)-1-procp<minimumdistance
%         disp('Error in minimum distance of intogl')
%             return
%         end
%     end    
    %calculate the acceptance
    
    for j=1:size(X,2)
        checkbox=find(cp{1,j}(1,:)<=procp+gltoinbox & cp{1,j}(1,:)>=procp-gltoinbox, 1);
        if checkbox==1 %start cp
            checkbox=[];
        elseif checkbox==model(j)+1 %end cp
            checkbox=[];
        end
        if isempty(checkbox)  %if there is no changepoint inside this region (find(x,1) just takes one value from the find)
            left(j)=find(cp{1,j}(1,:)<=procp+gltoinbox,1,'last'); %find the left segment below the proposed changepoint + the gltoinbox
            %we need to add a changepoint to this series
            addormove(j)=1; %1 means add
            %%%check if we are in minimum distance when we propose a new
            %%%global
            if procp-(cp{1,j}(left(j))+1)<minimumdistance||cp{1,j}(left(j)+1)-1-procp<minimumdistance
%                 disp('Error in minimum distance of intogl')
                return
            end
            %%% also check if we reach max model
            if model(j)==kmax
%                 disp('at max model in indeptogl')
                return
            end
        else %we move a cp 
            %%%check if we have two cps inside gltoinbox
            if numel(find(cp{1,j}(1,:)<=procp+gltoinbox & cp{1,j}(1,:)>=procp-gltoinbox))>1 %if we have two cps inside the gltoinbox
                %%%We reject move as we will breach minimum distance
                return
            end
            left(j)=find(cp{1,j}(1,:)<=procp+gltoinbox,1,'last')-1; %find the left segment below the proposed changepoint + the gltoinbox
            if procp==cp{1,j}(1,left(j)+1) %%No need to do anything
            else 
                %%% check if the cp to be absorb is left or right and if it is
                %%% a global changepoint
                if procp>cp{1,j}(1,left(j)+1) %procp-cp{1,j}(1,left(j)+1)<=gltoinbox %if the left is to be absorb
                    leftorright(j)=1; %keep track of if it is left or right, left is 1
                    %need to check if gcp or icp
                    if left(j)~=1 % we are not in boundary case                   
%                         if gcp(cp{1,j}(1,left(j))+1)==1 %%%%% Not sure about this %i.e the changepoint that we want to absorb is global, then we reject move
                        if gcp(cp{1,j}(1,left(j)+1))==1 %i.e the changepoint that we want to absorb is global, then we reject move
%                             disp('unable to absorb 1');
%                             cell(aw)
                            return
                        end
                        %%%also need to check if we now are in minimum distance; if
                        %%%we absorb left, then we bring it closer to right
                        if cp{1,j}(left(j)+2)-procp<minimumdistance
%                             disp('Error in minimum distance (right) of intogl')
%                             cell(aw)
                            return
                        end
                    else %if we are in boundary case, 0 is not a gcp or icp
                        %%%also need to check if we now are in minimum distance; if
                        %%%we absorb left, then we bring it closer to right
                        if cp{1,j}(left(j)+2)-procp<minimumdistance
%                             disp('Error in minimum distance (right) of intogl')
%                             cell(aw)
                            return
                        end
                    end
                else %the right is to be absorbed
                    %leftorright(j)=0 is right
                    %need to check if gcp or icp
                    if gcp(cp{1,j}(1,left(j)+1))==1 %i.e the changepoint that we want to absorb is global, then we reject move
%                         disp('unable to absorb 2');
%                         cell(aw)
                        return
                    end
                     %%%also need to check if we now are in minimum distance; if
                    %%%we absorb right, then we bring it closer to left
                    if procp-cp{1,j}(left(j))<minimumdistance
%                         disp('Error in minimum distance (left) of intogl')
                        return
                    end
                end
            end
            addormove(j)=0; %0 means move
        end
        if addormove(j)==0 %i.e we move cp
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
                
                post(j)=loglikelihood([newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)],X(cp{1,j}(left(j))+1:procp,j),likelihood)+loglikelihood([newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)],X(procp+1:cp{1,j}(left(j)+2),j),likelihood)+logprior([newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)],hyper,likelihood)+logprior([newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)],hyper,likelihood);
                cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j),likelihood)+loglikelihood(para{1,j}(1,left(j)+1,:),X(cp{1,j}(left(j)+1)+1:cp{1,j}(left(j)+2),j),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood)+logprior(para{1,j}(1,left(j)+1,:),hyper,likelihood)+log(rprior(j,procp)); %%%should be existing changepoints
        else %if we add a cp
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
%                 disp('Error in sigma in add cp')
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

            newtooldproposal(j)=log(fastnormalpdf(para{1,j}(1,left(j),1)-mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j)),0,adproposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j),2)^2)-std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1),j))^2,0,adproposalsd(2)));
            oldtonewproposal(j)=log(fastnormalpdf(q1(j),0,adproposalsd(1)))+log(fastnormalpdf(q2(j),0,adproposalsd(1)))+log(fastnormalpdf(q3(j),0,adproposalsd(2)))+log(fastnormalpdf(q4(j),0,adproposalsd(2)))+log(fastnormalpdf(q5(j),0,adproposalsd(3)));
        end
    end
        %%%reverse move is selecting the gltoin move, selecting the
        %%%proposed global cp, the correct proportions and the gltoinbox
        if sum(addormove)==0 %i.e all are move
            logjacobian=0;
        else
            logjacobian=zeros(1,size(X,2));
            indexadd=find(addormove==1);
            for c=1:sum(addormove)
                logjacobian(indexadd(c))=log(jacobian(indexadd(c))); %i.e add a jacobian ofr each add
            end
        end
        
        jumpnum=log(gjumpfrom(4)*(inglcl(2)))+log(1/(sum(gcp,2)+1))+(size(X,2)-sum(addormove))*log(1/(2*gltoinbox+1))+sum(addormove)*log(1-proportion1or0)+(size(X,2)-sum(addormove))*log(proportion1or0);
        jumpden=log(gjumpto(5)*(inglcl(2)))+log(1/size(X,2))+log(1/sum(icp(1,:,v)))+log(1/(2*gltoinbox+1));
        a=min(exp(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal+logjacobian)+jumpnum-jumpden+log(gprior(procp))),1);
        %%%chec if proposal is -inf
        if isinf(sum(post-cur+newmodelprior-oldmodelprior+newtooldproposal-oldtonewproposal+logjacobian)+jumpnum-jumpden+log(gprior(procp)))
            if isinf(sum(post))
%                 disp('indeptogl')
            else
            cell(aw)
            end
        end
        
        z=rand(1);
        if z<=a
            newi2gcount=i2gcount+1;
            newicp(1,curcp,v)=0;
             for j=1:size(X,2)
                 if addormove(j)==1  %if we add to get gcp
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
                else %if we move to get gcp
                    if leftorright(j)==1 %i.e we move left icp
                        newicp(1,cp{1,j}(1,left(j)+1),j)=0;
                        newcp{1,j}(1,left(j)+1)=procp;
                        newpara{1,j}(1,left(j),:)=[newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)];
                        newpara{1,j}(1,left(j)+1,:)=[newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)];
                    else %we move right icp
                        newicp(1,cp{1,j}(1,left(j)+1),j)=0;
                        newcp{1,j}(1,left(j)+1)=procp;
                        newpara{1,j}(1,left(j),:)=[newleftmu(j),sqrt(newleftsig(j)),newleftdf(j)];
                        newpara{1,j}(1,left(j)+1,:)=[newrightmu(j),sqrt(newrightsig(j)),newrightdf(j)];
                    end
                 end
             end
             newgcp(procp)=1;
        end
%         if addormove(j)==0; %i.e we only move cps
%             left(j)=find(cp{1,j}(1,:)<=procp+gltoinbox,1,'last')-1; %find the left segment below the proposed changepoint + the gltoinbox
%         else
%             left(j)=find(cp{1,j}(1,:)<=procp+gltoinbox,1,'last'); %find the left segment below the proposed changepoint + the gltoinbox
%         end
        %%%by construction, we can only have 
%         if procp(j)-cp{1,j}(left(j))<=gltoinbox
%             %absorb the cp inside i.e move cp
%             addormove(j)=0; %means moved
%         end
%         if cp{1,j}(left(j)+1)-procp(j)<=gltoinbox
%             %absorb the cp i.e move cp
%             addormove(j)=0; %means moved
%         end
        %%%%%need to check if the absorbed cp is independent or global
%         if procp-cp{1,j}(1,left(j))<=gltoinbox %if the proposed global cp is next to a cp on the left that is inside the global to indep box
%             if find(icp(1,cp{1,j}(left(j)),j))==1 %i.e if the cp to the left is an icp
%                 newicp(1,cp{1,j}(left(j)),j)=0; %set the left icp to zero
%             end
%         elseif procp-cp{1,j}(1,left(j)+1)<=gltoinbox %if the proposed global cp is next to a cp on the right that is inside the global to indep box
%             if find(icp(1,cp{1,j}(left(j)+1),j))==1 %i.e if the cp to the right is an icp
%                 newicp(1,cp{1,j}(left(j)+1),j)=0; %set the right icp to zero
%             end
%         end
        %ignore the absorption of cp, just add a cp
%         q1(j)=normrnd(0,proposalsd(1));
%         q2(j)=normrnd(0,proposalsd(1));
%         q3(j)=normrnd(0,proposalsd(2));
%         q4(j)=normrnd(0,proposalsd(2));
% 
%         newleftmu(j)=mean(X(cp{1,j}(left(j))+1:procp))+q1(j);
%         newrightmu(j)=mean(X(procp+1:cp{1,j}(left(j)+1)))+q2(j);
%         newleftsig(j)=std(X(cp{1,j}(left(j))+1:procp))^2+q3(j);
%         newrightsig(j)=std(X(procp+1:cp{1,j}(left(j)+1)))^2+q4(j);
%         
%         while newleftsig(j)<=0||newrightsig(j)<=0
%             q3(j)=normrnd(0,proposalsd(2));
%             q4(j)=normrnd(0,proposalsd(2));
%             newleftsig(j)=std(X(cp{1,j}(left(j))+1:procp))^2+q3(j);
%             newrightsig(j)=std(X(procp+1:cp{1,j}(left(j)+1)))^2+q4(j);
%             disp('Error in sigma in add cp for intogl')
%         end
%         
%         post(j)=loglikelihood([newleftmu(j),sqrt(newleftsig(j))],X(cp{1,j}(left(j))+1:procp),likelihood)+loglikelihood([newrightmu(j),sqrt(newrightsig(j))],X(procp+1:cp{1,j}(left(j)+1)),likelihood)+logprior([newleftmu(j),sqrt(newleftsig(j))],hyper,likelihood)+logprior([newrightmu(j),sqrt(newrightsig(j))],hyper,likelihood);
%         cur(j)=loglikelihood(para{1,j}(1,left(j),:),X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1)),likelihood)+logprior(para{1,j}(1,left(j),:),hyper,likelihood);
% 
%         newmodelprior(j)=log(mprior(model(j)+1));
%         oldmodelprior(j)=log(mprior(model(j)));
% 
%         newtooldproposal(j)=log(fastnormalpdf(para{1,j}(1,left(j),1)-mean(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1))),0,proposalsd(1)))+log(fastnormalpdf((para{1,j}(1,left(j),2)^2)-std(X(cp{1,j}(left(j))+1:cp{1,j}(left(j)+1)))^2,0,proposalsd(2)));
%         oldtonewproposal(j)=log(fastnormalpdf(q1(j),0,proposalsd(1)))+log(fastnormalpdf(q2(j),0,proposalsd(1)))+log(fastnormalpdf(q3(j),0,proposalsd(2)))+log(fastnormalpdf(q4(j),0,proposalsd(2)));
    for j=1:newmodel(1)
        if newcp{1,1}(1,j+1)-newcp{1,1}(1,j)<minimumdistance
            cell(aw);
        end
    end
else 
    disp('Need to change Likelihood')
    disp('Error in indeptogl')
end
end

