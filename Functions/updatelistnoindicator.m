function [modellist,cplist,clcplist,clusterlist,icplist,paralist] = updatelistnoindicator(X,model,cp,clcp,icp,para,u,intogl,typeofmove,v,vv,numofclusters)
%Updates the list with the current values of each segment and changepoint
%locations
modellist=zeros(1,size(X,2),1); 
cplist=cell(1,size(X,2),1);
clcplist=cell(1,numofclusters);%zeros(1,size(X,1)-1,1);
clusterlist=cell(1,size(X,2));
icplist=cell(1,size(X,2),1);
paralist=cell(1,size(X,2),1);

if typeofmove==1 %i.e. we choose a independent move %v~=0 %if we have picked a series
    %independant update
    modellist(1,v)=model(1,v); %update just the selected independent series
    
    cplist(1,v)=cp(1,v);%update just the selected independent seri
    
%     gcptemp=find(gcp==1);
%     gcplist={gcptemp};

    %%%Lets not save clustercp here
%     clcptemp=find(clcp(1,:,currentcluster(v))==1);
%     clcplist{1,:,currentcluster(v)}=clcptemp;
%     clcplist{
     clcplist=cell(1,numofclusters); %i.e. no update for independent moves
    clusterlist=zeros(1,size(X,2)); %i.e. no update for independent moves
    
    icptemp=find(icp(1,:,v)==1);
    icplist(1,v)={icptemp};%update just the selected independent seri
    
    paralist(1,v)=para(1,v);%update just the selected independent seri
elseif typeofmove==2 %we have done a cluster move % Need to know which cluster moved (we have vv selected)
    %global update
    modellist=model;
    cplist=cp;
    
%     gcptemp=find(gcp==1);
%     gcplist={gcptemp};
    clcptemp=find(clcp(1,:,vv)==1);
    clcplist{1,vv}=clcptemp; %update cluster changepoint
    clusterlist=zeros(1,size(X,2)); % no change for cluster since no reassignment
    
    for j=1:size(X,2) %for each series, find the changepoint vector and update
        icptemp=find(icp(1,:,j)==1);
        icplist(1,j)={icptemp};
    end
    
    paralist=para;
else %we do cluster reassign
    
end


end

