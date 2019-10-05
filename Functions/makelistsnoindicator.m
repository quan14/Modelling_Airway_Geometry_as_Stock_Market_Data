function [model,cpoint,gcplist,clcplist,clusterlist,icplist,para] = makelistsnoindicator(nameofdist,X,kmax,N,numofclusters)
%Checks the distrubtion used for likelihood and sets up the correct sized
%matrices and cells. model is the model matrix; cpoint is the changepoint
%cells and para is the parameters cells. nameofdistof a string of only two options: normal or
%locationscalet; X is data, kmax is maximum number of segments; N is the
%number of MCMC iterations

if strcmp(nameofdist,'normal')==1
    model=zeros(N,size(X,2));
    cpoint=cell(N,size(X,2),1);
    para=cell(N,size(X,2),1);
    gcplist=cell(N,1);%zeros(N,size(X,1)-1); %globalchangepoint positions
    clcplist=cell(N,numofclusters);%zeros(N,size(X,1)-1); %globalchangepoint positions
    clusterlist=zeros(N,size(X,2));%cell(N,size(X,2)); %the posterior cluster frequency of each series
    icplist=cell(N,size(X,2));
elseif strcmp(nameofdist,'studentt')==1
    model=zeros(N,size(X,2));
    cpoint=cell(N,size(X,2),1);
    para=cell(N,size(X,2),1);
    gcplist=cell(N,1);%zeros(N,size(X,1)-1); %globalchangepoint positions
    clcplist=cell(N,numofclusters);%zeros(N,size(X,1)-1); %globalchangepoint positions
    clusterlist=zeros(N,size(X,2));%cell(N,size(X,2)); %the posterior cluster frequency of each series
    icplist=cell(N,size(X,2));
else
    disp('Need to update likelihood')
    disp('makelist')
end

