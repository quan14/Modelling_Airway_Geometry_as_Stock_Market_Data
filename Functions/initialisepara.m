function [para,cp,gcp,clcp,icp,model] = initialisepara(X,numberofpara,kmax,likelihood,numofclusters)
%Initial a cell to hold the current values of the parameters, changepoints
%and model
model=ones(1,size(X,2),1); %i.e in the one segment model
gcp=zeros(1,size(X,1)-1); %i.e we have no global changepoints
clcp=zeros(1,size(X,1)-1,numofclusters);
icp=zeros(1,size(X,1)-1,size(X,2)); %i.e we have no independent changepoint
para=cell(1,size(X,2),1); %set up size(X,2) by numberofpara cell
[para{:}]=deal(zeros(1,1,numberofpara)); %set up empty matrices in each cell
cp=cell(1,size(X,2),1); 
[cp{:}]=deal([0,size(X,1)]);
if strcmp(likelihood,'normal')==1
    for j=1:size(X,2)
        %mu and sigma
        para{1,j}(1,1,1)=mean(X(:,j));
        para{1,j}(1,1,2)=std(X(:,j));
    end
elseif strcmp(likelihood,'studentt')==1
    for j=1:size(X,2)
        %mu and sigma
        para{1,j}(1,1,1)=mean(X(:,j));
        para{1,j}(1,1,2)=std(X(:,j));
        para{1,j}(1,1,3)=randi([59 60]); %arbitrary set the degrees of freedom
    end
else
    disp('Need to update likelihood');
    disp('ERROR initialisepara')
end
end

