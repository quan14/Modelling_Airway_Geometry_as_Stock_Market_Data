icpfreq=cell(1,size(X,2)); %store sorted frequency of icp for each series
seriesicplistfreq=cell(1,size(X,2)); %store tally of each location as a cp
l1=cell(1,size(X,2)); %store frequency of cp model for each series
icpall=cell(1,size(X,2)); %%store cp frequency across all time points for each series

for ts=1:size(X,2)
    l1{1,ts}=[unique(modellist(modellist(:,ts)>0,ts)),histc(modellist(modellist(:,ts)>0,ts),unique(modellist(modellist(:,ts)>0,ts)))];
    
    seriesicplist=cat(2,icplist{:,ts})'; %store all of the cps (including endpoints)
    seriesicplistfreq{1,ts}=[unique(seriesicplist),histc(seriesicplist,unique(seriesicplist))]; %get the frequency table of the changepoints (incluing endpoints)
    if ~isempty(seriesicplistfreq{1,ts})
        %dont need to delete first and last as we dont keep 0 and length(X)
%         seriesicplistfreq{1,ts}(1,:)=[];
%         seriesicplistfreq{1,ts}(size(seriesicplistfreq{1,ts},1),:)=[];
        icpfreq{1,ts}=sortrows(seriesicplistfreq{1,ts},2); %sort the rows by frequnecy
    else
    end
    icpall{1,ts}(:,1)=1:size(X,1);
    for m=1:size(icpfreq{1,ts},1)
        icpall{1,ts}(icpfreq{1,ts}(m,1),2)=icpfreq{1,ts}(m,2);
    end
end