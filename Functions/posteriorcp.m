l=cell(1,size(X,2),1); %store the posterior modelfor each series
cpfreq=cell(1,size(X,2),1); %store the frequency of changepoints in each series
modemodel=zeros(1,size(X,2)); % store the mode posterior model
seriescplistfreq=cell(1,size(X,2));

 %%%make a length(X) vector to store frequency of changepoints
cpall=cell(1,size(X,2));

modecps=cell(1,size(X,2),1); %store the most popular changepoints
for ts=1:size(X,2)
    l{1,ts}=[unique(modellist(modellist(:,ts)>0,ts)),histc(modellist(modellist(:,ts)>0,ts),unique(modellist(modellist(:,ts)>0,ts)))];

%     modemodelnumber=max(l{1,ts}(:,2));
    modemodel(1,ts)=l{1,ts}(l{1,ts}(:,2)==max(l{1,ts}(:,2)),1);
    
    seriescplist=cat(2,cplist{:,ts})'; %store all of the cps (including endpoints)
    seriescplistfreq{1,ts}=[unique(seriescplist),histc(seriescplist,unique(seriescplist))]; %get the frequency table of the changepoints (incluing endpoints)
    seriescplistfreq{1,ts}(1,:)=[];
    seriescplistfreq{1,ts}(size(seriescplistfreq{1,ts},1),:)=[];
    cpfreq{1,ts}=sortrows(seriescplistfreq{1,ts},2); %sort the rows by frequnecy
%     cpfreq{1,ts}(length(cpfreq{1,ts})-1:length(cpfreq{1,ts}),:)=[]; %delete the final two rows as they will be o and size(X,1)
    if modemodel(1,ts)==1; %if one segment model, i.e. no changepoints
        continue
    end
    tempmodecps=zeros(1,modemodel(1,ts)-1); %make a temp vector to store the mode changepoints
    tempmodecps(1)=cpfreq{1,ts}(size(cpfreq{1,ts},1),1); %save the most popular changepoint 
    if modemodel(1,ts)==2 %i.e. only one changepoint
        %move to plot graph
    else %recursively check for the next popular changepoints
        fillcp=2; %i.e we always fill from 2nd vector position
        cpcounter=modemodel(1,ts)-fillcp; %number of changepoints remaining to find
        indexcounter=length(cpfreq{1,ts})-1; %where the counter for the changepoints are (i.e. already taken the end one)
        while cpcounter>0 %if we still have a positive counter
            if sum(abs(tempmodecps-cpfreq{1,ts}(indexcounter))<=minimumdistance)>0 %if we have a minimum distance error
                indexcounter=indexcounter-1;
            else %its a cp location that is acceptable
                tempmodecps(fillcp)=cpfreq{1,ts}(indexcounter);
                fillcp=fillcp+1;
                indexcounter=indexcounter-1;
                cpcounter=cpcounter-1;
            end
        end
    end
    modecps{1,ts}=sort(tempmodecps);
    
    cpall{1,ts}(:,1)=1:size(X,1);
    for m=1:size(cpfreq{1,ts},1)
        cpall{1,ts}(cpfreq{1,ts}(m,1),2)=cpfreq{1,ts}(m,2);
    end
    
    figure
    
    subplot(3,1,1)
    plot(X(:,ts))
    if count==1
        title(strcat(editcompanies{ts,1}, ' - ',[' ','Log Returns with Mode Changepoints in Mode Model']))
    else
        title(strcat(editcompanies{ts,1},' - ',[' ',num2str(count)],' Day Average Log Returns with Mode Changepoints in Mode Model'))
    end
    xlim([0 length(X)])
    for m=1:modemodel(ts)-1
        line([modecps{1,ts}(m) modecps{1,ts}(m)],ylim,'color','r') %plot each changepoint location as a line onto of the DLR
    end
    
    subplot(3,1,2)
%     plot(seriescplistfreq{1,ts}(:,1),seriescplistfreq{1,ts}(:,2)/sum(seriescplistfreq{1,ts}(:,2)))
    plot(cpall{1,ts}(:,1),cpall{1,ts}(:,2)/sum(cpall{1,ts}(:,2)));
    xlim([0 length(X)])
    title('Posterior Density for Changepoint Locations');

    subplot(3,1,3)
    bar(l{1,ts}(:,1),l{1,ts}(:,2)/sum(l{1,ts}(:,2)))
    title('Posterior Distribution for Segment Models')
    
    
end
