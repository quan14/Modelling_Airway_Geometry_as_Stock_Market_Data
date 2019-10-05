function [ typeofmovelist,vlist,vvlist,reassignlist] = listofmoves(N)
%make lists to store which move was used

typeofmovelist=zeros(N,1); %zeros to store the move selected (1= independent, 2=cluster, 3=reassign)
vlist=zeros(N,1); %store which series was selected
vvlist=zeros(N,1); %store the clusters selected
reassignlist= zeros(N,1); %store which series was selected for reassign

end

