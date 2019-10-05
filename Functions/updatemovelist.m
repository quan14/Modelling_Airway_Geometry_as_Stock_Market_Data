function [ typeofmovelist,vlist,vvlist,reassignlist ] = updatemovelist( typeofmove,v,vv,reassign )
%Update the list with the move chosen.

%wont need this structure as only the one selected will change value. i.e.
%typeofmove will be non-zero, then either one of v, vv, reassign will be
%non zero.
% if typeofmove==1;%independent
%     
% elseif typeofmove==2; %cluster move
%     
% else %typeofmove==3; %reassign move
%     
% end

typeofmovelist=typeofmove;
vlist=v;
vvlist=vv;
reassignlist=reassign;

end

