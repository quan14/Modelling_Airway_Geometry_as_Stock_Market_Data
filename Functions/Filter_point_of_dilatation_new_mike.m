function [index_of_dilatation,transfromed_pos] = ...
    Filter_point_of_dilatation_new_mike(posterior_disturbution)
%Need to find correct change point from the distubution
%I - The posteriors disrubution - needs to be unmodifed - straght from the
%function
%O - The index of dilation and the the transfromed pos

%% Tranfromation 

transfromed_pos = posterior_disturbution;
% transfromed_pos = Recale_to_max(transfromed_pos,1);

%Need to get the max
[~,max_index] = max(transfromed_pos);
max_index = max_index + 5;

transfromed_pos(1:60) = 0;
[~,index_of_dilatation] = max(transfromed_pos);
% 
% figure
% plot(transfromed_pos)

end

