function [new_area_1,new_area_2] = ...
    Transfrom_first_array_and_crop(area_1,area_2,pertub_step)
%Need to transfrom the first area and crop to fit
%I - the area arrays and the pertub step
%O - areas that been shifted and crop

%% Need to tranfrom

if pertub_step < 0
    
    current_area_1 = area_1((-pertub_step)+1:end);
    current_area_2 = area_2;
    
elseif pertub_step == 0
    
    current_area_1 = area_1;
    current_area_2 = area_2;
    
else
    current_area_1 = area_1;
    current_area_2 = area_2((pertub_step+1):end);
    
end

new_area_1 = current_area_1;
new_area_2 = current_area_2;

%% Need to crop

max_length = min(length(new_area_1),length(new_area_2));
new_area_1 = new_area_1(1:max_length);
new_area_2 = new_area_2(1:max_length);

end

