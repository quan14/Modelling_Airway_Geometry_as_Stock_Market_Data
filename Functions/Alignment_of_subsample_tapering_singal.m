function [shifted_singal] = ...
    Alignment_of_subsample_tapering_singal(...
    area_1,area_2,pertubed_step)
%This foind a transformation of the singal to correct for for
%any misalignment.
%I - The area arrays and the perturbed steo
%O - The shifted singal

%% We will find the area

%Need to compare define the cost
collect_cost_vaules = [];

for i = 1:length(pertubed_step)
    
    %This is to get the paths - the perturm will only perturb area 1
    current_step = pertubed_step(i);
    
    if current_step < 0
        
        current_area_1 = area_1((-current_step)+1:end);
        current_area_2 = area_2(1:length(current_area_1));
        
    elseif current_step == 0
        
        current_area_1 = area_1;
        current_area_2 = area_2;
        
    else
        current_area_1 = area_1(1:end-current_step);
        current_area_2 = area_2((current_step+1):end);
        
    end
    
    %% The L2 norm cost
    
    differnce_vector = current_area_1 - current_area_2;
    diifernce_vector = norm(differnce_vector,2);
    
    collect_cost_vaules = cat(1,collect_cost_vaules,diifernce_vector);
    
    %Need to look at it
%     figure
%     plot(current_area_1)
%     hold on
%     plot(current_area_2)
%     xlim([-1 51])
%     title([num2str(current_step) ' ' num2str(diifernce_vector)])
%     
end

%Need to find the correct perturbation
[~,ind_arg_min] = min(collect_cost_vaules);
shifted_singal = pertubed_step(ind_arg_min);

end

