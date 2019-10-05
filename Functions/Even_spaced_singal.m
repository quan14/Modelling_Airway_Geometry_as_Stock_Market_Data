function [spaced_arclength,spaced_area] = ...
    Even_spaced_singal(grid_sturct)
%This get an even array so we can explore if there are any time series
%The input is:
% gird_sturct = struct;
% gird_sturct.spacing = 0.1;
% gird_sturct.int_method = 'PCHIP'; %needs to be the same text from the interpertation
% gird_sturct.arclength = arclength_array;
% gird_sturct.area = log_area_array;
%The output is the even spaced singal

%% Need to get the intorplation

%The points thhat need to be interplated
spaced_arclength = grid_sturct.arclength(1):grid_sturct.spacing:grid_sturct.arclength(end);

%placing the interpolation
spaced_area = interp1( grid_sturct.arclength , grid_sturct.area, ...
    spaced_arclength, grid_sturct.int_method );

end

