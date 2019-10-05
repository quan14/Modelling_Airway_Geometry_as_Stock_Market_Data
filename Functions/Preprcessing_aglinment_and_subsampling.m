function [series_differece] = ...
    Preprcessing_aglinment_and_subsampling(input_data)
%The function pre prcess the data into a 1D singal at the end

%% Doing the subsmapling

%This is the properties of the downsampling
gird_sturct = struct;
gird_sturct.spacing = 1;
gird_sturct.int_method = 'PCHIP';
gird_sturct.arclength = [];
gird_sturct.area = [];


%Need to set up the sturcts for the
baseline_airway = gird_sturct;
baseline_airway.area = input_data.baseline.area;
baseline_airway.arclength = input_data.baseline.arclength;

follow_up_airway = gird_sturct;
follow_up_airway.area = input_data.follow_up.area;
follow_up_airway.arclength = input_data.follow_up.arclength;

%perfoming alingment
[baseline_airway.sub_length,baseline_airway.sub_area] = ...
    Even_spaced_singal(baseline_airway);

[follow_up_airway.sub_length,follow_up_airway.sub_area] = ...
    Even_spaced_singal(follow_up_airway);

%% We need to crop and align the singals

%To only consider the first 50mm of the airways
croped_singal = 50;

croped_area_1 = log(baseline_airway.sub_area(1:croped_singal));
croped_area_2 = log(follow_up_airway.sub_area(1:croped_singal));

%Need to ailgn
moved_step = ...
    Alignment_of_subsample_tapering_singal(croped_area_1,...
    croped_area_2,-5:5);

%Need to shift - the shift
[baseline_airway.aligned,follow_up_airway.aligned] = ...
    Transfrom_first_array_and_crop(baseline_airway.sub_area,...
    follow_up_airway.sub_area,moved_step);

%We need to get the difference - for interpertation if a time singal a
%dilation is positive - NEED TO BE PERFROMED IN LOG SPACW
series_differece = log(follow_up_airway.aligned) - log(baseline_airway.aligned);

end

