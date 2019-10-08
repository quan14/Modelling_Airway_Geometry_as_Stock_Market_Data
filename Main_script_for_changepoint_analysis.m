%% Main script
%An implmentation of out MICCAI Paper

%% Getting the correct varaiables
clear all
addpath(genpath(pwd))
%Loading an example
input_data = Unpack_loaded_sturct('Input_aiway_data.mat');

%% Precessing data

%This alings and
[series_differece] = ...
    Preprcessing_aglinment_and_subsampling(input_data);

%% Getting the methods

[post_pt_tdis, ~] = ...
    Code_from_Micheal_RJ_MCMC_t_disturbution(series_differece,1);

posterior_dist = post_pt_tdis{1}(:,2)/sum(post_pt_tdis{1}(:,2));

%% Post prcessing 

%Finding the start of dilitation
[index_of_dilatation,~] = ...
    Filter_point_of_dilatation_new_mike(posterior_dist);

%% Ploting

%Ploting the results
figure
subplot(2,1,1)
plot(series_differece)
xlabel('Arclength mm')
ylabel('Log Area')
title('Log Cross Sectional Area Changes')

subplot(2,1,2)
plot(posterior_dist)
hold on
scatter(index_of_dilatation,posterior_dist(index_of_dilatation))
xlabel('Arclength mm')
ylabel('Probability Density')
title('Posterior p(\tau |y)')%,'Interpreter','latex')