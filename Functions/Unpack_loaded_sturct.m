function [ unpacked_sturct ] = ...
    Unpack_loaded_sturct( file_name_mat )
%Returns a sturct comataining all the varaibles
% I - the name of the file
% O - the sturct conating the saved variables

%% Load and unpack

loaded_sturct = load(file_name_mat);

%Unpack by finding the field name
field_name = fieldnames(loaded_sturct);
%I'm assuming that there is obly one fiels entry
field_name = field_name{1};

%Getting the sturct array
unpacked_sturct = getfield(loaded_sturct,field_name); %#ok<GFLD>
end

