function output_vector = Return_normalise_dis( scale_vector )
%This function returns the unit normal vector

%The input is the scale vector of the form [x y z]

%The output is the unit vector of the form [x y z]

%% Normalising the vector

output_vector = scale_vector/sum(scale_vector);

end

