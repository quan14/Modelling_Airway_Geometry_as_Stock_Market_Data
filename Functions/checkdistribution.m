function [numberofpara,paratype] = checkdistribution(likelihood)
%Check with likelihood is used and return the number of parameters

if strcmp(likelihood,'normal')==1
    disp('Number of parameters are 2')
    numberofpara=2;
    paratype={'Mu ','Sigma '};
elseif strcmp(likelihood,'studentt')==1
    disp('Number of parameters are 3')
    numberofpara=3;
    paratype={'Mu','Sigma','Degrees of freedom'};
else
    disp('Need to update likelihood')
    numberofpara=0;
end

end

