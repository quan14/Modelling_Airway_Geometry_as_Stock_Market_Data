function [posterior_pt_locations,posterior_segments] = ...
    Code_from_Micheal_RJ_MCMC_t_disturbution(Y,varargin)
%This is the code that Micheal used - the only things thats been apadetd is
%the I removed the log difference

%% KIN - the only pre - processing
Y = Y(:);

%% Pre processing of data
lrY=Y; %log rate of retrun
X=lrY; %save as X just as code is with X

%% Junk
%Not needed for you
count=1; %taking average at every 'count' i.e. count=5 means 5 day average
check=0; %check 0 means that we dont do an average portfolio step

%% Making the namiing
%Just for Naming purposes
editcompanies=cell(size(X,2),1); %create a blank cell for the titles of each series
[editcompanies{:}]=deal('Series'); %titling each series


%% Normalisation
% standardize the whole system%
X=zscore(X); %zscore standardises the matrix, by columns

%% The likelihood we use
likelihood='studentt'; %using a normal likelihood to descirbe the data
%  likelihood='studentt'; % using student t distribution

%% check the distribution and create the correct parameter types
[numberofpara,paratype]=checkdistribution(likelihood); %check the distribution and create the correct parameter types

%% No need to worry about this
%%Cluster assign %We begin with known number of clusters (With know changepoints)
%%%use gcplist method to do the clusters for now!!!!
%Known cluster, begin with 2
numofclusters=1;
currentcluster=zeros(1,size(X,2)); %allocate clusters for each series
currentcluster(1,1:1)=1;


%% set up the running parameters for MCMC
N=300000; %number of MCMC Samples
burnin=N/4;%N/10; %remove the first (N/10)/((N/10)+N)) percent
T=5; %Number of samples before we save, i.e thin from the total iterations of MCMC
apprxnum=N/T; %number of samples tosave
disp('Number of MCMC iterations')
disp(N+burnin);%samples to be saved
disp('Number of MCMC iterations used as burnin')
disp(burnin);%samples to be saved
disp('Number of Saved Samples')
disp(apprxnum);
kmax=1000; %max number of segments are set to k i.e k-1 changepoints

%% proportion of independent moves to global
%indep,global,cluster probability
unnorminglcl=[4,0,0]; %i.e 4,0,0 is only independent moves; 0,4,0 is only cluster
inglcl=unnorminglcl/sum(unnorminglcl);

%independent moves
resampleprob=3;
moveprob=3;
addprob=3;
deleteprob=3;
unnormijump=[resampleprob,moveprob,addprob,deleteprob];
ijump=unnormijump/(sum(unnormijump));

%% No need to worry as we ont have clusters
%now cluster moves
%global moves
clmoveprob=7;
claddprob=2;
cldeleteprob=2;
clustertoindep=7;
indeptocluster=7;
unnormcljump=[clmoveprob,claddprob,cldeleteprob,clustertoindep,indeptocluster];
cljump=unnormcljump/sum(unnormcljump);
%%%

%% proposal terms - split the proposals into dimension varying and within dimension (MCMC jump parameters)
%%%rm means resample and move (i.e. within mode moves)
%%%ad means add and delete (i.e dimension moving)
rmproposalsd=zeros(1,numberofpara);
rmproposaldefine=[0.1,0.2,1] %use normal distribution for degrees of freedom
adproposalsd=zeros(1,numberofpara);
adproposaldefine=[0.2,0.5,4] % use a larger variance for the adding and deleting
for i=1:numberofpara
    rmproposalsd(i)=rmproposaldefine(i); %just defining these standard deviations
    adproposalsd(i)=adproposaldefine(i); %just defining these standard deviations
end
proposalmove=3; %poisson parameter for moving a changepoint

minimumdistance=3; %mdl(minimum distance length) for choosing changepoint

cltoinbox=2; %the width of the box created around a global changepoint, pick this less than minimum distance to avoid overlapping the changepoints
proportion1or0=0.9; %proportion of global to independent moves that have changepoint

%% check the likelihood distrubtion to create the correct sized arrays
%set up empty matrix for model, cells for changepoints and parameters
[modellist,cplist,~,clcplist,clusterlist,icplist,paralist]=makelistsnoindicator(likelihood,X,kmax,(N/T),numofclusters);

%Set up the priors
rpriordist='uniform';locationr=0;%uniform prior for indepedent changepoint, no specific changepoint with weighting
% rpriordist='single';location=25; %single point prior

%Not needed for you
clpriordist='uniform';locationcl=0;% uniform prior for global changepoint, no specfic changepoint weighting
clhyper=0; %no hyper prior needed for uniform
clautoscale=1; %1 is to use global prior as rprior^(gscalefactor), !=1 is to use gscale (manual)
clscalefactor=0.2*size(X,2);%0.5 %the prior scaling on the global changepoints(normally scale with the number of series)
clscale=1; %manual scale factor for global changepoint

mpriordist='binomial'; %binomial prior on the models
mhyper=0.5; %number of expected changepoints
rhyper=0; %no hyper prior needed for uniform
rscale=0.5; %manual scale factor for independent changepoints

[rprior,clprior,mprior]=modelandcpprior(X,rpriordist,rhyper,rscale,...
    clpriordist,clhyper,clscale,clscalefactor,...
    clautoscale,locationr,locationcl,...
    mpriordist,mhyper,kmax); % define the priors for independent changepoint, cluster changepoint and model prior
hyper=[0,1,5,0.4,2,100]; %hyper(5) and hyper(6) are hyperparameters for unifrom prior on df
%%%we use the priors mu=norm(mu0,sig^2/ka0) and sig^2=scaled-inv-chi2(nu0,tau0)
%%%i.e to generate sig^2, use z=nu0 *tau0^2 / chi2rnd(nu0)
%%%i.e to generate mu, use mu=(mu0,sig^2/ka0)

%% Not needed for you
%%% Priors for Clusters
%%% for equal prob of selecting clusters
clusterprior=(1/numofclusters)*ones(1,numofclusters);

%% Just for debugging
dc=0; %check which was the last independent function called
i2clcount=0; %number of times we do indep2global moves
cl2icount=0; %number of times we do global2indep moves
%setting parameters

%% Not needed for you
jacobian=ones(1,size(X,2)); %jacobian from having a deterministic reversal jump
if strcmp(likelihood,'studentt')==1
    jacobian=2*jacobian; %using v1=v0+u, v2=v0-u
end

%% Initialise cells to hold MCMC samples
[para,cp,~,clcp,icp,model]=initialisepara(X,numberofpara,kmax,likelihood,numofclusters);
%%% need to add the starting point for the clusters
if numofclusters==2
    clcp(1,50,2)=0;
end

t=0; % counter for the number of thinning
b=0; %counter for the number of samples to burn
disp('Iteration number:')

%%%make a typeofmove list, with records of independent, cluster, and
%%%reassign chosen
[typeofmovelist,vlist,vvlist,reassignlist]=listofmoves(N/T);

%% Start of MCMC (No need to worry about elseif and else moves are they are cluster moves)
for i=1:N+burnin %from 1 to N+burnin
    if mod(i,20000)==0 %show iteration number at certain steps (in this case 20000)
        disp(i)
    end
    u=rand(1); %uniform random number between 0 and 1
    typeofmove=0;% keep track of the type of jump (i.e. indep, cluster or reassign
    v=0; %set v =0
    vv=0; %index of cluster selected
    reassign=0; %index of the series selected for cluster reassign
    if u<inglcl(1); %if less than the defined independent ratio, choose independant moves
        typeofmove=1;
        v=randi([1 size(X,2)]); %pick a series at random
        [model(1,v),cp{1,v},clcp(1,:,currentcluster(v)),icp(1,:,v),para{1,v}]=inRJMCMC(X(:,v),...
            likelihood,numberofpara,N,kmax,inglcl,ijump,cljump,adproposalsd,rmproposalsd,...
            proposalmove,minimumdistance,cltoinbox,proportion1or0,rprior,mprior,hyper,model(1,v),cp{1,v},clcp(1,:,currentcluster(v)),...
            icp(1,:,v),para{1,v},v,jacobian(v),numofclusters,currentcluster);%independant series
    elseif u<sum(inglcl(1:2))   %now cluster specific  %global series
        typeofmove=2;
        %pick a cluster
        vv=randi([1 max(currentcluster)]); % take the max as we may select a cluster that has no series in it
        %find the series in that cluster
        seriesincluster=find(currentcluster==vv);
        [model(1,seriesincluster),cp(1,seriesincluster),clcp(1,:,vv),icp(1,:,seriesincluster),para(1,seriesincluster),dc,i2clcount,cl2icount]=glRJMCMC(X(:,seriesincluster),likelihood,numberofpara,N,kmax,inglcl,ijump,cljump,adproposalsd,rmproposalsd,proposalmove,minimumdistance,cltoinbox,proportion1or0,rprior,clprior,mprior,hyper,model(1,seriesincluster),cp(1,seriesincluster),clcp(:,:,vv),icp(1,:,seriesincluster),para(1,seriesincluster),jacobian,i2clcount,cl2icount,numofclusters,currentcluster);%global series
    else %cluster reassign for a series
        typeofmove=3;
        reassign=randi([1 size(X,2)]); %pick a series at random
        %Cluster move
    end
    if b>burnin % if the counter is greater than the burnin, then we begin to save the thinned samples
        if mod(i,T)==0
            t=t+1; %add one to thinning counter
            [modellist(t,:),cplist(t,:),clcplist(t,:),clusterlist(t,:),icplist(t,:),paralist(t,:)]=updatelistnoindicator(X,model,cp,clcp,icp,para,u,inglcl,typeofmove,v,vv,numofclusters);
            [typeofmovelist(t,:),vlist(t,:),vvlist(t,:),reassignlist(t,:)]=updatemovelist(typeofmove,v,vv,reassign);
            %             clcplist(t,:)
        end
        %%
        %         %%% check minimum distance
        %         for j=1:model(1)
        %             if cp{1,1}(1,j+1)-cp{1,1}(1,j)<minimumdistance
        %                 cell(aw);
        %             end
        %         end
        
    end
    b=b+1; % add one to burnin counter
    
    if size(icp,3)==3;
        cell(aw)
    end
end

%%% can rewrite the follow scripts as we now have all the indexing


icptally % Tally of icp (independent changepoint) frequency

if isunix
    posteriorcp_cluster %plot icp posterior
else
    posteriorcp
    
end
%% KIN - Getting the output

if ~(isempty(varargin))
    close all
end

%This is the get the normalised disturbution
posterior_pt_locations = cpall;
posterior_segments = l;

end

