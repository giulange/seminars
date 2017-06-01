function scores = ga_pcr_fitness(x, perf, O, T)
% function scores = ga_pcr_fitness(x, perf, O, T)
% 
% |-------------------------|
% | CREATED on 30-May-2017  |
% |-------------------------|
%
% DESCRIPTION:
%  This function is a fitness function to be used with the genetic algorithm toolbox.
%  The 'x' input is used to select how many and which NN components should be used
%  to conform ensemble output aggregated by PCR with target data.
%  It is assumed that training was performed using Tr to calibrate weights and biases,
%  Va as a method of early stopping criterion and Te as the real world signals
%  evaluating the performance.
%
% INPUTS:
%	x      : population [PopulationSize x GenomeLength] if 'UseVectorized', OR
%            nvars		[GenomeLength x 1] otherwise
%            The more general case (population) was here implemented.
%	perf   : selects the performance indicator used in the GA/PCR method.
%                1: 'mse'
%                2: 'rmse'
%                3: 'mae'
%                4: 'mbe2'
%                5: 'r'
%                6: 'smape'
%                7: 'eff'
%                8: 'D'
%	O      : input  data, structure array containing {otr,ova,ote}.
%	T      : output data, structure array containing {ttr,tva,tte}.
%
% STEPS:
%  1.\ select which metrics to use;
%  2.\ use 'x' to select NN components (how many and which);
%  3.\ compute ensemble response aggregating results by PCR;
%       a. cross-validation weights on Tr       [PCRcv.m]
%       b. BAGAP response           on Te       [bagnet.m]
%       c. select N eigenvectors    on Te       [MSEg.m]
%       d. BAGAP response           on Va       [bagnet.m]
%       e. performance              on Va       [perfind.m]
%  4.\ extract selected metric whose value is assigned to scores.
%
% NOTES:
% -The mat file you load through uigetfile command should contain: {otr,ova,ttr,tva} variables in
%  order to complete successfully this function. These are training and validation data. If you splitted
%  randomly your Tr subset in Tr and Te subsets (as within sann.m function), and you don't know the exact
%  combination of the splitting of each trained NN component, don't worry because you can use the entire
%  Tr subset from which you splitted Tr and Te.

%% load data
%eval(['load ' mat_file]) 
ttr = T.ttr;
tva = T.tva;
tte = T.tte;
%% 1.\ metric selection
perf_Fnc = {'mse', 'rmse', 'mae', 'mbe2', 'r', 'smape', 'eff', 'D'};
%-------------------------------------------
if nargin > 1
    curr_perf = perf;
else
    curr_perf = 2;
end
%-------------------------------------------

%% 2.\ select the response of 'x' NNs
scores = zeros(size(x,1),1);
ttr = T.ttr;
tva = T.tva;
tte = T.tte;
for jj = 1:size(x,1)
	% select genome of current individual from population:
	genome = x{jj};
	% select the bootstrap replicates to be combined in the GA/PCR method:
	otr = O.otr(:,genome);
	ova = O.ova(:,genome);
	ote = O.ote(:,genome);
%% 3.\ compute ensemble response
	% %  IMPLEMENTED PROCEDURE
	%  a.\ cross-validation weights on Tr [PCRcv]
	[w,~] 			= PCRcv(otr,ttr);
	%  b.\ create BAGNET response on Va (the one used for early stopping criteria!!!!!!!!!!!!!!)
	Yb_va 			= bagnet(ova,w);
	%  c.\ recognize the best PC number on Va
	[~,w_min_va] 	= MSEg(Yb_va,tva);
	%  d.\ create BAGAP response on Te
	Yb_te 			= bagnet(ote,w(:,w_min_va)); 
	%  e.\ performance evaluation on Te
	Ind 			= perfind(Yb_te, tte);

%% 4.\ assign selected metric to scores
	scores(jj) 		= Ind.(perf_Fnc{curr_perf});
end

%% return
end