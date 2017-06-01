function pop = create_pop_rand(NVARS,FitnessFcn,options, nBoots)
%create_pop_rand Creates a population of randsamples.
%   POP = create_pop_rand(NVARS,FitnessFcn,options,...)
%	  creates a population of randsampled POP each with
%	  a length of NVARS. 
%
%   The arguments to the function are
%     NVARS: 		Number of variables 
%     FitnessFcn: 	Fitness function 
%     options: 		Options structure used by the GA
%     nBoots:       Number of bootstrap replicates used as experts
%                   to be processed within the fitness function.

totalPopulationSize = sum(options.PopulationSize);
n 					= NVARS;
pop 				= cell(totalPopulationSize,1);

for i = 1:totalPopulationSize
    %pop{i} = randperm(n);
    pop{i} = randsample(nBoots,n);
end

%% return
end