function mutationChildren = mutationpuntiform_cell(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation,mutationRate,nBoots)
%mutationpuntiform_cell Uniform multi-point mutation.
%   MUTATIONCHILDREN = mutationpuntiform_cell(PARENTS,OPTIONS,GENOMELENGTH,...
%                      FITNESSFCN,STATE,THISSCORE,THISPOPULATION, ...
%                      MUTATIONRATE, nBoots) Creates the mutated children using
%   puntiform mutations at multiple points. Mutated genes are randomly 
%   distributed over the range of genes (1 to nBoots).The new value is NOT a function
%   of the parents value for the gene.
%
%   Example:
%     options = optimoptions('ga','MutationFcn', @mutationpuntiform_cell); 
%
%   This will create an options structure specifying the mutation
%   function to be used is mutationpuntiform_cell.  Since the MUTATIONRATE is
%   not specified, the default value of 1 / length(parents) is used, which ensures
%   that at least 1 gene mutates passing from parent to child.
%
%     mutationRate = 0.05;
%     options = optimoptions('ga','MutationFcn', {@mutationpuntiform_cell, mutationRate});
%
%   This will create an options structure specifying the mutation
%   function to be used is mutationpuntiform_cell and the MUTATIONRATE is
%   user specified to be 0.05.


%% check
if nargin < 8 || isempty(mutationRate)
    mutationRate    = 1 / length(parents); % default mutation rate, at least 1 gene mutates
end

%% main
mutationChildren    = cell(length(parents),1); % zeros(length(parents),GenomeLength);
for i=1:length(parents)
    child           = thisPopulation{parents(i)}; % thisPopulation(parents(i),:);
    mutationPoints  = find(rand(1,length(child)) < mutationRate);
    solved          = false;
    while ~solved
        newGenes    = randsample( nBoots, min(nBoots,numel(mutationPoints)*10) );
        [~,iA,~]    = setxor( newGenes, child );
        if numel(mutationPoints) <= numel(iA)
            solved  = true;
        end
    end
    child(mutationPoints) = newGenes(iA(1:numel(mutationPoints)));
    mutationChildren{i}   = child; % mutationChildren(i,:)
end
    
%% return
end