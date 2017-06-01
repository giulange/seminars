function xoverKids  = crossoverscattered_cell(parents,options,GenomeLength,FitnessFcn,thisScore,thisPopulation)
%CROSSOVERSCATTERED_CELL Position independent crossover function.
%   XOVERKIDS = CROSSOVERSCATTERED_CELL(PARENTS,OPTIONS,GENOMELENGTH, ...
%   FITNESSFCN,SCORES,THISPOPULATION) creates the crossover children XOVERKIDS
%   of the given population THISPOPULATION using the available PARENTS.
%   In single or double point crossover, genomes that are near each other tend
%   to survive together, whereas genomes that are far apart tend to be
%   separated. The technique used here eliminates that effect. Each gene has an
%   equal chance of coming from either parent. This is sometimes called uniform
%   or random crossover.
%
%   Example:
%    Create an options structure using CROSSOVERSCATTERED_CELL as the crossover
%    function
%     options = optimoptions('ga', 'CrossoverFcn', @crossoverscattered);
%    (Note: If calling gamultiobj, replace 'ga' with 'gamultiobj') 

%   Copyright 2003-2015 The MathWorks, Inc.

% How many children to produce?
nKids = length(parents)/2;

% % Extract information about linear constraints, if any
%linCon = options.LinearConstr;
%constr = ~isequal(linCon.type,'unconstrained');

% Allocate space for the kids
xoverKids = cell(nKids,1); %zeros(nKids,GenomeLength);

% To move through the parents twice as fast as thekids are
% being produced, a separate index for the parents is needed
index = 1;
% for each kid...
for i=1:nKids
    % get parents
    r1 = parents(index);
    index = index + 1;
    r2 = parents(index);
    index = index + 1;
    % Randomly select half of the genes from each parent
    % This loop may seem like brute force, but it is twice as fast as the
    % vectorized version, because it does no allocation.
    for j = 1:GenomeLength
        if(rand > 0.5)
            %xoverKids(i,j) = thisPopulation(r1,j);
            xoverKids{i} = thisPopulation{r1};
        else
            %xoverKids(i,j) = thisPopulation(r2,j);
            xoverKids{i} = thisPopulation{r2};
        end
    end

%    % Make sure that offspring are feasible w.r.t. linear constraints
%    if constr
%        feasible = isTrialFeasible(xoverKids(i,:)',linCon.Aineq,linCon.bineq,linCon.Aeq, ...
%                            linCon.beq,linCon.lb,linCon.ub,options.TolCon);
%        if ~feasible % Kid is not feasible
%            % Children are arithmetic mean of two parents (feasible w.r.t
%            % linear constraints)
%            alpha = rand;
%            xoverKids(i,:) = alpha*thisPopulation(r1,:) + ...
%                (1-alpha)*thisPopulation(r2,:);
%        end
%    end
end

%% return
end