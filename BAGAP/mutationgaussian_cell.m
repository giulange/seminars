function mutationChildren = mutationgaussian_cell(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation,scale,shrink)
%mutationgaussian_cell Gaussian mutation.
%   MUTATIONCHILDREN = mutationgaussian_cell(PARENTS,OPTIONS,GENOMELENGTH,...
%   FITNESSFCN,STATE,THISSCORE,THISPOPULATION,SCALE,SHRINK) Creates the
%   mutated children using the Gaussian distribution.
%
%   SCALE controls what fraction of the gene's range is searched. A
%   value of 0 will result in no change, a SCALE of 1 will result in a
%   distribution whose standard deviation is equal to the range of this gene.
%   Intermediate values will produce ranges in between these extremes.
%
%   SHRINK controls how fast the SCALE is reduced as generations go by. A
%   SHRINK value of 0 will result in no shrinkage, yielding a constant
%   search size. A value of 1 will result in SCALE shrinking linearly to 0
%   as GA progresses to the number of generations specified by the options
%   structure. (See 'MaxGenerations' in OPTIMOPTIONS for more details).
%   Intermediate values of SHRINK will produce shrinkage between these
%   extremes. Note: SHRINK may be outside the interval (0,1), but this is
%   ill-advised.
%
%   The arguments to the function are:
%     -parents        : Row vector of parents chosen by the selection function
%     -options        : Options structure
%     -nvars          : Number of variables
%     -FitnessFcn     : Fitness function
%     -state          : Structure containing information about the current generation.
%                       The State Structure describes the fields of state.
%     -thisScore      : Vector of scores of the current population
%     -thisPopulation : Matrix of individuals in the current population

%   Example:
%     options = optimoptions('ga','MutationFcn',{@mutationgaussian});
%
%   This specifies that the mutation function used will be
%   mutationgaussian_cell, and since no values for SCALE or SHRINK are specified
%   the default values are used.
%
%     scale = 0.5; shrink = 0.75;
%     options = optimoptions('ga','MutationFcn',{@mutationgaussian,scale,shrink});
%
%   This specifies that the mutation function used will be
%   mutationgaussian_cell, and the values for SCALE or SHRINK are specified
%   as 0.5 and 0.75 respectively.
%
%   (Note: If calling gamultiobj, replace 'ga' with 'gamultiobj' in the
%   above examples)

% Use default parameters if the are not passed in.
% If these defaults are not what you prefer, you can pass in your own
% values when you set the mutation function:
%
% options.MutationFcn = { mutationgaussian, 0.3, 0} ;
%

%% main
if nargin < 9 || isempty(shrink)
    shrink = 1;
end
if nargin < 8 || isempty(scale)
    scale = 1;
end

if (shrink > 1) || (shrink < 0)
    warning(message('globaloptim:mutationgaussian:shrinkFactor'));
end

scale = scale - shrink * scale * state.Generation/options.Generations;

range = options.PopInitRange;
lower = range(1,:);
upper = range(2,:);
scale = scale * (upper - lower);

mutationChildren = cell(length(parents),1); % zeros(length(parents),GenomeLength);
for i=1:length(parents)
    parent = thisPopulation{parents(i)}; % thisPopulation(parents(i),:);
    mutationChildren{i} = parent  + scale .* randn(1,length(parent)); % mutationChildren(i,:) 
end

%% return
end