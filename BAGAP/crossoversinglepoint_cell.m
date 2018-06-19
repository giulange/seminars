function xoverKids  = crossoversinglepoint_cell(parents,options,GenomeLength,FitnessFcn,thisScore,thisPopulation,useSorting)
%crossoversinglepoint_cell Single point crossover.
%   XOVERKIDS = crossoversinglepoint_cell(PARENTS,OPTIONS,GENOMELENGTH, ...
%   FITNESSFCN,SCORES,THISPOPULATION) creates the crossover children XOVERKIDS
%   of the given population THISPOPULATION using the available parents
%   PARENTS. A single crossover point for each child is chosen at random.
%   The child has the genes of the first parent up to this point and the genes
%   of the other parent after this point.
%
% Note: I can use the crossover single point algorithm (adapted to cell array)
%       after sorting the Genome composition of a parent. So, this crossover
%       function can highlight good or bad bootstrap replicates because they
%       belong to low or high index position within the [1,1000] range, where
%       the first is the worst and the 1000th is the best ANN evaluated on Tr.
%       In fact when useSorting=TRUE, the xOverPoint value prunes out the potentially 
%       worst performing ANN replicates from both parents to build more performing
%       children in this way (parents are both sorted):
%         e.g. length(parent1)=50; xOverPoint = 41; compl = length(parent1)-xOverPoint = 9;
%          > parent1[ 1 ––> xOverPoint | xOverPoint+1 ––> end]
%          > parent2[ 1 ––> compl | compl+1 ––> end]
%          > child[ xOverPoint+1 ––> end | compl+1 ––> end ]
%
%       BUT REMEMBER TO SORT INPUT BOOTSTRAP ANN REPLICATES IN ASCENDING ORDER 
%       OF ERROR PERFORMANCE on Tr (from worst to best performance) BEFORE STARTING
%       Genetic Algorithm analysis.
%
%   Example:
%    Create an options structure using crossoversinglepoint_cell as the crossover
%    function
%     options = optimoptions('ga', 'CrossoverFcn', @crossoversinglepoint);
%    (Note: If calling gamultiobj, replace 'ga' with 'gamultiobj') 


%% Check input arguments
if nargin < 7
    useSorting = false;
end

%% main
% How many children to produce?
nKids = length(parents)/2;

% % Extract information about linear constraints, if any
% linCon = options.LinearConstr;
% constr = ~isequal(linCon.type,'unconstrained');

% Allocate space for the kids
xoverKids = cell(nKids,1); %zeros(nKids,GenomeLength);

% To move through the parents twice as fast as thekids are
% being produced, a separate index for the parents is needed
index = 1;

for i=1:nKids
    % get parents
    %parent1 = thisPopulation(parents(index),:);
    parent1 = thisPopulation{parents(index)};
    index = index + 1;
    %parent2 = thisPopulation(parents(index),:);
    parent2 = thisPopulation{parents(index)};
    index = index + 1;

    % cut point is AFTER this index.
    xOverPoint = ceil(rand * (length(parent1) - 1));

    if useSorting% see the note above to understand this feature!
        parent1 = sort(parent1);
        parent2 = sort(parent2);

        compl = length(parent1) - xOverPoint;

        % make one child
        xoverKids{i} = [ parent1(( xOverPoint + 1 ):  end ); parent2(compl+1:end) ];
    else
        % make one child
        %xoverKids(i,:) = [ parent1(1:xOverPoint),parent2(( xOverPoint + 1 ):  end )  ];
        xoverKids{i} = [ parent1(1:xOverPoint);parent2(( xOverPoint + 1 ):  end )  ];
    end

    % Here I have to avoid duplicates in current child.
    % The procedure takes one random gene from the parent(s) showing the
    % best (lowest) scores, avoiding duplicated values.
    if length(unique(xoverKids{i})) < GenomeLength
        %figure(3),plot(thisScore)
        [~,iScores]  = sort(thisScore,1,'ascend');
        [dupl,Ndupl,pos] = findduplicates( xoverKids{i} );
        for d = 1:numel(dupl)
            for Nd = 1:Ndupl-1
                Found = false;
                % loop each elite parent until duplicates are substituted
                for s = 1:numel(iScores)
                    % select current Elite Parent
                    currEliteParent = thisPopulation{iScores(s)};
                    % Gene selection criterion:
                    %   > randsample,
                    %   > order (because Gene having higher position
                    %            indexes should be more performing).
                    rS = randsample(GenomeLength,GenomeLength);
                    irS = 0;
                    % loop each elite gene untill duplicates are
                    % substituted
                    while ~Found
                        irS = irS + 1;
                        % change Elite Parent when current cannot provide
                        % all required substitutions:
                        if irS>GenomeLength, break, end
                        eliteGene = currEliteParent(rS(irS));
                        if isempty( find(xoverKids{i}==eliteGene,1) )
                            xoverKids{i}(pos{d}(Nd)) = eliteGene;
                            Found = true;
                        end
                    end
                    if Found == true, break, end % break
                end
            end
        end
    end
    if length(unique(xoverKids{i})) < GenomeLength
        error('The customized procedure didn''t work properly!')
    end
end        
%% return
end
%% identify duplicates
function [repeatedValues,numberOfAppearancesOfRepeatedValues,pos] = findduplicates( X )

    uniqueX = unique(X);
    countOfX = hist(X,uniqueX);
    indexToRepeatedValue = (countOfX~=1);
    repeatedValues = uniqueX(indexToRepeatedValue);
    numberOfAppearancesOfRepeatedValues = countOfX(indexToRepeatedValue);
    pos = cell(numel(repeatedValues,1));
    for i = 1:numel(repeatedValues)
        pos{i} = find(X==repeatedValues(i));
    end
    
end
