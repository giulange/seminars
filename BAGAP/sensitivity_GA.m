function sensitivity_GA(WDIR,name,save_ga_fil,useSorting,nBoots,nvars,popsize,maxgen,mutationRate,CrossFract,FitnessPerf)
%% P A R s
% WDIR                = '~/git/seminars/BAGAP';
% name                = 'ann1kh2_al05fe_3D__cl_o_r__corr.mat';
% save_ga_fil         = 'save_GA_conditions.mat';
% 
% useSorting          = true;
% nBoots              = 1000;
% nvars               = 20;  % 25 / 50 / 75 / 100 % == GenomeLength
% popsize             = 50;  % 100 : 100 : 500
% maxgen              = 200; % 250 : 250 : 100*nvars
% mutationRate        = 4 / nvars; % 0.08%
% CrossFract          = 0.8;
% FitnessPerf         = 2; % 2:rmse, ...
% 
% % migrint           = 5;   % migration interval [generations] | USELESS !!

%% L O A D
cd( WDIR )
load( fullfile(WDIR,name) )

%load( fullfile(WDIR,'ann1kh2_al05fe_3D__cl_o_r__corr__SMALL.mat') )
%load( fullfile(WDIR,'data4models.mat'), 'DATA2' )
%% -- prepare
sclorpt         = name( strfind(name,'__')+2 : strfind(name,'.')-1 );
annResamples    = numel(iPerf_tr);
[R,Q]           = size(x_in);
Ntr             = size(y_sim_tr,2);
Nva             = size(y_sim_va,2);
Nte             = size(y_sim_te,2);
Ntot            = Ntr + Nva + Nte;
%% S O R T  | by Tr performance [rmse]
if useSorting
    [~,iTr] = sort([iPerf_tr.rmse],2,'descend');
    figure(9),clf
    yyaxis left
    plot(1:numel(iTr),[iPerf_tr(iTr).rmse],'LineWidth',2),ylabel('RMSE')
    yyaxis right
    plot(1:numel(iTr),[iPerf_tr(iTr).r],'LineWidth',2),ylabel('r')
    legend('rmse','r')
else
    iTr = 1:P.annResamples;
end
%% DATA EXTRACTION | applying rmse sorting
% simulations:
O.otr           = y_sim_tr(iTr,:)';
O.ova           = y_sim_va(iTr,:)';
O.ote           = y_sim_te(iTr,:)';
% targets:
T.ttr           = ttr';
T.tva           = tva';
T.tte           = tte';
%% --- see
% size of subsets
fprintf(' DATA   [%4s,%4s]\n','rows','cols')
fprintf('--------------------\n')
fprintf(' otr    [%4d,%4d]\n',size(O.otr))
fprintf(' ova    [%4d,%4d]\n',size(O.ova))
fprintf(' ote    [%4d,%4d]\n',size(O.ote))
fprintf(' ttr    [%4d,%4d]\n',size(T.ttr))
fprintf(' tva    [%4d,%4d]\n',size(T.tva))
fprintf(' tte    [%4d,%4d]\n',size(T.tte))
fprintf('--------------------\n')
fprintf('        [%4s,%4s]\n','repl','boot')
%% Repeatability
seedN = randi( [1,10^6], 1, 1 );
rng( seedN, 'twister' )
save_rng = rng;
%% FUNCTIONS
%% -- create pop
CreateFnc = @(nvars,fitfnc,opt) create_pop_rand(nvars,fitfnc,opt,nBoots);
%% -- crossover
% https://it.mathworks.com/help/gads/vary-mutation-and-crossover.html?searchHighlight=ga%20crossover&s_tid=doc_srchtitle
CrossFnc = @(parents,options,GenomeLength,FitnessFcn,thisScore,thisPopulation) ...
            crossoversinglepoint_cell(parents,options,GenomeLength,FitnessFcn,thisScore,thisPopulation,useSorting); % *** advised!
%CrossFnc = @crossoverscattered_cell;

% Note: I can use the crossover single point algorithm (adapted to cell array)
%       after sorting the Genome composition of a parent. So, this crossover
%       function can highlight good or bad bootstrap replicates because they
%       belong to low or high index position within the [1,1000] range, where
%       the first is the worst ANN and the 1000th is the best evaluated on Tr.
%       In fact when useSorting=TRUE, the xOverPoint value prunes out the worst
%       performing ANN replicates from both parents to build more performing
%       children in this way (parents are both sorted):
%         e.g. length(parent1)=50; xOverPoint = 41; compl = length(parent1)-xOverPoint = 9;
%          > parent1[ 1 ??> xOverPoint | xOverPoint+1 ??> end]
%          > parent2[ 1 ??> compl | compl+1 ??> end]
%          > child[ xOverPoint+1 ??> end | compl+1 ??> end ]
%
%       BUT REMEMBER TO SORT INPUT BOOTSTRAP ANN REPLICATES IN ASCENDING ORDER 
%       OF ERROR PERFORMANCE (from worst to best performance) BEFORE STARTING
%       Genetic Algorithm analysis.
%% -- mutation
MutateFnc = @(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation)...
             mutationpuntiform_cell(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation,mutationRate,nBoots);
%% -- fitness
FitnessFnc = @(x) ga_pcr_fitness(x, FitnessPerf, O, T);
%% -- plot | incomplete find/develop something more interesting! | you can plot more fncs
% purge = 1;
% PltFnc = @(options,state,flag) gaplotgenealogy(options,state,flag,purge);

% PltFnc = @gaplotselection;
PltFnc = '';
%% G.A.
%% -- set options
% do not use gaoptimset but optimoptions
%   |NO| ga_opt = gaoptimset(@ga);
ga_opt = optimoptions(@ga, ...
                        'PopulationType',         'custom',...
                        'InitialPopulationRange', [1;popsize]);

ga_opt = optimoptions(ga_opt, ...
                        'CreationFcn',            CreateFnc, ...
                        'CrossoverFcn',           CrossFnc,...
                        'MutationFcn',            MutateFnc, ...
                        'PlotFcns',               PltFnc, ...
                        'PopulationSize',         popsize, ...
                        'MaxGenerations',         maxgen,...
                        'MaxStallGenerations',    200,...
                        'UseVectorized',          true,...
                        'UseParallel',            false...
                      );
%% -- run
[x,fval,reason,output] = ga(FitnessFnc,nvars,[],[],[],[],[],[],[],ga_opt);
%% ---- eval exitflag
% Possible values of EXITFLAG and the corresponding exit conditions are:
%
%     1 Average change in value of the fitness function over
%        options.MaxStallGenerations generations less than 
%        options.FunctionTolerance and constraint violation less than 
%        options.ConstraintTolerance.
%     3 The value of the fitness function did not change in
%        options.MaxStallGenerations generations and constraint violation 
%        less than options.ConstraintTolerance.
%     4 Magnitude of step smaller than machine precision and constraint
%        violation less than options.ConstraintTolerance. This exit 
%        condition applies only to nonlinear constraints.
%     5 Fitness limit reached and constraint violation less than
%        options.ConstraintTolerance. 
%     0 Maximum number of generations exceeded.
%    -1 Optimization terminated by the output or plot function.
%    -2 No feasible point found.
%    -4 Stall time limit exceeded.
%    -5 Time limit exceeded.
eFlags = {
    1, 'Average change in value of the fitness function over options.MaxStallGenerations generations less than options.FunctionTolerance and constraint violation less than options.ConstraintTolerance.';
    3, 'The value of the fitness function did not change in options.MaxStallGenerations generations and constraint violation less than options.ConstraintTolerance.';
    4, 'Magnitude of step smaller than machine precision and constraint violation less than options.ConstraintTolerance. This exit condition applies only to nonlinear constraints.';
    5, 'Fitness limit reached and constraint violation less than options.ConstraintTolerance.';
    0, 'Maximum number of generations exceeded.';
   -1, 'Optimization terminated by the output or plot function.';
   -2, 'No feasible point found.';
   -4, 'Stall time limit exceeded.';
   -5, 'Time limit exceeded.';
   };
%% -- eval results
pearson_r   = ga_pcr_fitness(x, 5, O, T);
rmse        = ga_pcr_fitness(x, 2, O, T);
eFlag_str   = eFlags( cell2mat(eFlags(:,1))==reason, 2 );
%% -- print
% PRINT ON SCREEN
fprintf('%s\n',repmat('_',1,110))
fprintf('Sorted iANN:\t')
fprintf('%d,', sort(x{1})')
fprintf('\n')
% useSorting,(nBoots),nvars,popsize,maxgen,mutationRate,CrossFract,(FitnessPerf)
fprintf('useSorting: %2d,  nvars: %3d,  popsize: %4d,  maxgen: %5d,  mutationRate: %.3f,  CrossFract: %.3f.\n', ...
        useSorting,nvars,popsize,maxgen,mutationRate,CrossFract)
fprintf('r: %.3f,  rmse: %.3f.\n\n',pearson_r,rmse)

% PRINT IN FILE:
fid = fopen( fullfile(WDIR,'ga_runs_results.txt'),'a');
fprintf(fid,'%s\n',repmat('_',1,110));
fprintf(fid,'Sorted iANN:\t');
fprintf(fid,'%d,', sort(x{1})');
fprintf(fid,'\n');
% useSorting,(nBoots),nvars,popsize,maxgen,mutationRate,CrossFract,(FitnessPerf)
fprintf(fid,'useSorting: %2d,  nvars: %3d,  popsize: %4d,  maxgen: %5d,  mutationRate: %.3f,  CrossFract: %.3f.\n', ...
        useSorting,nvars,popsize,maxgen,mutationRate,CrossFract);
fprintf(fid,'r: %.3f,  rmse: %.3f.\n\n',pearson_r,rmse);
fclose(fid);
%% -- save conditions
GA.name              = name;
GA.useSorting        = useSorting;
GA.nBoots            = nBoots;
GA.nvars             = nvars;
GA.popsize           = popsize;
GA.maxgen            = maxgen;
GA.mutationRate      = mutationRate;
GA.CrossFract        = CrossFract;
GA.FitnessPerf       = FitnessPerf;
GA.save_rng          = save_rng;
GA.pearson_r         = pearson_r;
GA.rmse              = rmse;
GA.x                 = x;
GA.fval              = fval;
GA.reason            = reason;
GA.eFlag_str         = eFlag_str;
GA.output            = output;
GA.ga_opt            = ga_opt;
GA.BAGAP             = P;

if exist( fullfile(WDIR,save_ga_fil), 'file' )
    load( fullfile(WDIR,save_ga_fil), 'GA' )
    GA(end+1) = GA;
    save( fullfile(WDIR,save_ga_fil), 'GA' )
else
    save( fullfile(WDIR,save_ga_fil), 'GA' )
end