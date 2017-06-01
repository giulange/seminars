%% PARS
WDIR                = '~/git/seminars/BAGAP';
name                = 'ann1kh2_al05fe_3D__cl_o_r__corr.mat';
save_ga_fil         = 'save_GA_conditions.mat';

% useSorting          = true;
nBoots              = 1000;
% nvars               = 20;  % 25 / 50 / 75 / 100 % == GenomeLength
% popsize             = 50;  % 100 : 100 : 500
% maxgen              = 200; % 250 : 250 : 100*nvars
% mutationRate        = 4 / nvars; % 0.08%
% CrossFract          = 0.8;
FitnessPerf         = 2; % 2:rmse, ...


% usage ??>     WDIR, name, save_ga_fil, useSorting, nBoots, nvars, popsize, maxgen, mutationRate, CrossFract, FitnessPerf)
sensitivity_GA( WDIR, name, save_ga_fil, true,       nBoots, 20,    50,      200,    0.10,         0.80,       FitnessPerf)

