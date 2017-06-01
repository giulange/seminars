%% log file for paper "clay_model_comparison" [created on 29/APR/2009]
% this log file is written after EGU conference in Vienna. Here the sann macro is modified in a more
% performing function in which the following parameters are implemented:
% early stopping criteria, activation functions selection, training function selection, network training
% parameters setting, normalization.

%% OLD --> SELECT PEDOLOGICAL SUPPORT, THEN LOAD AND PREPARE DATA || see sannd function
clear
%***************** SETTINGS *****************
%-select pedological support:
%--1:Topsoil, 2:Profile, 3:Horizon
curr_support = 3;
%-network settings
net_set.epochs = 60;
net_set.show = NaN;
%-number of simulations
nSim = 1000;
%***************** SETTINGS *****************

cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
% ---------------- SUPPORT ----------------
% ROW 1: tables; ROW 2: split; ROW 3: predictors; ROW 4: target; ROW 5: p1; ROW 6: p2; ROW 7: p3.
support = {'Topsoil', 'Profile', 'Horizon'; 'topsoil','profile','horizon'; 6:10,5:9,9:14; 5,4,8; ...
           1:3,1:3,2:4; 2:4,2:4,1:5; 1:4,1:4,1:6};
% ---------------- SUPPORT ----------------

% ---------------- LOAD ----------------
eval(['load Data.mat ' support{1,curr_support} ' F split'])
%-predictors...
eval(['p = cell2mat(' support{1,curr_support} '(2:end,support{3,curr_support}))'';'])
eval(['h_p = ' support{1,curr_support} '(1,support{3,curr_support});'])
%...and target
eval(['t = cell2mat(' support{1,curr_support} '(2:end,support{4,curr_support}))'';'])
eval(['h_t = ' support{1,curr_support} '(1,support{4,curr_support});'])
% ---------------- LOAD ----------------

% ---------------- SUBSETTING ----------------
% create Tr and Va subsets
ptr = p(support{5,curr_support}, split.(support{2,curr_support}){1});
h_p(support{5,curr_support})
ttr = t(1, split.(support{2,curr_support}){1});
pva = p(support{5,curr_support}, split.(support{2,curr_support}){2});
tva = t(1, split.(support{2,curr_support}){2});
% ---------------- SUBSETTING ----------------

%% SIMULATIONS
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
%% ---- SUPPORTS: top:1 -- pro:2 -- hor:3
% This file 
%% (1)-------- {p1, 3:1, VV='n', 'tt', 'lm', net_set=[60, NaN], norm='y');
clear;clc
%---------------------------------------------
curr_support = 3;   % top:1 -- pro:2 -- hor:3
%---------------------------------------------

curr_data = 1;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva,support,net_set,nSim] = sannd(curr_support, curr_data, 60, NaN, 1000); %(support,data,epochs,show,nSim)

for i=1:nSim
    progress_bar(i, nSim)
    %                                |_____data_______|_ES_|__AF__|_TF_|_netset_|_norm_| 
    [n{i},trainRec{i},PERF(i)] = sann(ptr,ttr, pva,tva, 'n', 'tt' ,'lm', net_set,  'y' );
end

eval(['save sann_output\' support{2,curr_support} '_#1.mat n trainRec PERF'])

%% (2)-------- {p1, 3:3:1, VV='n', 'ttt', 'lm', net_set=[60, NaN], norm='y');
clear;clc
%---------------------------------------------
curr_support = 3;   % top:1 -- pro:2 -- hor:3
%---------------------------------------------

curr_data = 1;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva,support,net_set,nSim] = sannd(curr_support, curr_data, 60, NaN, 1000); %(support,epochs,show,nSim)

for i=1:nSim
    progress_bar(i, nSim)
    %                                |_____data_______|_ES_|__AF__|_TF_|_netset_|_norm_| 
    [n{i},trainRec{i},PERF(i)] = sann(ptr,ttr, pva,tva, 'n', 'ttt' ,'lm', net_set,  'y' );
end

eval(['save sann_output\' support{2,curr_support} '_#2.mat n trainRec PERF'])

%% (3)-------- {p1, 3:1, VV='y', 'tt', 'lm', net_set=[60, NaN], norm='y');
clear;clc
%---------------------------------------------
curr_support = 3;   % top:1 -- pro:2 -- hor:3
%---------------------------------------------

curr_data = 1;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva,support,net_set,nSim] = sannd(curr_support, curr_data, 60, NaN, 1000); %(support,epochs,show,nSim)

for i=1:nSim
    progress_bar(i, nSim)
    %                                |_____data_______|_ES_|__AF__|_TF_|_netset_|_norm_| 
    [n{i},trainRec{i},PERF(i)] = sann(ptr,ttr, pva,tva, 'y', 'tt' ,'lm', net_set,  'y' );
end

eval(['save sann_output\' support{2,curr_support} '_#3.mat n trainRec PERF'])

%% (4)-------- {p1, 3:3:1, VV='y', 'ttt', 'lm', net_set=[60, NaN], norm='y');
clear;clc
%---------------------------------------------
curr_support = 3;   % top:1 -- pro:2 -- hor:3
%---------------------------------------------

curr_data = 1;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva,support,net_set,nSim] = sannd(curr_support, curr_data, 60, NaN, 1000); %(support,epochs,show,nSim)

for i=1:nSim
    progress_bar(i, nSim)
    %                                |_____data_______|_ES_|__AF__|_TF_|_netset_|_norm_| 
    [n{i},trainRec{i},PERF(i)] = sann(ptr,ttr, pva,tva, 'y', 'ttt' ,'lm', net_set,  'y' );
end

eval(['save sann_output\' support{2,curr_support} '_#4.mat n trainRec PERF'])

%% (5)-------- {p2, 3:1, VV='y', 'tt', 'lm', net_set=[60, NaN], norm='y');
clear;clc
%---------------------------------------------
curr_support = 3;   % top:1 -- pro:2 -- hor:3
%---------------------------------------------

curr_data = 2;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva,support,net_set,nSim] = sannd(curr_support, curr_data, 60, NaN, 1000); %(support,epochs,show,nSim)

for i=1:nSim
    progress_bar(i, nSim)
    %                                |_____data_______|_ES_|__AF__|_TF_|_netset_|_norm_| 
    [n{i},trainRec{i},PERF(i)] = sann(ptr,ttr, pva,tva, 'y', 'tt' ,'lm', net_set,  'y' );
end

eval(['save sann_output\' support{2,curr_support} '_#5.mat n trainRec PERF'])

%% (6)-------- {p2, 3:3:1, VV='y', 'ttt', 'lm', net_set=[60, NaN], norm='y');
clear;clc
%---------------------------------------------
curr_support = 3;   % top:1 -- pro:2 -- hor:3
%---------------------------------------------

curr_data = 2;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva,support,net_set,nSim] = sannd(curr_support, curr_data, 60, NaN, 1000); %(support,epochs,show,nSim)

for i=1:nSim
    progress_bar(i, nSim)
    %                                |_____data_______|_ES_|__AF__|_TF_|_netset_|_norm_| 
    [n{i},trainRec{i},PERF(i)] = sann(ptr,ttr, pva,tva, 'y', 'ttt' ,'lm', net_set,  'y' );
end

eval(['save sann_output\' support{2,curr_support} '_#6.mat n trainRec PERF'])

%% (7)-------- {p3, 3:1, VV='y', 'tt', 'lm', net_set=[60, NaN], norm='y');
clear;clc
%---------------------------------------------
curr_support = 3;   % top:1 -- pro:2 -- hor:3
%---------------------------------------------

curr_data = 3;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva,support,net_set,nSim] = sannd(curr_support, curr_data, 60, NaN, 1000); %(support,epochs,show,nSim)

for i=1:nSim
    progress_bar(i, nSim)
    %                                |_____data_______|_ES_|__AF__|_TF_|_netset_|_norm_| 
    [n{i},trainRec{i},PERF(i)] = sann(ptr,ttr, pva,tva, 'y', 'tt' ,'lm', net_set,  'y' );
end

eval(['save sann_output\' support{2,curr_support} '_#7.mat n trainRec PERF'])

%% (8)-------- {p3, 3:3:1, VV='y', 'ttt', 'lm', net_set=[60, NaN], norm='y');
clear;clc
%---------------------------------------------
curr_support = 3;   % top:1 -- pro:2 -- hor:3
%---------------------------------------------

curr_data = 3;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva,support,net_set,nSim] = sannd(curr_support, curr_data, 60, NaN, 1000); %(support,epochs,show,nSim)

for i=1:nSim
    progress_bar(i, nSim)
    %                                |_____data_______|_ES_|__AF__|_TF_|_netset_|_norm_| 
    [n{i},trainRec{i},PERF(i)] = sann(ptr,ttr, pva,tva, 'y', 'ttt' ,'lm', net_set,  'y' );
end

eval(['save sann_output\' support{2,curr_support} '_#8.mat n trainRec PERF'])

%% (9)-------- {p3, 3:1, VV='n', 'tt', 'lm', net_set=[60, NaN], norm='y');
clear;clc
%---------------------------------------------
curr_support = 3;   % top:1 -- pro:2 -- hor:3
%---------------------------------------------

curr_data = 3;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva,support,net_set,nSim] = sannd(curr_support, curr_data, 60, NaN, 1000); %(support,epochs,show,nSim)

for i=1:nSim
    progress_bar(i, nSim, 'sann.m')
    %                                |_____data_______|_ES_|__AF__|_TF_|_netset_|_norm_| 
    [n{i},trainRec{i},PERF(i)] = sann(ptr,ttr, pva,tva, 'n', 'tt' ,'lm', net_set,  'y' );
end

eval(['save sann_output\' support{2,curr_support} '_#9.mat n trainRec PERF'])

%% (10)-------- {p3, 3:3:1, VV='n', 'ttt', 'lm', net_set=[60, NaN], norm='y');
clear;clc
%---------------------------------------------
curr_support = 3;   % top:1 -- pro:2 -- hor:3
%---------------------------------------------

curr_data = 3;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva,support,net_set,nSim] = sannd(curr_support, curr_data, 60, NaN, 1000); %(support,epochs,show,nSim)

for i=1:nSim
    progress_bar(i, nSim)
    %                                |_____data_______|_ES_|__AF__|_TF_|_netset_|_norm_| 
    [n{i},trainRec{i},PERF(i)] = sann(ptr,ttr, pva,tva, 'n', 'ttt' ,'lm', net_set,  'y' );
end

eval(['save sann_output\' support{2,curr_support} '_#10.mat n trainRec PERF'])

%% ---- BOXPLOT
% Here I compute the boxplot of all points

%------------------------------------------------
curr_support = 3;   % top:1 -- pro:2 -- hor:3
curr_perf = 4;
fontsize = 14;
linewidth = 2;
rotation = 0; % {0, 90}
Hor_Align = 'Left';     % { 'Left', 'Center', 'Right' }
Ver_Align = 'Middle';   % { 'Middle', 'Top', 'Bottom', ... }
%------------------------------------------------
support = {'Topsoil', 'Profile', 'Horizon'; 'topsoil','profile','horizon'; 'top','prof','hor'; ...
           'REGR', 'REGR', 'GLM'; 'MOCOK', 'MOCOK', 'OK'; ; ...
           'regr','regr','glm_all'; 'mocok','mocok','ok'; '7','6','9'};
perf = {'rmse', 'r', 'mbe2', 'smape', 'eff', 'D'; ...
        'RMSE', 'Pearson''s r', 'MBE', 'SMAPE', 'Efficiency', 'Willmott''s D'};

cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab\sann_output')
bplot_str = '[';
for i = 1:10
    eval(['load ' support{2,curr_support} '_#' num2str(i) '.mat PERF' ])
    eval(['perf_' num2str(i) ' = PERF;'])
    clear PERF
    bplot_str = [bplot_str ' [perf_' num2str(i) '.' (perf{1,curr_perf}) ']'''];
end
bplot_str = [bplot_str ' ]'];
clear i

% lines: glm, ok
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
load results ind*
figure
hold on
    %lab = {{'p1';'tt'}, {'p1';'ttt'}, {'p1';'tt';'ES'}, {'p1';'ttt';'ES'}, {'p2';'tt';'ES'}, {'p2';'ttt';'ES'},  {'p3';'tt';'ES'}, {'p3';'ttt';'ES'} };
    %      1        2         3           4            5           6            7           8
    lab = {'p1-tt', 'p1-ttt', 'p1-tt-ES', 'p1-ttt-ES', 'p2-tt-ES', 'p2-ttt-ES', 'p3-tt-ES', 'p3-ttt-ES', 'p3-tt', 'p3-ttt'};
    eval(['boxplot(' bplot_str ', ''labels'',lab, ''notch'',''on'')'])
    ylabel(perf{2,curr_perf}, 'FontSize',fontsize, 'FontWeight','b')
    %GLM
    eval(['s(1:size(lab,2)+2) = ind_' support{3,curr_support} '_' support{6,curr_support} '.(perf{1,curr_perf});'])
    line([0:size(s,2)-1],s,'LineStyle','--', 'LineWidth',linewidth, 'Color','k')
    text(size(lab,2)+0.5,s(1), support{4,curr_support}, 'Color','k', 'FontSize',fontsize, 'Rotation',rotation, 'VerticalAlignment',Ver_Align, 'HorizontalAlignment',Hor_Align)
    %KRIGING
    eval(['s(1:size(lab,2)+2) = ind_' support{3,curr_support} '_' support{7,curr_support} '.(perf{1,curr_perf});'])
    line([0:size(s,2)-1],s,'LineStyle','--', 'LineWidth',linewidth, 'Color','b')
    text(size(lab,2)+0.5,s(1), support{5,curr_support}, 'Color','b', 'FontSize',fontsize, 'Rotation',rotation, 'VerticalAlignment',Ver_Align, 'HorizontalAlignment',Hor_Align)
    %SMU
    eval(['s(1:size(lab,2)+2) = ind_' support{3,curr_support} '_udp.(perf{1,curr_perf});'])
    line([0:size(s,2)-1],s,'LineStyle','--', 'LineWidth',linewidth, 'Color','r')
    text(size(lab,2)+0.5,s(1), 'SMU', 'Color','r', 'FontSize',fontsize, 'Rotation',rotation, 'VerticalAlignment',Ver_Align, 'HorizontalAlignment',Hor_Align)
    if curr_perf == 1 %because the metric selected for within the ga_sel_pcr function was RMSE for all supports
        eval( ['load GA\ga_' support{3,curr_support} '#' support{8,curr_support} '_bag100rand_RUN#4 fval'])
        eval(['s(1:size(lab,2)+2) = fval;'])
        line([0:size(s,2)-1],s,'LineStyle','--', 'LineWidth',linewidth, 'Color','m')
        text(size(lab,2)+0.5,s(1), 'BAGGING', 'Color','m', 'FontSize',fontsize, 'Rotation',rotation, 'VerticalAlignment',Ver_Align, 'HorizontalAlignment',Hor_Align)
    end
hold off

%% -------- SAVE boxplot
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
eval(['print -dtiff -r600 ' support{3,curr_support} '_' perf{1,curr_perf}])

%% ---- scatter ANN
clear, clc
%------------------------------------------------
curr_support = 1;   % top:1 -- pro:2 -- hor:3
protocol = 7;       % 1 to 10   {top=#7, prof=#6, hor=#9}
nn_boot = 'min';    % {'min' or 1 to 1000}
fontsize = 14;
linewidth = 2;
%------------------------------------------------

support = {'topsoil','profile','horizon'; 'top','prof','hor'};
% cd('J:\work\MyPapers\clay_model_comparison\MatLab')
eval(['load sann_output',filesep,support{1,curr_support} '_#' num2str(protocol) '.mat PERF'])
eval(['load results clay_' support{2,curr_support} '_mea'])
if nn_boot == 'min'
    col = find([PERF.rmse]==min([PERF.rmse]));
else
    col = nn_boot;
end
eval(['scatter(PERF(col).pred'', clay_' support{2,curr_support} '_mea)'])

%% ------- scatter mANN
clear, clc
%------------------------------------------------
curr_support = 1;   % top:1 -- pro:2 -- hor:3
protocol = 7;       % 1 to 10   {top=#7, prof=#6, hor=#9}
nn_boot = 'median'; % {'median'}
fontsize = 14;
linewidth = 2;
%------------------------------------------------

support = {'topsoil','profile','horizon'; 'top','prof','hor'};
% cd('J:\work\MyPapers\clay_model_comparison\MatLab')
eval(['load sann_output',filesep,support{1,curr_support} '_#' num2str(protocol) '.mat PERF'])
eval(['load results clay_' support{2,curr_support} '_mea'])
if nn_boot == 'median'
    col = find([PERF.rmse]==min([PERF.rmse]));
else
    col = nn_boot;
end
eval(['scatter(PERF(col).pred'', clay_' support{2,curr_support} '_mea)'])

%% ---- map of clay predictions: Topsoil #7 best ANN
% topsoil_#7 has the following characteristics:
% -topology         4:1
% -input type       p3
% -I/O data are of type (nvars,npoints)

load /media/Vista/work/MyPapers/clay_model_comparison/MatLab/sann_output/topsoil_#7.mat
rmse = [PERF.rmse];
% % plot(sort(rmse))
% % min(rmse)   % note that median rmse for topsoil #7 is about 8 (see top_rmse.tif)

% best ANN within 1000 NNs, based on rmse
Nb = n{ rmse==min(rmse) };
% % Nb.inputs{:}

% ================ GRID ================
% retrive input GRID data
%   In sannd I found that:
% % varargout{5} = {'Topsoil', 'Profile', 'Horizon', 'tables'; 'topsoil','profile','horizon','split'; 6:10,5:9,9:14,'predictors'; 5,4,8,'target'; ...
% %            1:3,1:3,2:4,'p1'; 2:4,2:4,1:5,'p2'; 1:4,1:4,1:6,'p3'};
%   ...so predictors are extracted from Topsoil table in Data.mat
load /media/Vista/work/MyPapers/clay_model_comparison/MatLab/Data.mat Topsoil
%   ...and p3 data are columns 6:9 of Topsoil table, that is
% Topsoil(1,6:9)  % {'d_tel','il_spi_tel','ndvi_sum_8-12','ndvi_maxdiff_1-16';}

% load GRIDs
% --d_tel
str_dem = '/media/DATI/work_b/Dottorato/PEDOMETRICS/DATA/UOT/RegrKr/Layers/Raster/d_tel.asc';
[dem,demH] = ImportAsciiRaster(str_dem);
% figure(1);subplot(221);mapshow(dem,header2R(demH),'displaytype','surface'); axis tight
% --il_spi_tel
str_spi = '/media/DATI/work_b/Dottorato/PEDOMETRICS/DATA/UOT/RegrKr/Layers/ILWIS/il_spi_tel.asc';
[spi,spiH] = ImportAsciiRaster(str_spi);
% subplot(222);mapshow(spi,header2R(spiH),'displaytype','surface');axis tight
% --NDVIS5
str_ndvis5 = '/media/DATI/work_b/Datasets/NDVI/Italia meridionale 2004/Modis/NDVIS5.asc';
[ndvis5,ndvis5H] = ImportAsciiRaster(str_ndvis5);
% subplot(223);mapshow(ndvis5,header2R(ndvis5H),'displaytype','surface');axis tight
% --NDVID16
str_ndvid16 = '/media/DATI/work_b/Datasets/NDVI/Italia meridionale 2004/Modis/NDVID16.asc';
[ndvid16,ndvid16H] = ImportAsciiRaster(str_ndvid16);
% subplot(224);mapshow(ndvid16,header2R(ndvid16H),'displaytype','surface');axis tight
% --compare headers for consistency
% [demH,spiH,ndvis5H,ndvid16H]
% NOTE: there is 1 meter discrepancy in spiH, but I ignore since error is
% very very small!!
% ================ GRID ================

% create ptr for simulation: that is to say "psi"
psi = [dem(:),spi(:),ndvis5(:),ndvid16(:)]';

% clear
clear PERF ans dem ndvis5 ndvid16 spi n rmse str* trainRec varargout

% normalize GRIDs (see 'normalization' section in sann.m)
% p = cell2mat(Topsoil(2:end,6:9));
% t = cell2mat(Topsoil(2:end,5));
% mm_p = [max(p)' min(p)'];     % I cannot use Topsoil as source of
%                                 interval range! So I use the entire study
%                                 area with psi!
mm_p = [max(psi'); min(psi')]';
cd /media/Vista/work/MyMatLab/ANN/  % to use my normalize function
psin = normalize(psi, mm_p,-1,+1);
cd /media/Vista/work/MyPapers/clay_model_comparison/MatLab/

% SIMULATION PHASE
map_annb = sim(Nb,psin);
mm_t = [max(cell2mat(Topsoil(2:end,5))), min(cell2mat(Topsoil(2:end,5)))];
map_annb = denormalize(map_annb, mm_t,-1,1);
map_annb = reshape(map_annb,demH(2),demH(1));
mapshow(map_annb,header2R(demH),'DisplayType','surface'), axis tight

% SAVE
annbH = demH;
% save /media/Vista/work/MyPapers/clay_model_comparison/MatLab/maps-Topsoil.mat map_annb -append
% save /media/Vista/work/MyPapers/clay_model_comparison/MatLab/maps-Topsoil.mat annbH -append
%% =================================================================
%% ==============           SELECTION           ====================
%% Genetic Algorithm: select proper NNS from the 1000 used
%% ---- topsoil
%% -------- prepare data || protocol ==> {p3, 'tt', ES}
cd('J:\work\MyPapers\clay_model_comparison\MatLab')
%-PREPARE tr AND va SUBSETS
curr_support = 1;   % top:1 -- pro:2 -- hor:3
curr_data = 3;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva] = sannd(curr_support, curr_data);
clear curr_support curr_data

%-NORMALIZE
mm_p = [max([ptr,pva]')' min([ptr,pva]')'];
mm_t = [max([ttr,tva]')' min([ttr,tva]')'];
ptr = normalize(ptr, mm_p,-1,+1);
ttr = normalize(ttr, mm_t,-1,+1);
pva = normalize(pva, mm_p,-1,+1);
tva = normalize(tva, mm_t,-1,+1);

% ============   LOOP   ============
nSim = 200; %number of simulations
perc_te = 0.125; %percentage to assign to Te subset
TF = 'lm'; %training function to be used
%-initialization of output simulation arrays
otr = zeros(round(size(ptr,2)*(1-perc_te)),nSim);
ote = zeros(round(size(ptr,2)*perc_te),nSim);
ova = zeros(size(pva,2),nSim);
%-loop
for i = 1:nSim
    progress_bar(i, nSim);
    %-Te
    F_te(1,:) = randsample(size(ptr,2), round(size(ptr,2)*perc_te));
    F_tr = setxor(1:size(ptr,2),F_te);
    pte = ptr(:,F_te);
    tte = ttr(:,F_te);
    ptrs = ptr(:,F_tr); %small ptr
    ttrs = ttr(:,F_tr); %small ttr
    %-net
    %     newff_stack(in_range, topo, transferF, initial_weights);
    net = newff_stack([-1 +1], [size(ptrs,1) 1], 'tt', 1);
    net.(TF).trainParam.epochs = 60;
    net.(TF).trainParam.show = NaN;
    %-training
    TV.P = pva; TV.T = tva; %for estimating net ability to generalize
    VV.P = pte; VV.T = tte; %for early stopping criteria
    [n{i},trainRec{i}] = train(net.(TF), ptrs,ttrs,[],[],VV,TV);
    %-simulation
    Otr = sim(n{i},ptrs);
    Ote = sim(n{i},pte);
    Ova = sim(n{i},pva);
    Otr = denormalize(Otr, mm_t, -1, 1);
    Ote = denormalize(Ote, mm_t, -1, 1);
    Ova = denormalize(Ova, mm_t, -1, 1);
    otr(:,i) = Otr;
    ote(:,i) = Ote;
    ova(:,i) = Ova;
end
ttr = denormalize(ttrs, mm_t,-1,+1);
tte = denormalize(tte, mm_t,-1,+1);
tva = denormalize(tva, mm_t,-1,+1);
ttr = ttr'; tte = tte'; tva = tva';
% ============   LOOP   ============

%-SAVE
save GA\ga_top#7_bagnet200.mat otr ote ova ttr tte tva n trainRec

%% ------------ REDUCE at RANDOM

%------------------------------
F_size = 100; % {25, 50, 100}
%------------------------------

cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
load GA\ga_top#7_bagnet200
F = randsample(size(otr,2), F_size);
otr = otr(:,F);
ote = ote(:,F);
ova = ova(:,F);
n = n(:,F);
trainRec = trainRec(:,F);
clear F
eval(['save GA\ga_top#7_bagnet' num2str(F_size) '_rand.mat otr ote ova ttr tte tva n trainRec'])

%% -------- run genetic algorithm computations
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab\GA')
%   I created a ga option structure array from gatool and exported to 'ga_opt_ann_sel_clay_top'.
load GA\ga_options 
save GA\ga_options 
clear

%% ------------ random 100
clear
clc
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab\GA')
tic

% 1.\
disp('FIRST run')
load GA\ga_options ga_opt_ann_sel_clay_top
[x fval exitflag output population1 scores] = ga(@ga_sel_pcr, 100, [],[],[],[],[],[],[], ga_opt_ann_sel_clay_top);
elapsed_time = toc;
save ga_top#7_bag100rand_RUN#1.mat
disp([ 'Time:     ' num2str(elapsed_time) ]);
disp([ 'Fitness:  ' num2str(fval) ]);
clear

% 2.\
disp('SECOND run')
load GA\ga_options ga_opt_ann_sel_clay_top
load ga_top#7_bag100rand_RUN#1 population1
ga_opt_ann_sel_clay_top = gaoptimset(ga_opt_ann_sel_clay_top, 'InitialPop',int8(population1));
[x fval exitflag output population2 scores] = ga(@ga_sel_pcr, 100, [],[],[],[],[],[],[], ga_opt_ann_sel_clay_top);
elapsed_time = toc;
save ga_top#7_bag100rand_RUN#2.mat
disp([ 'Time:     ' num2str(elapsed_time) ]);
disp([ 'Fitness:  ' num2str(fval) ]);
clear

% 3.\
disp('THIRD run')
load GA\ga_options ga_opt_ann_sel_clay_top
load ga_top#7_bag100rand_RUN#2 population2
ga_opt_ann_sel_clay_top = gaoptimset(ga_opt_ann_sel_clay_top, 'InitialPop',int8(population2));
[x fval exitflag output population3 scores] = ga(@ga_sel_pcr, 100, [],[],[],[],[],[],[], ga_opt_ann_sel_clay_top);
elapsed_time = toc;
save ga_top#7_bag100rand_RUN#3.mat
disp([ 'Time:     ' num2str(elapsed_time) ]);
disp([ 'Fitness:  ' num2str(fval) ]);
clear

% 4.\
disp('FOURTH run')
load GA\ga_options ga_opt_ann_sel_clay_top
load ga_top#7_bag100rand_RUN#3 population3
ga_opt_ann_sel_clay_top = gaoptimset(ga_opt_ann_sel_clay_top, 'InitialPop',int8(population3));
[x fval exitflag output population4 scores] = ga(@ga_sel_pcr, 100, [],[],[],[],[],[],[], ga_opt_ann_sel_clay_top);
elapsed_time = toc;
save ga_top#7_bag100rand_RUN#4.mat
disp([ 'Time:     ' num2str(elapsed_time) ]);
disp([ 'Fitness:  ' num2str(fval) ]);
clear


%% ---- profile
%% -------- prepare data || protocol ==> {p2, 'ttt', ES}
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
%-PREPARE tr AND va SUBSETS
curr_support = 2;   % top:1 -- pro:2 -- hor:3
curr_data = 2;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva] = sannd(curr_support, curr_data);
clear curr_support curr_data

%-NORMALIZE
mm_p = [max([ptr,pva]')' min([ptr,pva]')'];
mm_t = [max([ttr,tva]')' min([ttr,tva]')'];
ptr = normalize(ptr, mm_p,-1,+1);
ttr = normalize(ttr, mm_t,-1,+1);
pva = normalize(pva, mm_p,-1,+1);
tva = normalize(tva, mm_t,-1,+1);

% ============   LOOP   ============
nSim = 200; %number of simulations
perc_te = 0.125; %percentage to assign to Te subset
TF = 'lm'; %training function to be used
%-initialization of output simulation arrays
otr = zeros(round(size(ptr,2)*(1-perc_te)),nSim);
ote = zeros(round(size(ptr,2)*perc_te),nSim);
ova = zeros(size(pva,2),nSim);
%-loop
for i = 1:nSim
    progress_bar(i, nSim);
    %-Te
    F_te(1,:) = randsample(size(ptr,2), round(size(ptr,2)*perc_te));
    F_tr = setxor(1:size(ptr,2),F_te);
    pte = ptr(:,F_te);
    tte = ttr(:,F_te);
    ptrs = ptr(:,F_tr); %small ptr
    ttrs = ttr(:,F_tr); %small ttr
    %-net
    %     newff_stack(in_range, topo, transferF, initial_weights);
    net = newff_stack([-1 +1], [size(ptrs,1) size(ptrs,1) 1], 'ttt', 1);
    net.(TF).trainParam.epochs = 60;
    net.(TF).trainParam.show = NaN;
    %-training
    TV.P = pva; TV.T = tva; %for estimating net ability to generalize
    VV.P = pte; VV.T = tte; %for early stopping criteria
    [n{i},trainRec{i}] = train(net.(TF), ptrs,ttrs,[],[],VV,TV);
    %-simulation
    Otr = sim(n{i},ptrs);
    Ote = sim(n{i},pte);
    Ova = sim(n{i},pva);
    Otr = denormalize(Otr, mm_t, -1, 1);
    Ote = denormalize(Ote, mm_t, -1, 1);
    Ova = denormalize(Ova, mm_t, -1, 1);
    otr(:,i) = Otr;
    ote(:,i) = Ote;
    ova(:,i) = Ova;
end
ttr = denormalize(ttrs, mm_t,-1,+1);
tte = denormalize(tte, mm_t,-1,+1);
tva = denormalize(tva, mm_t,-1,+1);
ttr = ttr'; tte = tte'; tva = tva';
% ============   LOOP   ============

%-SAVE
save GA\ga_prof#6_bagnet200.mat otr ote ova ttr tte tva n trainRec

%% ------------ REDUCE at RANDOM

%------------------------------
F_size = 100; % {25, 50, 100}
%------------------------------

cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
load GA\ga_prof#6_bagnet200
F = randsample(size(otr,2), F_size);
otr = otr(:,F);
ote = ote(:,F);
ova = ova(:,F);
n = n(:,F);
trainRec = trainRec(:,F);
clear F
eval(['save GA\ga_prof#6_bagnet' num2str(F_size) '_rand.mat otr ote ova ttr tte tva n trainRec'])

%% -------- run genetic algorithm computations
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   I use here the ga option structure array created for topsoil 'ga_opt_ann_sel_clay_top'.

%% ------------ random 100
clear
clc
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab\GA')
tic

% 1.\
disp('FIRST run')
load GA\ga_options ga_opt_ann_sel_clay_top
[x fval exitflag output population1 scores] = ga(@ga_sel_pcr, 100, [],[],[],[],[],[],[], ga_opt_ann_sel_clay_top);
elapsed_time = toc;
save ga_prof#6_bag100rand_RUN#1.mat
disp([ 'Time:     ' num2str(elapsed_time) ]);
disp([ 'Fitness:  ' num2str(fval) ]);
clear

% 2.\
disp('SECOND run')
load GA\ga_options ga_opt_ann_sel_clay_top
load ga_prof#6_bag100rand_RUN#1 population1
ga_opt_ann_sel_clay_top = gaoptimset(ga_opt_ann_sel_clay_top, 'InitialPop',int8(population1));
[x fval exitflag output population2 scores] = ga(@ga_sel_pcr, 100, [],[],[],[],[],[],[], ga_opt_ann_sel_clay_top);
elapsed_time = toc;
save ga_prof#6_bag100rand_RUN#2.mat
disp([ 'Time:     ' num2str(elapsed_time) ]);
disp([ 'Fitness:  ' num2str(fval) ]);
clear

% 3.\
disp('THIRD run')
load GA\ga_options ga_opt_ann_sel_clay_top
load ga_prof#6_bag100rand_RUN#2 population2
ga_opt_ann_sel_clay_top = gaoptimset(ga_opt_ann_sel_clay_top, 'InitialPop',int8(population2));
[x fval exitflag output population3 scores] = ga(@ga_sel_pcr, 100, [],[],[],[],[],[],[], ga_opt_ann_sel_clay_top);
elapsed_time = toc;
save ga_prof#6_bag100rand_RUN#3.mat
disp([ 'Time:     ' num2str(elapsed_time) ]);
disp([ 'Fitness:  ' num2str(fval) ]);
clear

% 4.\
disp('FOURTH run')
load GA\ga_options ga_opt_ann_sel_clay_top
load ga_prof#6_bag100rand_RUN#3 population3
ga_opt_ann_sel_clay_top = gaoptimset(ga_opt_ann_sel_clay_top, 'InitialPop',int8(population3));
[x fval exitflag output population4 scores] = ga(@ga_sel_pcr, 100, [],[],[],[],[],[],[], ga_opt_ann_sel_clay_top);
elapsed_time = toc;
save ga_prof#6_bag100rand_RUN#4.mat
disp([ 'Time:     ' num2str(elapsed_time) ]);
disp([ 'Fitness:  ' num2str(fval) ]);
clear

%% ---- horizon
%% -------- prepare data - NOT GOOD, OLD
% From boxplot cell above I recognized protocol #9 {p3,tt} as the best one for horizon support.
% Here I use the outputs recorded within horizon_#9.mat file to start genetic algorithm computations.

%-LOAD DATA
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
load sann_output\horizon_#9 n

%-PREPARE tr AND va SUBSETS
curr_support = 3;   % top:1 -- pro:2 -- hor:3
curr_data = 3;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva] = sannd(curr_support, curr_data, 60, NaN, 1000); %(support,epochs,show,nSim)
clear curr_support curr_data

%-NORMALIZE
pte = [];
tte = [];
mm_p = [max([ptr,pva,pte]')' min([ptr,pva,pte]')'];
mm_t = [max([ttr,tva,tte]')' min([ttr,tva,tte]')'];
ptr = normalize(ptr, mm_p,-1,+1);
pva = normalize(pva, mm_p,-1,+1);
clear pte tte

%-SIMULATION PHASE
% select 200 random NN from n
rand_sel = randsample(1000,200);
for i = 1:size(rand_sel,1)
    progress_bar(i, size(rand_sel,1))
    otr(:,i) = sim(n{rand_sel(i)},ptr);
    ova(:,i) = sim(n{rand_sel(i)},pva);
%    Ind(i) = perfind(otr(:,i),ttr');
end

%-DENORMALIZE
otr = norm2real(-1, 1,mm_t(2),mm_t(1),otr);
ova = norm2real(-1, 1,mm_t(2),mm_t(1),ova);

%-SELECTION BASED ON rmse
% r = [Ind.rmse];
% prctl = prctile(r,[25 75]);
% r1=find(r>prctl(1));
% r2=find(r<prctl(2));
% r_int = intersect(r1,r2);

    
%-SAVE
ttr = ttr'; tva = tva';
save GA\ga_hor#9.mat otr ova ttr tva

%% -------- prepare data - GOOD
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
%-PREPARE tr AND va SUBSETS
curr_support = 3;   % top:1 -- pro:2 -- hor:3
curr_data = 3;      % p1:1  -- p2:2  -- p3:3
[ptr,ttr,pva,tva] = sannd(curr_support, curr_data, 60, NaN, 1000); %(support,epochs,show,nSim)
clear curr_support curr_data

%-NORMALIZE
mm_p = [max([ptr,pva]')' min([ptr,pva]')'];
mm_t = [max([ttr,tva]')' min([ttr,tva]')'];
ptr = normalize(ptr, mm_p,-1,+1);
ttr = normalize(ttr, mm_t,-1,+1);
pva = normalize(pva, mm_p,-1,+1);
tva = normalize(tva, mm_t,-1,+1);

% ============   LOOP   ============
nSim = 200; %number of simulations
perc_te = 0.125; %percentage to assign to Te subset
TF = 'lm'; %training function to be used
%-initialization of output simulation arrays
otr = zeros(round(size(ptr,2)*(1-perc_te)),nSim);
ote = zeros(round(size(ptr,2)*perc_te),nSim);
ova = zeros(size(pva,2),nSim);
%-loop
for i = 1:nSim
    progress_bar(i, nSim);
    %-Te
    F_te(1,:) = randsample(size(ptr,2), round(size(ptr,2)*perc_te));
    F_tr = setxor(1:size(ptr,2),F_te);
    pte = ptr(:,F_te);
    tte = ttr(:,F_te);
    ptrs = ptr(:,F_tr); %small ptr
    ttrs = ttr(:,F_tr); %small ttr
    %-net
    %     newff_stack(in_range, topo, transferF, initial_weights);
    net = newff_stack([-1 +1], [size(ptrs,1) 1], 'tt', 1);
    net.(TF).trainParam.epochs = 60;
    net.(TF).trainParam.show = NaN;
    %-training
    [n{i},trainRec{i}] = train(net.lm, ptrs,ttrs);
    %-simulation
    Otr = sim(n{i},ptrs);
    Ote = sim(n{i},pte);
    Ova = sim(n{i},pva);
    Otr = denormalize(Otr, mm_t, -1, 1);
    Ote = denormalize(Ote, mm_t, -1, 1);
    Ova = denormalize(Ova, mm_t, -1, 1);
    otr(:,i) = Otr;
    ote(:,i) = Ote;
    ova(:,i) = Ova;
end
ttr = denormalize(ttrs, mm_t,-1,+1);
tte = denormalize(tte, mm_t,-1,+1);
tva = denormalize(tva, mm_t,-1,+1);
ttr = ttr'; tte = tte'; tva = tva';
% ============   LOOP   ============

%-SAVE
save GA\ga_hor#9_bagnet200.mat otr ote ova ttr tte tva n trainRec

%% ------------ REDUCE with RMSE
%% ---------------- reduce size of BAGNET to 115 elements
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
load GA\ga_hor#9_bagnet200
for i = 1:200
    P(i) = perfind(otr(:,i),ttr);
end
rmse = [P.rmse];
F = find(rmse(rmse<16));
otr = otr(:,F);
ote = ote(:,F);
ova = ova(:,F);
n = n(:,F);
trainRec = trainRec(:,F);
clear i rmse F P
save GA\ga_hor#9_bagnet115.mat otr ote ova ttr tte tva n trainRec

%% ---------------- reduce size of BAGNET to 50 elements
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
load GA\ga_hor#9_bagnet200
for i = 1:200
    P(i) = perfind(otr(:,i),ttr);
end
rmse = [P.rmse];
F = find(rmse(rmse<15.4));
otr = otr(:,F);
ote = ote(:,F);
ova = ova(:,F);
n = n(:,F);
trainRec = trainRec(:,F);
clear i rmse F P
save GA\ga_hor#9_bagnet50.mat otr ote ova ttr tte tva n trainRec

%% ---------------- reduce size of BAGNET to 15 elements
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
load GA\ga_hor#9_bagnet200
for i = 1:200
    P(i) = perfind(otr(:,i),ttr);
end
rmse = [P.rmse];
F = find(rmse(rmse<15));
otr = otr(:,F);
ote = ote(:,F);
ova = ova(:,F);
n = n(:,F);
trainRec = trainRec(:,F);
clear i rmse F P
save GA\ga_hor#9_bagnet15.mat otr ote ova ttr tte tva n trainRec

%% ------------ REDUCE at RANDOM

%------------------------------
F_size = 25; % {25, 50, 100}
%------------------------------

cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
load GA\ga_hor#9_bagnet200
F = randsample(200, F_size);
otr = otr(:,F);
ote = ote(:,F);
ova = ova(:,F);
n = n(:,F);
trainRec = trainRec(:,F);
clear F
eval(['save GA\ga_hor#9_bagnet' num2str(F_size) '_rand.mat otr ote ova ttr tte tva n trainRec'])

%% -------- run genetic algorithm computations
%% ---------------- example
a(1:200) = 0;
b = randsample(200,25);
a(b)=1;
z = ga_sel_pcr(a);

%% ---------------- try code
%% -------------------- random 100
cd('J:\work\MyPapers\clay_model_comparison\MatLab\GA')

% 1.\
load GA\ga_options ga*2 gra*
ga_opt_ann_sel_clay2.StallTimeLimit = 900;
clear graphics
[x fval exitflag output population1 scores] = ga(@ga_sel_pcr, 100, [],[],[],[],[],[],[],ga_opt_ann_sel_clay2);
save ga_hor#9_bag100rand_RUN#1.mat

% 2.\
load GA\ga_options ga*2 gra*
ga_opt_ann_sel_clay2.PlotFcns = graphics([1 2 3 10]);
ga_opt_ann_sel_clay2.StallTimeLimit = 3600;
clear graphics
[x fval exitflag output population2 scores] = ga(@ga_sel_pcr, 100, [],[],[],[],[],[],[],ga_opt_ann_sel_clay2);
save ga_hor#9_bag100rand_RUN#2.mat

% 3.\
load ga_hor#9_bag100rand_RUN#2.mat ga_opt_ann_sel_clay2 population
% Now I capitalize the above run by starting a new run with initial
% population equal to the final population of previous run named 'population'
ga_opt_ann_sel_clay2 = gaoptimset(ga_opt_ann_sel_clay2, 'InitialPop',int8(population2));
ga_opt_ann_sel_clay2.SelectionFcn = @selectionstochunif;
clear population2
%ga_opt_ann_sel_clay2.InitialPopulation = int8(population);
[x fval exitflag output population3 scores] = ga(@ga_sel_pcr, 100, [],[],[],[],[],[],[], ga_opt_ann_sel_clay2);
save ga_hor#9_bag100rand_RUN#3.mat %!!!!!!!! IL FILE E' STATO PERSO !!!!!!!!

% 4.\
load ga_hor#9_bag100rand_RUN#3.mat ga_opt_ann_sel_clay2 population3
% Now I capitalize the above run by starting a new run with initial
% population equal to the final population of previous run named 'population'
ga_opt_ann_sel_clay2 = gaoptimset(ga_opt_ann_sel_clay2, 'InitialPop',int8(population3));
clear population3
%ga_opt_ann_sel_clay2.InitialPopulation = int8(population);
[x fval exitflag output population4 scores] = ga(@ga_sel_pcr, 100, [],[],[],[],[],[],[], ga_opt_ann_sel_clay2);
save ga_hor#9_bag100rand_RUN#4.mat


%=========================
% NOTES:
% ga_hor#9_bag100rand_RUN#1.mat was terminated by "stall time limit"; RMSE = 9.1;
% ga_hor#9_bag100rand_RUN#2.mat was terminated by "stall time limit"; RMSE = 8.9;
% ga_hor#9_bag100rand_RUN#3.mat was terminated by "stall time limit"; RMSE = 8.7;
% ga_hor#9_bag100rand_RUN#4.mat was terminated by "stall time limit"; RMSE = 8.7;
%=========================

%% -------------------- all 200
load GA\ga_options ga*2 gra*
%load GA\ga_hor#9_bagnet200
ga_opt_ann_sel_clay2.PlotFcns = graphics([1 2 3 10]);
ga_opt_ann_sel_clay2.StallTimeLimit = 3600;
[x fval exitflag output population scores] = ga(@ga_sel_pcr, 200, [],[],[],[],[],[],[],ga_opt_ann_sel_clay2);

%% ---------------- try tool
% I run gatool with default options and stall time limit = 200 sec; run was stopped at 22nd generation
% I imported garesults in workspace from gatool.
x = garesults.x;
load GA\ga_hor#9_bagnet100_rand

%-compare rmse of all nn with respect to nn selected by GA (selected nn are those with x=1) 
F = find(x==0);
for i = 1:size(otr,2)
    P(i) = perfind(otr(:,i),ttr);
end
rmse = [P.rmse];
figure;scatter(rmse,x)

%-build up overall response
otr = otr(:,x==1);
ote = ote(:,x==1);
ova = ova(:,x==1);

%  1.\ PCRcv
   [w,nan] = PCRcv(otr,ttr);
%  2.\ create BAGNET response on validation subset (the one used for early stopping criteria)
   Yb_te = bagnet(ote,w);
%  3.\ recognize the best PC number on validation subset
   [mse_te w_min_te] = MSEg(Yb_te,tte);
%  4.\ create BAGNET response on validation subset
   Yb_va = bagnet(ova,w(:,w_min_te)); 
%  5.\ performance evaluation on validation subset
   Ind = perfind(Yb_va, tva);



%% =================================================================
%% BOXPLOT {SMU, linear, geostatistics, ANN, bagging-ANN-GA-PCR}
% Here I compute the boxplot of all points

clear; clc;

%------------------------------------------------
curr_support = 3;   % top:1 -- pro:2 -- hor:3
curr_perf = 2;      % 'rmse', 'r', 'mbe2', 'smape', 'eff', 'D'
GA_perf = 5;        % {'mse','rmse','mae','mbe2','r','smape','eff','D'};
fontsize = 14;      % ylabel fontsize
minus = 6;          % text fontsize = ylabel fontsize - minus
fontweight = {'n','b'}; %fontweight of all FontWeights
fw_text = 2;            %fontweight of all text commands
linewidth = 0.5;
linestyle = {'-'};
%color of line and of corresponding text
linecolor = {'k','k','k','k'};  % {'k','b','r','m'};
rotation = 0; % {0, 90}
Hor_Align = 'Left';     % { 'Left', 'Center', 'Right' }
Ver_Align = 'Middle';   % { 'Middle', 'Top', 'Bottom', ... }
%------------------------------------------------
support = {'Topsoil', 'Profile', 'Horizon'; 'topsoil','profile','horizon'; 'top','prof','hor'; ...
           'REGR', 'REGR', 'GLM'; 'MOCOK', 'MOCOK', 'OK'; ; ...
           'regr','regr','glm_all'; 'mocok','mocok','ok'; '7','6','9'};
perf = {'rmse', 'r', 'mbe2', 'smape', 'eff', 'D'; ...
        'RMSE', 'Pearson''s r', 'MBE', 'SMAPE', 'Efficiency', 'Willmott''s D'};

cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab\sann_output')
bplot_str = '[';
for i = 1:10
    eval(['load ' support{2,curr_support} '_#' num2str(i) '.mat PERF' ])
    eval(['perf_' num2str(i) ' = PERF;'])
    clear PERF
    bplot_str = [bplot_str ' [perf_' num2str(i) '.' (perf{1,curr_perf}) ']'''];
end
bplot_str = [bplot_str ' ]'];
clear i

% lines: glm, ok
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
load results ind*
figure
hold on
    %lab = {{'p1';'tt'}, {'p1';'ttt'}, {'p1';'tt';'ES'}, {'p1';'ttt';'ES'}, {'p2';'tt';'ES'}, {'p2';'ttt';'ES'},  {'p3';'tt';'ES'}, {'p3';'ttt';'ES'} };
    %      1        2         3           4            5           6            7           8
    lab = {'p1-tt', 'p1-ttt', 'p1-tt-ES', 'p1-ttt-ES', 'p2-tt-ES', 'p2-ttt-ES', 'p3-tt-ES', 'p3-ttt-ES', 'p3-tt', 'p3-ttt'};
    eval(['boxplot(' bplot_str ', ''labels'',lab, ''notch'',''on'')'])
    ylabel(perf{2,curr_perf}, 'FontSize',fontsize, 'FontWeight',fontweight{2})
    xlabel('ANN', 'FontSize',fontsize-minus, 'FontWeight',fontweight{fw_text})
    %GLM
    eval(['s(1:size(lab,2)+2) = ind_' support{3,curr_support} '_' support{6,curr_support} '.(perf{1,curr_perf});'])
    line([0:size(s,2)-1],s,'LineStyle',linestyle{1}, 'LineWidth',linewidth, 'Color',linecolor{1})
    text(size(lab,2)+0.5,s(1), support{4,curr_support}, 'Color',linecolor{1}, 'FontSize',fontsize-minus, 'FontWeight',fontweight{fw_text}, 'Rotation',rotation, 'VerticalAlignment',Ver_Align, 'HorizontalAlignment',Hor_Align)
    %KRIGING
    eval(['s(1:size(lab,2)+2) = ind_' support{3,curr_support} '_' support{7,curr_support} '.(perf{1,curr_perf});'])
    line([0:size(s,2)-1],s,'LineStyle',linestyle{1}, 'LineWidth',linewidth, 'Color',linecolor{2})
    text(size(lab,2)+0.5,s(1), support{5,curr_support}, 'Color',linecolor{2}, 'FontSize',fontsize-minus, 'FontWeight',fontweight{fw_text}, 'Rotation',rotation, 'VerticalAlignment',Ver_Align, 'HorizontalAlignment',Hor_Align)
    %SMU
    eval(['s(1:size(lab,2)+2) = ind_' support{3,curr_support} '_udp.(perf{1,curr_perf});'])
    line([0:size(s,2)-1],s,'LineStyle',linestyle{1}, 'LineWidth',linewidth, 'Color',linecolor{3})
    text(size(lab,2)+0.5,s(1), 'SMU', 'Color',linecolor{3}, 'FontSize',fontsize-minus, 'FontWeight',fontweight{fw_text}, 'Rotation',rotation, 'VerticalAlignment',Ver_Align, 'HorizontalAlignment',Hor_Align)
    %bagging-ANN-GA-PCR
    eval( ['load GA\ga_' support{3,curr_support} '#' support{8,curr_support} '_bag100rand_RUN#4 x'])
    bagging = ga_sel_pcr(x,GA_perf);
    eval(['s(1:size(lab,2)+2) = bagging;'])
    line([0:size(s,2)-1],s,'LineStyle',linestyle{1}, 'LineWidth',linewidth, 'Color',linecolor{4})
    text(size(lab,2)+0.5,s(1), 'BAGGING', 'Color',linecolor{4}, 'FontSize',fontsize-minus, 'FontWeight',fontweight{fw_text}, 'Rotation',rotation, 'VerticalAlignment',Ver_Align, 'HorizontalAlignment',Hor_Align)
%     y_curly = ylim;
%     text(size(lab,2)/2,y_curly(1)-abs(y_curly(1))*0.1, '{')
hold off

%% -------- SAVE boxplot
%------------------------
save_dir = 'figures';
%------------------------
cd('C:\Dottorato\MyPapers\clay_model_comparison\MatLab')
if exist(save_dir,'dir') == 0
    mkdir(save_dir)
end
eval(['print -dtiff -r1200 ' save_dir filesep support{3,curr_support} '_' perf{1,curr_perf}])
