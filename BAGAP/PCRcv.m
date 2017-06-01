function [w,Yb] = PCRcv(O,T)
% [w,Yb] = PCRcv(O,T)
%
% DESCRIPTION
%   Questa funzione serve ad eseguire n=(pcE-pcS)+1 volte una PCR con pc 
%   componenti principali su b output forniti da b esperti bootstrap. 
%   Sono creati due matrici rilevanti:
%   ..:: w(pcE,B) e Yb(staz*dec,pcE) dove pcE = B ::..
%   -cv
%   "cv" sta per cross validation che è lo scopo di questa FUNCTION: creare
%   w e Yb da adoperare in un'analisi di sensibilità del tipo cross
%   validation su testing set (stazioni e decadi non partecipanti
%   all'addestramento) per capire quante componenti principali adoperare
%   sul simulation set (tutti i pixel interessati dall'interpolazione nel
%   dominio del tempo e dello spazio).
%
% INPUTS
%   O(staz*dec,B)            : ANN simulation output. This array is
%                              obtained simulating each iANN with
%   T(staz*dec,1)            : target, measured values
%
% OUTPUTS
%   w(B,B)                   : matrice dei pesi della PCR
%   Yb(staz*dec,pcE-pcS+1=B) : matrice di aggregazione su pc principal
%                              components
%
% WHERE
%   B   : number of bootstrap resampling (or "experts")
%   pcS : Start number of principal components to employ
%   pcE : End number of principal components to employ

%% account for possible NaNs
if sum(isnan(T(:))) || sum(isnan(O(:)))
    error('Detected NaNs in O/T arrays. NaNs are not allowed.');
end
% CodeBar     = find(not(isnan(T)));
% With_NaNs   = find(isnan(T));
% T2          = T(CodeBar,:);
% O2          = O(CodeBar,:);
%T2 = T;
%O2 = O;

%% run PCA

[z,t] = size(O);

% B=t, pcS=1 sempre, pcE=t. Diversamente vedi PCRcv_primordiale.

% pc=loadings | zscores=scores | pcvars=eigenvalues of the COV(O)
%[pc, zscores, pcvars] = pca(O);
[pc, ~, pcvars] = pca(O);
cumpercent = cumsum(pcvars/sum(pcvars));

% if pcvars is below tolerance, then exclude all components after the first
% occurrence:
if ~isempty( find(pcvars<1e-3,1) )
    t2 = find(pcvars<1e-3,1) -1;
else
    t2 = t;
end

%% PLOT, eventually

fig = 0;
if fig == 1
    figure
    subplot(1,2,1), plot(1:1:t,cumpercent,'xr'); xlabel('Number of Components'); ylabel('% Cumulated Variance Explained'); legend('% pcvars')
    subplot(1,2,2), plot(pcvars); xlabel('Number of Components'); ylabel('Eigenvalues of Covariance Matrix'); legend('pcvars')
end

%% cross validation on pc

w = zeros(t,t2);
Yb = zeros(z,t2);
for k = 1:t2
    % filtro della matrice pc rispetto alla scelta del numero di componenti
    Pk = pc(:,1:k);
    % calcolo della matrice dei pesi w necessaria per aggregare i bootstrap
    %   [Zhang Jie, Neurocomputing 25 (1999) 93-113]
    w(:,k) = Pk*((Pk'*(O'*O)*Pk)^ -1)*Pk'*O'*T;
    % calcolo della pioggia aggregando tutti i bootstrap
    Yb(:,k) = O*w(:,k);
end
%Yb(CodeBar,:) = Yb2;
%Yb(With_NaNs,:) = NaN;

%% end
return