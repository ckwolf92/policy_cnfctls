%% MONETARY POLICY SHOCK IRFs: STANDARD ID SCHEMES
% Alisdair McKay & Christian Wolf
% this version: 03/24/2023

%% HOUSEKEEPING
 
clearvars -except path task
cd([path task]);

%% DATA

% import data

data_table = readtable('_data.csv');
data = table2array(data_table);
date = data(:,1);

series = data_table.Properties.VariableNames;
series = series(2:end);

% macro outcomes

ygap   = data(:,2);
pgdp   = data(:,3);
lpgdp  = log(pgdp);
dlpgdp = [0;diff(lpgdp)];
infl   = 4*dlpgdp*100;
ffr    = data(:,4);
lpcom  = data(:,5);

% monetary shocks

rr_3       = data(:,6);
mp1_tc     = data(:,7);
ad_shock   = data(:,8);
mar_svariv = data(:,9);

% technology shock

bzk_ist_news   = data(:,10);

% dates

startdate = find(date == 1969);
enddate = find(date == 2007.75);

% collect VAR inputs

vardata = [mp1_tc ygap infl lpcom rr_3 ffr];
vardata = vardata(startdate:enddate,:);
vardata(isnan(vardata)) = 0;

% series names

series_names = {'GK Shock','Output','Inflation','Commodity Prices','RR Shock','Interest Rate'};

%% SETTINGS

% VAR specification

n_lags     = 4;                    % number of lags
constant   = 2;                    % constant?
IRF_hor    = 100;
n_draws    = 10000;
n_y        = size(vardata,2);

% identification

IS.shock_pos = [1 5];
n_shocks     = length(IS.shock_pos);

% outcomes of interest

IS.y_pos = [2 3 6];

% placeholders

IS.IRF     = NaN(IRF_hor,n_y,n_shocks,n_draws);
IS.IRF_med = NaN(IRF_hor,n_y,n_shocks);
IS.IRF_lb  = NaN(IRF_hor,n_y,n_shocks);
IS.IRF_ub  = NaN(IRF_hor,n_y,n_shocks);
IS.IRF_OLS = NaN(IRF_hor,n_y,n_shocks);

%% VAR ESTIMATION

%----------------------------------------------------------------
% Estimate Reduced-Form VAR
%----------------------------------------------------------------

T = size(vardata,1) - n_lags;
[B_draws,Sigma_draws,B_OLS,Sigma_OLS] = bvar_fn(vardata,n_lags,constant,n_draws);

%----------------------------------------------------------------
% OLS IRFs
%----------------------------------------------------------------

% extract VAR inputs
    
Sigma_u   = Sigma_OLS;
B         = B_OLS;

% benchmark rotation

bench_rot = chol(Sigma_u,'lower');

% Wold IRFs

IRF_Wold = zeros(n_y,n_y,IRF_hor); % row is variable, column is shock
IRF_Wold(:,:,1) = eye(n_y);

for l = 1:IRF_hor
    
    if l < IRF_hor
        for j=1:min(l,n_lags)
            IRF_Wold(:,:,l+1) = IRF_Wold(:,:,l+1) + B(1+(j-1)*n_y:j*n_y,:)'*IRF_Wold(:,:,l-j+1);
        end
    end
    
end

W = bench_rot;

% get IRFs

IRF_OLS = NaN(n_y,n_y,IRF_hor);
for i_hor = 1:IRF_hor
    IRF_OLS(:,:,i_hor) = IRF_Wold(:,:,i_hor) * W;
end

% collect results

for i_shock = 1:length(IS.shock_pos)
    IS.IRF_OLS(:,:,i_shock) = squeeze(IRF_OLS(:,IS.shock_pos(i_shock),:))';
end

%----------------------------------------------------------------
% Draw IRFs
%----------------------------------------------------------------

for i_draw = 1:n_draws
    
% extract VAR inputs
    
Sigma_u   = Sigma_draws(:,:,i_draw);
B         = B_draws(:,:,i_draw);

% benchmark rotation

bench_rot = chol(Sigma_u,'lower');

% Wold IRFs

IRF_Wold = zeros(n_y,n_y,IRF_hor); % row is variable, column is shock
IRF_Wold(:,:,1) = eye(n_y);

for l = 1:IRF_hor
    
    if l < IRF_hor
        for j=1:min(l,n_lags)
            IRF_Wold(:,:,l+1) = IRF_Wold(:,:,l+1) + B(1+(j-1)*n_y:j*n_y,:)'*IRF_Wold(:,:,l-j+1);
        end
    end
    
end

W = bench_rot;

% get IRFs

IRF_draw = NaN(n_y,n_y,IRF_hor);
for i_hor = 1:IRF_hor
    IRF_draw(:,:,i_hor) = IRF_Wold(:,:,i_hor) * W;
end

% collect results

for i_shock = 1:length(IS.shock_pos)
    IS.IRF(:,:,i_shock,i_draw) = squeeze(IRF_draw(:,IS.shock_pos(i_shock),:))';
end

end

for ii=1:IRF_hor
    for jj=1:n_y
        for kk = 1:n_shocks
            IS.IRF_med(ii,jj,kk) = quantile(IS.IRF(ii,jj,kk,:),0.5);
            IS.IRF_lb(ii,jj,kk) = quantile(IS.IRF(ii,jj,kk,:),0.16);
            IS.IRF_ub(ii,jj,kk) = quantile(IS.IRF(ii,jj,kk,:),0.84);
        end
    end
end

%% SAVE RESULTS

cd([path task '/_results'])

IS_1.IRF     = squeeze(IS.IRF(:,IS.y_pos,1,:));
IS_1.IRF_lb  = squeeze(IS.IRF_lb(:,IS.y_pos,1));
IS_1.IRF_ub  = squeeze(IS.IRF_ub(:,IS.y_pos,1));
IS_1.IRF_med = squeeze(IS.IRF_med(:,IS.y_pos,1));

IS_2.IRF     = squeeze(IS.IRF(:,IS.y_pos,2,:));
IS_2.IRF_lb  = squeeze(IS.IRF_lb(:,IS.y_pos,2));
IS_2.IRF_ub  = squeeze(IS.IRF_ub(:,IS.y_pos,2));
IS_2.IRF_med = squeeze(IS.IRF_med(:,IS.y_pos,2));

save IS_mp_base IS_1 IS_2

clearvars -except path task
cd([path task]);