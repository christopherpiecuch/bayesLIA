%function HP = set_hyperparameters(GPS_DATA)
function HP = set_hyperparameters(GPS_DATA)

% Set variance inflation parameters
var_infl=5^2;
var_infl2=10^2;
var0_infl=1;

% Define hyperparameters
% alpha
HP.eta_tilde_alpha = nanmean(GPS_DATA); % Mean of alpha prior
HP.zeta_tilde_alpha_2 = var_infl*nanvar(GPS_DATA); % Variance of alpha prior

% omega_2
HP.xi_tilde_omega_2 = 1/2; % Shape of tau_2 prior
HP.chi_tilde_omega_2 = 1/2*nanvar(GPS_DATA); % 

% epsilon_2
HP.xi_tilde_epsilon_2 = 1/2; % Shape of tau_2 prior
HP.chi_tilde_epsilon_2 = 1/2*(1e-3)^2; % Guess (1 mm/yr)^2 error variance

% rho (this one's constrained; 95% within 500,2000 km)
HP.eta_tilde_rho = -6.9; % "Mean" of phi prior
HP.zeta_tilde_rho_2 = (0.35)^2; % "Variance" of phi prior

% GIA
load('VLM_st_Berg.mat')
VLM_st([21,22],:) = []; % remove two locations
load('RF_correction_BERG.mat') % load reference framew correction
HP.G=1e-3*(VLM_st+Ax);
[~,K]=size(HP.G);
HP.mu=1/K;
return
