% set_initial_values

% mean parameters
alpha=[];
alpha=normrnd(HP.eta_tilde_alpha,sqrt(HP.zeta_tilde_alpha_2));

% variance parameters
omega_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_omega_2,HP.xi_tilde_omega_2], [1,1])]); % use min to prevent needlessly large values
epsilon_2=min([1 1/randraw('gamma', [0,1/HP.chi_tilde_epsilon_2,HP.xi_tilde_epsilon_2], [1,1])]); % use min to prevent needlessly large values

% inverse length scale parameters
rho=exp(normrnd(HP.eta_tilde_rho,sqrt(HP.zeta_tilde_rho_2)));

% spatial fields
u=(mvnrnd(alpha*ones(N,1),omega_2*exp(-rho*D)))';        
v=(mvnrnd(u,epsilon_2*eye(N)))';        
pi=[]; pi=drchrnd(HP.mu*ones(1,K),1)';
c=mnrnd(1,pi)';