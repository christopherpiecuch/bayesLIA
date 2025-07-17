%   function bayes_main_code(namStr,itrNum)
function bayes_main_code(namStr,itrNum)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Say hello
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause(0.1)
disp('Hello.  Things have started.')
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preliminary input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% name of experiment and file to be saved
%nameOfExperiment='Bergrates_55st_newRF_FINALVERSION'; % experiment name
nameOfExperiment=namStr; % experiment name
save_name=[date,'_',nameOfExperiment,'_iteration_',num2str(itrNum)]; % file name

% number of draws to perform
NN_burn=100000;             % 100,000 warm-up draws
NN_post=100000;             % 100,000 post-warm-up draws
thin_period=100;            % thin chains keeping 1 of 100
%NN_burn=100;             % 100,000 warm-up draws
%NN_post=100;             % 100,000 post-warm-up draws
%thin_period=1;            % thin chains keeping 1 of 100

% Define iteration parameters based on input
NN_burn_thin=NN_burn/thin_period; % Total number of burn-in to keep
NN_post_thin=NN_post/thin_period; % Total number of post-burn-in to keep
NN=NN_burn+NN_post;               % Total number of draws to take 
NN_thin=NN_burn_thin+NN_post_thin;% Total number of draws to keep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load GPS data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load GPS data    
[GPS_DATA,GPS_ERROR,GPS_LON,GPS_LAT]=prepare_gps_data;

% Define space and time parameters related to data
[N_gps]=numel(GPS_DATA');
N=N_gps;
LON=[GPS_LON'];
LAT=[GPS_LAT'];

% define useful spatial values
D=EarthDistances([LON' LAT']); % distances between locations
L=sum(~isnan(GPS_DATA)); % number of crustal motion data sites

% Set the seeds of the random number generator
rng(itrNum*sum(clock))

% Define hyperparameter values
HP=set_hyperparameters(GPS_DATA);
[~,K]=size(HP.G);

% Allocate space for the sample arrays
alpha=[];
initialize_output

% Set initial values
set_initial_values

% Set up selection matrices
%H=zeros(N,N);
%H(1:N,1:N)=eye(N);
E=zeros(L,N);
E(1:L,1:L)=eye(L);
x=GPS_DATA;
DeltaMat=diag(GPS_ERROR.^2);
invDeltaMat=inv(DeltaMat);

% Set up identity matrices and vectors of zeros or ones
I_N=eye(N);
ONE_N=ones(N,1);
ZERO_N=zeros(N,1);

% Loop through the Gibbs sampler with Metropolis steps
%tic
for nn=1:NN
    if mod(nn,100)==0 
        %toc
        disp([num2str(nn),' of ',num2str(NN),' iterations done.'])
        %tic
    end
    nn_thin=[]; nn_thin=ceil(nn/thin_period);
     
    % Define matrices to save time
    OmegaMat=omega_2*exp(-rho*D); invOmegaMat=inv(OmegaMat);
    TMat=exp(-rho*D); invTMat=inv(TMat);
  
    % Sample from p(epsilon_2|.)
    inside1=[]; inside2=[];
    inside1=N/2;
    inside2=1/2*(v-u)'*(v-u);
    epsilon_2=1/randraw('gamma', [0,1/(HP.chi_tilde_epsilon_2+inside2),...
      	(HP.xi_tilde_epsilon_2+inside1)], [1,1]);
   	clear inside*
    
    % Sample from p(omega_2|.)
    inside1=[]; inside2=[];
    inside1=N/2;
    inside2=1/2*(u-HP.G*c-alpha*ONE_N)'*invTMat*(u-HP.G*c-alpha*ONE_N);
    omega_2=1/randraw('gamma', [0,1/(HP.chi_tilde_omega_2+inside2),...
      	(HP.xi_tilde_omega_2+inside1)], [1,1]);
   	clear inside*
    
    % Sample from p(rho|.)
    Rho_now=log(rho);
    Rho_std=0.05;
    Rho_prp=normrnd(Rho_now,Rho_std);
    L_now=exp(-exp(Rho_now)*D);
    L_prp=exp(-exp(Rho_prp)*D);
    invL_now=inv(L_now);
    invL_prp=inv(L_prp);
    sumk_now=0;
    sumk_prp=0;
    sumk_now=(u-alpha*ONE_N-HP.G*c)'*invL_now*(u-alpha*ONE_N-HP.G*c);
    sumk_prp=(u-alpha*ONE_N-HP.G*c)'*invL_prp*(u-alpha*ONE_N-HP.G*c);
 	ins_now=-1/(2*HP.zeta_tilde_rho_2)*(Rho_now-HP.eta_tilde_rho)^2-1/(2*omega_2)*sumk_now;
   	ins_prp=-1/(2*HP.zeta_tilde_rho_2)*(Rho_prp-HP.eta_tilde_rho)^2-1/(2*omega_2)*sumk_prp;
  	MetFrac=det(L_prp*invL_now)^(-1/2)*exp(ins_prp-ins_now);
   	success_rate=min(1,MetFrac);
   	if rand(1)<=success_rate
     	Rho_now=Rho_prp; 
    end
  	rho=exp(Rho_now);
  	clear Rho_now Rho_std Rho_prp mat_now mat_prp ins_* sumk MetFrac success_rate L_*
    % redefine matrices because you just updated omega_2 and rho
    OmegaMat=omega_2*exp(-rho*D); invOmegaMat=inv(OmegaMat);
    TMat=exp(-rho*D); invTMat=inv(TMat);

    % Sample from p(alpha|.)
    V_ALPHA=[]; PSI_ALPHA=[];
    V_ALPHA=HP.eta_tilde_alpha/HP.zeta_tilde_alpha_2+ONE_N'*invOmegaMat*(u-HP.G*c);
   	PSI_ALPHA=inv(1/HP.zeta_tilde_alpha_2+ONE_N'*invOmegaMat*ONE_N);
    alpha=normrnd(PSI_ALPHA*V_ALPHA,sqrt(PSI_ALPHA));
    clear V_ALPHA PSI_ALPHA  
    
    % Sample from p(v|.)
   	V_V=[]; PSI_V=[]; 
    V_V=u/epsilon_2+E'*invDeltaMat*x;
    PSI_V=inv(I_N/epsilon_2+E'*invDeltaMat*E);
    v=mvnrnd(PSI_V*V_V,PSI_V)';
    clear V_V PSI_V
    
    % Sample from p(u|.)
   	V_U=[]; PSI_U=[]; 
    V_U=invOmegaMat*(HP.G*c+alpha*ONE_N)+v/epsilon_2;
    PSI_U=inv(invOmegaMat+I_N/epsilon_2);
    u=mvnrnd(PSI_U*V_U,PSI_U)';
    clear V_U PSI_U
    
    % Sample from p(c|.)
    c_now=c;
    c_prp=mnrnd(1,pi)';
 	ins_now=-0.5*(u-alpha*ONE_N-HP.G*c_now)'*invOmegaMat*(u-alpha*ONE_N-HP.G*c_now);
 	ins_prp=-0.5*(u-alpha*ONE_N-HP.G*c_prp)'*invOmegaMat*(u-alpha*ONE_N-HP.G*c_prp);
  	MetFrac=prod(pi.^(c_prp-c_now))*exp(ins_prp-ins_now);
   	success_rate=min(1,MetFrac);
   	if rand(1)<=success_rate
     	c_now=c_prp; 
    end
  	c=c_now;
  	clear c_now c_prp ins_* MetFrac success_rate

    % Sample from p(pi|.)
    pi=drchrnd(HP.mu*ones(1,K)+c',1)';

    % Now update arrays
    update_all_arrays

end
%toc

% delete the burn-in period values
delete_burn_in

% save output
if exist('bayes_model_solutions/')==0
    mkdir bayes_model_solutions/
end
save(['bayes_model_solutions/experiment_',save_name,'.mat'],...
    'U','C','PI','V','ALPHA','RHO','OMEGA_2','EPSILON_2','HP','*DATA','*LON','*LAT',...
    '*ERROR','N','L','D','nn')