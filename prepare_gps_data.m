% function [V_GPS,S_GPS,X_GPS,Y_GPS]=prepare_gps_data
function [V_GPS,S_GPS,X_GPS,Y_GPS]=prepare_gps_data

% load data and errors
load('deglacialRrate.mat')
load('Berg_Uplift.txt') 

% remove two sites
Berg_Uplift([21 22],:)=[];
Rrate_interp(21:22)=[];
Rrate_std_interp(21:22)=[];

% rename
GIA=Rrate_interp'.*1000; % convert to mm/yr
GIAsigma=Rrate_std_interp'.*1000; % convert to mm/yr
long=Berg_Uplift(:,2);
lat=Berg_Uplift(:,1);
OBS=Berg_Uplift(:,5);
OBSsigma=Berg_Uplift(:,9);
ela_gris=Berg_Uplift(:,6);
ela_grpg=Berg_Uplift(:,7);
ela_canpg=Berg_Uplift(:,8);
ela_gris_sigma=Berg_Uplift(:,10);
ela_grpg_sigma=Berg_Uplift(:,11);
ela_canpg_sigma=Berg_Uplift(:,12);
ELA=ela_gris+ela_grpg+ela_canpg;

% compute output
X_GPS=long;
Y_GPS=lat;
V_GPS=1e-3*(OBS-ELA-GIA); % convert to m/yr
S_GPS=1e-3*sqrt(OBSsigma.^2+ela_gris_sigma.^2+ela_grpg_sigma.^2+ela_canpg_sigma.^2+GIAsigma.^2); % convert to m/yr
return