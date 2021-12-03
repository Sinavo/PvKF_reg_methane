%%
%%% ======================================================================
%%% calc_gosat_bias_corrected.m
%%% Created by Sina Voshtani
%%% Created on 21/09/2021
%%% =======================================================================

%%% read satellite observations from bias corrected file (.mat)
lat_sat  =  nonzeros(lat_sat_all_bc(:,step_size*(i-1)+nstep));
lon_sat  =  nonzeros(lon_sat_all_bc(:,step_size*(i-1)+nstep));
ch4_sat  =  ch4_sat_all_bc(1:length(lat_sat),step_size*(i-1)+nstep);
sig_sat1 =  sig_sat_all_bc(1:length(lat_sat),step_size*(i-1)+nstep);
plv_sat  =  plv_sat_all_bc(:,1:length(lat_sat),step_size*(i-1)+nstep);
pwe_sat  =  pwe_sat_all_bc(:,1:length(lat_sat),step_size*(i-1)+nstep);
fgs_sat  =  fgs_sat_all_bc(:,1:length(lat_sat),step_size*(i-1)+nstep);
avk_sat  =  avk_sat_all_bc(:,1:length(lat_sat),step_size*(i-1)+nstep);

f_rep    = 0.4; % observation error factor (Part I Eq. 40)
sig_sat  = f_rep.*(sig_sat1); % multiplicitive obs. error

%%% bias corrections for different variables
%  ch4_sat  = ch4_sat_bi - (-0.46.*lat_sat - 13); % lattitude bias correction
%  ch4_sat  = ch4_sat_bi - (-0.55.*sol_sat - 12); % solar zenith angle bias correction
%  A_sat    = (1./cos(pi*sol_sat/180)) +  (1./cos(sen_sat*pi/180)); % AMF bais correction
%  ch4_sat  = ch4_sat_bi - (-21.*A_sat + 19); %AMF bais correction

%%% X-Y grided H-CMAQ on projection plane (ll2psn.m should be in current directory)
lat_sat_rad = deg2rad(lat_sat);
lon_sat_rad = deg2rad(lon_sat);
[x_sat,y_sat,z_sat] = sph2cart(lon_sat_rad,lat_sat_rad,r);
r_sat  = [x_sat(:),y_sat(:),z_sat(:)];
[xproj_sat,yproj_sat] = ll2psn(lat_sat,lon_sat,'TrueLat',45,'EarthRadius',6370000,...
    'Eccentricity',1e-30,'meridian',-98);

%%% =======================================================================
%%% END
%%% =======================================================================