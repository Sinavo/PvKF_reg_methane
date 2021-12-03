%% 
%%% ======================================================================
%%% main_loop_hemi.m
%%% Created by Sina Voshtani 
%%% Created on 21/09/2021
%%% =======================================================================

%%% comment "load" below if bias-corrected GOSAT data is not used
load('gosat_cmaq_biascorrected_qc_on_latitude.mat','lat_sat_all_bc'...
     ,'lon_sat_all_bc','tim_sat_all_bc','ch4_sat_all_bc',...
     'sig_sat_all_bc','plv_sat_all_bc','pwe_sat_all_bc',...
     'fgs_sat_all_bc','avk_sat_all_bc','xb_sat_all_bc')
%%% 
%%% PvKF for loop over days
for i = start_day:end_day
    tic
    
%%% extract daily values from netcdf files that are read 
    prep_day % check if this file exists in the current directory
    
%%% read satellite observation data and process for assimilation  
    %read_gosat % uncomment when not using bias corrected GOSAT data (.mat)     
    %read_sciamachy
    %read_iasi
    %read_tes

%%% to re-start from the middle
    if i == 1
       start_hour = 1;
       %load('BiasAMF_DA_i4f0.8q0_gos.mat') % uncomment to restart. Check
       %if .mat exists
       else
       start_hour = 1;
    end
    
%%% PvKF for loop over hours    
    for nstep=start_hour:1:step_size
        
%%% model forecast runs of concentration and error varaince        
%	run_hcmaq % check if this file exists in the current directory

%%% extract hourly values from netcdf files that are read 
	prep_hour % check if this file exists in the current directory
    
%%% process satellite observation data: make it ready for assimilation    
    calc_gosat_bias_corrected
    % calc_gosat_offline
    % calc_sciamachy
    
%%% Analysis step of PvKF        
    analysis_plus
    
%%% replace analysis in model concentrations for the next forecast    
    write_cgrid
    
%%% save important variables for the post-processing     
save('optDA_i2.5f4q10.mat','lat_sat_all','lon_sat_all','ch4_sat_all','sig_sat_all',...
     'plv_sat_all','pwe_sat_all','fgs_sat_all','avk_sat_all','xb_sat_all_bc',...
     'chi_sq_all','obs_num_all','innov_all','xf_o_avg_all','var_B_all')
    end
end

%%% =======================================================================
%%% END
%%% =======================================================================
