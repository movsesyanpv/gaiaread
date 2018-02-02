    module CustomDataTypes
    
        type GaiaEntry
            integer(8)    :: solution_id = 0
            integer(8)    :: source_id = 0
            integer(8)    :: random_index = 0
            real(8)       :: ref_epoch = 0
            real(8)       :: ra = 420
            real(8)       :: ra_error = 0
            real(8)       :: dec = 420
            real(8)       :: dec_error = 0
            real(8)       :: parallax = 0
            real(8)       :: parallax_error = 0
            real(8)       :: pmra = 1296000000
            real(8)       :: pmra_error = 0
            real(8)       :: pmdec = 1296000000
            real(8)       :: pmdec_error = 0
            real(4)       :: ra_dec_corr = 2
            real(4)       :: ra_parallax_corr = 2
            real(4)       :: ra_pmra_corr = 2
            real(4)       :: ra_pmdec_corr = 2
            real(4)       :: dec_parallax_corr = 2
            real(4)       :: dec_pmra_corr = 2
            real(4)       :: dec_pmdec_corr = 2
            real(4)       :: parallax_pmra_corr = 2
            real(4)       :: parallax_pmdec_corr = 2
            real(4)       :: pmra_pmdec_corr = 2
            integer(4)    :: astrometric_n_obs_al = -1
            integer(4)    :: astrometric_n_obs_ac = -1
            integer(4)    :: astrometric_n_good_obs_al = -1
            integer(4)    :: astrometric_n_good_obs_ac = -1
            integer(4)    :: astrometric_n_bad_obs_al = -1
            integer(4)    :: astrometric_n_bad_obs_ac = -1
            real(4)       :: astrometric_delta_q = -1
            real(8)       :: astrometric_excess_noise = 1296000000
            real(8)       :: astrometric_excess_noise_sig = -1
            character(16) :: astrometric_primary_flag = "false"
            real(4)       :: astrometric_relegation_factor = -1
            real(4)       :: astrometric_weight_al = -1
            real(4)       :: astrometric_weight_ac = -1
            integer(4)    :: astrometric_priors_used = 0
            integer(2)    :: matched_observations = -1
            character(16) :: duplicated_source = "false"
            real(4)       :: scan_direction_strength_k1 = 2
            real(4)       :: scan_direction_strength_k2 = 2
            real(4)       :: scan_direction_strength_k3 = 2
            real(4)       :: scan_direction_strength_k4 = 2
            real(4)       :: scan_direction_mean_k1 = 270
            real(4)       :: scan_direction_mean_k2 = 270
            real(4)       :: scan_direction_mean_k3 = 270
            real(4)       :: scan_direction_mean_k4 = 270
            integer(4)    :: phot_g_n_obs = -1
            real(8)       :: phot_g_mean_flux = -1
            real(8)       :: phot_g_mean_flux_error = -1
            real(8)       :: phot_g_mean_mag = 100
            character(16) :: phot_variable_flag = "not found"
            real(8)       :: l = 420
            real(8)       :: b = 420
            real(8)       :: ecl_lon = 420
            real(8)       :: ecl_lat = 420
        end type GaiaEntry

        type TgasEntry
            character(16) :: hip = 'no'
            character(32) :: tycho2_id = 'no'
            integer(8)    :: solution_id = 0
            integer(8)    :: source_id = 0
            integer(8)    :: random_index = 0
            real(8)       :: ref_epoch = 0
            real(8)       :: ra = 420
            real(8)       :: ra_error = 0
            real(8)       :: dec = 420
            real(8)       :: dec_error = 0
            real(8)       :: parallax = 0
            real(8)       :: parallax_error = 0
            real(8)       :: pmra = 1296000000
            real(8)       :: pmra_error = 0
            real(8)       :: pmdec = 1296000000
            real(8)       :: pmdec_error = 0
            real(4)       :: ra_dec_corr = 2
            real(4)       :: ra_parallax_corr = 2
            real(4)       :: ra_pmra_corr = 2
            real(4)       :: ra_pmdec_corr = 2
            real(4)       :: dec_parallax_corr = 2
            real(4)       :: dec_pmra_corr = 2
            real(4)       :: dec_pmdec_corr = 2
            real(4)       :: parallax_pmra_corr = 2
            real(4)       :: parallax_pmdec_corr = 2
            real(4)       :: pmra_pmdec_corr = 2
            integer(4)    :: astrometric_n_obs_al = -1
            integer(4)    :: astrometric_n_obs_ac = -1
            integer(4)    :: astrometric_n_good_obs_al = -1
            integer(4)    :: astrometric_n_good_obs_ac = -1
            integer(4)    :: astrometric_n_bad_obs_al = -1
            integer(4)    :: astrometric_n_bad_obs_ac = -1
            real(4)       :: astrometric_delta_q = -1
            real(8)       :: astrometric_excess_noise = 1296000000
            real(8)       :: astrometric_excess_noise_sig = -1
            character(16) :: astrometric_primary_flag = "false"
            real(4)       :: astrometric_relegation_factor = -1
            real(4)       :: astrometric_weight_al = -1
            real(4)       :: astrometric_weight_ac = -1
            integer(4)    :: astrometric_priors_used = 0
            integer(2)    :: matched_observations = -1
            character(16) :: duplicated_source = "false"
            real(4)       :: scan_direction_strength_k1 = 2
            real(4)       :: scan_direction_strength_k2 = 2
            real(4)       :: scan_direction_strength_k3 = 2
            real(4)       :: scan_direction_strength_k4 = 2
            real(4)       :: scan_direction_mean_k1 = 270
            real(4)       :: scan_direction_mean_k2 = 270
            real(4)       :: scan_direction_mean_k3 = 270
            real(4)       :: scan_direction_mean_k4 = 270
            integer(4)    :: phot_g_n_obs = -1
            real(8)       :: phot_g_mean_flux = -1
            real(8)       :: phot_g_mean_flux_error = -1
            real(8)       :: phot_g_mean_mag = 100
            character(16) :: phot_variable_flag = "not found"
            real(8)       :: l = 420
            real(8)       :: b = 420
            real(8)       :: ecl_lon = 420
            real(8)       :: ecl_lat = 420
        end type TgasEntry
    
        type DataString
            character(1024) :: threadString
        end type DataString
    
        type MedTable 
            integer :: ihealp = 0
            real(8) :: lcen = 0
            real(8) :: bcen = 0
            real(8) :: l = 0
            real(8) :: b = 0
            real(8) :: pml = 0
            real(8) :: pmb = 0
            real(8) :: parallax = 0
            integer :: nhealp = 0
            real(8) :: ra = 0
            real(8) :: dec = 0
        end type MedTable
        
    end module CustomDataTypes