
declare @cmd varchar(8000)
set @cmd = 'bcp "select '+
's.source_id,designation,ra,dec,ra_error,dec_error,ra_dec_corr,pmra,pmdec,phot_g_mean_mag,astrometric_chi2_al,astrometric_n_good_obs_al,'+
'astrometric_excess_noise,astrometric_excess_noise_sig,astrometric_matched_transits,'+
'phot_rp_mean_mag,phot_bp_mean_mag,ruwe,phot_bp_rp_excess_factor,bp_rp,'+
'phot_bp_mean_flux_over_error,phot_g_mean_flux_over_error,phot_rp_mean_flux_over_error,'+
'phot_bp_n_obs,phot_g_n_obs,phot_rp_n_obs,'+
'phot_bp_n_blended_transits,phot_rp_n_blended_transits,astrometric_params_solved,b,phot_bp_mean_flux_error,phot_g_mean_flux_error,phot_rp_mean_flux_error '+
'from gaia_source as s " ' +
' queryout G:\share\mar\gaiaSkinnyishDR3.csv -c '+
' -S ramses25 -d GAIAEDR3 -U sa -P XXXXX -t , -o G:\share\mar\gaiaDR3.log '
select @cmd
exec xp_cmdshell @cmd
