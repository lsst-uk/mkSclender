
declare @cmd varchar(8000)
set @cmd = 'bcp "select '+
'pts_key,''2MASS J''+designation as designation,ra,decl,err_maj,err_min,j_m,j_cmsig,h_m,h_cmsig,k_m,k_cmsig,bl_flg,rd_flg,cc_flg,'+
'ph_qual,jdate,j_snr,h_snr,k_snr '+
'from twomass_psc where ext_key is null " ' +
' queryout G:\share\mar\2MASS_PSC_all_SNR.csv -c '+
' -S ramses27 -d TWOMASS2 -U sa -P XXXXX -t , -o G:\share\mar\2MASS_PSC.log '
select @cmd
exec xp_cmdshell @cmd

declare @cmd varchar(8000)
set @cmd = 'bcp "select '+
'ext_key,''2MASX J''+designation as designation,sup_ra,sup_dec,0.3,0.3,j_m_k20fe,j_msig_k20fe,j_flg_k20fe,h_m_k20fe,h_flg_k20fe,h_msig_k20fe,k_m_k20fe,k_flg_k20fe,k_msig_k20fe,'+
'cc_flg,jdate '+
'from twomass_xsc " ' +
' queryout G:\share\mar\2MASS_XSC_all_SNR.csv -c '+
' -S ramses27 -d TWOMASS2 -U sa -P XXXXXX -t , -o G:\share\mar\2MASS_XSC.log '
select @cmd
exec xp_cmdshell @cmd
