declare @cmd varchar(8000)
set @cmd='bcp "select cntr,source_name,ra_pm, dec_pm, sigra_pm, sigdec_pm, sigradec_pm,ab_flags,w1mpro, w1sigmpro, w1sat,w1snr, w2mpro, w2sigmpro, w2sat, w2snr from catwise_2020   option (maxdop 1)"
  queryout G\:\share\mar\catwiseSkinnyishPM.csv -c  -S ramses25 -d WISE -U sa -P XXXXXX -t , -o G:\share\mar\catwisePM.log'

select @cmd
exec xp_cmdshell @cmd
