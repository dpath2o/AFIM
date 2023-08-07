#! /bin/csh -f
### qsubs -I -q copyq -l walltime=10:00:00 -l mem=200GB
############################################################################

cd /g/data/cj50/access-om2/raw-output/access-om2-01/01deg_jra55v140_iaf/
pwd
ls output188/ocean/ocean-3d-u-1-daily-mean-ym_2005*
set yr = 2005
set sim = 188
set sim2 = 189
set sim3 = 190
set sim4 = 191

echo $yr $sim $sim2 $sim3 $sim4

while ( $yr < 2007)
	
	echo $yr $sim $sim2 $sim3 $sim4
	
	echo ncrcat output${sim}/ocean/ocean-3d-u-1-daily-mean-ym_${yr}\* output${sim2}/ocean/ocean-3d-u-1-daily-mean-ym_${yr}\* output${sim3}/ocean/ocean-3d-u-1-daily-mean-ym_${yr}\* output${sim4}/ocean/ocean-3d-u-1-daily-mean-ym_${yr}\* /g/data/jk72/pas561/for_dan/ocean-3d-u-1-daily-mean-ym_${yr}.nc

	ncrcat -d st_ocean,0,1 output${sim}/ocean/ocean-3d-u-1-daily-mean-ym_${yr}* output${sim2}/ocean/ocean-3d-u-1-daily-mean-ym_${yr}* output${sim3}/ocean/ocean-3d-u-1-daily-mean-ym_${yr}* output${sim4}/ocean/ocean-3d-u-1-daily-mean-ym_${yr}* /g/data/jk72/pas561/for_dan/ocean-3d-u-1-daily-mean-ym_${yr}.nc
	
        @ sim =$sim4 + 1
	@ sim2 = $sim4 + 2
	@ sim3 = $sim4 + 3
	@ sim4 = $sim4 + 4

	@ yr = $yr + 1

	#echo $yr $sim $sim2 $sim3 $sim4

end


############################################################################
