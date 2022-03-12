#!/bin/bash
# -------------------------
# run python plotting script passing following parameters
# date, variable, pressure level, geographical domain
# -------------------------
# Usage: ./run_plot_multiple_timesteps.sh

## load modules
#source /etc/profile.d/modules.sh
#module 'load dyn_tools'
#module 'load miniconda3'
#source activate iacpy3_2019

# parameters
startdate=20200207_00
enddate=20200207_12 #20200208_00
VarVal=("virtpotT" "T" "QV" "dD" "dexc")
LevVal=(750. 850. 900. 950.)
DomVal=("domain3") #("domain1" "domain2")

# create plot for different parameters
currdate=$startdate
while [ "${currdate}" != "${enddate}" ];do
	for var in ${VarVal[@]}; do
	    for lev in ${LevVal[@]};do
		    for domain in ${DomVal[@]};do
			    echo ''
			    echo 'creating plot for' $date $var $lev $domain
			    # run plotting script
			    python plot_map_echam_cosmo10km_5km_1km.py $currdate $var $lev $domain
		done
	    done
	done
        # update date (add n hours)
        currdate=`/usr/local/bin/newtime ${currdate} +6`
done
