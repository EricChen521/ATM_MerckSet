#!/bin/bash

free_energy_dir=$1
start_frame_index=$2

SCRIPT_PATH=`realpath $0`
SCRIPT_DIR=`dirname $SCRIPT_PATH`
UWHAT_SCRIPT=$SCRIPT_DIR/uwham_analysis.R

cd $free_energy_dir
echo "name0,name1,mean,error,start_frame,end_frame"
for pair in *
do
	if [ -d $pair ]; then
		cd $pair
		R CMD BATCH -complex -$start_frame_index -10000 $UWHAT_SCRIPT
		result=($(grep -r "^DDGb =" uwham_analysis.Rout))
		name1=${pair%%~*}
		name2=${pair##*~}
		ddG=$(echo ${result[2]} | xargs printf "%.2f")
		uncertanity=$(echo ${result[4]} | xargs printf "%.2f")
		end_frame=$(echo ${result[7]})
		echo "$name1,$name2,$ddG,$uncertanity,$start_frame_index,$end_frame"
		cd ..
	fi
done
