
source /Users/ericchen/miniconda3/bin/activate atm_analysis

target_name=$1
python run_diffnet.py --csv atm_results.dat --atm-flag True --ref-file ref.dat --new-csv

python bfe_correlation.py --input_csv atm_diffnet_results.csv --x EXP_dG --y Predicted_dG --y_error Predicted_Error --title ${target_name}
