#/bin/bash

# In this example, the atm virtural environment was install at /projects2/insite/eric.chen/anaconda3/envs/atm/. You will need to change it to your own installation path 

# The $atom_build_path is the path to your local AToM_OpenMM repo

 
/projects2/insite/eric.chen/anaconda3/envs/atm/bin/python mintherm_atm.py || exit
/projects2/insite/eric.chen/anaconda3/envs/atm/bin/python mdlambda_atm.py || exit
/projects2/insite/eric.chen/anaconda3/envs/atm/bin/python mdlambda_equil_atm.py || exit



/projects2/insite/eric.chen/anaconda3/envs/atm/bin/python $atom_build_path/rbfe_explicit.py atom.cntl || exit

