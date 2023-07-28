#/bin/bash

# In this example, the atm virtural environment was install at /projects2/insite/eric.chen/anaconda3/envs/atm/. You will need to change it to your own installation path 

# The $atom_build_path is the path to your local AToM_OpenMM repo

if [ ! -f complex_0.xml ]; then
 
   /projects2/insite/eric.chen/anaconda3/envs/atm/bin/python mintherm_atm.py || exit
   /projects2/insite/eric.chen/anaconda3/envs/atm/bin/python mdlambda_atm.py || exit
   /projects2/insite/eric.chen/anaconda3/envs/atm/bin/python mdlambda_equil_atm.py || exit

fi


# if you only have one GPU.
#echo -e "localhost,0:0,1,CUDA,,/tmp" > $work_dir/nodefile

#if you have 4 GPUs
echo -e "localhost,0:0,1,CUDA,,/tmp\nlocalhost,0:1,1,CUDA,,/tmp\nlocalhost,0:2,1,CUDA,,/tmp\nlocalhost,0:3,1,CUDA,,/tmp"  > nodefile


/projects2/insite/eric.chen/anaconda3/envs/atm/bin/python $atom_build_path/rbfe_explicit.py atom.cntl || exit

