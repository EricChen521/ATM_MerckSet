Description
===========

The structural files and scripts to launch ATM calculation on [Schindler](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.0c00900) et al. dataset. 

The repository is organized as following:

Under the **targets** directory, the repo includes 8 different target systems: **cdk8**, **cmet**, **eg5**, **hif2a**, **pfkfb3**, **shp2**, **syk**, **tnks2**.

Inside each target, the repo provides:

* free_energy/_pairname_
	+ atom.cntl (ATM control file).
	+ mintherm_atm.py/mdlambda_equil_atm.py/mdlambda_atm.py (ATM structural equilibration protocol).
	+ complex.prmtop/complex.inpcrd (molecular system for each pair).
	+ nodefile (file that specifies available GPU devices).

Under the **scripts** directroy, the repo provides two useful scripts:

* run_atm.sh (script to run ATM calculation for each RBFE pair).
* uwhat_analysis.R (script to analyze perturbation with UWHAT).  
	

 


