Description
===========

The structural files and scripts to launch ATM calculation on [Schindler](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.0c00900) et al. dataset. 

The repository is organized as following:

Under the **targets** directory, the repo includes 8 different target systems: **cdk8**, **cmet**, **eg5**, **hif2a**, **pfkfb3**, **shp2**, **syk**, **tnks2**.

Inside each target, the repo provides:

* protein/protein.pdb (Prepared receptor structure after water placement).
* ligands/_ligandname_/ligand.sdf (Prepared ligand structure for each ligand).
* forcefield/_ligandname_: FFEngine forcefield parametrization for each ligand.
* free_energy/_pairname_
	+ atom.cntl (ATM control file).
	+ mintherm_atm.py/mdlambda_equil_atm.py/mdlambda_atm.py (ATM structural equilibration protocol).
	+ complex.prmtop/complex.inpcrd (ATM system for each pair).
	+ nodefile (Available GPU devices for ATM).

Under the **scripts** directroy, the repo provides useful scripts to launch ATM claculation and analyze ATM results:  
	
* run_atm.sh (To launch ATM calculation for the RBFE pair).
* run_uwhat.sh (To launch UWHAT analysis on trajectory energies to get ddG).
* run_diffnet.py (To do DiffNet analysis on ddG to get dG).
* bfe_correlation.py (Plot correlation between experimental and predicted dG).

 


