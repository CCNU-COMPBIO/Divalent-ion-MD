# parameter files for running simulations using AMBER forcefield and TIP4PEW water model.

## [Nucleosome_models](With tail and without tails)
Different nucleosome models used for simulations

## gen_nucl.leap 

Script to prepare the residue library and force field parameters for use with LEaP progarm.

To execuate the script:

tleap -f gen_nucl

## parmed_1264_na.in

Script to add 12-6-4 parmerater for divalent ions. 

To execuate the script:

parmed -i parmed_1264_na.in

## Min.in, Equil_v.in, Equil_pt.in, Prod.in

Simulation configration files for energy minimizations, equilibration, and production run.

## mist_sub_job_GPU.sh

Sample scripts to run the simulations on GPU.

## Required Progarms
Amber (https://ambermd.org)

Simulations were performed with the version of Amber22.

## MD trajectories

MD Trajectories are archived at: 
https://doi.org/10.5281/zenodo.4771269

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4771269.svg)](https://doi.org/10.5281/zenodo.4771269)
