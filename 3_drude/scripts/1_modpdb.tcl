# Load the original PDB
mol load pdb ../2_gromacs/3_mph30_npt.pdb

# Rename "PDB general atom name" to "CHARMM-specific atom name"
# HISE => HSE (protonated); SOL => TIP3; CL => CLA
[atomselect top "resname HISE"] set resname HSE
[atomselect top "resname SOL"] set resname TIP3
[atomselect top "resname CL"] set resname CLA

# OW, HW1, HW2 => OH2, H1, H2
[atomselect top "name OW"] set name OH2
[atomselect top "name HW1"] set name H1
[atomselect top "name HW2"] set name H2

# Select peptide and TFE
set pep [atomselect top "protein"]
set tfe [atomselect top "resname TFE"]
set tip [atomselect top "resname TIP3"]
set cl [atomselect top "resname CLA"]

# Combine peptide & TFE 
set sys [atomselect top "protein or resname TFE or resname TIP3 or resname CLA"]

# Write the modified PDB of the selected atoms (peptide & TFE)
$sys writepdb 1_mph30_solv.pdb
$pep writepdb 1_mph.pdb
$tfe writepdb 1_tfe.pdb
$tip writepdb 1_tip3.pdb
$cl writepdb 1_cl.pdb

exit
