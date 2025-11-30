# Load psfgen-plugin and CHARMM topology files
package require psfgen
resetpsf
topology ../../toppar/top_all36_prot.rtf
topology ../../toppar/toppar_all36_prot_modify_res.str
topology ../../toppar/top_all36_cgenff.rtf
topology ../../toppar/par_all36_lipid.prm
topology ../../toppar/toppar_water_ions.str

# Define the peptide segment "PROA"
segment PROA {
    pdb 1_mph.pdb
    first MYRN   ;# Patch N-terminus with MYRN
    last CT2     ;# Patch C-terminus with CT2
    #auto none
}

# Assign the coordinates of atoms to PROA using the PDB data
coordpdb 1_mph.pdb PROA

# Define the TFE segment
segment TFE {
    pdb 1_tfe.pdb
    #auto none
}

# Assign the coordinates for TFE
coordpdb 1_tfe.pdb TFE


# Define the SOLV segment
segment SOLV {
    pdb 1_tip3.pdb
    #auto none
}

# Assign the coordinates for SOLV
coordpdb 1_tip3.pdb SOLV

# Define the CLA segment
segment CLA {
    pdb 1_cl.pdb
    #auto none
}

# Assign the coordinates for CLA
coordpdb 1_cl.pdb CLA

# Guess the coordinates of missing atoms (mainly hydrogen)
#guesscoord 

# Generate PDB and PSF files
writepdb 2_mph30.pdb
writepsf 2_mph30.psf
exit

