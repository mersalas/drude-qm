# drude-qm
MD and QM files for the drude paper

Paper: [Length of the methylene side chain in lysine modulates non-covalent interactions in histidine-rich peptides]()

![Drude MD Workflow]()

## Files
- [drude_setup](): instructions on how to setup Drude MD input files and run simulations
- [drude_ana](): Gromacs scripts for MD analysis and ORCA scripts to run QM calculations
- [1_gromacs]():  
    - [coord](): initial structures of lipopeptides in pdb format
    - [paramtr](): mdp files
- [2_ffparam](): input & output files for Drude parameterization 
- [3_drude](): 
    - [coord](): processed pdb & psf files
    - [scripts](): Drude MD scripts 
- [4_openmm](): scripts to run Drude MD in OpenMM
- [5_ana](): inputs files & scripts to analyse trajectory files & calculate free energy surface
- [6_qm](): sample input files & script to perform QM calculations 

## How to cite:
```bibtex
@unpublished{Salas2025,
author = {Remmer L. Salas and Portia Mahal G. Sabido and Ricky B. Nellas},
title = {Length of the methylene side chain in lysine modulates non-covalent interactions in histidine-rich peptides},
journal = {},
volume = {},
pages = {},
year = {2026},
doi = {},
}