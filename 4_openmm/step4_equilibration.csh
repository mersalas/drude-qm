#!/bin/csh
#
# Equilibration

set init = step3_charmm2omm
set istep = step4_equilibration

python -u openmm_run.py -i ${istep}.inp -t toppar.str -p ${init}.psf -c ${init}.crd -orst ${istep}.rst -odcd ${istep}.dcd > ${istep}.out
