#!/bin/csh
#
# Production

# The OpenMM check point file (.chk) cannot be used in a different machine environment.
# So please make sure if you are using the same GPU and CUDA version of machine while doing additional
# production steps with the check point file.

set init1 = step3_charmm2omm
set init2 = step4_equilibration

set cnt = 1
set cntmax = 2

set input = step5_production

while ( ${cnt} <= ${cntmax} )
    @ pcnt = ${cnt} - 1
    set istep = step5_${cnt}
    set pstep = step5_${pcnt}

    if ( ${cnt} == 1 ) set pstep = step4_equilibration

    python -u openmm_run.py -i ${input}.inp -t toppar.str -p ${init1}.psf -c ${init2}.crd -irst ${pstep}.rst -orst ${istep}.rst -odcd ${istep}.dcd > ${istep}.out
    @ cnt += 1
end
