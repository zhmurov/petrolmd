#!/bin/bash
GMX=/usr/local/gromacs/bin/gmx

for filename in *.pdb; do
    name=$(basename $filename .pdb)
    # Run GROMACS to create top files for all the PDBs in the folder
    $GMX pdb2gmx -f $filename -o ${name}.gro -p ${name}.top -i ${name}_posre.itp <<< $'16\n3\n'
    # Remove data from `[ pairs ]` section as it is required by the TRAPPE-UA forcefield
    sed -n '1,/pairs/p;/angles/,$p' ${name}.top >> ${name}_nopairs.top
done