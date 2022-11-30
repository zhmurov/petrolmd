#!/bin/bash
GMX=/usr/local/gromacs/bin/gmx

for filename in *.pdb; do
    name=$(basename $filename .pdb)
    # Run GROMACS to create top files for all the PDBs in the folder
    $GMX pdb2gmx -f $filename -o ${name}.gro -p ${name}.top -i ${name}_posre.itp -ff trappeua -water tip4p
    # Create a copy of the topology that can be included
    cp ${name}.top ${name}.itp
    # Remove th1e header
    sed -i -n '/\[ moleculetype \]/,$p' ${name}.itp
    # Remove pairs section (needed by TraPPE forcefield)
    if grep -Fxq "[ pairs ]" ${name}.itp
    then
        sed -i -n '1,/pairs/p;/angles/,$p' ${name}.itp
        sed -i '\[ pairs \]/d' C8H18.itp
    else
        echo "Skipping ${name}.itp"
    fi
    # Remove the footer
    sed -i '/; Include Position restraint file/,$d' ${name}.itp
    # Rename the molecule
    sed -i "s/Other/${name}/g" ${name}.itp
done