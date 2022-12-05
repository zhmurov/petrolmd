#!/bin/bash
GMX=/usr/local/gromacs/bin/gmx

for filename in *.pdb; do
    name=$(basename $filename .pdb)
    # Run GROMACS to create top files for all the PDBs in the folder
    $GMX pdb2gmx -f $filename -o ${name}.gro -p ${name}.top -i ${name}_posre.itp -ff charmm36-jul2021 -water tip3p
    # Create a copy of the topology that can be included
    cp ${name}.top ${name}.itp
    # Remove th1e header
    sed -i -n '/\[ moleculetype \]/,$p' ${name}.itp
    # Remove the footer
    sed -i '/; Include Position restraint file/,$d' ${name}.itp
    # Rename the molecule
    sed -i "s/Other/${name}/g" ${name}.itp
    # Combine topologies into one file
    # cat ${name}.itp >> alkanes.itp
done